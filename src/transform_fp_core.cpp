#include <Rcpp.h>
#include <cmath>
#include <string>

using namespace Rcpp;


// Internal reusable fractional-polynomial transformation kernel.
//
// This function contains the actual numerical FP transformation logic. It is
// deliberately kept as a private C++ helper rather than an R-callable function:
// both exported Rcpp wrappers below call this same implementation, so the FP
// rules are written in exactly one place.
//
// R-level preprocessing is still responsible for:
//   * removing/sorting invalid or missing powers as needed,
//   * bypassing binary variables when check_binary = TRUE,
//   * supplying default shift/scale values,
//   * forcing shift = 0 when zero = TRUE,
//   * user-facing input validation,
//   * output column names.
//
// The FP rules implemented here are:
//   * p = 0      -> log(x)
//   * p != 0     -> x^p
//   * repeated p -> previous FP column multiplied by log(x)
//
// Example:
//   powers c(1, 1, 1) produce x, x * log(x), x * log(x)^2.
static NumericMatrix transform_fp_core_internal(const NumericVector& x_raw,
                                                const NumericVector& power,
                                                const double shift_val,
                                                const double scale_val,
                                                const bool zero) {
  
  const int n = x_raw.size();
  const int m = power.size();
  
  // An empty power vector gives an n x 0 matrix. The R wrapper normally avoids
  // this case, but returning a valid empty matrix makes the C++ helper safe.
  if (m == 0) {
    return NumericMatrix(n, 0);
  }
  
  // Guard the division below. This is intentionally kept in C++ even if the R
  // wrapper validates scale, because this kernel may be reused by other C++
  // helpers and should never silently divide by zero.
  if (!R_finite(scale_val) || scale_val == 0.0) {
    stop("`scale_val` must be finite and non-zero.");
  }
  
  NumericMatrix out(n, m);
  
  const double inv_scale = 1.0 / scale_val;
  
  for (int i = 0; i < n; ++i) {
    
    // In zero mode, non-positive raw values are assigned zero-valued
    // transformation columns. NumericMatrix initializes to zero, so continuing
    // here leaves the whole row at zero.
    if (zero && x_raw[i] <= 0.0) {
      continue;
    }
    
    // Apply the shift/scale convention used by transform_vector_fp(). For the
    // batch generator below, x has already been shifted/scaled in R, so that
    // caller passes shift_val = 0 and scale_val = 1.
    const double xi = (x_raw[i] + shift_val) * inv_scale;
    
    // In zero mode, values that remain non-positive after shift/scale are also
    // assigned zero-valued transformation columns.
    if (zero && xi <= 0.0) {
      continue;
    }
    
    const double log_xi = std::log(xi);
    
    // First FP column. Power 0 represents log(x); every other power represents
    // x^power.
    if (power[0] == 0.0) {
      out(i, 0) = log_xi;
    } else {
      out(i, 0) = std::pow(xi, power[0]);
    }
    
    // Remaining FP columns. Repeated powers follow the standard FP convention:
    // reuse the previous column and multiply by one additional log(x) factor.
    for (int j = 1; j < m; ++j) {
      if (power[j] == power[j - 1]) {
        out(i, j) = out(i, j - 1) * log_xi;
      } else {
        if (power[j] == 0.0) {
          out(i, j) = log_xi;
        } else {
          out(i, j) = std::pow(xi, power[j]);
        }
      }
    }
  }
  
  return out;
}


// R-callable wrapper used by the R transform_vector_fp() function.
//
// The Rcpp export creates an internal R wrapper in R/RcppExports.R. It does not
// make the function part of the public package API unless it is also exported
// from NAMESPACE. 
//
// [[Rcpp::export]]
NumericMatrix transform_fp_core(NumericVector x_raw,
                                NumericVector power,
                                double shift_val,
                                double scale_val,
                                bool zero) {
  return transform_fp_core_internal(
    x_raw,
    power,
    shift_val,
    scale_val,
    zero
  );
}


// Return TRUE if x has at most two distinct values.
//
// This helper mirrors the old R-side binary shortcut:
//
//   length(unique(x)) <= 2L
//
// That means missing values are counted as one distinct level. For example,
// c(0, 1) is binary, c(0, NA) is binary, but c(0, 1, NA) is not binary under
// the old R rule.
static bool is_binary_vector_cpp(const NumericVector& x) {
  const int n = x.size();
  
  bool seen_missing = false;
  bool seen_first = false;
  bool seen_second = false;
  
  double first = 0.0;
  double second = 0.0;
  
  int n_distinct = 0;
  
  for (int i = 0; i < n; ++i) {
    const double value = x[i];
    
    // Match unique(x): missing values count as one distinct value.
    if (NumericVector::is_na(value)) {
      if (!seen_missing) {
        seen_missing = true;
        ++n_distinct;
        
        if (n_distinct > 2) {
          return false;
        }
      }
      
      continue;
    }
    
    // Existing first non-missing value.
    if (seen_first && value == first) {
      continue;
    }
    
    // Existing second non-missing value.
    if (seen_second && value == second) {
      continue;
    }
    
    // New first non-missing value.
    if (!seen_first) {
      first = value;
      seen_first = true;
      ++n_distinct;
      
      if (n_distinct > 2) {
        return false;
      }
      
      continue;
    }
    
    // New second non-missing value.
    if (!seen_second) {
      second = value;
      seen_second = true;
      ++n_distinct;
      
      if (n_distinct > 2) {
        return false;
      }
      
      continue;
    }
    
    // Any additional distinct value means x is not binary.
    return false;
  }
  
  return true;
}


// Generate all requested FP transformation matrices for one variable.
//
// `powers` is the matrix returned by generate_powers_fp() in R. Each row is one
// candidate FP power combination. The return value is a list with one matrix per
// candidate row, matching the old R implementation conceptually:
//
//   lapply(seq_len(nrow(powers)), function(i) {
//     transform_vector_fp(x, power = powers[i, ], zero = zero)
//   })
//
// If `catzero` is provided, it is prepended as the first column of each returned
// matrix, matching the old R code:
//
//   mat <- cbind(catzero, mat)
//   colnames(mat)[1] <- "catzero"
//
// The speed gain comes from batching all candidate transformations into one
// R-to-C++ call instead of making one R call to transform_vector_fp() for every
// candidate power row.
//
// [[Rcpp::export]]
List generate_transformations_fp_cpp(const NumericVector& x,
                                     const NumericMatrix& powers,
                                     const bool zero,
                                     Nullable<NumericMatrix> catzero = R_NilValue) {
  
  const int n = x.size();
  const int n_candidates = powers.nrow();
  const int degree = powers.ncol();
  
  // `catzero` is optional. When present, it must be an n x 1 numeric matrix and
  // is copied as the first column of every candidate matrix.
  const bool use_catzero = catzero.isNotNull();
  
  NumericMatrix catzero_mat;
  
  if (use_catzero) {
    catzero_mat = NumericMatrix(catzero);
    
    if (catzero_mat.nrow() != n) {
      stop("`catzero` must have one row per observation in `x`.");
    }
    
    if (catzero_mat.ncol() != 1) {
      stop("`catzero` must be an n x 1 matrix.");
    }
  }
  
  // Match transform_vector_fp(..., check_binary = TRUE): if x is binary and
  // zero = FALSE, return x unchanged rather than applying FP powers. When
  // zero = TRUE, the zero/catzero logic is active and FP transformation should
  // still be applied to the positive part.
  const bool binary_x = (!zero && is_binary_vector_cpp(x));
  
  // Preallocate the output list. Each element is the design matrix for one
  // candidate power row.
  List out(n_candidates);
  
  for (int i = 0; i < n_candidates; ++i) {
    NumericMatrix fp_mat;
    
    if (binary_x) {
      // Binary variables are not FP-transformed. The same one-column matrix is
      // returned for every candidate row, matching transform_vector_fp().
      fp_mat = NumericMatrix(n, 1);
      
      for (int r = 0; r < n; ++r) {
        fp_mat(r, 0) = x[r];
      }
      
    } else {
      // Extract the i-th candidate power row as a vector for the shared FP
      // kernel. The candidate powers are generated and ordered by R before
      // entering this function.
      NumericVector candidate_power(degree);
      
      for (int j = 0; j < degree; ++j) {
        candidate_power[j] = powers(i, j);
      }
      
      // x is already shifted/scaled at the R level for ordinary FP candidate
      // generation. Therefore this batch helper uses shift = 0 and scale = 1,
      // matching the previous call:
      //   transform_vector_fp(x = x, power = powers[i, ], zero = zero)
      fp_mat = transform_fp_core_internal(
        x,
        candidate_power,
        0.0,
        1.0,
        zero
      );
    }
    
    if (!use_catzero) {
      out[i] = fp_mat;
      continue;
    }
    
    // Prepend catzero to the candidate FP matrix.
    const int fp_cols = fp_mat.ncol();
    
    NumericMatrix fp_mat_with_catzero(n, fp_cols + 1);
    
    // First column: structural-zero indicator.
    for (int r = 0; r < n; ++r) {
      fp_mat_with_catzero(r, 0) = catzero_mat(r, 0);
    }
    
    // Remaining columns: FP-transformed candidate terms.
    for (int c = 0; c < fp_cols; ++c) {
      for (int r = 0; r < n; ++r) {
        fp_mat_with_catzero(r, c + 1) = fp_mat(r, c);
      }
    }
    
    // Assign stable column names. The first column is the structural-zero
    // indicator; the FP columns are named V1, V2, ..., so the output has explicit
    // and predictable column names.
    CharacterVector col_names(fp_cols + 1);
    
    col_names[0] = "catzero";
    
    for (int c = 0; c < fp_cols; ++c) {
      col_names[c + 1] = "V" + std::to_string(c + 1);
    }
    
    colnames(fp_mat_with_catzero) = col_names;
    
    out[i] = fp_mat_with_catzero;
  }
  
  return out;
}
