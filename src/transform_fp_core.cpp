#include <Rcpp.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>

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
    
    // Preserve R-style missing/non-finite propagation.
    //
    // Without this guard, std::log() and std::pow() may turn R's NA_REAL into
    // NaN-like values at the C++ level. The R implementation propagates missing
    // values through the transformed columns, so the C++ kernel should explicitly
    // fill the whole output row with the original non-finite value.
    //
    // In normal fitting paths, x should already have been checked for missing and
    // non-finite values. This branch is nevertheless needed because this kernel is
    // also used by lower-level helpers and should preserve transform_vector_fp()
    // semantics when called directly.
    if (!R_finite(x_raw[i])) {
      for (int j = 0; j < m; ++j) {
        out(i, j) = x_raw[i];
      }
      continue;
    }
    
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
    
    // Preserve missing/non-finite values introduced by shift/scale.
    if (!R_finite(xi)) {
      for (int j = 0; j < m; ++j) {
        out(i, j) = xi;
      }
      continue;
    }
    
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


// CODES FOR BUILDING ADJUSTMENT VARIABLE STEP IN TRANSFORM_DATA_STEP()
// -----------------------------------------------------------------------------
// build_adjustment_step() C++ utility helpers ---------------------------------
// -----------------------------------------------------------------------------
//
// These helpers support the future C++ replacement of the hot per-variable loop
// inside build_adjustment_step().
//
// Scope:
//   * These helpers are internal C++ utilities.
//   * They do not change package behaviour by themselves.
//   * They assume higher-level R code has already validated/prepared inputs.
//   * They support FP transforms and stored-parameter ACD application.
//   * They do NOT fit/estimate ACD parameters.
//
// Important build_adjustment_step() invariant:
//   ACD adjustment variables should already have stored acd_parameter values.
//   Therefore this C++ path only needs to apply stored ACD parameters.


// Create an n x 0 numeric matrix.
//
// Used when an adjustment variable contributes no columns, for example when
// the variable is eliminated because all selected powers are NA.
static NumericMatrix mfp2_empty_matrix_cpp(const int n) {
  return NumericMatrix(n, 0);
}


// Convert a numeric vector into an n x 1 numeric matrix.
//
// This mirrors the R-side normalization:
//   if (is.null(dim(transformed))) matrix(transformed, ncol = 1L)
static NumericMatrix mfp2_vector_to_matrix_cpp(const NumericVector& x) {
  const int n = x.size();
  NumericMatrix out(n, 1);
  
  for (int i = 0; i < n; ++i) {
    out(i, 0) = x[i];
  }
  
  return out;
}


// Column-bind two numeric matrices.
//
// This is used for:
//   * cbind(catzero, transformed)
//   * cbind(x_fp, x_acd)
// It is intentionally small and explicit to avoid R-level cbind overhead in
// the hot adjustment-building path.
static NumericMatrix mfp2_cbind_two_cpp(const NumericMatrix& a,
                                        const NumericMatrix& b) {
  const int n = a.nrow();
  
  if (b.nrow() != n) {
    stop("Internal error: matrices to cbind have incompatible row counts.");
  }
  
  const int pa = a.ncol();
  const int pb = b.ncol();
  
  NumericMatrix out(n, pa + pb);
  
  for (int j = 0; j < pa; ++j) {
    for (int i = 0; i < n; ++i) {
      out(i, j) = a(i, j);
    }
  }
  
  for (int j = 0; j < pb; ++j) {
    for (int i = 0; i < n; ++i) {
      out(i, pa + j) = b(i, j);
    }
  }
  
  return out;
}


// Compare two normalized numeric power keys.
//
// This replaces the R-side:
//   identical(prev_power_key, current_power_key)
//
// for numeric vectors/lists that have already been normalized by
// normalize_powers_for_convergence() in R.
//
// NA values are considered equal only when both sides are NA in the same
// position. Length must match exactly.
static bool mfp2_numeric_keys_identical_cpp(SEXP a, SEXP b) {
  if (Rf_isNull(a) && Rf_isNull(b)) {
    return true;
  }
  
  if (Rf_isNull(a) || Rf_isNull(b)) {
    return false;
  }
  
  NumericVector x(a);
  NumericVector y(b);
  
  if (x.size() != y.size()) {
    return false;
  }
  
  for (int i = 0; i < x.size(); ++i) {
    const bool x_na = NumericVector::is_na(x[i]);
    const bool y_na = NumericVector::is_na(y[i]);
    
    if (x_na || y_na) {
      if (!(x_na && y_na)) {
        return false;
      }
      continue;
    }
    
    if (x[i] != y[i]) {
      return false;
    }
  }
  
  return true;
}


// Remove NA powers and sort finite/non-missing powers.
//
// This mirrors the behaviour needed before calling the FP core. In R,
// transform_vector_fp() effectively works with the non-missing selected powers.
// Eliminated variables are handled before this helper is called.
static NumericVector mfp2_clean_sorted_power_cpp(SEXP power_sexp) {
  if (Rf_isNull(power_sexp)) {
    return NumericVector(0);
  }
  
  NumericVector power(power_sexp);
  std::vector<double> clean;
  clean.reserve(power.size());
  
  for (int i = 0; i < power.size(); ++i) {
    if (!NumericVector::is_na(power[i])) {
      clean.push_back(power[i]);
    }
  }
  
  std::sort(clean.begin(), clean.end());
  
  NumericVector out(clean.size());
  
  for (int i = 0; i < static_cast<int>(clean.size()); ++i) {
    out[i] = clean[i];
  }
  
  return out;
}


// Transform one adjustment variable using ordinary FP powers.
//
// This supports the non-ACD branch of build_adjustment_step():
//   transform_vector_fp(x = xvec, power = power_current, zero = zero)
//
// It returns an n x k matrix, or n x 0 when there are no non-missing powers.
static NumericMatrix mfp2_transform_fp_adjustment_cpp(const NumericVector& x,
                                                      SEXP power_sexp,
                                                      const bool zero) {
  const int n = x.size();
  NumericVector power = mfp2_clean_sorted_power_cpp(power_sexp);
  
  if (power.size() == 0) {
    return mfp2_empty_matrix_cpp(n);
  }
  
  // Match transform_vector_fp(..., check_binary = TRUE) behaviour.
  // Binary variables are kept as a single untransformed column when zero
  // handling is not active.
  if (!zero && is_binary_vector_cpp(x)) {
    return mfp2_vector_to_matrix_cpp(x);
  }
  
  return transform_fp_core_internal(
    x,
    power,
    0.0,
    1.0,
    zero
  );
}


// Transform a vector by a single FP power.
//
// This helper is used by stored-ACD application:
//   1. transform x using the ACD model's stored power/shift/scale
//   2. transform acd(x) using the second selected ACD interaction power
//
// If power is NA, an n x 0 matrix is returned. This mirrors the R-side
// transform_vector_fp() behaviour where an NA power removes that component.
static NumericMatrix mfp2_transform_fp_one_cpp(const NumericVector& x,
                                               const double power,
                                               const double shift,
                                               const double scale,
                                               const bool zero,
                                               const bool check_binary) {
  const int n = x.size();
  
  if (NumericVector::is_na(power)) {
    return mfp2_empty_matrix_cpp(n);
  }
  
  if (check_binary && !zero && is_binary_vector_cpp(x)) {
    return mfp2_vector_to_matrix_cpp(x);
  }
  
  NumericVector power_one = NumericVector::create(power);
  
  return transform_fp_core_internal(
    x,
    power_one,
    shift,
    scale,
    zero
  );
}


// Apply stored ACD parameters to x.
//
// This implements the apply-only ACD path needed by build_adjustment_step().
// It does NOT fit or estimate ACD parameters.
//
// Equivalent R logic:
//
//   x_power <- transform_vector_fp(
//     x = x,
//     power = acd_parameter$power,
//     shift = acd_parameter$shift,
//     scale = acd_parameter$scale,
//     zero = zero,
//     check_binary = FALSE
//   )
//
//   pnorm(acd_parameter$beta0 + acd_parameter$beta1 * x_power[, 1L])
//
// Note:
//   transform_vector_fp() forces shift = 0 when zero = TRUE. We mirror that
//   behaviour here.
static NumericVector mfp2_apply_acd_stored_cpp(const NumericVector& x,
                                               const List& acd_parameter,
                                               const bool zero) {
  if (!acd_parameter.containsElementNamed("beta0") ||
      !acd_parameter.containsElementNamed("beta1") ||
      !acd_parameter.containsElementNamed("power") ||
      !acd_parameter.containsElementNamed("shift") ||
      !acd_parameter.containsElementNamed("scale")) {
      stop("Internal error: stored ACD parameters are incomplete.");
  }
  
  const double beta0 = as<double>(acd_parameter["beta0"]);
  const double beta1 = as<double>(acd_parameter["beta1"]);
  const double power = as<double>(acd_parameter["power"]);
  const double scale = as<double>(acd_parameter["scale"]);
  
  // Match transform_vector_fp(): when zero = TRUE, shift is ignored/forced to 0.
  const double shift = zero ? 0.0 : as<double>(acd_parameter["shift"]);
  
  NumericMatrix x_power = mfp2_transform_fp_one_cpp(
    x,
    power,
    shift,
    scale,
    zero,
    false
  );
  
  if (x_power.ncol() != 1) {
    stop("Internal error: stored ACD application must produce one FP column.");
  }
  
  const int n = x.size();
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    const double zhat = beta0 + beta1 * x_power(i, 0);
    out[i] = R::pnorm5(zhat, 0.0, 1.0, 1, 0);
  }
  
  return out;
}


// Transform one ACD adjustment variable using stored ACD parameters.
//
// This mirrors the apply-only branch of:
//
//   transform_vector_acd(..., acd_parameter = stored_parameter)
//
// for build_adjustment_step().
//
// For power = c(p_x, p_acd):
//   * p_x transforms the original x
//   * p_acd transforms acd(x)
//   * NA removes the corresponding component
//
// Returns:
//   cbind(FP(x, p_x), FP(acd(x), p_acd))
//
// or an n x 0 matrix if both powers are NA.
static NumericMatrix mfp2_transform_acd_adjustment_apply_cpp(
    const NumericVector& x,
    SEXP power_sexp,
    SEXP acd_parameter_sexp,
    const bool zero) {
  
  const int n = x.size();
  
  if (Rf_isNull(acd_parameter_sexp)) {
    stop(
      "Internal error: ACD adjustment variables require stored acd_parameter "
      "values in build_adjustment_step()."
    );
  }
  
  NumericVector power(power_sexp);
  
  if (power.size() != 2) {
    stop("Internal error: ACD adjustment powers must have length two.");
  }
  
  if (NumericVector::is_na(power[0]) && NumericVector::is_na(power[1])) {
    return mfp2_empty_matrix_cpp(n);
  }
  
  List acd_parameter(acd_parameter_sexp);
  
  NumericVector x_acd_raw = mfp2_apply_acd_stored_cpp(
    x,
    acd_parameter,
    zero
  );
  
  // First component: ordinary FP transform of x.
  NumericMatrix x_fp = mfp2_transform_fp_one_cpp(
    x,
    power[0],
         0.0,
         1.0,
         zero,
         true
  );
  
  // Second component: FP transform of acd(x).
  NumericMatrix x_acd = mfp2_transform_fp_one_cpp(
    x_acd_raw,
    power[1],
         0.0,
         1.0,
         false,
         true
  );
  
  return mfp2_cbind_two_cpp(x_fp, x_acd);
}

// Convert a catzero indicator to an n x 1 numeric matrix.
//
// In build_adjustment_step(), catzero[[varname]] is expected to be either:
//   * NULL, or
//   * an n x 1 matrix.
//
// The R object can be numeric, integer, or logical. The C++ adjustment loop
// works with NumericMatrix objects, so this helper validates the shape and
// coerces values to double.
static NumericMatrix mfp2_catzero_to_numeric_matrix_cpp(SEXP cz_sexp,
                                                        const int n,
                                                        const std::string& varname) {
  if (Rf_isNull(cz_sexp)) {
    return mfp2_empty_matrix_cpp(n);
  }
  
  if (!Rf_isMatrix(cz_sexp)) {
    stop(
      "Internal error: catzero for variable '" +
        varname +
        "' must be a matrix."
    );
  }
  
  IntegerVector dims = Rf_getAttrib(cz_sexp, R_DimSymbol);
  
  if (dims.size() != 2 || dims[0] != n || dims[1] != 1) {
    stop(
      "Internal error: catzero for variable '" +
        varname +
        "' must be an n x 1 matrix."
    );
  }
  
  NumericMatrix out(n, 1);
  
  switch (TYPEOF(cz_sexp)) {
  case REALSXP: {
    NumericMatrix cz(cz_sexp);
    
    for (int i = 0; i < n; ++i) {
      out(i, 0) = cz(i, 0);
    }
    
    break;
  }
  case INTSXP: {
    IntegerMatrix cz(cz_sexp);
    
    for (int i = 0; i < n; ++i) {
      const int value = cz(i, 0);
      out(i, 0) = IntegerVector::is_na(value) ?
      NA_REAL :
        static_cast<double>(value);
    }
    
    break;
  }
  case LGLSXP: {
    LogicalMatrix cz(cz_sexp);
    
    for (int i = 0; i < n; ++i) {
      const int value = cz(i, 0);
      out(i, 0) = LogicalVector::is_na(value) ?
      NA_REAL :
        static_cast<double>(value);
    }
    
    break;
  }
  default:
    stop(
      "Internal error: catzero for variable '" +
        varname +
        "' must be numeric, integer, or logical."
    );
  }
  
  return out;
}


// Assign adjustment-column names.
//
// This mirrors the R-side naming convention:
//
//   colnames(adj_mat) <- paste0(varname, "_adj", seq_len(ncol(adj_mat)))
//
// Variables that contribute zero columns are left unnamed.
static void mfp2_set_adjustment_colnames_cpp(NumericMatrix& mat,
                                             const std::string& varname) {
  const int p = mat.ncol();
  
  if (p == 0) {
    return;
  }
  
  CharacterVector names(p);
  
  for (int j = 0; j < p; ++j) {
    names[j] = varname + "_adj" + std::to_string(j + 1);
  }
  
  colnames(mat) = names;
}
// -----------------------------------------------------------------------------
// build_adjustment_step() exported C++ loop -----------------------------------
// -----------------------------------------------------------------------------
//
// This function is the C++ replacement for the hot R loop:
//
//   for (varname in vars_adj) {
//     ...
//     data_adj_list[[varname]] <- adj_mat
//   }
//
// plus:
//
//   data_adj <- do.call(cbind, data_adj_list)
//
// Scope:
//   * Used only by build_adjustment_step().
//   * Assumes R has already prepared aligned inputs.
//   * Applies stored ACD parameters; does not estimate/refit ACD.
//   * Handles cache reuse, eliminated variables, binary-only spike variables,
//     FP/ACD recomputation, catzero/spike assembly, column naming, and final
//     matrix binding.
//
// Input conventions:
//   * x_col_index is zero-based, because it is used directly in C++.
//   * vars_adj gives the output list names.
//   * powers_adj, acdx_adj, zero_adj, catzero_adj, spike_adj,
//     spike_decision_int_adj, acd_parameter_adj, eliminated,
//     spike_binary_only_flags, current_power_keys_adj, prev_power_keys_adj,
//     prev_data_adj_list, and prev_spike_decision_int_adj must all be aligned
//     to vars_adj.
//   * catzero_adj elements are either NULL or n x 1 matrices.
//   * prev_data_adj_list elements are either NULL or previously built matrices.
//
// Return:
//   list(
//     data_adj_list = named list of per-variable matrices,
//     data_adj      = final cbind-ed numeric matrix
//   )
//
// [[Rcpp::export]]
List build_adjustment_step_loop_cpp(const NumericMatrix& x,
                                    const IntegerVector& x_col_index,
                                    const CharacterVector& vars_adj,
                                    const List& powers_adj,
                                    const LogicalVector& acdx_adj,
                                    const LogicalVector& zero_adj,
                                    const List& catzero_adj,
                                    const LogicalVector& spike_adj,
                                    const IntegerVector& spike_decision_int_adj,
                                    const List& acd_parameter_adj,
                                    const LogicalVector& eliminated,
                                    const LogicalVector& spike_binary_only_flags,
                                    const List& current_power_keys_adj,
                                    const List& prev_power_keys_adj,
                                    const List& prev_data_adj_list,
                                    const IntegerVector& prev_spike_decision_int_adj,
                                    const bool has_prev) {
  
  const int n = x.nrow();
  const int p_x = x.ncol();
  const int n_vars = vars_adj.size();
  
  // --------------------------------------------------------------------------
  // Defensive alignment checks
  // --------------------------------------------------------------------------
  // These are internal checks. The R wrapper should already pass aligned inputs,
  // but failing here gives a clearer error than returning corrupted matrices.
  if (x_col_index.size() != n_vars ||
      powers_adj.size() != n_vars ||
      acdx_adj.size() != n_vars ||
      zero_adj.size() != n_vars ||
      catzero_adj.size() != n_vars ||
      spike_adj.size() != n_vars ||
      spike_decision_int_adj.size() != n_vars ||
      acd_parameter_adj.size() != n_vars ||
      eliminated.size() != n_vars ||
      spike_binary_only_flags.size() != n_vars ||
      current_power_keys_adj.size() != n_vars ||
      prev_power_keys_adj.size() != n_vars ||
      prev_data_adj_list.size() != n_vars ||
      prev_spike_decision_int_adj.size() != n_vars) {
    stop("Internal error: build_adjustment_step_loop_cpp() inputs are not aligned.");
  }
  
  // Per-variable matrices to return for cache reuse in later MFP iterations.
  List data_adj_list(n_vars);
  CharacterVector data_adj_names(n_vars);
  
  // Total number of columns across all per-variable adjustment matrices.
  // Needed to allocate final data_adj once.
  int total_cols = 0;
  
  // --------------------------------------------------------------------------
  // First pass: build or reuse each per-variable adjustment matrix
  // --------------------------------------------------------------------------
  for (int v = 0; v < n_vars; ++v) {
    const std::string varname = as<std::string>(vars_adj[v]);
    data_adj_names[v] = varname;
    
    const int col = x_col_index[v];
    
    if (col < 0 || col >= p_x) {
      stop(
        "Internal error: invalid column index for adjustment variable '" +
          varname +
          "'."
      );
    }
    
    const int spike_current = spike_decision_int_adj[v];
    
    // ------------------------------------------------------------------------
    // Previous cache entry
    // ------------------------------------------------------------------------
    // These correspond to the R-side:
    //
    //   prev_data_adj <- prev_xi$data_adj_list[[varname]]
    //   prev_spike_decision <- prev_spike_decision_int_adj[[varname]]
    //
    SEXP prev_data_sexp = R_NilValue;
    int prev_spike_decision = NA_INTEGER;
    
    if (has_prev) {
      prev_data_sexp = prev_data_adj_list[v];
      prev_spike_decision = prev_spike_decision_int_adj[v];
    }
    
    // ------------------------------------------------------------------------
    // Cache invalidation
    // ------------------------------------------------------------------------
    // Recompute only if either:
    //   1. normalized FP power key changed, or
    //   2. spike decision changed.
    //
    // This mirrors the current R logic using identical() for normalized power
    // keys and integer spike decisions.
    bool recompute = true;
    
    if (has_prev && !Rf_isNull(prev_data_sexp)) {
      const bool powers_same = mfp2_numeric_keys_identical_cpp(
        prev_power_keys_adj[v],
                           current_power_keys_adj[v]
      );
      
      const bool spike_decision_same =
        !IntegerVector::is_na(prev_spike_decision) &&
        prev_spike_decision == spike_current;
      
      recompute = !(powers_same && spike_decision_same);
    }
    
    NumericMatrix adj_mat;
    
    // ------------------------------------------------------------------------
    // Build or reuse this variable's adjustment matrix
    // ------------------------------------------------------------------------
    if (static_cast<bool>(spike_binary_only_flags[v])) {
      // Binary-only spike:
      // Use only the structural-zero indicator matrix. This branch is checked
      // before eliminated variables because powers are intentionally irrelevant
      // when spike_decision == 3.
      SEXP cz_sexp = catzero_adj[v];
      
      if (Rf_isNull(cz_sexp)) {
        stop(
          "Internal error: binary-only spike variable '" +
            varname +
            "' has no catzero indicator."
        );
      }
      
      adj_mat = mfp2_catzero_to_numeric_matrix_cpp(
        cz_sexp,
        n,
        varname
      );
      
    } else if (static_cast<bool>(eliminated[v])) {
      // Eliminated non-binary-only variable:
      // contributes no adjustment columns.
      adj_mat = mfp2_empty_matrix_cpp(n);
      
    } else if (!recompute) {
      // Cache hit:
      // reuse the previous transformed matrix for this variable.
      if (Rf_isNull(prev_data_sexp)) {
        stop(
          "Internal error: cache hit for '" +
            varname +
            "' has NULL cached adjustment data."
        );
      }
      
      adj_mat = NumericMatrix(prev_data_sexp);
      
      if (adj_mat.nrow() != n) {
        stop(
          "Internal error: cached adjustment matrix for '" +
            varname +
            "' has the wrong number of rows."
        );
      }
      
    } else {
      // Cache miss:
      // recompute the continuous/FP or stored-ACD part.
      NumericVector xvec(n);
      
      for (int i = 0; i < n; ++i) {
        xvec[i] = x(i, col);
      }
      
      if (static_cast<bool>(acdx_adj[v])) {
        // ACD adjustment variable:
        // build cbind(FP(x), FP(acd(x))) using stored ACD parameters.
        // This does not fit ACD parameters.
        adj_mat = mfp2_transform_acd_adjustment_apply_cpp(
          xvec,
          powers_adj[v],
                    acd_parameter_adj[v],
                                     static_cast<bool>(zero_adj[v])
        );
      } else {
        // Ordinary FP adjustment variable.
        adj_mat = mfp2_transform_fp_adjustment_cpp(
          xvec,
          powers_adj[v],
                    static_cast<bool>(zero_adj[v])
        );
      }
      
      // ----------------------------------------------------------------------
      // Add structural-zero indicator if needed
      // ----------------------------------------------------------------------
      // Matches the R-side catzero/spike rules:
      //
      //   spike_decision == 1: catzero + continuous/FP part
      //   spike_decision == 2: continuous/FP part only
      //   spike_decision == 3: catzero only
      //
      // Non-spike variables with catzero get catzero prepended.
      SEXP cz_sexp = catzero_adj[v];
      
      if (!Rf_isNull(cz_sexp)) {
        NumericMatrix cz_mat = mfp2_catzero_to_numeric_matrix_cpp(
          cz_sexp,
          n,
          varname
        );
        
        if (static_cast<bool>(spike_adj[v])) {
          if (spike_current == 1) {
            adj_mat = mfp2_cbind_two_cpp(cz_mat, adj_mat);
          } else if (spike_current == 2) {
            // Continuous/FP component only. No change.
          } else if (spike_current == 3) {
            adj_mat = cz_mat;
          }
        } else {
          adj_mat = mfp2_cbind_two_cpp(cz_mat, adj_mat);
        }
      }
    }
    
    // Assign informative column names:
    //   varname_adj1, varname_adj2, ...
    //
    // This preserves the current R naming convention before the final cbind.
    mfp2_set_adjustment_colnames_cpp(adj_mat, varname);
    
    data_adj_list[v] = adj_mat;
    total_cols += adj_mat.ncol();
  }
  
  data_adj_list.attr("names") = data_adj_names;
  
  // --------------------------------------------------------------------------
  // Second pass: cbind all per-variable matrices into final data_adj
  // --------------------------------------------------------------------------
  // This replaces:
  //
  //   data_adj <- do.call(cbind, data_adj_list)
  //
  // Empty n x 0 matrices contribute no columns.
  NumericMatrix data_adj(n, total_cols);
  CharacterVector data_adj_colnames(total_cols);
  
  int out_col = 0;
  
  for (int v = 0; v < n_vars; ++v) {
    NumericMatrix mat(data_adj_list[v]);
    CharacterVector cn = colnames(mat);
    
    for (int j = 0; j < mat.ncol(); ++j) {
      for (int i = 0; i < n; ++i) {
        data_adj(i, out_col) = mat(i, j);
      }
      
      if (cn.size() == mat.ncol()) {
        data_adj_colnames[out_col] = cn[j];
      } else {
        const std::string varname = as<std::string>(vars_adj[v]);
        data_adj_colnames[out_col] =
          varname + "_adj" + std::to_string(j + 1);
      }
      
      ++out_col;
    }
  }
  
  if (total_cols > 0) {
    colnames(data_adj) = data_adj_colnames;
  }
  
  return List::create(
    _["data_adj_list"] = data_adj_list,
    _["data_adj"]      = data_adj
  );
}
