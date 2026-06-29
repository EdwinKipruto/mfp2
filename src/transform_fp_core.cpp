#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Core fractional-polynomial transformation kernel.
// R-level preprocessing must handle:
//   * all-NA powers,
//   * removal/sorting of NA powers,
//   * binary-variable bypass,
//   * NULL shift/scale defaults,
//   * forcing shift = 0 when zero = TRUE,
//   * input validation,
//   * output column names.

// [[Rcpp::export]]
NumericMatrix transform_fp_core(NumericVector x_raw,
                                NumericVector power,
                                double shift_val,
                                double scale_val,
                                bool zero) {
  
  const int n = x_raw.size();
  const int m = power.size();
  
  if (m == 0) {
    return NumericMatrix(n, 0);
  }
  
  NumericMatrix out(n, m);
  
  const double inv_scale = 1.0 / scale_val;
  
  for (int i = 0; i < n; ++i) {
    
    if (zero && x_raw[i] <= 0.0) {
      continue;
    }
    
    double xi = (x_raw[i] + shift_val) * inv_scale;
    
    if (zero && xi <= 0.0) {
      continue;
    }
    
    const double log_xi = std::log(xi);
    
    out(i, 0) = (power[0] == 0.0) ? log_xi : std::pow(xi, power[0]);
    
    for (int j = 1; j < m; ++j) {
      if (power[j] == power[j - 1]) {
        out(i, j) = out(i, j - 1) * log_xi;
      } else {
        out(i, j) = (power[j] == 0.0) ? log_xi : std::pow(xi, power[j]);
      }
    }
  }
  
  return out;
}