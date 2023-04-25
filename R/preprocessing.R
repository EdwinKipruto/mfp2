#' Function that calculates an integer used to scale predictor
#' 
#' @details
#' For details on why scaling is useful see the corresponding section in the
#' documentation of [mfpa()].
#' 
#' The determination of the scaling factor is independent (i.e. not affected 
#' by) shifts in the input data, as it only depends on the range of the 
#' input data.
#' 
#' Note that the estimation of powers is unaffected by scaling, the same powers 
#' are found for scaled input data. In extreme cases scaling is necessary to 
#' preserve accuracy, see Royston and Sauerbrei (2008).
#' This formula uses the scaling formula from Section 4.11.1 of 
#' Royston and Sauerbrei (2008). Further information can also be found in the 
#' Stata manual for mfp at https://www.stata.com/manuals/rfp.pdf.
#'
#' @param x a numeric vector already shifted to positive values (see 
#' [find_shift_factor()]). This function requires at least 2 distinct values to 
#' work.
#' 
#' @return 
#' An integer that can be used to scale `x` to a reasonable range. For binary
#' variables 1 is returned.
#' 
#' @references 
#' Royston, P., and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @export
find_scale_factor <- function(x) {

  n_unique = length(unique(x))
  
  if (n_unique == 1)
    stop("! Input data must not be constant.", 
         "i All values of x are identical, hence log(max(x)-min(x)) = log(0) is not defined.")
  
  if (n_unique == 2)
    return(1)
    
  lrange <- log10(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  
  10^(sign(lrange) * floor(abs(lrange)))
}

#' Function that calculates a value used to shift predictor
#' 
#' @details
#' For details on why shifting is necessary see the corresponding section in the
#' documentation of [mfpa()].
#' 
#' This function implements the formula in Section 4.7 of Royston and 
#' Sauerbrei (2008).
#' 
#' @param x a numeric vector.
#' 
#' @return 
#' A numeric value that can be used to shift `x` to positive values. 
#' If all values are positive, or if `x` is binary then 0 is returned. 
#' 
#' @references 
#' Royston, P., and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @export
find_shift_factor <- function(x) {
  
  n_unique <- length(unique(x))
  
  if (all(x > 0) || n_unique <= 2) 
    return(0)
  
  difx <- diff(sort(x))
  eps <- min(difx[difx != 0])
  
  eps - min(x, na.rm = TRUE)
}

#' A function that is used to shift x values to positive values if it contains
#' negative or zero values.If all values of x are positive then the original
#' values of x is returned without shifting .
#' if x has already been shifted and scaled then the function does nothing
#' 
#' @param x A vector of predictor variable
#' @param scale scaling factors for x of interest. Must be positive integers.
#' Default is NULL and  scaling factors are automatically estimated using
#' find_scale_factor() function else it uses user supplied scaling factors. If no scaling
#' is needed just use scale = 1
#' @param shift adjustment factors required for shifting x to positive
#' values. Default is NULL and adjustment factors are estimated automatically
#' using find_shift_factor() function
#' 
#' @export
apply_shift_scale <- function(x, scale, shift) {
  # restrict x to be a vector not matrix
  if (is.matrix(x)) stop("x must be a vector not a matrix")
  N <- length(x)
  # If adjustment factors are NULL we use R&S formula to shift x to positive values
  if (is.null(shift)) {
    # Check whether all x values are positive. If true-No need of shifting
    if (all(x > 0)) {
      x <- x
    } else {
      # estimate adjustment factors
      shift <- find_shift_factor(x)
      # Shift x to positive
      x <- x + shift
    }
    # use the adjustment factors supplied by the user to shift x to positive values
  } else {
    # term + a
    x <- x + shift
    # check whether all x are now positive
    if (!all(x > 0)) stop("The minimum value of x after shifting x is ", min(x, na.rm = T), " which is not > 0. Check your adjustment factors")
  }
  # if scale is NULL then scale x for computational stability using R&S formula
  if (is.null(scale)) { # No scaling
    x <- x / find_scale_factor(x)
  } else {
    x <- x / scale
  }
  return(x)
}
