
#' A function that is used to shift x values to positive values if it contains
#' negative or zero values.If all values of x are positive then the original
#' values of x is returned without shifting.
#' if x has already been shifted and scaled then the function does nothing

#' @param x A vector of predictor variable
#' @param scale scaling factors for x of interest. Must be positive integers.
#' Default is NULL and  scaling factors are automatically estimated using
#' scalefn() function else it uses user suppled scaling factors. If no scaling
#' is needed just use scale = 1 because x/1 = x
#' @param shift adjustment factors required for shifting x to positive
#' values. Default is NULL and adjustment factors are estimated automatically
#' using Royston and Sauerbrei formula.
#' @export
# ===============================================================================
shift.scale <- function(x, scale, shift) {
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
      shift <- shift.factor(x)
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
    x <- x / scalefn(x)
  } else {
    x <- x / scale
  }
  return(x)
}
