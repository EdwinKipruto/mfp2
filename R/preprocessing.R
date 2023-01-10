#' A function that calculates an integer used to scale x. Estimation of powers
#' is unaffected by scaling, the same powers are found for x/scale. In extreme
#' cases scaling is necessary to preserve accuracy (Royston and Sauerbrei 2008).
#' We used the formula in section 4.11.1 of Royston and Sauerbrei (2008) book
#'
#' @param x a vector of predictor variable where x is already shifted to positive values.
#' @return An integer
#' @export
find_scale_factor <- function(x) {
  # see https://www.stata.com/manuals/rfp.pdf
  mx <- unique(x)
  if (length(mx) == 1) stop("all values of x are identical hence log(max(x)-min(x)) = log(0) not allowed")
  lrange <- log10(max(x, na.rm = T) - min(x, na.rm = T))
  # scales <- 10^(sign(lrange)*as.integer(lrange))
  scales <- 10^(sign(lrange) * floor(abs(lrange)))
  return(scales)
}

# estimates factors for shifting continuous variables to positive values. If x
# is a two level variable or all values are positive then a factor of 0 is return
# hence no shifting will be conducted
find_shift_factor <- function(x) {
  nxx <- length(unique(x))
  if (all(x > 0) || nxx <= 2) {
    adj.factors <- 0
  } else {
    # see section 4.7 in Royston and Sauerbrei (2008) book for the formula
    difx <- diff(sort(x))
    # get rid of 0 differences and find the minimum
    eps <- min(difx[difx != 0])
    # adjustment factor
    adj.factors <- eps - min(x, na.rm = T)
  }
  return(adj.factors)
}

#' A function that is used to shift x values to positive values if it contains
#' negative or zero values.If all values of x are positive then the original
#' values of x is returned without shifting.
#' if x has already been shifted and scaled then the function does nothing
#' 
#' @param x A vector of predictor variable
#' @param scale scaling factors for x of interest. Must be positive integers.
#' Default is NULL and  scaling factors are automatically estimated using
#' find_scale_factor() function else it uses user suppled scaling factors. If no scaling
#' is needed just use scale = 1 because x/1 = x
#' @param shift adjustment factors required for shifting x to positive
#' values. Default is NULL and adjustment factors are estimated automatically
#' using Royston and Sauerbrei formula.
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
