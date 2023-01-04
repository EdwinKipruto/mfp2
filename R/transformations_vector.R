#' Functions to transform a variable using fractional polynomial powers or acd
#' 
#' These functions generate fractional polynomials for a variable similar to
#' `fracgen` in Stata. `transform_acd_vector` generates the acd transformation
#' for a variable.
#' 
#' @details 
#' For example, this function  generates new variables for x given the power(s)
#' e.g., xnew = log(x) when power = 0 and scale = F and center = F.
#' 
#' @param x a vector of a predictor variable.
#' @param power The user-supplied FP power. Default is 1 (linear). Must be a 
#' vector of length 2 for acd transformation.
#' @param scale scaling factors for x of interest. Must be positive integers.
#' Default is `NULL` and scaling factors are automatically estimated by the
#' program. If no scaling is needed set `scale = 1`.
#' @param shift adjustment factors required for shifting x to positive
#' values. Default is `NULL` and adjustment factors are estimated automatically
#' using Royston and Sauerbrei formula iff any `x` <= 0. 
#' If no shifts are needed set `shift = 0`.
#' @param s passed to [acd()].
#' @param center Specification of centering for variable using
#' the mean i.e. `f(x) - mean(f(x))` for continuous variables and 
#' `x - min(x)` for binary variables. Default is no centering.
#' 
#' @return 
#' Returns a matrix of transformed variable(s). The number of columns
#' depends on the number of powers provided, the number of rows is equal to the
#' length of `x`.
#' 
#' @export
transform_fp_vector <- function(x, 
                                power = 1,
                                scale = NULL, 
                                center = FALSE, 
                                shift = NULL) {
  # number of observations
  N <- length(x)
  if (all(is.na(power))) { # important because omitted variables have power of NA
    x <- NULL
    xcenter <- NULL
  } else {
    # Do not transform x when it is a two level variable but center if necessary
    nx <- length(unique(x))
    if (nx <= 2) {
      # warning(paste0("x has 2 different values, original x ",if(center){paste0("(centered by min(x))")}, " values returned"))
      x <- as.matrix(x)
      # center using the lower(minimum) of the two distinct values of the covariates
      # as also done in stata mfp by Patrick Royston
      if (center) {
        xcenter <- min(x)
        x <- x - xcenter
      }
    } else {
      # shift x
      if (is.null(shift)) {
        shift <- find_shift_factor(x)
      }
      x <- x + shift
      # scale x
      if (is.null(scale)) {
        scale <- find_scale_factor(x)
      }
      x <- x / scale
      # sort the powers such that 1,0,1 is equal to 0,1,1. p1<p2<....<pm. see
      # Royston and Altman 1994 near equation 6 for explanation
      power <- sort(power) # if NA exist it will be removed
      # create a list that will hold the transformed variables
      H <- vector(mode = "list", length = length(power))
      # Deal with the first power
      H[[1]] <- ifelse(power[1] == rep(0, N), log(x), x^power[1])
      # check the length of the power. If length >=2 use a loop
      np <- length(power)
      # Deal with the remaining powers
      if (np > 1) {
        # check whether all powers are identical
        nu <- length(unique(power))
        # all powers are repeated
        if (nu == 1) {
          for (j in seq_len(np - 1)) {
            k <- j + 1
            H[[k]] <- H[[j]] * log(x)
          }
        } else {
          # Either a mix of powers like 1,2,2 or 1,2,3
          for (j in seq_len(np - 1)) { # check whether while loop is possible
            k <- j + 1
            # the subsequent power is repeated
            if (power[k] == power[j]) {
              H[[k]] <- H[[j]] * log(x)
              # the subsequent powers are not repeated e.g 1,2,4. so we are dealing with say 2,4
            } else {
              H[[k]] <- ifelse(power[k] == rep(0, N), log(x), x^power[k])
            }
          }
        }
      }
      # combine the transformed variables into a matrix
      x <- do.call(cbind, H)
      # center each transformed variable using: f(x) - mean(f(x))
      if (center) {
        x <- scale(x, scale = F)
      }
    }
  }
  
  x
}

#' @describeIn transform_fp_vector Function to generate acd transformation.
transform_acd_vector <- function(x, 
                                 power = c(1, 1),
                                 shift = NULL, 
                                 s = NULL, 
                                 scale = NULL, 
                                 center = FALSE) {
  # the length of powers must be equal to 2
  np <- length(power)
  if (np != 2) stop("The length of powers supplied is ", np, ". Two powers are needed")
  if (all(is.na(power))) {
    xx <- NULL
  } else {
    # Transform x and acdx using the supplied powers
    # Approximate cumulative distribution (ACD)
    xa <- acd(x, power = NULL, shift = shift, s = s, scale = scale)$acd
    xa <- transform_fp_vector(x = xa, power = power[2],
                              scale = scale, shift = shift, 
                              center = center)
    x2 <- transform_fp_vector(x = x, power = power[1],
                              scale = scale, shift = shift,
                              center = center)
    # return a matrix of two variables
    xx <- cbind(x2, xa)
  }
  
  xx
}
