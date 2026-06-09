#' Calculate centering values using Royston's formula
#'
#' This function computes centering parameters using Royston's formula.
#' Centering values are obtained by applying the function to the mean of x,
#' resulting in `f(mean(x))`. The function is then centered as `f(x) - f(mean(x))`.
#' This function is not currently used in `mfpi()`. The `mfpi()` function
#' calculates its centering values using `mean(f(x))`.
#'
#' @param x a continuous variable in which the mean will be calculated.
#' @param power a vector of powers to be used to transform the mean of x.
#' @returns a matrix of centering values of dimension 1 x length(power)
#' @references
#' Royston, P., & Sauerbrei, W. (2008). Multivariable model-building: a pragmatic approach to
#' regression analysis based on fractional polynomials for modelling continuous variables.
#' John Wiley & Sons.
#' @keywords internal
royston_centering_parameters <- function(x, power = 1) {
  
  # Reference link: https://www.stata.com/manuals/rfp.pdf
  # Calculate the mean of x
  meanx <- mean(x)
  
  # Calculate the number of powers
  num_powers <- length(power)
  
  # Initialize matrix to store centering parameters
  center_parameters <- matrix(NA, nrow = 1, ncol = num_powers)
  
  # Calculate centering parameters for each power
  center_parameters[, 1] <- transform_vector_fp(x = meanx, power = power[1], check_binary = FALSE)
  if (num_powers != 1) {
    for (j in seq_len(num_powers - 1)) {
      k <- j + 1
      if (power[k] == power[j]) {
        center_parameters[, k] <- transform_vector_fp(x = meanx, power = power[k], check_binary = FALSE) * log(meanx)
      } else {
        center_parameters[, k] <- transform_vector_fp(x = meanx, power = power[k], check_binary = FALSE)
      }
    }
  }
  
  return(center_parameters)
}

