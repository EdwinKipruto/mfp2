#' Function to estimate Approximate cumulative distribution (ACD)
#' 
#' @param x a numeric vector.
#' @param shift a numeric that is used to change the values of `x` to positive
#' values. The default value is `NULL`, and the program will estimate it
#' automatically.
#' @param s a vector of allowed FP powers. The default value is `NULL`,
#' meaning that the set {-2, -1, -0.5, 0, 0.5, 1, 2, 3} is used.
#' @param scale a numeric used to scale `x`. Default is `NULL` and the program
#' will estimate it automatically. 
acd <- function(x, power = NULL, shift = NULL, s = NULL, scale = NULL) {
  # # The user must supplier both power and shift or set both to NULL
  # if(is.null(power)){
  #   if(!is.null(shift)) stop("If shift is given then power must also be given")
  # }
  # if(is.null(shift)){
  #   if(!is.null(power)) stop("If power is given then shift must also be given")
  # }
  #
  # if(is.null(shift)){
  #   if(!is.null(power)) stop("If power is given then shift must also be given")
  # }
  # Shift x if the shifting factor is provided
  if (!is.null(shift)) {
    x <- x + shift
    # Estimate shift from the data when not provided
  } else {
    shift <- find_shift_factor(x)
    x <- x + shift
  }
  # check whether x > 0 after shifting.
  if (!all(x > 0)) stop("All x are not > 0 after shifting x by ", shift)
  # Scale x using user-supllied scaling factor. scale = 1 is noscaling
  if (!is.null(scale)) {
    x <- x / scale
    # Estimate scale from the data when not provided. Uses R&S 2008 formula. see the book
  } else {
    scale <- find_scale_factor(x)
    x <- x / scale
  }
  # Compute the rank of x. Ties are handled by finding the average value
  n <- length(x)
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.5) / n)
  # Estimate power from the data when not given
  if (is.null(power)) {
    # Estimate the best p in model E(z) = beta0 + beta1*x^p using the data
    fit <- fit_fp(
      y = z, x = x, family = "gaussian", method = NULL, shift = shift, scalex = NULL,
      weights = rep(1, n), offset = rep(0, n), strata = NULL, degree = 1,
      powers = s, control = NULL, rownames = NULL, nocenter = NULL
    )
    # The best FP Power estimated from the data
    power <- fit$bestpower
    # The regression coefficients of model E(z) = beta0 + beta1*x^p
    coefx <- fit$bestmodel$coefficients
    # fitted values of the model z_hat
    zhat <- fit$bestmodel$fitted.values
  } else {
    if (length(power) != 1) stop("Only one power term is supported")
    xnew <- ifelse(power == rep(0, n), log(x), x^power)
    # fit a linear model of z by x
    fit <- stats::lm.fit(x = cbind(rep(1, n), x), y = z)
    # coefficients
    coefx <- fit$coefficients
    # fitted values of the model z_hat
    zhat <- fit$fit$fitted.values
  }

  # Approximate cumulative distribution (ACD)
  acdx <- stats::pnorm(zhat)
  # Return acd, power, beta0, beta1, shift, scale
  return(list(acd = acdx, beta0 = coefx[1], beta1 = coefx[2], power = power, shift = shift, scale = scale))
}
