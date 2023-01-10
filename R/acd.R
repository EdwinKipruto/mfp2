#' Function to estimate Approximate cumulative distribution (ACD)
#' 
#' Fits acd transformation as outlined in Royston and Sauerbrei (2016).
#' 
#' @details 
#' Briefly, the estimation works as follows. First, the input data are shifted
#' to positive values and scaled as requested. Then 
#' \deqn{z = \Phi^{-1}(\frac{rank(x) - 0.5}{n}) }
#' is computed, where \eqn{n} is the number of elements in `x`, 
#' with ties in the ranks handled as averages. To approximate \eqn{z}, 
#' an FP1 model (least squares) is used, i.e. 
#' \eqn{E(z) = \beta_0 + \beta_1 (x)^p}, where \eqn{p} is chosen such that it 
#' provides the best fitting model among all possible FP1 models. 
#' The ACD transformation is then given as
#' \deqn{acd(x) = \Phi(\hat{z}),}
#' where the fitted values of the estimated model are used. 
#' 
#' @param x a numeric vector.
#' @param powers a vector of allowed FP powers. The default value is `NULL`,
#' meaning that the set {-2, -1, -0.5, 0, 0.5, 1, 2, 3} is used.
#' @param shift a numeric that is used to shift the values of `x` to positive
#' values. The default value is 0, meaning no shifting is conducted. 
#' If `NULL`, then the program will estimate an appropriate shift automatically
#' (see [find_shift_factor()]).
#' @param scale a numeric used to scale `x`. The default value is 1, meaning 
#' no scaling is conducted. If `NULL`, then the program will estimate 
#' an appropriate scaling factor automatically (see [find_scale_factor()]). 
#' 
#' @return 
#' A list is returned with components
#' 
#' * `acd`: the acd transformed input data.
#' * `beta0`: intercept of estimated model.
#' * `beta1`: coefficient of estimated model.
#' * `power`: estimated power.
#' * `shift`: shift value used for computations.
#' * `scale`: scaling factor used for computations. 
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2016. \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}
#' 
#' @export
fit_acd <- function(x, 
                    powers = NULL, 
                    shift = 0, 
                    scale = 1) {

  # preprocess data
  if (is.null(shift)) 
    shift <- find_shift_factor(x)
  
  if (is.null(scale)) 
    scale <- find_scale_factor(x)
  
  x <- (x + shift) / scale
  
  # check whether acd is estimable
  if (!all(x > 0)) 
    stop("! All x must be positive after shifting.", 
         "i Please use another shift value or let it be estimated by the programm by setting `shift = NULL`.")
  
  n <- length(x)
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.5) / n)

  # estimate the best p in model E(z) = beta0 + beta1*x^p using the data
  fit <- fit_fp(
    y = z, x = x, family = "gaussian", 
    method = NULL, shift = shift, scalex = NULL,
    weights = rep(1, n), offset = rep(0, n), 
    strata = NULL, degree = 1,
    powers = powers, 
    control = NULL, rownames = NULL, nocenter = NULL
  )
  
  coefx <- fit$bestmodel$coefficients
  zhat <- fit$bestmodel$fitted.values

  list(acd = stats::pnorm(zhat), 
       beta0 = coefx[1], 
       beta1 = coefx[2], 
       power = fit$bestpower, 
       shift = shift, 
       scale = scale)
}
