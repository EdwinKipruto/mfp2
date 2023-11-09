#' Function to estimate approximate cumulative distribution (ACD)
#' 
#' Fits ACD transformation as outlined in Royston (2014). The ACD transformation 
#' smoothly maps the observed distribution of a continuous covariate x onto one scale,
#' namely, that of an approximate uniform distribution on the interval (0, 1). 
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
#' If the relationship between a response Y and acd(x) is linear,
#' say, \eqn{E(Y) = \beta_0 + \beta_1 acd(x)}, the relationship between Y 
#' and x is nonlinear and is typically sigmoid in shape.  
#' The parameters \eqn{\beta_0} and \eqn{\beta_0 + \beta_1} in such a model are 
#' interpreted as the expected values of Y at the minimum and maximum of x,
#' that is, at acd(x) = 0 and 1, respectively. 
#' The parameter \eqn{\beta_1} represents the range of predictions of \eqn{E(Y)} 
#' across the whole observed distribution of x (Royston 2014).
#' 
#' @param x a numeric vector.
#' @param powers a vector of allowed FP powers. The default value is `NULL`,
#' meaning that the set \eqn{S = (-2, -1, -0.5, 0, 0.5, 1, 2, 3)} is used.
#' @param shift a numeric that is used to shift the values of `x` to positive
#' values. The default value is 0, meaning no shifting is conducted. 
#' If `NULL`, then the program will estimate an appropriate shift automatically
#' (see [find_shift_factor()]).
#' @param scale a numeric used to scale `x`. The default value is 1, meaning 
#' no scaling is conducted. If `NULL`, then the program will estimate 
#' an appropriate scaling factor automatically (see [find_scale_factor()]). 
#' @examples
#'
#'set.seed(42)
#' x = apply_shift_scale(rnorm(100))
#' y = rnorm(100)
#' fit_acd(x, y)
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
#' Royston, P. and Sauerbrei, W. (2016). \emph{mfpa: Extension of mfp using the
#' ACD covariate transformation for enhanced parametric multivariable modeling. 
#' The Stata Journal, 16(1), pp.72-87.}
#' 
#' Royston, P. (2014). \emph{A smooth covariate rank transformation for use in 
#' regression models with a sigmoid dose–response function. The Stata Journal,
#'  14(2), 329-341}. 
#' 
#' @export
fit_acd <- function(x, 
                    powers = NULL, 
                    shift = 0, 
                    scale = 1) {
  
  if (is.null(powers)) {
    # default FP powers proposed by Royston and Sauerbrei (2008)
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }

  # preprocess data
  if (is.null(shift)) {
    shift <- find_shift_factor(x)
  } 
  
  if (is.null(scale)) {
    scale <- find_scale_factor(x)
  }
  
  x <- (x + shift) / scale
  
  # check whether acd is estimable
  if (!all(x > 0)) 
    stop("! All x must be positive after shifting.", 
         "i Please use another shift value or let it be estimated by the program by setting `shift = NULL`.")
  
  # see here for details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5663339/pdf/emss-59479.pdf)
  n <- length(x)
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.5) / n)

  # estimate the best p in model E(z) = beta0 + beta1*x^p using the data
  fit <- find_best_fp1_for_acd(x = x, y = z, powers = powers)
  
  coefx <- fit$fit$coefficients
  zhat <- fit$fit$fitted.values

  list(acd = stats::pnorm(zhat), 
       beta0 = coefx[1], 
       beta1 = coefx[2], 
       power = fit$power, 
       shift = shift, 
       scale = scale)
}

#' Function to apply Approximate Cumulative Distribution (ACD)
#' 
#' Applies the acd transformation as outlined in Royston (2014) and Royston and 
#' Sauerbrei (2016). 
#' Designed to work with the output of [fit_acd()], Please refer to the corresponding
#' documentation for more details.
#' 
#' @param x a numeric vector.
#' @param beta0,beta1 each a numeric value, representing the coefficients of 
#' the FP1 model for the ACD transformation.
#' @param power a numeric value, estimated power to be used in the FP1 model for 
#' the ACD transformation.
#' @param shift a numeric value that is used to shift the values of `x` to
#' positive values. 
#' @param scale a numeric value used to scale `x`. 
#' @param ... not used.
#' 
#' @return 
#' The transformed input vector `x`.
apply_acd <- function(x, beta0, beta1, power, shift, scale, ...) {
  
  if (length(power) != 1) 
    stop("! `power` must be a single numeric value.")
  
  x <- (x + shift) / scale
  zhat <- beta0 + beta1 * transform_vector_power(x, power)
  
  stats::pnorm(zhat)
}

#' Function to fit univariable FP1 models for acd transformation
#' 
#' To be used in [fit_acd()].
#' 
#' @inheritParams fit_acd
#' @param y normal cdf of rank transform of `x`.
#' 
#' @return 
#' The best FP power with smallest deviance and the fitted model.
find_best_fp1_for_acd <- function(x, 
                                  y, 
                                  powers) {
  
  if (!is.null(dim(x))) 
    stop("! `x` must be a vector.")
  
  # generate all possible FP1 transformations
  trafo <- generate_transformations_fp(x = x, degree = 1, powers = powers)$data
  
  # Fit linear models for each FP1 function
  n_powers <- length(powers)
  
  # store deviance and model object
  devs <- vector("numeric", n_powers)
  fits <- vector("list", n_powers)
  
  for (i in seq_len(n_powers)) {
    # fit linear model
    fit <- fit_model(x = trafo[[i]], y = y, family = "gaussian")
    
    devs[i] <- -2 * fit$logl
    fits[[i]] <- fit$fit
  }
  
  # find best model
  index <- which.min(devs)
  
  list(
    power = powers[[index]], 
    fit = fits[[index]]
  )
}
