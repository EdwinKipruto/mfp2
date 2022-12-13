#' Function that fits generalized linear models including gaussian. It uses
#' glm.fit() function and returns loglikelihood and regression coefficients.
#'
#' @param x A design matrix of dimension n * (p+1) where 1 is the intercept
#' @param y  A vector of response variable of length n.
#' @param family a family function e.g gaussian().  A character string naming
#' a family function is not allowed since we are using internal function, glm.fit
#' @param weights a numeric vector of 'prior weights' to be used in the fitting process.
#' @param offset a numeric vector of of a priori known component to be included
#' in the linear predictor during fitting. This should be of length equal to the number of cases
#' @return The log likelihood of the fitted model, regression coefficients and the number of parameters
#' @import stats
#' @export
fit_glm <- function(x, y, family, weights, offset) {
  # x should include 1s for the intercept, family should be a function, weights and
  # offset should be specified
  fit <- stats::glm.fit(x = x, y = y, family = family, weights = weights, offset = offset)
  # number of parameters- we can as well use fit$coefficients
  p <- fit$rank
  # allow for estimated scale parameter for Gaussian family
  fam <- fit$family$family
  if (fam == "gaussian") p <- p + 1
  # loglikelihood of the glm model fitted idea borrowed from logLik.glm in github
  loglk <- p - fit$aic / 2
  # add residual sum of squares for calculating F statistic in Gaussian models
  SSE <- sum(fit$residuals^2, na.rm = T)
  # we are returning loglikelihood so that we can calculate AIC and BIC of
  # an MFP model which will consider the powers as additional parameters

  # we should also return fit for use with other existing functions like summary etc.
  return(list(logl = loglk, coefficients = fit$coefficients, df = p, SSE = SSE, fit = fit))
}
