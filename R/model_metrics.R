#' Deviance computations as used in mfp in stata
#' 
#' Note that this is not the usual formula of deviance used in R, but
#' uses the formula found here https://www.stata.com/manuals/rfp.pdf.
deviance_stata <- function(rss, weights, n) {
  
  # calculate lognormalized weights
  if (length(unique(weights)) == 1) {
    meanwts <- 0
  } else {
    logwts <- log(weights)
    mlogwts <- mean(logwts)
    # lognormalize weights
    lognormweights <- logwts / mlogwts
    # average of lognormalized weights
    meanwts <- mean(lognormweights)
  }
  k <- log((2 * pi * rss) / n)
  
  n * (1 - meanwts + k)
}

#' Function to compute model metrics to be used within `mfpa`
#' 
#' Mostly used within an mfp step to compare between the different fp models
#' of a variable. 
#' 
#' @param obj a list returned by [fit_model()] representing a glm or Cox model
#' fit.
#' @param n_obs a numeric value indicating the number of observations for the
#' data used to fit `obj`.
#' 
#' @return 
#' A list with the following entries:
#' 
#' * `df`: number of degrees of freedom of model (i.e. coefficients).
#' * `deviance_rs`: "deviance", i.e. minus twice the log likelihood. This is not
#' usual definition of deviance used by R, which is defined as twice the 
#' difference between the log likelihoods of the saturated model (one parameter 
#' per observation) and the null (or reduced) model. It is, however, the 
#' definition used in Royston and Sauerbrei (2008) and in `mfp`. For selection
#' of fps this does not really play a role, as the common factor would be 
#' cancelled anyway when comparing models based on deviances. 
#' * `sse`: sum of squared residuals as returned by [fit_model()].
#' * `deviance_stata`: deviance computed by [deviance_stata()].
#' * `aic`: Akaike information criterion, defined as `-2logL + 2df`.
#' * `bic`: Bayesian information criterion, defined as `-2logL + log(n_obs)df`.
#' * `df_resid`: residual degrees of freedom. For consistency with stata we 
#' subtract the scale parameter. 
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
calculate_model_metrics <- function(obj, n_obs) {

  res <- list(
    df = obj$df, 
    deviance_rs = -2 * obj$logl, 
    sse = obj$sse
  )
  
  res$deviance_stata <- deviance_stata(
    rss = res$sse, weights = obj$weights, n = n_obs
  )
  res$aic <- res$deviance_rs + 2 * res$df
  res$bic <- res$deviance_rs + log(n_obs) * res$df
  res$df_resid <- n_obs - (res$df - 1)
  
  res
}