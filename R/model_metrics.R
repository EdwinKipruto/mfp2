#' Deviance computations as used in mfp in stata
#' 
#' @details 
#' Note that this is not the usual formula of deviance used in R, but
#' uses the formula found here https://www.stata.com/manuals/rfp.pdf.
#' 
#' It can be applied for normal error models, but should not be used for other
#' kinds of glms.
deviance_gaussian <- function(rss, weights, n) {
  
  if (any(is.null(rss), is.null(weights), is.null(n))) return(NULL)
  
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

#' Function to compute model metrics to be used within `mfp2`
#' 
#' Mostly used within an mfp step to compare between the different fp models
#' of a variable. 
#' 
#' @param obj a list returned by [fit_model()] representing a glm or Cox model
#' fit.
#' @param n_obs a numeric value indicating the number of observations for the
#' data used to fit `obj`.
#' @param df_additional a numeric value indicating the number of additional
#' degrees of freedom to be accounted for in the computations of AIC and BIC. 
#' These may be necessary when a model uses FP terms, as these add another
#' degree of freedom per estimated power. 
#' 
#' @return 
#' A numeric vector with the following entries:
#' 
#' * `df`: number of degrees of freedom of model (i.e. coefficients plus 
#' `df_additional`).
#' * `deviance_rs`: "deviance", i.e. minus twice the log likelihood. 
#' This is not the usual definition of deviance used by R, which is defined as 
#' twice the difference between the log likelihoods of the saturated model
#' (one parameter per observation) and the null (or reduced) model. 
#' It is, however, the definition used in Royston and Sauerbrei (2008) and in 
#' `mfp`. For selection of fps this does not really play a role, as the common 
#' factor would be cancelled anyway when comparing models based on deviances. 
#' * `sse`: sum of squared residuals as returned by [fit_model()].
#' * `deviance_gaussian`: deviance computed by [deviance_gaussian()], 
#' applicable to Gaussian models and used for F-test computations.
#' * `aic`: Akaike information criterion, defined as
#' `-2logL + 2(df + df_additional)`.
#' * `bic`: Bayesian information criterion, defined as 
#' `-2logL + log(n_obs)(df + df_additional)`.
#' * `df_resid`: residual degrees of freedom, defined as `n_obs - df`. 
#' For consistency with stata we subtract the scale parameter from `df`. 
#' 
#' @references 
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
calculate_model_metrics <- function(obj, 
                                    n_obs, 
                                    df_additional = 0) {

  res <- c(
    logl = obj$logl,
    df = obj$df + df_additional, 
    deviance_rs = -2 * obj$logl, 
    sse = obj$sse
  )
  
  c(res, 
    deviance_gaussian = deviance_gaussian(
      rss = res[["sse"]], weights = obj$fit$weights, n = n_obs
    ), 
    aic = res[["deviance_rs"]] + 2 * res[["df"]],
    bic = res[["deviance_rs"]] + log(n_obs) * res[["df"]], 
    df_resid = n_obs - (res[["df"]] - 1)
  )
}