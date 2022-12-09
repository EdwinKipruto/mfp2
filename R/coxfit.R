
# Fits Cox proportional hazards regression model using coxph() function from
# survival package(). For now it fits cox model without strata. We will include
# it in later stages
#
# param x A design matrix of dimension n * p where n is the number of
# observations and p the number of predictors
# param y is a Surv object from the survival package generated using Surv()
# function
# param method A character string specifying the method for tie handling
# See coxph() for explanation of other parameters
# return The loglikelihood of the fitted model and regression coefficients
# import survival
# export
#' @import survival
coxfit <- function(x, y, strata, weights, offset, control, method, rownames, nocenter) {
  # Set default for control
  if (is.null(control)) control <- survival::coxph.control()
  # Fit cox model
  fit <- survival::coxph.fit(
    x = x, y = y, strata = strata, weights = weights, offset = offset,
    control = control, method = method, rownames = rownames, resid = T,
    nocenter = nocenter
  )
  # coxph.fit() returns loglikelihood for a null and a model with predictors.
  # If x is a null matrix, one loglikelihood for the null model is returned.
  loglik1 <- ifelse(!is.null(ncol(x)), fit$loglik[2], fit$loglik)
  # Degrees of freedom. sometimes coefficients can be NA like for example including
  # same variables in the model
  df <- length(fit$coefficients[!is.na(fit$coefficients)])
  # we will just output SSE for uniformity with glm() otherwise it's not used anywhere
  return(list(fit = fit, logl = loglik1, coefficients = fit$coefficients, df = df, SSE = sum(fit$residuals^2)))
}
