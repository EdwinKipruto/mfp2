#' Function that fits models supported by `mfpa`
#' 
#' Fits generalized linear models and Cox proportional hazard models. 
#' 
#' @details 
#' Computations rely on [fit_glm()] and [fit_cox()].
#'
#' @param x a matrix of predictors (excluding intercept) with column names.
#' If column names are not provided they are set according to
#' `colnames(x, do.NULL = FALSE)`.
#' @param y a vector for the outcome variable for glms, and a `Surv` object 
#' for Cox models.
#' @param method a character string specifying the method for tie handling. 
#' See [survival::coxph()].
#' @param family a character strong specifying glm family to be used, or "cox"
#' for Cox models.
#' @param strata,control,weights,offset,rownames,nocenter parameters for Cox 
#' model. See [survival::coxph()] for details.
#' @param fast passed to [fit_glm()] and [fit_cox()].
#' 
#'  @return 
#' A list with the following components: 
#' 
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `sse`: residual sum of squares.
#' * `fit`: the object returned by the fitting procedure.
#' 
#' @importFrom stats family
fit_model <- function(x,
                      y, 
                      family = "gaussian", 
                      weights = NULL,
                      offset = NULL, 
                      method = NULL, 
                      strata = NULL, 
                      control = NULL,
                      rownames = NULL,
                      nocenter = NULL, 
                      fast = TRUE) {
  
  if (is.null(colnames(x)))
    colnames(x) <- colnames(x, do.NULL = FALSE)
  
  if (family == "cox") {
    # cox needs more work especially on how to handle strata
    fit <- fit_cox(
      x = x, y = y, strata = strata, weights = weights, offset = offset,
      control = control, method = method, rownames = rownames,
      nocenter = nocenter, fast = fast
    )
  } else {
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) family <- family()
  
    fit <- fit_glm(
      y = y, x = x, family = family, weights = weights, offset = offset, 
      fast = fast
    )
  }
  
  fit
}

#' Function that fits generalized linear models 
#'
#' @param x a matrix of predictors including intercept with nobs observations.
#' @param y a vector for the outcome variable.
#' @param family a family function e.g. `stats::gaussian()`.  
#' @param weights a numeric vector of length nobs of 'prior weights' to be used 
#' in the fitting process.
#' @param offset a numeric vector of length nobs of of a priori known component 
#' to be included in the linear predictor during fitting. 
#' @param fast a logical which determines how the model is fitted. The default
#' `TRUE` uses fast fitting routines (i.e. [stats::glm.fit()]), while `FALSE`
#' uses the normal fitting routines (used for the final output of `mfpa`).
#' 
#' @return 
#' A list with the following components: 
#' 
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `sse`: residual sum of squares.
#' * `fit`: the fitted model object.
#' 
#' @import stats
fit_glm <- function(x,
                    y, 
                    family, 
                    weights, 
                    offset, 
                    fast = TRUE) {

  if (fast) {
    # add 1s for intercept
    xx <- cbind(rep(1, length(y)), x)
    # Note that x can be NULL when fitting a NULL model
    # (using only adjustment variable)
    # when all FP powers estimated are NA-all variables removed
    colnames(xx) <- if (is.null(ncol(x))) {
      "(Intercept)"
    } else {
      c("(Intercept)", colnames(x))
    }
    
    fit <- stats::glm.fit(
      x = xx, y = y, family = family, weights = weights, offset = offset
    )  
  } else {
    fit <- glm(y ~ ., data = data.frame(x, y), 
               family = family, weights = weights, offset = offset)
  }

  # account for estimation of variance parameter in gaussian models
  # computation as in logLik.glm using rank
  df <- if (fit$family$family == "gaussian") fit$rank + 1 else fit$rank
  
  list(
    fit = fit,
    # loglikelihood computed as in stats::logLik.glm
    logl = df - fit$aic / 2,
    coefficients = fit$coefficients,
    df = df,
    sse = sum(fit$residuals^2, na.rm = TRUE)
  )
}

#' Function that fits Cox proportional hazards models
#'
#' @param x a matrix of predictors excluding intercept with nobs observations.
#' @param y a `Surv` object.
#' @param weights a numeric vector of length nobs of 'prior weights' to be used 
#' in the fitting process.
#' @param offset a numeric vector of length nobs of of a priori known component 
#' to be included in the linear predictor during fitting. 
#' @param method a character string specifying the method for tie handling. 
#' See [survival::coxph()].
#' @param fast a logical which determines how the model is fitted. The default
#' `TRUE` uses fast fitting routines (i.e. [survival::coxph.fit()]), while
#' `FALSE`uses the normal fitting routines (used for the final output of 
#' `mfpa`).
#' 
#' @return 
#' A list with the following components: 
#' 
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `sse`: residual sum of squares (not used).
#' * `fit`: the fitted model object.
#' 
#' @import survival
fit_cox <- function(x, 
                    y, 
                    strata, 
                    weights, 
                    offset, 
                    control, 
                    method, 
                    rownames, 
                    nocenter, 
                    fast = TRUE) {
  
  # Set default for control
  if (is.null(control)) control <- survival::coxph.control()

  if (fast) {
    fit <- survival::coxph.fit(
      x = x, y = y, strata = strata, weights = weights, offset = offset,
      control = control, method = method, rownames = rownames, resid = TRUE,
      nocenter = nocenter
    )  
  } else {
    fit <- survival::coxph(
      x = x, y = y, strata = strata, weights = weights, offset = offset,
      control = control, method = method, rownames = rownames, resid = TRUE,
      nocenter = nocenter
    )
  }
  
  list(
    fit = fit, 
    # coxph.fit() returns loglikelihood for a null and a model with predictors.
    # If x is a null matrix, one loglikelihood for the null model is returned.
    logl = ifelse(!is.null(ncol(x)), fit$loglik[2], fit$loglik), 
    coefficients = fit$coefficients, 
    # sometimes coefficients can be NA
    # for example when including same variables in the model
    df = length(fit$coefficients[!is.na(fit$coefficients)]), 
    sse = sum(fit$residuals^2)
  )
}
