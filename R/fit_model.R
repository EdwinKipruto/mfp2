#' Function that fits models supported by `mfpa`
#' 
#' Fits generalized linear models and cox proportional hazard models. 
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
                      nocenter = NULL) {
  
  if (is.null(colnames(x)))
    colnames(x) <- colnames(x, do.NULL = FALSE)
  
  if (family == "cox") {
    # cox needs more work especially on how to handle strata
    fit <- fit_cox(
      x = x, y = y, strata = strata, weights = weights, offset = offset,
      control = control, method = method, rownames = rownames,
      nocenter = nocenter
    )
  } else {
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) family <- family()
    # add 1s for intercept
    xx <- cbind(rep(1, length(y)), x)
    # Note that x can be NULL when fitting a NULL model
    # (using only adjustment variable)
    # when all FP powers estimated are NA-all variables removed
    colnames(xx) <- if (is.null(ncol(x))) {
      "Intercept"
    } else {
      c("Intercept", colnames(x))
    }
    fit <- fit_glm(
      y = y, x = xx, family = family, weights = weights, offset = offset
    )
  }
  
  fit
}

#' Function that fits generalized linear models 
#' 
#' @details
#' Uses [stats::glm.fit()] to fit the model.
#'
#' @param x a matrix of predictors including intercept with nobs observations.
#' @param y a vector for the outcome variable.
#' @param family a family function e.g. `stats::gaussian()`.  
#' @param weights a numeric vector of length nobs of 'prior weights' to be used 
#' in the fitting process.
#' @param offset a numeric vector of length nobs of of a priori known component 
#' to be included in the linear predictor during fitting. 
#' 
#' @return 
#' A list with the following components: 
#' 
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `sse`: residual sum of squares.
#' * `fit`: the object returned by [stats::glm.fit()].
#' 
#' @import stats
fit_glm <- function(x,
                    y, 
                    family, 
                    weights, 
                    offset) {

  fit <- stats::glm.fit(
    x = x, y = y, family = family, weights = weights, offset = offset
  )

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
#' @details
#' Uses [survival::coxph.fit()] to fit the model.
#'
#' @param x a matrix of predictors excluding intercept with nobs observations.
#' @param y a `Surv` object.
#' @param weights a numeric vector of length nobs of 'prior weights' to be used 
#' in the fitting process.
#' @param offset a numeric vector of length nobs of of a priori known component 
#' to be included in the linear predictor during fitting. 
#' @param method a character string specifying the method for tie handling. 
#' See [survival::coxph()].
#' 
#' @return 
#' A list with the following components: 
#' 
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `sse`: residual sum of squares (not used).
#' * `fit`: the object returned by [survival::coxph.fit()].
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
                    nocenter) {
  
  # Set default for control
  if (is.null(control)) control <- survival::coxph.control()

  fit <- survival::coxph.fit(
    x = x, y = y, strata = strata, weights = weights, offset = offset,
    control = control, method = method, rownames = rownames, resid = TRUE,
    nocenter = nocenter
  )
  
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
