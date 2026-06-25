#' Function that fits models supported by `mfp2`
#' 
#' Fits generalized linear models and Cox proportional hazard models. 
#' 
#' @details 
#' Computations rely on \code{fit_glm()} and \code{fit_cox()}.
#'
#' @param x a matrix of predictors (excluding intercept) with column names.
#' If column names are not provided they are set according to
#' `colnames(x, do.NULL = FALSE)`.
#' @param y Response variable. For GLMs, this may be a numeric vector, a factor
#' response accepted by [stats::glm()], or for binomial models a two-column
#' matrix of grouped counts `cbind(successes, failures)`. For Cox models, this
#' must be a [survival::Surv()] object.
#' @param method a character string specifying the method for tie handling. 
#' See [survival::coxph()].
#' @param family a character strong specifying glm family to be used, or "cox"
#' for Cox models. The default family is set to 'Gaussian'.
#' @param strata,control,weights,offset,rownames,nocenter parameters for Cox 
#' or glm. See [survival::coxph()] or [stats::glm()] for details.
#' @param fast passed to \code{fit_glm()} and \code{fit_cox()}.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
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
#' @keywords internal
#' @noRd
fit_model <- function(x,
                      y, 
                      family, 
                      weights = NULL,
                      offset = NULL, 
                      method = NULL, 
                      strata = NULL, 
                      control = NULL,
                      rownames = NULL,
                      nocenter = NULL, 
                      fast = TRUE,
                      has_offset = FALSE) {
  
  # Set column names if not provided
  if (!is.null(dim(x)) && is.null(colnames(x))) {
    colnames(x) <- colnames(x, do.NULL = FALSE)
  }
  
  # Extract family string and convert to family object if needed
  # Extract family string. `family` is expected to have been validated by
  # mfp2.default() before reaching this internal fitting helper.
  if (is.character(family)) {
    family_string <- family
    
    if (family != "cox") {
      family <- switch(
        family,
        gaussian = stats::gaussian(),
        binomial = stats::binomial(),
        poisson  = stats::poisson()
      )
    }
  } else if (is.function(family)) {
    family <- family()
    family_string <- family$family
  } else {
    family_string <- family$family
  }
  
  if (family_string == "cox") {
    # cox needs more work especially on how to handle strata
    fit <- fit_cox(
      x = x, 
      y = y, 
      strata = strata, 
      weights = weights, 
      offset = offset,
      control = control, 
      method = method, 
      rownames = rownames,
      nocenter = nocenter, 
      fast = fast,
      has_offset = has_offset
    )
  } else {
    fit <- fit_glm(
      y = y,
      x = x,
      family = family, 
      weights = weights, 
      offset = offset, 
      fast = fast,
      has_offset = has_offset
    )
  }
  
  fit
}

#' Function that fits generalized linear models 
#'
#' @param x a matrix of predictors with nobs observations.
#' @param y Response variable. For GLMs, this may be a numeric vector, a factor
#' response accepted by [stats::glm()], or for binomial models a two-column
#' matrix of grouped counts `cbind(successes, failures)`. 
#' @param family a family function e.g. `stats::gaussian()`.  
#' @param weights a numeric vector of length nobs of 'prior weights' to be used 
#' in the fitting process. see [stats::glm()] for details.
#' @param offset a numeric vector of length nobs of of a priori known component 
#' to be included in the linear predictor during fitting. 
#' @param fast a logical which determines how the model is fitted. The default
#' `TRUE` uses fast fitting routines (i.e. [stats::glm.fit()]), while `FALSE`
#' uses the normal fitting routines ([stats::glm()]) (used for the final output
#' of `mfp2`). 
#' The difference is mainly due to the fact that normal fitting routines have
#' to handle data.frames, which is a lot slower than using the model matrix
#' and outcome vectors directly. 
#' @param has_offset logical indicating whether `offset` should be included in
#' the final formula-based fit when `fast = FALSE`.
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
#' @keywords internal
#' @noRd
fit_glm <- function(x,
                    y, 
                    family, 
                    weights, 
                    offset, 
                    fast = TRUE,
                    has_offset = FALSE) {
  
  nobs <- NROW(y)
  
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("! `offset` must be numeric.", call. = FALSE)
    }
    
    if (length(offset) != nobs) {
      stop("! `offset` must have one value per observation.", call. = FALSE)
    }
  }
  
  has_predictors <- !is.null(x) && NCOL(x) > 0L
  
  if (fast) {

    if (has_predictors) {
      xx <- cbind("(Intercept)" = rep.int(1, nobs), x)
    } else {
      xx <- matrix(
        rep.int(1, nobs),
        ncol = 1L,
        dimnames = list(NULL, "(Intercept)")
      )
    }
    
    fit <- stats::glm.fit(
      x = xx,
      y = y,
      family = family,
      weights = weights,
      offset = offset
    )
  } else {
    
    if (is.null(x) || NCOL(x) == 0) {
      data <- data.frame(y = y)
      
      if (isTRUE(has_offset)) {
        data$offset_ <- offset
        formula <- y ~ offset(offset_)
      } else {
        formula <- y ~ 1
      }
      
    } else {
      if (is.null(colnames(x)) || any(colnames(x) == "")) {
        stop("! Internal error: x must have non-empty column names.", call. = FALSE)
      }
      
      data <- data.frame(x, y = y, check.names = FALSE)
      
      rhs <- paste(sprintf("`%s`", colnames(x)), collapse = " + ")
      
      if (isTRUE(has_offset)) {
        data$offset_ <- offset
        rhs <- paste(rhs, "+ offset(offset_)")
      }
      
      formula <- stats::as.formula(paste("y ~", rhs))
    }

    
    fit <- stats::glm(
      formula = formula,
      data = data,
      family = family,
      weights = weights,
      x = TRUE,
      y = TRUE
    )
  }

  # account for estimation of variance parameter in gaussian models
  # computation as in logLik.glm using rank
  df <- if (fit$family$family == "gaussian") fit$rank + 1 else fit$rank
  
  
  # we need weighted rss for gaussian
  fit_weights <- fit$prior.weights
  if (is.null(fit_weights)) {
    fit_weights <- weights
  }
  
  if (length(fit_weights) != length(fit$residuals)) {
    stop(
      "Internal error: fitted residuals and weights have different lengths.",
      call. = FALSE
    )
  }
  
  list(
    fit = fit,
    # loglikelihood computed as in stats::logLik.glm
    logl = df - fit$aic / 2,
    coefficients = fit$coefficients,
    df = df,
    sse =  sum(fit_weights * fit$residuals^2, na.rm = TRUE),
    residuals = fit$residuals,
    weights = fit_weights,
    has_scale_parameter = fit$family$family == "gaussian"
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
#' `FALSE`uses the normal fitting routines ([survival::coxph()]) (used for
#'  the final output of `mfp2`).
#' @param has_offset logical indicating whether `offset` should be included in
#' the final formula-based fit when `fast = FALSE`.
#' @param strata,control,rownames,nocenter passed to [survival::coxph.fit()].
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
#' @keywords internal
#' @noRd
fit_cox <- function(x, 
                    y, 
                    strata, 
                    weights, 
                    offset, 
                    control, 
                    method, 
                    rownames, 
                    nocenter, 
                    fast = TRUE,
                    has_offset = FALSE) {
  
  # Set default for control
  if (is.null(control)) { 
    control <- survival::coxph.control()
  }
  
  has_predictors <- !is.null(x) && NCOL(x) > 0
  
  if (fast) {
    fit <- survival::coxph.fit(
      x = x,
      y = y, 
      strata = strata,
      weights = weights,
      offset = offset,
      control = control, 
      method = method, 
      rownames = rownames, 
      resid = TRUE,
      nocenter = nocenter
    )  
  } else {
    # construct appropriate formula incorporating offset and strata terms
    # cbinding y will lead to two variables: time and status
    
    if (!has_predictors) {
      d <- data.frame(y = y)
      rhs <- "1"
    } else {
      d <- data.frame(x, y = y, check.names = FALSE)
      
      if (is.null(colnames(x)) || any(colnames(x) == "")) {
        stop("! Internal error: x must have non-empty column names.", call. = FALSE)
      }
      
      rhs <- paste(sprintf("`%s`", colnames(x)), collapse = " + ")
    }
    
    # Add offset only when the model was structurally specified with one.
    # This distinguishes no-offset models from all-zero offset models.
    if (isTRUE(has_offset)) {
      d$offset_ <- offset
      rhs <- paste(rhs, "+ offset(offset_)")
    }
    
    if (!is.null(strata)) {
      d$strata_ <- strata
      rhs <- paste(rhs, "+ strata(strata_)")
    }
    
    ff <- stats::as.formula(paste("y ~", rhs))
    
    fit <- survival::coxph(
      ff,
      data = d, 
      weights = weights, 
      control = control,
      method = method, 
      nocenter = nocenter, 
      x = TRUE,
      y = TRUE
    )
  }
  
  logl <- if (length(fit$loglik) >= 2) {
    fit$loglik[2]
  } else {
    fit$loglik[1]
  }
  
  fit_weights <- fit$prior.weights
  
  if (is.null(fit_weights)) {
    fit_weights <- weights
  }
  
  list(
    fit = fit, 
    logl = logl,
    coefficients = fit$coefficients, 
    # sometimes coefficients can be NA
    # for example when including same variables in the model
    df = length(fit$coefficients[!is.na(fit$coefficients)]), 
    weights = fit_weights,
    sse = sum(fit_weights * fit$residuals^2, na.rm = TRUE),
    residuals = fit$residuals,
    has_scale_parameter = FALSE
  )
}