#' Function that fits models supported by `mfp2`
#'
#' Fits generalized linear models and Cox proportional-hazards models.
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
#' @param family a character string specifying the GLM family to be used, or
#' "cox" for Cox models. The default family is set to 'Gaussian'.
#' @param strata,control,weights,offset,rownames,nocenter parameters for Cox
#' or glm. See [survival::coxph()] or [stats::glm()] for details.
#' @param fast passed to \code{fit_glm()} and \code{fit_cox()}.
#' @param calculate_fit_statistics logical. If `TRUE`, return the
#' family-specific null and fitted-model statistics used for final reporting.
#' Candidate-selection fits leave this `FALSE`.
#' @param calculate_gaussian_deviance logical. If `TRUE` for a Gaussian GLM,
#' compute and retain the scalar Stata-style Gaussian deviance used by F-tests.
#' Leave this `FALSE` for likelihood-ratio, AIC, and BIC selection so candidate
#' fits do not copy residual and weight vectors into the lightweight wrapper.
#' @param keep_fit logical. If `TRUE`, retain the underlying fitted object.
#' Defaults to `!fast`, so ordinary final fits are retained and fast candidate
#' fits are discarded unless a caller explicitly needs the fitted object.
#' @param keep_fitted_values logical. If `TRUE` for a GLM, retain only the
#' fitted-values vector in the lightweight wrapper. This is useful for internal
#' candidate searches that need predictions but not the complete backend fit.
#' @param fitter GLM fitting backend; ignored for Cox models.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param x_has_intercept internal logical. If TRUE, GLM fast fitting treats
#' `x` as an already-intercepted model matrix. Must be FALSE for Cox models.
#'
#' @return
#' A list with the following components:
#'
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `rank`: fitted regression rank.
#' * `deviance_gaussian`: when `calculate_gaussian_deviance = TRUE` for a
#'   Gaussian GLM, the scalar Stata-style Gaussian deviance.
#' * `null_deviance`: when `calculate_fit_statistics = TRUE`, family-specific
#'   deviance of the null model. For GLMs,
#'   this is the `null.deviance` returned by [stats::glm.fit()] or
#'   [stats::glm()]. For Cox models, it is minus twice the null partial
#'   log-likelihood.
#' * `model_deviance`: when `calculate_fit_statistics = TRUE`, family-specific
#'   deviance of the fitted model. For GLMs,
#'   this is the `deviance` returned by [stats::glm.fit()] or [stats::glm()].
#'   For Cox models, it is minus twice the fitted partial log-likelihood.
#' * `null_logl`: when `calculate_fit_statistics = TRUE` for a Cox model, the
#'   null partial log-likelihood.
#' * `fitted_values`: for GLMs, when `keep_fitted_values = TRUE`, the backend
#'   fitted-values vector without retaining the complete fitted object.
#' * `fit`: when `keep_fit = TRUE`, the object returned by the fitting
#'   procedure.
#' * `transformed_to_model_columns`: for full fits, a named character vector
#'   mapping source design columns to the exact fitted coefficient names.
#'
#' @importFrom stats family
#' @keywords internal
#' @noRd
fit_model <- function(x,
                      y,
                      family,
                      family_string,
                      weights = NULL,
                      offset = NULL,
                      method = NULL,
                      strata = NULL,
                      control = NULL,
                      rownames = NULL,
                      nocenter = NULL,
                      fast = TRUE,
                      calculate_fit_statistics = FALSE,
                      calculate_gaussian_deviance = FALSE,
                      keep_fit = !fast,
                      has_offset = FALSE,
                      x_has_intercept = FALSE,
                      fitter = "base",
                      keep_fitted_values = FALSE) {
  # Set column names if not provided
  if (!is.null(dim(x)) && is.null(colnames(x))) {
    colnames(x) <- colnames(x, do.NULL = FALSE)
  }

  if (identical(family_string, "cox")) {
    if (isTRUE(keep_fitted_values)) {
      stop(
        "Internal error: keep_fitted_values is available only for GLMs.",
        call. = FALSE
      )
    }

    if (isTRUE(x_has_intercept)) {
      stop(
        "Internal error: Cox model matrix must not include an intercept.",
        call. = FALSE
      )
    }

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
      calculate_fit_statistics = calculate_fit_statistics,
      keep_fit = keep_fit,
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
      calculate_fit_statistics = calculate_fit_statistics,
      calculate_gaussian_deviance = calculate_gaussian_deviance,
      keep_fit = keep_fit,
      keep_fitted_values = keep_fitted_values,
      fitter = fitter,
      has_offset = has_offset,
      x_has_intercept = x_has_intercept,
      family_string = family_string
    )
  }


  # Full formula-based fits may quote non-syntactic source-column names in
  # their coefficient vectors. Store the exact positional relationship once at
  # fit time so prediction code never has to add or remove backticks.
  if (!isTRUE(fast)) {
    source_columns <- if (is.null(x) || NCOL(x) == 0L) {
      character(0L)
    } else {
      colnames(x)
    }

    if (isTRUE(x_has_intercept) && length(source_columns) > 0L &&
        identical(source_columns[[1L]], "(Intercept)")) {
      source_columns <- source_columns[-1L]
    }

    fitted_columns <- names(fit$coefficients)
    if (is.null(fitted_columns)) {
      fitted_columns <- character(0L)
    }
    fitted_columns <- setdiff(fitted_columns, "(Intercept)")

    if (length(source_columns) != length(fitted_columns)) {
      stop(
        "Internal error: fitted coefficient names do not align with source columns.",
        call. = FALSE
      )
    }

    fit$transformed_to_model_columns <- stats::setNames(
      fitted_columns,
      source_columns
    )
  }

  fit
}

#' Assemble a Design Matrix from Reusable Column Blocks
#'
#' Allocates the destination matrix once and fills contiguous column ranges
#' from the supplied blocks. This avoids constructing intermediate matrices
#' when an intercept and multiple predictor blocks must be combined.
#'
#' @param blocks list of matrix-like predictor blocks. NULL and zero-column
#'   blocks are ignored.
#' @param nobs number of rows in the resulting matrix.
#' @param intercept logical; prepend a column named `(Intercept)` when TRUE.
#'
#' @return A numeric matrix containing the optional intercept and all blocks in
#'   list order.
#' @keywords internal
#' @noRd
assemble_design_matrix <- function(blocks, nobs, intercept = FALSE) {
  keep <- vapply(
    blocks,
    function(block) !is.null(block) && NCOL(block) > 0L,
    logical(1L)
  )
  blocks <- blocks[keep]

  if (length(blocks) > 0L) {
    block_nrows <- vapply(blocks, NROW, integer(1L))
    if (any(block_nrows != nobs)) {
      stop("Internal error: design-matrix blocks have inconsistent rows.", call. = FALSE)
    }
  }

  block_ncols <- if (length(blocks) > 0L) {
    vapply(blocks, NCOL, integer(1L))
  } else {
    integer(0L)
  }

  block_names <- unlist(
    lapply(blocks, function(block) {
      names <- colnames(block)
      if (is.null(names)) {
        names <- colnames(block, do.NULL = FALSE)
      }
      names
    }),
    use.names = FALSE
  )

  row_names <- NULL
  if (length(blocks) > 0L) {
    row_name_index <- which(vapply(
      blocks,
      function(block) !is.null(rownames(block)),
      logical(1L)
    ))
    if (length(row_name_index) > 0L) {
      row_names <- rownames(blocks[[row_name_index[[1L]]]])
    }
  }

  intercept_width <- as.integer(isTRUE(intercept))
  result <- matrix(
    0,
    nrow = nobs,
    ncol = intercept_width + sum(block_ncols),
    dimnames = list(
      row_names,
      c(if (isTRUE(intercept)) "(Intercept)", block_names)
    )
  )

  next_column <- 1L
  if (isTRUE(intercept)) {
    result[, 1L] <- 1
    next_column <- 2L
  }

  if (length(blocks) > 0L) {
    for (block_index in seq_along(blocks)) {
      columns <- seq.int(
        from = next_column,
        length.out = block_ncols[[block_index]]
      )
      result[, columns] <- blocks[[block_index]]
      next_column <- next_column + block_ncols[[block_index]]
    }
  }

  result
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
#' @param calculate_fit_statistics logical. If `TRUE`, return the null and
#' fitted-model deviances required for final reporting.
#' @param calculate_gaussian_deviance logical. If `TRUE` for a Gaussian GLM,
#' compute and retain the scalar Stata-style Gaussian deviance required by
#' F-tests. Residual and weight vectors are not retained in the wrapper.
#' @param keep_fit logical. If `TRUE`, retain the underlying fitted object.
#' @param keep_fitted_values logical. If `TRUE`, retain only the backend
#' fitted-values vector in the lightweight return object.
#' @param fitter GLM fitting backend for the matrix fast path.
#' @param family_string Normalized family name supplied by `fit_model()`.
#' @param has_offset logical indicating whether `offset` should be included in
#' the final formula-based fit when `fast = FALSE`.
#' @param x_has_intercept internal logical. If TRUE, `x` is already the full
#' model matrix including an intercept column named `"(Intercept)"`.
#'
#' @return
#' A list with the following components:
#'
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients.
#' * `df`: number of parameters (degrees of freedom).
#' * `rank`: fitted regression rank.
#' * `deviance_gaussian`: when `calculate_gaussian_deviance = TRUE`, the
#'   scalar Stata-style Gaussian deviance.
#' * `null_deviance`: when `calculate_fit_statistics = TRUE`, the null-model
#'   deviance returned by the fitted GLM.
#' * `model_deviance`: when `calculate_fit_statistics = TRUE`, the residual
#'   deviance returned by the fitted GLM.
#' * `fitted_values`: when `keep_fitted_values = TRUE`, the backend
#'   fitted-values vector without retaining that full object.
#' * `fit`: when `keep_fit = TRUE`, the fitted model object.
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
                    calculate_fit_statistics = FALSE,
                    calculate_gaussian_deviance = FALSE,
                    keep_fit = !fast,
                    has_offset = FALSE,
                    x_has_intercept = FALSE,
                    fitter = "base",
                    family_string = NULL,
                    keep_fitted_values = FALSE) {
  if (is.null(family_string)) {
    family_string <- if (is.character(family) && length(family) == 1L) {
      family
    } else {
      family$family
    }
  }
  nobs <- NROW(y)

  has_predictors <- !is.null(x) && NCOL(x) > 0L

  if (isTRUE(x_has_intercept)) {
    if (!has_predictors) {
      stop(
        "Internal error: x_has_intercept = TRUE requires a non-empty matrix.",
        call. = FALSE
      )
    }

    if (is.null(colnames(x)) || colnames(x)[1L] != "(Intercept)") {
      stop(
        "Internal error: x_has_intercept = TRUE requires first column '(Intercept)'.",
        call. = FALSE
      )
    }
  }

  if (fast) {

    if (isTRUE(x_has_intercept)) {
      xx <- x
    } else {
      xx <- assemble_design_matrix(
        blocks = list(x),
        nobs = nobs,
        intercept = TRUE
      )
    }

    fit <- if (identical(family_string, "negbin")) {
      fit_glm_fastglm_nb(xx, y, weights, offset)
    } else if (identical(fitter, "fastglm")) {
      fit_glm_fastglm(xx, y, family, weights, offset)
    } else {
      stats::glm.fit(
        x = xx,
        y = y,
        family = family,
        weights = weights,
        offset = offset
      )
    }

  } else if (identical(family_string, "negbin")) {
    # Negative-binomial fitting has no stats::glm() formula equivalent because
    # theta is estimated jointly. Use the same rank-revealing matrix fitter for
    # the single final fit while leaving all other families on stats::glm().
    if (isTRUE(x_has_intercept)) {
      xx <- x
    } else {
      xx <- assemble_design_matrix(
        blocks = list(x),
        nobs = nobs,
        intercept = TRUE
      )
    }
    fit <- fit_glm_fastglm_nb(xx, y, weights, offset)

  } else {

    x_formula <- if (isTRUE(x_has_intercept)) {
      x[, -1L, drop = FALSE]
    } else {
      x
    }

    has_formula_predictors <- !is.null(x_formula) && NCOL(x_formula) > 0L

    if (has_formula_predictors) {
      if (is.null(colnames(x_formula)) || any(colnames(x_formula) == "")) {
        stop(
          "Internal error: x must have non-empty column names.",
          call. = FALSE
        )
      }

      data <- data.frame(x_formula, check.names = FALSE)
      rhs <- paste(sprintf("`%s`", colnames(x_formula)), collapse = " + ")

    } else {
      data <- data.frame(row.names = seq_len(nobs))
      rhs <- "1"
    }

    if (is.matrix(y)) {
      response_cols <- c("..mfp2_successes", "..mfp2_failures")

      if (any(response_cols %in% names(data))) {
        stop(
          "Internal error: response column names conflict with design matrix names.",
          call. = FALSE
        )
      }

      data[[response_cols[1L]]] <- y[, 1L]
      data[[response_cols[2L]]] <- y[, 2L]

      lhs <- sprintf(
        "cbind(%s, %s)",
        response_cols[1L],
        response_cols[2L]
      )

    } else {
      response_col <- "..mfp2_y"

      if (response_col %in% names(data)) {
        stop(
          "Internal error: response column name conflicts with design matrix names.",
          call. = FALSE
        )
      }

      data[[response_col]] <- y
      lhs <- response_col
    }

    if (isTRUE(has_offset)) {
      data$offset_ <- offset

      rhs <- if (identical(rhs, "1")) {
        "offset(offset_)"
      } else {
        paste(rhs, "+ offset(offset_)")
      }
    }

    formula <- stats::as.formula(paste(lhs, "~", rhs))

    fit <- stats::glm(
      formula = formula,
      data = data,
      family = family,
      weights = weights,
      x = TRUE,
      y = TRUE
    )
  }

  # Gaussian scale and negative-binomial theta are estimated nuisance
  # parameters in addition to the regression rank.
  is_gaussian <- identical(family_string, "gaussian")
  is_negbin <- identical(family_string, "negbin")
  df <- fit$rank + as.integer(is_gaussian || is_negbin)

  # The Stata-style Gaussian deviance is needed only for Gaussian F-tests.
  # Compute it here while the fitted object is available, then retain only the
  # scalar. This avoids copying observation-length residual and weight vectors
  # into every candidate fit used by LR, AIC, or BIC selection.
  gaussian_deviance <- NA_real_

  if (is_gaussian && isTRUE(calculate_gaussian_deviance)) {
    fit_residuals <- fit$residuals
    fit_weights <- fit$prior.weights
    if (is.null(fit_weights)) {
      fit_weights <- weights
    }

    if (!is.numeric(fit_residuals) || !is.numeric(fit_weights) ||
        length(fit_weights) != length(fit_residuals)) {
      stop(
        "Internal error: Gaussian residuals and weights are unavailable or ",
        "have different lengths.",
        call. = FALSE
      )
    }

    gaussian_deviance <- deviance_gaussian(
      residuals = fit_residuals,
      weights = fit_weights
    )

    if (is.null(gaussian_deviance)) {
      gaussian_deviance <- NA_real_
    }
  }

  result <- list(
    # fastglm_nb() reports the maximized twice-log-likelihood directly. Its
    # stored AIC does not include theta, so deriving logL from AIC would be off
    # by one parameter. Other GLMs retain stats::logLik.glm()'s convention.
    logl = if (is_negbin) fit$twologlik / 2 else df - fit$aic / 2,
    coefficients = fit$coefficients,
    rank = unname(fit$rank),
    df = df
  )

  if (is_gaussian && isTRUE(calculate_gaussian_deviance)) {
    result$deviance_gaussian <- gaussian_deviance
  }

  if (isTRUE(keep_fitted_values)) {
    # Some internal searches need only the n-length fitted-value vector.
    # Retain that lightweight result before the backend fit goes out of scope
    # instead of forcing callers to retain the complete fit object.
    fitted_values <- fit$fitted.values
    if (!is.numeric(fitted_values) || length(fitted_values) != nobs) {
      stop(
        "Internal error: GLM backend did not return usable fitted values.",
        call. = FALSE
      )
    }
    result$fitted_values <- fitted_values
  }

  if (isTRUE(keep_fit)) {
    result$fit <- fit
  }

  if (isTRUE(calculate_fit_statistics)) {
    # Use the family-specific deviances already computed by glm.fit()/glm().
    # These are not, in general, equal to minus twice the log-likelihood.
    result$null_deviance <- unname(fit$null.deviance)
    result$model_deviance <- unname(fit$deviance)
  }

  result
}

#' Check whether the optional fastglm package is available
#' @keywords internal
#' @noRd
fastglm_available <- function() {
  requireNamespace("fastglm", quietly = TRUE)
}

#' Return functions exported by the installed fastglm package
#' @keywords internal
#' @noRd
fastglm_exports <- function() {
  getNamespaceExports("fastglm")
}

#' Resolve the effective GLM fitting backend once per top-level fit
#'
#' Ordinary GLMs fall back to the base fitter with one warning when fastglm is
#' unavailable or does not expose the required matrix fitter. Negative-binomial
#' models have no base fallback and therefore fail before candidate fitting.
#' Cox models do not use a GLM fitter and always resolve to "base".
#'
#' @param fitter Requested fitter, "base" or "fastglm".
#' @param family_string Normalized family name.
#'
#' @return The effective fitter, either "base" or "fastglm".
#' @keywords internal
#' @noRd
resolve_fitter <- function(fitter, family_string) {
  fitter <- match.arg(fitter, c("base", "fastglm"))

  if (identical(family_string, "cox")) {
    return("base")
  }

  if (identical(family_string, "negbin") &&
      !identical(fitter, "fastglm")) {
    stop(
      "`family = \"negbin\"` is available only with `fitter = \"fastglm\"`; ",
      "stats::glm.fit() cannot estimate the negative-binomial dispersion parameter.",
      call. = FALSE
    )
  }

  if (identical(fitter, "base")) {
    return("base")
  }

  if (!fastglm_available()) {
    if (identical(family_string, "negbin")) {
      stop(
        "`family = \"negbin\"` requires the optional 'fastglm' package. ",
        "Install it with install.packages(\"fastglm\").",
        call. = FALSE
      )
    }

    warning(
      "`fitter = \"fastglm\"` was requested but the 'fastglm' package is ",
      "not installed; using `stats::glm.fit()` instead.",
      call. = FALSE
    )
    return("base")
  }

  exports <- fastglm_exports()
  has_glm <- "fastglm" %in% exports
  has_nb <- "fastglm_nb" %in% exports

  if (!has_glm) {
    if (identical(family_string, "negbin")) {
      stop(
        "The installed 'fastglm' package does not export `fastglm()`; ",
        "update fastglm before using `family = \"negbin\"`.",
        call. = FALSE
      )
    }

    warning(
      "The installed 'fastglm' package does not export `fastglm()`; ",
      "using `stats::glm.fit()` instead.",
      call. = FALSE
    )
    return("base")
  }

  if (identical(family_string, "negbin") && !has_nb) {
    stop(
      "The installed 'fastglm' package does not export `fastglm_nb()`; ",
      "update fastglm before using `family = \"negbin\"`.",
      call. = FALSE
    )
  }

  "fastglm"
}

#' Test whether a value is one finite numeric scalar
#' @keywords internal
#' @noRd
is_finite_numeric_scalar <- function(value) {
  is.numeric(value) && length(value) == 1L && !is.na(value) && is.finite(value)
}

#' Fit a GLM with fastglm and safely fall back to stats::glm.fit
#'
#' Package availability is resolved once by `resolve_fitter()` before entering
#' the candidate-search path. This helper retains only candidate-specific error
#' and return-value checks.
#'
#' @keywords internal
#' @noRd
fit_glm_fastglm <- function(x, y, family, weights, offset) {
  fallback <- function() {
    stats::glm.fit(
      x = x,
      y = y,
      family = family,
      weights = weights,
      offset = offset
    )
  }

  fit_error <- NULL
  fit <- tryCatch(
    fastglm::fastglm(
      x = x,
      y = y,
      family = family,
      weights = weights,
      offset = offset,
      method = 0L
    ),
    error = function(e) {
      fit_error <<- e
      NULL
    }
  )

  family_name <- if (!is.null(fit) && is.list(fit$family)) {
    fit$family$family
  } else {
    NULL
  }
  is_gaussian <- is.character(family_name) &&
    length(family_name) == 1L && identical(family_name, "gaussian")

  fit_is_usable <- !is.null(fit) &&
    all(vapply(
      list(fit$deviance, fit$aic, fit$rank, fit$null.deviance),
      is_finite_numeric_scalar,
      logical(1L)
    )) &&
    is.character(family_name) &&
    length(family_name) == 1L &&
    is.numeric(fit$coefficients)

  # Only Gaussian candidate comparisons require residuals and weights for the
  # package's F-test deviance. Other GLM families are compared by likelihood.
  if (fit_is_usable && is_gaussian) {
    fit_weights <- fit$prior.weights
    if (is.null(fit_weights)) fit_weights <- fit$weights

    fit_is_usable <- is.numeric(fit$residuals) &&
      length(fit$residuals) == NROW(y) &&
      is.numeric(fit_weights) &&
      length(fit_weights) == NROW(y)
  }

  if (!fit_is_usable) {
    if (!is.null(fit_error)) {
      warning(
        "fastglm failed for this candidate (", conditionMessage(fit_error),
        "); falling back to stats::glm.fit().",
        call. = FALSE
      )
    } else {
      warning(
        "fastglm returned a non-finite fit for this candidate; falling back ",
        "to stats::glm.fit().",
        call. = FALSE
      )
    }
    return(fallback())
  }

  fit
}

#' Fit a negative-binomial GLM with fastglm
#'
#' Package and export availability are resolved once by `resolve_fitter()`.
#' Negative-binomial fitting has no base fallback, so candidate-specific errors
#' remain fatal.
#'
#' @keywords internal
#' @noRd
fit_glm_fastglm_nb <- function(x, y, weights, offset) {
  fit <- tryCatch(
    fastglm::fastglm_nb(
      x = x,
      y = y,
      weights = weights,
      offset = offset,
      method = 0L
    ),
    error = function(e) {
      stop(
        "fastglm could not fit the requested negative-binomial model: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # fastglm_nb() does not return glm.fit()-style working residuals. They
  # are not needed for negative-binomial candidate selection, which is based on
  # the maximized likelihood, rank, deviance, and theta, so do not fabricate or
  # store them in the fit_model() wrapper.

  # fastglm_nb()'s stored aic counts the regression rank but not the jointly
  # estimated dispersion parameter. Recalculate it before validating the fit.
  if (is_finite_numeric_scalar(fit$twologlik) &&
      is_finite_numeric_scalar(fit$rank)) {
    fit$aic <- -fit$twologlik + 2 * (fit$rank + 1L)
  }

  scalar_fields_are_finite <- all(vapply(
    list(
      fit$deviance, fit$aic, fit$rank,
      fit$theta, fit$twologlik
    ),
    is_finite_numeric_scalar,
    logical(1L)
  ))

  vector_fields_are_usable <-
    is.numeric(fit$coefficients) &&
    is.numeric(fit$fitted.values) &&
    length(fit$fitted.values) == NROW(y) &&
    all(is.finite(fit$fitted.values)) &&
    is.numeric(fit$linear.predictors) &&
    length(fit$linear.predictors) == NROW(y) &&
    all(is.finite(fit$linear.predictors)) &&
    is.numeric(fit$prior.weights) &&
    length(fit$prior.weights) == NROW(y) &&
    all(is.finite(fit$prior.weights))

  family_is_usable <-
    is.list(fit$family) &&
    is.character(fit$family$family) &&
    length(fit$family$family) == 1L

  if (!(scalar_fields_are_finite && vector_fields_are_usable &&
        family_is_usable)) {
    stop(
      "fastglm returned an incomplete or non-finite negative-binomial fit; ",
      "no base fitter fallback is available for `family = \"negbin\"`.",
      call. = FALSE
    )
  }

  # The rank + 1 AIC convention above matches MASS::glm.nb() and counts
  # theta without over-counting aliased regression columns.
  fit
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
#' @param calculate_fit_statistics logical. If `TRUE`, return the null
#' partial log-likelihood and the null and fitted-model statistics required for
#' final reporting.
#' @param keep_fit logical. If `TRUE`, retain the underlying fitted object.
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
#' * `rank`: fitted regression rank.
#' * `null_logl`: when `calculate_fit_statistics = TRUE`, the null partial
#'   log-likelihood.
#' * `null_deviance`: when `calculate_fit_statistics = TRUE`, minus twice the
#'   null partial log-likelihood.
#' * `model_deviance`: when `calculate_fit_statistics = TRUE`, minus twice the
#'   fitted partial log-likelihood.
#' * `fit`: when `keep_fit = TRUE`, the fitted model object.
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
                    calculate_fit_statistics = FALSE,
                    keep_fit = !fast,
                    has_offset = FALSE) {

  # Set default for control
  if (is.null(control)) {
    control <- survival::coxph.control()
  }
  has_predictors <- !is.null(x) && NCOL(x) > 0

  # coxph.fit() requires integer stratum identifiers. Many internal fast-fit
  # calls already pass an integer vector, so avoid repeatedly coercing an
  # observation-length object and creating unnecessary allocation/GC pressure.
  istrata <- if (is.null(strata)) {
    NULL
  } else if (is.integer(strata)) {
    strata
  } else {
    as.integer(strata)
  }

  if (fast) {
    fit <- survival::coxph.fit(
      x = x,
      y = y,
      strata = istrata,
      offset = offset,
      control = control,
      weights = weights,
      method = method,
      rownames = rownames,
      resid = FALSE,
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

  # coxph.fit()/coxph() normally return the null and fitted partial
  # log-likelihoods in positions 1 and 2, respectively. Preserve both so the
  # caller can report null and fitted-model deviances without another fit.
  if (length(fit$loglik) >= 2L) {
    null_logl <- fit$loglik[1L]
    model_logl <- fit$loglik[2L]
  } else {
    # Defensive fallback for an irregular fitted object. The available value is
    # treated as the fitted partial log-likelihood; the null value is unknown.
    null_logl <- NA_real_
    model_logl <- fit$loglik[1L]
  }

  model_df <- length(fit$coefficients[!is.na(fit$coefficients)])
  result <- list(
    logl = model_logl,
    coefficients = fit$coefficients,
    # Sometimes coefficients can be NA, for example when duplicate or
    # linearly dependent predictors are included in the model.
    rank = model_df,
    df = model_df,
    # Duplicated residual/weight vectors are required only by Gaussian F-tests,
    # never by Cox likelihood comparisons.
    weights = NULL,
    residuals = NULL
  )

  if (isTRUE(keep_fit)) {
    result$fit <- fit
  }

  if (isTRUE(calculate_fit_statistics)) {
    # The null partial log-likelihood is computed by coxph at the start of
    # Newton-Raphson and can therefore be reported without another fit.
    result$null_logl <- if (is.finite(null_logl)) {
      unname(null_logl)
    } else {
      NA_real_
    }
    result$null_deviance <- if (is.finite(null_logl)) {
      unname(-2 * null_logl)
    } else {
      NA_real_
    }
    result$model_deviance <- unname(-2 * model_logl)
  }

  result
}