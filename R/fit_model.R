#' Function that fits models supported by `mfp2`
#'
#' Fits likelihood-based generalized linear, multinomial logistic, Cox
#' proportional-hazards, parametric survival, and Fine--Gray models.
#'
#' @details
#' Computations dispatch to the corresponding matrix-level candidate fitter and
#' retain the native high-level model class for the final fit.
#'
#' @param x a matrix of predictors (excluding intercept) with column names.
#' If column names are not provided they are set according to
#' `colnames(x, do.NULL = FALSE)`.
#' @param y Response variable. For GLMs, this may be a numeric vector, a factor
#' response accepted by [stats::glm()], or for binomial models a two-column
#' matrix of grouped integer counts `cbind(successes, failures)`. Poisson and
#' negative-binomial responses are nonnegative integer counts. Multinomial
#' responses are factors, character/numeric class-label vectors, or class-count
#' matrices. Survival families use
#' a family-appropriate [survival::Surv()] object.
#' @param method a character string specifying the method for tie handling.
#' See [survival::coxph()].
#' @param family Resolved GLM, multinomial, Cox, survreg, or Fine--Gray family.
#' @param strata,weights,offset,rownames,nocenter family-specific fitting
#' parameters. See [stats::glm()], [survival::coxph()], and
#' [survival::survreg()] for details.
#' @param control Required family-specific control object, normalized once by
#' the public fitting boundary before repeated candidate fits begin.
#' @param fast logical selecting the matrix-level candidate path rather than the
#' retained native final fit.
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
#' @param keep_fitted_values logical. If `TRUE` for a GLM or multinomial model,
#' retain only the fitted-values vector in the lightweight wrapper. This is
#' useful for internal candidate searches that need predictions but not the
#' complete backend fit.
#' @param keep_coefficients logical. For multinomial candidate fits, controls
#' whether the per-logit coefficient matrix and flattened coefficient vector are
#' materialized. Candidate searches that use only likelihood and degrees of
#' freedom set this to `FALSE`.
#' @param fitter GLM fitting backend; ignored for survival models.
#' @param has_offset logical indicating whether an offset was specified before
#' missing offsets were replaced by zeros internally.
#' @param x_has_intercept internal logical. If TRUE, GLM fast fitting treats
#' `x` as an already-intercepted model matrix. Must be FALSE for Cox models.
#' @param reserved_names Optional character vector of user-facing names that
#' must not be reused for package-created formula columns in a final refit.
#' Ignored by fast matrix fits.
#' @param multinomial_optimizer Optional precomputed mask/weight structure for
#' a batch of same-shaped multinomial candidates.
#'
#' @return
#' A list with the following components:
#'
#' * `logl`: the log likelihood of the fitted model.
#' * `coefficients`: regression coefficients. Omitted for multinomial candidate
#'   fits when `keep_coefficients = FALSE`.
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
#' * `fitted_values`: for GLMs and multinomial models, when
#'   `keep_fitted_values = TRUE`, the backend fitted-values vector or matrix
#'   without retaining the complete fitted object.
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
                      control,
                      rownames = NULL,
                      nocenter = c(-1, 0, 1),
                      fast = TRUE,
                      calculate_fit_statistics = FALSE,
                      calculate_gaussian_deviance = FALSE,
                      keep_fit = !fast,
                      has_offset = FALSE,
                      x_has_intercept = FALSE,
                      fitter = "base",
                      keep_fitted_values = FALSE,
                      keep_coefficients = TRUE,
                      multinomial_optimizer = NULL,
                      reserved_names = character()) {
  # Match coxph() by default: 0/1 and -1/0/1 design columns are not
  # internally recentered. Public callers pass this value explicitly;
  # keeping the same internal default also protects direct helper calls.
  # Set column names if not provided
  if (!is.null(dim(x)) && is.null(colnames(x))) {
    colnames(x) <- colnames(x, do.NULL = FALSE)
  }

  if (identical(family_string, "multinomial")) {
    if (isTRUE(calculate_gaussian_deviance)) {
      stop("Internal error: Gaussian deviance is unavailable for multinomial models.", call. = FALSE)
    }
    fit <- fit_multinomial(
      x = x,
      family = family,
      control = control,
      fast = fast,
      calculate_fit_statistics = calculate_fit_statistics,
      keep_fit = keep_fit,
      keep_fitted_values = keep_fitted_values,
      keep_coefficients = keep_coefficients,
      x_has_intercept = x_has_intercept,
      optimizer_structure = multinomial_optimizer
    )
  } else if (identical(family_string, "cox")) {
    if (isTRUE(keep_fitted_values)) {
      stop(
        "Internal error: keep_fitted_values is unavailable for Cox models.",
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
      has_offset = has_offset,
      reserved_names = reserved_names
    )
  } else if (identical(family_string, "finegray")) {
    if (isTRUE(keep_fitted_values)) {
      stop(
        "Internal error: keep_fitted_values is unavailable for Fine--Gray models.",
        call. = FALSE
      )
    }
    if (isTRUE(x_has_intercept)) {
      stop(
        "Internal error: Fine--Gray model matrix must not include an intercept.",
        call. = FALSE
      )
    }

    fit <- fit_finegray(
      x = x,
      family = family,
      offset = offset,
      control = control,
      method = method,
      nocenter = nocenter,
      fast = fast,
      calculate_fit_statistics = calculate_fit_statistics,
      keep_fit = keep_fit,
      has_offset = has_offset,
      reserved_names = reserved_names
    )
  } else if (identical(family_string, "ordinal")) {
    if (isTRUE(keep_fitted_values)) {
      stop(
        "Internal error: keep_fitted_values is unavailable for ordinal models.",
        call. = FALSE
      )
    }

    fit <- fit_ordinal(
      x = x,
      family = family,
      control = control,
      fast = fast,
      calculate_fit_statistics = calculate_fit_statistics,
      keep_fit = keep_fit,
      x_has_intercept = x_has_intercept,
      reserved_names = reserved_names
    )
  } else if (identical(family_string, "survreg")) {
    if (isTRUE(keep_fitted_values)) {
      stop(
        "Internal error: keep_fitted_values is unavailable for survreg models.",
        call. = FALSE
      )
    }

    fit <- fit_survreg(
      x = x,
      y = y,
      family = family,
      weights = weights,
      offset = offset,
      control = control,
      fast = fast,
      calculate_fit_statistics = calculate_fit_statistics,
      keep_fit = keep_fit,
      has_offset = has_offset,
      x_has_intercept = x_has_intercept,
      reserved_names = reserved_names
    )
  } else {
    fit <- fit_glm(
      y = y,
      x = x,
      family = family,
      weights = weights,
      offset = offset,
      control = control,
      fast = fast,
      calculate_fit_statistics = calculate_fit_statistics,
      calculate_gaussian_deviance = calculate_gaussian_deviance,
      keep_fit = keep_fit,
      keep_fitted_values = keep_fitted_values,
      fitter = fitter,
      has_offset = has_offset,
      x_has_intercept = x_has_intercept,
      family_string = family_string,
      reserved_names = reserved_names
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

    fitted_columns <- if (family_string %in% c("multinomial", "ordinal")) {
      # Multinomial stores a per-logit matrix and ordinal carries k-1 intercepts
      # ahead of the slopes; in both cases the slope coefficients are named by
      # the design columns, so map source columns to themselves.
      source_columns
    } else {
      names(fit$coefficients)
    }
    if (is.null(fitted_columns)) fitted_columns <- character(0L)
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


#' Validate a Fitted Model Before Its Criteria Enter MFP Selection
#'
#' Candidate models are never silently discarded or automatically refitted.
#' Every model that participates in an MFP decision must have a finite
#' likelihood/df pair, and a backend that explicitly reports non-convergence
#' makes that fit unusable. The checks are scalar and therefore negligible
#' compared with the model fit itself.
#'
#' @param logl Fitted model log-likelihood.
#' @param df Model degrees of freedom used by MFP comparisons.
#' @param converged Optional logical convergence flag returned by the backend.
#' @param gaussian_deviance Optional Gaussian deviance used by F-tests.
#' @param require_gaussian_deviance Whether a finite Gaussian deviance is
#' required for the current comparison path.
#' @param family_string Normalized model-family label used in error messages.
#' @param fast Whether this is a repeated candidate fit (`TRUE`) or the final
#' model refit (`FALSE`).
#'
#' @return Invisibly `TRUE`; otherwise stops with a targeted error.
#' @keywords internal
#' @noRd
validate_mfp_fit_result <- function(logl,
                                    df,
                                    converged = NULL,
                                    gaussian_deviance = NULL,
                                    require_gaussian_deviance = FALSE,
                                    family_string = NULL,
                                    fast = TRUE) {
  stage <- if (isTRUE(fast)) "candidate" else "final"
  family_label <- if (is.character(family_string) &&
                      length(family_string) == 1L &&
                      !is.na(family_string) && nzchar(family_string)) {
    paste0(" ", family_string)
  } else {
    ""
  }

  # Respect an explicit convergence flag from glm/glm.fit/fastglm rather than
  # attempting a second fit with different controls. Re-fitting would alter
  # runtime and could change the candidate search requested by the user.
  if (!is.null(converged)) {
    # Some compiled backends may expose a scalar 0/1 flag rather than a native
    # R logical. Accept that representation, but reject any ambiguous value.
    if (is.numeric(converged) && length(converged) == 1L &&
        !is.na(converged) && is.finite(converged) &&
        converged %in% c(0, 1)) {
      converged <- as.logical(converged)
    }

    if (!is.logical(converged) || length(converged) != 1L || is.na(converged)) {
      stop(
        "Internal error: model backend returned an invalid convergence flag.",
        call. = FALSE
      )
    }
    if (!isTRUE(converged)) {
      stop(
        "The ", stage, family_label,
        " model did not converge; MFP fitting cannot use an unconverged fit. ",
        "Adjust the fitting controls or inspect the data.",
        call. = FALSE
      )
    }
  }

  criterion_label <- if (isTRUE(fast)) "candidate criterion" else "fit criterion"

  if (!is.numeric(logl) || length(logl) != 1L || !is.finite(logl)) {
    stop(
      "The ", stage, family_label,
      " model produced a non-finite log-likelihood; MFP fitting cannot ",
      "continue with an invalid ", criterion_label, ".",
      call. = FALSE
    )
  }

  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 0) {
    stop(
      "The ", stage, family_label,
      " model produced invalid degrees of freedom; MFP fitting cannot ",
      "continue with an invalid ", criterion_label, ".",
      call. = FALSE
    )
  }

  if (isTRUE(require_gaussian_deviance) &&
      (!is.numeric(gaussian_deviance) ||
       length(gaussian_deviance) != 1L ||
       !is.finite(gaussian_deviance))) {
    stop(
      "The ", stage,
      " Gaussian model produced a non-finite deviance required for the F-test; ",
      "MFP fitting cannot continue.",
      call. = FALSE
    )
  }

  invisible(TRUE)
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

# -----------------------------------------------------------------------------
# Generalized linear models ---------------------------------------------------
# -----------------------------------------------------------------------------

#' Muffle Only the Base-GLM Non-Convergence Warning
#'
#' Evaluates a base-GLM fit while suppressing the specific backend
#' non-convergence warning that MFP immediately replaces with a targeted
#' fail-fast error in `validate_mfp_fit_result()`. Every other `glm()` or
#' `glm.fit()` warning (for example fitted probabilities near 0 or 1) is left
#' visible to the user. Suppression is condition-only; fitting itself and
#' the reported convergence status are unaffected.
#'
#' @param expr Unevaluated expression that fits a base GLM.
#'
#' @return The value of `expr`.
#'
#' @keywords internal
#' @noRd
mfp2_muffle_glm_nonconvergence_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (identical(conditionMessage(w), "glm.fit: algorithm did not converge")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Function that fits generalized linear models
#'
#' @param x a matrix of predictors with nobs observations.
#' @param y Response variable. For GLMs, this may be a numeric vector, a factor
#' response accepted by [stats::glm()], or for binomial models a two-column
#' matrix of grouped integer counts `cbind(successes, failures)`.
#' @param family a family function e.g. `stats::gaussian()`.
#' @param weights a numeric vector of length nobs of 'prior weights' to be used
#' in the fitting process. see [stats::glm()] for details.
#' @param offset a numeric vector of length nobs of of a priori known component
#' to be included in the linear predictor during fitting.
#' @param control A normalized GLM control list. For the fastglm backend,
#' `epsilon` is mapped to `tol` and `maxit` is passed through;
#' `trace = TRUE` is not supported by fastglm.
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
#' @param reserved_names Optional character vector of names that package-created
#' response/offset columns must avoid in the formula-based final refit.
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
                    control,
                    fast = TRUE,
                    calculate_fit_statistics = FALSE,
                    calculate_gaussian_deviance = FALSE,
                    keep_fit = !fast,
                    has_offset = FALSE,
                    x_has_intercept = FALSE,
                    fitter = "base",
                    family_string,
                    keep_fitted_values = FALSE,
                    reserved_names = character()) {
  nobs <- NROW(y)

  fit_flags <- attr(family, "mfp2_fit_flags", exact = TRUE)
  if (is.null(fit_flags)) {
    # Defensive support for direct internal helper calls. Public mfp2()/mfpi()
    # paths attach these invariant flags in prepare_family_for_fit().
    fit_flags <- list(
      is_gaussian = identical(family_string, "gaussian"),
      is_negbin = identical(family_string, "negbin"),
      estimates_dispersion = mfp2_glm_estimates_dispersion(
        family = family,
        family_string = family_string
      )
    )
  }
  is_gaussian <- isTRUE(fit_flags$is_gaussian)
  is_negbin <- isTRUE(fit_flags$is_negbin)
  estimates_dispersion <- isTRUE(fit_flags$estimates_dispersion)

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

    fit <- if (is_negbin) {
      fit_glm_fastglm_nb(xx, y, weights, offset, control)
    } else if (identical(fitter, "fastglm")) {
      fit_glm_fastglm(xx, y, family, weights, offset, control)
    } else {
      mfp2_muffle_glm_nonconvergence_warning(
        stats::glm.fit(
          x = xx,
          y = y,
          family = family,
          weights = weights,
          offset = offset,
          control = control
        )
      )
    }

  } else if (is_negbin) {
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
    fit <- fit_glm_fastglm_nb(xx, y, weights, offset, control)

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

    # Allocate package-created columns only for this single formula-based final
    # refit. Predictor names are never renamed, and the repeated candidate-fit
    # path above remains unchanged. Collision checks therefore cost only a few
    # string comparisons once per fitted model.
    # Reserve both final transformed-column names and the original user-facing
    # predictor names. A selected continuous predictor named, for example,
    # `offset_` is transformed to `offset_.1`; without the original name in this
    # set the temporary offset column could still reuse `offset_` and later
    # collide with prediction data. This check runs only once for the final fit.
    used_names <- unique(c(names(data), as.character(reserved_names)))
    internal_names <- list(response = NULL, offset = NULL, strata = NULL)

    if (is.matrix(y)) {
      response_cols <- character(2L)
      response_cols[1L] <- mfp2_internal_name("successes", used_names, preferred = "..mfp2_successes")
      used_names <- c(used_names, response_cols[1L])
      response_cols[2L] <- mfp2_internal_name("failures", used_names, preferred = "..mfp2_failures")
      used_names <- c(used_names, response_cols[2L])

      data[[response_cols[1L]]] <- y[, 1L]
      data[[response_cols[2L]]] <- y[, 2L]
      internal_names$response <- response_cols

      lhs <- sprintf(
        "cbind(%s, %s)",
        response_cols[1L],
        response_cols[2L]
      )

    } else {
      response_col <- mfp2_internal_name("response", used_names, preferred = "..mfp2_y")
      used_names <- c(used_names, response_col)

      data[[response_col]] <- y
      internal_names$response <- response_col
      lhs <- response_col
    }

    if (isTRUE(has_offset)) {
      offset_col <- mfp2_internal_name("offset", used_names, preferred = "offset_")
      used_names <- c(used_names, offset_col)
      data[[offset_col]] <- offset
      internal_names$offset <- offset_col

      offset_term <- paste0("offset(", offset_col, ")")
      rhs <- if (identical(rhs, "1")) {
        offset_term
      } else {
        paste(rhs, "+", offset_term)
      }
    }

    formula <- stats::as.formula(paste(lhs, "~", rhs))

    fit <- mfp2_muffle_glm_nonconvergence_warning(
      stats::glm(
        formula = formula,
        data = data,
        family = family,
        weights = weights,
        control = control,
        x = TRUE,
        y = TRUE
      )
    )

    # Prediction reuses these exact names so the stored formula never depends
    # on a fixed helper name such as `offset_`.
    fit$mfp2_internal_names <- internal_names
  }

  # Gaussian, Gamma, and inverse-Gaussian dispersion and negative-binomial
  # theta are estimated nuisance parameters in addition to regression rank.
  # fastglm_nb() computes its glm.fit-style null deviance from a constant
  # weighted mean. MASS::glm.nb() performs one additional intercept-only IRLS
  # fit when an offset is present, holding the full model's theta fixed. Do the
  # same only when fit statistics have explicitly been requested. Candidate
  # fits leave calculate_fit_statistics = FALSE, so this diagnostic refit stays
  # outside the repeated model-selection hot path.
  if (is_negbin && isTRUE(calculate_fit_statistics) &&
      !is.null(offset) && any(offset != 0)) {
    fit$null.deviance <- mfp2_negbin_offset_null_deviance(
      y = y,
      weights = weights,
      offset = offset,
      family = fit$family,
      control = control
    )
  }

  df <- fit$rank + as.integer(estimates_dispersion || is_negbin)

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

  # Validate immediately after fitting, before this model can enter an MFP
  # comparison. Invalid candidates are neither removed from the search space
  # nor retried with different fitting controls.
  convergence_flag <- if (is.null(fit$converged)) NULL else fit$converged
  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    converged = convergence_flag,
    gaussian_deviance = gaussian_deviance,
    require_gaussian_deviance = is_gaussian && isTRUE(calculate_gaussian_deviance),
    family_string = family_string,
    fast = fast
  )

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

#' Normalize GLM fitting controls
#'
#' Converts `NULL` or a user-supplied control list into a validated
#' [stats::glm.control()] list. Normalizing once at the fitting boundary gives
#' all GLM backends the same `epsilon`, `maxit`, and `trace` semantics.
#'
#' @param control `NULL` or a list accepted by [stats::glm.control()].
#'
#' @return A validated list returned by [stats::glm.control()].
#' @keywords internal
#' @noRd
normalize_glm_control <- function(control = NULL) {
  if (is.null(control)) {
    return(stats::glm.control())
  }

  if (!is.list(control)) {
    stop(
      "For GLMs, `control` must be `NULL` or a list accepted by ",
      "`stats::glm.control()`.",
      call. = FALSE
    )
  }

  tryCatch(
    do.call(stats::glm.control, control),
    error = function(e) {
      stop(
        "Invalid GLM `control`: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

#' Resolve the effective GLM fitting backend once per top-level fit
#'
#' Ordinary GLMs fall back to the base fitter with one warning when fastglm is
#' unavailable or does not expose the required matrix fitter. Negative-binomial
#' models have no base fallback and therefore fail before candidate fitting.
#' Cox models do not use a GLM fitter and always resolve to "base".
#'
#' @param fitter Requested fitter, "base" or "fastglm". Negative-binomial
#'   models always resolve silently to "fastglm".
#' @param family_string Normalized family name.
#'
#' @return The effective fitter, either "base" or "fastglm".
#' @keywords internal
#' @noRd
resolve_fitter <- function(fitter, family_string) {
  fitter <- match.arg(fitter, c("base", "fastglm"))

  if (!mfp2_family_is_glm(family_string)) {
    return("base")
  }

  # Negative-binomial fitting has only one valid backend. Resolve it here,
  # before package/export checks and before any response or candidate-model
  # work, so callers do not have to repeat `fitter = "fastglm"`. This override
  # is intentionally silent: it is a family requirement rather than a fallback.
  if (identical(family_string, "negbin")) fitter <- "fastglm"

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

#' Fit a GLM with fastglm Without Candidate-level Refitting
#'
#' Package availability is resolved once by `resolve_fitter()` before entering
#' the candidate-search path. Once fastglm is selected, an error or unusable
#' result is fatal for that MFP fit: the candidate is not silently retried with
#' `stats::glm.fit()`, because changing backend after seeing a failed candidate
#' would make the search path backend-dependent.
#'
#' @keywords internal
#' @noRd
fit_glm_fastglm <- function(x, y, family, weights, offset, control) {
  fit <- tryCatch(
    fastglm::fastglm(
      x = x,
      y = y,
      family = family,
      weights = weights,
      offset = offset,
      method = 0L,
      # fastglm uses `tol` for glm.control()'s convergence `epsilon`.
      tol = control$epsilon,
      maxit = control$maxit
    ),
    error = function(e) {
      stop(
        "fastglm could not fit an MFP candidate model: ",
        conditionMessage(e),
        ". The candidate was not refit with another backend.",
        call. = FALSE
      )
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
    stop(
      "fastglm returned an incomplete or non-finite MFP candidate fit; ",
      "the candidate was not discarded or refit with another backend.",
      call. = FALSE
    )
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
fit_glm_fastglm_nb <- function(x, y, weights, offset, control) {
  fit <- tryCatch(
    fastglm::fastglm_nb(
      x = x,
      y = y,
      weights = weights,
      offset = offset,
      method = 0L,
      # Apply glm.control() to the inner IRLS fit. fastglm_nb-specific outer
      # optimization controls keep their documented defaults.
      tol = control$epsilon,
      maxit = control$maxit
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


#' Offset-Aware Null Deviance for Negative-Binomial Models
#'
#' Reproduces the offset-aware null deviance that [MASS::glm.nb()] would
#' report, without re-estimating the dispersion parameter theta. Only the
#' intercept-only model deviance is required, so theta is reused from the
#' full fit. Called only when the reporting statistics of the full-reference
#' or retained-final fit are requested; never invoked for ordinary MFP
#' candidate models.
#'
#' @param y Numeric response vector of finite nonnegative integer counts.
#' @param weights Numeric vector of strictly positive observation weights.
#' @param offset Numeric offset vector aligned with `y`.
#' @param family Prepared negative-binomial family object.
#' @param control Fitting control list used for the null-model iteration.
#'
#' @return Numeric scalar with the null deviance for an intercept-only
#'   negative-binomial model with the supplied offset.
#'
#' @keywords internal
#' @noRd
mfp2_negbin_offset_null_deviance <- function(y, weights, offset, family,
                                             control) {
  nobs <- NROW(y)
  null_x <- matrix(
    1,
    nrow = nobs,
    ncol = 1L,
    dimnames = list(names(y), "(Intercept)")
  )

  null_fit <- mfp2_muffle_glm_nonconvergence_warning(
    stats::glm.fit(
      x = null_x,
      y = y,
      weights = weights,
      offset = offset,
      family = family,
      control = control,
      intercept = TRUE
    )
  )

  if (!isTRUE(null_fit$converged) ||
      !is_finite_numeric_scalar(null_fit$deviance)) {
    stop(
      "The offset-aware negative-binomial null model did not converge; ",
      "adjust `control` or inspect the data.",
      call. = FALSE
    )
  }

  unname(null_fit$deviance)
}


# -----------------------------------------------------------------------------
# Shared control dispatch -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Validate a Control List Against the Selected Fitter Backend
#'
#' Cross-checks the user-supplied `control` list against the chosen `fitter`
#' backend once, at the public fitting boundary, so that repeated candidate
#' fits do not have to revalidate it. Family-specific rules are enforced
#' here: for example, `trace = TRUE` is rejected under
#' `fitter = "fastglm"`, and the negative-binomial family requires the
#' `"fastglm"` backend.
#'
#' @param control List of control settings, typically produced by
#'   [stats::glm.control()], [survival::coxph.control()], or
#'   [survival::survreg.control()]. May be `NULL`.
#' @param family_string Canonical family name.
#' @param fitter Fitting backend for repeated GLM candidate models. Either
#'   `"base"` or `"fastglm"`.
#'
#' @return Invisibly returns `TRUE` on success. Invalid combinations raise
#'   an error identifying the offending setting.
#'
#' @keywords internal
#' @noRd
validate_fit_control_for_fitter <- function(control, family_string,
                                            fitter = "base") {
  uses_fastglm <- identical(family_string, "negbin") ||
    (mfp2_family_is_glm(family_string) && identical(fitter, "fastglm"))

  if (uses_fastglm && isTRUE(control$trace)) {
    if (identical(family_string, "negbin")) {
      stop(
        "`control$trace = TRUE` is not supported for `family = \"negbin\"`; ",
        "the required fastglm backend does not provide GLM iteration tracing.",
        call. = FALSE
      )
    }

    stop(
      "`control$trace = TRUE` is not supported with `fitter = \"fastglm\"`; ",
      "use `fitter = \"base\"` to enable GLM iteration tracing.",
      call. = FALSE
    )
  }

  invisible(control)
}


#' Normalise Family-Specific Fitting Controls
#'
#' Resolves the user-supplied `control` list to the concrete object expected
#' by the low-level fitter for a given family. This normalisation happens
#' once at the public fitting boundary; every subsequent MFP candidate fit
#' reuses the resolved object unchanged.
#'
#' @param control User-supplied control list, or `NULL` for family defaults.
#' @param family_string Canonical family name.
#' @param fitter Backend selector for GLM families, one of `"base"` or
#'   `"fastglm"`. Ignored for non-GLM families.
#'
#' @return A control object appropriate for the selected family's fitter:
#'   a `glm.control()` list for GLM families, a `coxph.control()` list for
#'   Cox and Fine--Gray, and a `survreg.control()` list for parametric
#'   survival.
#'
#' @keywords internal
#' @noRd
normalize_fit_control <- function(control = NULL, family_string,
                                  fitter = "base") {
  normalized <- if (family_string %in% c("cox", "finegray")) {
    normalize_cox_control(control)
  } else if (identical(family_string, "survreg")) {
    normalize_survreg_control(control)
  } else if (identical(family_string, "multinomial")) {
    normalize_multinomial_control(control)
  } else if (identical(family_string, "ordinal")) {
    normalize_ordinal_control(control)
  } else {
    normalize_glm_control(control)
  }

  validate_fit_control_for_fitter(
    control = normalized,
    family_string = family_string,
    fitter = fitter
  )
  normalized
}


# -----------------------------------------------------------------------------
# Ordinal models --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Normalise a Control List for the Ordinal Fitter
#'
#' Resolves the user-supplied control list to explicit values for the
#' [rms::orm.fit()]-based ordinal fitter. Using explicit defaults here means
#' that MFP candidate-model behaviour does not silently drift with changes to
#' the defaults of successive `rms` releases. The GLM-style `epsilon`
#' argument is accepted as an alias for the `rms`-style `eps`.
#'
#' @param control User-supplied control list, or `NULL` to use MFP defaults.
#'
#' @return A named list carrying `maxit`, `eps`, and any additional entries
#'   accepted by [rms::orm.fit()], with unspecified fields filled from the
#'   MFP defaults.
#'
#' @keywords internal
#' @noRd
normalize_ordinal_control <- function(control = NULL) {
  defaults <- list(
    maxit = 30L,
    eps   = 0.005,
    tol   = 1.0e-7,
    trace = FALSE
  )
  if (is.null(control)) return(defaults)
  if (!is.list(control)) {
    stop(
      "For ordinal models, `control` must be `NULL` or a list of ordinal ",
      "control values (`maxit`, `eps`, `tol`, `trace`).",
      call. = FALSE
    )
  }

  # Use exact [[ extraction: `$` partial-matches, so control$eps would spuriously
  # match a supplied `epsilon`.
  if (!is.null(control[["epsilon"]]) && is.null(control[["eps"]])) {
    control[["eps"]] <- control[["epsilon"]]
  }
  control[["epsilon"]] <- NULL

  allowed <- names(defaults)
  unknown <- setdiff(names(control), allowed)
  if (length(unknown) > 0L) {
    stop(
      "Unknown ordinal control field(s): ", paste(unknown, collapse = ", "),
      ". Allowed fields are: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  out <- utils::modifyList(defaults, control)

  if (!is.numeric(out$maxit) || length(out$maxit) != 1L || anyNA(out$maxit) ||
      !is.finite(out$maxit) || out$maxit < 1 ||
      out$maxit > .Machine$integer.max || out$maxit != floor(out$maxit)) {
    stop("Ordinal `control$maxit` must be a single positive integer.",
         call. = FALSE)
  }
  out$maxit <- as.integer(out$maxit)
  for (nm in c("eps", "tol")) {
    if (!is.numeric(out[[nm]]) || length(out[[nm]]) != 1L || anyNA(out[[nm]]) ||
        !is.finite(out[[nm]]) || out[[nm]] <= 0) {
      stop("Ordinal `control$", nm, "` must be a single positive number.",
           call. = FALSE)
    }
  }
  if (!is.logical(out$trace) || length(out$trace) != 1L || anyNA(out$trace)) {
    stop("Ordinal `control$trace` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  out
}


#' Convert Ordinal-Fit Non-Convergence Warnings into Fail-Fast Errors
#'
#' Wraps an expression that fits an ordinal model with [rms::orm()] or
#' [rms::orm.fit()], catching the non-convergence warnings those functions
#' emit and rethrowing them as a targeted error. This mirrors the survreg
#' convergence guard so that candidate and final-fit paths surface fitter
#' failures immediately instead of silently producing suspect fits.
#'
#' @param expr Unevaluated expression that fits the ordinal model.
#' @param fast Logical. `TRUE` (default) indicates a candidate-model fit,
#'   which is reported as such in the error; `FALSE` indicates the retained
#'   final fit.
#'
#' @return The value of `expr` when it converges.
#'
#' @keywords internal
#' @noRd
mfp2_with_ordinal_convergence_guard <- function(expr, fast = TRUE) {
  stage <- if (isTRUE(fast)) "candidate" else "final"
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("did not converge|singular|Non-convergence", msg, ignore.case = TRUE)) {
        stop(
          "The ", stage, " ordinal model did not converge; adjust `control` ",
          "(for example raise `maxit`) or inspect the data.",
          call. = FALSE
        )
      }
    }
  )
}


#' Fit a Proportional-Odds Ordinal Model
#'
#' Fits an ordinal regression model with a cumulative-link (proportional-odds)
#' representation. Candidate MFP fits call [rms::orm.fit()] directly on the
#' model matrix and integer-coded response for speed; the retained final
#' model is returned as a native [rms::orm()] object so that the standard
#' `rms` methods (`coef()`, `vcov()`, `logLik()`, `AIC()`, `predict()`,
#' `summary()`) work without translation. The integer coding of the response
#' is prepared once in `prepare_ordinal_family()`.
#'
#' @param x Numeric model matrix.
#' @param family Prepared ordinal family object (see
#'   `prepare_ordinal_family()`).
#' @param control Normalised ordinal-fit control list (see
#'   `normalize_ordinal_control()`).
#' @param fast Logical. `TRUE` fits the fast candidate path through
#'   [rms::orm.fit()]; `FALSE` refits with [rms::orm()] to produce the
#'   retained final model.
#' @param calculate_fit_statistics Logical. If `TRUE`, compute log-likelihood
#'   and information criteria for reporting.
#' @param keep_fit Logical. If `TRUE`, retain the fitted-model object in the
#'   result. Defaults to `!fast`.
#' @param keep_fitted_values Logical. If `TRUE`, retain fitted values in the
#'   result.
#' @param ... Additional arguments forwarded to the underlying fitter.
#'
#' @return A list of fit results with components used by the MFP engine and,
#'   when `keep_fit` is `TRUE`, the fitted-model object.
#'
#' @keywords internal
#' @noRd
fit_ordinal <- function(x,
                        family,
                        control,
                        fast = TRUE,
                        calculate_fit_statistics = FALSE,
                        keep_fit = !fast,
                        x_has_intercept = FALSE,
                        reserved_names = character()) {
  prepared <- family$prepared
  ycodes <- prepared$y
  link <- prepared$link
  k <- prepared$n_classes
  n_int <- k - 1L
  nobs <- prepared$n_original
  has_offset <- isTRUE(prepared$has_offset)
  offset <- prepared$offset

  # rms::orm.fit() adds the k-1 intercepts itself, so the design must not carry
  # an intercept column. Strip a leading "(Intercept)" when the caller supplies
  # an intercept-augmented matrix (the convention for has-intercept families).
  if (isTRUE(x_has_intercept)) {
    if (is.null(x) || NCOL(x) == 0L || is.null(colnames(x)) ||
        !identical(colnames(x)[[1L]], "(Intercept)")) {
      stop(
        "Internal error: x_has_intercept = TRUE requires a leading '(Intercept)' column.",
        call. = FALSE
      )
    }
    x <- x[, -1L, drop = FALSE]
  } else if (!is.null(x) && NCOL(x) > 0L && !is.null(colnames(x)) &&
             identical(colnames(x)[[1L]], "(Intercept)")) {
    stop(
      "Internal error: ordinal predictor matrix contains an undeclared intercept.",
      call. = FALSE
    )
  }

  has_predictors <- !is.null(x) && NCOL(x) > 0L
  xx <- if (has_predictors) as.matrix(x) else NULL
  raw_args <- list(
    x = xx,
    y = ycodes,
    family = link,
    maxit = control$maxit,
    eps = control$eps,
    tol = control$tol,
    trace = isTRUE(control$trace)
  )
  if (has_offset) raw_args$offset <- offset

  # orm.fit() natively supports an intercept-only ordinal model when x is NULL.
  # Use the same backend for every candidate so thresholds, covariance metadata,
  # offsets, likelihoods, and convergence handling have one implementation.
  raw <- mfp2_with_ordinal_convergence_guard(
    do.call(rms::orm.fit, raw_args),
    fast = fast
  )
  if (isTRUE(raw$fail)) {
    validate_mfp_fit_result(
      logl = NA_real_, df = n_int, converged = FALSE,
      family_string = "ordinal", fast = fast
    )
  }
  dev <- raw$deviance
  model_dev <- unname(dev[length(dev)])
  null_dev <- if (length(dev) >= 2L) unname(dev[1L]) else NA_real_
  logl <- -model_dev / 2
  all_coef <- raw$coefficients
  betas <- if (length(all_coef) > n_int) {
    all_coef[seq.int(n_int + 1L, length(all_coef))]
  } else {
    numeric(0L)
  }
  if (length(betas) > 0L) names(betas) <- colnames(xx)

  converged <- TRUE
  p <- length(betas)
  rank <- p
  df <- n_int + rank

  result <- list(
    logl = unname(logl),
    coefficients = betas,
    rank = rank,
    df = df,
    weights = NULL,
    residuals = NULL
  )
  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    converged = converged,
    family_string = "ordinal",
    fast = fast
  )

  if (isTRUE(keep_fit)) {
    # Retain a native rms::orm object so downstream accessors and prediction use
    # the standard rms methods. Fit once, on the transformed design.
    if (has_predictors) {
      dfm <- as.data.frame(xx, check.names = FALSE)
    } else {
      dfm <- data.frame(row.names = seq_len(nobs))
    }
    used_names <- unique(c(names(dfm), as.character(reserved_names)))
    response_col <- mfp2_internal_name(
      "response", used_names, preferred = "..mfp2_ordinal_y"
    )
    used_names <- c(used_names, response_col)
    dfm[[response_col]] <- ycodes
    rhs <- if (has_predictors) {
      paste(sprintf("`%s`", colnames(xx)), collapse = " + ")
    } else {
      "1"
    }
    if (has_offset) {
      offset_col <- mfp2_internal_name("offset", used_names, preferred = "offset_")
      dfm[[offset_col]] <- offset
      rhs <- if (identical(rhs, "1")) {
        paste0("offset(", offset_col, ")")
      } else {
        paste(rhs, "+", paste0("offset(", offset_col, ")"))
      }
    }
    formula <- stats::as.formula(paste(response_col, "~", rhs))
    fitobj <- mfp2_with_ordinal_convergence_guard(
      rms::orm(formula, data = dfm, family = link, x = TRUE, y = TRUE,
               maxit = control$maxit, eps = control$eps, tol = control$tol,
               trace = isTRUE(control$trace)),
      fast = FALSE
    )
    if (isTRUE(fitobj$fail)) {
      stop(
        "The final ordinal model did not converge; adjust `control` (for ",
        "example raise `maxit`) or inspect the data.",
        call. = FALSE
      )
    }
    fitobj$mfp2_ordinal_levels <- prepared$levels
    fitobj$mfp2_ordinal_link <- link
    # The public family component is intentionally stripped before returning
    # the mfp2 object. Keep the prepared family separately for prediction and
    # nonlinear-term reduced refits.
    fitobj$mfp2_family <- family
    retained_coef <- fitobj$coefficients
    if (length(retained_coef) < n_int) {
      stop("Internal error: retained ordinal thresholds are incomplete.",
           call. = FALSE)
    }
    fitobj$mfp2_ordinal_intercepts <- retained_coef[seq_len(n_int)]
    result$fit <- fitobj
  }

  if (isTRUE(calculate_fit_statistics)) {
    result$null_logl <- if (is.finite(null_dev)) -null_dev / 2 else NA_real_
    result$null_deviance <- null_dev
    result$model_deviance <- model_dev
  }
  result
}


# -----------------------------------------------------------------------------
# Multinomial models ----------------------------------------------------------
# -----------------------------------------------------------------------------

# Normalize a control list for the nnet-based multinomial fitter.
#
#' Normalise a Control List for the Multinomial Fitter
#'
#' Returns a clean, `nnet`-only control list for
#' [nnet::multinom()]/[nnet::nnet.default()], whose native control knobs
#' differ from those of [stats::glm.control()]. Routing a multinomial fit
#' through `glm.control()` would impose `maxit = 25`, far below `nnet`'s own
#' default of 100, which is easily too tight for FP-transformed candidate
#' designs. The default `maxit` used here is raised to 300 to give BFGS
#' iterations enough head-room for MFP candidate models.
#'
#' @param control User-supplied control list, or `NULL` to use MFP defaults.
#'
#' @return A named list containing the `nnet`-recognised control fields
#'   (`maxit`, `abstol`, `reltol`, `trace`, `MaxNWts`, and related values),
#'   with unspecified fields taking MFP defaults.
#'
#' @keywords internal
#' @noRd
normalize_multinomial_control <- function(control = NULL) {
  defaults <- list(
    maxit   = 300L,      # nnet default is 100; raised for FP candidate fits
    reltol  = 1.0e-8,    # nnet default
    abstol  = 1.0e-4,    # nnet default
    trace   = FALSE,     # nnet default is TRUE; silenced for MFP candidate fits
    MaxNWts = 1000L      # nnet default
  )
  if (is.null(control)) return(defaults)
  if (!is.list(control)) {
    stop(
      "For multinomial models, `control` must be `NULL` or a list of nnet ",
      "control values (`maxit`, `reltol`, `abstol`, `trace`, `MaxNWts`).",
      call. = FALSE
    )
  }

  # Accept the glm-style `epsilon` as an alias for nnet's `reltol` so that
  # control lists shared with GLM families keep working.
  if (!is.null(control$epsilon) && is.null(control$reltol)) {
    control$reltol <- control$epsilon
  }
  control$epsilon <- NULL

  allowed <- names(defaults)
  unknown <- setdiff(names(control), allowed)
  if (length(unknown) > 0L) {
    stop(
      "Unknown multinomial control field(s): ", paste(unknown, collapse = ", "),
      ". Allowed fields are: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  out <- utils::modifyList(defaults, control)

  if (!is.numeric(out$maxit) || length(out$maxit) != 1L || anyNA(out$maxit) ||
      !is.finite(out$maxit) || out$maxit < 1 ||
      out$maxit > .Machine$integer.max || out$maxit != floor(out$maxit)) {
    stop("Multinomial `control$maxit` must be a single positive integer.",
         call. = FALSE)
  }
  out$maxit <- as.integer(out$maxit)
  for (nm in c("reltol", "abstol")) {
    if (!is.numeric(out[[nm]]) || length(out[[nm]]) != 1L || anyNA(out[[nm]]) ||
        !is.finite(out[[nm]]) || out[[nm]] <= 0) {
      stop("Multinomial `control$", nm, "` must be a single positive number.",
           call. = FALSE)
    }
  }
  if (!is.logical(out$trace) || length(out$trace) != 1L || anyNA(out$trace)) {
    stop("Multinomial `control$trace` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.numeric(out$MaxNWts) || length(out$MaxNWts) != 1L ||
      anyNA(out$MaxNWts) || !is.finite(out$MaxNWts) || out$MaxNWts < 1 ||
      out$MaxNWts > .Machine$integer.max || out$MaxNWts != floor(out$MaxNWts)) {
    stop("Multinomial `control$MaxNWts` must be a single positive integer.",
         call. = FALSE)
  }
  out$MaxNWts <- as.integer(out$MaxNWts)
  out
}


#' Flatten a Multinomial Coefficient Matrix to a Named Vector
#'
#' Converts a coefficient matrix with rows indexed by non-reference outcome
#' and columns indexed by design term into the flat named vector used by
#' [nnet::multinom()] and by MFP's internal candidate-fit reporting. The
#' order matches `coef.multinom()`[nnet::coef.multinom]: all design
#' coefficients for one non-reference outcome, followed by the next outcome.
#'
#' @param coefficient_matrix Numeric matrix of coefficients with dimensions
#'   `Q` (non-reference outcomes) by `p` (design terms).
#'
#' @return Numeric vector of length `Q * p` with element names of the form
#'   `"<outcome>:<term>"`.
#'
#' @keywords internal
#' @noRd
mfp2_multinomial_flatten <- function(coefficient_matrix) {
  values <- as.vector(t(coefficient_matrix))
  names(values) <- paste0(
    rep(rownames(coefficient_matrix), each = ncol(coefficient_matrix)),
    "::",
    rep(colnames(coefficient_matrix), times = nrow(coefficient_matrix))
  )
  values
}


#' Fit a Multinomial Model via `nnet::nnet.default()`
#'
#' Runs a single native `nnet` fit for a supplied design matrix `xx`. This
#' is the fast candidate path used by the MFP multinomial engine and is also
#' used to construct the retained final model when a callable
#' [nnet::multinom()] wrapper is not needed. Response, weights, and class
#' offsets are read from the prepared family object.
#'
#' @param xx Numeric design matrix already carrying any intercept column.
#' @param family Prepared multinomial family object (see
#'   `prepare_multinomial_family()`).
#' @param control Normalised multinomial control list (see
#'   `normalize_multinomial_control()`).
#' @param Hess Logical. If `TRUE`, request the observed-information Hessian
#'   from `nnet.default()`.
#'
#' @return A list with fitted coefficients, log-likelihood, convergence
#'   status, and, when requested, the Hessian.
#'
#' @keywords internal
#' @noRd
mfp2_fit_multinom_native <- function(xx, family, control, Hess = FALSE) {
  prepared <- family$prepared
  offset_matrix <- prepared$offset_matrix
  include_offset <- isTRUE(prepared$has_offset)
  predictor_names <- colnames(xx)
  if (is.null(predictor_names) || !identical(predictor_names[[1L]], "(Intercept)")) {
    stop("Internal error: multinomial design must start with an intercept.", call. = FALSE)
  }

  x_no_intercept <- xx[, -1L, drop = FALSE]
  data <- if (ncol(x_no_intercept) == 0L) {
    data.frame(row.names = seq_len(nrow(xx)))
  } else {
    as.data.frame(x_no_intercept, check.names = FALSE)
  }
  rhs <- if (ncol(x_no_intercept) == 0L) {
    "1"
  } else {
    paste(sprintf("`%s`", colnames(x_no_intercept)), collapse = " + ")
  }
  fit_env <- new.env(parent = parent.frame())
  fit_env$.mfp2_multinomial_response <- prepared$y
  if (isTRUE(include_offset)) {
    fit_env$.mfp2_multinomial_offset <- offset_matrix
    rhs <- paste(rhs, "+ offset(.mfp2_multinomial_offset)")
  }
  formula <- stats::as.formula(
    paste(".mfp2_multinomial_response ~", rhs),
    env = fit_env
  )

  max_weights <- max(
    if (is.null(control$MaxNWts)) 1000L else control$MaxNWts,
    as.integer((ncol(xx) + prepared$n_classes + 1L) * prepared$n_classes)
  )
  # nnet::multinom() appends fixed class-offset columns to its internal design
  # before calling multinomHess(). multinomHess() cannot align those fixed
  # columns with the freely estimated coefficient vector. Fit natively, but
  # calculate the information matrix below from the estimable design when a
  # matrix-valued offset is present.
  fit <- nnet::multinom(
    formula = formula,
    data = data,
    weights = prepared$case_weights,
    Hess = isTRUE(Hess) && !isTRUE(include_offset),
    model = TRUE,
    trace = isTRUE(control$trace),
    maxit = control$maxit,
    reltol = control$reltol,
    abstol = control$abstol,
    MaxNWts = max_weights
  )

  if (isTRUE(Hess) && isTRUE(include_offset)) {
    fit$Hessian <- mfp2_multinomial_information(
      x = xx,
      probabilities = fit$fitted.values,
      weights = fit$weights,
      outcomes = prepared$nonreference,
      class_levels = prepared$levels
    )
  }

  fit
}


# Expected information for the freely estimated baseline-category logit
#' Observed Information Matrix for a Multinomial Fit
#'
#' Assembles the weighted observed-information matrix for a baseline-category
#' multinomial model, evaluated at the current coefficient estimates. Fixed
#' class offsets shift probabilities but have no derivative with respect to
#' the regression coefficients, so they must not appear as design columns.
#' The parameter order matches [nnet::coef.multinom()] and
#' [nnet::vcov.multinom()]: all design coefficients for one non-reference
#' outcome, followed by the next outcome.
#'
#' @param x Numeric design matrix.
#' @param probabilities Numeric matrix of predicted class probabilities with
#'   one row per observation and one column per class.
#' @param weights Numeric vector of case weights.
#' @param outcomes Numeric matrix of observed class outcomes (indicator
#'   representation).
#' @param class_levels Character vector of class labels in fitted order; the
#'   first level is the reference class.
#'
#' @return Numeric observed-information matrix of order `Q p` by `Q p`,
#'   where `Q` is the number of non-reference classes and `p` is the number
#'   of estimable design columns.
#'
#' @keywords internal
#' @noRd
mfp2_multinomial_information <- function(x, probabilities, weights, outcomes,
                                         class_levels) {
  probabilities <- as.matrix(probabilities)
  if (!is.null(colnames(probabilities)) &&
      all(class_levels %in% colnames(probabilities))) {
    probabilities <- probabilities[, class_levels, drop = FALSE]
  }
  if (ncol(probabilities) != length(class_levels)) {
    stop(
      "Internal error: multinomial fitted probabilities are incomplete.",
      call. = FALSE
    )
  }

  q_logits <- length(outcomes)
  p_columns <- ncol(x)
  parameter_names <- paste0(
    rep(outcomes, each = p_columns),
    ":",
    rep(colnames(x), times = q_logits)
  )
  information <- matrix(
    0,
    nrow = q_logits * p_columns,
    ncol = q_logits * p_columns,
    dimnames = list(parameter_names, parameter_names)
  )
  nonreference_probabilities <- probabilities[, -1L, drop = FALSE]
  weights <- as.numeric(weights)

  for (j in seq_len(q_logits)) {
    rows <- (j - 1L) * p_columns + seq_len(p_columns)
    for (k in seq_len(q_logits)) {
      columns <- (k - 1L) * p_columns + seq_len(p_columns)
      working_weights <- weights * nonreference_probabilities[, j] *
        ((j == k) - nonreference_probabilities[, k])
      information[rows, columns] <- crossprod(
        x,
        x * working_weights
      )
    }
  }

  information
}


#' Build the Shared `nnet` Mask and Offset Weights for Multinomial Candidates
#'
#' Constructs the `nnet` parameter mask and fixed offset weights shared by
#' every multinomial candidate that has the same number of estimable design
#' columns. Building this once per candidate class avoids rebuilding the
#' mask for every candidate fit within the class, which is one of the
#' hottest paths of the multinomial MFP engine.
#'
#' @param p Integer number of estimable design columns.
#' @param c_classes Integer number of response classes.
#' @param has_offset Logical flag indicating whether the model uses a class
#'   offset.
#' @param MaxNWts Integer maximum number of weights allowed by `nnet`.
#'
#' @return A list with the parameter mask, offset-weight vector, and helper
#'   indices consumed by [nnet::nnet.default()] during candidate fitting.
#'
#' @keywords internal
#' @noRd
mfp2_multinomial_optimizer_structure <- function(p, c_classes, has_offset,
                                                  MaxNWts = 1000L) {
  if (isTRUE(has_offset)) {
    mask <- c(
      rep(FALSE, p + 1L + c_classes),
      rep(
        c(FALSE, rep(TRUE, p), rep(FALSE, c_classes)),
        c_classes - 1L
      )
    )
    fixed_offset_weights <- as.vector(
      rbind(matrix(0, p + 1L, c_classes), diag(c_classes))
    )
  } else {
    mask <- c(
      rep(FALSE, p + 1L),
      rep(c(FALSE, rep(TRUE, p)), c_classes - 1L)
    )
    fixed_offset_weights <- NULL
  }

  list(
    p = p,
    n_classes = c_classes,
    has_offset = isTRUE(has_offset),
    mask = mask,
    fixed_offset_weights = fixed_offset_weights,
    MaxNWts = max(MaxNWts, length(mask))
  )
}


#' Matrix-Level Multinomial Candidate Fitter
#'
#' Fast candidate fitter used by the MFP multinomial engine. The response
#' normalisation, effective case weights, and class offsets are prepared
#' once in `prepare_multinomial_family()` and then shared across every
#' candidate model. Candidate fits call [nnet::nnet.default()] directly;
#' only the retained final model is refitted through the formula-based
#' [nnet::multinom()] wrapper so that user-facing methods keep working.
#'
#' @param x Numeric model matrix.
#' @param family Prepared multinomial family object.
#' @param control Normalised multinomial-fit control list.
#' @param fast Logical. `TRUE` selects the candidate hot path;
#'   `FALSE` refits with [nnet::multinom()] to build the retained final
#'   model.
#' @param calculate_fit_statistics Logical. If `TRUE`, compute
#'   log-likelihood and information criteria for reporting.
#' @param keep_fit Logical. If `TRUE`, retain the fitted-model object in the
#'   result. Defaults to `!fast`.
#' @param keep_fitted_values Logical. If `TRUE`, retain the fitted class
#'   probabilities.
#' @param ... Additional arguments forwarded to the underlying fitter.
#'
#' @return A list of fit results used by the MFP engine; when `keep_fit` is
#'   `TRUE`, it also contains the fitted-model object.
#'
#' @keywords internal
#' @noRd
fit_multinomial <- function(x, family, control,
                            fast = TRUE, calculate_fit_statistics = FALSE,
                            keep_fit = !fast, keep_fitted_values = FALSE,
                            keep_coefficients = TRUE,
                            x_has_intercept = FALSE,
                            optimizer_structure = NULL) {
  prepared <- family$prepared
  nobs <- nrow(prepared$y_matrix)
  xx <- if (isTRUE(x_has_intercept)) {
    x
  } else {
    assemble_design_matrix(list(x), nobs = nobs, intercept = TRUE)
  }
  if (is.null(colnames(xx)) || !identical(colnames(xx)[1L], "(Intercept)")) {
    stop("Internal error: multinomial design must start with '(Intercept)'.", call. = FALSE)
  }

  use_direct <- isTRUE(fast)
  if (use_direct) {
    # nnet.default() does not report the design rank. Candidate transforms can
    # become aliased even when the original input matrix was full rank, so the
    # direct path must calculate it before assigning likelihood-test df.
    design_rank <- qr(xx)$rank
    if (design_rank != ncol(xx)) {
      stop("The multinomial candidate design is rank deficient.", call. = FALSE)
    }

    p <- ncol(xx)
    c_classes <- prepared$n_classes
    has_offset <- isTRUE(prepared$has_offset)
    if (is.null(optimizer_structure)) {
      optimizer_structure <- mfp2_multinomial_optimizer_structure(
        p = p,
        c_classes = c_classes,
        has_offset = has_offset,
        MaxNWts = control$MaxNWts
      )
    }
    structure_is_valid <-
      identical(optimizer_structure$p, p) &&
      identical(optimizer_structure$n_classes, c_classes) &&
      identical(optimizer_structure$has_offset, has_offset)
    if (!isTRUE(structure_is_valid)) {
      stop("Internal error: multinomial optimizer structure does not match the candidate design.",
           call. = FALSE)
    }
    max_weights <- max(
      optimizer_structure$MaxNWts,
      control$MaxNWts,
      length(optimizer_structure$mask)
    )

    x_fit <- if (has_offset) {
      cbind(xx, prepared$offset_matrix)
    } else {
      xx
    }
    fit_args <- list(
      x = x_fit,
      y = prepared$y_fit,
      weights = prepared$effective_weights,
      size = 0L, skip = TRUE, softmax = TRUE, rang = 0, decay = 0,
      mask = optimizer_structure$mask,
      trace = isTRUE(control$trace), maxit = control$maxit,
      reltol = control$reltol, abstol = control$abstol,
      MaxNWts = max_weights
    )
    # Wts is supplied only to encode the fixed identity map from class offsets;
    # it never carries coefficients forward from an earlier candidate.
    if (!is.null(optimizer_structure$fixed_offset_weights)) {
      fit_args$Wts <- optimizer_structure$fixed_offset_weights
    }
    raw_fit <- do.call(nnet::nnet.default, fit_args)
    logl <- -unname(raw_fit$value)
    converged <- length(raw_fit$convergence) == 1L && raw_fit$convergence == 0
    native_fit <- raw_fit
  } else {
    native_fit <- mfp2_fit_multinom_native(
      xx = xx,
      family = family,
      control = control,
      Hess = TRUE
    )
    # multinom() computes qr(X)$rank and its corresponding effective parameter
    # count internally. Reuse those values instead of repeating the QR here.
    design_rank <- native_fit$rank
    if (!is.numeric(design_rank) || length(design_rank) != 1L ||
        !is.finite(design_rank) || design_rank != ncol(xx)) {
      stop("The final multinomial design is rank deficient.", call. = FALSE)
    }
    logl <- -unname(native_fit$value)
    converged <- length(native_fit$convergence) == 1L && native_fit$convergence == 0
  }

  regression_df <- if (use_direct) {
    prepared$n_logits * design_rank
  } else {
    native_fit$edf
  }
  result <- list(
    logl = logl,
    rank = regression_df,
    df = regression_df,
    class_levels = prepared$levels,
    reference_class = prepared$reference,
    n_logits = prepared$n_logits
  )
  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    converged = converged,
    family_string = "multinomial",
    fast = fast
  )

  need_coefficients <- isTRUE(keep_coefficients) || isTRUE(keep_fit)
  coefficient_matrix <- NULL
  if (need_coefficients) {
    if (use_direct) {
      all_coefficients <- matrix(
        raw_fit$wts,
        nrow = prepared$n_classes,
        byrow = TRUE
      )[, 1L + seq_len(ncol(xx)), drop = FALSE]
      coefficient_matrix <- all_coefficients[-1L, , drop = FALSE]
      dimnames(coefficient_matrix) <- list(
        prepared$nonreference, colnames(xx)
      )
    } else {
      coefficient_matrix <- stats::coef(native_fit)
      if (is.null(dim(coefficient_matrix))) {
        coefficient_matrix <- matrix(
          coefficient_matrix,
          nrow = prepared$n_logits,
          dimnames = list(prepared$nonreference, names(coefficient_matrix))
        )
      }
    }
    result$coefficients <- mfp2_multinomial_flatten(coefficient_matrix)
    result$coefficient_matrix <- coefficient_matrix
  }

  if (isTRUE(keep_fitted_values)) {
    result$fitted_values <- native_fit$fitted.values
  }
  if (isTRUE(keep_fit)) {
    # Candidate preparation stores NULL when no offset was supplied so repeated
    # fits do not allocate an n x C zero matrix. A retained multinom object,
    # however, exposes the complete class-offset matrix to prediction and
    # downstream inspection. Materialize that neutral matrix once here, outside
    # the candidate hot path, and preserve its class column names.
    retained_offset_matrix <- prepared$offset_matrix
    if (is.null(retained_offset_matrix)) {
      retained_offset_matrix <- mfp2_multinomial_offset(
        offset = NULL,
        family = family,
        nobs = nobs,
        has_offset = FALSE
      )
    }
    native_fit$coefficients <- coefficient_matrix
    native_fit$mfp2_coefficient_matrix <- coefficient_matrix
    native_fit$mfp2_class_levels <- prepared$levels
    native_fit$mfp2_reference_class <- prepared$reference
    native_fit$mfp2_offset_matrix <- retained_offset_matrix
    native_fit$mfp2_design <- xx
    native_fit$mfp2_family <- family
    native_fit$has_offset <- isTRUE(prepared$has_offset)
    result$fit <- native_fit
  }

  if (isTRUE(calculate_fit_statistics)) {
    null_fit <- mfp2_fit_multinom_native(
      xx = matrix(1, nrow = nobs, ncol = 1L,
                  dimnames = list(NULL, "(Intercept)")),
      family = family,
      control = control,
      Hess = FALSE
    )
    result$null_logl <- -unname(null_fit$value)
    result$null_deviance <- -2 * result$null_logl
    result$model_deviance <- -2 * result$logl
  }
  result
}


# -----------------------------------------------------------------------------
# Survival models -------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Normalise a Control List for the `survreg` Fitter
#'
#' Resolves the user-supplied control list to the concrete object expected
#' by [survival::survreg.fit()] and [survival::survreg()], filling in
#' defaults from [survival::survreg.control()] when values are absent. This
#' is called once per MFP analysis; the resolved list is then reused for
#' every candidate model.
#'
#' @param control User-supplied control list, or `NULL` to use survreg
#'   defaults.
#'
#' @return A `survreg.control()` list carrying the merged fields.
#'
#' @keywords internal
#' @noRd
normalize_survreg_control <- function(control = NULL) {
  if (is.null(control)) return(survival::survreg.control())
  if (!is.list(control)) {
    stop(
      "For `survreg` models, `control` must be `NULL` or a list accepted by `survival::survreg.control()`.",
      call. = FALSE
    )
  }
  tryCatch(
    do.call(survival::survreg.control, control),
    error = function(e) {
      stop("Invalid `survreg` control: ", conditionMessage(e), call. = FALSE)
    }
  )
}


#' Normalise Controls Shared by Cox and Fine--Gray Fitters
#'
#' Resolves the user-supplied control list to a [survival::coxph.control()]
#' object used by both Cox and Fine--Gray subdistribution-hazard candidate
#' fits. Any fields not supplied by the user take defaults from
#' `survival::coxph.control()`.
#'
#' @param control User-supplied control list, or `NULL` to use coxph
#'   defaults.
#'
#' @return A `coxph.control()` list carrying the merged fields.
#'
#' @keywords internal
#' @noRd
normalize_cox_control <- function(control = NULL) {
  if (is.null(control)) return(survival::coxph.control())
  if (!is.list(control)) {
    stop(
      "For Cox and Fine--Gray models, `control` must be `NULL` or a list accepted by `survival::coxph.control()`.",
      call. = FALSE
    )
  }
  tryCatch(
    do.call(survival::coxph.control, control),
    error = function(e) {
      stop("Invalid Cox/Fine--Gray `control`: ", conditionMessage(e), call. = FALSE)
    }
  )
}


#' Convert `survreg` Non-Convergence Warnings into Fail-Fast Errors
#'
#' Wraps an expression that fits a parametric survival model with
#' [survival::survreg()] or [survival::survreg.fit()], catching
#' non-convergence warnings and rethrowing them as targeted errors. This
#' ensures that candidate and final-fit paths surface fitter failures
#' immediately rather than producing suspect fits.
#'
#' @param expr Unevaluated expression that fits the parametric survival
#'   model.
#' @param fast Logical. `TRUE` (default) marks the message as a
#'   candidate-model failure; `FALSE` marks it as a final-model failure.
#'
#' @return The value of `expr` when it converges.
#'
#' @keywords internal
#' @noRd
mfp2_with_survreg_convergence_guard <- function(expr, fast = TRUE) {
  stage <- if (isTRUE(fast)) "candidate" else "final"
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Ran out of iterations and did not converge", conditionMessage(w), fixed = TRUE)) {
        stop(
          "The ", stage, " survreg model did not converge; adjust `control` or inspect the data.",
          call. = FALSE
        )
      }
    }
  )
}


#' Evaluate a Cox Fit While Treating Explicit Non-convergence as Fatal
#'
#' `survival::coxph.fit()` reports iteration exhaustion by warning rather than
#' by a returned `converged` flag. Convert only that specific warning into an
#' error. All other Cox warnings retain survival's normal behavior.
#'
#' @param expr Cox fitting expression to evaluate.
#' @param fast Whether this is a repeated candidate fit or the final refit.
#'
#' @return The evaluated Cox fit.
#' @keywords internal
#' @noRd
mfp2_with_cox_convergence_guard <- function(expr, fast = TRUE) {
  stage <- if (isTRUE(fast)) "candidate" else "final"

  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("Ran out of iterations and did not converge", msg, fixed = TRUE)) {
        stop(
          "The ", stage,
          " Cox model did not converge within the configured iteration limit; ",
          "MFP fitting will not refit or silently discard the model. ",
          "Adjust `control$iter.max` or inspect the data.",
          call. = FALSE
        )
      }
    }
  )
}

#' Fit a Parametric Survival Model
#'
#' Fits an accelerated-failure-time parametric survival model. Candidate MFP
#' fits use the matrix-level hot path through [survival::survreg.fit()] for
#' speed. The retained final model is refitted with [survival::survreg()] so
#' that the standard `survival` methods work through the usual formula
#' interface. The transformed response is prepared once by
#' `prepare_survreg_family()`.
#'
#' @param x Numeric model matrix.
#' @param y Response, typically a [survival::Surv()] object.
#' @param family Prepared survreg family object.
#' @param control Normalised survreg-fit control list.
#' @param fast Logical. `TRUE` selects the candidate hot path;
#'   `FALSE` refits with [survival::survreg()] to build the retained final
#'   model.
#' @param calculate_fit_statistics Logical. If `TRUE`, compute
#'   log-likelihood and information criteria for reporting.
#' @param keep_fit Logical. If `TRUE`, retain the fitted-model object in the
#'   result. Defaults to `!fast`.
#' @param keep_fitted_values Logical. If `TRUE`, retain fitted linear
#'   predictors.
#' @param ... Additional arguments forwarded to the underlying fitter.
#'
#' @return A list of fit results used by the MFP engine; when `keep_fit` is
#'   `TRUE`, also carries the fitted-model object.
#'
#' @keywords internal
#' @noRd
fit_survreg <- function(x,
                        y,
                        family,
                        weights,
                        offset,
                        control,
                        fast = TRUE,
                        calculate_fit_statistics = FALSE,
                        keep_fit = !fast,
                        has_offset = FALSE,
                        x_has_intercept = FALSE,
                        reserved_names = character()) {
  prepared <- family$prepared
  nobs <- prepared$n_original
  has_predictors <- !is.null(x) && NCOL(x) > 0L

  if (isTRUE(x_has_intercept)) {
    if (!has_predictors || is.null(colnames(x)) ||
        !identical(colnames(x)[1L], "(Intercept)")) {
      stop("Internal error: survreg intercept template is invalid.", call. = FALSE)
    }
    xx <- x
  } else {
    xx <- assemble_design_matrix(list(x), nobs = nobs, intercept = TRUE)
  }

  if (fast) {
    raw_fit <- mfp2_with_survreg_convergence_guard(
      survival::survreg.fit(
        x = xx,
        y = prepared$y,
        weights = weights,
        offset = offset,
        init = NULL,
        controlvals = control,
        dist = prepared$dist,
        scale = prepared$scale,
        nstrat = prepared$nstrat,
        strata = prepared$strata,
        parms = prepared$parms,
        # No penalized terms are fitted. Passing NULL satisfies survreg.fit()'s
        # public low-level interface without constructing per-candidate term
        # assignment metadata that this unpenalized routine does not inspect.
        assign = NULL
      ),
      fast = TRUE
    )

    corrected_loglik <- raw_fit$loglik + prepared$logcorrect
    nreg <- NCOL(xx)
    regression_indices <- seq_len(nreg)
    coefficients <- raw_fit$coefficients[regression_indices]
    singular <- diag(raw_fit$var)[regression_indices] == 0
    coefficients[singular] <- NA_real_
    names(coefficients) <- colnames(xx)
    rank <- sum(!singular)
    nuisance_df <- if (prepared$scale == 0) prepared$nstrat else 0L

    result <- list(
      logl = unname(corrected_loglik[length(corrected_loglik)]),
      coefficients = coefficients,
      rank = rank,
      df = rank + nuisance_df,
      weights = NULL,
      residuals = NULL
    )
    validate_mfp_fit_result(
      logl = result$logl,
      df = result$df,
      family_string = "survreg",
      fast = TRUE
    )

    # Candidate selection normally discards the backend fit. Materialize the
    # covariance submatrix and survreg-compatible fields only for the uncommon
    # diagnostic path that explicitly asks to retain it.
    if (isTRUE(keep_fit)) {
      candidate_fit <- raw_fit
      candidate_fit$loglik <- corrected_loglik
      candidate_fit$coefficients <- coefficients
      candidate_fit$var <- raw_fit$var[
        regression_indices, regression_indices, drop = FALSE
      ]
      candidate_fit$scale <- if (prepared$scale == 0) {
        exp(raw_fit$coefficients[nreg + seq_len(prepared$nstrat)])
      } else {
        prepared$scale
      }
      candidate_fit$df <- result$df
      class(candidate_fit) <- "survreg"
      result$fit <- candidate_fit
    }
    if (isTRUE(calculate_fit_statistics)) {
      null_logl <- if (length(corrected_loglik) >= 2L) {
        unname(corrected_loglik[1L])
      } else {
        NA_real_
      }
      result$null_logl <- null_logl
      result$null_deviance <- if (is.finite(null_logl)) -2 * null_logl else NA_real_
      result$model_deviance <- -2 * result$logl
    }
    return(result)
  }

  x_formula <- if (isTRUE(x_has_intercept)) {
    x[, -1L, drop = FALSE]
  } else {
    x
  }
  has_formula_predictors <- !is.null(x_formula) && NCOL(x_formula) > 0L
  if (has_formula_predictors) {
    if (is.null(colnames(x_formula)) || any(colnames(x_formula) == "")) {
      stop("Internal error: x must have non-empty column names.", call. = FALSE)
    }
    data <- data.frame(x_formula, check.names = FALSE)
    rhs <- paste(sprintf("`%s`", colnames(x_formula)), collapse = " + ")
  } else {
    data <- data.frame(row.names = seq_len(nobs))
    rhs <- "1"
  }

  used_names <- unique(c(names(data), as.character(reserved_names)))
  internal_names <- list(response = NULL, offset = NULL, strata = NULL)
  response_col <- mfp2_internal_name("response", used_names, preferred = "..mfp2_survreg_y")
  used_names <- c(used_names, response_col)
  data[[response_col]] <- prepared$y_original
  internal_names$response <- response_col

  if (isTRUE(has_offset)) {
    offset_col <- mfp2_internal_name("offset", used_names, preferred = "offset_")
    used_names <- c(used_names, offset_col)
    data[[offset_col]] <- offset
    internal_names$offset <- offset_col
    rhs <- if (identical(rhs, "1")) {
      paste0("offset(", offset_col, ")")
    } else {
      paste(rhs, "+", paste0("offset(", offset_col, ")"))
    }
  }

  if (isTRUE(prepared$strata_supplied)) {
    strata_col <- mfp2_internal_name("strata", used_names, preferred = "strata_")
    data[[strata_col]] <- prepared$strata_factor
    internal_names$strata <- strata_col
    rhs <- paste(rhs, "+", paste0("strata(", strata_col, ")"))
  }

  formula <- stats::as.formula(paste(response_col, "~", rhs))
  fit_args <- list(
    formula = formula,
    data = data,
    weights = weights,
    dist = family$dist,
    control = control,
    x = TRUE,
    y = TRUE
  )
  if (isTRUE(family$scale_supplied) &&
      !isTRUE(prepared$distribution_has_fixed_scale)) {
    fit_args$scale <- family$scale
  }
  # Use the same resolved parameter vector as the matrix-level candidate fits.
  # This also preserves defaults supplied by a custom distribution.
  if (!is.null(prepared$parms)) fit_args$parms <- prepared$parms

  fit <- mfp2_with_survreg_convergence_guard(
    do.call(survival::survreg, fit_args),
    fast = FALSE
  )
  fit$mfp2_internal_names <- internal_names
  fit$mfp2_strata_levels <- levels(prepared$strata_factor)
  fit$mfp2_survreg_strata <- if (isTRUE(prepared$strata_supplied)) {
    prepared$strata_factor
  } else {
    NULL
  }
  # Retain resolved model metadata for stable, model-specific print and summary
  # output. In particular, `survreg`'s numeric scale alone does not reveal
  # whether it was estimated, supplied by the user, or fixed by the selected
  # distribution.
  requested_dist <- prepared$dist_requested
  fit$mfp2_survreg_distribution <- if (is.character(requested_dist)) {
    requested_dist
  } else if (is.list(requested_dist) &&
             is.character(requested_dist$name) &&
             length(requested_dist$name) == 1L) {
    requested_dist$name
  } else {
    prepared$dist$name
  }
  fit$mfp2_survreg_scale_fixed <- isTRUE(prepared$scale > 0)
  fit$mfp2_survreg_parms <- prepared$parms

  loglik <- fit$loglik
  null_logl <- if (length(loglik) >= 2L) unname(loglik[1L]) else NA_real_
  model_logl <- unname(loglik[length(loglik)])
  rank <- sum(!is.na(fit$coefficients))
  total_df <- rank + if (prepared$scale == 0) prepared$nstrat else 0L

  result <- list(
    logl = model_logl,
    coefficients = fit$coefficients,
    rank = rank,
    df = total_df,
    weights = NULL,
    residuals = NULL
  )
  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    family_string = "survreg",
    fast = fast
  )

  if (isTRUE(keep_fit)) result$fit <- fit
  if (isTRUE(calculate_fit_statistics)) {
    result$null_logl <- null_logl
    result$null_deviance <- if (is.finite(null_logl)) -2 * null_logl else NA_real_
    result$model_deviance <- -2 * model_logl
  }
  result
}


#' Fit a Fine--Gray Subdistribution-Hazard Model
#'
#' Fits the weighted counting-process Cox representation produced once by
#' `prepare_finegray_family()`. The expanded response, expansion weights,
#' and row-index mapping are read from the family object so that the
#' expansion is not repeated for every MFP candidate model.
#'
#' @param x Numeric model matrix aligned with the expanded weighted response.
#' @param family Prepared Fine--Gray family object.
#' @param offset Numeric offset vector aligned with the expanded response,
#'   or `NULL`.
#' @param strata Optional stratification vector aligned with the expanded
#'   response.
#' @param control Normalised Cox-fit control list (see
#'   `normalize_cox_control()`).
#' @param fast Logical. `TRUE` selects the candidate hot path;
#'   `FALSE` refits with [survival::coxph()] to build the retained final
#'   model.
#' @param calculate_fit_statistics Logical. If `TRUE`, compute
#'   log-likelihood and information criteria for reporting.
#' @param keep_fit Logical. If `TRUE`, retain the fitted-model object.
#'   Defaults to `!fast`.
#' @param keep_fitted_values Logical. If `TRUE`, retain fitted linear
#'   predictors on the expanded data.
#' @param ... Additional arguments forwarded to the underlying fitter.
#'
#' @return A list of fit results used by the MFP engine.
#'
#' @keywords internal
#' @noRd
fit_finegray <- function(x,
                         family,
                         offset,
                         control,
                         method,
                         nocenter,
                         fast = TRUE,
                         calculate_fit_statistics = FALSE,
                         keep_fit = !fast,
                         has_offset = FALSE,
                         reserved_names = character()) {
  prepared <- family$prepared

  row_map <- prepared$row_map
  x_expanded <- if (is.null(x) || NCOL(x) == 0L) {
    matrix(numeric(0L), nrow = length(row_map), ncol = 0L)
  } else {
    x[row_map, , drop = FALSE]
  }
  offset_expanded <- prepared$offset_expanded

  # Baseline subdistribution hazard stratification (Zhou et al. 2011).
  # When strata_action is "both" or "baseline", prepared$strata_expanded
  # contains the expanded strata factor aligned to the pseudo-observations.
  fg_strata <- prepared$strata_expanded
  istrata_fg <- prepared$strata_expanded_codes

  if (fast) {
    fit <- mfp2_with_cox_convergence_guard(
      survival::agreg.fit(
        x = x_expanded,
        y = prepared$y,
        strata = istrata_fg,
        offset = offset_expanded,
        init = NULL,
        control = control,
        weights = prepared$weights,
        method = method,
        # agreg.fit() consults row names only when residuals are requested.
        # Candidate fits use resid = FALSE, so expanding names here would be
        # pure O(n_expanded) allocation on every candidate.
        rownames = NULL,
        resid = FALSE,
        nocenter = nocenter
      ),
      fast = TRUE
    )
  } else {
    has_predictors <- NCOL(x_expanded) > 0L
    if (has_predictors) {
      if (is.null(colnames(x_expanded)) || any(colnames(x_expanded) == "")) {
        stop("Internal error: x must have non-empty column names.", call. = FALSE)
      }
      data <- data.frame(x_expanded, check.names = FALSE)
      rhs <- paste(sprintf("`%s`", colnames(x_expanded)), collapse = " + ")
    } else {
      data <- data.frame(row.names = seq_len(length(row_map)))
      rhs <- "1"
    }

    used_names <- unique(c(names(data), as.character(reserved_names)))
    internal_names <- list(response = NULL, offset = NULL, strata = NULL, cluster = NULL)
    response_col <- mfp2_internal_name("response", used_names, preferred = "..mfp2_finegray_y")
    used_names <- c(used_names, response_col)
    data[[response_col]] <- prepared$y
    internal_names$response <- response_col

    if (isTRUE(has_offset)) {
      offset_col <- mfp2_internal_name("offset", used_names, preferred = "offset_")
      used_names <- c(used_names, offset_col)
      data[[offset_col]] <- offset_expanded
      internal_names$offset <- offset_col
      rhs <- if (identical(rhs, "1")) {
        paste0("offset(", offset_col, ")")
      } else {
        paste(rhs, "+", paste0("offset(", offset_col, ")"))
      }
    }

    # Add strata to the formula for baseline hazard stratification
    if (!is.null(fg_strata)) {
      strata_col <- mfp2_internal_name("strata", used_names, preferred = "..mfp2_fg_strata")
      used_names <- c(used_names, strata_col)
      data[[strata_col]] <- fg_strata
      internal_names$strata <- strata_col
      rhs <- paste(rhs, "+", paste0("strata(", strata_col, ")"))
    }

    cluster_col <- mfp2_internal_name("cluster", used_names, preferred = "..mfp2_subject")
    data[[cluster_col]] <- prepared$subject_id
    internal_names$cluster <- cluster_col
    rhs <- paste(rhs, "+", paste0("cluster(", cluster_col, ")"))
    formula <- stats::as.formula(paste(response_col, "~", rhs))

    fit <- mfp2_with_cox_convergence_guard(
      survival::coxph(
        formula,
        data = data,
        weights = prepared$weights,
        control = control,
        method = method,
        nocenter = nocenter,
        robust = TRUE,
        x = TRUE,
        y = TRUE
      ),
      fast = FALSE
    )
    fit$mfp2_internal_names <- internal_names
    fit$mfp2_finegray_event <- prepared$event
    fit$mfp2_finegray_row_map <- row_map
    fit$mfp2_finegray_n_original <- prepared$n_original
    # coxph stores a centered offset on the expanded pseudo-observation rows.
    # Retain the user-supplied original-row offset as separate metadata so
    # MFPI can reconstruct training-row predictions without first expanding
    # and then ambiguously collapsing that offset.
    fit$mfp2_original_offset <- offset
    # Keep one formula-compatible prediction row per original observation. The
    # expanded weighted Cox data cannot be used directly for subject-level CIFs,
    # and reconstructing it after fitting would lose rows omitted by finegray().
    training_prediction_data <- if (!is.null(x) && NCOL(x) > 0L) {
      as.data.frame(x, check.names = FALSE)
    } else {
      data.frame(row.names = seq_len(prepared$n_original))
    }
    if (isTRUE(has_offset)) {
      training_prediction_data[[internal_names$offset]] <- as.numeric(offset)
    }
    if (!is.null(prepared$strata_original)) {
      training_prediction_data[[internal_names$strata]] <- prepared$strata_original
    }
    training_prediction_data[[internal_names$cluster]] <- prepared$subject_id_original
    fit$mfp2_finegray_training_newdata <- training_prediction_data
  }

  if (length(fit$loglik) >= 2L) {
    null_logl <- fit$loglik[1L]
    model_logl <- fit$loglik[2L]
  } else {
    null_logl <- NA_real_
    model_logl <- fit$loglik[1L]
  }
  model_df <- length(fit$coefficients[!is.na(fit$coefficients)])
  result <- list(
    logl = unname(model_logl),
    coefficients = fit$coefficients,
    rank = model_df,
    df = model_df,
    weights = NULL,
    residuals = NULL
  )
  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    family_string = "finegray",
    fast = fast
  )
  if (isTRUE(keep_fit)) result$fit <- fit
  if (isTRUE(calculate_fit_statistics)) {
    result$null_logl <- if (is.finite(null_logl)) unname(null_logl) else NA_real_
    result$null_deviance <- if (is.finite(null_logl)) unname(-2 * null_logl) else NA_real_
    result$model_deviance <- unname(-2 * model_logl)
  }
  result
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
#' @param strata `NULL`, a normalized factor, or internal integer stratum codes.
#' @param reserved_names Optional character vector of names that package-created
#' response/offset/strata columns must avoid in the formula-based final refit.
#' @param control,rownames,nocenter passed to [survival::coxph.fit()].
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
                    has_offset = FALSE,
                    reserved_names = character()) {

  has_predictors <- !is.null(x) && NCOL(x) > 0

  istrata <- if (is.null(strata)) {
    NULL
  } else {
    cached <- attr(strata, "mfp2_integer_codes", exact = TRUE)
    if (!is.null(cached)) cached else as.integer(strata)
  }

  if (fast) {
    fit <- mfp2_with_cox_convergence_guard(
      survival::coxph.fit(
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
      ),
      fast = TRUE
    )
  } else {
    # Build the single formula-based final Cox refit with collision-safe names
    # for package-created columns. The Surv response itself is unchanged; only
    # the temporary data-frame name used by coxph() is allocated dynamically.
    if (!has_predictors) {
      d <- data.frame(row.names = seq_len(NROW(y)))
      rhs <- "1"
    } else {
      if (is.null(colnames(x)) || any(colnames(x) == "")) {
        stop("! Internal error: x must have non-empty column names.", call. = FALSE)
      }

      d <- data.frame(x, check.names = FALSE)
      rhs <- paste(sprintf("`%s`", colnames(x)), collapse = " + ")
    }

    # Include original user-facing predictor names as reserved names as well
    # as the transformed columns present in this final design. This prevents a
    # helper such as `y`, `offset_`, or `strata_` from reusing a real predictor
    # name even when that predictor was transformed to a different model column.
    used_names <- unique(c(names(d), as.character(reserved_names)))
    internal_names <- list(response = NULL, offset = NULL, strata = NULL)

    response_col <- mfp2_internal_name("response", used_names, preferred = "y")
    used_names <- c(used_names, response_col)
    # Assign the Surv object directly. Wrapping it in I() adds an AsIs layer
    # that prevents coxph() from recognizing the model-frame response as Surv.
    d[[response_col]] <- y
    internal_names$response <- response_col

    # Add offset only when the model was structurally specified with one.
    # This distinguishes no-offset models from all-zero offset models.
    if (isTRUE(has_offset)) {
      offset_col <- mfp2_internal_name("offset", used_names, preferred = "offset_")
      used_names <- c(used_names, offset_col)
      d[[offset_col]] <- offset
      internal_names$offset <- offset_col
      rhs <- paste(rhs, "+", paste0("offset(", offset_col, ")"))
    }

    if (!is.null(strata)) {
      strata_col <- mfp2_internal_name("strata", used_names, preferred = "strata_")
      used_names <- c(used_names, strata_col)
      d[[strata_col]] <- strata
      internal_names$strata <- strata_col
      rhs <- paste(rhs, "+", paste0("strata(", strata_col, ")"))
    }

    ff <- stats::as.formula(paste(response_col, "~", rhs))

    fit <- mfp2_with_cox_convergence_guard(
      survival::coxph(
        ff,
        data = d,
        weights = weights,
        control = control,
        method = method,
        nocenter = nocenter,
        x = TRUE,
        y = TRUE
      ),
      fast = FALSE
    )

    # Prediction must reconstruct offset/strata columns under the exact names
    # embedded in this stored formula. Retaining the metadata avoids formula
    # parsing and prevents package-created columns from colliding with predictors.
    fit$mfp2_internal_names <- internal_names
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

  # Cox does not expose a glm-style convergence flag. Iteration exhaustion is
  # handled above at warning time; the finite likelihood/df check here protects
  # every remaining candidate before its metric is used for selection.
  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    family_string = "cox",
    fast = fast
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
