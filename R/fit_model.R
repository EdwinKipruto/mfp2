#' Fit a model supported by mfp2
#'
#' Fits generalized linear, multinomial, ordinal, Cox, parametric survival,
#' Fine--Gray, and generalized estimating equation (GEE) models.
#'
#' @details
#' Candidate models use the fitting method for their family. Final fits use
#' the corresponding standard model class when the fitted object is retained.
#'
#' @param x A predictor matrix with column names. If names are missing,
#' they are generated. An intercept column may be present when
#' `x_has_intercept = TRUE`.
#' @param y The response. Generalized linear models accept responses supported
#' by their family, including two-column grouped binomial counts. Multinomial
#' responses may be class labels or a matrix of class counts. Ordinal
#' responses are ordered categories. Survival models use an appropriate
#' [survival::Surv()] response. For multinomial, ordinal, parametric survival,
#' Fine--Gray, and GEE models, the response is taken from `family` rather
#' than `y`.
#' @param family The prepared specification for the selected model family.
#' @param family_string The name of the selected model family, used to choose
#' the fitting method.
#' @param weights Case weights, where supported by the family.
#' @param offset An offset, where supported by the family.
#' @param method The method for handling tied event times in Cox and
#' Fine--Gray models. See [survival::coxph()].
#' @param strata Strata for a Cox model. Other families obtain any required
#' stratification information from their prepared family specification.
#' @param control The family-specific fitting controls, prepared before
#' candidate models are fitted.
#' @param rownames Row names passed to the Cox candidate fitter.
#' @param nocenter Values that identify Cox and Fine--Gray predictor columns
#' that should not be centered internally.
#' @param fast Logical. If `TRUE` (the default), fit candidate models with
#' lower-level functions such as [stats::glm.fit()],
#' [survival::coxph.fit()], [survival::survreg.fit()],
#' [nnet::nnet.default()], [rms::orm.fit()], or
#' [geepack::geese.fit()], according to the family. Fine--Gray
#' candidates use [survival::agreg.fit()]. If `FALSE`, fit the final
#' model with the corresponding full-model function. Generalized linear
#' models may use `fastglm`; negative-binomial models use `fastglm_nb`
#' for both candidate and final fits.
#' @param calculate_fit_statistics Whether to include null-model and
#' fitted-model statistics for reporting.
#' @param calculate_gaussian_deviance Whether to calculate the Gaussian
#' deviance used for F-tests. Applies only to Gaussian generalized linear
#' models.
#' @param keep_fit Whether to include the fitted model object in the result.
#' Defaults to `!fast`.
#' @param has_offset Whether an offset was supplied before missing offsets
#' were replaced internally by zeros.
#' @param x_has_intercept Whether `x` already has a leading intercept column.
#' This must be `FALSE` for Cox and Fine--Gray models.
#' @param fitter The fitting method for generalized linear models. Ignored
#' by the other families.
#' @param keep_fitted_values Whether to return fitted values separately from
#' the fitted object. Supported for generalized linear, multinomial, and
#' GEE models.
#' @param keep_coefficients Whether to calculate and return coefficients for
#' a multinomial candidate fit. Coefficients are still calculated when
#' `keep_fit = TRUE`.
#' @param multinomial_optimizer Previously prepared information used across
#' multinomial candidates with the same design dimensions.
#' @param gee_alpha_start An optional starting value for the working
#' correlation in a GEE candidate fit.
#' @param gee_column_scales Optional column scaling factors for a GEE
#' candidate whose predictor matrix has already been scaled.
#' @param gee_selection_criterion The selection criterion for a GEE
#' candidate. AIC and BIC candidates omit covariance calculations they
#' do not need.
#' @param gee_allow_failed_candidate Whether an unsuccessful GEE candidate
#' may be returned as an ineligible result during a fractional polynomial
#' power search. This does not apply to the final fit.
#' @param reserved_names Names that cannot be used for temporary columns
#' created during a final fit.
#'
#' @return A list containing `logl`, `df`, and `rank`, along with
#' family-specific results. `logl` is the fitted log likelihood, except
#' for GEE models, where it contains a quasi-likelihood measure. An
#' unsuccessful GEE candidate permitted by `gee_allow_failed_candidate`
#' has `logl = NA`.
#'
#' `df` is the number of parameters used for model comparisons; `rank`
#' describes the fitted regression parameters according to the family.
#' These can differ, for example, when a survival model has additional
#' scale parameters or an ordinal model has threshold parameters.
#'
#' The list may also contain:
#'
#' * `coefficients`: fitted regression coefficients. These may be omitted
#'   for a multinomial candidate when `keep_coefficients = FALSE` and
#'   `keep_fit = FALSE`.
#' * `deviance_gaussian`: Gaussian deviance when requested.
#' * `null_deviance` and `model_deviance`: family-specific deviances when
#'   `calculate_fit_statistics = TRUE`. For a fast GEE fit,
#'   `null_deviance` is `NA`.
#' * `null_logl`: null-model log likelihood when fit statistics are
#'   requested for a multinomial, ordinal, Cox, Fine--Gray, or parametric
#'   survival model.
#' * `fitted_values`: fitted values when requested for a supported family.
#' * `fit`: the fitted object when `keep_fit = TRUE`.
#' * `transformed_to_model_columns`: for a final fit, the mapping from
#'   predictor column names to fitted coefficient names.
#'
#' GEE results also include information used for model selection and
#' working-correlation calculations.
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
                      gee_alpha_start = NULL,
                      gee_column_scales = NULL,
                      gee_selection_criterion = NULL,
                      gee_allow_failed_candidate = FALSE,
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
  } else if (identical(family_string, "gee")) {
    if (isTRUE(calculate_gaussian_deviance)) {
      stop(
        "Internal error: Gaussian F-test deviance is unavailable for GEE models.",
        call. = FALSE
      )
    }

    fit <- fit_gee(
      x = x,
      y = y,
      family = family,
      weights = weights,
      offset = offset,
      control = control,
      fast = fast,
      calculate_fit_statistics = calculate_fit_statistics,
      keep_fit = keep_fit,
      keep_fitted_values = keep_fitted_values,
      has_offset = has_offset,
      x_has_intercept = x_has_intercept,
      alpha_start = gee_alpha_start,
      column_scales = gee_column_scales,
      selection_criterion = gee_selection_criterion,
      allow_failed_candidate = gee_allow_failed_candidate,
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
#' Signal a Model-Fitting Failure
#'
#' Candidate fits (`fast = TRUE`) signal a classed condition,
#' `mfp2_candidate_failure`, so that the FP power search can skip that one
#' candidate. Final fits signal an ordinary error.
#'
#' @param message Error message.
#' @param fast Whether the failing fit is a repeated candidate fit.
#' @keywords internal
#' @noRd
mfp2_stop_fit <- function(message, fast = TRUE) {
  classes <- if (isTRUE(fast)) {
    c("mfp2_candidate_failure", "error", "condition")
  } else {
    c("error", "condition")
  }
  stop(structure(class = classes, list(message = message, call = NULL)))
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
      mfp2_stop_fit(
        paste0(
          "The ", stage, family_label,
          " model did not converge; MFP fitting cannot use an unconverged fit. ",
          "Adjust the fitting controls or inspect the data."
        ),
        fast = fast
      )
    }
  }

  criterion_label <- if (isTRUE(fast)) "candidate criterion" else "fit criterion"

  if (!is.numeric(logl) || length(logl) != 1L || !is.finite(logl)) {
    mfp2_stop_fit(
      paste0(
        "The ", stage, family_label,
        " model produced a non-finite log-likelihood; MFP fitting cannot ",
        "continue with an invalid ", criterion_label, "."
      ),
      fast = fast
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

#' Fit a generalized linear model
#'
#' Fits a candidate model directly from a predictor matrix, or fits the
#' final model using [stats::glm()]. Negative-binomial models use
#' `fastglm::fastglm_nb()` in both cases.
#'
#' @param x A predictor matrix with one row per observation. If
#'   `x_has_intercept = FALSE`, an intercept column is added for fitting.
#' @param y The response accepted by the selected family. Binomial
#'   responses may be a two-column matrix of success and failure counts.
#' @param family A GLM family object, such as `stats::gaussian()`.
#' @param weights Optional prior weights, one per observation.
#' @param offset Optional known values added to the linear predictor,
#'   one per observation.
#' @param control GLM fitting controls. With `fastglm`, `epsilon` is passed
#'   as `tol` and `maxit` is passed through. `trace = TRUE` is not supported
#'   with `fastglm`.
#' @param fast Logical. If `TRUE` (the default), fit a candidate directly
#'   from its predictor matrix using [stats::glm.fit()] or, when requested,
#'   `fastglm::fastglm()`. If `FALSE`, fit the final model with
#'   [stats::glm()]. Negative-binomial models use `fastglm::fastglm_nb()`
#'   for both candidate and final fits.
#' @param calculate_fit_statistics Whether to return the null-model and
#'   fitted-model deviances.
#' @param calculate_gaussian_deviance Whether to calculate the Gaussian
#'   deviance used for F-tests. Applies only to Gaussian models.
#' @param keep_fit Whether to include the fitted object in the result.
#'   Defaults to `!fast`.
#' @param has_offset Whether an offset was supplied and should appear in
#'   the formula for a final fit.
#' @param x_has_intercept Whether `x` already begins with a column named
#'   `"(Intercept)"`.
#' @param fitter The candidate fitting method: `"base"` or `"fastglm"`.
#'   Negative-binomial models always use `fastglm::fastglm_nb()`.
#' @param family_string The family name supplied by `fit_model()`.
#' @param keep_fitted_values Whether to return fitted values separately
#'   from the fitted object.
#' @param reserved_names Names that cannot be used for temporary response
#'   or offset columns created for a final fit.
#'
#' @return A list containing:
#'
#' * `logl`: the fitted log likelihood.
#' * `coefficients`: the regression coefficients.
#' * `rank`: the fitted regression rank.
#' * `df`: the number of parameters used for model comparisons. This
#'   includes an estimated dispersion parameter where applicable.
#'
#' When requested, the list also contains `deviance_gaussian`, the
#' Gaussian deviance used for F-tests; `null_deviance` and
#' `model_deviance`, the deviances reported by the fitting method;
#' `fitted_values`, the fitted values; and `fit`, the fitted object.
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
  negbin_link <- if (is_negbin) {
    if (is.null(fit_flags$negbin_link)) {
      mfp2_negbin_link(family)
    } else {
      fit_flags$negbin_link
    }
  } else {
    NULL
  }
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
      fit_glm_fastglm_nb(xx, y, weights, offset, control,
                         link = negbin_link, fast = TRUE)
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
    fit <- fit_glm_fastglm_nb(xx, y, weights, offset, control,
                              link = negbin_link, fast = FALSE)

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
        data    = data,
        family  = family,
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
  if (identical(family_string, "negbin")) {
    fitter <- "fastglm"
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

#' Fit a candidate generalized linear model with fastglm
#'
#' Fits the model directly from a predictor matrix using
#' [fastglm::fastglm()]. Checks that the result contains the statistics
#' needed for model selection. Gaussian fits must also contain residuals
#' and weights needed for F-tests.
#'
#' If fitting fails or the result is incomplete, fitting stops. The
#' candidate is not fitted again with [stats::glm.fit()].
#'
#' @param x Predictor matrix, including an intercept column.
#' @param y Response vector or matrix.
#' @param family GLM family object.
#' @param weights Prior weights.
#' @param offset Values added to the linear predictor.
#' @param control GLM controls. `epsilon` is passed to fastglm as `tol`;
#'   `maxit` is passed through.
#'
#' @return The fitted `fastglm` object.
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

#' Fit a negative-binomial model with fastglm
#'
#' Fits candidate and final models with [fastglm::fastglm_nb()] using the
#' requested link. Checks that the result contains finite values needed for
#' model selection and counts the estimated dispersion parameter in the AIC.
#' There is no alternative fitting method for negative-binomial models.
#'
#' The square-root and identity links have a restricted feasible domain:
#' `sqrt` requires a positive linear predictor and `identity` a positive mean
#' for every observation. A fit outside that domain, or a backend error for
#' one of these links, is reported through `mfp2_stop_fit()`, so a single FP
#' power candidate is skipped while null, linear, and final fits stop. The
#' link is never changed and fitted values are never clipped. Log-link
#' backend errors keep their original strict behaviour.
#'
#' @param link Negative-binomial link: `"log"`, `"sqrt"`, or `"identity"`.
#' @param fast Whether this is a repeated candidate fit.
#'
#' @keywords internal
#' @noRd
fit_glm_fastglm_nb <- function(x, y, weights, offset, control,
                               link = "log", fast = TRUE) {
  link <- mfp2_validate_negbin_link(link)
  restricted_link <- !identical(link, "log")
  stage <- if (isTRUE(fast)) "candidate" else "final"

  fit <- tryCatch(
    fastglm::fastglm_nb(
      x = x,
      y = y,
      weights = weights,
      offset = offset,
      link = link,
      method = 0L,
      # Apply glm.control() to the inner IRLS fit. fastglm_nb-specific outer
      # optimization controls keep their documented defaults.
      tol = control$epsilon,
      maxit = control$maxit
    ),
    error = function(e) {
      message <- paste0(
        "fastglm could not fit the requested negative-binomial model (",
        link, " link): ", conditionMessage(e)
      )
      if (restricted_link) {
        mfp2_stop_fit(message, fast = fast)
      }
      stop(message, call. = FALSE)
    }
  )

  # Guard against a backend that ignores the requested link.
  fitted_link <- if (is.list(fit$family)) fit$family$link else NULL
  if (!identical(fitted_link, link)) {
    stop(
      "fastglm returned a negative-binomial fit with link '",
      if (is.character(fitted_link)) fitted_link else "unknown",
      "' instead of the requested '", link, "' link. ",
      "Update fastglm to a version whose `fastglm_nb()` supports `link`.",
      call. = FALSE
    )
  }

  # The square-root and identity links can leave their feasible domain. Check
  # it before any likelihood-based validation, because an infeasible mean also
  # makes the log-likelihood non-finite.
  if (restricted_link) {
    eta <- fit$linear.predictors
    mu <- fit$fitted.values
    domain_ok <- is.numeric(eta) && length(eta) == NROW(y) &&
      all(is.finite(eta)) &&
      is.numeric(mu) && length(mu) == NROW(y) && all(is.finite(mu)) &&
      all(mu > 0) &&
      (!identical(link, "sqrt") || all(eta > 0))
    if (!domain_ok) {
      requirement <- if (identical(link, "sqrt")) {
        "a positive linear predictor"
      } else {
        "a positive mean"
      }
      mfp2_stop_fit(
        paste0(
          "The ", stage, " negative-binomial model with the ", link,
          " link has fitted values outside the link's valid domain; the ",
          link, " link requires ", requirement, " for every observation. ",
          "The link was not changed and fitted values were not adjusted."
        ),
        fast = fast
      )
    }
  }

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
      "no base fitter fallback is available for negative-binomial models.",
      call. = FALSE
    )
  }

  # The rank + 1 AIC convention above matches MASS::glm.nb() and counts
  # theta without over-counting aliased regression columns.
  fit
}

#' Calculate negative-binomial null deviance with an offset
#'
#' Uses [stats::glm.fit()] to fit an intercept-only model with the supplied
#' offset and weights. The `family` argument comes from the full
#' negative-binomial fit: it contains that fit's estimated dispersion
#' parameter, theta, in its functions. This calculation holds theta fixed
#' rather than estimating it again.
#'
#' This helper is called only when fit statistics are requested and the
#' offset contains a nonzero value. It returns the intercept-only fit's
#' deviance and stops if that fit does not converge or has a non-finite
#' deviance.
#'
#' @param y Numeric vector of finite, nonnegative integer counts.
#' @param weights Numeric vector of positive observation weights.
#' @param offset Numeric vector aligned with `y`.
#' @param family Negative-binomial family from the full fit, constructed
#'   using its estimated theta.
#' @param control Controls passed to [stats::glm.fit()].
#'
#' @return A numeric value: the deviance of the intercept-only model
#'   fitted with the offset and the full model's theta.
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

#' Check whether the selected fitter supports GLM tracing
#'
#' Rejects `control$trace = TRUE` when a GLM uses fastglm. This includes
#' negative-binomial models, which always use fastglm. Called after fitting
#' controls and the fitter have been resolved, before candidate fitting begins.
#'
#' @param control Family-specific fitting controls. May be `NULL`.
#' @param family_string Name of the selected model family.
#' @param fitter Resolved GLM fitter, `"base"` or `"fastglm"`.
#'
#' @return The supplied `control` object, invisibly. An unsupported trace
#'   setting raises an error.
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


#' Prepare fitting controls for the selected model family
#'
#' Fills in family-specific defaults and checks the supplied control settings
#' before candidate models are fitted. For GLM and negative-binomial models,
#' it also rejects iteration tracing when the fitter is fastglm. The prepared
#' controls are reused during model selection.
#'
#' @param control A list of fitting controls, or `NULL` to use the defaults
#'   for the selected family.
#' @param family_string Name of the selected model family.
#' @param fitter Resolved GLM fitter, `"base"` or `"fastglm"`. Negative-binomial
#'   models use fastglm regardless of this argument.
#'
#' @return A validated control list for the selected family: GLM and
#'   negative-binomial, Cox and Fine--Gray, parametric survival, multinomial,
#'   ordinal, or GEE.
#'
#' @keywords internal
#' @noRd
normalize_fit_control <- function(control = NULL, family_string,
                                  fitter = "base") {
  normalized <- if (identical(family_string, "gee")) {
    normalize_gee_control(control)
  } else if (family_string %in% c("cox", "finegray")) {
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

#' Prepare controls for ordinal model fitting
#'
#' Fills in and validates the controls used for ordinal fits. Accepts
#' `epsilon` as an alternative name for `eps`; if both are supplied, `eps`
#' takes precedence. Other control names are rejected.
#'
#' @param control A list containing `maxit`, `eps`, `tol`, or `trace`, or
#'   `NULL` to use the defaults. `epsilon` may be supplied in place of `eps`.
#'
#' @return A list containing `maxit`, `eps`, `tol`, and `trace`. Their defaults
#'   are `30`, `0.005`, `1e-7`, and `FALSE`, respectively.
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

  if (is.null(control)) {
    return(defaults)
  }

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

#' Stop when an ordinal fit warns of convergence failure
#'
#' Evaluates an ordinal fitting expression and turns warnings about
#' non-convergence or a singular fit into errors. Other warnings retain
#' their usual behavior. Used around both [rms::orm.fit()] candidate fits
#' and [rms::orm()] final fits.
#'
#' @param expr Expression that fits the ordinal model.
#' @param fast Whether this is a candidate fit. Defaults to `TRUE`;
#'   `FALSE` identifies a final fit in the error.
#'
#' @return The value of `expr` if no matching warning occurs.
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
        mfp2_stop_fit(
          paste0(
            "The ", stage, " ordinal model did not converge; adjust `control` ",
            "(for example raise `maxit`) or inspect the data."
          ),
          fast = fast
        )
      }
    }
  )
}

#' Fit an ordinal regression model
#'
#' Fits ordered outcomes using the link specified in `family`. The response
#' and its category order are prepared before this function is called.
#' Every fit uses [rms::orm.fit()] on the predictor matrix. When `keep_fit`
#' is `TRUE`, the function also fits and retains an [rms::orm()] object.
#'
#' @param x Numeric predictor matrix, or `NULL` for an intercept-only model.
#' @param family Prepared ordinal family containing the response, link, and
#'   any offset.
#' @param control Validated ordinal fitting controls.
#' @param fast Whether this is a candidate fit. Defaults to `TRUE`.
#' @param calculate_fit_statistics Whether to return the null-model log
#'   likelihood and the null and fitted-model deviances.
#' @param keep_fit Whether to include an `rms::orm()` fitted object in the
#'   result. Defaults to `!fast`.
#' @param x_has_intercept Whether `x` begins with a column named
#'   `"(Intercept)"`. That column is removed because `orm.fit()` fits
#'   the ordinal intercepts itself.
#' @param reserved_names Names that temporary response and offset columns
#'   must avoid when creating a retained fit.
#'
#' @return A list containing `logl`, the fitted log likelihood;
#'   `coefficients`, the predictor coefficients without the ordinal
#'   intercepts; `rank`, the number of predictor coefficients; and `df`,
#'   the number of predictor coefficients plus ordinal intercepts.
#'   `weights` contains the supplied case weights; `residuals` is `NULL`.
#'
#'   When `keep_fit = TRUE`, the list also contains `fit`, an `rms::orm()`
#'   object. When `calculate_fit_statistics = TRUE`, it contains
#'   `null_logl`, `null_deviance`, and `model_deviance`.
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
  weights <- prepared$weights
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
    weights = weights,
    normwt = FALSE,
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
  # With an offset, orm.fit() records intercepts-only, intercepts+offset,
  # and (when predictors exist) intercepts+offset+predictors deviances.
  # Compare the selected model with the null model that retains its offset.
  null_dev <- if (has_offset) {
    if (has_predictors) unname(dev[length(dev) - 1L]) else model_dev
  } else if (length(dev) >= 2L) {
    unname(dev[1L])
  } else {
    NA_real_
  }
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
    weights = weights,
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
    # Pass weights as values so a predictor named "weights" cannot mask them
    # during the formula fitter's model-frame evaluation.
    fit_args <- list(
      formula = formula, data = dfm, family = link, x = TRUE, y = TRUE,
      normwt = FALSE,
      maxit = control$maxit, eps = control$eps, tol = control$tol,
      trace = isTRUE(control$trace)
    )
    # Explicit unit weights trigger rms's weighted-validation warning even
    # though the fit is unweighted. Omit that argument for this equivalent case.
    if (any(weights != 1)) fit_args$weights <- weights
    fitobj <- mfp2_with_ordinal_convergence_guard(
      do.call(rms::orm, fit_args),
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

#' Prepare controls for multinomial model fitting
#'
#' Fills in and validates the controls passed to the nnet fitting functions.
#' Accepts `epsilon` as an alternative name for `reltol`; if both are
#' supplied, `reltol` takes precedence. Other control names are rejected.
#'
#' @param control A list containing `maxit`, `reltol`, `abstol`, `trace`,
#'   or `MaxNWts`, or `NULL` to use the defaults. `epsilon` may be supplied
#'   in place of `reltol`.
#'
#' @return A list containing `maxit`, `reltol`, `abstol`, `trace`, and
#'   `MaxNWts`. Their defaults are `300`, `1e-8`, `1e-4`, `FALSE`, and
#'   `1000`, respectively.
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
  if (is.null(control)) {
    return(defaults)
    }

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


#' Flatten multinomial coefficients
#'
#' Converts a coefficient matrix into the named vector returned by mfp2's
#' multinomial fitter. Coefficients are grouped by non-reference outcome,
#' with all predictor coefficients for one outcome followed by those for
#' the next.
#'
#' @param coefficient_matrix Numeric matrix with one row per non-reference
#'   outcome and one column per model coefficient. Rows and columns must
#'   have names.
#'
#' @return A numeric vector of length
#'   `nrow(coefficient_matrix) * ncol(coefficient_matrix)`, with names
#'   in the form `"<outcome>::<coefficient>"`.
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


#' Fit a multinomial model with nnet::multinom()
#'
#' Builds a formula from the supplied design matrix and fits the model with
#' [nnet::multinom()]. The response, case weights, and any class-specific
#' offsets come from the prepared `family` object. This helper is used for
#' retained fits and for null-model statistics.
#'
#' @param xx Numeric design matrix whose first column is named
#'   `"(Intercept)"`.
#' @param family Prepared multinomial family containing the response,
#'   weights, and any offsets.
#' @param control Validated multinomial fitting controls.
#' @param Hess Whether to include a coefficient information matrix. With
#'   class-specific offsets, the helper calculates it from the fitted
#'   probabilities; otherwise it requests it from `nnet::multinom()`.
#'
#' @return A fitted `nnet::multinom()` object, with an information matrix
#'   when `Hess = TRUE`.
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


#' Calculate the information matrix for a multinomial fit
#'
#' Calculates the coefficient information matrix from the design matrix,
#' fitted class probabilities, and case weights. Class-specific offsets
#' affect the fitted probabilities but are not treated as estimated
#' coefficients.
#'
#' Parameters are ordered by non-reference class: all design coefficients
#' for one class, followed by those for the next.
#'
#' @param x Numeric design matrix with one row per observation.
#' @param probabilities Matrix of fitted probabilities, with one row per
#'   observation and one column per class.
#' @param weights Numeric vector of case weights.
#' @param outcomes Character vector naming the non-reference classes.
#' @param class_levels Character vector of all class names in fitted order,
#'   with the reference class first.
#'
#' @return A square information matrix with
#'   `length(outcomes) * ncol(x)` rows and columns. Its row and column names
#'   have the form `"<class>:<design column>"`.
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


#' Prepare nnet weights for multinomial candidate fits
#'
#' Builds the parameter mask used by [nnet::nnet.default()] to identify
#' which coefficients may be estimated. When class-specific offsets are
#' present, it also builds fixed weights that carry those offsets into
#' the fit. Candidates with the same design width, number of classes,
#' and offset setting can reuse this structure.
#'
#' @param p Number of design columns, including the intercept.
#' @param c_classes Number of response classes.
#' @param has_offset Whether class-specific offsets are present.
#' @param MaxNWts Requested maximum number of nnet weights.
#'
#' @return A list containing `p`, `n_classes`, `has_offset`, the logical
#'   parameter `mask`, `fixed_offset_weights` (`NULL` when there is no
#'   offset), and `MaxNWts`. The returned `MaxNWts` is at least the
#'   length of the mask.
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

#' Fit a multinomial model
#'
#' Fits candidate models directly from a design matrix with
#' [nnet::nnet.default()]. When `fast = FALSE`, fits through
#' [nnet::multinom()]. The response, case weights, and any class-specific
#' offsets come from the prepared `family` object.
#'
#' @param x Numeric predictor matrix.
#' @param family Prepared multinomial family containing the response,
#'   weights, and any offsets.
#' @param control Validated multinomial fitting controls.
#' @param fast Whether to use `nnet.default()` for a candidate fit.
#'   If `FALSE`, use `nnet::multinom()`. Defaults to `TRUE`.
#' @param calculate_fit_statistics Whether to fit a null model and return
#'   its log likelihood and the null and fitted-model deviances.
#' @param keep_fit Whether to include the fitted object in the result.
#'   Defaults to `!fast`.
#' @param keep_fitted_values Whether to return fitted class probabilities.
#' @param keep_coefficients Whether to return the coefficient vector and
#'   matrix. Coefficients are always returned when `keep_fit = TRUE`.
#' @param x_has_intercept Whether `x` already begins with a column named
#'   `"(Intercept)"`; otherwise an intercept is added.
#' @param optimizer_structure Optional prepared parameter mask and offset
#'   weights for a candidate fit with matching design dimensions.
#'
#' @return A list containing `logl`, `rank`, `df`, `class_levels`,
#'   `reference_class`, and `n_logits`. It also contains `coefficients`
#'   and `coefficient_matrix` when requested; `fitted_values` and `fit`
#'   when requested; and `null_logl`, `null_deviance`, and
#'   `model_deviance` when fit statistics are requested.
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
    # nnet.default() does not report the design rank. A transformed candidate
    # may have aliased columns, so check its rank before fitting and reject it
    # if the design is rank deficient.
    design_rank <- qr(xx)$rank
    if (design_rank != ncol(xx)) {
      mfp2_stop_fit("The multinomial candidate design is rank deficient.", fast = TRUE)
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
    # nnet stores fitted values and residuals in the internal class order
    # (reference first). Report them in the original response order, as
    # predict(type = "response") does, so fitted() and predict() agree.
    original_levels <- prepared$original_levels
    for (component in c("fitted.values", "residuals")) {
      values <- native_fit[[component]]
      if (is.matrix(values) && all(original_levels %in% colnames(values))) {
        native_fit[[component]] <- values[, original_levels, drop = FALSE]
      }
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

#' Prepare controls for parametric survival fitting
#'
#' Applies the defaults and validates the supplied settings using
#' [survival::survreg.control()]. The resulting list is reused when fitting
#' candidate and final models.
#'
#' @param control A list of settings accepted by
#'   [survival::survreg.control()], or `NULL` to use its defaults.
#'
#' @return A validated list of `survreg` fitting controls.
#'
#' @keywords internal
#' @noRd
normalize_survreg_control <- function(control = NULL) {
  if (is.null(control)) {
    return(survival::survreg.control())
  }

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

#' Prepare controls for Cox and Fine--Gray fitting
#'
#' Applies the defaults and validates the supplied settings using
#' [survival::coxph.control()]. The resulting list is used for Cox and
#' Fine--Gray candidate and final fits.
#'
#' @param control A list of settings accepted by
#'   [survival::coxph.control()], or `NULL` to use its defaults.
#'
#' @return A validated list of Cox fitting controls.
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

#' Stop when a survreg fit runs out of iterations
#'
#' Evaluates a parametric survival fitting expression. If
#' [survival::survreg()] or [survival::survreg.fit()] warns that it ran out
#' of iterations without converging, turns that warning into an error.
#' Other warnings retain their usual behavior.
#'
#' @param expr Expression that fits the parametric survival model.
#' @param fast Whether this is a candidate fit. Defaults to `TRUE`;
#'   `FALSE` identifies a final fit in the error.
#'
#' @return The value of `expr` if the specified warning does not occur.
#'
#' @keywords internal
#' @noRd
mfp2_with_survreg_convergence_guard <- function(expr, fast = TRUE) {
  stage <- if (isTRUE(fast)) "candidate" else "final"
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Ran out of iterations and did not converge", conditionMessage(w), fixed = TRUE)) {
        mfp2_stop_fit(
          paste0("The ", stage, " survreg model did not converge; adjust `control` or inspect the data."),
          fast = fast
        )
      }
    }
  )
}

#' Stop when a Cox or Fine--Gray fit runs out of iterations
#'
#' Turns an iteration-exhaustion warning into an error that identifies
#' the model being fitted. Other warnings retain their usual behavior.
#'
#' @param expr Fitting expression to evaluate.
#' @param fast Whether this is a candidate fit. Defaults to `TRUE`.
#' @param model Model name to use in the error message.
#'
#' @return The value of `expr` if the specified warning does not occur.
#' @keywords internal
#' @noRd
mfp2_with_cox_convergence_guard <- function(expr, fast = TRUE, model = "Cox") {
  stage <- if (isTRUE(fast)) "candidate" else "final"

  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("Ran out of iterations and did not converge",
                msg, fixed = TRUE)) {
        mfp2_stop_fit(
          paste0(
            "The ", stage, " ", model,
            " model did not converge within the configured iteration limit. ",
            "Adjust `control$iter.max` or inspect the data."
          ),
          fast = fast
        )
      }
    }
  )
}

#' Fit a parametric survival model
#'
#' Fits candidate models from a design matrix with [survival::survreg.fit()].
#' When `fast = FALSE`, fits the final model with [survival::survreg()].
#' The prepared `family` contains the response, distribution, scale, and
#' any stratification information.
#'
#' @param x Numeric predictor matrix, or `NULL` for an intercept-only model.
#' @param y Response argument retained for the common fitting interface.
#'   This function uses the response stored in `family`.
#' @param family Prepared parametric survival family.
#' @param weights Observation weights.
#' @param offset Values added to the linear predictor.
#' @param control Validated `survreg` fitting controls.
#' @param fast Whether to use `survreg.fit()` for a candidate fit.
#'   If `FALSE`, use `survreg()`. Defaults to `TRUE`.
#' @param calculate_fit_statistics Whether to return the null-model log
#'   likelihood and the null and fitted-model deviances.
#' @param keep_fit Whether to include the fitted object in the result.
#'   Defaults to `!fast`.
#' @param has_offset Whether an offset was supplied and should appear in
#'   the formula for a final fit.
#' @param x_has_intercept Whether `x` already begins with a column named
#'   `"(Intercept)"`; otherwise an intercept is added for fitting.
#' @param reserved_names Names that temporary response, offset, and strata
#'   columns must avoid when creating a final fit.
#'
#' @return A list containing `logl`, `coefficients`, `rank`, and `df`.
#'   `df` includes estimated scale parameters where applicable.
#'   `weights` and `residuals` are `NULL`.
#'
#'   When `keep_fit = TRUE`, the list also contains `fit`. When
#'   `calculate_fit_statistics = TRUE`, it contains `null_logl`,
#'   `null_deviance`, and `model_deviance`.
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

#' Fit a Fine--Gray subdistribution hazard model
#'
#' Fits the weighted Cox representation prepared for a Fine--Gray model.
#' The prepared `family` contains the expanded response, weights, offset,
#' strata, and mapping from expanded rows to original observations.
#' Candidate fits use [survival::agreg.fit()]; final fits use
#' [survival::coxph()].
#'
#' @param x Numeric predictor matrix with one row per original observation,
#'   or `NULL` for a model without predictors.
#' @param family Prepared Fine--Gray family.
#' @param offset Offset aligned with the original observations, or `NULL`.
#'   The expanded offset used for fitting comes from `family`.
#' @param control Validated Cox fitting controls.
#' @param method Method for handling tied event times.
#' @param nocenter Values identifying columns that the Cox fitter should
#'   leave uncentered.
#' @param fast Whether to use `agreg.fit()` for a candidate fit. If `FALSE`,
#'   use `coxph()` for the final fit. Defaults to `TRUE`.
#' @param calculate_fit_statistics Whether to return the null partial log
#'   likelihood and the null and fitted-model deviances.
#' @param keep_fit Whether to include the fitted object in the result.
#'   Defaults to `!fast`.
#' @param has_offset Whether an offset was supplied and should appear in
#'   the final fit's formula.
#' @param reserved_names Names that temporary formula columns must avoid
#'   when creating a final fit.
#'
#' @return A list containing `logl`, `coefficients`, `rank`, and `df`.
#'   `weights` and `residuals` are `NULL`. When requested, the list also
#'   contains `fit`, or `null_logl`, `null_deviance`, and
#'   `model_deviance`.
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
      fast = TRUE,
      model = "Fine-Gray"
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
      fast = FALSE,
      model = "Fine-Gray"
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

#' Fit a Cox proportional hazards model
#'
#' Fits candidate models directly from a predictor matrix with
#' [survival::coxph.fit()]. When `fast = FALSE`, fits the final model
#' with [survival::coxph()].
#'
#' @param x Numeric predictor matrix without an intercept, or `NULL`
#'   for a model without predictors.
#' @param y A [survival::Surv()] response.
#' @param strata Optional strata. Candidate fits use integer stratum codes;
#'   final fits include strata in the model formula.
#' @param weights Observation weights.
#' @param offset Values added to the linear predictor.
#' @param control Validated Cox fitting controls.
#' @param method Method for handling tied event times. See
#'   [survival::coxph()].
#' @param rownames Row names passed to `coxph.fit()` for candidate fits.
#' @param nocenter Values identifying predictor columns that the Cox fitter
#'   should leave uncentered.
#' @param fast Whether to use `coxph.fit()` for a candidate fit. If `FALSE`,
#'   use `coxph()` for the final fit. Defaults to `TRUE`.
#' @param calculate_fit_statistics Whether to return the null partial log
#'   likelihood and the null and fitted-model statistics.
#' @param keep_fit Whether to include the fitted object in the result.
#'   Defaults to `!fast`.
#' @param has_offset Whether an offset was supplied and should appear in
#'   the final fit's formula.
#' @param reserved_names Names that temporary response, offset, and strata
#'   columns must avoid when creating a final fit.
#'
#' @return A list containing `logl`, the fitted partial log likelihood;
#'   `coefficients`; and `rank` and `df`, both equal to the number of
#'   non-missing coefficients. `weights` and `residuals` are `NULL`.
#'
#'   When `keep_fit = TRUE`, the list also contains `fit`. When
#'   `calculate_fit_statistics = TRUE`, it contains `null_logl` and
#'   `null_deviance` and `model_deviance`, calculated as minus twice
#'   the respective partial log likelihoods.
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
