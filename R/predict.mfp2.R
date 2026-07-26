#' Predict from an `mfp2` Model
#'
#' Obtain predictions from a model fitted with [mfp2()]. The function can
#' return predictions from the complete model, term-specific fitted values,
#' or differences from chosen reference values.
#'
#' @details
#' Supply predictors using their original values. The fitted shifts, FP or ACD
#' functions, centering, zero handling, and spike-at-zero decisions are applied
#' automatically.
#'
#' Do not create fractional-polynomial basis columns, centered columns, or
#' structural-zero indicator columns yourself.
#'
#' Standard errors and confidence intervals are conditional on the final fitted
#' model. They do not include uncertainty from variable selection or from
#' choosing the fitted function.
#'
#' @section Supplying prediction data:
#' Supply `newdata` in the same form used when fitting the model.
#'
#' For a model fitted with a formula, normally supply the original variables.
#' For example, supply the factor `race`, not manually created contrast columns
#' such as `raceB` and `raceC`. Use factor levels that were available when the
#' model was fitted. The required formula terms and factor contrasts are
#' recreated automatically.
#'
#' For a model fitted with a numeric matrix, supply the original numeric columns
#' with the same names and coding used in the fitting matrix.
#'
#' If several matrix columns were grouped using `term_groups`, supply every
#' member column. Use the group name, rather than the member-column names, in
#' the `terms` argument.
#'
#' For complete-model prediction, normally supply only the variables needed by
#' terms included in the final model. For term or contrast prediction, only the
#' requested terms are required. Additional variables can be needed to evaluate
#' a formula offset, Cox strata, or an absolute Cox prediction.
#'
#' Required values must not be missing. No numeric column in `newdata` may
#' contain `Inf` or `-Inf`, including unused extra columns. Prediction stops if
#' a supplied value is outside the valid range of the fitted transformation.
#' For example, some logarithmic or negative-power functions require a
#' positive value after the fitted shift is applied.
#'
#' @section Choosing the prediction type:
#' For Gaussian, binomial, Poisson, and negative-binomial models:
#'
#' - `type = "link"` returns predictions on the link scale.
#' - `type = "response"` returns predictions on the response scale.
#' - `type = "lp"` is accepted as another name for `"link"`.
#'
#' For Cox models:
#'
#' - `type = "lp"` returns the relative log-hazard.
#' - `type = "risk"` returns the relative hazard, calculated as `exp(lp)`.
#' - `type = "expected"` returns the predicted cumulative hazard up to the
#'   supplied follow-up time.
#' - `type = "survival"` returns the estimated survival probability at the
#'   supplied follow-up time.
#' - `type = "link"` is accepted as another name for `"lp"`.
#'
#' For all supported models:
#'
#' - `type = "terms"` returns a fitted value for each requested term on the
#'   linear-predictor scale. For Gaussian, binomial, Poisson, and negative-binomial models, the
#'   intercept is included by default.
#' - `type = "contrasts"` compares each fitted term value with a reference
#'   value.
#'
#' When `type` is not supplied, the default is `"link"` for Gaussian,
#' binomial, and Poisson models and `"lp"` for Cox models.
#'
#' @section Term predictions:
#' A term can contain more than one fitted column. For example, one term can
#' represent an FP2 function, a factor with several contrast columns, a
#' zero-handled variable, or several matrix columns grouped with `term_groups`.
#' These columns are combined and returned as one term result.
#'
#' Use the original variable or group name in `terms`. Use
#' [get_selected_variable_names()] to see the terms included in the final model.
#' A requested term that is not in the final model is omitted with a warning.
#'
#' For a single numeric term:
#'
#' - `terms_seq = "equidistant"` uses `nseq` equally spaced values between the
#'   smallest and largest available values.
#' - `terms_seq = "data"` uses the values in `newdata`, or the fitted values
#'   when `newdata = NULL`.
#'
#' Factors and other grouped terms are evaluated row by row.
#'
#' For Gaussian, binomial, Poisson, and negative-binomial models, `add_intercept = TRUE` includes
#' the fitted intercept in each term result. Set `add_intercept = FALSE` to
#' return only the contribution from the requested term. This argument does not
#' affect Cox models or contrasts.
#'
#' @section Contrasts and reference values:
#' With `type = "contrasts"`, the returned value is:
#'
#' `fitted term value at the evaluated value - fitted term value at the reference`.
#'
#' Supply references as a named list, for example:
#'
#' `ref = list(age = 50, race = "A")`.
#'
#' Supply numeric references on the original variable scale. If a variable was
#' shifted during fitting, the same shift is added to the reference automatically.
#' The fitted FP or ACD function and centering are then applied automatically.
#' Do not shift, transform, or center the reference yourself.
#'
#' For example, if the fitted function uses `log(age + 10)`, supply
#' `ref = list(age = 50)`, not 60 and not `log(60)`.
#'
#' For a formula factor, supply a fitted level label. For another grouped term,
#' such as a formula interaction or a matrix term created with `term_groups`,
#' supply one numeric value for every member column. A named vector is
#' recommended, for example:
#'
#' `ref = list(race = c(raceB = 0, raceC = 0))`.
#' The member names are the raw model-matrix columns used for that grouped term.
#'
#' If no reference is supplied, the function uses:
#'
#' - the smaller value for a numeric term with two distinct values;
#' - the mean for another single numeric term;
#' - the first fitted level for a formula factor;
#' - the first available row for another grouped term.
#'
#' The numeric and grouped-row defaults are based on `newdata` when it is
#' supplied and on the fitted data otherwise. A formula factor always uses its
#' first fitted level.
#'
#' The model intercept is not included in contrasts because it cancels from the
#' comparison.
#'
#' @section Zero and spike-at-zero variables:
#' Supply zero-handled and spike-at-zero variables on their original scale.
#' Do not create a separate zero-indicator column.
#'
#' The zero handling selected during fitting is reused during prediction. A
#' nonpositive value is treated in the same way as a nonpositive value in the
#' fitted model, and a positive value is passed to the selected continuous
#' function when applicable.
#'
#' If the final spike-at-zero model includes both a continuous component and a
#' zero indicator, both are used automatically. If it includes only one
#' component, only that component is used.
#'
#' The same rule is applied to a supplied reference. For a term that includes
#' only the zero indicator, a nonpositive reference is treated as indicator 1
#' and a positive reference as indicator 0.
#'
#' @section Offsets:
#' Offsets are used only for complete-model predictions. They are not included
#' in term or contrast results.
#'
#' When `newdata = NULL`, the offset values used for fitting are used.
#'
#' When `newdata` is supplied:
#'
#' - If the model formula contains `offset(...)`, include the variable or
#'   variables needed by that expression in `newdata`. The offset is evaluated
#'   automatically.
#' - If the offset was supplied through the separate `offset` argument, supply
#'   the offset for the new rows through `newoffset`.
#' - A supplied `newoffset` takes precedence over a formula offset.
#'
#' Supply one finite `newoffset` value for each row of `newdata`.
#'
#' @section Cox predictions:
#' Cox predictions with `type = "lp"` and `type = "risk"` are relative
#' quantities. They are not absolute hazards or survival probabilities.
#'
#' Use `cox_reference` to choose the covariate reference for these relative
#' predictions:
#'
#' - `"zero"` uses zero on the fitted predictor scale and is the default.
#' - `"sample"` uses the fitted-sample predictor means.
#' - `"strata"` uses fitted-sample means within each stratum.
#'
#' Changing `cox_reference` changes the displayed linear predictors and risks,
#' but comparisons made using the same reference remain unchanged. With
#' `cox_reference = "strata"`, comparisons should be made within the same
#' stratum.
#'
#' For the right-censored Cox models supported by `mfp2`,
#' `type = "expected"` or `type = "survival"` requires `newdata` to contain
#' follow-up information. For a formula model, include the variables used in the
#' fitted `Surv()` response. Alternatively, include exactly one column that is a
#' `Surv` object. Wrap it in `I()` when placing it in a data frame, for example:
#'
#' `data.frame(x = xnew, followup = I(survival::Surv(time, event)))`.
#'
#' When `newdata = NULL`, the follow-up values used to fit the model are used.
#'
#' Do not use `cox_reference` with `type = "expected"` or
#' `type = "survival"`.
#'
#' For a stratified Cox model, each prediction row must have a valid stratum.
#' If `strata()` was included in the formula, include its original variable or
#' variables in `newdata`. Otherwise, supply the new strata through `strata`.
#'
#' @param object A fitted object of class `"mfp2"`.
#'
#' @param newdata An optional data frame or matrix containing observations to
#'   predict. See **Supplying prediction data**.
#'
#' @param type A character value selecting the prediction type. See
#'   **Choosing the prediction type**.
#'
#' @param se.fit A single `TRUE` or `FALSE` value. For complete-model
#'   predictions, `TRUE` requests standard errors from the underlying model
#'   method. Term and contrast results always include standard errors and
#'   confidence limits, so this argument is ignored for those prediction types.
#'
#' @param terms An optional character vector naming terms to return with
#'   `type = "terms"` or `type = "contrasts"`. Use original variable or group
#'   names. If `NULL`, all terms included in the final model are used.
#'
#' @param terms_seq A character value controlling the evaluation values for a
#'   single numeric term. Use `"equidistant"` or `"data"`. The default is
#'   `"equidistant"`. Grouped terms are always evaluated row by row.
#'
#' @param alpha A number between 0 and 1 used to calculate confidence intervals
#'   for term and contrast results. The confidence level is `1 - alpha`. The
#'   default `alpha = 0.05` gives 95 percent confidence intervals.
#'
#' @param ref A named list of reference values used with
#'   `type = "contrasts"`. Numeric references must be supplied on the original
#'   variable scale. See **Contrasts and reference values**.
#'
#' @param strata Optional stratum information for complete-model prediction
#'   from a stratified Cox model when `newdata` is supplied. For one
#'   stratification variable, supply one value per row. For several
#'   stratification variables, supply a
#'   matrix or data frame with one row per prediction row.
#'
#'   This argument can be supplied only for complete-model prediction from a
#'   stratified Cox model when `newdata` is supplied.
#'
#' @param newoffset An optional finite numeric vector containing one offset
#'   value for each row of `newdata`. See **Offsets**.
#'
#'   This argument can be supplied only for a complete-model prediction when the
#'   fitted model used an offset and `newdata` is supplied.
#'
#' @param nseq A positive integer giving the number of equally spaced values
#'   used when `terms_seq = "equidistant"`. The default is 100.
#'
#' @param add_intercept A single `TRUE` or `FALSE` value. For term predictions
#'   from Gaussian, binomial, Poisson, and negative-binomial models, `TRUE` includes the fitted
#'   intercept in each result. The default is `TRUE`. It has no effect for Cox
#'   models or contrasts.
#'
#' @param cox_reference A character value choosing the covariate reference for
#'   Cox predictions with `type = "lp"` or `type = "risk"`. Use `"zero"`,
#'   `"sample"`, or `"strata"`. The default is `"zero"`. It cannot be
#'   used for non-Cox models, term or contrast predictions, or Cox predictions
#'   with `type = "expected"` or `type = "survival"`.
#'
#' @param ... Further arguments passed to [stats::predict.glm()] or
#'   [survival::predict.coxph()] for complete-model predictions.
#'
#' @return
#' For complete-model predictions, the result follows [stats::predict.glm()] or
#' [survival::predict.coxph()]. It is usually a numeric vector. When
#' `se.fit = TRUE`, it is usually a list containing fitted values and standard
#' errors.
#'
#' For `type = "terms"` or `type = "contrasts"`, the result is a named list
#' with one data frame for each requested term. Each data frame contains:
#'
#' - `variable`: the evaluated value or factor level;
#' - `variable_pre`: the value after the fitted shift and before the final
#'   transformation for a single numeric term. For an indicator-only spike
#'   term, this is the zero indicator;
#' - `value`: the term-specific fitted value or contrast;
#' - `se`: the standard error;
#' - `lower`: the lower confidence limit;
#' - `upper`: the upper confidence limit.
#'
#' Grouped terms can also include their member columns so that each evaluated
#' row can be identified.
#'
#' @examples
#' data("prostate")
#'
#' fit <- mfp2(
#'   lpsa ~ fp(age) + fp(cavol) + fp(weight) + svi,
#'   data = prostate,
#'   keep = c("age", "cavol"),
#'   verbose = FALSE
#' )
#'
#' # Predict for the data used to fit the model.
#' predict(fit)
#'
#' # Predict for five observations on the response scale.
#' predict(fit, newdata = prostate[1:5, ], type = "response")
#'
#' # Link-scale predictions with standard errors.
#' predict(fit, newdata = prostate[1:5, ], se.fit = TRUE)
#'
#' # Fitted contributions of age and cancer volume.
#' predict(
#'   fit,
#'   type = "terms",
#'   terms = c("age", "cavol"),
#'   add_intercept = FALSE
#' )
#'
#' # Compare age values with age 50 and cancer volume values with 1.
#' # The references are supplied on the original variable scale.
#' predict(
#'   fit,
#'   type = "contrasts",
#'   terms = c("age", "cavol"),
#'   ref = list(age = 50, cavol = 1)
#' )
#'
#' \dontrun{
#' # A formula offset is evaluated from newdata.
#' set.seed(1)
#' d <- data.frame(
#'   y = rpois(80, 2),
#'   age = runif(80, 30, 75),
#'   exposure = runif(80, 0.5, 3)
#' )
#'
#' fit_offset <- mfp2(
#'   y ~ fp(age) + offset(log(exposure)),
#'   data = d,
#'   family = "poisson",
#'   verbose = FALSE
#' )
#'
#' predict(
#'   fit_offset,
#'   newdata = d[1:4, c("age", "exposure")],
#'   type = "response"
#' )
#'
#' # An offset supplied separately must be supplied again for new rows.
#' fit_external_offset <- mfp2(
#'   y ~ fp(age),
#'   data = d,
#'   family = "poisson",
#'   offset = log(d$exposure),
#'   verbose = FALSE
#' )
#'
#' predict(
#'   fit_external_offset,
#'   newdata = d[1:4, "age", drop = FALSE],
#'   newoffset = log(d$exposure[1:4]),
#'   type = "response"
#' )
#'
#' # Spike-at-zero references are supplied on the original scale.
#' d_spike <- data.frame(
#'   y = rnorm(120),
#'   exposure = c(rep(0, 35), runif(85, 0.1, 8))
#' )
#'
#' fit_spike <- mfp2(
#'   y ~ fp(exposure, spike = TRUE),
#'   data = d_spike,
#'   keep = "exposure",
#'   verbose = FALSE
#' )
#'
#' predict(
#'   fit_spike,
#'   newdata = data.frame(exposure = c(0, 1, 4)),
#'   type = "contrasts",
#'   terms = "exposure",
#'   ref = list(exposure = 0)
#' )
#'
#' # Cox predictions.
#' data("gbsg")
#'
#' fit_cox <- mfp2(
#'   survival::Surv(rectime, censrec) ~ fp(age) + fp(nodes),
#'   data = gbsg,
#'   family = "cox",
#'   df = 1,
#'   select = 1,
#'   alpha = 1,
#'   verbose = FALSE
#' )
#'
#' relative_profiles <- gbsg[1:3, c("age", "nodes")]
#' predict(fit_cox, relative_profiles, type = "lp")
#' predict(fit_cox, relative_profiles, type = "risk")
#'
#' # Survival predictions also need the original follow-up variables.
#' absolute_profiles <- gbsg[
#'   1:3,
#'   c("rectime", "censrec", "age", "nodes")
#' ]
#' predict(fit_cox, absolute_profiles, type = "survival")
#' predict(fit_cox, absolute_profiles, type = "expected")
#' }
#'
#' @seealso
#' [mfp2()], [get_selected_variable_names()], [stats::predict.glm()],
#' [survival::predict.coxph()]
#'
#' @method predict mfp2
#' @export
predict.mfp2 <- function(object,
                         newdata = NULL,
                         type = NULL,
                         se.fit = FALSE,
                         terms = NULL,
                         terms_seq = c("equidistant", "data"),
                         alpha = 0.05,
                         ref = NULL,
                         strata = NULL,
                         newoffset = NULL,
                         nseq = 100,
                         add_intercept = TRUE,
                         cox_reference = NULL,
                         ...) {

  # Record whether the caller supplied a non-NULL Cox covariate reference.
  # NULL follows predict.mfpi() semantics and means "use the method default".
  # The explicit-supply distinction is essential for expected/survival
  # predictions, where the underlying Cox method fixes its own sample reference.
  cox_reference_supplied <- !missing(cox_reference) && !is.null(cox_reference)

  # `newy` was considered during development but is intentionally not part of
  # the public API. Absolute Cox prediction data belong in `newdata`, matching
  # the data-oriented convention used by predict.coxph(). Catch legacy draft
  # calls here so the argument cannot disappear silently into `...`.
  dots_call <- match.call(expand.dots = FALSE)[["..."]]
  dots_names <- names(as.list(dots_call))
  if (!is.null(dots_names) && "reference" %in% dots_names) {
    stop(
      "'reference' has been renamed to 'cox_reference' in predict.mfp2().",
      call. = FALSE
    )
  }
  if (!is.null(dots_names) && "newy" %in% dots_names) {
    stop(
      "'newy' has been removed. Supply the Cox prediction response in ",
      "'newdata': include the original formula response variables or one ",
      "Surv column.",
      call. = FALSE
    )
  }

  terms_seq <- match.arg(terms_seq)

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(nseq) || nseq <= 0) {
    stop("'nseq' must be a positive integer.", call. = FALSE)
  }

  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be a single TRUE or FALSE value.", call. = FALSE)
  }

  # Prediction metadata and fitted-model methods are defined only for mfp2 objects.
  if (!inherits(object, "mfp2")) {
    stop("The object is not an mfp2 object.", call. = FALSE)
  }

  # Resolve the prediction scale. Formula fits store the same conceptual term
  # names used by selection and printing, including source-variable names for
  # simple wrappers such as factor(x9).
  if (is.null(type)) {
    type <- ifelse(object$family_string == "cox", "lp", "link")
  }

  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single non-missing character value.", call. = FALSE)
  }

  # `link` and `lp` are family-specific names for the same linear-predictor
  # quantity. Normalize that unambiguous alias before method-specific type
  # matching so users receive the intended prediction rather than a downstream
  # match.arg() error from the base model method.
  type <- normalize_prediction_type(
    type = type,
    family_string = object$family_string
  )

  # Cox prediction types are resolved here rather than left to
  # predict.coxph(). This lets mfp2 distinguish its own term/contrast paths from
  # the relative and absolute full-model paths and report argument conflicts
  # before any transformed data are constructed.
  resolved_cox_reference <- NULL
  if (identical(object$family_string, "cox")) {
    type <- mfp2_match_cox_prediction_type(type)
    resolved_cox_reference <- mfp2_validate_cox_reference(
      type = type,
      cox_reference = cox_reference,
      cox_reference_supplied = cox_reference_supplied
    )
  } else {
    type <- mfp2_match_glm_prediction_type(type)
    if (cox_reference_supplied) {
      stop(
        "'cox_reference' is only available for Cox models.",
        call. = FALSE
      )
    }
  }

  # `strata` and `newoffset` alter only native full-model predictions on
  # supplied rows. Reject irrelevant combinations before formula evaluation or
  # transformation so that an accepted argument is never silently ignored.
  full_model_prediction <- !type %in% c("terms", "contrasts")

  if (!is.null(strata)) {
    if (!identical(object$family_string, "cox")) {
      stop("'strata' is available only for Cox predictions.", call. = FALSE)
    }
    if (!full_model_prediction) {
      stop(
        "'strata' is not used for mfp2 term or contrast predictions.",
        call. = FALSE
      )
    }
    if (is.null(newdata)) {
      stop(
        "'strata' can be supplied only together with 'newdata'.",
        call. = FALSE
      )
    }
    if (!prediction_fit_has_strata(object)) {
      stop(
        "'strata' was supplied, but the retained Cox model is not stratified.",
        call. = FALSE
      )
    }
  }

  if (!is.null(newoffset)) {
    if (!full_model_prediction) {
      stop(
        "'newoffset' is not used for mfp2 term or contrast predictions.",
        call. = FALSE
      )
    }
    if (is.null(newdata)) {
      stop(
        "'newoffset' can be supplied only together with 'newdata'.",
        call. = FALSE
      )
    }
    if (!isTRUE(object$has_offset)) {
      stop(
        "'newoffset' can be supplied only when the fitted model used an offset.",
        call. = FALSE
      )
    }
    if (!is.numeric(newoffset) || anyNA(newoffset) ||
        any(!is.finite(newoffset))) {
      stop(
        "'newoffset' must be a finite numeric vector with one value per prediction row.",
        call. = FALSE
      )
    }
    if (length(newoffset) != NROW(newdata)) {
      stop(
        "The length of 'newoffset' must equal the number of rows in 'newdata'.",
        call. = FALSE
      )
    }
  }

  selected_internal_terms <- get_selected_variable_names(object)

  if (is.null(terms)) {
    terms <- selected_internal_terms
  }

  if (is.null(ref)) {
    ref <- stats::setNames(lapply(terms, function(v) NULL), terms)
  }

  # Preserve the user's original formula-style data for factor labels, offsets,
  # and strata. Formula fits are then rebuilt on their fitted model-matrix scale.
  if (!is.null(newdata) && is.null(colnames(newdata))) {
    stop("Newdata must have column names", call. = FALSE)
  }

  newdata_raw <- newdata

  # Reject infinite numeric values before formula reconstruction can discard
  # unselected columns or model.matrix() can propagate them into derived terms.
  # Missing-value validation remains limited to the reconstructed predictors
  # actually required by the fitted model.
  if (!is.null(newdata_raw)) {
    mfp2_validate_prediction_newdata(
      newdata = newdata_raw,
      check_missing = FALSE
    )
  }

  # Absolute Cox predictions with newdata need a Surv response because
  # predict.coxph() evaluates the cumulative baseline hazard at the supplied
  # time or interval. Reconstruct it before formula newdata are reduced to the
  # active predictor design columns.
  cox_prediction_response <- if (
    !is.null(newdata_raw) &&
    identical(object$family_string, "cox") &&
    type %in% c("expected", "survival")
  ) {
    reconstruct_cox_prediction_response(
      object = object,
      fit_obj = object,
      newdata = newdata_raw
    )
  } else {
    NULL
  }

  prediction_terms <- if (type %in% c("terms", "contrasts")) {
    intersect(terms, selected_internal_terms)
  } else {
    selected_internal_terms
  }

  if (!is.null(newdata)) {
    newdata <- reconstruct_formula_newdata(
      object,
      newdata,
      terms = prediction_terms
    )
  }

  if (!is.null(newdata)) {
    # Validate only the reconstructed predictors used by the requested
    # prediction. Matrix-interface callers may include irrelevant extra columns
    # containing missing values; those columns are intentionally ignored.
    validation_lookup <- prediction_term_to_columns(
      object,
      terms = prediction_terms
    )
    validation_columns <- unique(unlist(validation_lookup, use.names = FALSE))
    validation_columns <- intersect(validation_columns, colnames(newdata))
    mfp2_validate_prediction_newdata(
      newdata[, validation_columns, drop = FALSE]
    )
  }

  if (type == "contrasts" && length(ref) != sum(names(ref) != "", na.rm = TRUE)) {
    warning(
      paste0(
        "i The supplied reference values (ref) must all be named.\n",
        "i predict() continues with the default reference: the mean for ",
        "continuous terms, the minimum for binary terms, or the first fitted ",
        "level for formula factors."
      ),
      call. = FALSE
    )
  }

  # Term and contrast prediction is handled term by term. A conceptual term
  # may correspond to one transformed continuous predictor or to a block of
  # dummy/contrast columns for a categorical term.
  if (type %in% c("terms", "contrasts")) {
    requested_terms <- terms
    internal_terms <- requested_terms

    selected <- internal_terms %in% selected_internal_terms
    if (!any(selected)) {
      warning(
        "i All the terms supplied are not in the final model.\n",
        "i predict() continues but returns an empty list.",
        call. = FALSE
      )
    } else if (!all(selected)) {
      warning(
        paste0(
          "i Some terms supplied are not in the final model.\n",
          "i predict() continues but omits terms not in the model."
        ),
        call. = FALSE
      )
    }

    requested_terms <- requested_terms[selected]
    internal_terms <- internal_terms[selected]

    valid_reference_names <- requested_terms
    if (!all(vapply(ref, is.null, logical(1L))) &&
        any(!names(ref) %in% valid_reference_names)) {
      warning(
        "i Some reference names are not selected prediction terms and are ignored.",
        call. = FALSE
      )
    }

    # The intercept contributes to GLM term predictions only when requested.
    # It cancels from contrasts and is not present in Cox partial predictors.
    cf <- coef(object)
    intercept <- if ("(Intercept)" %in% names(cf)) cf[["(Intercept)"]] else 0
    if (is.na(intercept) || !add_intercept || type == "contrasts") {
      intercept <- 0
    }

    # Resolve conceptual prediction terms once. Each element names the raw
    # design-matrix columns that must be transformed and combined for that term.
    lookup <- prediction_term_to_columns(object)
    res_list <- list()

    for (term_index in seq_along(internal_terms)) {
      t <- internal_terms[[term_index]]
      display_term <- requested_terms[[term_index]]
      term_columns <- lookup[[t]]
      factor_info <- if (!is.null(object$formula_factor_info)) {
        object$formula_factor_info[[t]]
      } else {
        NULL
      }
      block_term <- term_uses_column_mapping(t, term_columns) || !is.null(factor_info)

      if (block_term) {
        # Grouped terms are predicted as complete design blocks. This includes
        # multi-column factors and one-column formula factors carrying factor
        # metadata for user-facing level labels.
        # Reconstructed `newdata` is on the raw model-matrix scale and still
        # needs the fitted shift. In contrast, `object$x_original` is stored on
        # the shifted-but-not-scaled final-fitting scale, so its shift must not
        # be applied a second time.
        source_needs_preprocessing <- !is.null(newdata)

        source_block <- if (!is.null(newdata)) {
          missing_columns <- setdiff(term_columns, colnames(newdata))
          if (length(missing_columns) > 0L) {
            stop(
              "Prediction data are missing grouped-term column(s): ",
              paste(missing_columns, collapse = ", "),
              call. = FALSE
            )
          }
          newdata[, term_columns, drop = FALSE]
        } else {
          missing_columns <- setdiff(term_columns, colnames(object$x_original))
          if (length(missing_columns) > 0L) {
            stop(
              "Fitted data do not contain grouped-term column(s): ",
              paste(missing_columns, collapse = ", "),
              call. = FALSE
            )
          }
          object$x_original[, term_columns, drop = FALSE]
        }

        source_block <- as.matrix(source_block)
        x_trafo <- as.matrix(prepare_newdata_for_predict(
          object,
          source_block,
          terms = t,
          apply_pre = source_needs_preprocessing,
          allow_missing_predictors = TRUE,
          check_binary = FALSE
        ))

        model_columns <- prediction_model_column_names(
          object,
          colnames(x_trafo)
        )
        term_coef <- stats::coef(object)[model_columns]
        term_coef[is.na(term_coef)] <- 0

        labels <- decode_grouped_term_values(
          object,
          term = t,
          block = source_block
        )

        res <- data.frame(
          variable = labels,
          variable_pre = labels,
          source_block,
          value = as.numeric(x_trafo %*% term_coef + intercept),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )

        x_ref_trafo <- NULL
        if (type == "contrasts") {
          reference_block <- resolve_grouped_reference(
            object,
            term = t,
            reference = if (display_term %in% names(ref)) ref[[display_term]] else ref[[t]],
            block = source_block
          )
          x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            reference_block,
            terms = t,
            apply_pre = TRUE,
            allow_missing_predictors = TRUE,
            check_binary = FALSE
          ))
          res$value <- res$value - as.numeric(x_ref_trafo %*% term_coef)
        }
      } else {
        # Singleton numeric terms retain the ordinary FP prediction path.
        # Evaluation values are represented on the shifted, unscaled scale used
        # by the final fitted transformation.
        if (terms_seq == "equidistant") {
          if (!is.null(newdata)) {
            x_range <- range(newdata[, t]) + object$transformations[t, "shift"]
          } else {
            x_range <- range(object$x_original[, t])
          }

          x_seq <- matrix(
            seq(x_range[1L], x_range[2L], length.out = nseq),
            ncol = 1L,
            dimnames = list(NULL, t)
          )

          x_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_seq,
            terms = t,
            apply_pre = FALSE,
            allow_missing_predictors = TRUE
          ))
        } else {
          x_seq <- if (!is.null(newdata)) {
            matrix(
              newdata[[t]] + object$transformations[t, "shift"],
              ncol = 1L,
              dimnames = list(NULL, t)
            )
          } else {
            object$x_original[, t, drop = FALSE]
          }

          x_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_seq,
            terms = t,
            apply_pre = FALSE,
            allow_missing_predictors = TRUE
          ))
        }

        model_columns <- prediction_model_column_names(
          object,
          colnames(x_trafo)
        )
        term_coef <- stats::coef(object)[model_columns]
        term_coef[is.na(term_coef)] <- 0

        variable <- as.numeric(x_seq) - object$transformations[t, "shift"]
        variable_pre <- as.numeric(x_seq)

        if (object$spike_dec[t] == saz_decision_codes[["binary_only"]]) {
          variable <- as.integer(x_seq[, t] <= 0)
          variable_pre <- variable
        }

        res <- data.frame(
          variable = variable,
          variable_pre = variable_pre,
          value = as.numeric(x_trafo %*% term_coef + intercept)
        )

        x_ref_trafo <- NULL
        if (type == "contrasts") {
          x_ref <- if (display_term %in% names(ref)) ref[[display_term]] else ref[[t]]
          if (is.null(x_ref)) {
            v <- if (!is.null(newdata)) {
              newdata[, t] + object$transformations[t, "shift"]
            } else {
              object$x_original[, t]
            }
            x_ref <- if (length(unique(stats::na.omit(v))) == 2L) {
              min(v, na.rm = TRUE)
            } else {
              mean(v, na.rm = TRUE)
            }
          } else {
            x_ref <- x_ref + object$transformations[t, "shift"]
          }

          x_ref <- matrix(x_ref, nrow = 1L, ncol = 1L, dimnames = list(NULL, t))
          x_ref_trafo <- as.matrix(prepare_newdata_for_predict(
            object,
            x_ref,
            terms = t,
            apply_pre = FALSE,
            check_binary = FALSE,
            reset_zero = FALSE,
            allow_missing_predictors = TRUE
          ))
          res$value <- res$value - as.numeric(x_ref_trafo %*% term_coef)
        }
      }

      # Standard errors use the coefficient covariance block corresponding to
      # the transformed columns for this conceptual term.
      res$se <- calculate_standard_error(
        object,
        x_trafo,
        x_ref_trafo,
        include_intercept = add_intercept && type == "terms"
      )
      mult <- stats::qnorm(1 - alpha / 2)
      res$lower <- res$value - mult * res$se
      res$upper <- res$value + mult * res$se
      res_list[[display_term]] <- res
    }

    return(res_list)
  }

  # Full-model prediction reconstructs every required raw design column,
  # applies the fitted transformations, and then delegates to the underlying
  # glm or coxph prediction method.

  if (!is.null(newdata)) {

    if (
      is.null(newoffset) &&
      !is.null(newdata_raw) &&
      !is.null(object$formula_offset_terms)
    ) {
      newoffset <- reconstruct_formula_offset_newdata(object, newdata_raw)
    }

    # Formula offsets are reconstructed automatically when possible; otherwise
    # an offset must be supplied explicitly for models fitted with one.
    has_offset <- isTRUE(object$has_offset)
    if (has_offset) {
      if (is.null(newoffset)) {
        stop(
          "No newoffset provided for prediction, yet offset was used in mfp2",
          call. = FALSE
        )
      }

      if (!is.numeric(newoffset) || anyNA(newoffset) ||
          any(!is.finite(newoffset))) {
        stop(
          "! newoffset must be a finite numeric vector.",
          call. = FALSE
        )
      }

      if (length(newoffset) != nrow(newdata)) {
        stop(
          "The length of newoffset must be equal to the number of rows of newdata",
          call. = FALSE
        )
      }

    } else {
      newoffset <- NULL
    }

    if (
      object$family_string == "cox" &&
      is.null(strata) &&
      !is.null(newdata_raw) &&
      !is.null(object$formula_strata_terms)
    ) {
      strata <- reconstruct_formula_strata_newdata(object, newdata_raw)
    }

    newdata <- prepare_newdata_for_predict(
      object,
      newdata,
      terms = selected_internal_terms,
      strata = strata,
      offset = newoffset,
      check_binary = FALSE
    )

    if (!is.null(cox_prediction_response)) {
      newdata <- attach_cox_prediction_response(
        fit_obj = object,
        newdata = newdata,
        response = cox_prediction_response
      )
    }

    if (identical(object$family_string, "cox")) {
      return(
        mfp2_predict_cox_base(
          object = object,
          newdata = newdata,
          type = type,
          se.fit = se.fit,
          cox_reference = resolved_cox_reference,
          ...
        )
      )
    }

    if (inherits(object, "fastglm")) {
      prediction_data <- as.data.frame(newdata, check.names = FALSE)
      prediction_offset <- if ("offset_" %in% names(prediction_data)) {
        prediction_data$offset_
      } else {
        rep(0, nrow(prediction_data))
      }
      prediction_data$offset_ <- NULL
      prediction_matrix <- as.matrix(prediction_data)
      model_columns <- prediction_model_column_names(
        object,
        colnames(prediction_matrix)
      )
      colnames(prediction_matrix) <- model_columns

      return(
        mfp2_predict_fastglm_matrix(
          object = object,
          newx = prediction_matrix,
          offset = prediction_offset,
          type = type,
          se.fit = se.fit,
          ...
        )
      )
    }

    # Strip "mfp2" so S3 dispatch reaches predict.glm() rather than recursively
    # entering predict.mfp2().
    obj_base <- object
    class(obj_base) <- setdiff(class(obj_base), "mfp2")

    return(
      stats::predict(
        obj_base,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        ...
      )
    )
  }

  # With no newdata, the fitted base object already contains the training
  # response, transformed design, strata, and offsets required by its native
  # prediction method.
  if (identical(object$family_string, "cox")) {
    return(
      mfp2_predict_cox_base(
        object = object,
        type = type,
        se.fit = se.fit,
        cox_reference = resolved_cox_reference,
        ...
      )
    )
  }

  if (inherits(object, "fastglm")) {
    return(
      mfp2_predict_fastglm_matrix(
        object = object,
        type = type,
        se.fit = se.fit,
        ...
      )
    )
  }

  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "mfp2")
  stats::predict(
    obj_base,
    type = type,
    se.fit = se.fit,
    ...
  )
}


#' Match a Cox Prediction Type Used by `predict.mfp2()`
#'
#' Resolves exact or unambiguous partial Cox prediction-type names before the
#' method chooses a prediction path. `mfp2` implements `"terms"` and
#' `"contrasts"` itself, delegates `"lp"` and `"risk"` to the relative
#' prediction branch of `predict.coxph()`, and delegates `"expected"` and
#' `"survival"` to its baseline-hazard branch.
#'
#' Performing this match centrally prevents a partially matched type from being
#' routed to the wrong branch and makes invalid values fail before prediction
#' data are transformed.
#'
#' @param type Character scalar supplied to `predict.mfp2()`.
#'
#' @return One of `"lp"`, `"risk"`, `"expected"`, `"survival"`,
#'   `"terms"`, or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_cox_prediction_type <- function(type) {
  choices <- c("lp", "risk", "expected", "survival", "terms", "contrasts")

  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For Cox models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  )
}


#' Match a GLM Prediction Type Used by `predict.mfp2()`
#'
#' Resolves exact or unambiguous partial GLM prediction-type names before the
#' method constructs transformed prediction data. Linear-predictor aliases have
#' already been normalized by `normalize_prediction_type()`.
#'
#' @param type Character scalar supplied to `predict.mfp2()` after alias
#'   normalization.
#'
#' @return One of `"link"`, `"response"`, `"terms"`, or `"contrasts"`.
#'
#' @keywords internal
#' @noRd
mfp2_match_glm_prediction_type <- function(type) {
  choices <- c("link", "response", "terms", "contrasts")

  if (identical(type, "risk")) {
    stop(
      "For GLM models, type = 'risk' is not supported. Use type = ",
      "'response' only when prediction on the model response scale is ",
      "intended.",
      call. = FALSE
    )
  }

  tryCatch(
    match.arg(type, choices),
    error = function(e) {
      stop(
        "For GLM models, 'type' must be one of: ",
        paste(shQuote(choices), collapse = ", "),
        ". The alias 'lp' is also accepted for 'link'.",
        call. = FALSE
      )
    }
  )
}


#' Validate Prediction Data
#'
#' Checks raw or reconstructed prediction data before ranges,
#' fractional-polynomial transformations, or model matrices are calculated.
#' Missing values retain their dedicated diagnostic when requested, while
#' positive and negative infinity are always rejected and reported by column.
#'
#' @param newdata Prediction data used by `predict.mfp2()`.
#' @param check_missing Logical scalar. When `TRUE`, reject missing values in
#'   addition to non-finite numeric values. Raw user data are checked with
#'   `FALSE` before formula reconstruction; reconstructed required predictors
#'   are checked with `TRUE` before transformation.
#'
#' @return Invisibly `TRUE`.
#'
#' @keywords internal
#' @noRd
mfp2_validate_prediction_newdata <- function(newdata, check_missing = TRUE) {
  if (!is.logical(check_missing) || length(check_missing) != 1L ||
      is.na(check_missing)) {
    stop("Internal error: `check_missing` must be TRUE or FALSE.",
         call. = FALSE)
  }

  newdata <- as.data.frame(newdata, check.names = FALSE)

  if (isTRUE(check_missing) && anyNA(newdata)) {
    stop(
      "! newdata must not contain any NA (missing data).\n",
      "i Please remove any missing data before passing newdata to this function.",
      call. = FALSE
    )
  }

  non_finite_columns <- names(newdata)[vapply(
    newdata,
    function(column) {
      if (!is.numeric(column)) {
        return(FALSE)
      }
      bad <- if (isTRUE(check_missing)) {
        !is.finite(column)
      } else {
        is.infinite(column)
      }
      any(bad)
    },
    logical(1L)
  )]

  if (length(non_finite_columns) > 0L) {
    stop(
      "! newdata must contain only finite numeric values.\n",
      "i Infinite value(s) were found in numeric column(s): ",
      paste(non_finite_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate the Cox Covariate Reference for the Requested Prediction Type
#'
#' `survival::predict.coxph()` uses `reference` only for relative
#' linear-predictor, risk, and native term predictions. `mfp2` computes its own
#' term and contrast output, while the Cox expected-event/survival path
#' unconditionally uses the sample reference internally. This helper enforces
#' those distinctions and returns a value only when the delegated Cox method can
#' use it.
#'
#' The explicit-supply flag is kept separate from the argument value so that
#' `NULL` means "use the method default", matching `predict.mfpi()`, while a
#' non-NULL value can be rejected on prediction paths where references do not
#' apply.
#'
#' @param type Resolved Cox prediction type.
#' @param cox_reference Candidate Cox reference argument.
#' @param cox_reference_supplied Logical scalar indicating whether the caller
#'   explicitly supplied `cox_reference` to `predict.mfp2()`.
#'
#' @return A single validated reference string for `type = "lp"` or
#'   `type = "risk"`; otherwise `NULL`.
#'
#' @keywords internal
#' @noRd
mfp2_validate_cox_reference <- function(type,
                                        cox_reference,
                                        cox_reference_supplied) {
  if (!is.logical(cox_reference_supplied) ||
      length(cox_reference_supplied) != 1L ||
      is.na(cox_reference_supplied)) {
    stop("Internal error: invalid cox_reference-supply flag.", call. = FALSE)
  }

  if (type %in% c("terms", "contrasts")) {
    if (cox_reference_supplied) {
      stop(
        "'cox_reference' is not used for mfp2 term or contrast predictions. ",
        "Use 'ref' to choose variable-specific contrast reference values.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  if (type %in% c("expected", "survival")) {
    if (cox_reference_supplied) {
      stop(
        "'cox_reference' does not apply to Cox predictions with type = '",
        type,
        "'. The baseline-hazard calculation uses the fitted sample ",
        "reference internally; remove the 'cox_reference' argument.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  # The public NULL default resolves to "zero" for relative Cox predictions.
  # A non-NULL user value must be scalar and is validated below.
  if (!cox_reference_supplied) {
    return("zero")
  }

  match_cox_reference(
    value = cox_reference,
    default = "zero",
    argument = "cox_reference"
  )
}


#' Delegate a Full-Model Prediction to the Stored Cox Fit
#'
#' Removes only the `"mfp2"` class and then calls the native Cox prediction
#' method through `stats::predict()`. The underlying object remains a complete
#' `coxph` fit with its fitted response, transformed design matrix, residuals,
#' means, covariance matrix, strata metadata, and tie method intact.
#'
#' Relative predictions (`"lp"` and `"risk"`) receive the validated covariate
#' reference. Absolute predictions (`"expected"` and `"survival"`) deliberately
#' omit that argument because `predict.coxph()` resets it to `"sample"` before
#' calculating the baseline cumulative hazard. Omitting it also avoids implying
#' that an ignored user value controls absolute predictions.
#'
#' When the fitted model contains an offset, the stored training-offset origin
#' is validated before delegation. Offset centering itself remains the
#' responsibility of `predict.coxph()` and is independent of the covariate
#' reference.
#'
#' @param object Fitted `mfp2` Cox object.
#' @param newdata Optional fully transformed prediction data frame. For absolute
#'   predictions it must also contain the internal `Surv` response column.
#' @param type Resolved Cox prediction type.
#' @param se.fit Logical scalar passed to `predict.coxph()`.
#' @param cox_reference Validated reference string for relative predictions, or
#'   `NULL` for absolute predictions.
#' @param ... Further arguments passed to the native Cox prediction method.
#'
#' @return The unmodified result returned by `predict.coxph()`.
#'
#' @keywords internal
#' @noRd
mfp2_predict_cox_base <- function(object,
                                  newdata = NULL,
                                  type,
                                  se.fit,
                                  cox_reference = NULL,
                                  ...) {
  if (!inherits(object, "mfp2") ||
      !identical(object$family_string, "cox")) {
    stop("Internal error: a fitted mfp2 Cox object is required.", call. = FALSE)
  }

  if (!type %in% c("lp", "risk", "expected", "survival")) {
    stop("Internal error: unsupported delegated Cox prediction type.",
         call. = FALSE)
  }

  if (isTRUE(object$has_offset)) {
    mfp2_cox_offset_reference(object)
  }

  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "mfp2")

  # Keep the native Cox method authoritative for both covariate-reference and
  # offset origins. In particular, do not add or subtract
  # object$cox_offset_reference here: predict.coxph() reconstructs the fitted
  # model frame when an offset is present, subtracts the mean training offset,
  # and applies the requested covariate reference in the same calculation.
  # Reimplementing either adjustment in mfp2 would double-center one component
  # and can diverge from the exact fitted coxph design.
  if (type %in% c("lp", "risk")) {
    if (is.null(cox_reference)) {
      stop("Internal error: relative Cox prediction lacks a reference.",
           call. = FALSE)
    }

    if (is.null(newdata)) {
      return(stats::predict(
        obj_base,
        type = type,
        se.fit = se.fit,
        reference = cox_reference,
        ...
      ))
    }

    return(stats::predict(
      obj_base,
      newdata = newdata,
      type = type,
      se.fit = se.fit,
      reference = cox_reference,
      ...
    ))
  }

  if (is.null(newdata)) {
    return(stats::predict(
      obj_base,
      type = type,
      se.fit = se.fit,
      ...
    ))
  }

  stats::predict(
    obj_base,
    newdata = newdata,
    type = type,
    se.fit = se.fit,
    ...
  )
}


#' Resolve Conceptual Terms to Raw Prediction Columns
#'
#' Returns the term lookup stored by grouped-term fits. The names are conceptual
#' terms and each value contains the corresponding raw design-matrix columns in
#' fitted term order. The lookup can be restricted to active prediction terms.
#'
#' @param object A fitted \code{mfp2} object.
#' @param terms Optional character vector of conceptual terms to retain.
#'
#' @return Named list mapping conceptual terms to raw input columns.
#'
#' @keywords internal
#' @noRd
prediction_term_to_columns <- function(object, terms = NULL) {
  fitted_terms <- rownames(object$transformations)

  if (is.null(fitted_terms) || length(fitted_terms) == 0L) {
    stop("The fitted object does not contain transformation term names.",
         call. = FALSE)
  }

  lookup <- object$term_to_columns

  if (is.null(lookup)) {
    # Objects created before grouped categorical support did not store an
    # explicit conceptual-term lookup. Formula fits from newer releases may
    # still carry the equivalent formula-level lookup, so prefer it when
    # available. Otherwise a continuous-only legacy object has an identity
    # relationship between conceptual terms and raw predictor columns.
    lookup <- object$formula_term_to_columns

    if (is.null(lookup)) {
      selected_terms <- tryCatch(
        get_selected_variable_names(object),
        error = function(e) character(0L)
      )
      selected_raw_columns <- if (is.null(object$x_original)) {
        character(0L)
      } else {
        colnames(object$x_original)
      }

      # A retained grouped term cannot be reconstructed safely without its
      # stored member-column mapping. Fail explicitly rather than silently
      # treating the conceptual name as a raw design column.
      unresolved_selected <- setdiff(selected_terms, selected_raw_columns)
      if (length(unresolved_selected) > 0L) {
        stop(
          "The fitted object lacks grouped-term mapping metadata for selected term(s): ",
          paste(unresolved_selected, collapse = ", "),
          ". Refit the model with the current mfp2 version.",
          call. = FALSE
        )
      }

      lookup <- stats::setNames(as.list(fitted_terms), fitted_terms)
    }
  }

  if (!is.list(lookup) || is.null(names(lookup)) || anyDuplicated(names(lookup))) {
    stop("The fitted grouped-term lookup is malformed.", call. = FALSE)
  }

  missing_terms <- setdiff(fitted_terms, names(lookup))
  if (length(missing_terms) > 0L) {
    stop(
      "The fitted grouped-term lookup is missing term(s): ",
      paste(missing_terms, collapse = ", "),
      call. = FALSE
    )
  }

  lookup <- lookup[fitted_terms]

  invalid <- vapply(
    lookup,
    function(columns) {
      !is.character(columns) || length(columns) == 0L || anyNA(columns) ||
        any(!nzchar(columns))
    },
    logical(1L)
  )

  if (any(invalid) || anyDuplicated(unlist(lookup, use.names = FALSE))) {
    stop("The fitted grouped-term lookup contains invalid raw columns.",
         call. = FALSE)
  }

  if (!is.null(terms)) {
    terms <- unique(as.character(terms))
    unknown <- setdiff(terms, names(lookup))
    if (length(unknown) > 0L) {
      stop(
        "Unknown fitted prediction term(s): ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    lookup <- lookup[terms]
  }

  lookup
}

# Resolve transformed design columns to the exact coefficient names stored by
# the final formula-based fit. The map is created once in fit_mfp() from the
# transformed matrix and fitted coefficient order; prediction therefore does
# not need to quote, unquote, sanitize, or otherwise guess model column names.
prediction_model_column_names <- function(object, transformed_columns) {
  column_map <- object$transformed_to_model_columns

  if (is.null(column_map)) {
    # Legacy continuous-only objects predate the explicit map. Recreate the
    # complete transformed source-column order from the stored training data,
    # then align it positionally with the fitted coefficient vector. Final
    # glm()/coxph() fits preserve formula-column order, including aliased
    # coefficients, so this recovers exact quoted or non-syntactic names.
    selected_terms <- get_selected_variable_names(object)

    if (length(selected_terms) == 0L) {
      if (length(transformed_columns) > 0L) {
        stop(
          "Cannot map transformed columns for an intercept-only fitted object.",
          call. = FALSE
        )
      }
      return(character(0L))
    }

    if (is.null(object$x_original)) {
      stop(
        "The legacy fitted object lacks training predictors needed to reconstruct coefficient-column metadata.",
        call. = FALSE
      )
    }

    reconstructed <- prepare_newdata_for_predict(
      object,
      object$x_original,
      terms = selected_terms,
      apply_pre = FALSE,
      allow_missing_predictors = FALSE
    )
    source_columns <- colnames(reconstructed)
    if (is.null(source_columns)) {
      source_columns <- character(0L)
    }

    fitted_columns <- names(stats::coef(object))
    if (is.null(fitted_columns)) {
      fitted_columns <- character(0L)
    }
    fitted_columns <- fitted_columns[fitted_columns != "(Intercept)"]

    if (length(source_columns) != length(fitted_columns)) {
      stop(
        "The legacy fitted object does not contain enough metadata to align transformed columns with fitted coefficients.",
        call. = FALSE
      )
    }

    column_map <- stats::setNames(fitted_columns, source_columns)
  }

  if (!is.character(column_map) || is.null(names(column_map)) ||
      anyDuplicated(names(column_map))) {
    stop(
      "The fitted transformed-to-model column mapping is malformed.",
      call. = FALSE
    )
  }

  model_columns <- unname(column_map[transformed_columns])
  if (anyNA(model_columns)) {
    stop(
      "The fitted object lacks coefficient mappings for transformed column(s): ",
      paste(transformed_columns[is.na(model_columns)], collapse = ", "),
      call. = FALSE
    )
  }

  model_columns
}



#' Expand Term-Level Metadata for Raw Prediction Columns
#'
#' Grouped fits store selection and transformation metadata once per conceptual
#' term, while model fitting and prediction use the underlying raw design
#' columns. This helper recreates the raw-column metadata needed by
#' \code{transform_matrix()}.
#'
#' @param object A fitted \code{mfp2} object.
#' @param raw_columns Character vector of raw prediction columns.
#'
#' @return List of raw-column metadata aligned with \code{raw_columns}.
#'
#' @keywords internal
#' @noRd
expand_prediction_metadata <- function(object, raw_columns) {
  lookup <- prediction_term_to_columns(object)
  spike_source <- if (!is.null(object$spike)) {
    object$spike
  } else {
    object$fp_terms[, "spike"]
  }

  expanded <- expand_term_metadata_to_columns(
    term_to_columns = lookup,
    powers = object$fp_powers,
    raw_columns = raw_columns,
    acdx = object$acd,
    zero = object$zero,
    catzero = object$catzero,
    spike = spike_source,
    spike_decision = object$spike_dec,
    acd_parameter = object$acd_parameter
  )

  raw_terms <- unname(expanded$terms)
  shifts <- stats::setNames(
    vapply(
      raw_terms,
      function(term) as.numeric(object$transformations[term, "shift"]),
      numeric(1L)
    ),
    raw_columns
  )

  list(
    terms = expanded$terms,
    powers = expanded$powers,
    shift = shifts,
    acdx = expanded$acdx,
    zero = expanded$zero,
    catzero = expanded$catzero,
    spike = expanded$spike,
    spike_decision = expanded$spike_decision,
    acd_parameter = expanded$acd_parameter
  )
}


#' Decode Formula Factor Levels from a Raw Design Block
#'
#' @param object A fitted \code{mfp2} object.
#' @param term Conceptual factor term name.
#' @param block Raw design block for the term.
#' @return Character vector of level labels or design signatures.
#'
#' @keywords internal
#' @noRd
decode_grouped_term_values <- function(object, term, block) {
  # Match each reconstructed design row to the fitted level-to-design mapping.
  # This remains correct when an inline factor() call changes labels, level
  # order, or contrast coding relative to the raw source variable.
  info <- if (!is.null(object$formula_factor_info)) {
    object$formula_factor_info[[term]]
  } else {
    NULL
  }

  if (!is.null(info) && !is.null(info$design_by_level)) {
    design <- as.matrix(info$design_by_level)
    out <- rep(NA_character_, nrow(block))

    for (i in seq_len(nrow(block))) {
      matches <- which(vapply(
        seq_len(nrow(design)),
        function(j) {
          all(is.finite(design[j, ])) &&
            isTRUE(all.equal(
              as.numeric(block[i, ]),
              as.numeric(design[j, ]),
              tolerance = 1e-10,
              check.attributes = FALSE
            ))
        },
        logical(1L)
      ))

      if (length(matches) > 0L) {
        out[[i]] <- rownames(design)[matches[[1L]]]
      }
    }

    if (all(!is.na(out))) {
      return(out)
    }
  }

  # Non-factor grouped terms have no level labels; expose their design row as
  # a compact signature so the returned value remains identifiable.
  apply(
    block,
    1L,
    function(row) paste(format(row, trim = TRUE), collapse = ",")
  )
}


#' Resolve a Reference Design Row for a Grouped Term
#'
#' @param object A fitted \code{mfp2} object.
#' @param term Conceptual term name.
#' @param reference User-supplied reference value or \code{NULL}.
#' @param block Observed raw design block used as a fallback reference source.
#'
#' @return One-row numeric matrix with columns matching \code{block}.
#'
#' @keywords internal
#' @noRd
resolve_grouped_reference <- function(object, term, reference, block) {
  # Formula factors accept level labels. Manually grouped matrix terms accept a
  # complete numeric design row with one value per raw member column.
  columns <- colnames(block)
  info <- if (!is.null(object$formula_factor_info)) {
    object$formula_factor_info[[term]]
  } else {
    NULL
  }

  if (!is.null(info)) {
    level <- if (is.null(reference)) info$levels[[1L]] else as.character(reference)

    if (length(level) == 1L && level %in% rownames(info$design_by_level)) {
      values <- info$design_by_level[level, columns, drop = FALSE]
      if (anyNA(values)) {
        stop(
          "No fitted design row is available for reference level '", level,
          "' of term '", term, "'.",
          call. = FALSE
        )
      }
      return(values)
    }

    if (!is.numeric(reference)) {
      stop(
        "Reference for factor term '", term, "' must be one of: ",
        paste(info$levels, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # For manually grouped terms, default to the first observed complete design
  # row when no explicit reference is supplied.
  if (is.null(reference)) {
    reference <- block[1L, , drop = FALSE]
  }

  if (is.data.frame(reference) || is.matrix(reference)) {
    reference <- as.matrix(reference)
    if (nrow(reference) != 1L) {
      stop("Grouped-term reference must contain exactly one row.", call. = FALSE)
    }
    if (!is.null(colnames(reference))) {
      missing <- setdiff(columns, colnames(reference))
      if (length(missing) > 0L) {
        stop(
          "Grouped-term reference is missing column(s): ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      reference <- reference[, columns, drop = FALSE]
    }
  } else {
    reference_names <- names(reference)
    reference <- as.numeric(reference)

    if (!is.null(reference_names) && all(nzchar(reference_names))) {
      missing <- setdiff(columns, reference_names)
      if (length(missing) > 0L) {
        stop(
          "Grouped-term reference is missing column(s): ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      names(reference) <- reference_names
      reference <- reference[columns]
    }

    if (length(reference) != length(columns)) {
      stop(
        "Reference for grouped term '", term, "' must contain ",
        length(columns), " value(s), one per raw column.",
        call. = FALSE
      )
    }
    reference <- matrix(reference, nrow = 1L, dimnames = list(NULL, columns))
  }

  storage.mode(reference) <- "double"
  reference
}


#' Transform Linear Predictor to Response or Risk
#'
#' Converts linear predictors (`nfit`) from a model to the appropriate scale
#' for interpretation or prediction, depending on the model family, link function,
#' and type of prediction. This is an internal helper function for GLMs, survival models,
#' used by mfp2 prediction helpers.
#'
#' @param nfit Numeric vector of linear predictors \eqn{X\beta}.
#' @param family Character string specifying the model family. Supported families include:
#'   - `"gaussian"`: linear regression
#'   - `"binomial"`: binary outcomes or proportions
#'   - `"poisson"`: count data
#'   - `"cox"`: proportional hazards model
#'   - other families fallback to returning `nfit`
#' @param link Character string specifying the link function. Defaults:
#'   - Gaussian: `"identity"`
#'   - Binomial: `"logit"`
#'   - Poisson: `"log"`
#'   Supported links:
#'   - Gaussian: `"identity"`, `"log"`, `"inverse"`
#'   - Binomial: `"logit"`, `"probit"`, `"cloglog"`, `"cauchit"`, `"identity"`, `"log"`
#'   - Poisson: `"log"`, `"identity"`, `"sqrt"`
#' @param type Character string specifying the type of prediction:
#'   - `"response"` or `"risk"`: returns predictions on the response/probability scale
#'   - `NULL` (default): returns the linear predictor itself
#'
#' @return Numeric vector of predictions on the requested scale.
#'
#' @examples
#' \dontrun{
#' # Binomial example
#' lp_bin <- c(-1, 0, 1)
#' transform_linear_predictor(lp_bin, family = "binomial", link = "logit", type = "response")
#'
#' # Poisson example
#' lp_pois <- c(0.5, 1, 1.5)
#' transform_linear_predictor(lp_pois, family = "poisson", link = "log", type = "response")
#'
#' # Gaussian example
#' lp_gauss <- c(1, 2, 3)
#' transform_linear_predictor(lp_gauss, family = "gaussian", link = "log", type = "response")
#' }
#' @keywords internal
#' @noRd
transform_linear_predictor <- function(nfit, family, link = NULL, type = NULL) {

  # Use the canonical link for supported families when the fitted link is absent.
  if (is.null(link)) {
    link <- switch(family,
                   gaussian = "identity",
                   binomial = "logit",
                   poisson  = "log",
                   negbin   = "log",
                   NULL)
  }

  switch(family,

         # Gaussian family with multiple links
         gaussian = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    identity = nfit,
                    log      = exp(nfit),
                    inverse  = 1 / nfit,
                    stop("Unknown Gaussian link"))
           } else nfit
         },

         # Binomial family with multiple links
         binomial = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    logit   = 1 / (1 + exp(-nfit)),
                    probit  = pnorm(nfit),
                    cloglog = 1 - exp(-exp(nfit)),
                    cauchit = pcauchy(nfit),
                    identity = nfit,
                    log     = exp(nfit),
                    stop("Unknown binomial link"))
           } else nfit
         },

         # Poisson and negative-binomial count models use the same links.
         poisson = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    log      = exp(nfit),
                    identity = nfit,
                    sqrt     = nfit^2,
                    stop("Unknown Poisson link"))
           } else nfit
         },

         negbin = {
           if (!is.null(type) && type == "response") {
             switch(link,
                    log      = exp(nfit),
                    identity = nfit,
                    sqrt     = nfit^2,
                    stop("Unknown negative-binomial link"))
           } else nfit
         },

         # Cox proportional hazards: exponentiate for risk
         cox = {
           if (!is.null(type) && type %in% c("response", "risk")) exp(nfit) else nfit
         },

         # Unknown families remain on the linear-predictor scale.
         nfit
  )
}

#' Does a fitted term require positive raw prediction input?
#'
#' Internal helper used by `predict.mfp2()`.
#'
#' The prediction method applies the shift values learned during fitting before
#' constructing transformed variables. This helper determines whether the
#' shifted raw covariate for variable `v` must be strictly positive, based on
#' the final fitted model metadata.
#'
#' For ordinary FP terms, all selected powers in `object$fp_powers[[v]]` apply
#' directly to the raw covariate.
#'
#' For ACD terms, `object$fp_powers[[v]][1]` applies to the raw covariate and
#' `object$fp_powers[[v]][2]` applies to the ACD component A(x). Therefore only
#' the first power is checked against the raw covariate here.
#'
#' Variables fitted with zero-handling are exempt because nonpositive values are
#' represented structurally. Variables with `spike_decision = 3` are also exempt
#' because the continuous component is dropped and only the spike/catzero
#' indicator is retained.
#'
#' @param object A fitted `mfp2` object.
#' @param v Character scalar giving the conceptual term name.
#'
#' @return A single logical value.
#'
#' @keywords internal
#' @noRd
requires_positive_raw_input <- function(object, v) {
  powers <- object$fp_powers[[v]]

  if (is.null(powers) || all(is.na(powers))) {
    return(FALSE)
  }

  if (!is.null(object$zero) && isTRUE(object$zero[[v]])) {
    return(FALSE)
  }

  if (!is.null(object$spike_dec) &&
      v %in% names(object$spike_dec) &&
      !is.na(object$spike_dec[[v]]) &&
      as.integer(object$spike_dec[[v]]) == saz_decision_codes[["binary_only"]]) {
    return(FALSE)
  }

  is_acd <- !is.null(object$acd) &&
    v %in% names(object$acd) &&
    isTRUE(object$acd[[v]])

  raw_powers <- if (is_acd) {
    # For ACD terms:
    #   powers[1] applies to raw x
    #   powers[2] applies to A(x)
    #
    # Only the raw-x power determines whether shifted raw newdata must be
    # strictly positive at this stage.
    powers[1L]
  } else {
    powers
  }

  fp_power_requires_positive_input(raw_powers)
}

#' Resolve Formula Labels for Active Prediction Terms
#'
#' Formula fits retain original formula labels for evaluation and conceptual
#' term names for MFP selection. This helper links the two representations so
#' prediction can evaluate only the terms present in the final model.
#'
#' @keywords internal
#' @noRd
prediction_formula_source_name <- function(term_label) {
  factor_source <- formula_factor_source_name(term_label)
  if (!is.null(factor_source)) {
    return(factor_source)
  }

  formula_fp_source_name(term_label)
}


prediction_formula_term_map <- function(object) {
  map <- object$formula_prediction_term_names
  labels <- attr(object$formula_terms, "term.labels")

  if (!is.null(map)) {
    if (!is.character(map) || is.null(names(map)) || anyDuplicated(names(map))) {
      stop("The fitted formula prediction-term map is malformed.",
           call. = FALSE)
    }
    missing_labels <- setdiff(labels, names(map))
    if (length(missing_labels) > 0L) {
      stop(
        "The fitted formula prediction-term map is missing formula term(s): ",
        paste(missing_labels, collapse = ", "),
        call. = FALSE
      )
    }
    return(map[labels])
  }

  # Backward-compatible fallback for objects fitted before the explicit map was
  # stored. Bare terms and simple factor wrappers can be resolved safely. More
  # complex legacy formula objects retain the previous complete-formula replay.
  map <- stats::setNames(labels, labels)
  lookup <- prediction_term_to_columns(object)
  for (label in labels) {
    source_name <- prediction_formula_source_name(label)
    if (!is.null(source_name) && source_name %in% names(lookup)) {
      map[[label]] <- source_name
    }
  }
  map
}


#' Restrict Stored Formula Terms to Active Conceptual Terms
#'
#' @keywords internal
#' @noRd
prediction_formula_terms <- function(object, terms) {
  terms_object <- object$formula_terms
  term_map <- prediction_formula_term_map(object)
  labels <- attr(terms_object, "term.labels")
  active_labels <- names(term_map)[unname(term_map) %in% terms]
  keep_labels <- active_labels

  # Preserve lower-order terms that use only variables from an active term.
  # model.matrix() uses the term hierarchy when deciding factor contrast coding;
  # retaining these support terms reproduces the fit-time columns for active
  # interactions without reintroducing unrelated predictor dependencies.
  factors <- attr(terms_object, "factors")
  if (!is.null(factors) && length(active_labels) > 0L) {
    for (active_label in active_labels) {
      active_variables <- rownames(factors)[factors[, active_label] != 0]
      support_labels <- colnames(factors)[
        vapply(
          seq_len(ncol(factors)),
          function(index) {
            term_variables <- rownames(factors)[factors[, index] != 0]
            length(term_variables) > 0L &&
              all(term_variables %in% active_variables)
          },
          logical(1L)
        )
      ]
      keep_labels <- union(keep_labels, support_labels)
    }
  }

  drop_indices <- which(!labels %in% keep_labels)
  if (length(drop_indices) > 0L) {
    terms_object <- terms_object[-drop_indices]
  }

  terms_object
}


#' Rebuild formula-interface newdata with model.matrix()
#'
#' Formula-fitted mfp2 objects are fitted on the expanded numeric design matrix
#' returned by model.matrix(). This helper lets users pass ordinary newdata with
#' the original formula variables, including factors, and reconstructs the same
#' expanded columns used during fitting. If all fitted design columns are already
#' present, no formula reconstruction is necessary.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_newdata <- function(object, newdata, terms = NULL) {
  # Matrix-interface fits already receive raw design columns and need no formula
  # evaluation. Formula fits are reconstructed with stored levels and contrasts.
  if (!isTRUE(object$formula_interface)) {
    return(newdata)
  }

  if (is.null(terms)) {
    terms <- get_selected_variable_names(object)
  }
  lookup <- prediction_term_to_columns(object, terms = terms)
  expected <- unlist(lookup, use.names = FALSE)

  # A final intercept-only model has no predictor dependencies. Preserve the
  # requested row count without evaluating any original formula expression.
  if (length(expected) == 0L) {
    return(data.frame(row.names = seq_len(nrow(newdata))))
  }

  # Accept callers that already supplied the exact active model-matrix columns.
  if (all(expected %in% colnames(newdata))) {
    return(newdata[, expected, drop = FALSE])
  }

  if (is.null(object$formula_terms)) {
    stop(
      "! This object was fitted using the formula interface, but formula terms ",
      "needed for prediction were not stored.",
      call. = FALSE
    )
  }

  newdata_df <- as.data.frame(newdata)
  active_formula_terms <- prediction_formula_terms(object, terms)
  active_factors <- attr(active_formula_terms, "factors")
  active_frame_variables <- if (is.null(active_factors)) {
    character(0L)
  } else {
    rownames(active_factors)[rowSums(active_factors != 0) > 0L]
  }
  active_xlevels <- object$formula_xlevels[
    intersect(names(object$formula_xlevels), active_frame_variables)
  ]

  # model.frame() evaluates only active formula expressions and enforces their
  # fitted factor levels; model.matrix() recreates the fitted contrast coding.
  mf <- stats::model.frame(
    active_formula_terms,
    data = newdata_df,
    na.action = stats::na.pass,
    xlev = active_xlevels
  )
  active_contrasts <- object$formula_contrasts[
    intersect(names(object$formula_contrasts), names(mf))
  ]
  mm <- stats::model.matrix(
    active_formula_terms,
    data = mf,
    contrasts.arg = active_contrasts
  )

  keep_cols <- colnames(mm) != "(Intercept)"
  mm <- mm[, keep_cols, drop = FALSE]

  # Restore the exact column names used by the fitted object after formula
  # sanitisation or internal fp()/fp2() renaming.
  column_map <- object$formula_column_map
  if (!is.null(column_map)) {
    mapped <- unname(column_map[colnames(mm)])
    colnames(mm) <- ifelse(is.na(mapped), colnames(mm), mapped)
  }

  missing <- setdiff(expected, colnames(mm))
  if (length(missing) > 0L) {
    stop(
      "! Could not reconstruct required model-matrix column(s) from newdata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  mm[, expected, drop = FALSE]
}

#' Reconstruct formula-level Cox strata from prediction data
#'
#' Reconstructs the original formula-level Cox stratification term from
#' `newdata` for `mfp2` objects fitted through the formula interface.
#'
#' This helper deliberately mimics the high-level `survival::coxph()` strata
#' handling:
#'
#' - for a single strata variable, it returns the evaluated column from the
#'   model frame;
#' - for multiple strata variables, it combines them with
#'   `survival::strata(..., shortlabel = TRUE)`;
#' - it does **not** convert the result to integer codes.
#'
#' Integer conversion belongs only in the low-level Cox fitting path immediately
#' before calling `survival::coxph.fit()`, analogous to the `istrat` object used
#' internally by `survival::coxph()`.
#'
#' @param object A fitted `mfp2` object.
#' @param newdata A user-supplied prediction data frame before formula
#'   reconstruction has been applied.
#'
#' @return Either `NULL`, if no formula-level strata metadata are available, or
#'   a Cox strata object suitable for passing as the `strata` argument to
#'   `prepare_newdata_for_predict()`.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_strata_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_strata_terms)) {
    return(NULL)
  }

  newdata_df <- as.data.frame(newdata)

  mf <- tryCatch(
    stats::model.frame(
      object$formula_strata_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_strata_xlevels
    ),
    error = function(e) {
      stop(
        "! This `mfp2` object was fitted with formula-level Cox strata, ",
        "but the strata term could not be reconstructed from `newdata`.\n",
        "i Include the original strata variable(s) in `newdata`, or supply ",
        "`strata = ...` explicitly to `predict()`.\n",
        "i Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (NROW(mf) != NROW(newdata_df)) {
    stop(
      "! Formula-level Cox strata reconstruction returned the wrong number of rows.",
      call. = FALSE
    )
  }

  if (anyNA(mf)) {
    stop(
      "! Reconstructed Cox strata contain missing values.\n",
      "i Please remove missing strata values from `newdata` before prediction.",
      call. = FALSE
    )
  }

  if (NCOL(mf) == 1L) {
    mf[[1L]]
  } else {
    do.call(
      survival::strata,
      c(as.list(mf), list(shortlabel = TRUE))
    )
  }
}

#' Retrieve and Validate the Cox Offset Origin Stored at Fit Time
#'
#' `survival::predict.coxph()` handles offsets independently of the covariate
#' `reference`. For both training and new-data predictions it subtracts the mean
#' training offset, so a supplied prediction offset contributes as
#' `new_offset - mean(training_offset)`.
#'
#' The final `fit_mfp()` call stores that mean once as `cox_offset_reference`.
#' This helper validates the metadata before delegation. It intentionally does
#' not alter prediction offsets or combine them with covariate centering; the
#' native Cox method remains the single implementation of that arithmetic.
#'
#' A malformed value normally indicates a Cox object created by an older package
#' version or an object whose internals were modified after fitting. Failing here
#' gives a direct remediation message rather than allowing an unexplained shift
#' in the predicted linear predictor.
#'
#' @param object A fitted `mfp2` Cox object containing
#'   `cox_offset_reference`.
#'
#' @return Unnamed finite numeric scalar equal to the mean offset used in the
#'   final Cox fitting data.
#'
#' @keywords internal
#' @noRd
mfp2_cox_offset_reference <- function(object) {
  reference <- object$cox_offset_reference

  if (!is.numeric(reference) || length(reference) != 1L ||
      is.na(reference) || !is.finite(reference)) {
    stop(
      paste0(
        "The fitted Cox model lacks valid offset-reference metadata. ",
        "Refit the mfp2 object with the current package version."
      ),
      call. = FALSE
    )
  }

  unname(reference)
}


#' Rebuild a Formula-Level Offset from Prediction Newdata
#'
#' Formula-interface fits store a response-free terms object containing only the
#' original `offset()` expression and the factor levels required to evaluate it.
#' This helper evaluates that expression in ordinary user-facing `newdata` and
#' returns the raw offset vector expected by the internal fitted Cox or GLM
#' formula.
#'
#' No centering is performed here. For Cox models the raw vector is attached as
#' `offset_`, after which `predict.coxph()` subtracts the mean training offset.
#' Keeping reconstruction separate from centering prevents the offset origin
#' from being applied twice and keeps GLM and Cox delegation consistent with
#' their native prediction methods.
#'
#' The helper returns `NULL` for matrix-interface fits and for formula fits that
#' did not contain an offset. Evaluation errors identify the missing original
#' variable and direct callers to the explicit `newoffset` alternative.
#'
#' @param object Fitted `mfp2` object.
#' @param newdata User-facing prediction data before formula reconstruction.
#'
#' @return `NULL` when no formula offset is stored, otherwise a finite numeric
#'   vector with one value per prediction row.
#'
#' @keywords internal
#' @noRd
reconstruct_formula_offset_newdata <- function(object, newdata) {
  if (!isTRUE(object$formula_interface) ||
      is.null(object$formula_offset_terms)) {
    return(NULL)
  }

  newdata_df <- as.data.frame(newdata)

  # Validate the raw numeric variables used by the stored offset expression
  # before model.frame() evaluates transformations such as log(). This gives
  # missing and non-finite inputs a deterministic prediction error and prevents
  # transformation warnings (for example, "NaNs produced") from escaping.
  offset_variables_call <- attr(object$formula_offset_terms, "variables")
  offset_source_variables <- if (is.null(offset_variables_call)) {
    character()
  } else {
    unique(all.vars(offset_variables_call))
  }
  offset_source_columns <- intersect(offset_source_variables, names(newdata_df))
  if (length(offset_source_columns) > 0L) {
    mfp2_validate_prediction_newdata(
      newdata_df[, offset_source_columns, drop = FALSE]
    )
  }

  mf <- tryCatch(
    stats::model.frame(
      object$formula_offset_terms,
      data = newdata_df,
      na.action = stats::na.pass,
      xlev = object$formula_offset_xlevels
    ),
    error = function(e) {
      stop(
        "! This `mfp2` object was fitted with a formula-level offset, ",
        "but the offset could not be reconstructed from `newdata`.\n",
        "i Include the original offset variable(s) in `newdata`, or supply ",
        "`newoffset = ...` explicitly to `predict()`.\n",
        "i Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  out <- stats::model.offset(mf)

  if (is.null(out)) {
    stop(
      "! Formula-level offset could not be evaluated from `newdata`.",
      call. = FALSE
    )
  }

  out <- as.vector(out)

  if (length(out) != NROW(newdata_df)) {
    stop(
      "! Formula-level offset reconstruction returned the wrong number of rows.",
      call. = FALSE
    )
  }

  if (anyNA(out) || any(!is.finite(out))) {
    stop(
      "! Reconstructed offset contains missing or non-finite values.",
      call. = FALSE
    )
  }

  out
}

#' Helper function to prepare newdata for predict function
#'
#' To be used in \code{predict.mfp2()}.
#'
#' @param object fitted `mfp2` model object.
#' @param newdata Numeric design data on the raw model-matrix scale. Full-model
#'   prediction requires every fitted term. Internal term prediction may supply
#'   a subset consisting of complete conceptual terms.
#' @param strata,offset passed from \code{predict.mfp2()}. For Cox
#'   prediction, \code{strata} must be kept as a high-level vector/factor or
#'   combined multi-column strata object; integer conversion is not performed
#'   here because the stored Cox formula evaluates \code{strata(strata_)}.
#' @param terms Character vector of conceptual terms required for this prediction.
#'   Defaults to the terms selected in the final model.
#' @param apply_pre logical indicating whether the fitted pre-transformation
#' is applied or not.
#' @param apply_center logical indicating whether the fitted centers are applied
#' after transformation or not.
#' @param check_binary passed to \code{transform_vector_fp()}.
#' @param reset_zero Logical. If `TRUE`, variables marked as `zero = TRUE`
#' but containing only positive values are reset to `FALSE` before transformation.
#' The prediction helper defaults to `FALSE` because prediction must preserve
#' the zero/catzero/spike structure learned during model fitting. Adaptive
#' resetting based on the contents of `newdata` can create a design matrix that
#' is inconsistent with the fitted coefficients. Passed to
#' \code{transform_matrix()}.
#' @param allow_missing_predictors Logical. If `FALSE`, the default, `newdata`
#'   must contain every raw input column required by the fitted term lookup.
#'   For a grouped term, all member columns are required together. If `TRUE`,
#'   entire terms may be omitted and only complete terms present in `newdata`
#'   are transformed. This is intended for internal term-specific prediction,
#'   not for ordinary full-model prediction.
#' @return A data frame containing the transformed prediction design and any
#'   required Cox strata or offset columns.
#' @keywords internal
#' @noRd
prepare_newdata_for_predict <- function(object,
                                        newdata,
                                        terms = NULL,
                                        strata = NULL,
                                        offset = NULL,
                                        apply_pre = TRUE,
                                        apply_center = TRUE,
                                        check_binary = TRUE,
                                        reset_zero = FALSE,
                                        allow_missing_predictors = FALSE) {
  # Keep formula-level strata and offset variables out of the numeric predictor
  # matrix. A data frame containing a factor stratum would otherwise be coerced
  # wholesale to a character matrix before the fitted predictor columns are
  # selected. Predictor isolation therefore precedes matrix conversion.
  newdata <- as.data.frame(newdata, check.names = FALSE)
  n_newdata <- nrow(newdata)

  if (is.null(colnames(newdata)) || any(colnames(newdata) == "")) {
    stop("! newdata must have non-empty column names.", call. = FALSE)
  }

  # Resolve the conceptual-term lookup once. Full prediction requires every raw
  # member column; term-specific prediction may include complete term subsets.
  if (is.null(terms)) {
    terms <- get_selected_variable_names(object)
  }
  lookup <- prediction_term_to_columns(object, terms = terms)
  expected <- unlist(lookup, use.names = FALSE)
  missing <- setdiff(expected, colnames(newdata))

  if (length(missing) > 0L && !isTRUE(allow_missing_predictors)) {
    stop(
      "! Missing required predictor column(s) in newdata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Count supplied raw columns per conceptual term. Partial categorical blocks
  # are invalid because they would change the fitted parameterisation.
  present_by_term <- vapply(
    lookup,
    function(columns) sum(columns %in% colnames(newdata)),
    integer(1L)
  )
  partial_terms <- names(lookup)[
    present_by_term > 0L & present_by_term < lengths(lookup)
  ]

  if (length(partial_terms) > 0L) {
    details <- vapply(
      partial_terms,
      function(term) {
        absent <- setdiff(lookup[[term]], colnames(newdata))
        paste0(term, " (missing: ", paste(absent, collapse = ", "), ")")
      },
      character(1L)
    )
    stop(
      "! Grouped terms must be supplied with all member columns. Incomplete term(s): ",
      paste(details, collapse = "; "),
      call. = FALSE
    )
  }

  active_terms <- names(lookup)[present_by_term == lengths(lookup)]

  # Expand complete conceptual terms in lookup order. This keeps every grouped
  # block contiguous and gives transform_matrix() raw-column metadata in the
  # same order as the supplied design block.
  raw_columns <- unlist(lookup[active_terms], use.names = FALSE)

  if (length(raw_columns) == 0L) {
    newdata <- data.frame(row.names = seq_len(n_newdata))
  } else {
    newdata <- newdata[, raw_columns, drop = FALSE]

    non_numeric <- raw_columns[
      !vapply(newdata, is.numeric, logical(1L))
    ]
    if (length(non_numeric) > 0L) {
      stop(
        "! Reconstructed predictor column(s) must be numeric: ",
        paste(non_numeric, collapse = ", "),
        call. = FALSE
      )
    }

    # Keep this guard immediately before shifting and transformation. It also
    # protects internal callers that bypass the public method's earlier check.
    mfp2_validate_prediction_newdata(newdata)

    newdata <- as.matrix(newdata)
    storage.mode(newdata) <- "double"
    metadata <- expand_prediction_metadata(object, raw_columns)

    # Apply the fitted shift only to raw user data. Stored x_original data and
    # internally generated term grids are already on the shifted prediction scale.
    if (apply_pre) {
      newdata <- sweep(newdata, 2L, metadata$shift[raw_columns], "+")

      positive_terms <- active_terms[
        vapply(
          active_terms,
          function(term) {
            length(lookup[[term]]) == 1L && requires_positive_raw_input(object, term)
          },
          logical(1L)
        )
      ]

      if (length(positive_terms) > 0L) {
        bad_terms <- positive_terms[
          vapply(
            positive_terms,
            function(term) {
              column <- lookup[[term]][[1L]]
              any(!is.na(newdata[, column]) & newdata[, column] <= 0)
            },
            logical(1L)
          )
        ]

        if (length(bad_terms) > 0L) {
          bad_summary <- vapply(
            bad_terms,
            function(term) {
              column <- lookup[[term]][[1L]]
              count <- sum(!is.na(newdata[, column]) & newdata[, column] <= 0)
              paste0(term, " (", count, if (count == 1L) " row)" else " rows)")
            },
            character(1L)
          )

          stop(
            "After applying the shift values learned during fitting, some values in ",
            "`newdata` remain non-positive for terms whose fitted transformation ",
            "requires strictly positive input.\n",
            "i Problematic term(s): ", paste(bad_summary, collapse = ", "), ".",
            call. = FALSE
          )
        }
      }
    }

    # ACD prediction needs the parameters estimated at fit time for each active
    # singleton ACD term.
    active_acd <- raw_columns[
      metadata$acdx &
        !vapply(metadata$powers, function(p) is.null(p) || all(is.na(p)), logical(1L))
    ]
    missing_acd <- active_acd[
      vapply(
        active_acd,
        function(column) is.null(metadata$acd_parameter[[column]]),
        logical(1L)
      )
    ]

    if (length(missing_acd) > 0L) {
      stop(
        "Missing stored ACD parameters for prediction column(s): ",
        paste(missing_acd, collapse = ", "),
        call. = FALSE
      )
    }

    # Recreate the exact uncentred transformed columns used by the final model.
    x_trans <- transform_matrix(
      newdata,
      power_list = metadata$powers,
      center = stats::setNames(rep(FALSE, length(raw_columns)), raw_columns),
      keep_x_order = TRUE,
      acdx = metadata$acdx,
      acd_parameter_list = metadata$acd_parameter,
      check_binary = check_binary,
      zero = metadata$zero,
      catzero = metadata$catzero,
      spike = metadata$spike,
      spike_decision = metadata$spike_decision,
      reset_zero = reset_zero
    )

    if (is.null(x_trans)) {
      newdata <- matrix(numeric(0L), nrow = n_newdata, ncol = 0L)
    } else {
      newdata <- x_trans$x_transformed

      # Apply the stored final-model centres after transformation.
      if (apply_center && !is.null(object$centers)) {
        newdata <- center_matrix(
          newdata,
          centers = object$centers[colnames(newdata)],
          zero = x_trans$zero_expanded
        )
      }
    }

    if (NCOL(newdata) == 0L) {
      newdata <- data.frame(row.names = seq_len(n_newdata))
    } else {
      newdata <- data.frame(newdata, check.names = FALSE)
    }
  }

  # Attach high-level strata and offset variables expected by the stored model
  # formula after the predictor design has been reconstructed.
  if (object$family_string == "cox" && !is.null(strata)) {
    strata_n <- if (is.vector(strata) || is.factor(strata)) length(strata) else NROW(strata)
    if (strata_n != nrow(newdata)) {
      stop("! `strata` must have one value or row per prediction row.", call. = FALSE)
    }
    if (anyNA(strata)) {
      stop("! `strata` must not contain missing values.", call. = FALSE)
    }
    strata_frame <- as.data.frame(strata, check.names = FALSE)
    bad_numeric_strata <- names(strata_frame)[vapply(
      strata_frame,
      function(column) is.numeric(column) && any(!is.finite(column)),
      logical(1L)
    )]
    if (length(bad_numeric_strata) > 0L) {
      stop("! Numeric `strata` values must be finite.", call. = FALSE)
    }

    newdata$strata_ <- if (is.matrix(strata) || is.data.frame(strata)) {
      do.call(
        survival::strata,
        c(as.list(as.data.frame(strata)), list(shortlabel = TRUE))
      )
    } else {
      strata
    }
  }

  if (!is.null(offset)) {
    if (!is.numeric(offset) || anyNA(offset) || any(!is.finite(offset))) {
      stop("! offset must be a finite numeric vector.", call. = FALSE)
    }
    if (length(offset) != nrow(newdata)) {
      stop("! offset must have one value per observation.", call. = FALSE)
    }
    newdata$offset_ <- offset
  }

  newdata
}
#' Helper function to compute standard error of a partial predictor
#'
#' To be used in \code{predict.mfp2()}.
#'
#' @param model fitted `mfp2` object.
#' @param X transformed input matrix with variables of interest for partial predictor.
#' @param xref transformed reference value for variable of interest. Default is
#'  `NULL`, in which case this function computes standard errors without reference
#' values.
#' @param include_intercept logical indicating whether intercept variance and
#' covariance should be included when no reference value is supplied. This must
#' match the `add_intercept` choice used for term predictions.
#'
#' @details
#' See pages 91-92 and following in the book by Royston and Sauerbrei 2008
#' for the formulas and mathematical details.
#'
#' @return
#' Numeric vector of conditional standard errors.
#'
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building:
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' @keywords internal
#' @noRd
calculate_standard_error <- function(model,
                                     X,
                                     xref = NULL,
                                     include_intercept = TRUE) {

  vcovx <- vcov(object = model)

  # Small or rank-deficient fitted samples can yield undefined covariance entries.
  if (any(is.nan(vcovx)))
    warning("i NaN detected in the covariance matrix of the model.",
            "i Standard errors for calculation of confidence intervals may not exist")

  # Use the exact fit-time mapping between transformed matrix columns and
  # fitted coefficient names. This handles non-syntactic formula terms without
  # modifying either the transformed matrix names or the model covariance names.
  xnames <- colnames(X)
  model_columns <- prediction_model_column_names(model, xnames)
  ind <- match(model_columns, colnames(vcovx))

  if (anyNA(ind)) {
    stop(
      "Cannot compute SE; covariance matrix lacks mapped columns: ",
      paste(model_columns[is.na(ind)], collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(xref)) {
    # Subtract the reference value: f(x)-f(xref)
    X <- sweep(X, 2, xref, "-")
  }

  # Augment X by the intercept only when the reported term value includes
  # the intercept. Contrasts and Cox partial predictors do not include intercept
  # variance.
  if (isTRUE(include_intercept) && model$family_string != "cox" && is.null(xref)) {
    intercept_ind <- match("(Intercept)", colnames(vcovx))
    if (is.na(intercept_ind)) {
      stop(
        "Cannot compute SE with intercept; covariance matrix lacks `(Intercept)`.",
        call. = FALSE
      )
    }
    X <- cbind("(Intercept)" = 1, X)
    ind <- c(intercept_ind, ind)
  }

  # Compute diag(X V X') without constructing the full n-by-n product.
  vcovx <- vcovx[ind, ind, drop = FALSE]
  v <- rowSums((X %*% vcovx) * X)

  sqrt(v)
}
