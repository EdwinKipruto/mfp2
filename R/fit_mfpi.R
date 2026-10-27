# Internal implementation functions for mfpi()
#
# None of these functions are exported. They are called exclusively by
# mfpi.default() after all argument checking and pre-processing has been
# completed. Parameters therefore arrive already validated, expanded to full
# length, shifted, and scaled.

# -----------------------------------------------------------------------------
# Adjustment-power extraction -------------------------------------------------
# -----------------------------------------------------------------------------

#' Extract Canonical Powers for Selected MFPI Adjustment Terms
#'
#' Uses the fitted adjustment model's `fp_powers` list as the primary source of
#' selected powers. Unlike the display-oriented `fp_terms` table, this list
#' preserves structural `NA` positions required by ACD terms. For example, an
#' ACD-only linear term is represented as `c(NA, 1)`, where the first slot is
#' the inactive transformation of x and the second slot is the active
#' transformation of A(x).
#'
#' @param adjustment_model Fitted `mfp2` adjustment-model object.
#' @param selected_vars Character vector of selected conceptual adjustment
#'   terms.
#'
#' @return A named list of numeric power vectors aligned to `selected_vars`.
#'
#' @keywords internal
#' @noRd
mfpi_extract_adjustment_powers <- function(adjustment_model, selected_vars) {
  if (length(selected_vars) == 0L) {
    return(stats::setNames(vector("list", 0L), character(0L)))
  }

  canonical <- adjustment_model$fp_powers

  if (!is.null(canonical)) {
    if (is.null(names(canonical))) {
      stop(
        "Internal error: adjustment-model `fp_powers` must be named.",
        call. = FALSE
      )
    }

    missing_vars <- setdiff(selected_vars, names(canonical))
    if (length(missing_vars) > 0L) {
      stop(
        paste0(
          "Internal error: adjustment-model powers are missing selected term(s): ",
          paste(missing_vars, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    powers <- canonical[selected_vars]
  } else if (!is.null(adjustment_model$fp_terms)) {
    # Compatibility fallback for older fitted objects. This is adequate for
    # ordinary FP terms, but display tables may not preserve inactive ACD power
    # slots; the ACD validation below therefore remains mandatory.
    powers <- get_fp_powers(selected_vars, adjustment_model$fp_terms)
  } else {
    stop(
      paste0(
        "Internal error: the adjustment model has neither `fp_powers` nor ",
        "`fp_terms` metadata."
      ),
      call. = FALSE
    )
  }

  powers <- lapply(powers, function(p) {
    if (is.null(p)) numeric(0L) else unname(as.numeric(p))
  })
  names(powers) <- selected_vars

  missing_power_values <- selected_vars[
    vapply(powers, length, integer(1L)) == 0L
  ]
  if (length(missing_power_values) > 0L) {
    stop(
      paste0(
        "Internal error: selected adjustment term(s) have no stored powers: ",
        paste(missing_power_values, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  acd_flags <- adjustment_model$acd
  acd_selected <- selected_vars[
    vapply(selected_vars, function(v) {
      !is.null(acd_flags) &&
        !is.null(names(acd_flags)) &&
        v %in% names(acd_flags) &&
        isTRUE(acd_flags[[v]])
    }, logical(1L))
  ]

  invalid_acd <- acd_selected[
    vapply(powers[acd_selected], length, integer(1L)) != 2L
  ]
  if (length(invalid_acd) > 0L) {
    stop(
      paste0(
        "Internal error: selected ACD adjustment term(s) must retain two ",
        "power positions (for x and A(x)): ",
        paste(invalid_acd, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  powers
}

# -----------------------------------------------------------------------------
# fit_mfpi() ------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit MFPI Adjustment and Interaction Models
#'
#' Selects an adjustment model and evaluates each specified interaction with
#' the grouping variable separately. Uses inputs prepared by \code{mfpi()}.
#'
#' The selected adjustment terms and transformations are reused in each
#' interaction analysis, excluding the variable being tested. Each interaction
#' uses its resolved functional form; the function does not select between
#' linear, FP1, and FP2 forms.
#'
#' @param x Numeric predictor matrix, including \code{group_var}, with shifts
#'   and scales already applied and no intercept column.
#' @param y Response vector, response matrix, or \code{survival::Surv()} object.
#' @param family Resolved model-family specification.
#' @param family_string Character string giving the canonical family name.
#' @param weights Numeric vector of observation weights.
#' @param offset Linear-predictor offsets: a numeric vector or, for multinomial
#'   models, a matrix with one column per non-reference outcome category.
#' @param has_offset Logical indicating whether an offset was supplied.
#' @param cycles Maximum number of MFP backfitting cycles.
#' @param center Named logical vector controlling predictor centring.
#' @param criterion Selection criterion: \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}. Used for adjustment selection and interaction comparisons.
#' @param select Named numeric vector of significance levels for adjustment
#'   variable selection under \code{criterion = "pvalue"}.
#' @param alpha Named numeric vector of significance levels for adjustment
#'   FP degree selection under \code{criterion = "pvalue"}.
#' @param keep Character vector naming adjustment terms or columns to retain.
#'   Naming a column of a categorical term retains the complete term.
#' @param force_max_fp Named logical vector identifying adjustment variables
#'   whose maximum permitted FP degree must be retained. Powers are still
#'   selected within that degree.
#' @param df Named integer vector of adjustment degrees of freedom, after
#'   restrictions based on the number of distinct values.
#' @param xorder Order of adjustment terms during MFP backfitting:
#'   \code{"ascending"}, \code{"descending"}, or \code{"original"}.
#' @param fp_powers Named list of candidate FP powers for each predictor.
#' @param ties Method for handling tied event times in Cox models.
#' @param strata Optional prepared survival-model stratification information.
#' @param nocenter Numeric values excluded from automatic centring by
#'   \code{survival::coxph()}.
#' @param acd_vars Named logical vector identifying predictors to receive an
#'   approximate cumulative distribution transformation.
#' @param zero_vars Named logical vector identifying predictors whose exact
#'   zeros are handled separately during transformation.
#' @param catzero_vars Named logical vector identifying predictors requiring
#'   a zero indicator alongside their transformed positive values.
#' @param spike_vars Named logical vector identifying predictors eligible for
#'   spike-at-zero modelling.
#' @param min_saz_prop Minimum proportion required in both the zero and
#'   positive components of a spike-at-zero adjustment variable.
#' @param use_ftest Whether to use F-tests instead of likelihood-ratio
#'   chi-square tests for Gaussian p-value comparisons in both fitting stages.
#'   Ignored for other families.
#' @param control List of model-fitting control parameters.
#' @param group_var Character string naming the grouping variable.
#' @param include_group_var Whether to include group main effects during
#'   adjustment-model selection.
#' @param flex Interaction procedure: \code{"flex1"}, \code{"flex2"},
#'   \code{"flex3"}, or \code{"flex4"}.
#' @param interaction_vars Character vector naming terms to test for
#'   interaction with \code{group_var}.
#' @param interaction_forms Named character vector of resolved interaction
#'   forms: \code{"linear"}, \code{"fp1"}, or \code{"fp2"}.
#' @param interaction_specs Named list describing each interaction term's
#'   type and predictor columns.
#' @param p_interact Significance threshold for retaining interactions under
#'   \code{criterion = "pvalue"}.
#' @param min_improvement Minimum AIC or BIC improvement required to retain
#'   an interaction.
#' @param show_models Whether to print coefficient tables for selected
#'   interaction models when \code{verbose = TRUE}.
#' @param verbose Whether to print progress messages.
#' @param digits Number of decimal places for printed numerical results.
#' @param scale Named numeric vector of predictor scale factors.
#' @param shift Named numeric vector of predictor shifts.
#' @param quiet Internal verbosity flag, reset according to \code{verbose}.
#' @param center_type Centring method for interaction transformations:
#'   \code{"grand"} or \code{"group"}.
#' @param p_adjust_method Method passed to \code{stats::p.adjust()} for
#'   interaction p-values under \code{criterion = "pvalue"}.
#' @param group_input_levels Numeric group values in the input matrix.
#' @param group_levels_original Original labels corresponding to
#'   \code{group_input_levels}.
#' @param term_to_columns Named list mapping adjustment terms to predictor
#'   columns. Categorical terms are selected jointly.
#' @param fitter Model-fitting backend.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{best_model_metrics}, \code{all_model_metrics}}{Data frames
#'     of comparison metrics for selected and evaluated interactions,
#'     respectively.}
#'   \item{\code{best_interaction_model},
#'     \code{all_interaction_models}}{Named lists of fitted models for selected
#'     and evaluated interactions, respectively.}
#'   \item{\code{center_vals_list}}{Centring constants for selected interactions.}
#'   \item{\code{adjustment_model}, \code{adjust_terms}}{The adjustment-model
#'     results and its selected-term table.}
#'   \item{\code{adjustment_term_to_columns}}{Mapping of adjustment terms to
#'     predictor columns.}
#'   \item{\code{var_winners}, \code{univariable_interactions}}{Per-variable
#'     results and the complete interaction-analysis results.}
#'   \item{\code{x_train_internal}}{Prepared training predictors, including
#'     shifted and scaled values and recoded groups.}
#' }
#' The list also contains the family, fitting settings, group labels and
#' mappings, interaction specifications, and transformation information.
#'
#' @seealso \code{mfpi()}
#'
#' @keywords internal
#' @noRd
fit_mfpi <- function(x, y, family, family_string, weights, offset, cycles,
                     center, criterion, select, alpha, keep,
                     force_max_fp,
                     df, xorder,
                     fp_powers, ties, strata, nocenter, acd_vars, zero_vars,
                     catzero_vars, spike_vars,  min_saz_prop, use_ftest,
                     control, group_var, include_group_var, flex, interaction_vars,
                     interaction_forms, interaction_specs,
                     p_interact, min_improvement, show_models, verbose,
                     digits, scale, shift = NULL, quiet = FALSE, has_offset,
                     shift_automatic = NULL,
                     center_type = c("grand", "group"),
                     p_adjust_method = "none",
                     group_input_levels = NULL,
                     group_levels_original = NULL,
                     term_to_columns,
                     fitter = "base") {

  center_type <- match.arg(center_type)
  quiet  <- !isTRUE(verbose)
  # Resolve the likelihood-GLM or survival family once for all repeated
  # internal model fits. Prepared survival specifications retain their cached
  # response representation throughout both MFPI stages.
  family_fit <- resolve_fit_model_family(family)

  # ---------------------------------------------------------------------------
  # Pre-processing: synchronise parameter vectors and remap group levels
  # ---------------------------------------------------------------------------
  processed_data <- preprocess_data(
    x                 = x,
    group_var         = group_var,
    include_group_var = include_group_var,
    select            = select,
    alpha             = alpha,
    df                = df,
    center            = center,
    acd_vars          = acd_vars,
    fp_powers         = fp_powers,
    zero_vars         = zero_vars,
    catzero_vars      = catzero_vars,
    spike_vars        = spike_vars,
    keep        = keep,
    force_max_fp      = force_max_fp,
    group_input_levels    = group_input_levels,
    group_levels_original = group_levels_original,
    term_to_columns       = term_to_columns
  )

  # ---------------------------------------------------------------------------
  # Step 1: Fit adjustment model via MFP
  # ---------------------------------------------------------------------------
  if (verbose) {
    rule_thick <- strrep("=", 70)
    rule_thin  <- strrep("-", 70)
    mfp2_message("")
    mfp2_message(rule_thick)
    mfp2_message("  STEP 1: Selected Adjustment Model (MFP)")
    mfp2_message(rule_thick)
  }

  adjustment_model <- fit_adjustment_model(
    x               = processed_data$x,
    y               = y,
    weights         = weights,
    offset          = offset,
    cycles          = cycles,
    family          = family_fit,
    family_string   = family_string,
    fitter          = fitter,
    criterion       = criterion,
    updated_params  = processed_data$updated_params,
    xorder          = xorder,
    ties            = ties,
    strata          = strata,
    nocenter        = nocenter,
    min_saz_prop    = min_saz_prop,
    use_ftest       = use_ftest,
    control         = control,
    force_max_fp    = processed_data$updated_params$force_max_fp,
    scale           = scale,
    has_offset      = has_offset,
    term_to_columns = processed_data$term_to_columns,
    verbose         = FALSE
  )

  # `mfpi.default()` applies the covariate-wide shifts before `fit_mfpi()` is
  # called. Consequently, `fit_adjustment_model()` must pass shift = 0 to
  # `fit_mfp()` so that those shifts are not applied a second time during Stage 1.
  # The zero values stored by `fit_mfp()`, however, describe only that internal
  # fitting step. The returned adjustment model is an ordinary `mfp2` object and
  # must describe preprocessing relative to raw covariate values, because
  # `print.mfp2()` reports `transformations` and `predict.mfp2(newdata = ...)`
  # uses the stored shifts to transform raw new data. Restore the actual MFPI
  # shifts now that fitting is complete. This does not alter fitted coefficients
  # or the training design; it only makes the returned model metadata consistent
  # with the shifted, unscaled scale on which its final FP terms were fitted.
  adjustment_model <- restore_mfpi_adjustment_shifts(
    adjustment_model = adjustment_model,
    shift = shift,
    term_to_columns = processed_data$term_to_columns
  )
  if (!identical(family_string, "negbin")) {
    # fit_mfpi() calls fit_mfp() directly rather than through mfp2.default(),
    # so attach the same public, unprepared family metadata that a standalone
    # mfp2 fit receives. This keeps summary/prediction on the nested adjustment
    # model complete for survreg families.
    adjustment_model$family <- mfp2_strip_prepared_family(family_fit)
  }
  # fit_mfpi() calls fit_mfp() below the public mfp2.default() boundary, so
  # mirror the observation-level offset metadata that mfp2.default() normally
  # adds.
  adjustment_model$offset <- offset
  adjustment_model$has_offset <- has_offset

  # Identify variables retained by the adjustment model
  selected_vars <- get_selected_variables(adjustment_model)
  if (include_group_var) {
    # Group dummies form one conceptual term named by group_var. Exclude that
    # term from the ordinary adjustment set used in interaction models.
    selected_vars <- setdiff(
      selected_vars,
      c(group_var, processed_data$dummy_names)
    )
  }

  if (verbose) {
    adjustment_info <- mfpi_prepare_adjustment_display(
      list(
        criterion = criterion,
        p_interact = p_interact,
        min_improvement = min_improvement,
        p_adjust_method = p_adjust_method,
        adjust_terms = adjustment_model$fp_terms
      ),
      digits = digits
    )

    mfp2_message("")
    mfp2_message(
      "  Model: ",
      mfpi_print_model_label(list(
        family_string = family_string,
        family = mfp2_strip_prepared_family(family_fit)
      ))
    )
    mfp2_message("  Criterion: ", adjustment_info$settings$criterion_label)
    mfp2_message("")

    if (!is.null(adjustment_info$display) &&
        nrow(adjustment_info$display) > 0L) {
      mfp2_message_capture(print(adjustment_info$display))
      mfp2_message("")
      mfp2_message(
        "  Selected adjustment variables (",
        length(adjustment_info$selected_names),
        "): ",
        paste(adjustment_info$selected_names, collapse = ", ")
      )
    } else {
      mfp2_message("  No adjustment variables selected.")
    }
  }

  # Use the adjustment model's canonical power metadata. In particular,
  # ACD terms must retain both structural power positions (x and A(x)); the
  # display-oriented fp_terms table can collapse an inactive NA slot and is
  # therefore not a safe source for rebuilding the adjustment design matrix.
  adj_fp_powers <- mfpi_extract_adjustment_powers(
    adjustment_model = adjustment_model,
    selected_vars = selected_vars
  )

  # Extract spike decisions from the adjustment model for use when
  # re-transforming adjustment variables inside evaluate_interactions().
  # spike_decision is a named numeric vector (1 = FP+binary, 2 = FP only,
  # 3 = binary only) stored in the mfp2 object by fit_mfp().
  adj_spike_decision <- if (!is.null(adjustment_model$spike_dec)) {
    adjustment_model$spike_dec
  } else {
    setNames(rep(2L, length(selected_vars)), selected_vars)
  }

  adj_spike_decision <- adj_spike_decision[selected_vars]

  # Developer note:
  # fit_mfp() may reset zero/catzero/spike flags during adjustment-model fitting,
  # especially through the spike-at-zero eligibility checks. When the selected
  # adjustment variables are re-transformed inside evaluate_interactions(), we
  # must use these final post-fit flags, not the original pre-fit flags from
  # processed_data$updated_params.
  adj_zero <- mfpi_extract_adjustment_flag(
    adjustment_model,
    "zero",
    selected_vars,
    default = FALSE
  )

  adj_catzero <- mfpi_extract_adjustment_flag(
    adjustment_model,
    "catzero",
    selected_vars,
    default = FALSE
  )

  adj_spike <- mfpi_extract_adjustment_flag(
    adjustment_model,
    "spike",
    selected_vars,
    default = FALSE
  )

  # ---------------------------------------------------------------------------
  # Step 2: Evaluate univariable interactions
  # ---------------------------------------------------------------------------
  if (verbose) {
    rule_thick <- strrep("=", 70)
    interaction_settings <- mfpi_interaction_settings(
      criterion = criterion,
      p_interact = p_interact,
      min_improvement = min_improvement,
      p_adjust_method = p_adjust_method,
      digits = digits
    )

    mfp2_message("")
    mfp2_message(rule_thick)
    mfp2_message(
      sprintf(
        "  STEP 2: Evaluating Interactions (flex = %s)",
        flex
      )
    )
    mfp2_message(rule_thick)
    mfp2_message("")
    mfp2_message("  Group variable: ", group_var)
    mfp2_message("  Selection rule: ", interaction_settings$selection)
    if (criterion == "pvalue") {
      mfp2_message(
        "  P-value adjustment: ",
        interaction_settings$p_adjust_method
      )
    }
  }

  univ_results <- evaluate_interactions(
    y                  = y,
    processed_data     = processed_data,
    selected_vars      = selected_vars,
    adj_fp_powers      = adj_fp_powers,
    adj_zero           = adj_zero,
    adj_catzero        = adj_catzero,
    adj_spike          = adj_spike,
    adj_spike_decision = adj_spike_decision,
    interaction_vars   = interaction_vars,
    interaction_forms  = interaction_forms,
    interaction_specs  = interaction_specs,
    flex               = flex,
    weights            = weights,
    offset             = offset,
    has_offset         = has_offset,
    xorder             = xorder,
    ties               = ties,
    strata             = strata,
    use_ftest          = use_ftest,
    control            = control,
    nocenter           = nocenter,
    family             = family_fit,
    family_string      = family_string,
    fitter             = fitter,
    fp_powers          = fp_powers,
    cycles             = cycles,
    criterion          = criterion,
    group_var          = group_var,
    show_models        = show_models,
    p_interact         = p_interact,
    min_improvement    = min_improvement,
    min_saz_prop       = min_saz_prop,
    center_type        = center_type,
    scale              = scale,
    shift              = shift,
    shift_automatic    = shift_automatic,
    adj_acd_parameter  = adjustment_model$acd_parameter,
    quiet              = quiet,
    p_adjust_method    = p_adjust_method,
    digits             = digits,
    verbose            = verbose
  )

  # Prepared family data are a fitting cache, not public model metadata. Strip
  # that cache before storing the full Stage-2 result.
  univ_results$family <- mfp2_strip_prepared_family(univ_results$family)

  n_selected <- nrow(univ_results$best_model_metrics)

  if (verbose) {
    mfp2_message(
      "Done: ",
      n_selected,
      " of ",
      length(interaction_vars),
      " interactions selected."
    )
    mfp2_message("")
  }

  # Every proportional-hazards interaction candidate is fitted on the same
  # working rows with the same training offset as the final adjustment model.
  # fit_mfp() has already calculated predict.coxph()'s offset origin once, so
  # expose that scalar at the top level rather than recomputing it per term.
  cox_offset_reference <- if (identical(family_string, "cox")) {
    adjustment_model$cox_offset_reference
  } else {
    NULL
  }

  list(
    best_model_metrics       = univ_results$best_model_metrics,
    all_model_metrics        = univ_results$all_model_metrics,
    best_interaction_model   = univ_results$best_interaction_model,
    all_interaction_models   = univ_results$all_interaction_models,
    center_vals_list         = univ_results$center_vals_list,
    adjust_terms             = adjustment_model$fp_terms,
    group_var                = univ_results$group_var,
    include_group_var        = include_group_var,
    show_models              = univ_results$show_models,
    flex                     = univ_results$flex,
    criterion                = univ_results$criterion,
    group_levels_new         = univ_results$group_levels_new,
    group_levels_original    = univ_results$group_levels_original,
    group_level_map          = univ_results$group_level_map,
    family                   = univ_results$family,
    fitter                   = fitter,
    nobs                     = univ_results$nobs,
    cox_offset_reference     = cox_offset_reference,
    p_adjust_method          = univ_results$p_adjust_method,
    scale                    = scale,
    shift                    = shift,
    # Training matrix after MFPI preprocessing.
    #
    # This matrix contains shifted/scaled predictors and the internally remapped
    # group_var codes produced by preprocess_data(). It is the correct matrix for
    # manual training-data prediction reconstruction, especially when
    # predict.mfpi(newdata = NULL, offset = ...) cannot delegate to the stored fit.
    x_train_internal           = processed_data$original_x,

    var_winners                = univ_results$var_winners,
    adjustment_model           = adjustment_model,
    adjustment_term_to_columns = processed_data$term_to_columns,
    interaction_specs          = interaction_specs,
    univariable_interactions   = univ_results,
    digits                     = digits
  )
}


# -----------------------------------------------------------------------------
# preprocess_data() -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Prepare Predictors for the MFPI Adjustment Model
#'
#' Removes the grouping variable from the adjustment predictors and optionally
#' replaces it with dummy variables for its non-reference levels. Aligns the
#' modelling settings and term-to-column mapping with the resulting matrix.
#'
#' Group values are recoded as consecutive integers starting at zero.
#' Original labels and their mapping to these codes are retained.
#'
#' @param x Numeric predictor matrix including \code{group_var}.
#' @param group_var Character string naming the grouping variable.
#' @param include_group_var Whether to include group dummy variables in the
#'   adjustment model.
#' @param select Named numeric vector of variable-selection significance levels.
#' @param alpha Named numeric vector of FP-degree significance levels.
#' @param df Named numeric vector of adjustment degrees of freedom.
#' @param center Named logical vector controlling predictor centring.
#' @param acd_vars Named logical vector identifying predictors to receive an
#'   approximate cumulative distribution transformation.
#' @param fp_powers Named list of candidate FP powers for each predictor.
#' @param zero_vars Named logical vector identifying predictors whose exact
#'   zeros are handled separately during transformation.
#' @param catzero_vars Named logical vector identifying predictors requiring
#'   a zero indicator alongside their transformed positive values.
#' @param spike_vars Named logical vector identifying predictors eligible for
#'   spike-at-zero modelling.
#' @param keep Character vector naming adjustment terms or columns to retain.
#' @param force_max_fp Named logical vector identifying predictors whose
#'   maximum permitted FP degree must be retained.
#' @param group_input_levels Numeric group values in \code{x}, in group order.
#'   If \code{NULL}, uses the sorted observed values.
#' @param group_levels_original Original labels corresponding to
#'   \code{group_input_levels}. If \code{NULL}, uses the numeric group values
#'   as labels.
#' @param term_to_columns Named list mapping adjustment terms to predictor
#'   columns.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{x}}{Adjustment predictor matrix, with optional group dummy
#'     columns.}
#'   \item{\code{original_x}}{Input predictor matrix with recoded group values.}
#'   \item{\code{cat_info}}{Group values, labels, level mappings, and dummy
#'     matrix.}
#'   \item{\code{dummy_names}}{Names of appended group dummy columns, or
#'     \code{NULL} when group main effects are excluded.}
#'   \item{\code{term_to_columns}}{Updated mapping of adjustment terms to
#'     predictor columns.}
#'   \item{\code{updated_params}}{Modelling settings aligned with the adjustment
#'     predictors.}
#' }
#'
#' @keywords internal
#' @noRd
preprocess_data <- function(x, group_var, include_group_var,
                            select, alpha, df, center, acd_vars,
                            fp_powers, zero_vars, catzero_vars, spike_vars,
                            keep, force_max_fp, group_input_levels = NULL,
                            group_levels_original = NULL,
                            term_to_columns) {

  original_names <- colnames(x)

  # Developer note:
  # Internally remap group_var to 0, 1, ..., K - 1. Several downstream helpers
  # build and parse column names using the group value as a suffix, so arbitrary
  # original levels such as 10/20, 0.5/1.5, -1/1, or character labels are unsafe.
  #
  # `input_levels` are the numeric values currently stored in x[, group_var].
  # `original_levels` are the user-facing labels used for printing, summaries,
  # plotting, and prediction metadata.
  group_values_original <- x[, group_var, drop = FALSE]

  input_levels <- if (!is.null(group_input_levels)) {
    as.numeric(group_input_levels)
  } else {
    sort(unique(drop(group_values_original)))
  }

  original_levels <- if (!is.null(group_levels_original)) {
    as.character(group_levels_original)
  } else {
    as.character(input_levels)
  }

  if (length(original_levels) != length(input_levels)) {
    stop(
      "Internal error: `group_levels_original` and `group_input_levels` ",
      "must have the same length.",
      call. = FALSE
    )
  }

  # Keep only group levels that are still present after any row-level filtering
  # such as `subset`. `drop()` removes matrix/data-frame dimensions if present.
  observed_input_levels <- sort(unique(drop(group_values_original)))

  # Identify metadata levels that were known before filtering but are no longer
  # represented in the data. Keeping these would create all-zero group dummies.
  unused_input_levels <- setdiff(input_levels, observed_input_levels)

  # If filtering removed one or more groups, remove the corresponding entries
  # from both the numeric/internal level vector and the original-label vector.
  # This keeps the two vectors aligned before new internal levels are assigned.
  if (length(unused_input_levels) > 0L) {
    keep <- input_levels %in% observed_input_levels
    input_levels <- input_levels[keep]
    original_levels <- original_levels[keep]
  }

  new_levels <- seq_along(input_levels) - 1L

  n_group_levels <- length(new_levels)

  if (isTRUE(include_group_var) && n_group_levels > 6L) {
    warning(
      paste0(
        "`include_group_var = TRUE` will add ",
        n_group_levels - 1L,
        " forced group dummy variables to the adjustment model. ",
        "This may be unstable when group levels are sparse."
      ),
      call. = FALSE
    )
  }

  group_mapped <- new_levels[
    match(drop(group_values_original), input_levels)
  ]

  if (anyNA(group_mapped)) {
    bad <- unique(drop(group_values_original)[is.na(group_mapped)])
    stop(
      paste0(
        "Internal error: group values could not be mapped to stored group levels. ",
        "Unmatched value(s): ",
        paste(bad, collapse = ", "),
        ". This usually indicates that `group_var` was modified before ",
        "group-level remapping."
      ),
      call. = FALSE
    )
  }

  group_values <- matrix(
    group_mapped,
    ncol = 1L,
    dimnames = dimnames(group_values_original)
  )

  x_internal <- x
  x_internal[, group_var] <- group_mapped

  # Developer note:
  # original_x is "original" with respect to later filtering and dummy expansion,
  # but it intentionally contains the internally remapped group_var. The raw
  # group values are stored separately in cat_info$original_values.
  original_x <- x_internal

  # Remove group_var from the adjustment predictor set. Interaction-stage code
  # still has access to the remapped group variable through original_x/cat_info.
  predictor_names <- setdiff(original_names, group_var)
  x_filtered      <- x_internal[, predictor_names, drop = FALSE]

  # Remove the dedicated grouping variable from the adjustment-term lookup but
  # retain every other conceptual term and its surviving raw columns.
  adjustment_term_to_columns <- lapply(term_to_columns, function(cols) {
    intersect(cols, predictor_names)
  })
  adjustment_term_to_columns <- adjustment_term_to_columns[
    lengths(adjustment_term_to_columns) > 0L
  ]

  adjustment_column_to_term <- stats::setNames(
    rep(names(adjustment_term_to_columns), lengths(adjustment_term_to_columns)),
    unlist(adjustment_term_to_columns, use.names = FALSE)
  )

  # Convert keep entries to conceptual adjustment terms. A raw member column of
  # a grouped categorical term therefore retains the complete block.
  keep_filtered <- character(0L)
  if (!is.null(keep)) {
    keep_filtered <- unique(vapply(keep, function(value) {
      if (value %in% names(adjustment_term_to_columns)) {
        value
      } else if (value %in% names(adjustment_column_to_term)) {
        unname(adjustment_column_to_term[[value]])
      } else {
        NA_character_
      }
    }, character(1L)))
    keep_filtered <- keep_filtered[!is.na(keep_filtered)]
  }

  # Developer note:
  # Slice per-variable parameters by name, not by position. This is safer after
  # group_var removal and after formula-interface processing. Required modelling
  # parameters must never be NULL here; a NULL value means the upstream interface
  # failed to expand defaults correctly before calling fit_mfpi().
  slice_param <- function(v, param_name) {
    if (is.null(v)) {
      stop(
        paste0(
          "! Internal error in preprocess_data(): `", param_name, "` is NULL.\n",
          "i This parameter should have been expanded to one value per column ",
          "before calling fit_mfpi()."
        ),
        call. = FALSE
      )
    }

    if (is.null(names(v))) {
      if (length(v) != length(original_names)) {
        stop(
          paste0(
            "! Internal error in preprocess_data(): `", param_name,
            "` is unnamed and has length ", length(v),
            ", but expected length ", length(original_names), "."
          ),
          call. = FALSE
        )
      }
      names(v) <- original_names
    }

    missing_names <- setdiff(predictor_names, names(v))
    if (length(missing_names) > 0L) {
      stop(
        paste0(
          "! Internal error in preprocess_data(): `", param_name,
          "` is missing entries for: ",
          paste(missing_names, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }

    out <- v[predictor_names]
    names(out) <- predictor_names
    out
  }

  updated_params <- list(
    select       = slice_param(select,       "select"),
    alpha        = slice_param(alpha,        "alpha"),
    df           = slice_param(df,           "df"),
    center       = slice_param(center,       "center"),
    acd_vars     = slice_param(acd_vars,     "acd_vars"),
    fp_powers    = slice_param(fp_powers,    "fp_powers"),
    zero_vars    = slice_param(zero_vars,    "zero_vars"),
    catzero_vars = slice_param(catzero_vars, "catzero_vars"),
    spike_vars   = slice_param(spike_vars,   "spike_vars"),
    force_max_fp = slice_param(force_max_fp, "force_max_fp"),
    keep         = keep_filtered
  )

  dummy_names <- NULL

  # Compute group dummies once. They are always stored in cat_info and are also
  # appended to the adjustment matrix when include_group_var = TRUE.
  group_dummies_mat <- create_group_dummies(
    group_values,
    levels = new_levels,
    quiet = TRUE
  )

  # Group dummy names are generated from the pattern paste0(group_var, level),
  # for example "trt1". A user-supplied predictor can legitimately have the same
  # name. If that happens, duplicate column names would corrupt model fitting,
  # adjustment-variable selection, coefficient lookup, and prediction
  # reconstruction. Stop early with a clear error instead of allowing silent
  # name collisions.
  dummy_names_candidate <- colnames(group_dummies_mat)
  dummy_name_collisions <- intersect(dummy_names_candidate, colnames(x_filtered))

  if (length(dummy_name_collisions) > 0L) {
    stop(
      paste0(
        "Generated group dummy column name(s) collide with existing predictor ",
        "column name(s): ",
        paste(dummy_name_collisions, collapse = ", "),
        ". Rename the existing predictor column(s) or the grouping variable."
      ),
      call. = FALSE
    )
  }

  #
  if (include_group_var) {
    dummy_names <- colnames(group_dummies_mat)
    x_filtered  <- cbind(x_filtered, group_dummies_mat)
    n_dummies   <- length(dummy_names)

    make_dummy_param <- function(val) {
      setNames(rep(val, n_dummies), dummy_names)
    }

    dummy_params <- list(
      select       = make_dummy_param(1),     # force group dummies into MFP
      alpha        = make_dummy_param(0),     # no FP degree selection
      df           = make_dummy_param(1L),    # linear binary dummy terms
      center       = make_dummy_param(FALSE), # do not centre dummy terms
      acd_vars     = make_dummy_param(FALSE),
      zero_vars    = make_dummy_param(FALSE),
      catzero_vars = make_dummy_param(FALSE),
      spike_vars   = make_dummy_param(FALSE),
      force_max_fp = make_dummy_param(FALSE),
      fp_powers    = setNames(
        replicate(n_dummies, 1, simplify = FALSE),
        dummy_names
      )
    )

    for (param in names(dummy_params)) {
      updated_params[[param]] <- c(
        updated_params[[param]],
        dummy_params[[param]]
      )
    }

    # Treat all group dummy columns as one conceptual term. They are forced
    # into the adjustment model together and removed together before stage 2.
    adjustment_term_to_columns[[group_var]] <- dummy_names
    updated_params$keep <- unique(c(keep_filtered, group_var))
  }

  list(
    x              = x_filtered,
    original_x     = original_x,
    cat_info       = list(
      values          = group_values,
      original_values = group_values_original,
      original_levels = original_levels,
      input_levels    = input_levels,
      new_levels      = new_levels,
      level_map       = data.frame(
        original = original_levels,
        input    = input_levels,
        internal = new_levels,
        stringsAsFactors = FALSE
      ),
      group_var      = group_var,
      dummies        = group_dummies_mat
    ),
    dummy_names    = dummy_names,
    term_to_columns = adjustment_term_to_columns,
    updated_params = updated_params
  )
}

# -----------------------------------------------------------------------------
# restore_mfpi_adjustment_shifts()
# -----------------------------------------------------------------------------

#' Restore Predictor Shifts in the MFPI Adjustment Model
#'
#' Records the original predictor shifts in the fitted adjustment model so
#' that reported transformations and predictions from new data use the correct
#' values. The predictors are already shifted before fitting, so the adjustment
#' model is initially fitted with zero shifts.
#'
#' @param adjustment_model Fitted \code{mfp2} adjustment model.
#' @param shift Named numeric vector of predictor shifts.
#' @param term_to_columns Named list mapping adjustment terms to predictor
#'   columns.
#'
#' @return The adjustment model with the \code{shift} column of
#'   \code{transformations} updated to the original predictor shifts.
#'
#' @keywords internal
#' @noRd
restore_mfpi_adjustment_shifts <- function(adjustment_model,
                                           shift,
                                           term_to_columns) {
  if (is.null(adjustment_model$transformations) ||
      is.null(rownames(adjustment_model$transformations)) ||
      !"shift" %in% colnames(adjustment_model$transformations) ||
      is.null(shift) || length(shift) == 0L) {
    return(adjustment_model)
  }

  if (is.null(names(shift)) || anyNA(names(shift)) || any(!nzchar(names(shift)))) {
    stop(
      "Internal error: MFPI shift values must be named before restoring adjustment-model metadata.",
      call. = FALSE
    )
  }

  transformation_terms <- rownames(adjustment_model$transformations)

  for (term in transformation_terms) {
    # A singleton conceptual term normally has the same name as its raw column.
    # Prefer that direct lookup because it also covers formula-created singleton
    # terms whose mapping is intentionally simple.
    if (term %in% names(shift)) {
      term_shift <- unname(shift[[term]])
    } else {
      cols <- term_to_columns[[term]]
      cols <- intersect(cols, names(shift))

      # The MFPI grouping term can be represented by internally created dummy
      # columns that have no entry in the original shift vector. Its shift must
      # remain the Stage-1 value (zero), so simply skip terms with no source shift.
      if (length(cols) == 0L) next

      values <- unname(shift[cols])
      if (length(values) > 1L && !all(values == values[[1L]])) {
        stop(
          sprintf(
            "Internal error: grouped adjustment term '%s' has inconsistent MFPI shift values.",
            term
          ),
          call. = FALSE
        )
      }
      term_shift <- values[[1L]]
    }

    if (length(term_shift) != 1L || is.na(term_shift) || !is.finite(term_shift)) {
      stop(
        sprintf(
          "Internal error: adjustment term '%s' does not have a finite MFPI shift.",
          term
        ),
        call. = FALSE
      )
    }

    adjustment_model$transformations[term, "shift"] <- term_shift
  }

  adjustment_model
}


# -----------------------------------------------------------------------------
# fit_adjustment_model() ------------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the MFP Adjustment Model
#'
#' Selects adjustment terms and their fractional-polynomial transformations
#' using \code{fit_mfp()}. Categorical terms are selected jointly.
#'
#' @param x Prepared numeric adjustment predictor matrix, optionally including
#'   group dummy variables.
#' @param y Response vector, response matrix, or \code{survival::Surv()} object.
#' @param weights Numeric vector of observation weights.
#' @param offset Linear-predictor offsets: a numeric vector or, for multinomial
#'   models, a matrix with one column per non-reference outcome category.
#' @param cycles Maximum number of MFP backfitting cycles.
#' @param family Resolved model-family specification.
#' @param family_string Character string giving the canonical family name.
#' @param criterion Adjustment selection criterion: \code{"pvalue"},
#'   \code{"aic"}, or \code{"bic"}.
#' @param updated_params Named list of modelling settings returned by
#'   \code{preprocess_data()}, including selection levels, degrees of freedom,
#'   candidate powers, transformation options, and retention settings.
#' @param xorder Order of adjustment terms during MFP backfitting:
#'   \code{"ascending"}, \code{"descending"}, or \code{"original"}.
#' @param ties Method for handling tied event times in Cox models.
#' @param strata Optional prepared survival-model stratification information.
#' @param nocenter Numeric values excluded from automatic centring by
#'   \code{survival::coxph()}.
#' @param min_saz_prop Minimum proportion required in both the zero and
#'   positive components of a spike-at-zero adjustment variable.
#' @param use_ftest Whether to use F-tests instead of likelihood-ratio
#'   chi-square tests for Gaussian p-value comparisons. Ignored for other
#'   families.
#' @param control List of model-fitting control parameters.
#' @param force_max_fp Retained internal argument. The settings used for
#'   fitting are read from \code{updated_params$force_max_fp}.
#' @param scale Named numeric vector of predictor scale factors.
#' @param has_offset Logical indicating whether an offset was supplied.
#' @param term_to_columns Named list mapping adjustment terms to predictor
#'   columns.
#' @param verbose Whether to print fitting progress.
#' @param fitter Model-fitting backend.
#'
#' @details
#' Predictors have already been shifted and scaled. Zero shifts are passed to
#' \code{fit_mfp()} to avoid applying shifts twice. Scale factors are retained
#' so that selected terms are constructed on the shifted, unscaled predictor
#' scale. The caller restores the original shifts in the fitted model.
#'
#' Variables marked in \code{updated_params$force_max_fp} retain their maximum
#' permitted FP degree, with powers selected within that degree. For eligible
#' spike-at-zero terms, both the zero indicator and the maximum-degree
#' continuous component are retained.
#'
#' @return The adjustment-model results returned by \code{fit_mfp()}, including
#'   the selected-term table \code{fp_terms} and the underlying fitted model.
#'
#' @keywords internal
#' @noRd
fit_adjustment_model <- function(x, y, weights, offset, cycles, family,
                                 family_string, criterion, updated_params,
                                 xorder, ties, strata, nocenter,
                                 min_saz_prop, use_ftest, control,
                                 force_max_fp, scale, has_offset,
                                 term_to_columns,
                                 verbose = FALSE, fitter = "base") {
  n_vars    <- ncol(x)
  x_names   <- colnames(x)

  # Build raw-column scale values first, then collapse all modelling settings
  # to one value per conceptual term for the grouped MFP core. Group dummy and
  # categorical contrast columns use scale 1.
  scale_vec <- setNames(rep(1, n_vars), x_names)
  if (!is.null(scale)) {
    shared <- intersect(x_names, names(scale))
    scale_vec[shared] <- scale[shared]
  }

  collapse <- function(values, setting) {
    collapse_option_to_terms(values, term_to_columns, setting)
  }

  df_term <- collapse(updated_params$df, "df")
  center_term <- collapse(updated_params$center, "center")
  select_term <- collapse(updated_params$select, "select")
  alpha_term <- collapse(updated_params$alpha, "alpha")
  acd_term <- collapse(updated_params$acd_vars, "acd_vars")
  zero_term <- collapse(updated_params$zero_vars, "zero_vars")
  catzero_term <- collapse(updated_params$catzero_vars, "catzero_vars")
  spike_term <- collapse(updated_params$spike_vars, "spike_vars")
  force_term <- collapse(updated_params$force_max_fp, "force_max_fp")
  scale_term <- collapse(scale_vec, "scale")

  powers_term <- lapply(names(term_to_columns), function(term) {
    cols <- term_to_columns[[term]]
    if (term_uses_column_mapping(term, cols)) {
      1
    } else {
      updated_params$fp_powers[[cols[[1L]]]]
    }
  })
  names(powers_term) <- names(term_to_columns)

  fit_mfp(
    x                = x,
    y                = y,
    weights          = weights,
    offset           = offset,
    cycles           = cycles,
    method           = ties,
    strata           = strata,
    nocenter         = nocenter,
    scale            = scale_term,
    shift            = stats::setNames(rep(0, length(term_to_columns)), names(term_to_columns)),
    ftest            = use_ftest,
    control          = control,
    family           = family,
    family_string    = family_string,
    fitter           = fitter,
    criterion        = criterion,
    xorder           = xorder,
    df               = df_term,
    center           = center_term,
    select           = select_term,
    alpha            = alpha_term,
    keep             = updated_params$keep,
    powers           = powers_term,
    acdx             = acd_term,
    zero             = zero_term,
    catzero          = catzero_term,
    spike            = spike_term,
    min_saz_prop     = min_saz_prop,
    saz_pre_resolved = TRUE,
    force_max_fp     = force_term,
    term_to_columns  = term_to_columns,
    has_offset       = has_offset,
    verbose          = verbose
  )
}

# -----------------------------------------------------------------------------
# format_type_label() ---------------------------------------------------------
# -----------------------------------------------------------------------------

#' Format Functional Form Labels
#'
#' Converts \code{"linear"}, \code{"fp1"}, and \code{"fp2"} to
#' \code{"Linear"}, \code{"FP1"}, and \code{"FP2"}, respectively.
#'
#' @param type Character vector of functional form names.
#'
#' @return Character vector of display labels.
#'
#' @keywords internal
#' @noRd
format_type_label <- function(type) {
  ifelse(
    tolower(type) == "linear",
    "Linear",
    toupper(type)   # fp1 -> FP1, fp2 -> FP2
  )
}

# -----------------------------------------------------------------------------
# format_candidate_table() ----------------------------------------------------
# -----------------------------------------------------------------------------

#' Print an Interaction Comparison Table
#'
#' Prints comparison metrics for one interaction variable. Displayed columns
#' depend on the selection criterion:
#' \itemize{
#'   \item \code{"pvalue"}: interaction-model deviance, deviance difference,
#'     interaction degrees of freedom, and p-value.
#'   \item \code{"aic"}: main-effects and interaction-model AIC values and
#'     their difference.
#'   \item \code{"bic"}: main-effects and interaction-model BIC values and
#'     their difference.
#' }
#'
#' Also displays FP powers for the main-effects and interaction models.
#' For FLEX4, interaction powers are shown separately for each group.
#' AIC and BIC differences are calculated as main-effects minus interaction
#' values, so positive differences favour the interaction model.
#'
#' @param rows Named list of one-row data frames containing comparison metrics
#'   and FP powers for the evaluated interaction form.
#' @param criterion Selection criterion: \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}.
#' @param digits Number of decimal places for displayed statistics.
#'
#' @return Invisibly returns \code{NULL} after printing the table.
#'
#' @keywords internal
#' @noRd
format_candidate_table <- function(rows, criterion, digits,  best_type = NULL) {

  # Format a single power entry as a parenthesised string. Always returns a
  # length-1 character string ("." for missing/empty).
  # For flex1-3, `pow` is a numeric vector e.g. c(1, 2).
  # For flex4, `pow` is a list of per-group vectors e.g. list(c(1, 2), c(0.5)).
  format_powers <- function(pow) {
    if (is.null(pow) || length(pow) == 0L) return(".")
    if (is.atomic(pow) && all(is.na(pow))) return(".")
    if (is.list(pow)) {
      paste(
        vapply(pow, function(p) {
          if (is.null(p) || length(p) == 0L) "()"
          else paste0("(", paste(p, collapse = ", "), ")")
        }, character(1L)),
        collapse = ", "
      )
    } else {
      paste0("(", paste(pow, collapse = ", "), ")")
    }
  }

  format_row <- function(metrics) {
    # Safe extraction helpers
    safe_decimal <- function(val) {
      if (is.null(val) || length(val) == 0L || is.na(val[1L])) return("NA")
      format_print_decimal(val[1L], digits)
    }
    safe_pvalue <- function(val) {
      if (is.null(val) || length(val) == 0L || is.na(val[1L])) return("NA")
      format_print_pvalue(val[1L], digits)
    }
    safe_df <- function(val) {
      if (is.null(val) || length(val) == 0L || is.na(val[1L])) return("NA")
      format(val[1L], trim = TRUE, scientific = FALSE)
    }
    safe_chr <- function(val) {
      if (is.null(val) || length(val) == 0L) return("")
      as.character(val[1L])
    }

    type_raw   <- safe_chr(metrics$type)
    type_label <- if (nzchar(type_raw)) format_type_label(type_raw) else "?"

    main_pow <- if (length(metrics$fp_powers_main) > 0L)
      metrics$fp_powers_main[[1L]] else NULL
    int_pow  <- if (length(metrics$fp_powers_int)  > 0L)
      metrics$fp_powers_int[[1L]]  else NULL

    # Base columns: Type and powers always shown
    base_cols <- data.frame(
      Type       = type_label,
      main_power = format_powers(main_pow),
      int_power  = format_powers(int_pow),
      stringsAsFactors = FALSE
    )

    # Criterion-specific columns inserted between Type and powers
    crit_cols <- switch(
      criterion,
      "pvalue" = data.frame(
        deviance_int  = safe_decimal(metrics$deviance_int),
        deviance_diff = safe_decimal(metrics$deviance_diff),
        df            = safe_df(metrics$df_int),
        pvalue        = safe_pvalue(metrics$pvalue),
        stringsAsFactors = FALSE
      ),
      "aic"    = data.frame(
        AIC_main        = safe_decimal(metrics$AIC_main),
        AIC_interaction = safe_decimal(metrics$AIC_interaction),
        dAIC            = safe_decimal(metrics$AIC_main_minus_int),
        stringsAsFactors = FALSE
      ),
      "bic"    = data.frame(
        BIC_main        = safe_decimal(metrics$BIC_main),
        BIC_interaction = safe_decimal(metrics$BIC_interaction),
        dBIC            = safe_decimal(metrics$BIC_main_minus_int),
        stringsAsFactors = FALSE
      )
    )

    # Final column order: Type, criterion metrics, then powers
    cbind(base_cols[, "Type", drop = FALSE],
          crit_cols,
          base_cols[, c("main_power", "int_power"), drop = FALSE])
  }

  tbl <- do.call(rbind, lapply(rows, format_row))
  print(tbl, row.names = FALSE, right = FALSE)
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# evaluate_interactions() -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Evaluate Separate MFPI Interactions
#'
#' Evaluates each specified interaction with the grouping variable, adjusting
#' for the selected adjustment terms other than the variable being tested.
#' Continuous, binary, and categorical interaction terms are supported.
#'
#' Each variable uses its resolved interaction form; no selection between
#' linear, FP1, and FP2 forms is performed. Results are retained for all
#' evaluated interactions, with selected interactions also stored separately.
#'
#' @param y Response vector, response matrix, or \code{survival::Surv()} object.
#' @param processed_data List returned by \code{preprocess_data()}.
#' @param selected_vars Character vector naming selected adjustment terms.
#' @param adj_fp_powers Named list of selected adjustment powers, including
#'   structural \code{NA} entries for ACD terms.
#' @param adj_zero Named logical vector of fitted zero-handling flags.
#' @param adj_catzero Named logical vector of fitted zero-indicator flags.
#' @param adj_spike Named logical vector of fitted spike-at-zero flags.
#' @param adj_spike_decision Named numeric vector of fitted spike-at-zero
#'   decisions: \code{1} for continuous and binary components, \code{2} for
#'   the continuous component only, and \code{3} for the binary component only.
#' @param adj_acd_parameter Fitted ACD transformation parameters for adjustment
#'   terms.
#' @param interaction_vars Character vector naming terms to test for
#'   interaction with the grouping variable.
#' @param interaction_forms Named character vector of resolved interaction
#'   forms: \code{"linear"}, \code{"fp1"}, or \code{"fp2"}.
#' @param interaction_specs Named list describing each interaction term's type
#'   and predictor columns.
#' @param flex Interaction procedure: \code{"flex1"}, \code{"flex2"},
#'   \code{"flex3"}, or \code{"flex4"}.
#' @param weights Numeric vector of observation weights.
#' @param offset Linear-predictor offsets: a numeric vector or, for multinomial
#'   models, a matrix with one column per non-reference outcome category.
#' @param xorder Order of terms during MFP backfitting:
#'   \code{"ascending"}, \code{"descending"}, or \code{"original"}.
#' @param ties Method for handling tied event times in Cox models.
#' @param strata Optional survival-model stratification information.
#' @param use_ftest Whether to use F-tests instead of likelihood-ratio
#'   chi-square tests for Gaussian p-value comparisons. Ignored for other
#'   families.
#' @param control List of model-fitting control parameters.
#' @param nocenter Numeric values excluded from automatic centring by
#'   \code{survival::coxph()}.
#' @param family Resolved model-family specification.
#' @param family_string Character string giving the canonical family name.
#' @param fp_powers Named list of candidate FP powers for each predictor.
#' @param cycles Maximum number of MFP backfitting cycles.
#' @param criterion Interaction selection criterion: \code{"pvalue"},
#'   \code{"aic"}, or \code{"bic"}.
#' @param digits Number of decimal places for printed numerical results.
#' @param scale Named numeric vector of predictor scale factors.
#' @param shift Named numeric vector of predictor shifts.
#' @param center_type Centring method for interaction transformations:
#'   \code{"grand"} or \code{"group"}.
#' @param show_models Whether to print coefficient tables for selected
#'   interaction models when \code{verbose = TRUE}.
#' @param group_var Character string naming the grouping variable.
#' @param p_interact Significance threshold for p-value selection.
#' @param min_improvement Minimum AIC or BIC improvement required for selection.
#' @param min_saz_prop Minimum proportion required in both components of a
#'   spike-at-zero variable.
#' @param xadj Optional matrix of transformed adjustment predictors, used when
#'   \code{skip_adjustment = TRUE}.
#' @param skip_adjustment Whether to use \code{xadj} without rebuilding the
#'   adjustment transformations.
#' @param p_adjust_method Method passed to \code{stats::p.adjust()} for
#'   interaction p-values. Used only under \code{criterion = "pvalue"}.
#' @param has_offset Logical indicating whether an offset was supplied.
#' @param quiet Internal diagnostic-output control.
#' @param verbose Whether to print progress messages.
#' @param fitter Model-fitting backend.
#'
#' @section Selection:
#' Each interaction is compared with its corresponding main-effects model.
#' Under \code{criterion = "pvalue"}, selection requires a p-value strictly
#' below \code{p_interact}. If requested, p-values are adjusted across all
#' interaction variables before selection.
#'
#' Under \code{"aic"} or \code{"bic"}, selection requires the main-effects
#' criterion minus the interaction criterion to exceed \code{min_improvement}.
#' P-value adjustment does not affect AIC or BIC selection.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{best_model_metrics}, \code{all_model_metrics}}{Data frames of
#'     comparison metrics for selected and evaluated interactions, respectively.}
#'   \item{\code{best_interaction_model},
#'     \code{all_interaction_models}}{Named lists of fitted models for selected
#'     and evaluated interactions, respectively.}
#'   \item{\code{center_vals_list}}{Centring constants for selected interactions.}
#'   \item{\code{var_winners}}{Fitted results for each interaction variable,
#'     regardless of selection.}
#' }
#' The list also contains group labels and mappings, the family, observation
#' count, selection settings, and predictor shifts and scales.
#'
#' @keywords internal
#' @noRd
evaluate_interactions <- function(y, processed_data, selected_vars,
                                  adj_fp_powers,
                                  adj_zero, adj_catzero, adj_spike,
                                  adj_spike_decision,
                                  interaction_vars, interaction_forms,
                                  interaction_specs, flex,
                                  weights, offset, xorder, ties, strata,
                                  use_ftest, control, nocenter, family,
                                  family_string, fp_powers, cycles, criterion,
                                  digits, group_var, show_models, p_interact,
                                  min_improvement, min_saz_prop,
                                  center_type = c("grand", "group"),
                                  scale = NULL, shift = NULL,
                                  shift_automatic = NULL,
                                  adj_acd_parameter = NULL,
                                  xadj = NULL, skip_adjustment = FALSE,
                                  p_adjust_method = "none", has_offset,
                                  quiet = FALSE, verbose = FALSE, fitter = "base") {

  # Normalise the centering strategy. Interaction FP bases can either use one
  # grand centering vector shared across groups or group-specific centering.
  center_type <- match.arg(center_type)

  # ---------------------------------------------------------------------------
  # Validate adjustment-matrix control arguments
  # ---------------------------------------------------------------------------
  # `skip_adjustment = TRUE` means an already prepared adjustment matrix may be
  # supplied via `xadj`; otherwise the adjustment matrix is built from selected
  # adjustment variables below.
  if (!is.logical(skip_adjustment) || length(skip_adjustment) != 1L) {
    stop("`skip_adjustment` must be a single logical value.", call. = FALSE)
  }

  if (skip_adjustment && !is.null(xadj)) {
    if (!is.matrix(xadj) || !is.numeric(xadj)) {
      stop(
        "`xadj` must be a numeric matrix when `skip_adjustment = TRUE`.",
        call. = FALSE
      )
    }
  }

  # Main processed predictor matrix and group metadata used throughout the
  # interaction-evaluation stage.
  x        <- processed_data$original_x
  nobs     <- nrow(x)
  cat_info <- processed_data$cat_info

  # Map the user-facing interaction form to the FP degree expected by flex_fit().
  # linear -> degree 0, fp1 -> degree 1, fp2 -> degree 2.
  degree_lookup <- c(linear = 0L, fp1 = 1L, fp2 = 2L)

  # ---------------------------------------------------------------------------
  # Pre-allocate result containers
  # ---------------------------------------------------------------------------
  # Each variable contributes exactly one candidate form because MFPI now uses
  # pre-specified `interaction_forms`; there is no data-driven search over linear,
  # FP1, and FP2 inside this loop.
  n_vars  <- length(interaction_vars)
  n_types <- 1L

  best_metrics_list      <- vector("list", n_vars)
  all_metrics_list       <- vector("list", n_vars * n_types)
  best_interaction_model <- vector("list", n_vars)
  all_interaction_models <- vector("list", n_vars)
  center_vals_list       <- vector("list", n_vars)

  names(best_interaction_model) <- interaction_vars
  names(all_interaction_models) <- interaction_vars
  names(center_vals_list)       <- interaction_vars

  # Store the best candidate for every variable, even if it is not retained.
  # This is needed later for multiplicity-adjusted p-value decisions.
  var_winners <- vector("list", n_vars)
  names(var_winners) <- interaction_vars

  # P-value selection is deferred only when a real multiplicity adjustment is
  # requested. If p_adjust_method = "none", raw p-values are final immediately.
  adjusting_pvals <- criterion == "pvalue" && p_adjust_method != "none"

  # Shared AIC/BIC metadata used by the per-variable winner logic, verbose
  # selected/not-selected messages, and the final AIC/BIC summary. Keeping this
  # mapping in one place prevents the metric column names and display labels from
  # drifting across branches.
  ic_col <- NULL
  ic_label <- NULL

  if (criterion == "aic") {
    ic_col   <- "AIC_main_minus_int"
    ic_label <- "dAIC"
  } else if (criterion == "bic") {
    ic_col   <- "BIC_main_minus_int"
    ic_label <- "dBIC"
  }

  # Write counters for the pre-allocated result lists.
  best_idx <- 0L
  all_idx  <- 0L

  # Filled only when p-value adjustment is requested. In that case it becomes
  # the authoritative all-model metric table returned by this function.
  all_candidate_metrics <- NULL

  # ---------------------------------------------------------------------------
  # Build the adjustment-block cache once
  # ---------------------------------------------------------------------------
  # The adjustment set is fixed across the interaction stage except that the
  # current tested variable must be removed from its own adjustment set. Instead
  # of rebuilding the full transformed adjustment matrix inside every cont_var
  # iteration, transform each selected adjustment variable once and reuse the
  # cached source-variable blocks.
  group_dummy_names <- colnames(processed_data$cat_info$dummies)

  adjustment_cache <- mfpi_build_adjustment_block_cache(
    x                  = x,
    selected_vars      = selected_vars,
    group_dummy_names  = group_dummy_names,
    group_dummy_term   = cat_info$group_var,
    term_to_columns    = processed_data$term_to_columns,
    skip_adjustment    = skip_adjustment,
    scale              = scale,
    adj_fp_powers      = adj_fp_powers,
    center             = processed_data$updated_params$center,
    acd_vars           = processed_data$updated_params$acd_vars,
    adj_zero           = adj_zero,
    adj_catzero        = adj_catzero,
    adj_spike          = adj_spike,
    adj_spike_decision = adj_spike_decision,
    adj_acd_parameter  = adj_acd_parameter
  )

  # Source-variable names for selected adjustment variables after group dummies
  # have been removed. This is used both for verbose display and xadj assembly.
  selected_adj_vars <- adjustment_cache$selected_adj_vars

  # ===========================================================================
  # Main loop: evaluate one pre-specified interaction candidate per cont_var
  # ===========================================================================
  for (var_name in interaction_vars) {
    var_idx <- match(var_name, interaction_vars)

    res <- evaluate_interaction_for_variable(
      var_name               = var_name,
      var_idx                = var_idx,
      n_vars                 = n_vars,
      x                      = x,
      y                      = y,
      cat_info               = cat_info,
      adjustment_cache       = adjustment_cache,
      selected_adj_vars      = selected_adj_vars,
      skip_adjustment        = skip_adjustment,
      xadj                   = xadj,
      interaction_forms      = interaction_forms,
      interaction_specs      = interaction_specs,
      degree_lookup          = degree_lookup,
      processed_data         = processed_data,
      criterion              = criterion,
      ties                   = ties,
      family                 = family,
      family_string          = family_string,
      fitter                 = fitter,
      use_ftest              = use_ftest,
      center_type            = center_type,
      xorder                 = xorder,
      weights                = weights,
      offset                 = offset,
      strata                 = strata,
      control                = control,
      nocenter               = nocenter,
      cycles                 = cycles,
      min_saz_prop           = min_saz_prop,
      flex                   = flex,
      scale                  = scale,
      shift                  = shift,
      shift_automatic        = shift_automatic,
      has_offset             = has_offset,
      ic_col                 = ic_col,
      ic_label               = ic_label,
      adjusting_pvals        = adjusting_pvals,
      p_interact             = p_interact,
      min_improvement        = min_improvement,
      verbose                = verbose,
      quiet                  = quiet,
      show_models            = show_models,
      digits                 = digits
    )

    # Store this variable's candidate metrics rows in the full all-model
    # metrics list, mirroring the original per-candidate all_idx bookkeeping
    # (currently always one row per variable, since there is exactly one
    # pre-specified candidate form per variable).
    for (interaction_type in names(res$candidate_metrics)) {
      all_idx <- all_idx + 1L
      all_metrics_list[[all_idx]] <- res$candidate_metrics[[interaction_type]]
    }

    # Store the per-variable winner for later summaries and, when needed,
    # multiplicity-adjusted final selection.
    var_winners[[var_name]] <- list(
      fit         = res$best_fit,
      metric      = res$best_metric,
      type        = res$best_type,
      score       = res$best_score,
      center_vals = if (!is.null(res$best_fit)) res$best_fit$center_vals else NULL
    )

    # Immediate final selection when p-values are not being adjusted. The
    # decision is stored here but is not printed until the Step 3 summary.
    if (!adjusting_pvals && isTRUE(res$retained)) {
      best_idx <- best_idx + 1L
      best_metrics_list[[best_idx]]      <- res$best_metric
      best_interaction_model[[var_name]] <-
        res$best_fit$test_results$interaction_model
      center_vals_list[[var_name]]       <- res$best_fit$center_vals
    }

    # Store the fitted interaction model for this variable directly, regardless
    # of whether it was selected. Each variable has exactly one pre-specified
    # form, so exposing the legacy candidate-list layer would make the public
    # object one level deeper than documented.
    interaction_type <- interaction_forms[[var_name]]
    interaction_model <- res$candidate_models[[interaction_type]]

    if (is.null(interaction_model)) {
      stop(
        sprintf(
          paste0(
            "Internal error: no interaction model was stored for '%s' ",
            "using form '%s'."
          ),
          var_name,
          interaction_type
        ),
        call. = FALSE
      )
    }

    all_interaction_models[[var_name]] <- interaction_model

    # Replace the all-metrics row with the canonical best_metric object where
    # available, preserving normalized reporting columns.
    if (!is.null(res$best_metric)) {
      all_metrics_list[[all_idx]] <- res$best_metric
    }
  }

  # ===========================================================================
  # Step 3 final interaction summary
  # ===========================================================================
  # Final selection decisions are reported only after all variables have been
  # evaluated. Step 2 deliberately reports observed statistics without a
  # provisional selected/not-selected label.
  if (verbose && criterion %in% c("aic", "bic") && !adjusting_pvals) {
    print_interaction_step3_summary(
      var_winners            = var_winners,
      interaction_vars       = interaction_vars,
      mode                   = "ic",
      ic_col                 = ic_col,
      ic_label               = ic_label,
      min_improvement        = min_improvement,
      digits                 = digits,
      show_models            = show_models,
      best_interaction_model = best_interaction_model
    )
  }

  if (verbose && criterion == "pvalue" && !adjusting_pvals) {
    print_interaction_step3_summary(
      var_winners            = var_winners,
      interaction_vars       = interaction_vars,
      mode                   = "pvalue",
      p_interact             = p_interact,
      p_adjust_method        = p_adjust_method,
      digits                 = digits,
      show_models            = show_models,
      best_interaction_model = best_interaction_model
    )
  }

  # ===========================================================================
  # Phase 2: p-value adjustment and final p-value decisions
  # ===========================================================================
  # When multiplicity adjustment is requested, raw p-values from all planned
  # tests must be collected before any final retain/drop decisions are made.
  if (adjusting_pvals) {
    finalized <- finalize_adjusted_pvalue_selection(
      var_winners      = var_winners,
      all_metrics_list = all_metrics_list,
      all_idx          = all_idx,
      interaction_vars = interaction_vars,
      n_vars           = n_vars,
      n_types          = n_types,
      p_adjust_method  = p_adjust_method,
      p_interact       = p_interact
    )

    var_winners            <- finalized$var_winners
    all_candidate_metrics  <- finalized$all_candidate_metrics
    best_idx               <- finalized$best_idx
    best_metrics_list      <- finalized$best_metrics_list
    best_interaction_model <- finalized$best_interaction_model
    center_vals_list       <- finalized$center_vals_list

    # Final adjusted p-value interaction summary.
    if (verbose) {
      print_interaction_step3_summary(
        var_winners            = var_winners,
        interaction_vars       = interaction_vars,
        mode                   = "pvalue",
        p_interact             = p_interact,
        p_adjust_method        = p_adjust_method,
        digits                 = digits,
        show_models            = show_models,
        best_interaction_model = best_interaction_model
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Trim pre-allocated containers
  # ---------------------------------------------------------------------------
  # When p-value adjustment was used, all_candidate_metrics is already the final
  # all-model metric table, so all_metrics_list does not need to be rebound for
  # the returned all_model_metrics.
  best_metrics_list <- best_metrics_list[seq_len(best_idx)]

  if (!adjusting_pvals) {
    all_metrics_list <- all_metrics_list[seq_len(all_idx)]
  }

  # Remove variables that were not finally retained.
  best_interaction_model <- Filter(Negate(is.null), best_interaction_model)

  # ---------------------------------------------------------------------------
  # Build returned metric tables
  # ---------------------------------------------------------------------------
  combined_metrics <- bind_metric_rows(best_metrics_list)

  if (nrow(combined_metrics) > 0L && "type" %in% names(combined_metrics)) {
    combined_metrics <- combined_metrics[
      ,
      c("type", setdiff(names(combined_metrics), "type")),
      drop = FALSE
    ]
  }

  all_model_metrics <- if (adjusting_pvals && !is.null(all_candidate_metrics)) {
    all_candidate_metrics
  } else {
    bind_metric_rows(all_metrics_list)
  }

  # Keep a stable display/storage order for all-model metric columns.
  if (nrow(all_model_metrics) > 0L) {
    metric_cols <- setdiff(names(all_model_metrics), c("variable", "type"))
    all_model_metrics <- all_model_metrics[
      ,
      c("type", "variable", metric_cols),
      drop = FALSE
    ]
  }

  # ---------------------------------------------------------------------------
  # Return MFPI interaction-evaluation result object
  # ---------------------------------------------------------------------------
  list(
    best_model_metrics       = combined_metrics,
    all_model_metrics        = all_model_metrics,
    best_interaction_model   = best_interaction_model,
    all_interaction_models   = all_interaction_models,
    center_vals_list         = center_vals_list,
    group_var                = group_var,
    show_models              = show_models,
    group_levels_new         = cat_info$new_levels,
    group_levels_original    = cat_info$original_levels,
    group_level_map          = cat_info$level_map,
    flex                     = flex,
    criterion                = criterion,
    family                   = family,
    nobs                     = nobs,
    p_adjust_method          = p_adjust_method,
    scale                    = scale,
    shift                    = shift,
    var_winners              = var_winners
  )
}

# -----------------------------------------------------------------------------
# evaluate_interactions() helpers --------------------------------------------
# -----------------------------------------------------------------------------

#' Evaluate One MFPI Interaction
#'
#' Fits the resolved interaction form for one continuous, binary, or
#' categorical term and compares it with the corresponding main-effects model.
#' Makes the selection decision immediately unless p-value adjustment is
#' requested, in which case selection is deferred.
#'
#' @inheritParams evaluate_interactions
#'
#' @param var_name Character string naming the interaction term.
#' @param var_idx,n_vars Current variable index and total number of interaction
#'   variables, used in progress messages.
#' @param x Prepared numeric predictor matrix.
#' @param cat_info Group metadata from \code{preprocess_data()}.
#' @param adjustment_cache Cached transformed adjustment terms.
#' @param selected_adj_vars Character vector naming selected adjustment terms,
#'   excluding group dummy variables.
#' @param degree_lookup Named integer vector mapping interaction forms to
#'   FP degrees: \code{linear = 0}, \code{fp1 = 1}, and \code{fp2 = 2}.
#' @param ic_col,ic_label Metric column name and display label for the active
#'   information criterion, or \code{NULL} for p-value selection.
#' @param adjusting_pvals Whether selection is deferred until interaction
#'   p-values have been adjusted across variables.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{candidate_metrics}}{Named list containing the comparison
#'     metrics for the evaluated form.}
#'   \item{\code{candidate_models}}{Named list containing the fitted
#'     interaction model.}
#'   \item{\code{best_fit}, \code{best_metric}, \code{best_type},
#'     \code{best_score}}{The fitted result, comparison metrics, form, and
#'     comparison score for the evaluated interaction.}
#'   \item{\code{retained}}{Selection decision, or \code{NA} when selection
#'     is deferred for p-value adjustment.}
#' }
#'
#' @keywords internal
#' @noRd
evaluate_interaction_for_variable <- function(var_name,
                                              var_idx,
                                              n_vars,
                                              x,
                                              y,
                                              cat_info,
                                              adjustment_cache,
                                              selected_adj_vars,
                                              skip_adjustment,
                                              xadj,
                                              interaction_forms,
                                              interaction_specs,
                                              degree_lookup,
                                              processed_data,
                                              criterion,
                                              ties,
                                              family,
                                              family_string,
                                              use_ftest,
                                              center_type,
                                              xorder,
                                              weights,
                                              offset,
                                              strata,
                                              control,
                                              nocenter,
                                              cycles,
                                              min_saz_prop,
                                              flex,
                                              scale,
                                              shift,
                                              shift_automatic = NULL,
                                              has_offset,
                                              ic_col,
                                              ic_label,
                                              adjusting_pvals,
                                              p_interact,
                                              min_improvement,
                                              verbose,
                                              quiet,
                                              show_models,
                                              digits,
                                              fitter = "base") {

  # -----------------------------------------------------------------------
  # Verbose per-variable header
  # -----------------------------------------------------------------------
  if (verbose) {
    rule_thin <- strrep("-", 70)

    mfp2_message("")
    mfp2_message(rule_thin)
    mfp2_message(
      sprintf(
        "  [%d/%d] Variable: '%s' x '%s'",
        var_idx, n_vars, var_name, cat_info$group_var
      )
    )

    # Adjustment variables for this variable-specific test. The group dummies
    # have already been excluded in the cache; here we only exclude var_name
    # itself if it was selected by the adjustment model.
    adj_vars <- setdiff(selected_adj_vars, var_name)

    if (length(adj_vars) > 0L) {
      mfp2_message("  Adjusting for: ", paste(adj_vars, collapse = ", "))
    } else {
      mfp2_message("  Adjusting for: (none)")
    }
    mfp2_message(rule_thin)
  }

  # Track the winning candidate for the current variable. Since there is only
  # one pre-specified form per variable, this mostly records successful fit
  # state and the relevant score.
  best_fit    <- NULL
  best_metric <- NULL
  best_type   <- NULL
  best_score  <- Inf

  # The form for this variable is pre-specified by interaction_forms.
  active_type <- interaction_forms[[var_name]]
  interaction_spec <- interaction_specs[[var_name]]
  if (is.null(interaction_spec)) {
    stop(
      paste0("Internal error: no interaction specification is stored for `",
             var_name, "`."),
      call. = FALSE
    )
  }

  # Candidate containers are still list-based for consistency with the older
  # multi-form structure and with downstream storage conventions.
  n_types <- 1L
  candidate_models    <- vector("list", n_types)
  candidate_metrics   <- vector("list", n_types)
  candidate_full_fits <- vector("list", n_types)

  if (verbose) {
    verbose_rows <- vector("list", n_types)
    names(verbose_rows) <- active_type
  }

  # Assemble the transformed adjustment matrix for this variable by combining
  # cached blocks and excluding the current source-variable block if needed.
  xadj_current <- mfpi_make_xadj_current(
    var_name         = var_name,
    adjustment_cache = adjustment_cache,
    skip_adjustment  = skip_adjustment,
    xadj             = xadj
  )

  # -------------------------------------------------------------------------
  # Fit the single pre-specified interaction form
  # -------------------------------------------------------------------------
  for (interaction_type in active_type) {
    degree <- degree_lookup[[interaction_type]]

    fit_result <- mfpi_with_rank_hint(flex_fit(
      x                     = x,
      y                     = y,
      cont_var              = var_name,
      interaction_columns   = interaction_spec$columns,
      interaction_kind      = interaction_spec$kind,
      group_var             = cat_info$group_var,
      group_dummies         = cat_info$dummies,
      xadj                  = xadj_current,
      criterion             = criterion,
      ties                  = ties,
      degree                = degree,
      family                = family,
      family_string         = family_string,
      fitter                = fitter,
      fp_cand               = if (isTRUE(interaction_spec$discrete)) 1 else
        processed_data$updated_params$fp_powers[[var_name]],
      use_ftest             = use_ftest,
      center                = if (isTRUE(interaction_spec$discrete)) FALSE else
        processed_data$updated_params$center[var_name],
      center_type           = center_type,
      xorder                = xorder,
      weights               = weights,
      offset                = offset,
      strata                = strata,
      control               = control,
      nocenter              = nocenter,
      cycles                = cycles,
      zero_var              = if (isTRUE(interaction_spec$discrete)) FALSE else
        processed_data$updated_params$zero_vars[var_name],
      spike_var             = FALSE,
      min_saz_prop          = min_saz_prop,
      flex          = flex,
      scale_var     = if (isTRUE(interaction_spec$discrete)) 1 else if (!is.null(scale)) unname(scale[var_name]) else 1,
      shift_var     = if (isTRUE(interaction_spec$discrete)) 0 else if (!is.null(shift)) unname(shift[var_name]) else 0,
      run_test      = TRUE,
      has_offset    = has_offset
    ),
      var_name         = var_name,
      interaction_spec = interaction_spec,
      interaction_form = interaction_type,
      x                = x,
      scale            = scale,
      shift            = shift,
      shift_automatic  = shift_automatic,
      zero_var         = isTRUE(processed_data$updated_params$zero_vars[var_name])
    )

    # -----------------------------------------------------------------------
    # Extract and annotate model-comparison metrics
    # -----------------------------------------------------------------------
    metrics          <- fit_result$test_results$evaluation_metrics
    metrics$variable <- var_name
    metrics$type     <- interaction_type

    # Normalise FP power columns for reporting. This keeps metric tables
    # readable without changing the fitted interaction model object.
    metrics <- normalise_metric_fp_powers(
      metrics      = metrics,
      variable     = var_name,
      group_levels = cat_info$original_levels
    )

    candidate_metrics[[interaction_type]]   <- metrics
    candidate_full_fits[[interaction_type]] <- fit_result

    if (verbose) {
      verbose_rows[[interaction_type]] <- metrics
    }

    # Store the fitted interaction model regardless of whether it is retained.
    candidate_models[[interaction_type]] <-
      fit_result$test_results$interaction_model

    # For p-value selection, the only candidate is the per-variable winner if
    # its p-value is finite/non-missing.
    if (criterion == "pvalue") {
      pval <- metrics$pvalue[1L]
      if (!is.na(pval)) {
        best_fit    <- fit_result
        best_metric <- metrics
        best_type   <- interaction_type
        best_score  <- pval
      }
    } else if (!criterion %in% c("aic", "bic")) {
      stop(
        paste0(
          "! `criterion` must be one of 'pvalue', 'aic', or 'bic'; got '",
          criterion,
          "'."
        ),
        call. = FALSE
      )
    }
  }

  # -------------------------------------------------------------------------
  # Resolve AIC/BIC winner for the current variable
  # -------------------------------------------------------------------------
  # There is only one active candidate, so this checks whether its improvement
  # statistic is finite and records it as the winner.
  if (criterion %in% c("aic", "bic")) {
    m <- candidate_metrics[[active_type]]

    if (!is.null(m)) {
      if (!ic_col %in% names(m)) {
        stop(
          paste0(
            "Internal error: interaction metrics are missing `",
            ic_col,
            "`."
          ),
          call. = FALSE
        )
      }

      dIC <- m[[ic_col]][1L]

      if (is.finite(dIC)) {
        best_type   <- active_type
        best_fit    <- candidate_full_fits[[active_type]]
        best_metric <- m
        best_score  <- dIC
      }
    }
  }

  # -------------------------------------------------------------------------
  # Verbose per-variable candidate table
  # -------------------------------------------------------------------------
  if (verbose) {
    mfp2_message_capture(
      format_candidate_table(
        verbose_rows,
        criterion = criterion,
        digits = digits,
        best_type = best_type
      )
    )
  }

  # -------------------------------------------------------------------------
  # Immediate final selection when p-values are not being adjusted
  # -------------------------------------------------------------------------
  retained <- NA

  if (!adjusting_pvals && !is.null(best_fit)) {
    retained <- if (criterion == "pvalue") {
      best_score < p_interact
    } else {
      best_score > min_improvement
    }
  } else if (!adjusting_pvals && is.null(best_fit)) {
    retained <- FALSE
    if (verbose) {
      mfp2_message("        >> model fit failed")
    }
  }

  list(
    candidate_metrics = candidate_metrics,
    candidate_models  = candidate_models,
    best_fit          = best_fit,
    best_metric       = best_metric,
    best_type         = best_type,
    best_score        = best_score,
    retained          = retained
  )
}

#' Select Interactions Using Adjusted P-values
#'
#' Adjusts interaction p-values using \code{stats::p.adjust()} and selects
#' interactions whose adjusted p-values are strictly below \code{p_interact}.
#' Updates the per-variable results and selected-model lists.
#'
#' @param var_winners Named list of per-variable fitted results and comparison
#'   metrics.
#' @param all_metrics_list List of comparison-metric data frames.
#' @param all_idx Number of populated entries in \code{all_metrics_list}.
#' @param interaction_vars Character vector naming interaction terms in
#'   evaluation order.
#' @param n_vars Number of interaction variables.
#' @param n_types Number of evaluated forms per variable.
#' @param p_adjust_method Adjustment method passed to \code{stats::p.adjust()}.
#' @param p_interact Significance threshold for adjusted p-values.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{var_winners}}{Per-variable results updated with adjusted
#'     p-values and scores.}
#'   \item{\code{all_candidate_metrics}}{Complete comparison-metric table,
#'     including a \code{p_adjusted} column.}
#'   \item{\code{best_idx}}{Number of selected interactions.}
#'   \item{\code{best_metrics_list}, \code{best_interaction_model},
#'     \code{center_vals_list}}{Comparison metrics, fitted models, and centring
#'     constants for selected interactions.}
#' }
#'
#' @keywords internal
#' @noRd
finalize_adjusted_pvalue_selection <- function(var_winners,
                                               all_metrics_list,
                                               all_idx,
                                               interaction_vars,
                                               n_vars,
                                               n_types,
                                               p_adjust_method,
                                               p_interact) {

  all_candidate_metrics <- bind_metric_rows(all_metrics_list[seq_len(all_idx)])

  if (nrow(all_candidate_metrics) > 0L &&
      "pvalue" %in% names(all_candidate_metrics)) {
    all_candidate_metrics$p_adjusted <- NA_real_

    valid_p <- is.finite(all_candidate_metrics$pvalue)
    n_planned_tests <- length(interaction_vars) * n_types

    if (any(valid_p)) {
      all_candidate_metrics$p_adjusted[valid_p] <- if (p_adjust_method == "none") {
        all_candidate_metrics$pvalue[valid_p]
      } else {
        stats::p.adjust(
          all_candidate_metrics$pvalue[valid_p],
          method = p_adjust_method,
          n      = n_planned_tests
        )
      }
    }

    # Copy adjusted p-values back into the per-variable winners so the final
    # selection and verbose summary use adjusted scores.
    for (vn in interaction_vars) {
      rows_idx <- which(all_candidate_metrics$variable == vn &
                          is.finite(all_candidate_metrics$p_adjusted))
      if (length(rows_idx) == 0L) next

      win_idx <- rows_idx[1L]
      w <- var_winners[[vn]]

      if (!is.null(w$fit)) {
        w$metric <- all_candidate_metrics[win_idx, , drop = FALSE]
        w$score  <- all_candidate_metrics$p_adjusted[win_idx]
        var_winners[[vn]] <- w
      }
    }
  }

  # Rebuild retained-model containers from adjusted p-value decisions.
  best_idx <- 0L
  best_metrics_list <- vector("list", n_vars)
  best_interaction_model <- stats::setNames(vector("list", n_vars), interaction_vars)
  center_vals_list       <- stats::setNames(vector("list", n_vars), interaction_vars)

  for (vn in interaction_vars) {
    w <- var_winners[[vn]]
    if (is.null(w$fit) || is.null(w$metric)) next

    padj_val <- w$metric$p_adjusted[1L]

    if (!is.na(padj_val) && padj_val < p_interact) {
      best_idx <- best_idx + 1L
      best_metrics_list[[best_idx]]      <- w$metric
      best_interaction_model[[vn]]       <- w$fit$test_results$interaction_model
      center_vals_list[[vn]]             <- w$center_vals
    }
  }

  list(
    var_winners            = var_winners,
    all_candidate_metrics  = all_candidate_metrics,
    best_idx               = best_idx,
    best_metrics_list      = best_metrics_list,
    best_interaction_model = best_interaction_model,
    center_vals_list       = center_vals_list
  )
}

#' Print the Final Interaction Summary
#'
#' Prints a summary of interaction comparisons. For AIC or BIC selection,
#' shows criterion improvements for selected interactions. For p-value
#' selection, shows p-values and selection decisions, including adjusted
#' p-values when applicable.
#'
#' @param var_winners Named list of per-variable fitted results and comparison
#'   metrics.
#' @param interaction_vars Character vector of tested terms in display order.
#' @param mode Summary format: \code{"ic"} or \code{"pvalue"}.
#' @param ic_col,ic_label Metric column name and display label for AIC or BIC
#'   improvements. Used when \code{mode = "ic"}.
#' @param min_improvement Selection threshold for AIC or BIC improvements.
#' @param p_interact Significance threshold for p-value selection.
#' @param p_adjust_method P-value adjustment method.
#' @param digits Number of decimal places for displayed statistics.
#' @param show_models Whether to print selected-model summaries after the table.
#' @param best_interaction_model Named list of selected interaction models.
#'
#' @return Invisibly returns \code{NULL} after printing the summary.
#'
#' @keywords internal
#' @noRd
print_interaction_step3_summary <- function(var_winners,
                                            interaction_vars,
                                            mode = c("ic", "pvalue"),
                                            ic_col = NULL,
                                            ic_label = NULL,
                                            min_improvement = NULL,
                                            p_interact = NULL,
                                            p_adjust_method = NULL,
                                            digits = 3,
                                            show_models = FALSE,
                                            best_interaction_model = NULL) {
  mode <- match.arg(mode)

  rule_thick <- strrep("=", 70)
  rule_thin <- strrep("-", 70)
  show_adjusted <- mode == "pvalue" &&
    is.character(p_adjust_method) &&
    length(p_adjust_method) == 1L &&
    !is.na(p_adjust_method) &&
    !identical(tolower(p_adjust_method), "none")

  mfp2_message("")
  mfp2_message(rule_thick)
  mfp2_message("  STEP 3: Interaction Summary")
  mfp2_message(rule_thick)

  if (mode == "ic") {
    mfp2_message(sprintf(
      "  %-12s %-8s %12s %10s",
      "Variable", "Type", ic_label, "selected"
    ))
  } else if (show_adjusted) {
    mfp2_message(sprintf(
      "  %-12s %-8s %12s %12s %10s",
      "Variable", "Type", "p_raw", "p_adjusted", "selected"
    ))
  } else {
    mfp2_message(sprintf(
      "  %-12s %-8s %12s %10s",
      "Variable", "Type", "p-value", "selected"
    ))
  }
  mfp2_message(rule_thin)

  for (vn in interaction_vars) {
    w <- var_winners[[vn]]

    if (is.null(w$fit) || is.null(w$metric)) {
      if (mode == "ic") {
        mfp2_message(sprintf(
          "  %-12s %-8s %12s %10s",
          vn, "---", "NA", "NA"
        ))
      } else if (show_adjusted) {
        mfp2_message(sprintf(
          "  %-12s %-8s %12s %12s %10s",
          vn, "---", "NA", "NA", "NA"
        ))
      } else {
        mfp2_message(sprintf(
          "  %-12s %-8s %12s %10s",
          vn, "---", "NA", "NA"
        ))
      }
      next
    }

    type_label <- format_type_label(w$type)

    if (mode == "ic") {
      ic_val <- if (ic_col %in% names(w$metric)) {
        w$metric[[ic_col]][1L]
      } else {
        NA_real_
      }
      selected <- if (is.finite(ic_val) && ic_val > min_improvement) {
        "Yes"
      } else {
        "No"
      }

      mfp2_message(sprintf(
        "  %-12s %-8s %12s %10s",
        vn,
        type_label,
        format_print_decimal(ic_val, digits),
        selected
      ))
    } else {
      p_raw <- w$metric$pvalue[1L]
      p_adjusted <- if ("p_adjusted" %in% names(w$metric)) {
        w$metric$p_adjusted[1L]
      } else {
        p_raw
      }
      p_for_selection <- if (show_adjusted) p_adjusted else p_raw
      selected <- if (is.finite(p_for_selection) &&
                      p_for_selection < p_interact) {
        "Yes"
      } else {
        "No"
      }

      if (show_adjusted) {
        mfp2_message(sprintf(
          "  %-12s %-8s %12s %12s %10s",
          vn,
          type_label,
          format_print_pvalue(p_raw, digits),
          format_print_pvalue(p_adjusted, digits),
          selected
        ))
      } else {
        mfp2_message(sprintf(
          "  %-12s %-8s %12s %10s",
          vn,
          type_label,
          format_print_pvalue(p_raw, digits),
          selected
        ))
      }
    }
  }

  mfp2_message(rule_thin)
  mfp2_message("")

  if (show_models) {
    for (vn in names(Filter(Negate(is.null), best_interaction_model))) {
      w <- var_winners[[vn]]
      fit_obj <- w$fit$test_results$interaction_model$fit
      if (!is.null(fit_obj)) {
        type_label <- format_type_label(w$type)
        mfp2_message("  Interaction model (", vn, ", ", type_label, "):")
        mfp2_message(strrep("-", 50))
        mfp2_message_capture(print(summary(fit_obj)))
        mfp2_message("")
      }
    }
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# MFPI adjustment-block cache helpers -----------------------------------------
# -----------------------------------------------------------------------------

#' Cache Transformed MFPI Adjustment Terms
#'
#' Transforms each selected adjustment term once and stores its complete
#' design block by term name. Group dummy terms are excluded because
#' interaction models add group main effects separately.
#'
#' The cached blocks allow \code{mfpi_make_xadj_current()} to assemble each
#' interaction's adjustment matrix while excluding the term being tested.
#' Categorical contrasts and multiple transformed columns belonging to one
#' term are kept together.
#'
#' @param x Numeric predictor matrix with shifts and scales already applied.
#' @param selected_vars Character vector naming selected adjustment terms.
#' @param group_dummy_names Names of group dummy columns to exclude, or
#'   \code{NULL}.
#' @param group_dummy_term Character string naming the group dummy term.
#' @param term_to_columns Named list mapping adjustment terms to predictor
#'   columns.
#' @param skip_adjustment Whether to skip transformation and return an empty
#'   cache.
#' @param scale Named numeric vector of scale factors used to restore shifted,
#'   unscaled predictor values before transformation, or \code{NULL}.
#' @param adj_fp_powers Named list of selected adjustment powers.
#' @param center Named logical vector controlling adjustment-term centring.
#' @param acd_vars Named logical vector identifying ACD adjustment terms.
#' @param adj_zero Named logical vector of fitted zero-handling flags.
#' @param adj_catzero Fitted zero-indicator settings for adjustment terms.
#' @param adj_spike Named logical vector of fitted spike-at-zero flags.
#' @param adj_spike_decision Named integer vector of fitted spike-at-zero
#'   decisions.
#' @param adj_acd_parameter Named list of fitted ACD transformation parameters,
#'   or \code{NULL}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{selected_adj_vars}}{Selected adjustment term names after
#'     excluding group dummy terms.}
#'   \item{\code{adj_blocks}}{Named list containing the complete transformed
#'     predictor block for each adjustment term.}
#' }
#'
#' @keywords internal
#' @noRd
mfpi_build_adjustment_block_cache <- function(x,
                                              selected_vars,
                                              group_dummy_names,
                                              group_dummy_term,
                                              term_to_columns,
                                              skip_adjustment,
                                              scale,
                                              adj_fp_powers,
                                              center,
                                              acd_vars,
                                              adj_zero,
                                              adj_catzero,
                                              adj_spike,
                                              adj_spike_decision,
                                              adj_acd_parameter) {
  # Normalize the group dummy names. In some edge cases there may be no dummy
  # matrix or no column names; treating that as an empty exclusion set keeps the
  # downstream setdiff() simple.
  if (is.null(group_dummy_names)) {
    group_dummy_names <- character(0L)
  }

  # Adjustment variables are selected adjustment-model variables excluding the
  # group dummies. The current cont_var is not removed here; it is removed later
  # by mfpi_make_xadj_current(), because that removal changes per loop iteration.
  selected_adj_vars <- setdiff(
    selected_vars,
    c(group_dummy_term, group_dummy_names)
  )
  selected_adj_vars <- unique(
    selected_adj_vars[!is.na(selected_adj_vars) & nzchar(selected_adj_vars)]
  )

  # Always return the same cache structure, even when adjustment is skipped.
  # This keeps evaluate_interactions() simple and avoids special cases inside
  # the interaction_vars loop.
  adj_blocks <- stats::setNames(
    vector("list", length(selected_adj_vars)),
    selected_adj_vars
  )

  if (isTRUE(skip_adjustment) || length(selected_adj_vars) == 0L) {
    return(list(
      selected_adj_vars = selected_adj_vars,
      adj_blocks        = adj_blocks
    ))
  }

  for (v in selected_adj_vars) {
    cols <- term_to_columns[[v]]

    # Work one conceptual term at a time. A categorical term contributes all of
    # its raw contrast columns as one cached block; a continuous singleton uses
    # its selected FP transformation and optional extensions.
    x_one <- x[, cols, drop = FALSE]

    scale_cols <- stats::setNames(rep(1, length(cols)), cols)
    if (!is.null(scale)) {
      shared <- intersect(cols, names(scale))
      scale_cols[shared] <- scale[shared]
    }
    x_one <- sweep(x_one, 2L, scale_cols, "*")

    center_term <- collapse_option_to_terms(
      if (!is.null(center)) center else {
        stats::setNames(rep(FALSE, ncol(x)), colnames(x))
      },
      term_to_columns[v],
      "center"
    )
    expanded_adjustment <- expand_term_metadata_to_columns(
      term_to_columns = term_to_columns[v],
      powers = stats::setNames(list(adj_fp_powers[[v]]), v),
      raw_columns = cols,
      center = center_term,
      acdx = acd_vars,
      zero = adj_zero,
      catzero = adj_catzero,
      spike = adj_spike,
      spike_decision = adj_spike_decision,
      acd_parameter = adj_acd_parameter
    )
    power_cols <- expanded_adjustment$powers
    center_cols <- expanded_adjustment$center
    acd_cols <- expanded_adjustment$acdx
    zero_cols <- expanded_adjustment$zero
    catzero_cols <- expanded_adjustment$catzero
    spike_cols <- expanded_adjustment$spike
    spike_decision_cols <- expanded_adjustment$spike_decision
    acd_parameter_cols <- expanded_adjustment$acd_parameter

    transformed <- transform_matrix(
      x                  = x_one,
      power_list         = power_cols,
      center             = center_cols,
      acdx               = acd_cols,
      zero               = zero_cols,
      catzero            = catzero_cols,
      spike              = spike_cols,
      spike_decision     = spike_decision_cols,
      keep_x_order       = FALSE,
      acd_parameter_list = acd_parameter_cols,
      reset_zero         = FALSE,
      check_binary       = TRUE
    )

    block <- if (is.null(transformed)) NULL else transformed$x_transformed

    # A variable may theoretically produce no fitted columns, for example if it
    # was selected but represented in a degenerate way after transformation.
    # Store NULL explicitly so the assembler can drop it cleanly.
    if (is.null(block) || ncol(block) == 0L) {
      adj_blocks[[v]] <- NULL
      next
    }

    block <- as.matrix(block)
    storage.mode(block) <- "double"

    if (is.null(colnames(block)) ||
        anyNA(colnames(block)) ||
        any(!nzchar(colnames(block)))) {
      stop(
        paste0(
          "! Internal error: transformed adjustment block for `", v,
          "` has missing or empty column names."
        ),
        call. = FALSE
      )
    }

    if (anyNA(block) || any(!is.finite(block))) {
      stop(
        paste0(
          "! Internal error: non-finite values in transformed adjustment ",
          "block for `", v, "`."
        ),
        call. = FALSE
      )
    }

    adj_blocks[[v]] <- block
  }

  # Column names must be unique once blocks are cbind()ed. Duplicates would make
  # coefficient lookup and downstream prediction reconstruction ambiguous.
  adj_colnames <- unlist(
    lapply(adj_blocks, function(z) {
      if (is.null(z)) character(0L) else colnames(z)
    }),
    use.names = FALSE
  )

  if (length(adj_colnames) > 0L && anyDuplicated(adj_colnames)) {
    stop(
      "! Internal error: duplicated transformed adjustment-column names.",
      call. = FALSE
    )
  }

  list(
    selected_adj_vars = selected_adj_vars,
    adj_blocks        = adj_blocks
  )
}

#' Assemble the MFPI Adjustment Matrix for One Interaction
#'
#' Combines cached adjustment blocks, excluding the complete block for the
#' term being tested.
#'
#' @param var_name Character string naming the interaction term.
#' @param adjustment_cache List returned by
#'   \code{mfpi_build_adjustment_block_cache()}.
#' @param skip_adjustment Whether to return \code{xadj} unchanged instead of
#'   assembling cached blocks.
#' @param xadj Existing adjustment matrix or \code{NULL}, used when
#'   \code{skip_adjustment = TRUE}.
#'
#' @return A numeric adjustment matrix, or \code{NULL} if no adjustment terms
#'   remain. When \code{skip_adjustment = TRUE}, returns \code{xadj} unchanged.
#'
#' @keywords internal
#' @noRd
mfpi_make_xadj_current <- function(var_name,
                                   adjustment_cache,
                                   skip_adjustment,
                                   xadj = NULL) {
  if (isTRUE(skip_adjustment)) {
    return(xadj)
  }

  if (!is.character(var_name) || length(var_name) != 1L ||
      is.na(var_name) || !nzchar(var_name)) {
    stop("`var_name` must be a single non-empty character string.",
         call. = FALSE)
  }

  if (!is.list(adjustment_cache) ||
      !all(c("selected_adj_vars", "adj_blocks") %in% names(adjustment_cache))) {
    stop(
      "`adjustment_cache` must be returned by `mfpi_build_adjustment_block_cache()`.",
      call. = FALSE
    )
  }

  selected_adj_vars <- adjustment_cache$selected_adj_vars
  adj_blocks        <- adjustment_cache$adj_blocks

  # Remove the current tested variable from its own adjustment set. This drops
  # all derived columns for that source variable, not just columns matching a
  # transformed-name pattern.
  adj_vars_current <- setdiff(selected_adj_vars, var_name)

  if (length(adj_vars_current) == 0L) {
    return(NULL)
  }

  blocks <- adj_blocks[adj_vars_current]
  blocks <- blocks[!vapply(blocks, is.null, logical(1L))]

  if (length(blocks) == 0L) {
    return(NULL)
  }

  out <- do.call(cbind, blocks)

  if (is.null(out) || ncol(out) == 0L) {
    return(NULL)
  }

  out <- as.matrix(out)
  storage.mode(out) <- "double"

  if (is.null(colnames(out)) ||
      anyNA(colnames(out)) ||
      any(!nzchar(colnames(out)))) {
    stop(
      "! Internal error: current adjustment matrix has missing column names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(colnames(out))) {
    stop(
      "! Internal error: current adjustment matrix has duplicated columns.",
      call. = FALSE
    )
  }

  if (anyNA(out) || any(!is.finite(out))) {
    stop(
      paste0(
        "! Internal error: current adjustment matrix for `", var_name,
        "` contains non-finite values."
      ),
      call. = FALSE
    )
  }

  out
}


#' Add Cause-and-Fix Guidance to a Rank-Deficient MFPI Error
#'
#' Evaluates `expr`; if it signals `mfpi_rank_deficient_error`, re-signals the
#' same error with guidance appended when the failure has the typical
#' signature of an automatic shift: the variable was shifted automatically,
#' its smallest shifted value is very close to 0 relative to its range (below
#' 1%), and the selected powers include a negative or log (0) power. Otherwise
#' the error is re-signalled unchanged. No shift value is computed or
#' suggested; the guidance only names the cause and the available options.
#'
#' @param expr Expression to evaluate (the `flex_fit()` call).
#' @param var_name Interaction variable name.
#' @param interaction_spec Interaction specification for `var_name`.
#' @param interaction_form Resolved interaction form (`"fp1"`, `"fp2"`, ...).
#' @param x Internal (shifted and scaled) training matrix.
#' @param scale,shift Named scale and shift vectors used at fit time.
#' @param shift_automatic Named logical vector: `TRUE` where the shift was
#'   estimated automatically rather than supplied by the user.
#' @param zero_var Logical; whether `var_name` uses zero handling.
#' @keywords internal
#' @noRd
mfpi_with_rank_hint <- function(expr, var_name, interaction_spec,
                                interaction_form, x, scale, shift,
                                shift_automatic, zero_var = FALSE) {
  tryCatch(
    expr,
    mfpi_rank_deficient_error = function(e) {
      hint <- mfpi_rank_deficient_hint(
        powers = e$selected_powers,
        var_name = var_name,
        interaction_spec = interaction_spec,
        interaction_form = interaction_form,
        x = x,
        scale = scale,
        shift = shift,
        shift_automatic = shift_automatic,
        zero_var = zero_var
      )
      if (!is.null(hint)) e$message <- paste0(e$message, "\n", hint)
      stop(e)
    }
  )
}


#' Build the Automatic-Shift Guidance for a Rank-Deficient MFPI Design
#'
#' @inheritParams mfpi_with_rank_hint
#' @param powers Selected FP powers attached to the error, or `NULL`.
#' @return A character scalar, or `NULL` when the failure does not match the
#'   automatic-shift signature.
#' @keywords internal
#' @noRd
mfpi_rank_deficient_hint <- function(powers, var_name, interaction_spec,
                                     interaction_form, x, scale, shift,
                                     shift_automatic, zero_var = FALSE) {
  if (is.null(interaction_spec) || !identical(interaction_spec$kind, "continuous")) {
    return(NULL)
  }
  if (is.null(powers) || !any(powers <= 0)) return(NULL)
  if (is.null(shift_automatic) || !isTRUE(shift_automatic[var_name])) return(NULL)

  shift_v <- if (!is.null(shift) && var_name %in% names(shift)) unname(shift[[var_name]]) else 0
  scale_v <- if (!is.null(scale) && var_name %in% names(scale)) unname(scale[[var_name]]) else 1
  if (!is.finite(shift_v) || shift_v <= 0 || !var_name %in% colnames(x)) return(NULL)

  x_plus <- as.numeric(x[, var_name]) * scale_v
  if (isTRUE(zero_var)) x_plus <- x_plus[x_plus != 0]
  x_plus <- x_plus[is.finite(x_plus)]
  if (length(x_plus) < 2L) return(NULL)
  value_range <- diff(range(x_plus))
  smallest <- min(x_plus)
  if (!(value_range > 0) || smallest >= 0.01 * value_range) return(NULL)

  fmt <- function(v) format(signif(v, 3), scientific = abs(v) < 1e-3 && v != 0)
  simpler <- if (identical(interaction_form, "fp2")) {
    paste0("`interaction_forms = c(", var_name, " = \"fp1\")` or `\"linear\"`")
  } else {
    paste0("`interaction_forms = c(", var_name, " = \"linear\")`")
  }
  paste0(
    "i Likely cause: `", var_name, "` was shifted automatically (shift = ",
    fmt(shift_v), "), which places its smallest value at ", fmt(smallest),
    ", very close to 0 relative to its range (", fmt(value_range), "). ",
    "The selected power(s) (", paste(powers, collapse = ", "), ") include a ",
    "negative or log power, which gives the observations near 0 extreme values.\n",
    "i Possible fixes:\n",
    "  * supply your own shift for `", var_name, "` so that its smallest value ",
    "is not close to 0: `fp(", var_name, ", shift = ...)` in a formula, or ",
    "`shift = c(", var_name, " = ...)` with the matrix interface;\n",
    "  * use a simpler interaction form, e.g. ", simpler, ";\n",
    "  * try another `flex` level."
  )
}
