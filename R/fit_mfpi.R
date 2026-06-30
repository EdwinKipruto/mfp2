# Internal implementation functions for mfpi()
#
# None of these functions are exported. They are called exclusively by
# mfpi.default() after all argument checking and pre-processing has been
# completed. Parameters therefore arrive already validated, expanded to full
# length, shifted, and scaled.

# -----------------------------------------------------------------------------
# fit_mfpi() ------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the Full MFPI Model
#'
#' Core workhorse called by \code{mfpi()} after argument validation and
#' pre-processing. Runs the three-step MFPI algorithm: data pre-processing,
#' adjustment-model selection via MFP, and univariable interaction evaluation.
#'
#' This function performs **no argument checking**. All inputs are assumed to
#' have been validated, expanded to full-length named vectors, shifted, and
#' scaled by the calling \code{mfpi()} function.
#'
#' @param x Numeric matrix (\eqn{n \times p}) of predictors including
#'   \code{group_var}. Shift and scale have already been applied; no intercept
#'   column.
#' @param y Response vector or \code{survival::Surv()} object.
#' @param family Character string; one of \code{"gaussian"}, \code{"binomial"},
#'   \code{"poisson"}, or \code{"cox"}.
#' @param family_string Same as \code{family} but always a plain character
#'   string. Required separately by \code{fit_mfp()} for internal
#'   branching.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param center Named logical vector of length \eqn{p}. Whether to center
#'   each predictor before fitting.
#' @param criterion Character string; \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}. Governs two distinct selection steps:
#'   \enumerate{
#'     \item \strong{Adjustment-model selection} (Step 1): controls variable
#'       elimination and FP degree selection in \code{fit_mfp()}.
#'     \item \strong{Interaction test} (Step 2): tests the pre-specified
#'       functional form from \code{cont_var_forms} against its main-effects
#'       counterpart in \code{evaluate_interactions()}.
#'   }
#'   In the internal \code{fit_mfp()} calls within \code{flex1()} and
#'   \code{flex4()} - used only to estimate FP powers for \code{cont_var} -
#'   the user criterion is passed alongside \code{force_max_fp = TRUE}, which
#'   prevents AIC/BIC from simplifying the functional form below the requested
#'   degree. The best power combination within that degree is still selected by
#'   the criterion (equivalently, by deviance minimization at fixed df).
#' @param select Named numeric vector of length \eqn{p}. Nominal significance
#'   levels for backward elimination of each predictor.
#' @param alpha Named numeric vector of length \eqn{p}. Significance levels
#'   for FP degree selection. Under \code{criterion = "pvalue"}, \code{alpha = 1}
#'   guarantees the most complex FP degree is always accepted. Ignored under
#'   \code{criterion = "aic"} or \code{"bic"}.
#' @param keep Character vector of variable names forced into the
#'   adjustment model regardless of selection.
#' @param force_max_fp Named logical vector of length \eqn{p}. For each
#'   variable, if \code{TRUE}, forces \code{select_ic()} to select the most
#'   complex functional form at the degree specified by \code{df}, bypassing
#'   AIC/BIC competition against simpler forms. Expanded from a scalar or
#'   validated as a named vector of length \eqn{p} by \code{mfpi.default()}
#'   before being passed here. Passed directly to \code{fit_adjustment_model()};
#'   not forwarded to the flex functions, which construct their own internal
#'   \code{force_max_fp} vectors from \code{vnames}.
#' @param df Named integer vector of length \eqn{p}. Degrees of freedom per
#'   predictor (1 = linear, 2m = FP of degree m), after cardinality-based
#'   overrides by \code{assign_df()}.
#' @param xorder Character string; order of covariate entry into MFP
#'   backfitting - \code{"ascending"}, \code{"descending"}, or
#'   \code{"original"}.
#' @param fp_powers Named list of candidate FP power sets, one per predictor.
#' @param ties Character string; tie-handling method for Cox models.
#' @param strata Integer vector of stratum memberships for stratified Cox
#'   models, or \code{NULL}.
#' @param nocenter Numeric vector passed to \code{survival::coxph()} to
#'   suppress centring for specific predictors.
#' @param acd_vars Named logical vector of length \eqn{p}. Whether each
#'   predictor undergoes the ACD transformation. Must be \code{FALSE} for all
#'   variables in \code{cont_vars}.
#' @param zero_vars Named logical vector of length \eqn{p}. Whether each
#'   predictor treats non-positive values as zero before transformation.
#' @param catzero_vars Named logical vector of length \eqn{p}. Whether each
#'   predictor is semi-continuous and requires a binary indicator alongside
#'   its FP transformation.
#' @param spike_vars Named logical vector of length \eqn{p}. Whether each
#'   predictor is subject to spike-at-zero (SAZ) modelling.
#' @param min_prop Numeric in \eqn{[0,1]}. Minimum proportion of zeros
#'   required for SAZ modelling.
#' @param max_prop Numeric in \eqn{[0,1]}. Maximum proportion of zeros for
#'   SAZ modelling.
#' @param use_ftest Logical. If \code{TRUE} and \code{family = "gaussian"},
#'   use an F-test (via \code{calculate_f_test()}) rather than a
#'   chi-square likelihood-ratio test for the interaction p-value. Also
#'   passed to \code{fit_mfp()} for adjustment-variable selection.
#'   Ignored for non-Gaussian families.
#' @param control List of fitting control parameters.
#' @param group_var Character string. Name of the grouping variable in \code{x}.
#' @param include_group_var Logical. Whether \code{group_var} is included as a
#'   covariate in the adjustment model.
#' @param flex Character string; \code{"flex1"}, \code{"flex2"},
#'   \code{"flex3"}, or \code{"flex4"}.
#' @param cont_vars Character vector of continuous variables to test for
#'   interaction with \code{group_var}.
#' @param cont_var_forms Named character vector with names matching
#'   \code{cont_vars} and values in \code{c("linear", "fp1", "fp2")}.
#'   Specifies the functional form to use for each variable. Constructed and
#'   validated by \code{mfpi.default()} before being passed here.
#' @param p_interact Numeric. Significance threshold for
#'   \code{criterion = "pvalue"}.
#' @param min_improvement Numeric. Minimum criterion improvement required to
#'   retain an interaction.
#' @param show_models Logical. If \code{TRUE} and \code{verbose = TRUE},
#'   prints the full regression coefficient table for each selected
#'   interaction model. Has no effect when \code{verbose = FALSE}.
#'   Nothing is printed when no variable is selected.
#' @param verbose Logical. Whether to print progress messages.
#' @param quiet Logical. Suppresses per-variable messages when no significant
#'   interaction is found.
#' @param center_type Character string; \code{"grand"} (default) or
#'   \code{"group"}. Passed to \code{create_z_variables()} to control the
#'   centering strategy when \code{center = TRUE}. The constants are stored
#'   in \code{center_vals_list} on the returned list.
#' @param scale Named numeric vector of per-variable scale factors (one per
#'   column of \code{x}), as computed and applied in \code{mfpi.default()}.
#'   Passed to \code{flex_fit()} as \code{scale_var} so that
#'   \code{create_z_variables()} and \code{transform_z_variables()} can
#'   backscale \code{cont_var} (multiply by the scale factor) before FP
#'   transformation. This ensures interaction model coefficients are on the
#'   \eqn{\phi(x + \text{shift})} scale, matching the adjustment model and
#'   standalone \pkg{mfp2}. \code{NULL} or 1 means no backscaling.
#' @param shift Named numeric vector of per-variable shift factors (one per
#'   column of \code{x}), as computed and applied in \code{mfpi.default()}.
#'   Passed through to fitted-function generation so displayed x-values can be
#'   returned on the original raw scale.
#' @param digits Positive integer. Significant digits for printed output.
#' @param p_adjust_method Character string. Method for adjusting p-values
#'   across all variables in \code{cont_vars} when \code{criterion = "pvalue"}.
#'   Passed to \code{\link[stats]{p.adjust}}. Default \code{"none"}.
#' @return A list with top-level fields (accessible as \code{fit$field}) and
#'   two retained sub-objects. Top-level fields:
#' \describe{
#'   \item{\code{best_model_metrics}}{Data frame of metrics for the best
#'     interaction model per significant variable.}
#'   \item{\code{all_model_metrics}}{Tibble of metrics for the interaction
#'     model evaluated for each variable in \code{cont_vars}, using the form
#'     specified by \code{cont_var_forms}.}
#'   \item{\code{best_interaction_model}}{Named list of fitted interaction model
#'     objects for the **winning** model per significant variable.}
#'   \item{\code{all_interaction_models}}{Named list. Names are \code{cont_vars};
#'     each element is the fitted model object for that variable using the form
#'     specified in \code{cont_var_forms}, regardless of whether it was selected
#'     or significant.}
#'   \item{\code{center_vals_list}}{Named list (one element per
#'     \code{cont_var}) of centering constants used when fitting the interaction
#'     model. \code{NULL} when \code{center = FALSE}.}
#'   \item{\code{adjust_terms}}{The \code{fp_terms} table from the adjustment
#'     model.}
#'   \item{\code{group_var}, \code{flex}, \code{family}, \code{nobs},
#'     \code{show_models}, \code{group_levels_new},
#'     \code{group_levels_original}}{Metadata fields.}
#'   \item{\code{adjustment_model}}{Full adjustment model object.}
#'   \item{\code{p_adjust_method}}{The multiplicity adjustment method used.}
#'   \item{\code{var_winners}}{Named list of per-variable best candidates
#'     (regardless of selection). See \code{evaluate_interactions()}.}
#'   \item{\code{univariable_interactions}}{Full list from
#'     \code{evaluate_interactions()}.}
#' }
#'
#' @section Algorithm:
#' \strong{Pre-processing.} Levels of \code{group_var} are remapped to
#' consecutive integers. Group dummies are computed once and stored in
#' \code{cat_info$dummies}. When \code{include_group_var = TRUE}, dummies are
#' appended to \code{x}. All parameter vectors are synchronized silently.
#'
#' \strong{Step 1 - Adjustment model.} \code{fit_mfp()} selects
#' adjustment variables and FP transformations using \code{criterion}.
#' Selected variables, their FP powers, and spike decisions are fixed for
#' the remainder of the algorithm.
#'
#' \strong{Step 2 - Univariable interactions.} For each variable in
#' \code{cont_vars}, the interaction model is fitted using the form
#' pre-specified in \code{cont_var_forms} and compared to its main-effects
#' model via \code{evaluate_interactions()}. No data-driven selection among
#' functional forms is performed.
#'
#' @seealso \code{mfpi()}, \code{preprocess_data()},
#'   \code{fit_adjustment_model()}, \code{evaluate_interactions()}
#'
#' @keywords internal
#' @noRd
fit_mfpi <- function(x, y, family, family_string, weights, offset, cycles,
                     center, criterion, select, alpha, keep,
                     force_max_fp,
                     df, xorder,
                     fp_powers, ties, strata, nocenter, acd_vars, zero_vars,
                     catzero_vars, spike_vars, min_prop, max_prop, use_ftest,
                     control, group_var, include_group_var, flex, cont_vars,
                     cont_var_forms,
                     p_interact, min_improvement, show_models, verbose,
                     digits, scale, shift = NULL, quiet = FALSE, has_offset,
                     center_type = c("grand", "group"),
                     p_adjust_method = "none",
                     group_input_levels = NULL,
                     group_levels_original = NULL) {
  
  center_type <- match.arg(center_type)
  
  # Resolve GLM family objects once for all repeated internal model fits.
  # Public mfp2.default() already does this, but keeping it here makes direct
  # internal calls to fit_mfp() avoid repeated stats::gaussian()/binomial()/
  # poisson() construction as well. Cox remains the character string "cox".
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
    group_levels_original = group_levels_original
  )
  
  # ---------------------------------------------------------------------------
  # Step 1: Fit adjustment model via MFP
  # ---------------------------------------------------------------------------
  if (verbose) {
    rule_thick <- strrep("=", 70)
    rule_thin  <- strrep("-", 70)
    cat("\n", rule_thick, "\n", sep = "")
    cat("  STEP 1: Adjustment Model (MFP selection)\n")
    cat(rule_thick, "\n", sep = "")
  }
  
  adjustment_model <- fit_adjustment_model(
    x              = processed_data$x,
    y              = y,
    weights        = weights,
    offset         = offset,
    cycles         = cycles,
    family         = family_fit,
    family_string  = family_string,
    criterion      = criterion,
    updated_params = processed_data$updated_params,
    xorder         = xorder,
    ties           = ties,
    strata         = strata,
    nocenter       = nocenter,
    min_prop       = min_prop,
    max_prop       = max_prop,
    use_ftest      = use_ftest,
    control        = control,
    force_max_fp   = processed_data$updated_params$force_max_fp,
    scale          = scale,
    has_offset     = has_offset,
    verbose        = FALSE
  )
  
  # Identify variables retained by the adjustment model
  selected_vars <- get_selected_variables(adjustment_model)
  if (include_group_var) {
    selected_vars <- setdiff(selected_vars, processed_data$dummy_names)
  }
  
  if (verbose) {
    print(adjustment_model$fp_terms)
    
    # Print the names of selected adjustment variables for quick reference,
    # excluding the group dummies (which are structural, not adjustment).
    group_dummy_names <- colnames(processed_data$cat_info$dummies)
    adj_selected     <- setdiff(selected_vars, group_dummy_names)
    cat("\n")
    if (length(adj_selected) > 0L) {
      cat(sprintf(
        "  Selected adjustment variables (%d): %s\n",
        length(adj_selected),
        paste(adj_selected, collapse = ", ")
      ))
    } else {
      cat("  Selected adjustment variables: (none)\n")
    }
  }
  
  # `get_fp_powers()` is public and intentionally requires a non-empty
  # variable vector. Internally, however, an empty selected adjustment model
  # is valid: interaction models are then fitted without adjustment covariates.
  adj_fp_powers <- if (length(selected_vars) > 0L) {
    get_fp_powers(
      selected_vars,
      adjustment_model$fp_terms
    )
  } else {
    stats::setNames(vector("list", 0L), character(0L))
  }
  
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
  extract_adjustment_flag <- function(object, field, selected_vars, default = FALSE) {
    flag <- object[[field]]
    
    if (is.null(flag)) {
      return(stats::setNames(rep(default, length(selected_vars)), selected_vars))
    }
    
    if (is.null(names(flag))) {
      stop(
        paste0(
          "! Internal error: `adjustment_model$", field,
          "` must be named."
        ),
        call. = FALSE
      )
    }
    
    missing_flag <- setdiff(selected_vars, names(flag))
    if (length(missing_flag) > 0L) {
      stop(
        paste0(
          "! Internal error: `adjustment_model$", field,
          "` is missing selected variable(s): ",
          paste(missing_flag, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
    
    out <- flag[selected_vars]
    names(out) <- selected_vars
    out
  }
  
  adj_zero <- extract_adjustment_flag(
    adjustment_model,
    "zero",
    selected_vars,
    default = FALSE
  )
  
  adj_catzero <- extract_adjustment_flag(
    adjustment_model,
    "catzero",
    selected_vars,
    default = FALSE
  )
  
  adj_spike <- extract_adjustment_flag(
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
    cat("\n", rule_thick, "\n", sep = "")
    cat(sprintf(
      "  STEP 2: Evaluating Interactions  (flex = %s, criterion = '%s')\n",
      flex, criterion
    ))
    cat(rule_thick, "\n", sep = "")
    cat(sprintf(
      "  Testing %d continuous variable(s) against group '%s'\n",
      length(cont_vars), group_var
    ))
    if (criterion == "pvalue") {
      cat(sprintf("  Threshold: p_interact = %g\n", p_interact))
    } else {
      cat(sprintf("  Threshold: min_improvement = %g\n", min_improvement))
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
    cont_vars          = cont_vars,
    cont_var_forms     = cont_var_forms,
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
    fp_powers          = fp_powers,
    cycles             = cycles,
    criterion          = criterion,
    group_var          = group_var,
    show_models        = show_models,
    p_interact         = p_interact,
    min_improvement    = min_improvement,
    min_prop           = min_prop,
    max_prop           = max_prop,
    center_type        = center_type,
    scale              = scale,
    shift              = shift,
    adj_acd_parameter = adjustment_model$acd_parameter,
    quiet              = quiet,
    p_adjust_method    = p_adjust_method,
    digits             = digits,
    verbose            = verbose
  )
  
  n_significant <- nrow(univ_results$best_model_metrics)
  if (n_significant == 0L && !quiet) {
    message("! No significant interactions were found.", call. = FALSE)
  }
  
  if (verbose) {
    rule_thick <- strrep("=", 70)
    cat("\n", rule_thick, "\n", sep = "")
    cat(sprintf(
      "  DONE: %d of %d variable(s) show significant interaction\n",
      n_significant, length(cont_vars)
    ))
    cat(rule_thick, "\n\n", sep = "")
  }
  
  list(
    best_model_metrics       = univ_results$best_model_metrics,
    all_model_metrics        = univ_results$all_model_metrics,
    best_interaction_model   = univ_results$best_interaction_model,
    all_interaction_models   = univ_results$all_interaction_models,
    center_vals_list         = univ_results$center_vals_list,
    adjust_terms             = adjustment_model$fp_terms,
    group_var                = univ_results$group_var,
    show_models              = univ_results$show_models,
    flex                     = univ_results$flex,
    criterion                = univ_results$criterion,
    group_levels_new         = univ_results$group_levels_new,
    group_levels_original    = univ_results$group_levels_original,
    group_level_map          = univ_results$group_level_map,
    family                   = univ_results$family,
    nobs                     = univ_results$nobs,
    p_adjust_method          = univ_results$p_adjust_method,
    scale                    = scale,
    shift                    = shift,
    var_winners              = univ_results$var_winners,
    adjustment_model         = adjustment_model,
    univariable_interactions = univ_results,
    digits                   = digits
  )
}


# -----------------------------------------------------------------------------
# preprocess_data() -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Pre-process Predictor Data for the MFPI Adjustment Model
#'
#' Removes \code{group_var} from the predictor matrix and, when
#' \code{include_group_var = TRUE}, appends dummy variables for the non-reference
#' levels of \code{group_var}. All modelling-parameter vectors and lists
#' (\code{select}, \code{alpha}, \code{df}, \code{center}, \code{fp_powers},
#' etc.) are synchronised with the modified column set so that downstream
#' functions receive consistently aligned inputs.
#'
#' The grouping variable is internally remapped to consecutive integer levels
#' \code{0, 1, ..., K - 1}. This is important because later MFPI helper
#' functions build and parse column names using the group value as a suffix.
#' The original grouping values and the internal mapping are retained in
#' \code{cat_info} for reporting and interpretation.
#'
#' This function is called at the start of \code{fit_mfpi()} and its output is
#' threaded through the adjustment-model and interaction-model stages. It assumes
#' that public-facing validation and scalar expansion have already been handled
#' by \code{mfpi.default()}.
#'
#' @param x Numeric matrix of predictors including \code{group_var}.
#' @param group_var Character string. Name of the categorical grouping variable.
#' @param include_group_var Logical. If \code{TRUE}, dummy variables derived
#'   from \code{group_var} are appended to the filtered predictor matrix so that
#'   group membership can be adjusted for in the MFP step.
#' @param select Named numeric vector of selection levels, one per column of
#'   \code{x}.
#' @param alpha Named numeric vector of FP-degree significance levels, one per
#'   column of \code{x}.
#' @param df Named numeric vector of degrees of freedom, one per column of
#'   \code{x}.
#' @param center Named logical vector. Whether to centre each predictor.
#' @param acd_vars Named logical vector. Whether each predictor undergoes the
#'   ACD transformation.
#' @param fp_powers Named list of candidate FP power sets, one per predictor.
#' @param zero_vars Named logical vector. Whether each predictor should treat
#'   non-positive values as zero.
#' @param catzero_vars Named logical vector. Whether each predictor should add a
#'   structural-zero indicator.
#' @param spike_vars Named logical vector. Whether each predictor should be
#'   assessed for spike-at-zero handling.
#' @param keep Character vector of variable names always retained in the
#'   adjustment model.
#' @param force_max_fp Named logical vector of length \eqn{p}. Sliced to
#'   adjustment columns and stored in \code{updated_params$force_max_fp}. When
#'   \code{include_group_var = TRUE}, dummy columns receive \code{FALSE}
#'   because they have \code{df = 1} and no FP form to force.
#'
#' @return A list with five components:
#' \describe{
#'   \item{\code{x}}{Numeric matrix with \code{group_var} removed and,
#'     optionally, group dummy columns appended.}
#'   \item{\code{original_x}}{Internal predictor matrix with \code{group_var}
#'     remapped to consecutive integer levels \code{0, 1, ..., K - 1}.}
#'   \item{\code{cat_info}}{List of metadata for \code{group_var}, including
#'     remapped values, original values, original levels, internal levels, the
#'     level map, the group variable name, and the dummy matrix.}
#'   \item{\code{dummy_names}}{Character vector of appended dummy column names,
#'     or \code{NULL} if \code{include_group_var = FALSE}.}
#'   \item{\code{updated_params}}{Named list of parameter vectors/lists aligned
#'     with the columns of the returned \code{x}.}
#' }
#'
#' @keywords internal
#' @noRd
preprocess_data <- function(x, group_var, include_group_var,
                            select, alpha, df, center, acd_vars,
                            fp_powers, zero_vars, catzero_vars, spike_vars,
                            keep, force_max_fp, group_input_levels = NULL,
                            group_levels_original = NULL) {
  
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
  
  # Developer note:
  # group_var is removed from the adjustment matrix. If the user listed group_var
  # in keep, remove it here; otherwise fit_mfp() would receive a keep variable
  # that is no longer a column of x.
  keep_filtered <- intersect(
    if (is.null(keep)) character(0L) else keep,
    predictor_names
  )
  
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
    
    updated_params$keep <- c(keep_filtered, dummy_names)
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
    updated_params = updated_params
  )
}

# -----------------------------------------------------------------------------
# fit_adjustment_model() ------------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the MFP Adjustment Model
#'
#' A thin wrapper around `fit_mfp()` that selects adjustment variables
#' and their FP transformations. Scale and shift are set to 1 and 0
#' respectively because \code{mfpi()} has already applied them to `x`.
#'
#' @param x Numeric matrix of predictor variables passed to
#'   \code{fit_mfp()} for adjustment-variable selection. It is the
#'   filtered matrix produced by \code{preprocess_data()}, which means:
#'   \itemize{
#'     \item \code{group_var} has been \strong{removed} - it is not adjusted
#'       for here because it is the variable whose interaction with
#'       \code{cont_vars} is being tested.
#'     \item If \code{include_group_var = TRUE} in \code{mfpi()},
#'       \code{preprocess_data()} has already re-added \code{group_var} as a
#'       set of binary dummy columns (one per non-reference level), with
#'       \code{select = 1} to force them into the adjustment model.
#'     \item All columns have already been shifted and scaled by
#'       \code{mfpi.default()} - shift and scale are therefore passed as 0
#'       and 1 respectively to \code{fit_mfp()} inside this function.
#'   }
#' @param y Response vector or Surv object.
#' @param weights Numeric vector of observation weights.
#' @param offset Numeric vector of linear-predictor offsets.
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param family GLM family object (e.g. \code{gaussian()}) or the character
#'   string \code{"cox"}. Passed directly to \code{fit_mfp()} as its
#'   \code{family} argument.
#' @param family_string Character string of the family name (e.g.
#'   \code{"gaussian"}, \code{"cox"}). Required separately by
#'   \code{fit_mfp()} for internal branching.
#' @param criterion Character string; \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}. Governs full MFP selection: variable elimination, FP
#'   degree selection, and functional form for adjustment variables. Unlike
#'   the internal \code{fit_mfp()} calls in \code{flex1()} and
#'   \code{flex4()}, \code{force_max_fp} is not passed here - it defaults to
#'   \code{FALSE} so the user criterion fully controls degree selection, which
#'   is the correct behaviour for adjustment-model fitting.
#' @param updated_params Named list of per-variable modelling parameters as
#'   returned by \code{preprocess_data()}. Must contain \code{df},
#'   \code{center}, \code{select}, \code{alpha}, \code{keep},
#'   \code{fp_powers}, \code{acd_vars}, \code{zero_vars},
#'   \code{catzero_vars}, and \code{spike_vars}.
#' @param xorder Character string; entry order for MFP backfitting.
#' @param ties Character string; tie-handling for Cox models.
#' @param strata Integer stratum vector for stratified Cox models, or
#'   \code{NULL}.
#' @param nocenter Numeric vector passed to \code{survival::coxph()}.
#' @param min_prop Numeric. Minimum proportion of zeros for SAZ modelling.
#' @param max_prop Numeric. Maximum proportion of zeros for SAZ modelling.
#' @param use_ftest Logical. If \code{TRUE} and \code{family = "gaussian"}, use an F-test rather than a chi-square test. Applied to both adjustment-variable selection and the interaction test.
#' @param control Fitting control list.
#' @param verbose Logical. Passed directly to \code{fit_mfp()}.
#'
#' @return The fitted model object returned by \code{fit_mfp()}, which
#'   includes \code{fp_terms} (a data frame of selected variables and their FP
#'   powers) and the underlying model fit.
#'
#' @details
#' The real per-variable \code{scale} factors are now passed to
#' \code{fit_mfp()} so that it can backscale \code{x} before the final model
#' fit (step 4 of \code{fit_mfp()}), giving adjustment model coefficients on
#' the \eqn{\phi(x + \text{shift})} scale, matching standalone \code{mfp2()}.
#' \code{shift} is set to 0 because shifting was already applied upstream in
#' \code{mfpi.default()}.
#'
#' \code{force_max_fp} is a named logical vector already subsetted by
#' \code{fit_mfpi()} to match the columns of \code{x} (i.e. with
#' \code{group_var} removed). When the user supplies \code{force_max_fp = FALSE}
#' (the default), all elements are \code{FALSE} and \code{select_ic()} runs its
#' normal degree competition - correct behaviour for full adjustment-model
#' selection. \code{force_max_fp = TRUE} for specific variables forces the most
#' complex functional form for those variables in the adjustment model.
#' @keywords internal
#' @noRd
fit_adjustment_model <- function(x, y, weights, offset, cycles, family,
                                 family_string, criterion, updated_params,
                                 xorder, ties, strata, nocenter, min_prop,
                                 max_prop, use_ftest, control,
                                 force_max_fp, scale,has_offset,
                                 verbose = FALSE) {
  n_vars    <- ncol(x)
  x_names   <- colnames(x)
  
  # Build scale_vec: real scale factor for predictor columns, 1 for group
  # dummy columns. Dummy column names (e.g. "rx1", "svi1") are not present
  # in names(scale) since scale is keyed by the original variable names from
  # mfpi.default. The intersect() therefore naturally assigns scale = 1 to
  # dummies, which is correct -- dummies are binary and should not be scaled.
  scale_vec <- setNames(rep(1, n_vars), x_names)
  if (!is.null(scale)) {
    shared <- intersect(x_names, names(scale))
    scale_vec[shared] <- scale[shared]
  }
  
  fit_mfp(
    x             = x,
    y             = y,
    weights       = weights,
    offset        = offset,
    cycles        = cycles,
    method        = ties,
    strata        = strata,
    nocenter      = nocenter,
    scale         = scale_vec,
    shift         = rep(0, n_vars),
    ftest         = use_ftest,
    control       = control,
    family        = family,
    family_string = family_string,
    criterion     = criterion,
    xorder        = xorder,
    df            = updated_params$df,
    center        = updated_params$center,
    select        = updated_params$select,
    alpha         = updated_params$alpha,
    keep          = updated_params$keep,
    powers        = updated_params$fp_powers,
    acdx          = updated_params$acd_vars,
    zero          = updated_params$zero_vars,
    catzero       = updated_params$catzero_vars,
    spike         = updated_params$spike_vars,
    min_prop      = min_prop,
    max_prop      = max_prop,
    force_max_fp  = force_max_fp,
    has_offset    = has_offset,
    verbose       = verbose
  )
}

# -----------------------------------------------------------------------------
# format_type_label() ---------------------------------------------------------
# -----------------------------------------------------------------------------

#' Format a functional form type string for display
#'
#' Converts the internal lowercase type string to a display label:
#' \code{"linear"} -> \code{"Linear"}, \code{"fp1"} -> \code{"FP1"},
#' \code{"fp2"} -> \code{"FP2"}. Vectorised over character vectors.
#'
#' @param type Character vector of type strings.
#' @return Character vector of display labels.
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

#' Format and Print a Table of Interaction Candidates for One Variable
#'
#' Assembles a compact summary table from the per-candidate metric rows
#' collected during the interaction evaluation loop and prints it to the
#' console. The displayed columns adapt to the active \code{criterion}:
#'
#' \itemize{
#'   \item \code{criterion = "pvalue"}: shows \code{deviance_int}
#'     (deviance of the interaction model, \eqn{-2\ell_{\mathrm{int}}}),
#'     \code{deviance_diff} (likelihood-ratio statistic
#'     \eqn{-2\ell_{\mathrm{main}} + 2\ell_{\mathrm{int}}}),
#'     \code{df} (interaction degrees of freedom), and \code{pvalue}.
#'   \item \code{criterion = "aic"}: shows \code{AIC_main},
#'     \code{AIC_interaction}, and \code{dAIC}.
#'   \item \code{criterion = "bic"}: shows \code{BIC_main},
#'     \code{BIC_interaction}, and \code{dBIC}.
#' }
#'
#' The table also includes two power columns to help the user distinguish
#' the FP powers from the pooled main-effects model (\code{main_power}) and
#' the FP powers used in the interaction model (\code{int_power}). For
#' flex1-3, interaction powers are shown as a single tuple, e.g.
#' \code{(1, 2)}. For flex4, per-group powers are shown separately, e.g.
#' \code{(1, 2), (0.5)}.
#'
#' @param rows Named list of one-row data frames, one per fitted type. With
#'   pre-specified forms this will always contain exactly one entry. Each row
#'   must contain fields \code{pvalue},
#'   \code{AIC_main}, \code{AIC_interaction}, \code{AIC_main_minus_int},
#'   \code{BIC_main}, \code{BIC_interaction}, \code{BIC_main_minus_int},
#'   \code{fp_powers_main}, and \code{fp_powers_int}.
#' @param criterion Character string: \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}.
#' @param digits Integer; number of significant digits for p-values.
#'
#' @return Invisible \code{NULL}. The function prints to the console as a
#'   side effect.
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
    safe_fmt <- function(val, fmt = "g", digits = digits) {
      if (is.null(val) || length(val) == 0L || is.na(val[1L])) return("NA")
      formatC(val[1L], format = fmt, digits = digits)
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
        deviance_int  = safe_fmt(metrics$deviance_int,        fmt = "f", digits = digits),
        deviance_diff = safe_fmt(metrics$deviance_diff,       fmt = "f", digits = digits),
        df            = safe_fmt(metrics$df_int,      fmt = "d", digits = 0),
        pvalue        = safe_fmt(metrics$pvalue,              fmt = "g", digits = digits),
        stringsAsFactors = FALSE
      ),
      "aic"    = data.frame(
        AIC_main        = safe_fmt(metrics$AIC_main,           fmt = "f", digits = digits),
        AIC_interaction = safe_fmt(metrics$AIC_interaction,    fmt = "f", digits = digits),
        dAIC            = safe_fmt(metrics$AIC_main_minus_int, fmt = "f", digits = digits),
        stringsAsFactors = FALSE
      ),
      "bic"    = data.frame(
        BIC_main        = safe_fmt(metrics$BIC_main,           fmt = "f", digits = digits),
        BIC_interaction = safe_fmt(metrics$BIC_interaction,    fmt = "f", digits = digits),
        dBIC            = safe_fmt(metrics$BIC_main_minus_int, fmt = "f", digits = digits),
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

#' Evaluate Univariable Interactions Between a Grouping and Continuous Variables
#'
#' For each variable in `cont_vars`, fits the interaction model using the
#' pre-specified functional form from `cont_var_forms` and tests it against
#' the corresponding main-effects model according to `criterion`. Variables
#' whose model does not meet the selection threshold are excluded from the
#' output. No data-driven selection among functional forms is performed.
#'
#' Optionally, a precomputed adjustment matrix can be supplied via `xadj` with
#' `skip_adjustment = TRUE` to avoid redundant transformation of the adjustment
#' predictors when the function is called repeatedly with the same adjustment set.
#'
#' @param y Response vector or Surv object.
#' @param processed_data List returned by \code{preprocess_data()}.
#' @param selected_vars Character vector of variable names retained by the
#'   adjustment model.
#' @param adj_fp_powers Named list of FP powers for the selected adjustment
#'   variables, as returned by \code{get_fp_powers()}.
#' @param adj_zero Named logical vector for selected adjustment variables,
#'   giving the final post-fit zero-handling flags from
#'   \code{adjustment_model$zero}. These are used when rebuilding the selected
#'   adjustment matrix inside \code{evaluate_interactions()}.
#' @param adj_catzero Named logical vector for selected adjustment variables,
#'   giving the final post-fit structural-zero indicator flags from
#'   \code{adjustment_model$catzero}.
#' @param adj_spike Named logical vector for selected adjustment variables,
#'   giving the final post-fit spike-at-zero flags from
#'   \code{adjustment_model$spike}.
#' @param adj_spike_decision Named numeric vector of spike-at-zero decisions
#'   for the selected adjustment variables, as stored in
#'   \code{adjustment_model$spike_dec} by \code{fit_mfp()}. Values:
#'   \code{1} = FP term + binary indicator, \code{2} = FP term only,
#'   \code{3} = binary indicator only. Used together with \code{adj_spike}
#'   when re-transforming adjustment variables via \code{transform_matrix()}.
#' @param cont_vars Character vector of continuous variables to test.
#' @param cont_var_forms Named character vector of length \code{length(cont_vars)}.
#'   Names are variable names from \code{cont_vars}; values are one of
#'   \code{"linear"}, \code{"fp1"}, or \code{"fp2"}. Each variable is tested
#'   using exactly its pre-specified form; no selection among forms is performed.
#'   Constructed and validated by \code{mfpi.default()} before being passed here.
#' @param flex Character string; `"flex1"`, `"flex2"`, `"flex3"`, or
#'   `"flex4"`.
#' @param weights Numeric vector of observation weights.
#' @param offset Numeric vector of linear-predictor offsets.
#' @param xorder Character string; entry order for MFP backfitting.
#' @param ties Character string; tie-handling for Cox models.
#' @param strata Integer stratum vector, or `NULL`.
#' @param use_ftest Logical. If \code{TRUE} and \code{family = "gaussian"}, use an F-test rather than a chi-square test. Applied to both adjustment-variable selection and the interaction test.
#' @param control Fitting control list.
#' @param nocenter Numeric vector for Cox centring suppression.
#' @param family Character string; regression family.
#' @param fp_powers Named list of candidate FP power sets for all predictors.
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param criterion Character string; `"pvalue"`, `"aic"`, or `"bic"`.
#' @param digits Positive integer. Number of digits for printed output.
#' @param scale Named numeric vector of per-variable scale factors. These are
#'   the same scale factors applied to \code{x} upstream in
#'   \code{mfpi.default()}. Passed to \code{fit_adjustment_model()} so that
#'   \code{fit_mfp()} can backscale before the final model fit (giving
#'   coefficients on the \eqn{\phi(x + \text{shift})} scale), and to
#'   \code{evaluate_interactions()} for the same purpose in the interaction
#'   model. Default \code{NULL} (no backscaling).
#' @param shift Named numeric vector of per-variable shift factors. Passed to
#'   fitted-function generation so x coordinates can be returned/displayed on
#'   the original raw scale. Default \code{NULL} (no shift correction).
#' @param center_type Character string; \code{"grand"} (default) or
#'   \code{"group"}. Centering strategy for FP-transformed interaction
#' @param show_models Logical. If \code{TRUE} and \code{verbose = TRUE},
#'   prints the full regression coefficient table for each selected
#'   interaction model. Has no effect when \code{verbose = FALSE}.
#'   Nothing is printed when no variable is selected.
#' @param group_var Character string. Name of the grouping variable.
#' @param p_interact Numeric. Significance threshold for `criterion = "pvalue"`.
#' @param min_improvement Numeric. Minimum criterion improvement to retain an
#'   interaction.
#' @param xadj Optional numeric matrix of pre-transformed adjustment predictors.
#'   Used when `skip_adjustment = TRUE`.
#' @param skip_adjustment Logical. If `TRUE`, `xadj` is used directly and
#'   adjustment variables are not re-transformed. Default is `FALSE`.
#' @param p_adjust_method Character string. Method for adjusting p-values,
#'   passed to [stats::p.adjust()]. Only applied when `criterion = "pvalue"`.
#'   Default `"none"`. Adjusted p-values are applied across all variables in
#'   \code{cont_vars} (one test per variable).
#' @param quiet Logical. If `TRUE`, suppresses the per-variable message when no
#'   significant interaction is found. Default is `FALSE`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`best_model_metrics`}{Data frame of evaluation metrics for the
#'     selected interaction model per selected variable. Empty if none selected.}
#'   \item{`all_model_metrics`}{Data frame of metrics for the interaction model
#'     evaluated for each variable in `cont_vars`, using the form specified
#'     by `cont_var_forms`.}
#'   \item{`best_interaction_model`}{Named list of fitted model objects for
#'     selected variables.}
#'   \item{`all_interaction_models`}{Named list (one element per `cont_var`)
#'     of model objects, each fitted using the form from `cont_var_forms`.}
#'   \item{`center_vals_list`}{Named list of centering constants per
#'     selected variable.}
#'   \item{`var_winners`}{Named list (one element per `cont_var`) storing
#'     the fitted result for each variable regardless of selection.}
#'   \item{`group_var`}{Name of the grouping variable.}
#'   \item{`show_models`}{Value of the `show_models` argument.}
#'   \item{`group_levels_new`}{Integer-coded levels used internally.}
#'   \item{`group_levels_original`}{Original levels of `group_var`.}
#'   \item{`flex`}{The chosen flexibility level.}
#'   \item{`criterion`}{The selection criterion used.}
#'   \item{`family`}{The regression family.}
#'   \item{`nobs`}{Number of observations.}
#'   \item{`p_adjust_method`}{The multiplicity adjustment method used.}
#' }
#'
#' @section Selection logic:
#' For each variable in \code{cont_vars}, exactly one interaction model is
#' fitted using the form pre-specified in \code{cont_var_forms}. The model is
#' tested against its main-effects counterpart:
#'
#' \describe{
#'   \item{\code{"pvalue"}}{The interaction is selected if its (possibly
#'     adjusted) p-value is strictly below \code{p_interact}.}
#'   \item{\code{"aic"}}{Interaction is selected if
#'     \code{AIC_main - AIC_interaction > min_improvement}.}
#'   \item{\code{"bic"}}{Interaction is selected if
#'     \code{BIC_main - BIC_interaction > min_improvement}.}
#' }
#'
#' When \code{p_adjust_method != "none"} and \code{criterion = "pvalue"},
#' selection is deferred until all variables are processed: raw p-values are
#' collected, adjusted via \code{\link[stats]{p.adjust}}, and the adjusted
#' p-values are compared to \code{p_interact}. For AIC/BIC criteria,
#' multiplicity adjustment does not affect the selection decision.
#'
#' @keywords internal
#' @noRd
evaluate_interactions <- function(y, processed_data, selected_vars,
                                  adj_fp_powers, 
                                  adj_zero, adj_catzero, adj_spike,
                                  adj_spike_decision,
                                  cont_vars, cont_var_forms, flex,
                                  weights, offset, xorder, ties, strata,
                                  use_ftest, control, nocenter, family,
                                  family_string, fp_powers, cycles, criterion,
                                  digits, group_var, show_models, p_interact,
                                  min_improvement, min_prop, max_prop,
                                  center_type = c("grand", "group"),
                                  scale = NULL, shift = NULL,
                                  adj_acd_parameter = NULL,
                                  xadj = NULL, skip_adjustment = FALSE,
                                  p_adjust_method = "none", has_offset,
                                  quiet = FALSE, verbose = FALSE) {
  
  center_type <- match.arg(center_type)
  
  # Input validation for skip_adjustment / xadj --------------------------------
  if (!is.logical(skip_adjustment) || length(skip_adjustment) != 1L) {
    stop("`skip_adjustment` must be a single logical value.", call. = FALSE)
  }
  if (skip_adjustment && !is.null(xadj)) {
    if (!is.matrix(xadj) || !is.numeric(xadj)) {
      stop("`xadj` must be a numeric matrix when `skip_adjustment = TRUE`.",
           call. = FALSE)
    }
  }
  
  x        <- processed_data$original_x
  nobs     <- nrow(x)
  cat_info <- processed_data$cat_info
  
  # Full reference table mapping form name to FP degree
  degree_lookup <- c(linear = 0L, fp1 = 1L, fp2 = 2L)
  
  # Pre-allocate result containers with known maximum sizes.
  # Each variable contributes exactly one form (pre-specified via cont_var_forms),
  # so the total number of candidate rows equals n_vars.
  n_vars            <- length(cont_vars)
  n_types           <- 1L   # one pre-specified form per variable
  
  best_metrics_list      <- vector("list", n_vars)
  all_metrics_list       <- vector("list", n_vars * n_types)
  best_interaction_model <- vector("list", n_vars)   # winner only
  all_interaction_models <- vector("list", n_vars)   # one model per variable
  center_vals_list       <- vector("list", n_vars)   # centering constants per var
  
  names(best_interaction_model) <- cont_vars
  names(all_interaction_models) <- cont_vars
  names(center_vals_list)       <- cont_vars
  
  # Per-variable winners: stores the best candidate for EVERY variable
  # (regardless of significance) for use in phase-2 multiplicity adjustment.
  var_winners <- vector("list", n_vars)
  names(var_winners) <- cont_vars
  
  # Developer note:
  # Only defer p-value selection when an actual multiplicity adjustment is
  # requested. With p_adjust_method = "none", raw p-values are final and the
  # usual per-variable threshold can be applied immediately.
  adjusting_pvals <- criterion == "pvalue" && p_adjust_method != "none"
  
  best_idx     <- 0L   # counter for significant variables
  all_idx      <- 0L   # write counter for all candidate metric rows
  
  # Used only when p-value adjustment is requested. In that case the adjusted
  # candidate table becomes the authoritative all-model metric table returned
  # by evaluate_interactions().
  all_candidate_metrics <- NULL
  
  for (var_name in cont_vars) {
    
    if (verbose) {
      var_idx <- match(var_name, cont_vars)
      rule_thin <- strrep("-", 70)
      cat("\n", rule_thin, "\n", sep = "")
      cat(sprintf(
        "  [%d/%d] Variable: '%s' x '%s'\n",
        var_idx, n_vars, var_name, processed_data$cat_info$group_var
      ))
      # Adjustment variables for this evaluation: selected adjustment set
      # minus the current variable (main effect of the interaction) and
      # minus the group dummies (structural part of the interaction model).
      group_dummy_names <- colnames(processed_data$cat_info$dummies)
      adj_vars <- setdiff(selected_vars, c(var_name, group_dummy_names))
      if (length(adj_vars) > 0L) {
        cat(sprintf("        Adjusting for: %s\n",
                    paste(adj_vars, collapse = ", ")))
      } else {
        cat("        Adjusting for: (none)\n")
      }
      cat(rule_thin, "\n", sep = "")
    }
    
    best_fit    <- NULL
    best_metric <- NULL
    best_type   <- NULL
    best_score  <- Inf
    
    # The form for this variable is pre-specified; exactly one candidate.
    active_type <- cont_var_forms[[var_name]]
    
    # Collect the single candidate model and metrics for this variable
    candidate_models <- vector("list", n_types)
    candidate_metrics <- vector("list", n_types)
    candidate_full_fits <- vector("list", n_types)
    
    # Verbose: one row for the single candidate.
    if (verbose) {
      verbose_rows <- vector("list", n_types)
      names(verbose_rows) <- active_type
    }
    
    # Build adjustment matrix once per cont_var. The adjustment set changes only
    # by removing var_name from selected_vars, so it is identical regardless of
    # functional form and only needs to be built once per variable.
    if (!skip_adjustment) {
      #adj_vars     <- setdiff(selected_vars, var_name)
      group_dummy_names <- colnames(processed_data$cat_info$dummies)
      
      adj_vars <- setdiff(
        selected_vars,
        c(var_name, group_dummy_names)
      )
      
      xadj_current <- NULL
      if (length(adj_vars) > 0L) {
        # Backscale adj_vars to match fit_mfp convention: x was pre-scaled in
        # mfpi.default. Multiplying by scale restores x + shift so that
        # adjustment variable columns are on the phi(x + shift) scale,
        # consistent with the interaction variable and the adjustment model.
        x_adj <- x[, adj_vars, drop = FALSE]
        if (!is.null(scale)) {
          shared_adj <- intersect(adj_vars, names(scale))
          if (length(shared_adj) > 0L) {
            scale_adj <- scale[shared_adj]
            if (any(scale_adj != 1)) {
              x_adj[, shared_adj] <- sweep(
                x_adj[, shared_adj, drop = FALSE], 2L, scale_adj, "*"
              )
            }
          }
        }
        xadj_current <- transform_matrix(
          x                  = x_adj,
          power_list         = adj_fp_powers[adj_vars],
          center             = processed_data$updated_params$center[adj_vars],
          acdx               = processed_data$updated_params$acd_vars[adj_vars],
          zero               = adj_zero[adj_vars],
          catzero            = adj_catzero[adj_vars],
          spike              = adj_spike[adj_vars],
          spike_decision     = adj_spike_decision[adj_vars],
          keep_x_order       = FALSE,
          acd_parameter_list = if (!is.null(adj_acd_parameter)) adj_acd_parameter[adj_vars] else NULL,
          reset_zero         = FALSE,
          check_binary       = TRUE
        )$x_transformed
      }
    } else {
      xadj_current <- xadj
    }
    
    for (interaction_type in active_type) {
      degree <- degree_lookup[[interaction_type]]
      
      # Fit interaction model --------------------------------------------------
      fit_result <- flex_fit(
        x             = x,
        y             = y,
        cont_var      = var_name,
        group_var     = cat_info$group_var,
        group_dummies = cat_info$dummies,
        xadj          = xadj_current,
        criterion     = criterion,
        ties          = ties,
        degree        = degree,
        family        = family,
        family_string = family_string,
        fp_cand       = processed_data$updated_params$fp_powers[[var_name]],
        use_ftest     = use_ftest,
        center        = processed_data$updated_params$center[var_name],
        center_type   = center_type,
        xorder        = xorder,
        weights       = weights,
        offset        = offset,
        strata        = strata,
        control       = control,
        nocenter      = nocenter,
        cycles        = cycles,
        zero_var      = processed_data$updated_params$zero_vars[var_name],
        spike_var     = FALSE,
        min_prop      = min_prop,
        max_prop      = max_prop,
        flex          = flex,
        scale_var     = if (!is.null(scale)) unname(scale[var_name]) else 1,
        shift_var     = if (!is.null(shift)) unname(shift[var_name]) else 0,
        run_test      = TRUE,
        has_offset    = has_offset
      )
      
      # Annotate metrics -------------------------------------------------------
      metrics          <- fit_result$test_results$evaluation_metrics
      metrics$variable <- var_name
      metrics$type     <- interaction_type
      
      # Normalise reporting-oriented FP power columns. This keeps the metric object
      # readable and machine-readable without changing the fitted interaction model
      # or prediction internals.
      metrics <- normalise_metric_fp_powers(
        metrics      = metrics,
        variable     = var_name,
        group_levels = cat_info$original_levels
      )
      
      candidate_metrics[[interaction_type]] <- metrics
      candidate_full_fits[[interaction_type]] <- fit_result
      
      # Collect verbose row for this candidate (printed as table after inner loop)
      if (verbose) {
        verbose_rows[[interaction_type]] <- metrics
      }
      
      all_idx <- all_idx + 1L
      all_metrics_list[[all_idx]] <- metrics
      
      # Store candidate model regardless of significance
      candidate_models[[interaction_type]] <- fit_result$test_results$interaction_model
      
      # With one pre-specified form per variable, assign directly.
      # For AIC/BIC the winner is resolved in the block below after the
      # loop; for pvalue it is set here.
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
          paste0("! `criterion` must be one of 'pvalue', 'aic', or 'bic'; got '",
                 criterion, "'."),
          call. = FALSE
        )
      }
    }  # end interaction_type loop
    
    if (criterion %in% c("aic", "bic")) {
      m <- candidate_metrics[[active_type]]
      
      if (!is.null(m)) {
        ic_col <- if (criterion == "aic") {
          "AIC_main_minus_int"
        } else {
          "BIC_main_minus_int"
        }
        
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
    
    # Print the per-variable metrics table and per-variable result line
    if (verbose) {
      format_candidate_table(verbose_rows, criterion = criterion,
                             digits = digits, best_type = best_type)
      if (criterion == "pvalue" && !is.null(best_metric)) {
        pval_str <- formatC(best_metric$pvalue[1L], format = "g", digits = digits)
        if (adjusting_pvals) {
          # p_adjust_method != "none": final decision deferred to Phase 2
          cat(sprintf(
            "        >> provisional: p_raw = %s  (final after adjustment)\n",
            pval_str
          ))
        } else {
          # p_adjust_method = "none": this is the final result
          passes <- best_metric$pvalue[1L] < p_interact
          cat(sprintf(
            "        >> p_raw = %s  \u2192  %s\n",
            pval_str,
            if (passes) "Selected" else "Not selected"
          ))
        }
      }
    }
    
    # Store the per-variable winner (may be NULL if the model failed to fit)
    var_winners[[var_name]] <- list(
      fit         = best_fit,
      metric      = best_metric,
      type        = best_type,
      score       = best_score,
      center_vals = if (!is.null(best_fit)) best_fit$center_vals else NULL
    )
    
    # ------------------------------------------------------------------
    # Selection decision (when NOT deferring to Phase 2)
    # ------------------------------------------------------------------
    if (!adjusting_pvals && !is.null(best_fit)) {
      # Check whether the best candidate passes the threshold
      passes_threshold <- if (criterion == "pvalue") {
        best_score < p_interact
      } else {
        best_score > min_improvement
      }
      
      if (passes_threshold) {
        best_idx <- best_idx + 1L
        best_metrics_list[[best_idx]]      <- best_metric
        best_interaction_model[[var_name]] <- best_fit$test_results$interaction_model
        center_vals_list[[var_name]]       <- best_fit$center_vals
        
        if (verbose && criterion %in% c("aic", "bic")) {
          ic_col <- if (criterion == "aic") {
            "AIC_main_minus_int"
          } else {
            "BIC_main_minus_int"
          }
          
          ic_label <- if (criterion == "aic") "dAIC" else "dBIC"
          
          ic_val <- if (ic_col %in% names(best_metric)) {
            best_metric[[ic_col]][1L]
          } else {
            NA_real_
          }
          
          cat(sprintf(
            "        >> Selected  (%s = %s)\n",
            ic_label,
            formatC(ic_val, format = "f", digits = 2)
          ))
        }
        
        if (verbose && show_models) {
          fit_obj <- best_fit$test_results$interaction_model$fit
          if (!is.null(fit_obj)) {
            type_label <- format_type_label(best_type)
            cat(sprintf(
              "\n  Interaction model (%s, %s):\n",
              var_name, type_label
            ))
            cat(strrep("-", 50), "\n")
            print(summary(fit_obj))
            cat("\n")
          }
        }
      } else {
        if (verbose && criterion %in% c("aic", "bic")) {
          ic_col <- if (criterion == "aic") {
            "AIC_main_minus_int"
          } else {
            "BIC_main_minus_int"
          }
          
          ic_label <- if (criterion == "aic") "dAIC" else "dBIC"
          
          ic_val <- if (ic_col %in% names(best_metric)) {
            best_metric[[ic_col]][1L]
          } else {
            NA_real_
          }
          
          cat(sprintf(
            "        >> Not selected  (%s = %s)\n",
            ic_label,
            formatC(ic_val, format = "f", digits = 2)
          ))
        }
        
        if (!quiet) {
          message(sprintf(
            "i No significant interaction retained for '%s'.", var_name
          ))
        }
      }
    } else if (!adjusting_pvals && is.null(best_fit)) {
      if (verbose) {
        cat("        >> Not selected: model failed to fit\n")
      }
      if (!quiet) {
        message(sprintf(
          "i No significant interaction retained for '%s'.", var_name
        ))
      }
    }
    
    # Store the fitted interaction model for this variable regardless of significance
    all_interaction_models[[var_name]] <- candidate_models
    
    # Update all_metrics_list with the selected best_metric. For AIC/BIC,
    # the canonical improvement columns remain `AIC_main_minus_int` and
    # `BIC_main_minus_int`; `dAIC` and `dBIC` are display labels only.
    if (!is.null(best_metric)) {
      all_metrics_list[[all_idx]] <- best_metric
    }
  }  # end cont_vars loop
  
  # ============================================================================
  # Step 3 summary for AIC/BIC (printed after all variables are processed)
  # ============================================================================
  if (verbose && criterion %in% c("aic", "bic") && !adjusting_pvals) {
    rule_thick <- strrep("=", 70)
    rule_thin  <- strrep("-", 70)
    
    ic_col <- if (criterion == "aic") {
      "AIC_main_minus_int"
    } else {
      "BIC_main_minus_int"
    }
    
    ic_label <- if (criterion == "aic") "dAIC" else "dBIC"
    
    cat("\n", rule_thick, "\n", sep = "")
    cat("  STEP 3: Interaction Summary\n")
    cat(rule_thick, "\n", sep = "")
    cat(sprintf("  %-12s %-8s %12s  %s\n", "Variable", "Type", ic_label, ""))
    cat(rule_thin, "\n", sep = "")
    
    for (vn in cont_vars) {
      w <- var_winners[[vn]]
      
      if (is.null(w$fit) || is.null(w$metric)) {
        cat(sprintf("  %-12s %-8s %12s\n", vn, "---", "NA"))
      } else {
        type_label <- format_type_label(w$type)
        
        ic_val <- if (ic_col %in% names(w$metric)) {
          w$metric[[ic_col]][1L]
        } else {
          NA_real_
        }
        
        sel_flag <- if (!is.na(w$score) && w$score > min_improvement) " *" else ""
        
        cat(sprintf(
          "  %-12s %-8s %12s%s\n",
          vn,
          type_label,
          formatC(ic_val, format = "f", digits = 2),
          sel_flag
        ))
      }
    }
    
    cat(rule_thin, "\n", sep = "")
    cat(sprintf("  * = selected (%s > %g)\n\n", ic_label, min_improvement))
    
    if (show_models) {
      for (vn in names(Filter(Negate(is.null), best_interaction_model))) {
        w <- var_winners[[vn]]
        fit_obj <- w$fit$test_results$interaction_model$fit
        if (!is.null(fit_obj)) {
          type_label <- format_type_label(w$type)
          cat(sprintf("  Interaction model (%s, %s):\n", vn, type_label))
          cat(strrep("-", 50), "\n")
          print(summary(fit_obj))
          cat("\n")
        }
      }
    }
  }
  
  # ============================================================================
  # Phase 2: p-value adjustment and final p-value decisions
  # ============================================================================
  if (adjusting_pvals) {
    #all_candidate_metrics <- dplyr::bind_rows(all_metrics_list[seq_len(all_idx)])
    all_candidate_metrics <- bind_metric_rows(all_metrics_list[seq_len(all_idx)])
    if (nrow(all_candidate_metrics) > 0L && "pvalue" %in% names(all_candidate_metrics)) {
      all_candidate_metrics$p_adjusted <- NA_real_
      
      valid_p <- is.finite(all_candidate_metrics$pvalue)
      n_planned_tests <- length(cont_vars) * n_types
      
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
      
      # Update each variable's winner with the adjusted p-value
      for (vn in cont_vars) {
        rows_idx <- which(all_candidate_metrics$variable == vn &
                            is.finite(all_candidate_metrics$p_adjusted))
        if (length(rows_idx) == 0L) next
        # Only one row per variable; rows_idx has length 1
        win_idx <- rows_idx[1L]
        w <- var_winners[[vn]]
        if (!is.null(w$fit)) {
          w$metric <- all_candidate_metrics[win_idx, , drop = FALSE]
          w$score  <- all_candidate_metrics$p_adjusted[win_idx]
          var_winners[[vn]] <- w
        }
      }
    }
    
    # Reset final selected objects and fill from the adjusted variable winners.
    best_idx <- 0L
    best_metrics_list <- vector("list", n_vars)
    best_interaction_model <- stats::setNames(vector("list", n_vars), cont_vars)
    center_vals_list       <- stats::setNames(vector("list", n_vars), cont_vars)
    
    for (vn in cont_vars) {
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
    
    
    # Verbose: after p-value adjustment, print one final interaction summary.
    # The per-variable output printed during fitting is provisional because
    # adjusted p-values are only available after all variables are fitted.
    if (verbose) {
      rule_thick <- strrep("=", 70)
      rule_thin  <- strrep("-", 70)
      
      cat("\n", rule_thick, "\n", sep = "")
      cat("  STEP 3: Interaction Summary\n")
      cat(rule_thick, "\n", sep = "")
      cat(sprintf("  %-12s %-8s %12s %12s  %s\n",
                  "Variable", "Type", "p_raw", "p_adjusted", ""))
      cat(rule_thin, "\n", sep = "")
      
      for (vn in cont_vars) {
        w <- var_winners[[vn]]
        if (is.null(w$fit) || is.null(w$metric)) {
          cat(sprintf("  %-12s %-8s %12s %12s\n", vn, "---", "NA", "NA"))
        } else {
          type_label <- format_type_label(w$type)
          sel_flag <- if (!is.na(w$metric$p_adjusted[1L]) &&
                          w$metric$p_adjusted[1L] < p_interact) " *" else ""
          cat(sprintf("  %-12s %-8s %12s %12s%s\n",
                      vn, type_label,
                      formatC(w$metric$pvalue[1L], format = "g", digits = digits),
                      formatC(w$metric$p_adjusted[1L], format = "g", digits = digits),
                      sel_flag))
        }
      }
      cat(rule_thin, "\n", sep = "")
      cat(sprintf("  p_adjust_method = %s\n", p_adjust_method))
      cat(sprintf("  * = selected at p_interact = %g\n\n", p_interact))
      
      if (show_models) {
        for (vn in names(Filter(Negate(is.null), best_interaction_model))) {
          w <- var_winners[[vn]]
          fit_obj <- w$fit$test_results$interaction_model$fit
          if (!is.null(fit_obj)) {
            type_label <- format_type_label(w$type)
            cat(sprintf("  Interaction model (%s, %s):\n", vn, type_label))
            cat(strrep("-", 50), "\n")
            print(summary(fit_obj))
            cat("\n")
          }
        }
      }
    }
    
  }
  
  # Trim pre-allocated lists to actual fill level. When p-value adjustment was
  # used, all_candidate_metrics is the authoritative all-model metric table, so
  # all_metrics_list no longer needs to be trimmed for final binding.
  best_metrics_list <- best_metrics_list[seq_len(best_idx)]
  
  if (!adjusting_pvals) {
    all_metrics_list <- all_metrics_list[seq_len(all_idx)]
  }
  
  # Drop NULL slots from best_interaction_model
  best_interaction_model <- Filter(Negate(is.null), best_interaction_model)
  
  # Combine results -----------------------------------------------------------
  combined_metrics <- bind_metric_rows(best_metrics_list)
  
  if (nrow(combined_metrics) > 0L && "type" %in% names(combined_metrics)) {
    # Move type to first column using base R subsetting.
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
  
  # type first, variable second, then remaining metric columns
  if (nrow(all_model_metrics) > 0L) {
    metric_cols       <- setdiff(names(all_model_metrics), c("variable", "type"))
    all_model_metrics <- all_model_metrics[, c("type", "variable", metric_cols),
                                           drop = FALSE]
  }
  
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