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
#'     \item \strong{Interaction functional form selection} (Step 2): selects
#'       among linear, FP1, and FP2 interaction candidates in
#'       \code{evaluate_interactions()}.
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
#' @param force_keep Character vector of variable names forced into the
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
#' @param digits Positive integer. Significant digits for printed output.
#'
#' @return A list with top-level fields (accessible as \code{fit$field}) and
#'   two retained sub-objects. Top-level fields:
#' \describe{
#'   \item{\code{best_model_metrics}}{Tibble of metrics for the best
#'     interaction model per significant variable.}
#'   \item{\code{all_model_metrics}}{Tibble of metrics for all three candidates
#'     (linear, FP1, FP2) for every variable in \code{cont_vars}.}
#'   \item{\code{best_interaction_model}}{Named list of fitted interaction model
#'     objects for the **winning** model per significant variable.}
#'   \item{\code{all_interaction_models}}{Named list of lists. Outer names are
#'     \code{cont_vars}; inner names are \code{"linear"}, \code{"fp1"},
#'     \code{"fp2"}. Each leaf is the fitted model object for that candidate,
#'     regardless of whether it was selected or significant. Allows the user
#'     to inspect or compare all three functional forms for every variable.}
#'   \item{\code{best_fitted_functions}}{Named list of group-specific fitted values.}
#'   \item{\code{adjust_terms}}{The \code{fp_terms} table from the adjustment
#'     model.}
#'   \item{\code{group_var}, \code{flex}, \code{family}, \code{nobs},
#'     \code{show_models}, \code{group_levels_new},
#'     \code{group_levels_original}}{Metadata fields.}
#'   \item{\code{adjustment_model}}{Full adjustment model object.}
#'   \item{\code{univariable_interactions}}{Full list from
#'     \code{evaluate_interactions()}.}
#' }
#'
#' @section Algorithm:
#' \strong{Pre-processing.} Levels of \code{group_var} are remapped to
#' consecutive integers. Group dummies are computed once and stored in
#' \code{cat_info$dummies}. When \code{include_group_var = TRUE}, dummies are
#' appended to \code{x}. All parameter vectors are synchronized Silently.
#'
#' \strong{Step 1 - Adjustment model.} \code{fit_mfp()} selects
#' adjustment variables and FP transformations using \code{criterion}.
#' Selected variables, their FP powers, and spike decisions are fixed for
#' the remainder of the algorithm.
#'
#' \strong{Step 2 - Univariable interactions.} For each variable in
#' \code{cont_vars}, linear, FP1, and FP2 interaction models are fitted and
#' compared to their main-effects models via \code{evaluate_interactions()}.
#' Within \code{flex1()} and \code{flex4()}, FP powers are estimated using
#' \code{fit_mfp()} with \code{force_max_fp = TRUE} to prevent AIC/BIC
#' from simplifying the functional form below the requested degree. Functional
#' form selection (linear/FP1/FP2) is performed by \code{evaluate_interactions()}
#' using the user \code{criterion}.
#'
#' @seealso \code{mfpi()}, \code{preprocess_data()},
#'   \code{fit_adjustment_model()}, \code{evaluate_interactions()}
#'
#' @keywords internal
#' @noRd
fit_mfpi <- function(x, y, family, family_string, weights, offset, cycles,
                     center, criterion, select, alpha, force_keep,
                     force_max_fp,
                     df, xorder,
                     fp_powers, ties, strata, nocenter, acd_vars, zero_vars,
                     catzero_vars, spike_vars, min_prop, max_prop, use_ftest,
                     control, group_var, include_group_var, flex, cont_vars,
                     p_interact, min_improvement, show_models, verbose,
                     digits, scale,
                     center_type = c("grand", "group")) {
  
  center_type <- match.arg(center_type)
  
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
    force_keep        = force_keep,
    force_max_fp      = force_max_fp
  )
  
  # ---------------------------------------------------------------------------
  # Step 1: Fit adjustment model via MFP
  # ---------------------------------------------------------------------------
  if (verbose) {
    rule_thick <- strrep("=", 70)
    #rule_thin  <- strrep("-", 70)
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
    family         = family,
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
  
  adj_fp_powers <- get_fp_powers(
    selected_vars,
    adjustment_model$fp_terms
  )
  
  # Extract spike decisions from the adjustment model for use when
  # re-transforming adjustment variables inside evaluate_interactions().
  # spike_decision is a named numeric vector (1 = FP+binary, 2 = FP only,
  # 3 = binary only) stored in the mfp2 object by mfp2:::fit_mfp().
  adj_spike_decision <- if (!is.null(adjustment_model$spike_dec)) {
    adjustment_model$spike_decision
   } else {
    setNames(rep(2L, length(selected_vars)), selected_vars)
   }
  
  # ---------------------------------------------------------------------------
  # Step 2: Evaluate univariable interactions
  # ---------------------------------------------------------------------------
  if (verbose) {
    #rule_thick <- strrep("=", 70)
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
    adj_spike_decision = adj_spike_decision,
    cont_vars          = cont_vars,
    flex               = flex,
    weights            = weights,
    offset             = offset,
    xorder             = xorder,
    ties               = ties,
    strata             = strata,
    use_ftest          = use_ftest,
    control            = control,
    nocenter           = nocenter,
    family             = family,
    family_string      = family_string,
    fp_powers          = fp_powers,
    cycles             = cycles,
    criterion          = criterion,
    digits             = digits,
    group_var          = group_var,
    show_models        = show_models,
    p_interact         = p_interact,
    min_improvement    = min_improvement,
    min_prop           = min_prop,
    max_prop           = max_prop,
    center_type        = center_type,
    scale              = scale,
    verbose            = verbose
  )
  
  n_significant <- nrow(univ_results$best_model_metrics)
  if (n_significant == 0L) {
    warning("! No significant interactions were found.", call. = FALSE)
  }
  
  if (verbose) {
    #rule_thick <- strrep("=", 70)
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
    best_fitted_functions    = univ_results$best_fitted_functions,
    all_fitted_functions     = univ_results$all_fitted_functions,
    center_vals_list         = univ_results$center_vals_list,
    adjust_terms             = adjustment_model$fp_terms,
    group_var                = univ_results$group_var,
    show_models              = univ_results$show_models,
    flex                     = univ_results$flex,
    criterion                = univ_results$criterion,
    group_levels_new         = univ_results$group_levels_new,
    group_levels_original    = univ_results$group_levels_original,
    family                   = univ_results$family,
    nobs                     = univ_results$nobs,
    adjustment_model         = adjustment_model,
    univariable_interactions = univ_results
  )
}


# -----------------------------------------------------------------------------
# preprocess_data() -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Pre-process Predictor Data for the Adjustment Model
#'
#' Removes `group_var` from the predictor matrix and, when
#' `include_group_var = TRUE`, appends orthogonal dummy variables for that
#' variable. All modeling-parameter vectors (`select`, `alpha`, `df`, etc.) are
#' synchronised with the modified column set so that every downstream function
#' receives consistently aligned inputs.
#'
#' This function is called at the start of \code{fit_mfpi()} and its output is
#' threaded through the entire algorithm. It does not perform argument checking.
#'
#' @param x Numeric matrix of predictors including `group_var`.
#' @param group_var Character string. Name of the categorical grouping variable.
#' @param include_group_var Logical. If `TRUE`, dummy variables derived from
#'   `group_var` are appended to the filtered predictor matrix so that group
#'   membership can be adjusted for in the MFP step.
#' @param select Named numeric vector of selection levels, one per column of
#'   `x`.
#' @param alpha Named numeric vector of FP-degree significance levels, one per
#'   column of `x`.
#' @param df Named numeric vector of degrees of freedom, one per column of `x`.
#' @param center Named logical vector. Whether to centre each predictor.
#' @param acd_vars Named logical vector. Whether each predictor undergoes the
#'   ACD transformation.
#' @param fp_powers Named list of candidate FP power sets, one per predictor.
#' @param zero_vars Named logical vector. Whether each predictor should treat
#'   non-positive values as zero.
#' @param force_keep Character vector of variable names always retained in the
#'   adjustment model.
#' @param force_max_fp Named logical vector of length \eqn{p}. Sliced to
#'   adjustment columns and stored in \code{updated_params$force_max_fp}.
#'   When \code{include_group_var = TRUE}, dummy columns receive
#'   \code{FALSE} (they have \code{df = 1} and no functional form to force).
#'
#' @return A list with five components:
#' \describe{
#'   \item{`x`}{Numeric matrix with `group_var` removed and, optionally, dummy
#'     columns appended.}
#'   \item{`original_x`}{The unmodified input matrix.}
#'   \item{`cat_info`}{List of metadata for `group_var`:
#'     `values` (original column), `original_levels` (sorted unique values),
#'     `new_levels` (0-based integer codes), and `group_var` (the name).}
#'   \item{`dummy_names`}{Character vector of appended dummy column names, or
#'     `NULL` if `include_group_var = FALSE`.}
#'   \item{`updated_params`}{Named list of parameter vectors aligned with the
#'     columns of the returned `x`.}
#' }
#'
#' @keywords internal
#' @noRd
preprocess_data <- function(x, group_var, include_group_var,
                            select, alpha, df, center, acd_vars,
                            fp_powers, zero_vars, catzero_vars, spike_vars,
                            force_keep, force_max_fp) {
  
  original_x     <- x
  original_names <- colnames(x)
  
  # Remap group_var levels to 0-based integers
  group_values    <- x[, group_var, drop = FALSE]
  original_levels <- sort(unique(drop(group_values)))
  new_levels      <- seq_along(original_levels) - 1L
  
  # Remove group_var from the predictor set
  predictor_names <- setdiff(original_names, group_var)
  x_filtered      <- x[, predictor_names, drop = FALSE]
  
  # Helper: extract and rename a parameter vector to match the filtered columns
  pred_idx    <- which(original_names %in% predictor_names)
  slice_param <- function(v) setNames(v[pred_idx], predictor_names)
  
  updated_params <- list(
    select       = slice_param(select),
    alpha        = slice_param(alpha),
    df           = slice_param(df),
    center       = slice_param(center),
    acd_vars     = slice_param(acd_vars),
    fp_powers    = slice_param(fp_powers),
    zero_vars    = slice_param(zero_vars),
    catzero_vars = slice_param(catzero_vars),
    spike_vars   = slice_param(spike_vars),
    force_max_fp = slice_param(force_max_fp),
    force_keep   = force_keep
  )
  
  dummy_names <- NULL
  
  # Compute group dummies once - reused for cat_info and (if needed) appended
  # to x_filtered when include_group_var = TRUE.
  group_dummies_mat <- create_group_dummies(group_values)
  
  if (include_group_var) {
    dummy_names <- colnames(group_dummies_mat)
    x_filtered  <- cbind(x_filtered, group_dummies_mat)
    n_dummies   <- length(dummy_names)
    
    make_dummy_param <- function(val) setNames(rep(val, n_dummies), dummy_names)
    
    dummy_params <- list(
      select       = make_dummy_param(1),    # always selected
      alpha        = make_dummy_param(0),    # no FP transformation
      df           = make_dummy_param(1L),   # linear (binary)
      center       = make_dummy_param(FALSE),
      acd_vars     = make_dummy_param(FALSE),
      zero_vars    = make_dummy_param(FALSE),
      catzero_vars = make_dummy_param(FALSE),
      spike_vars   = make_dummy_param(FALSE),
      force_max_fp = make_dummy_param(FALSE),  # dummies: df=1, nothing to force
      fp_powers    = as.list(setNames(rep(1, n_dummies), dummy_names))
    )
    
    for (param in names(dummy_params)) {
      updated_params[[param]] <- c(updated_params[[param]], dummy_params[[param]])
    }
    updated_params$force_keep <- c(force_keep, dummy_names)
  }
  
  list(
    x           = x_filtered,
    original_x  = original_x,
    cat_info    = list(
      values          = group_values,
      original_levels = original_levels,
      new_levels      = new_levels,
      group_var       = group_var,
      dummies         = group_dummies_mat   # computed once above, reused here
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
#'   \code{center}, \code{select}, \code{alpha}, \code{force_keep},
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
                                 force_max_fp, scale,
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
    keep          = updated_params$force_keep,
    powers        = updated_params$fp_powers,
    acdx          = updated_params$acd_vars,
    zero          = updated_params$zero_vars,
    catzero       = updated_params$catzero_vars,
    spike         = updated_params$spike_vars,
    min_prop      = min_prop,
    max_prop      = max_prop,
    force_max_fp  = force_max_fp,
    verbose       = verbose
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
#'     \code{df} (interaction degrees of freedom), \code{pvalue}, and
#'     \code{dBIC}. \code{dBIC} (\eqn{\mathrm{BIC}_{\mathrm{main}} -
#'     \mathrm{BIC}_{\mathrm{int}}}) is the tiebreaker when two candidates
#'     have identical p-values; the candidate with the larger \code{dBIC}
#'     wins, favouring the simpler functional form.
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
#' @param rows Named list of one-row data frames, one per candidate type
#'   (linear, fp1, fp2). Each row must contain fields \code{pvalue},
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
format_candidate_table <- function(rows, criterion, digits) {
  
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
    safe_fmt <- function(val, fmt = "g", digits = 3) {
      if (is.null(val) || length(val) == 0L || is.na(val[1L])) return("NA")
      formatC(val[1L], format = fmt, digits = digits)
    }
    safe_chr <- function(val) {
      if (is.null(val) || length(val) == 0L) return("")
      as.character(val[1L])
    }
    
    type_raw   <- safe_chr(metrics$type)
    type_label <- if (nzchar(type_raw)) toupper(gsub("fp", "FP", type_raw)) else "?"
    
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
        deviance_int  = safe_fmt(metrics$deviance_int,        fmt = "f", digits = 2),
        deviance_diff = safe_fmt(metrics$deviance_diff,       fmt = "f", digits = 2),
        df            = safe_fmt(metrics$df_interaction,      fmt = "d", digits = 0),
        pvalue        = safe_fmt(metrics$pvalue,              fmt = "g", digits = digits),
        dBIC          = safe_fmt(metrics$BIC_main_minus_int,  fmt = "f", digits = 2),
        stringsAsFactors = FALSE
      ),
      "aic"    = data.frame(
        AIC_main        = safe_fmt(metrics$AIC_main,           fmt = "f", digits = 2),
        AIC_interaction = safe_fmt(metrics$AIC_interaction,    fmt = "f", digits = 2),
        dAIC            = safe_fmt(metrics$AIC_main_minus_int, fmt = "f", digits = 2),
        stringsAsFactors = FALSE
      ),
      "bic"    = data.frame(
        BIC_main        = safe_fmt(metrics$BIC_main,           fmt = "f", digits = 2),
        BIC_interaction = safe_fmt(metrics$BIC_interaction,    fmt = "f", digits = 2),
        dBIC            = safe_fmt(metrics$BIC_main_minus_int, fmt = "f", digits = 2),
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
#' For each variable in `cont_vars`, fits linear, FP1, and FP2 interaction
#' models with `group_var` and selects the best functional form according to
#' `criterion` and `min_improvement`. Variables whose best model does not meet
#' the selection threshold are excluded from the output.
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
#' @param adj_spike_decision Named numeric vector of spike-at-zero decisions
#'   for the adjustment variables, as stored in \code{adjustment_model$spike_dec}
#'   by \code{fit_mfp()}. Values: \code{1} = FP term + binary indicator,
#'   \code{2} = FP term only (default), \code{3} = binary indicator only.
#'   Used when re-transforming adjustment variables via \code{transform_matrix()}
#'   to ensure binary indicators for semi-continuous adjustment variables are
#'   included or suppressed consistently with the adjustment model.
#' @param cont_vars Character vector of continuous variables to test.
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
#' @param digits Positive integer. Significant digits for printed output.
#' @param scale Named numeric vector of per-variable scale factors. These are
#'   the same scale factors applied to \code{x} upstream in
#'   \code{mfpi.default()}. Passed to \code{fit_adjustment_model()} so that
#'   \code{fit_mfp()} can backscale before the final model fit (giving
#'   coefficients on the \eqn{\phi(x + \text{shift})} scale), and to
#'   \code{evaluate_interactions()} for the same purpose in the interaction
#'   model. Default \code{NULL} (no backscaling).
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
#' @param quiet Logical. If `TRUE`, suppresses the per-variable message when no
#'   significant interaction is found. Default is `FALSE`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`best_model_metrics`}{Data frame of evaluation metrics for the
#'     best interaction model per variable. Empty if no interactions are
#'     significant.}
#'   \item{`all_model_metrics`}{Data frame of metrics for all models tested
#'     (linear, FP1, FP2) for every variable in `cont_vars`, regardless of
#'     significance.}
#'   \item{`best_interaction_model`}{Named list of fitted model objects for the best
#'     interaction type per significant variable.}
#'   \item{`best_fitted_functions`}{Named list of fitted FP function values per
#'     significant variable, structured for plotting.}
#'   \item{`group_var`}{Name of the grouping variable.}
#'   \item{`show_models`}{Value of the `show_models` argument.}
#'   \item{`group_levels_new`}{Integer-coded levels used internally.}
#'   \item{`group_levels_original`}{Original levels of `group_var`.}
#'   \item{`flex`}{The chosen flexibility level.}
#'   \item{`family`}{The regression family.}
#'   \item{`nobs`}{Number of observations.}
#' }
#'
#' @section Selection logic:
#' For each variable in \code{cont_vars}, three candidate interaction models
#' are fitted in order: linear (degree 0), FP1 (degree 1), FP2 (degree 2).
#' The models are evaluated using \code{criterion}:
#'
#' \describe{
#'   \item{\code{"pvalue"}}{Retain candidates with p-value strictly below
#'     \code{p_interact}. Among retained candidates, select the one with the
#'     \strong{smallest p-value}. Ties are broken by the \strong{largest
#'     \code{BIC_main_minus_int}}: when two candidates achieve identical
#'     p-values, the one with the larger BIC improvement wins. BIC penalizes
#'     complexity more heavily than AIC (penalty grows as
#'     \eqn{k\log n} rather than \eqn{2k}), so this rule favours the
#'     \strong{simpler} functional form when statistical evidence is equally
#'     strong - aligning with the parsimony principle of FP selection.}
#'   \item{\code{"aic"}}{Retain candidates with \code{AIC_main_minus_int}
#'     strictly above \code{min_improvement}. Select the one with the
#'     \strong{largest \code{AIC_main_minus_int}}.}
#'   \item{\code{"bic"}}{Same as \code{"aic"} using \code{BIC_main_minus_int}.}
#' }
#'
#' \code{AIC_main_minus_int} and \code{BIC_main_minus_int} are defined as
#' \eqn{\mathrm{AIC}_\text{main} - \mathrm{AIC}_\text{interaction}} and
#' \eqn{\mathrm{BIC}_\text{main} - \mathrm{BIC}_\text{interaction}}
#' respectively. A \strong{positive} value indicates that the interaction
#' model fits better (lower information criterion) than the main-effects
#' model.
#'
#' If no candidate clears its threshold, no interaction is retained for that
#' variable and it is excluded from \code{best_model_metrics}.
#'
#' @keywords internal
#' @importFrom dplyr bind_rows
#' @noRd
evaluate_interactions <- function(y, processed_data, selected_vars,
                                  adj_fp_powers, adj_spike_decision,
                                  cont_vars, flex,
                                  weights, offset, xorder, ties, strata,
                                  use_ftest, control, nocenter, family,
                                  family_string, fp_powers, cycles, criterion,
                                  digits, group_var, show_models, p_interact,
                                  min_improvement, min_prop, max_prop,
                                  center_type = c("grand", "group"),
                                  scale = NULL,
                                  xadj = NULL, skip_adjustment = FALSE,
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
  
  # Map interaction type to FP degree
  degree_lookup <- c(linear = 0L, fp1 = 1L, fp2 = 2L)
  
  # Pre-allocate result containers with known maximum sizes.
  # Indexed by position during the loop then trimmed at the end.
  n_vars            <- length(cont_vars)
  n_types           <- length(degree_lookup)   # 3: linear, fp1, fp2
  
  best_metrics_list      <- vector("list", n_vars)
  all_metrics_list       <- vector("list", n_vars * n_types)
  best_interaction_model <- vector("list", n_vars)   # winner only
  all_interaction_models <- vector("list", n_vars)   # all three candidates
  best_fitted_functions       <- vector("list", n_vars)   # winner fitted values
  all_fitted_functions   <- vector("list", n_vars)   # all three candidates
  center_vals_list       <- vector("list", n_vars)   # centering constants per var
  
  names(best_interaction_model) <- cont_vars
  names(all_interaction_models) <- cont_vars
  names(best_fitted_functions)       <- cont_vars
  names(all_fitted_functions)   <- cont_vars
  names(center_vals_list)       <- cont_vars
  
  best_idx     <- 0L   # counter for significant variables
  all_idx      <- 0L   # counter for all models
  
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
    best_score  <- if (criterion == "pvalue") Inf else -Inf
    
    # Collect all three candidate models and fitted functions for this variable
    candidate_models <- vector("list", n_types)
    candidate_fitted <- vector("list", n_types)
    names(candidate_models) <- names(degree_lookup)
    names(candidate_fitted) <- names(degree_lookup)
    
    # Verbose: collect one row per candidate, printed as table after inner loop.
    # Names required so format_candidate_table can preserve the type ordering.
    if (verbose) {
      verbose_rows <- vector("list", n_types)
      names(verbose_rows) <- names(degree_lookup)
    }
    
    # Build adjustment matrix once per cont_var, reused across all three
    # degrees (linear / FP1 / FP2). The adjustment set changes only by
    # removing var_name from selected_vars, so it is identical for all
    # three flex_fit calls for this variable.
    if (!skip_adjustment) {
      adj_vars     <- setdiff(selected_vars, var_name)
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
          zero               = processed_data$updated_params$zero_vars[adj_vars],
          catzero            = processed_data$updated_params$catzero_vars[adj_vars],
          spike              = processed_data$updated_params$spike_vars[adj_vars],
          spike_decision     = adj_spike_decision[adj_vars],
          keep_x_order       = FALSE,
          acd_parameter_list = NULL,
          check_binary       = TRUE
        )$x_transformed
      }
    } else {
      xadj_current <- xadj
    }
    
    for (interaction_type in names(degree_lookup)) {
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
        digits        = digits,
        scale_var     = if (!is.null(scale)) unname(scale[var_name]) else 1,
        run_test      = TRUE,
        compute_fitted = TRUE
      )
      
      # Annotate metrics -------------------------------------------------------
      metrics           <- fit_result$test_results$evaluation_metrics
      metrics$variable  <- var_name
      metrics$type      <- interaction_type
      
      # Collect verbose row for this candidate (printed as table after inner loop)
      if (verbose) {
        verbose_rows[[interaction_type]] <- metrics
      }
      
      all_idx <- all_idx + 1L
      all_metrics_list[[all_idx]] <- metrics
      
      # Store candidate model and fitted functions regardless of significance
      candidate_models[[interaction_type]]  <- fit_result$test_results$interaction_model
      candidate_fitted[[interaction_type]]  <- fit_result$fitted_functions
      
      # Select best model according to criterion --------------------------------
      if (criterion == "pvalue") {
        pval      <- metrics$pvalue[1L]
        # Tiebreaker: when two candidates have identical p-values, prefer the
        # one with the larger BIC improvement (dBIC = BIC_main - BIC_int).
        # BIC penalises complexity more heavily than AIC, so it favours the
        # simpler functional form when both are equally significant.
        bic_delta <- metrics$BIC_main_minus_int[1L]
        if (!is.na(pval) && pval < p_interact) {
          best_bic <- if (!is.null(best_fit))
            best_fit$test_results$evaluation_metrics$BIC_main_minus_int[1L]
          else -Inf
          if (pval < best_score ||
              (pval == best_score && !is.na(bic_delta) && bic_delta > best_bic)) {
            best_fit    <- fit_result
            best_metric <- metrics
            best_type   <- interaction_type
            best_score  <- pval
          }
        }
      } else if (criterion == "aic") {
        delta <- metrics$AIC_main_minus_int[1L]
        if (!is.na(delta) && delta > min_improvement && delta > best_score) {
          best_fit    <- fit_result
          best_metric <- metrics
          best_type   <- interaction_type
          best_score  <- delta
        }
      } else if (criterion == "bic") {
        delta <- metrics$BIC_main_minus_int[1L]
        if (!is.na(delta) && delta > min_improvement && delta > best_score) {
          best_fit    <- fit_result
          best_metric <- metrics
          best_type   <- interaction_type
          best_score  <- delta
        }
      } else {
        stop(
          paste0("! `criterion` must be one of 'pvalue', 'aic', or 'bic'; got '",
                 criterion, "'."),
          call. = FALSE
        )
      }
    }  # end interaction_type loop
    
    # Print the candidate table for this variable
    if (verbose) {
      format_candidate_table(verbose_rows, criterion = criterion, digits = digits)
    }
    
    if (!is.null(best_fit)) {
      best_idx <- best_idx + 1L
      best_metrics_list[[best_idx]]         <- best_metric
      best_interaction_model[[var_name]]    <- best_fit$test_results$interaction_model
      best_fitted_functions[[var_name]]     <- best_fit$fitted_functions
      center_vals_list[[var_name]]          <- best_fit$center_vals
      
      if (verbose) {
        type_label <- toupper(gsub("fp", "FP", best_type))
        if (criterion == "pvalue") {
          cat(sprintf(
            "        >> SELECTED: %s  (p = %s)\n",
            type_label,
            formatC(best_metric$pvalue[1L], format = "g", digits = digits)
          ))
        } else {
          ic_col <- if (criterion == "aic") "AIC_main_minus_int" else "BIC_main_minus_int"
          cat(sprintf(
            "        >> SELECTED: %s  (d%s = %s)\n",
            type_label, toupper(criterion),
            formatC(best_metric[[ic_col]][1L], format = "f", digits = 2)
          ))
        }
      }
      
      # Print regression table for the winning interaction model when requested.
      # show_models = TRUE only has effect when verbose = TRUE; suppressing
      # verbose output means suppressing all output including model summaries.
      # When no variable is selected, nothing is printed (best_fit is NULL).
      if (verbose && show_models) {
        fit_obj <- best_fit$test_results$interaction_model$fit
        if (!is.null(fit_obj)) {
          type_label <- toupper(gsub("fp", "FP", best_type))
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
      if (verbose) {
        cat("        >> NOT SELECTED: no candidate met the threshold.\n")
      }
      if (!quiet) {
        message(sprintf(
          "i No significant interaction retained for '%s'.", var_name
        ))
      }
    }
    
    # Store all three candidates for every variable regardless of significance
    all_interaction_models[[var_name]] <- candidate_models
    all_fitted_functions[[var_name]]   <- candidate_fitted
  }  # end cont_vars loop
  
  # Trim pre-allocated lists to actual fill level
  best_metrics_list <- best_metrics_list[seq_len(best_idx)]
  all_metrics_list  <- all_metrics_list[seq_len(all_idx)]
  
  # Drop NULL slots from best_interaction_model and best_fitted_functions
  best_interaction_model <- Filter(Negate(is.null), best_interaction_model)
  best_fitted_functions       <- Filter(Negate(is.null), best_fitted_functions)
  
  # Combine results -----------------------------------------------------------
  if (length(best_metrics_list) > 0L) {
    combined_metrics <- dplyr::bind_rows(best_metrics_list)
    # Move type to first column using base R subsetting - avoids rlang dependency
    combined_metrics <- combined_metrics[, c("type", setdiff(names(combined_metrics), "type")),
                                         drop = FALSE]
  } else {
    combined_metrics <- data.frame()
  }
  
  all_model_metrics <- dplyr::bind_rows(all_metrics_list)
  
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
    best_fitted_functions    = best_fitted_functions,
    all_fitted_functions     = all_fitted_functions,
    center_vals_list         = center_vals_list,
    group_var                = group_var,
    show_models              = show_models,
    group_levels_new         = cat_info$new_levels,
    group_levels_original    = cat_info$original_levels,
    flex                     = flex,
    criterion                = criterion,
    family                   = family,
    nobs                     = nobs
  )
}