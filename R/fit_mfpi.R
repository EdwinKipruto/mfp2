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
#' have been validated, expanded to full-length vectors, shifted, and scaled by
#' the calling \code{mfpi()} function.
#'
#' @param x Numeric matrix (\eqn{n \times p}) of predictors including
#'   `group_var`. Shift and scale have already been applied; no intercept
#'   column.
#' @param y Response vector or [survival::Surv()] object. For
#'   `family = "binomial"`, a numeric vector with exactly two distinct values.
#'   For `family = "cox"`, a two-column right-censored Surv object.
#' @param family Character string; one of `"gaussian"`, `"binomial"`,
#'   `"poisson"`, or `"cox"`.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param center Named logical vector of length \eqn{p}. Whether to centre
#'   each predictor before fitting.
#' @param criterion Character string; `"pvalue"`, `"aic"`, or `"bic"`.
#'   Governs adjustment-variable selection in Step 2.
#' @param select Named numeric vector of length \eqn{p}. Nominal significance
#'   levels for backward elimination of each predictor.
#' @param alpha Named numeric vector of length \eqn{p}. Significance levels
#'   for choosing between FP degrees.
#' @param force_keep Character vector of variable names forced into the
#'   adjustment model regardless of selection.
#' @param df Named numeric vector of length \eqn{p}. Degrees of freedom per
#'   predictor (1 = linear, 2m = FP of degree m).
#' @param xorder Character string; order of covariate entry into MFP
#'   backfitting — `"ascending"`, `"descending"`, or `"original"`.
#' @param fp_powers Named list of candidate FP power sets, one per predictor.
#' @param ties Character string; tie-handling method for Cox models —
#'   `"breslow"`, `"efron"`, or `"exact"`.
#' @param strata Integer vector of stratum memberships for stratified Cox
#'   models, or `NULL`.
#' @param nocenter Numeric vector passed to [survival::coxph()] to suppress
#'   centring for specific predictors. Cox models only.
#' @param acd_vars Named logical vector of length \eqn{p}. Whether each
#'   predictor undergoes the approximate cumulative distribution (ACD)
#'   transformation. Must be `FALSE` for all variables in `cont_vars`.
#' @param zero_vars Named logical vector of length \eqn{p}. Whether each
#'   predictor should treat non-positive values as zero before transformation.
#' @param use_ftest Logical. Use F-test rather than chi-square test for
#'   Gaussian models. Ignored for other families.
#' @param control List of fitting control parameters from
#'   [stats::glm.control()] or [survival::coxph.control()].
#' @param group_var Character string. Name of the categorical grouping variable
#'   in `x`.
#' @param include_group_var Logical. Whether `group_var` is included as a
#'   covariate in the MFP adjustment model.
#' @param flex Character string; `"flex1"`, `"flex2"`, `"flex3"`, or
#'   `"flex4"`. Controls how FP powers are estimated and constrained across
#'   groups.
#' @param cont_vars Character vector of continuous variables to test for
#'   interaction with `group_var`.
#' @param p_interact Numeric. Nominal significance level for the interaction
#'   test when `criterion = "pvalue"`.
#' @param min_improvement Numeric. Minimum required improvement in the
#'   selection criterion for retaining an interaction term.
#' @param show_models Logical. Whether to print model summaries for the final
#'   interaction models.
#' @param verbose Logical. Whether to print progress messages.
#' @param digits Positive integer. Significant digits for printed output.
#'
#' @return A list with two components:
#' \describe{
#'   \item{`adjustment_model`}{Fitted MFP adjustment model from [mfp2::mfp2()],
#'     including selected variables and their FP terms.}
#'   \item{`univariable_interactions`}{List of results from the univariable
#'     interaction evaluation; see \code{evaluate_interactions()} for structure.}
#' }
#'
#' @section Algorithm:
#' **Pre-processing (internal).** The levels of `group_var` are remapped to
#' consecutive integers starting at 0. Dummy variables are appended to `x`
#' when `include_group_var = TRUE`. All modelling-parameter vectors are
#' synchronised with the modified column set. This step is silent.
#'
#' **Step 1 — Adjustment model.** \code{mfp2:::fit_mfp()} selects adjustment
#' variables and their FP transformations. If `criterion = "pvalue"`, variables
#' may be dropped by backward elimination; if `alpha > 0` and `df > 1`,
#' nonlinear transforms may be estimated. Selected variables and their FP
#' powers are fixed for the remainder of the algorithm.
#'
#' **Step 2 — Univariable interactions.** For each variable in `cont_vars`,
#' linear, FP1, and FP2 interaction models are fitted and compared to their
#' respective main-effects models. The best functional form is selected
#' according to `criterion` and `min_improvement`; variables whose best model
#' does not meet the threshold are dropped.
#'
#' @seealso \code{mfpi()}, \code{preprocess_data()}, \code{fit_adjustment_model()},
#'   \code{evaluate_interactions()}
#'
#' @keywords internal
#' @noRd
fit_mfpi <- function(x, y, family, family_string, weights, offset, cycles,
                     center, criterion, select, alpha, force_keep, df, xorder,
                     fp_powers, ties, strata, nocenter, acd_vars, zero_vars,
                     catzero_vars, spike_vars, min_prop, max_prop, use_ftest,
                     control, group_var, include_group_var, flex, cont_vars,
                     p_interact, min_improvement, show_models, verbose,
                     digits) {
  
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
    force_keep        = force_keep
  )
  
  # ---------------------------------------------------------------------------
  # Step 1: Fit adjustment model via MFP
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("\ni === Fitting Adjustment Model via MFP ===\n")
    cat("i Step 1: Adjustment model\n")
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
    verbose        = FALSE
  )
  
  if (verbose) {
    print(adjustment_model$fp_terms)
  }
  
  # Identify variables retained by the adjustment model
  selected_vars <- get_selected_variables(adjustment_model)
  if (include_group_var) {
    selected_vars <- setdiff(selected_vars, processed_data$dummy_names)
  }
  
  adj_fp_powers <- get_fp_powers(
    rownames(adjustment_model$fp_terms),
    adjustment_model$fp_terms
  )
  
  # ---------------------------------------------------------------------------
  # Step 2: Evaluate univariable interactions
  # ---------------------------------------------------------------------------
  univ_results <- evaluate_interactions(
    y                  = y,
    processed_data     = processed_data,
    selected_vars      = selected_vars,
    adj_fp_powers      = adj_fp_powers,
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
    min_improvement    = min_improvement
  )
  
  n_significant <- nrow(univ_results$model_evaluation_metrics)
  if (n_significant == 0L) {
    warning("! No significant interactions were found.", call. = FALSE)
  }
  
  list(
    # Top-level fields: directly accessible as fit$field by print/summary methods
    model_evaluation_metrics = univ_results$model_evaluation_metrics,
    all_model_metrics        = univ_results$all_model_metrics,
    interaction_models       = univ_results$interaction_models,
    fitted_functions         = univ_results$fitted_functions,
    adjust_terms             = adjustment_model$fp_terms,
    group_var                = univ_results$group_var,
    show_models              = univ_results$show_models,
    flex                     = univ_results$flex,
    group_levels_new         = univ_results$group_levels_new,
    group_levels_original    = univ_results$group_levels_original,
    family                   = univ_results$family,
    nobs                     = univ_results$nobs,
    # Full sub-objects retained for programmatic access
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
                            force_keep) {
  
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
    force_keep   = force_keep
  )
  
  dummy_names <- NULL
  
  if (include_group_var) {
    # Create dummy variables (reference category dropped automatically)
    catvar_dummies <- create_dummies(group_values)
    dummy_names    <- colnames(catvar_dummies)
    x_filtered     <- cbind(x_filtered, catvar_dummies)
    n_dummies      <- length(dummy_names)
    
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
      group_var       = group_var
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
#' A thin wrapper around `mfp2:::fit_mfp()` that selects adjustment variables
#' and their FP transformations. Scale and shift are set to 1 and 0
#' respectively because \code{mfpi()} has already applied them to `x`.
#'
#' @param x Numeric matrix of adjustment predictors (group_var already
#'   removed).
#' @param y Response vector or Surv object.
#' @param weights Numeric vector of observation weights.
#' @param offset Numeric vector of linear-predictor offsets.
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param family GLM family object (e.g. \code{gaussian()}) or the character
#'   string \code{"cox"}. Passed directly to \code{mfp2:::fit_mfp()} as its
#'   \code{family} argument.
#' @param family_string Character string of the family name (e.g.
#'   \code{"gaussian"}, \code{"cox"}). Required separately by
#'   \code{mfp2:::fit_mfp()} for internal branching.
#' @param criterion Character string; \code{"pvalue"}, \code{"aic"}, or
#'   \code{"bic"}.
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
#' @param use_ftest Logical. Use F-test for Gaussian models.
#' @param control Fitting control list.
#' @param verbose Logical. Passed directly to \code{mfp2:::fit_mfp()}.
#'
#' @return The fitted model object returned by \code{mfp2:::fit_mfp()}, which
#'   includes \code{fp_terms} (a data frame of selected variables and their FP
#'   powers) and the underlying model fit.
#'
#' @details
#' This function calls \code{mfp2:::fit_mfp()} directly via \code{:::} because
#' \code{fit_mfp()} is not currently exported by the \pkg{mfp2} package. Once
#' \pkg{mfp2} exports \code{fit_mfp()}, this call should be updated to use
#' \code{::}. The \code{scale} and \code{shift} arguments are set to 1 and 0
#' respectively because \code{mfpi()} has already applied them to \code{x}.
#' @keywords internal
#' @noRd
fit_adjustment_model <- function(x, y, weights, offset, cycles, family,
                                 family_string, criterion, updated_params,
                                 xorder, ties, strata, nocenter, min_prop,
                                 max_prop, use_ftest, control,
                                 verbose = FALSE) {
  n_vars <- ncol(x)
  
  mfp2:::fit_mfp(
    x             = x,
    y             = y,
    weights       = weights,
    offset        = offset,
    cycles        = cycles,
    method        = ties,
    strata        = strata,
    nocenter      = nocenter,
    scale         = rep(1, n_vars),   # already applied by mfpi()
    shift         = rep(0, n_vars),   # already applied by mfpi()
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
    verbose       = verbose
  )
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
#'   variables, as returned by `get_fp_powers()`.
#' @param cont_vars Character vector of continuous variables to test.
#' @param flex Character string; `"flex1"`, `"flex2"`, `"flex3"`, or
#'   `"flex4"`.
#' @param weights Numeric vector of observation weights.
#' @param offset Numeric vector of linear-predictor offsets.
#' @param xorder Character string; entry order for MFP backfitting.
#' @param ties Character string; tie-handling for Cox models.
#' @param strata Integer stratum vector, or `NULL`.
#' @param use_ftest Logical. Use F-test for Gaussian models.
#' @param control Fitting control list.
#' @param nocenter Numeric vector for Cox centring suppression.
#' @param family Character string; regression family.
#' @param fp_powers Named list of candidate FP power sets for all predictors.
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param criterion Character string; `"pvalue"`, `"aic"`, or `"bic"`.
#' @param digits Positive integer. Significant digits for printed output.
#' @param group_var Character string. Name of the grouping variable.
#' @param show_models Logical. Whether to print model summaries.
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
#'   \item{`model_evaluation_metrics`}{Data frame of evaluation metrics for the
#'     best interaction model per variable. Empty if no interactions are
#'     significant.}
#'   \item{`all_model_metrics`}{Data frame of metrics for all models tested
#'     (linear, FP1, FP2) for every variable in `cont_vars`, regardless of
#'     significance.}
#'   \item{`interaction_models`}{Named list of fitted model objects for the best
#'     interaction type per significant variable.}
#'   \item{`fitted_functions`}{Named list of fitted FP function values per
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
#'     **smallest p-value**. Ties in p-value (e.g. both rounded to the same
#'     value) are broken by selecting the candidate with the **largest
#'     \code{AIC_diff}**, favouring the more parsimonious fit.}
#'   \item{\code{"aic"}}{Retain candidates with \code{AIC_diff} strictly above
#'     \code{min_improvement}. Among retained candidates, select the one with
#'     the **largest \code{AIC_diff}}.}
#'   \item{\code{"bic"}}{Same as \code{"aic"} using \code{BIC_diff}.}
#' }
#'
#' \code{AIC_diff} and \code{BIC_diff} are defined as
#' \eqn{\mathrm{AIC}_\text{main} - \mathrm{AIC}_\text{int}} and
#' \eqn{\mathrm{BIC}_\text{main} - \mathrm{BIC}_\text{int}} respectively.
#' A **positive** value indicates that the interaction model fits better
#' (lower information criterion) than the main-effects model.
#'
#' If no candidate clears its threshold, no interaction is retained for that
#' variable and it is excluded from \code{model_evaluation_metrics}.
#'
#' @keywords internal
#' @importFrom dplyr bind_rows relocate
#' @importFrom rlang .data
#' @noRd
evaluate_interactions <- function(y, processed_data, selected_vars,
                                  adj_fp_powers, cont_vars, flex,
                                  weights, offset, xorder, ties, strata,
                                  use_ftest, control, nocenter, family,
                                  family_string, fp_powers, cycles, criterion,
                                  digits, group_var, show_models, p_interact,
                                  min_improvement,
                                  xadj = NULL, skip_adjustment = FALSE,
                                  quiet = FALSE) {
  
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
  
  # Result containers
  best_metrics_list <- list()
  all_metrics_list  <- list()
  interaction_models <- list()
  fitted_functions   <- list()
  
  for (var_name in cont_vars) {
    
    best_fit   <- NULL
    best_metric <- NULL
    best_type  <- NULL
    # Initialise best_score so any valid model will beat it:
    #   pvalue criterion: seek the smallest value -> initialise to +Inf
    #   aic/bic criterion: seek the largest improvement -> initialise to -Inf
    best_score <- if (criterion == "pvalue") Inf else -Inf
    
    for (interaction_type in names(degree_lookup)) {
      degree <- degree_lookup[[interaction_type]]
      
      # Build or reuse the adjustment matrix ---------------------------------
      if (!skip_adjustment) {
        # Exclude the variable currently being tested from adjustment
        adj_vars <- setdiff(selected_vars, var_name)
        xadj_current <- NULL
        if (length(adj_vars) > 0L) {
          xadj_current <- mfp2::transform_matrix(
            x             = x[, adj_vars, drop = FALSE],
            power_list    = adj_fp_powers[adj_vars],
            center        = processed_data$updated_params$center[adj_vars],
            acdx          = processed_data$updated_params$acd_vars[adj_vars],
            zero          = processed_data$updated_params$zero_vars[adj_vars],
            catzero       = setNames(rep(FALSE, length(adj_vars)), adj_vars),
            keep_x_order  = FALSE,
            acd_parameter_list = NULL,
            check_binary  = TRUE
          )$x_transformed
        }
      } else {
        xadj_current <- xadj
      }
      
      # Fit and test the interaction model -----------------------------------
      fit_result <- flex_fit(
        x             = x,
        y             = y,
        cont_var      = var_name,
        group_var     = cat_info$group_var,
        xadj          = xadj_current,
        criterion     = criterion,
        ties          = ties,
        degree        = degree,
        family        = family,
        family_string = family_string,
        fp_cand       = processed_data$updated_params$fp_powers[[var_name]],
        use_ftest     = use_ftest,
        center        = processed_data$updated_params$center[var_name],
        xorder        = xorder,
        weights       = weights,
        offset        = offset,
        strata        = strata,
        control       = control,
        nocenter      = nocenter,
        cycles        = cycles,
        zero_var      = processed_data$updated_params$zero_vars[var_name],
        spike_var     = FALSE,
        flex          = flex,
        digits        = digits
      )
      
      # Annotate metrics with variable and interaction type ------------------
      metrics           <- fit_result$test_results$evaluation_metrics
      metrics$variable  <- var_name
      metrics$type      <- interaction_type
      
      all_metrics_list[[paste(var_name, interaction_type, sep = "_")]] <- metrics
      
      # Select best model according to criterion ----------------------------
      if (criterion == "pvalue") {
        pval      <- metrics$pvalue[1L]
        aic_delta <- metrics$AIC_diff[1L]
        if (!is.na(pval) && pval < p_interact) {
          # Primary: smallest p-value. Tie-break: largest AIC_diff.
          best_aic <- if (!is.null(best_fit))
            best_fit$test_results$evaluation_metrics$AIC_diff[1L]
          else -Inf
          if (pval < best_score ||
              (pval == best_score && !is.na(aic_delta) && aic_delta > best_aic)) {
            best_fit    <- fit_result
            best_metric <- metrics
            best_type   <- interaction_type
            best_score  <- pval
          }
        }
      } else if (criterion == "aic") {
        delta <- metrics$AIC_diff[1L]
        if (!is.na(delta) && delta > min_improvement && delta > best_score) {
          best_fit    <- fit_result
          best_metric <- metrics
          best_type   <- interaction_type
          best_score  <- delta
        }
      } else if (criterion == "bic") {
        delta <- metrics$BIC_diff[1L]
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
    
    if (!is.null(best_fit)) {
      best_metrics_list[[length(best_metrics_list) + 1L]] <- best_metric
      interaction_models[[var_name]] <- best_fit$test_results$interaction_model
      fitted_functions[[var_name]]   <- best_fit$fitted_functions
    } else if (!quiet) {
      message(sprintf(
        "i No significant interaction retained for '%s'.", var_name
      ))
    }
  }  # end cont_vars loop
  
  # Combine results -----------------------------------------------------------
  if (length(best_metrics_list) > 0L) {
    combined_metrics <- dplyr::bind_rows(best_metrics_list)
    combined_metrics <- dplyr::relocate(combined_metrics, .data$type, .before = 1L)
  } else {
    combined_metrics <- data.frame()
  }
  
  all_model_metrics <- dplyr::bind_rows(all_metrics_list, .id = "var_model")
  
  list(
    model_evaluation_metrics = combined_metrics,
    all_model_metrics        = all_model_metrics,
    interaction_models       = interaction_models,
    fitted_functions         = fitted_functions,
    group_var                = group_var,
    show_models              = show_models,
    group_levels_new         = cat_info$new_levels,
    group_levels_original    = cat_info$original_levels,
    flex                     = flex,
    family                   = family,
    nobs                     = nobs
  )
}