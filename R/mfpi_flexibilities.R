# Flex strategy functions for MFPI
#
# These internal functions estimate FP powers, build main-effects and
# interaction design matrices, and test for interaction between a continuous
# variable (`cont_var`) and a grouping variable (`group_var`). None are
# exported; they are called exclusively via `flex_fit()`.
#
# Naming convention used throughout:
#   cont_var   - the continuous variable being tested for interaction
#   group_var  - the categorical grouping variable (e.g. treatment)
#   xadj       - pre-transformed, pre-centered adjustment matrix (or NULL)
#   ties       - Cox tie-handling method (was `method` in original)
#   use_ftest  - F-test flag for Gaussian models (was `ftest`)
#   zero_var   - logical; treat non-positive values of cont_var as zero
#   fp_cand    - candidate FP powers for cont_var (was `powers`)


# -----------------------------------------------------------------------------
# flex_fit() - dispatcher -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Dispatch to the Appropriate Flex Implementation
#'
#' Routes to `flex0`, `flex1`, `flex2`, `flex3`, or `flex4` depending on the
#' value of `flex` and `degree`. When `degree < 1` the call is always routed
#' to `flex0` (linear), regardless of the `flex` argument.
#'
#' This function performs light input validation (checking `flex`, `cont_var`,
#' and `group_var`) before dispatching; the individual flex functions assume
#' these checks have already been done.
#'
#' @param x Numeric matrix (\eqn{n \times p}) of all predictors, including
#'   `group_var`. Shift and scale have already been applied. No intercept
#'   column.
#' @param y Response vector or [survival::Surv()] object.
#' @param cont_var Character string. Name of the continuous variable in `x` to
#'   test for interaction with `group_var`.
#' @param group_var Character string. Name of the categorical grouping variable
#'   in `x`. The lowest unique value is the reference group.
#' @param xadj Numeric matrix of pre-transformed, pre-centered adjustment
#'   predictors, or `NULL` if there are no adjustment variables.
#' @param criterion Character string; `"pvalue"`, `"aic"`, or `"bic"`.
#' @param ties Character string; Cox tie-handling - `"breslow"`, `"efron"`,
#'   or `"exact"`. Ignored for non-Cox families.
#' @param degree Integer; FP degree - `0` (linear), `1` (FP1), or `2` (FP2).
#'   Values below 1 are silently treated as `0` and routed to `flex0`.
#' @param family Character string; `"gaussian"`, `"binomial"`, `"poisson"`,
#'   or `"cox"`.
#' @param fp_cand Numeric vector of candidate FP powers for `cont_var`.
#'   Corresponds to one element of the `fp_powers` list in [mfp2::mfpi()]. When
#'   `degree == 1` (FP1), power = 1 is automatically excluded from `fp_cand`
#'   inside `flex_fit()` because linear is handled separately by `flex0`;
#'   keeping power = 1 in the FP1 candidate set would allow FP1 to collapse
#'   to a duplicate linear fit. Higher degrees retain power = 1 since
#'   combinations such as `(1, 2)` are distinct from linear.
#' @param use_ftest Logical. Use F-test rather than chi-square for Gaussian
#'   models. Currently applied only to adjustment-variable selection; not yet
#'   implemented for the interaction test itself.
#' @param center Logical scalar. Whether to centre `cont_var` before fitting.
#'   Adjustment variables are assumed already centred.
#' @param xorder Character string; entry order for MFP backfitting -
#'   `"ascending"`, `"descending"`, or `"original"`. Retained for API
#'   compatibility; has no effect when variable selection is disabled.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param strata Integer stratum vector for stratified Cox models, or `NULL`.
#' @param control Fitting control list from [stats::glm.control()] or
#'   [survival::coxph.control()].
#' @param nocenter Numeric vector for Cox centring suppression; see
#'   [survival::coxph()].
#' @param cycles Positive integer. Maximum MFP backfitting iterations.
#' @param zero_var Logical scalar. Whether non-positive values of `cont_var`
#'   should be treated as structural zeros before FP transformation.
#' @param spike_var Logical scalar. Whether `cont_var` is subject to
#'   spike-at-zero handling. Always `FALSE` for interaction testing; see
#'   \code{evaluate_interactions()}.
#' @param min_prop Numeric. Minimum proportion of zeros for SAZ modelling.
#' @param max_prop Numeric. Maximum proportion of zeros for SAZ modelling.
#' @param flex Character string; `"flex0"`, `"flex1"`, `"flex2"`, `"flex3"`,
#'   or `"flex4"`. Controls how FP powers are estimated and constrained across
#'   groups. See the *Flex levels* section in \code{mfpi()} for details.
#' @param digits Positive integer. Significant digits for printed output.
#' @param run_test Logical. Whether to perform the interaction test. Default
#'   `TRUE`.
#' @param compute_fitted Logical. Whether to compute group-specific fitted
#'   functions. Requires `run_test = TRUE`. Default `TRUE`.
#'
#' @note \code{force_max_fp} is constructed internally inside \code{flex1()}
#'   and \code{flex4()} and passed to \code{fit_mfp()}. It is a named
#'   logical vector of length \code{nvars} that is \code{TRUE} only for the
#'   \code{cont_var} columns, telling \code{select_ic()} to use the most
#'   complex functional form at the requested degree without AIC/BIC
#'   simplification. It is not a parameter of \code{flex_fit()} or the
#'   individual flex functions.
#'
#' @return A list with components:
#' \describe{
#'   \item{`bestfp_main`}{Numeric vector. Selected FP powers for the
#'     main-effects model.}
#'   \item{`bestfp_interaction`}{Named list of FP power vectors, one per group
#'     of `group_var`.}
#'   \item{`xmain`}{Design matrix for the main-effects model.}
#'   \item{`xinteraction`}{Design matrix for the interaction model.}
#'   \item{`znames`}{Character vector of interaction-term column names.}
#'   \item{`test_results`}{Interaction test output from `test_interaction()`,
#'     or `NULL` if `run_test = FALSE`.}
#'   \item{`fitted_functions`}{Group-specific fitted values from
#'     `gen_fitted_values_per_group()`, or `NULL` if `compute_fitted = FALSE`.}
#' }
#'
#' @keywords internal
#' @noRd
flex_fit <- function(x, y, cont_var, group_var, group_dummies, xadj,
                     criterion, ties, degree, family, family_string, fp_cand,
                     use_ftest, center, xorder, weights, offset, strata,
                     control, nocenter, cycles, zero_var, spike_var = FALSE,
                     min_prop = 0.05, max_prop = 0.95,
                     flex, digits,
                     run_test = TRUE, compute_fitted = TRUE) {
  
  # Input validation -----------------------------------------------------------
  if (!is.character(flex) || length(flex) != 1L) {
    stop("`flex` must be a single character string.", call. = FALSE)
  }
  valid_flex <- c("flex0", "flex1", "flex2", "flex3", "flex4")
  if (!flex %in% valid_flex) {
    stop(
      paste0("`flex` must be one of: ", paste(valid_flex, collapse = ", "), "."),
      call. = FALSE
    )
  }
  if (!is.character(cont_var) || length(cont_var) != 1L || !cont_var %in% colnames(x)) {
    stop("`cont_var` must be a single character string naming a column of `x`.",
         call. = FALSE)
  }
  if (!is.character(group_var) || length(group_var) != 1L || !group_var %in% colnames(x)) {
    stop("`group_var` must be a single character string naming a column of `x`.",
         call. = FALSE)
  }
  
  # Linear degree always uses flex0, regardless of the flex setting
  if (degree < 1L) flex <- "flex0"
  
  # For FP1 (degree = 1), exclude power = 1 from the candidate set so that
  # the FP1 search cannot collapse to a linear fit. Linear is already handled
  # separately by flex0 (degree = 0), so testing it again as FP1 with power 1
  # would produce a vacuous duplicate of the LINEAR candidate.
  # Higher degrees (FP2+) retain power 1 since combinations like (1, 2),
  # (1, 3) etc. are genuinely distinct from linear.
  if (degree == 1L) {
    fp_cand <- setdiff(fp_cand, 1)
    if (length(fp_cand) < 1L) {
      stop(
        "! `fp_cand` is empty for FP1 after excluding power = 1. ",
        "Provide a candidate set containing at least one non-unity power.",
        call. = FALSE
      )
    }
  }
  
  # Build shared argument list and dispatch ------------------------------------
  common_args <- list(
    x              = x,
    y              = y,
    cont_var       = cont_var,
    group_var      = group_var,
    group_dummies  = group_dummies,
    xadj           = xadj,
    criterion      = criterion,
    ties           = ties,
    degree         = degree,
    family         = family,
    family_string  = family_string,
    fp_cand        = fp_cand,
    use_ftest      = use_ftest,
    center         = center,
    xorder         = xorder,
    weights        = weights,
    offset         = offset,
    strata         = strata,
    control        = control,
    nocenter       = nocenter,
    cycles         = cycles,
    zero_var       = zero_var,
    spike_var      = spike_var,
    min_prop       = min_prop,
    max_prop       = max_prop,
    digits         = digits,
    run_test       = run_test,
    compute_fitted = compute_fitted
  )
  
  do.call(get(flex, mode = "function"), common_args)
}


# -----------------------------------------------------------------------------
# flex0() - linear interaction ------------------------------------------------
# -----------------------------------------------------------------------------

#' Linear Interaction Model (flex0)
#'
#' Fits main-effects and interaction models using `cont_var` as a plain linear
#' term (no FP transformation). This is the baseline case used whenever
#' `degree = 0`, or when `degree < 1` is passed to \code{flex_fit()}.
#'
#' @inheritParams flex_fit
#' @return See \code{flex_fit()} for the return structure. `bestfp_main` and all
#'   elements of `bestfp_interaction` are set to `1` (i.e. linear power).
#' @keywords internal
#' @noRd
flex0 <- function(x, y, cont_var, group_var, xadj, criterion, ties,
                  degree, family, family_string, fp_cand, use_ftest,
                  center, xorder, weights, offset, strata, control,
                  nocenter, cycles, zero_var, spike_var = FALSE,
                  min_prop = 0.05, max_prop = 0.95,
                  group_dummies, digits,
                  run_test = TRUE, compute_fitted = FALSE) {
  
  if (compute_fitted && !run_test) {
    stop("! `compute_fitted = TRUE` requires `run_test = TRUE`.", call. = FALSE)
  }
  if (!is.null(xadj)) {
    if (!is.matrix(xadj)) stop("`xadj` must be a matrix.", call. = FALSE)
    if (is.null(colnames(xadj))) stop("`xadj` must have column names.", call. = FALSE)
  }
  
  # Extract column vectors -----------------------------------------------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]
  group_levels <- sort(unique(groupvar_vec))
  k            <- length(group_levels)
  
  # Design matrices ------------------------------------------------------------
  
  z_vars <- create_z_variables(
    cont_var  = contvar_vec,
    group_var = groupvar_vec,
    power     = 1L,
    scale     = 1,
    shift     = 0,
    center    = center,
    zero      = zero_var
  )
  
  x_main        <- cbind(group_dummies, z_vars$xtransformed, xadj)
  x_interaction <- cbind(group_dummies, z_vars$z,            xadj)
  
  # Power bookkeeping ----------------------------------------------------------
  znames              <- sprintf("%s%d%d", cont_var, group_levels, 1L)
  bestfp_main         <- 1L
  bestfp_interaction  <- setNames(replicate(k, list(1L), simplify = TRUE), znames)
  
  # Interaction test -----------------------------------------------------------
  test_results <- NULL
  if (run_test) {
    test_results <- test_interaction(
      y                  = y,
      cont_var       = contvar_vec,
      group_var      = groupvar_vec,
      xmain              = x_main,
      xinteraction       = x_interaction,
      degree             = 0L,
      bestfp_main        = bestfp_main,
      bestfp_interaction = bestfp_interaction,
      flex               = "flex0",
      use_ftest          = use_ftest,
      family             = family,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      family_string      = family_string,
      digits             = digits
    )
  }
  
  # Fitted functions -----------------------------------------------------------
  fitted_functions <- NULL
  if (compute_fitted) {
    fitted_functions <- gen_fitted_values_per_group(
      cont_var       = contvar_vec,
      group_fp_powers   = bestfp_interaction,
      interaction_model = test_results$interaction_model,
      group_var      = groupvar_vec,
      family            = family,
      family_string     = family_string
    )
  }
  
  list(
    bestfp_main        = bestfp_main,
    bestfp_interaction = bestfp_interaction,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    test_results       = test_results,
    fitted_functions   = fitted_functions
  )
}


# -----------------------------------------------------------------------------
# flex1() - main-effect FP powers applied within groups ----------------------
# -----------------------------------------------------------------------------

#' FP Powers from the Pooled Main-Effects Model (flex1)
#'
#' Estimates FP powers from a pooled model that includes only the main effect
#' of `cont_var` (without interaction terms). The same powers are then used to
#' transform `cont_var` within each group for both the main-effects and
#' interaction models. This is the default and least flexible strategy.
#'
#' Using identical powers across groups reduces the risk of overfitting and
#' keeps the main-effects and interaction models nested, which validates the
#' likelihood-ratio interaction test (Royston and Sauerbrei 2004).
#'
#' @inheritParams flex_fit
#' @return See \code{flex_fit()} for the return structure.
#' @keywords internal
#' @noRd
flex1 <- function(x, y, cont_var, group_var, xadj, criterion, ties,
                  degree, family, family_string, fp_cand, use_ftest,
                  center, xorder, weights, offset, strata, control,
                  nocenter, cycles, zero_var, spike_var = FALSE,
                  min_prop = 0.05, max_prop = 0.95,
                  group_dummies, digits,
                  run_test = TRUE, compute_fitted = FALSE) {
  
  if (compute_fitted && !run_test) {
    stop("! `compute_fitted = TRUE` requires `run_test = TRUE`.", call. = FALSE)
  }
  if (degree < 1L) {
    stop("! `degree` must be >= 1 for flex1; use flex0 for linear models.",
         call. = FALSE)
  }
  
  # Extract column vectors -----------------------------------------------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]
  group_levels <- sort(unique(groupvar_vec))
  k            <- length(group_levels)
  
  # Step 1: Fit pooled main-effects MFP model to find best FP powers ----------
  # All variables are forced in (select = alpha = 1, keep = vnames).
  # force_max_fp: named logical vector, TRUE for cont_var only.
  # Tells select_ic() to select the most complex functional form at the degree
  # fixed by df_vec, without competing against simpler forms under AIC/BIC.
  # The best power combination within that degree is still found by
  # find_best_fpm_step() via deviance minimisation (equivalent to AIC/BIC
  # at fixed df). Functional form selection (linear/FP1/FP2) happens in
  # evaluate_interactions() via the user criterion.
  xnew    <- cbind(contvar_vec, group_dummies, xadj)
  vnames  <- colnames(xnew)
  n_adj   <- if (is.null(xadj)) 0L else ncol(xadj)
  
  df_vec <- c(
    2L * degree,               # cont_var: FP of requested degree
    rep(1L, k - 1L),           # group dummies: linear
    if (n_adj > 0L) rep(1L, n_adj)
  )
  
  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  pw_list <- setNames(replicate(length(vnames), default_powers, simplify = FALSE),
                      vnames)
  pw_list[[cont_var]] <- fp_cand
  
  n_total      <- length(vnames)
  center_vec   <- c(center,    rep(FALSE, n_total - 1L))
  zero_vec     <- c(zero_var,  rep(FALSE, n_total - 1L))
  spike_vec    <- c(spike_var, rep(FALSE, n_total - 1L))
  acd_vec      <- setNames(rep(FALSE, n_total), vnames)
  catzero      <- setNames(rep(FALSE, n_total), vnames)
  # force_max_fp: TRUE only for cont_var; dummies and adj have df = 1
  # so there is no simpler form to downgrade to.
  force_max_fp <- setNames(c(TRUE, rep(FALSE, n_total - 1L)), vnames)
  
  fit_main_pool <- fit_mfp(
    x             = xnew,
    y             = y,
    df            = df_vec,
    criterion     = criterion,
    select        = rep(1, n_total),
    alpha         = rep(1, n_total),
    keep          = vnames,
    force_max_fp  = force_max_fp,
    method        = ties,
    family        = family,
    family_string = family_string,
    powers        = pw_list,
    ftest         = use_ftest,
    center        = center_vec,
    shift         = rep(0, n_total),
    scale         = rep(1, n_total),
    acdx          = acd_vec,
    xorder        = xorder,
    weights       = weights,
    offset        = offset,
    strata        = strata,
    control       = control,
    nocenter      = nocenter,
    cycles        = cycles,
    catzero       = catzero,
    zero          = zero_vec,
    spike         = spike_vec,
    min_prop      = min_prop,
    max_prop      = max_prop,
    verbose       = FALSE
  )
  
  # Step 2: Extract best FP powers for cont_var --------------------------------
  bestfp <- unlist(get_fp_powers(cont_var, fit_main_pool$fp_terms))
  
  # Defensive guard: if mfp2 returned an unexpected linear or empty result for
  # an FP1 search (despite force_max_fp = TRUE and power = 1 excluded from
  # fp_cand by flex_fit()), perform a direct deviance-minimisation search over
  # the non-unity candidates. This ensures FP1 is always a non-linear functional
  # form, since linear is handled separately by flex0.
  needs_fallback <- degree == 1L && (
    length(bestfp) == 0L ||
      (length(bestfp) == 1L && isTRUE(bestfp == 1))
  )
  if (needs_fallback) {
    nonunit_cand <- setdiff(fp_cand, 1)
    if (length(nonunit_cand) == 0L) {
      stop(
        "! FP1 fallback failed: `fp_cand` contains no non-unity powers.",
        call. = FALSE
      )
    }
    # Fit a one-term FP transformation at each candidate power and choose
    # the lowest-deviance fit. Adjustment columns and group dummies are
    # included as covariates so the comparison is on equal footing with the
    # MFP-fitted model.
    fixed_x <- cbind(group_dummies, xadj)
    dev_vec <- vapply(nonunit_cand, function(p) {
      z_p <- transform_vector_fp(
        x     = contvar_vec,
        power = p,
        shift = 0,
        scale = 1,
        zero  = zero_var,
        name  = cont_var
      )
      fit_p <- fit_model(
        x        = cbind(z_p, fixed_x),
        y        = y,
        family   = family,
        weights  = weights,
        offset   = offset,
        method   = ties,
        strata   = strata,
        control  = control,
        nocenter = nocenter,
        rownames = NULL,
        fast     = TRUE
      )
      -2 * fit_p$logl
    }, numeric(1L))
    bestfp <- nonunit_cand[which.min(dev_vec)]
  }
  
  # Step 3: Build interaction terms using the pooled FP powers ----------------
  znames             <- sprintf("%s%d1", cont_var, group_levels)
  bestfp_interaction <- setNames(
    replicate(length(znames), list(bestfp), simplify = TRUE), znames
  )
  
  z_vars <- create_z_variables(
    cont_var  = contvar_vec,
    group_var = groupvar_vec,
    power     = bestfp,
    shift     = 0,
    scale     = 1,
    center    = center,
    zero      = zero_var
  )
  
  x_main        <- cbind(group_dummies, z_vars$xtransformed, xadj)
  x_interaction <- cbind(group_dummies, z_vars$z,            xadj)
  
  # Step 4: Interaction test ---------------------------------------------------
  test_results <- NULL
  if (run_test) {
    test_results <- test_interaction(
      y                  = y,
      cont_var       = contvar_vec,
      group_var      = groupvar_vec,
      xmain              = x_main,
      xinteraction       = x_interaction,
      degree             = degree,
      bestfp_main        = bestfp,
      bestfp_interaction = bestfp_interaction,
      flex               = "flex1",
      use_ftest          = use_ftest,
      family             = family,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      family_string      = family_string,
      digits             = digits
    )
  }
  
  # Step 5: Fitted functions ---------------------------------------------------
  fitted_functions <- NULL
  if (compute_fitted) {
    fitted_functions <- gen_fitted_values_per_group(
      cont_var       = contvar_vec,
      group_fp_powers   = bestfp_interaction,
      interaction_model = test_results$interaction_model,
      group_var      = groupvar_vec,
      family            = family,
      family_string     = family_string
    )
  }
  
  list(
    bestfp_main        = bestfp,
    bestfp_interaction = bestfp_interaction,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    test_results       = test_results,
    fitted_functions   = fitted_functions
  )
}


# -----------------------------------------------------------------------------
# flex2() - within-group FP powers, constrained equal across groups -----------
# -----------------------------------------------------------------------------

#' FP Powers Estimated Within Groups, Constrained Equal (flex2)
#'
#' Estimates FP powers by fitting `cont_var` within each group simultaneously,
#' subject to the constraint that the powers are the same across groups. The
#' best power set (lowest deviance) is then used for both the main-effects and
#' interaction design matrices. The degrees of freedom for the interaction test
#' are the same as for `flex1`, but because power estimation consumes additional
#' degrees of freedom that are not accounted for, p-values are slightly
#' anti-conservative and should be treated as indicative.
#'
#' @inheritParams flex_fit
#' @return See \code{flex_fit()} for the return structure.
#' @keywords internal
#' @noRd
flex2 <- function(x, y, cont_var, group_var, xadj, criterion, ties,
                  degree, family, family_string, fp_cand, use_ftest,
                  center, xorder, weights, offset, strata, control,
                  nocenter, cycles, zero_var, spike_var = FALSE,
                  min_prop = 0.05, max_prop = 0.95,
                  group_dummies, digits,
                  run_test = TRUE, compute_fitted = FALSE) {
  
  if (compute_fitted && !run_test) {
    stop("! `compute_fitted = TRUE` requires `run_test = TRUE`.", call. = FALSE)
  }
  
  # Extract column vectors -----------------------------------------------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]
  
  # Step 1: Generate all candidate within-group transformed variables ----------
  # transform_z_variables() tries every combination of powers from fp_cand
  # and returns a named list of matrices plus a power matrix.
  transformed <- transform_z_variables(
    cont_var  = contvar_vec,
    group_var = groupvar_vec,
    shift     = 0,
    scale     = 1,
    fp_cand   = fp_cand,
    fp_degree = degree,
    acdx      = FALSE,
    center    = FALSE    # centering deferred to Step 3
  )
  
  power_matrix      <- transformed$powers_matrix
  transformed_groups <- transformed$z_transformed
  znames            <- transformed$znames
  
  # Step 2: Select the power combination with the lowest deviance -------------
  # Deviance is the same whether we use AIC/BIC/pvalue because df is fixed
  # for a given degree, so it is sufficient to minimise -2*logLik.
  fixed_x <- cbind(group_dummies, xadj)
  
  deviance_vec <- vapply(transformed_groups, function(z) {
    z[is.na(z)] <- 0    # structural zeros
    fit <- fit_model(
      x        = cbind(z, fixed_x),
      y        = y,
      family   = family,
      weights  = weights,
      offset   = offset,
      method   = ties,
      strata   = strata,
      control  = control,
      nocenter = nocenter,
      rownames = NULL,
      fast     = TRUE
    )
    -2 * fit$logl
  }, numeric(1L))
  
  best_idx <- which.min(deviance_vec)
  bestfp <- power_matrix[best_idx, ]
  
  # Step 3: Transform cont_var using best powers; optionally centre -----------
  contvar_transformed <- mfp2::transform_vector_fp(
    x     = contvar_vec,
    power = bestfp,
    name  = cont_var,
    scale = 1,
    shift = 0
  )
  z_best <- transformed_groups[[best_idx]]
  
  if (center) {
    # contvar_transformed: use mean of positive values only when zero_var = TRUE
    # (spike at zero); grand mean otherwise. center_matrix() handles both cases.
    ct_names <- colnames(contvar_transformed)
    if (is.null(ct_names))
      ct_names <- paste0(cont_var, seq_len(ncol(contvar_transformed)))
    contvar_transformed <- center_matrix(
      mat    = contvar_transformed,
      centers = NULL,
      zero   = setNames(rep(zero_var, ncol(contvar_transformed)), ct_names)
    )
    
    # z_best: each column is non-zero only for one group (structural zeros
    # elsewhere). zero = TRUE tells center_matrix() to compute the mean only
    # over non-zero rows, leaving structural zeros at zero.
    z_best <- center_matrix(
      mat    = z_best,
      centers = NULL,
      zero   = setNames(rep(TRUE, ncol(z_best)), colnames(z_best))
    )
  }
  
  # Step 4: Assemble design matrices -------------------------------------------
  x_main        <- cbind(group_dummies, contvar_transformed, xadj)
  x_interaction <- cbind(group_dummies, z_best,              xadj)
  
  bestfp_interaction <- setNames(
    replicate(length(znames), list(bestfp), simplify = TRUE), znames
  )
  
  # Step 5: Interaction test ---------------------------------------------------
  test_results <- NULL
  if (run_test) {
    test_results <- test_interaction(
      y                  = y,
      cont_var       = contvar_vec,
      group_var      = groupvar_vec,
      xmain              = x_main,
      xinteraction       = x_interaction,
      degree             = degree,
      bestfp_main        = bestfp,
      bestfp_interaction = bestfp_interaction,
      flex               = "flex2",
      use_ftest          = use_ftest,
      family             = family,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      family_string      = family_string,
      digits             = digits
    )
  }
  
  # Step 6: Fitted functions ---------------------------------------------------
  fitted_functions <- NULL
  if (compute_fitted) {
    fitted_functions <- gen_fitted_values_per_group(
      cont_var       = contvar_vec,
      group_fp_powers   = bestfp_interaction,
      interaction_model = test_results$interaction_model,
      group_var      = groupvar_vec,
      family            = family,
      family_string     = family_string
    )
  }
  
  list(
    bestfp_main        = bestfp,
    bestfp_interaction = bestfp_interaction,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    test_results       = test_results,
    fitted_functions   = fitted_functions
  )
}


# -----------------------------------------------------------------------------
# flex3() - separate FP powers for main and interaction, equal across groups --
# -----------------------------------------------------------------------------

#' Separate FP Powers for Main-Effects and Interaction Models (flex3)
#'
#' Combines `flex1` and `flex2`: the main-effects model uses FP powers
#' estimated from the pooled data (as in `flex1`), while the interaction model
#' uses powers estimated within groups subject to the equal-powers constraint
#' (as in `flex2`). Because the two models may use different FP families they
#' are non-nested, so the likelihood-ratio p-value is indicative rather than
#' exact. Simulation evidence suggests this approach recovers within-group
#' functional forms more accurately than `flex1` when a true interaction is
#' present (Royston and Sauerbrei 2014).
#'
#' @inheritParams flex_fit
#' @return See \code{flex_fit()} for the return structure.
#' @keywords internal
#' @noRd
flex3 <- function(x, y, cont_var, group_var, xadj, criterion, ties,
                  degree, family, family_string, fp_cand, use_ftest,
                  center, xorder, weights, offset, strata, control,
                  nocenter, cycles, zero_var, spike_var = FALSE,
                  min_prop = 0.05, max_prop = 0.95,
                  group_dummies, digits,
                  run_test = TRUE, compute_fitted = FALSE) {
  
  if (compute_fitted && !run_test) {
    stop("! `compute_fitted = TRUE` requires `run_test = TRUE`.", call. = FALSE)
  }
  
  # Shared arguments for flex1 and flex2 calls ---------------------------------
  shared <- list(
    x             = x, y = y, cont_var = cont_var, group_var = group_var,
    group_dummies = group_dummies,
    xadj          = xadj, criterion = criterion, ties = ties, degree = degree,
    family        = family, family_string = family_string, fp_cand = fp_cand,
    use_ftest     = use_ftest, center = center, xorder = xorder,
    weights       = weights, offset = offset, strata = strata,
    control       = control, nocenter = nocenter, cycles = cycles,
    zero_var      = zero_var, spike_var = spike_var,
    min_prop      = min_prop, max_prop = max_prop,
    digits        = digits,
    run_test      = FALSE, compute_fitted = FALSE
  )
  
  # Step 1: Main-effects design matrix from flex1 (pooled FP powers) ----------
  fit_flex1   <- do.call(flex1, shared)
  x_main      <- fit_flex1$xmain
  bestfp_main <- fit_flex1$bestfp_main
  
  # Step 2: Interaction design matrix from flex2 (within-group FP powers) -----
  fit_flex2          <- do.call(flex2, shared)
  x_interaction      <- fit_flex2$xinteraction
  bestfp_interaction <- fit_flex2$bestfp_interaction
  znames             <- fit_flex2$znames
  
  # Step 3: Interaction test using the two separate design matrices ------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]
  
  test_results <- NULL
  if (run_test) {
    test_results <- test_interaction(
      y                  = y,
      cont_var       = contvar_vec,
      group_var      = groupvar_vec,
      xmain              = x_main,
      xinteraction       = x_interaction,
      degree             = degree,
      bestfp_main        = bestfp_main,
      bestfp_interaction = bestfp_interaction,
      flex               = "flex3",
      use_ftest          = use_ftest,
      family             = family,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      family_string      = family_string,
      digits             = digits
    )
  }
  
  # Step 4: Fitted functions ---------------------------------------------------
  fitted_functions <- NULL
  if (compute_fitted) {
    fitted_functions <- gen_fitted_values_per_group(
      cont_var       = contvar_vec,
      group_fp_powers   = bestfp_interaction,
      interaction_model = test_results$interaction_model,
      group_var      = groupvar_vec,
      family            = family,
      family_string     = family_string
    )
  }
  
  list(
    bestfp_main        = bestfp_main,
    bestfp_interaction = bestfp_interaction,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    test_results       = test_results,
    fitted_functions   = fitted_functions
  )
}


# -----------------------------------------------------------------------------
# flex4() - group-specific FP powers (most flexible) -------------------------
# -----------------------------------------------------------------------------

#' Group-Specific FP Powers for Main-Effects and Interaction Models (flex4)
#'
#' The most flexible strategy. FP powers for the interaction model are
#' estimated independently within each group, so the functional form of
#' `cont_var` can differ across groups. The main-effects model uses pooled
#' powers estimated by `flex1`. The additional degrees of freedom reduce power
#' to detect interaction relative to `flex1`-`flex3`.
#'
#' @inheritParams flex_fit
#' @return See \code{flex_fit()} for the return structure.
#' @keywords internal
#' @noRd
flex4 <- function(x, y, cont_var, group_var, xadj, criterion, ties,
                  degree, family, family_string, fp_cand, use_ftest,
                  center, xorder, weights, offset, strata, control,
                  nocenter, cycles, zero_var, spike_var = FALSE,
                  min_prop = 0.05, max_prop = 0.95,
                  group_dummies, digits,
                  run_test = TRUE, compute_fitted = FALSE) {
  
  if (compute_fitted && !run_test) {
    stop("! `compute_fitted = TRUE` requires `run_test = TRUE`.", call. = FALSE)
  }
  
  # Extract column vectors -----------------------------------------------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]
  k            <- length(unique(groupvar_vec))
  
  # Step 1: Build untransformed group-wise variables (linear placeholder) ------
  # create_z_variables() with power = 1 produces one column per group;
  # these will be replaced by the FP-transformed versions in Step 3.
  z_linear <- create_z_variables(
    cont_var  = contvar_vec,
    group_var = groupvar_vec,
    power     = 1L,
    shift     = 0,
    scale     = 1,
    center    = FALSE
  )$z
  znames <- colnames(z_linear)
  
  # Step 2: Fit group-specific FP model via MFP --------------------------------
  # Each group's column is allowed its own FP powers (no equality constraint).
  # Non-positive values within groups are handled via zero = TRUE.
  xnew   <- cbind(z_linear, group_dummies, xadj)
  vnames <- colnames(xnew)
  n_adj  <- if (is.null(xadj)) 0L else ncol(xadj)
  n_total <- length(vnames)
  
  df_vec <- c(
    rep(2L * degree, k),     # one FP term per group
    rep(1L, k - 1L),         # group dummies
    if (n_adj > 0L) rep(1L, n_adj)
  )
  
  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  pw_list <- setNames(replicate(n_total, default_powers, simplify = FALSE), vnames)
  pw_list[znames] <- replicate(k, fp_cand, simplify = FALSE)
  
  # Centre only the within-group cont_var columns; adjustment already centred
  center_vec  <- c(rep(center,    k), rep(FALSE, n_total - k))
  zero_vec    <- c(rep(TRUE,      k), rep(FALSE, n_total - k))
  spike_vec   <- c(rep(spike_var, k), rep(FALSE, n_total - k))
  acd_vec      <- setNames(rep(FALSE, n_total), vnames)
  catzero_vec  <- setNames(rep(FALSE, n_total), vnames)
  # force_max_fp: named logical vector, TRUE for the k per-group cont_var
  # columns, FALSE for group dummies and adjustment variables (df = 1,
  # nothing to force). Same rationale as flex1.
  force_max_fp <- setNames(c(rep(TRUE, k), rep(FALSE, n_total - k)), vnames)
  
  
  fit_int <- fit_mfp(
    x             = xnew,
    y             = y,
    df            = df_vec,
    criterion     = criterion,
    select        = rep(1, n_total),
    alpha         = rep(1, n_total),
    keep          = vnames,
    force_max_fp  = force_max_fp,
    method        = ties,
    family        = family,
    family_string = family_string,
    powers        = pw_list,
    shift         = rep(0, n_total),
    scale         = rep(1, n_total),
    ftest         = use_ftest,
    acdx          = acd_vec,
    center        = center_vec,
    xorder        = xorder,
    weights       = weights,
    offset        = offset,
    strata        = strata,
    control       = control,
    nocenter      = nocenter,
    cycles        = cycles,
    zero          = zero_vec,
    catzero       = catzero_vec,
    spike         = spike_vec,
    min_prop      = min_prop,
    max_prop      = max_prop,
    verbose       = FALSE
  )
  
  bestfp_interaction <- get_fp_powers(znames, fit_int$fp_terms)
  
  # Step 3: Apply group-specific FP powers to build the interaction matrix -----
  transformed <- transform_matrix(
    x          = z_linear,
    power_list = bestfp_interaction,
    acdx       = setNames(rep(FALSE, k), znames),
    center     = setNames(rep(FALSE, k), znames),   # centring done below
    zero       = setNames(rep(TRUE,  k), znames),
    catzero    = setNames(rep(FALSE, k), znames)
  )$x_transformed
  
  # Infinite values arise when zero is log- or negative-power transformed;
  # Non-finite values arise from structural zeros transformed by log/negative
  # powers. Mark as NA temporarily so they do not distort centering.
  transformed[!is.finite(transformed)] <- NA
  
  if (center) {
    # zero = TRUE: center_matrix computes means only over non-NA (non-zero)
    # rows per column, leaving structural zeros (NA) unaffected.
    transformed <- center_matrix(
      mat    = transformed,
      centers = NULL,
      zero   = setNames(rep(TRUE, ncol(transformed)), colnames(transformed))
    )
  }
  transformed[is.na(transformed)] <- 0    # restore structural zeros
  
  x_interaction <- cbind(group_dummies, transformed, xadj)
  
  # Step 4: Main-effects design matrix via flex1 (pooled FP powers) -----------
  fit_main <- flex1(
    x             = x, y = y, cont_var = cont_var, group_var = group_var,
    group_dummies = group_dummies,
    xadj          = xadj, criterion = criterion, ties = ties, degree = degree,
    family        = family, family_string = family_string,
    fp_cand       = fp_cand, use_ftest = use_ftest,
    center        = center, xorder = xorder, weights = weights, offset = offset,
    strata        = strata, control = control, nocenter = nocenter,
    cycles        = cycles, zero_var = zero_var, spike_var = spike_var,
    min_prop      = min_prop, max_prop = max_prop,
    digits        = digits, run_test = FALSE, compute_fitted = FALSE
  )
  
  x_main      <- fit_main$xmain
  bestfp_main <- fit_main$bestfp_main
  
  # Step 5: Interaction test ---------------------------------------------------
  test_results <- NULL
  if (run_test) {
    test_results <- test_interaction(
      y                  = y,
      cont_var       = contvar_vec,
      group_var      = groupvar_vec,
      xmain              = x_main,
      xinteraction       = x_interaction,
      degree             = degree,
      bestfp_main        = bestfp_main,
      bestfp_interaction = bestfp_interaction,
      flex               = "flex4",
      use_ftest          = use_ftest,
      family             = family,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      family_string      = family_string,
      digits             = digits
    )
  }
  
  # Step 6: Fitted functions ---------------------------------------------------
  fitted_functions <- NULL
  if (compute_fitted) {
    fitted_functions <- gen_fitted_values_per_group(
      cont_var       = contvar_vec,
      group_fp_powers   = bestfp_interaction,
      interaction_model = test_results$interaction_model,
      group_var      = groupvar_vec,
      family            = family,
      family_string     = family_string
    )
  }
  
  list(
    bestfp_main        = bestfp_main,
    bestfp_interaction = bestfp_interaction,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    test_results       = test_results,
    fitted_functions   = fitted_functions
  )
}