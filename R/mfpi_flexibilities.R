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
#'   `"negbin"`, or `"cox"`.
#' @param fp_cand Numeric vector of candidate FP powers for `cont_var`.
#'   Corresponds to one element of the `fp_powers` list in [mfp2::mfpi()]. When
#'   `degree == 1` (FP1), power = 1 is automatically excluded from `fp_cand`
#'   inside `flex_fit()` because linear is handled separately by `flex0`;
#'   keeping power = 1 in the FP1 candidate set would allow FP1 to collapse
#'   to a duplicate linear fit. Higher degrees retain power = 1 since
#'   combinations such as `(1, 2)` are distinct from linear.
#' @param use_ftest Logical. If \code{TRUE} and \code{family = "gaussian"},
#'   use an F-test rather than a chi-square likelihood-ratio test when
#'   computing p-values. Applied to both the adjustment-variable selection
#'   (via \code{fit_mfp()}) and the interaction test
#'   (via \code{calculate_f_test()}). Ignored for non-Gaussian
#'   families. Default \code{FALSE}.
#' @param center Logical scalar. Whether to centre `cont_var` before fitting.
#'   Adjustment variables are assumed already centred.
#' @param center_type Character string; \code{"grand"} (default) or
#'   \code{"group"}. Passed to \code{create_z_variables()} to control the
#'   centering strategy when \code{center = TRUE}. See
#'   \code{create_z_variables()} for details.
#' @param scale_var Numeric scalar. The scale factor for \code{cont_var},
#'   as computed and applied in \code{mfpi.default()}. Used to backscale
#'   \code{cont_var} (multiply by \code{scale_var}) inside
#'   \code{create_z_variables()} and \code{transform_z_variables()} before
#'   FP transformation, so that interaction model coefficients are on the
#'   \eqn{\phi(x + \text{shift})} scale matching the adjustment model and
#'   standalone \pkg{mfp2}. Default \code{1} (no backscaling).
#' @param shift_var Numeric scalar. The shift factor for \code{cont_var}.
#'   Passed to \code{gen_fitted_values_per_group()} so fitted-function x
#'   coordinates are returned/displayed on the original raw scale. Default
#'   \code{0}.
#' @param xorder Character string; entry order for MFP backfitting -
#'   `"ascending"`, `"descending"`, or `"original"`. Retained for API
#'   compatibility; has no effect when variable selection is disabled.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param strata Optional high-level Cox stratification object, or `NULL`.
#'   Passed through to Cox model fits without converting ordinary vector/factor
#'   strata to integer codes.
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
#' @param min_saz_component_prop Numeric in \eqn{(0, 0.5)}. Minimum required
#'   proportion in each component of a spike-at-zero covariate. Passed through
#'   to \code{fit_mfp()} when an internal flex model uses spike-at-zero handling.
#' @param flex Character string; `"flex0"`, `"flex1"`, `"flex2"`, `"flex3"`,
#'   or `"flex4"`. Controls how FP powers are estimated and constrained across
#'   groups. See the *Flex levels* section in \code{mfpi()} for details.
#' @param run_test Logical. Whether to perform the interaction test. Default
#'   `TRUE`.
#'
#' @note \code{force_max_fp} is constructed internally inside \code{flex1()}
#'   and \code{flex4()} and passed to \code{mfp2:::fit_mfp()}. It is a named
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
#' }
#'
#' @keywords internal
#' @noRd
flex_fit <- function(x, y, cont_var, group_var, group_dummies, xadj,
                     criterion, ties, degree, family, family_string, fp_cand,
                     use_ftest, center, xorder, weights, offset, strata,
                     control, nocenter, cycles, zero_var, spike_var = FALSE,
                     min_saz_component_prop = 0.10,
                     flex, scale_var = 1, shift_var = 0,
                     center_type = c("grand", "group"), has_offset,
                     run_test = TRUE, fitter = "base") {

  center_type <- match.arg(center_type)

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
    fitter         = fitter,
    fp_cand        = fp_cand,
    use_ftest      = use_ftest,
    center         = center,
    center_type    = center_type,
    scale_var      = scale_var,
    shift_var      = shift_var,
    xorder         = xorder,
    weights        = weights,
    offset         = offset,
    strata         = strata,
    control        = control,
    nocenter       = nocenter,
    cycles         = cycles,
    zero_var       = zero_var,
    spike_var      = spike_var,
    min_saz_component_prop = min_saz_component_prop,
    run_test       = run_test,
    has_offset     = has_offset
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
                  min_saz_component_prop = 0.10,
                  group_dummies,
                  center_type = c("grand", "group"),
                  scale_var = 1, shift_var = 0, has_offset,
                  run_test = TRUE, fitter = "base") {

  center_type <- match.arg(center_type)


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
    zero        = zero_var,
    center_type = center_type,
    scale_var   = scale_var
  )

  x_main        <- cbind(group_dummies, z_vars$xtransformed, xadj)
  x_interaction <- cbind(group_dummies, z_vars$z,            xadj)
  center_vals   <- z_vars$center_vals
  coefficient_groups <- z_vars$column_groups
  # Power bookkeeping ----------------------------------------------------------
  znames              <- sprintf("%s%d%d", cont_var, group_levels, 1L)
  bestfp_main         <- 1L
  bestfp_interaction  <- setNames(replicate(k, 1L, simplify = FALSE), znames)

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
      family_string      = family_string,
      fitter             = fitter,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      has_offset         = has_offset
    )
  }


  list(
    bestfp_main        = bestfp_main,
    bestfp_interaction = bestfp_interaction,
    center_vals        = center_vals,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    coefficient_groups = coefficient_groups,
    test_results       = test_results
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
                  min_saz_component_prop = 0.10,
                  group_dummies,
                  center_type = c("grand", "group"),
                  scale_var = 1, shift_var = 0, has_offset,
                  run_test = TRUE, fitter = "base") {

  center_type <- match.arg(center_type)


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

  # x is already shifted upstream so shift is
  # intentionally set to 0. This call is for power
  # selection only; coefficients are not used.
  # Scaling does not affect FP power selection.
  # The final model with correct backscaling is
  # fitted later in test_interaction via fit_model.
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
    fitter        = fitter,
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
    min_saz_component_prop = min_saz_component_prop,
    saz_pre_resolved = TRUE,
    has_offset    = has_offset,
    verbose       = FALSE
  )

  # Step 2: Extract best FP powers for cont_var --------------------------------
  bestfp <- unlist(get_fp_powers(cont_var, fit_main_pool$fp_terms))

  # Defensive guard: if fit_mfp returned an unexpected linear or empty result for
  # an FP1 search (despite force_max_fp = TRUE and power = 1 excluded from
  # fp_cand by flex_fit()), perform a direct deviance-minimization search over
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

    # Fit a one-term FP transformation at each candidate power and choose the
    # lowest-deviance fit. Adjustment columns and group dummies are included as
    # covariates so the comparison is on equal footing with the MFP-fitted model.
    #
    # Although this deviance search uses fast = TRUE, where `offset` is passed
    # directly to the underlying fitting routine, forward `has_offset` as well.
    # This keeps offset handling consistent with the rest of the MFPI fitting
    # path and avoids future inconsistencies if this call is changed to use a
    # formula-based fit.
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
        x             = cbind(z_p, fixed_x),
        y             = y,
        family        = family,
        family_string = family_string,
        fitter        = fitter,
        weights       = weights,
        offset        = offset,
        method        = ties,
        strata        = strata,
        control       = control,
        nocenter      = nocenter,
        rownames      = NULL,
        fast          = TRUE,
        has_offset    = has_offset
      )
      -2 * fit_p$logl
    }, numeric(1L))
    bestfp <- nonunit_cand[which.min(dev_vec)]
  }

  # Step 3: Build interaction terms using the pooled FP powers ----------------
  znames             <- sprintf("%s%d1", cont_var, group_levels)
  bestfp_interaction <- setNames(
    replicate(length(znames), bestfp, simplify = FALSE), znames
  )

  z_vars <- create_z_variables(
    cont_var    = contvar_vec,
    group_var   = groupvar_vec,
    power       = bestfp,
    shift       = 0,
    scale       = 1,
    center      = center,
    zero        = zero_var,
    center_type = center_type,
    scale_var   = scale_var
  )

  x_main        <- cbind(group_dummies, z_vars$xtransformed, xadj)
  x_interaction <- cbind(group_dummies, z_vars$z,            xadj)
  center_vals   <- z_vars$center_vals
  coefficient_groups <- z_vars$column_groups
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
      family_string      = family_string,
      fitter             = fitter,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      has_offset         = has_offset
    )
  }

  list(
    bestfp_main        = bestfp,
    bestfp_interaction = bestfp_interaction,
    center_vals        = center_vals,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    coefficient_groups = coefficient_groups,
    test_results       = test_results
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
                  min_saz_component_prop = 0.10,
                  group_dummies,
                  center_type = c("grand", "group"),
                  scale_var = 1, shift_var = 0, has_offset,
                  run_test = TRUE, fitter = "base") {

  center_type <- match.arg(center_type)


  # Extract column vectors -----------------------------------------------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]

  # Step 1: Generate all candidate within-group transformed variables ----------
  # transform_z_variables() tries every combination of powers from fp_cand,
  # applies the chosen centering (center_type) using the group membership mask
  # correctly, and returns center_vals_list alongside z_transformed.
  transformed <- transform_z_variables(
    cont_var    = contvar_vec,
    group_var   = groupvar_vec,
    shift       = 0,
    scale       = 1,
    fp_cand     = fp_cand,
    fp_degree   = degree,
    acdx        = FALSE,
    center      = center,
    center_type = center_type,
    zero        = zero_var,
    scale_var   = scale_var
  )

  power_matrix       <- transformed$powers_matrix
  transformed_groups <- transformed$z_transformed
  znames             <- transformed$znames
  center_vals_list   <- transformed$center_vals_list   # NULL when center = FALSE
  column_groups_list  <- transformed$column_groups_list

  # Step 2: Select the power combination with the lowest deviance -------------
  # Each candidate within-group FP basis is fitted with the same adjustment
  # columns and group dummies. The candidate with the smallest deviance supplies
  # the common FP powers used for the interaction-side design.
  #
  # The search uses fast = TRUE, so `offset` is passed directly to the underlying
  # fitting routine. We still forward `has_offset` to keep fit_model() calls
  # consistent across the MFPI codebase and to preserve correct behaviour if the
  # implementation later switches to a formula-based fit.
  fixed_x <- cbind(group_dummies, xadj)

  deviance_vec <- vapply(transformed_groups, function(z) {
    if (anyNA(z)) {
      stop("! Internal error: NA values in flex2 candidate design matrix.", call. = FALSE)
    }
    fit <- fit_model(
      x             = cbind(z, fixed_x),
      y             = y,
      family        = family,
      family_string = family_string,
      fitter        = fitter,
      weights       = weights,
      offset        = offset,
      method        = ties,
      strata        = strata,
      control       = control,
      nocenter      = nocenter,
      rownames      = NULL,
      fast          = TRUE,
      has_offset = has_offset
    )
    -2 * fit$logl
  }, numeric(1L))

  best_idx <- which.min(deviance_vec)
  bestfp   <- power_matrix[best_idx, ]
  coefficient_groups <- column_groups_list[[best_idx]]

  # Extract the centering constants for the selected combination --------------
  center_vals <- if (center) center_vals_list[[best_idx]] else NULL

  # Step 3: Build contvar_transformed for the main-effects model --------------
  # Developer note:
  # flex2 selects FP powers from the interaction model, then reuses those same
  # powers in the pooled main-effects model. The pooled main-effect FP basis must
  # therefore use the same scale/zero convention as the interaction-side basis.
  #
  # contvar_vec arrives from mfpi.default() already shifted and divided by the
  # variable-specific scale. Multiplying by scale_var restores the shifted
  # original scale, x + shift, before applying the FP transformation.
  #
  # Do not call create_z_variables() here: flex2 only needs the pooled main-effect
  # basis, not the full group-specific interaction design.
  contvar_main <- contvar_vec

  if (!is.null(scale_var) && length(scale_var) == 1L &&
      is.finite(scale_var) && scale_var != 1) {
    contvar_main <- contvar_main * scale_var
  }

  contvar_transformed <- transform_vector_fp(
    x            = contvar_main,
    power        = bestfp,
    name         = cont_var,
    scale        = 1,
    shift        = 0,
    zero         = zero_var,
    check_binary = FALSE
  )

  contvar_transformed <- as.matrix(contvar_transformed)
  storage.mode(contvar_transformed) <- "double"

  ct_names <- colnames(contvar_transformed)
  if (is.null(ct_names)) {
    ct_names <- paste0(cont_var, seq_len(ncol(contvar_transformed)))
  }
  colnames(contvar_transformed) <- ct_names

  if (center) {
    if (isTRUE(zero_var)) {
      # Developer note:
      # For zero-handled variables, compute centering constants from the positive
      # part only. Otherwise, a large structural-zero mass can dominate the mean
      # and distort the centered FP function for the positive distribution.
      positive_rows <- as.vector(contvar_main) > 0

      if (!any(positive_rows)) {
        stop(
          "! No positive values are available to centre the flex2 main-effect FP transformation.",
          call. = FALSE
        )
      }

      ct_means <- colMeans(
        contvar_transformed[positive_rows, , drop = FALSE],
        na.rm = TRUE
      )
    } else {
      ct_means <- colMeans(contvar_transformed, na.rm = TRUE)
    }

    if (any(!is.finite(ct_means))) {
      stop(
        "! Could not compute finite centering constants for the flex2 main-effect FP transformation.",
        call. = FALSE
      )
    }

    contvar_transformed <- sweep(
      contvar_transformed,
      2L,
      ct_means,
      "-",
      check.margin = FALSE
    )
  }

  # Developer note:
  # Structural-zero rows represent absence of the positive-part FP contribution.
  # Restore them after centering; otherwise centering would assign non-zero
  # artificial FP contributions to structural-zero observations.
  if (isTRUE(zero_var)) {
    zero_rows <- as.vector(contvar_main) <= 0
    if (any(zero_rows)) {
      contvar_transformed[zero_rows, ] <- 0
    }
  }

  z_best <- transformed_groups[[best_idx]]

  if (
    is.null(coefficient_groups) ||
    !identical(
      unname(unlist(coefficient_groups, use.names = FALSE)),
      colnames(z_best)
    )
  ) {
    stop(
      "! Internal flex2 error: coefficient_groups do not match the selected interaction design.",
      call. = FALSE
    )
  }

  # Developer note:
  # transform_z_variables() should return a finite design matrix: inactive group
  # blocks and structural-zero rows are coded as 0. NA, Inf, or -Inf here indicate
  # an internal transformation problem and should not be silently repaired.
  if (any(!is.finite(z_best))) {
    stop(
      "! Non-finite values were produced in the selected flex2 interaction design.",
      call. = FALSE
    )
  }

  # Step 4: Assemble design matrices -------------------------------------------
  x_main        <- cbind(group_dummies, contvar_transformed, xadj)
  x_interaction <- cbind(group_dummies, z_best,              xadj)

  bestfp_interaction <- setNames(
    replicate(length(znames), bestfp, simplify = FALSE), znames
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
      family_string      = family_string,
      fitter             = fitter,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      has_offset         = has_offset
    )
  }

  list(
    bestfp_main        = bestfp,
    bestfp_interaction = bestfp_interaction,
    center_vals        = center_vals,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    coefficient_groups = coefficient_groups,
    test_results       = test_results
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
                  min_saz_component_prop = 0.10,
                  group_dummies,
                  center_type = c("grand", "group"),
                  scale_var = 1, shift_var = 0, has_offset,
                  run_test = TRUE, fitter = "base") {

  center_type <- match.arg(center_type)

  # Shared arguments for flex1 and flex2 calls ---------------------------------
  shared <- list(
    x             = x, y = y, cont_var = cont_var, group_var = group_var,
    group_dummies = group_dummies,
    xadj          = xadj, criterion = criterion, ties = ties, degree = degree,
    family        = family, family_string = family_string, fitter = fitter,
    fp_cand       = fp_cand,
    use_ftest     = use_ftest, center = center, center_type = center_type,
    scale_var     = scale_var, shift_var = shift_var,
    xorder        = xorder,
    weights       = weights, offset = offset, strata = strata,
    control       = control, nocenter = nocenter, cycles = cycles,
    zero_var      = zero_var, spike_var = spike_var,
    min_saz_component_prop = min_saz_component_prop,
    has_offset    = has_offset,
    run_test      = FALSE
  )

  # Step 1: Main-effects design matrix from flex1 (pooled FP powers) ----------
  fit_flex1   <- do.call(flex1, shared)
  x_main      <- fit_flex1$xmain
  bestfp_main <- fit_flex1$bestfp_main

  # Step 2: Interaction design matrix from flex2 (within-group FP powers) -----
  fit_flex2          <- do.call(flex2, shared)
  x_interaction      <- fit_flex2$xinteraction
  bestfp_interaction <- fit_flex2$bestfp_interaction
  center_vals        <- fit_flex2$center_vals
  znames             <- fit_flex2$znames
  coefficient_groups <- fit_flex2$coefficient_groups

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
      family_string      = family_string,
      fitter             = fitter,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      has_offset         = has_offset
    )
  }

  list(
    bestfp_main        = bestfp_main,
    bestfp_interaction = bestfp_interaction,
    center_vals        = center_vals,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    coefficient_groups = coefficient_groups,
    test_results       = test_results
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
#' powers estimated by `flex1`.
#'
#' @inheritParams flex_fit
#' @return See \code{flex_fit()} for the return structure.
#' @keywords internal
#' @noRd
flex4 <- function(x, y, cont_var, group_var, xadj, criterion, ties,
                  degree, family, family_string, fp_cand, use_ftest,
                  center, xorder, weights, offset, strata, control,
                  nocenter, cycles, zero_var, spike_var = FALSE,
                  min_saz_component_prop = 0.10,
                  group_dummies,
                  center_type = c("grand", "group"),
                  scale_var = 1, shift_var = 0, has_offset,
                  run_test = TRUE, fitter = "base") {

  center_type <- match.arg(center_type)


  # Extract column vectors -----------------------------------------------------
  contvar_vec  <- x[, cont_var,  drop = FALSE]
  groupvar_vec <- x[, group_var, drop = FALSE]
  k            <- length(unique(groupvar_vec))

  # Step 1: Build untransformed group-wise variables (linear placeholder) ------
  # Two versions are needed:
  #   z_linear         - built from scaled contvar_vec (no backscaling) for
  #                      power selection in Step 2. Scaling improves numerical
  #                      stability during MFP backfitting cycles.
  #   z_linear_bs      - built from backscaled contvar_vec (* scale_var) for
  #                      the final transform_matrix call in Step 3. This ensures
  #                      interaction model coefficients are on phi(x + shift)
  #                      scale matching mfp2 and the adjustment model.
  z_linear <- create_z_variables(
    cont_var  = contvar_vec,
    group_var = groupvar_vec,
    power     = 1L,
    shift     = 0,
    scale     = 1,
    center    = FALSE,
    scale_var = 1          # no backscaling for power selection
  )$z
  znames <- colnames(z_linear)

  # Backscaled version for final model fit and centering
  z_linear_bs <- create_z_variables(
    cont_var  = contvar_vec,
    group_var = groupvar_vec,
    power     = 1L,
    shift     = 0,
    scale     = 1,
    center    = FALSE,
    scale_var = scale_var  # backscale: cont_var * scale_var = x + shift
  )$z

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
  # Developer note: Non-positive values within groups are handled via zero = TRUE"),
  # this is unrelated to the user's zero_var flag.
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
    fitter        = fitter,
    powers        = pw_list,
    shift         = rep(0, n_total),   # x already shifted upstream
    scale         = rep(1, n_total),   # intentionally 1: power selection only;
    ftest         = use_ftest,         # coefficients not used from this call.
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
    min_saz_component_prop = min_saz_component_prop,
    saz_pre_resolved = TRUE,
    has_offset    = has_offset,
    verbose       = FALSE
  )

  bestfp_interaction <- get_fp_powers(znames, fit_int$fp_terms)

  # Step 3: Apply group-specific FP powers to build the interaction matrix -----
  # Use z_linear_bs (backscaled) so the final model is on phi(x + shift) scale.
  transformed <- transform_matrix(
    x          = z_linear_bs,
    power_list = bestfp_interaction,
    acdx       = setNames(rep(FALSE, k), znames),
    center     = setNames(rep(FALSE, k), znames),   # centring done below
    zero       = setNames(rep(TRUE,  k), znames),
    catzero    = setNames(rep(FALSE, k), znames)
  )$x_transformed

  # transform_matrix() should already return a finite matrix here.  In
  # particular, structural-zero rows for zero-handled variables are coded as 0,
  # not as Inf/NA.  Treat any remaining non-finite value as an internal
  # transformation error rather than silently repairing it.
  if (any(!is.finite(transformed))) {
    stop(
      "! Internal error. Non-finite values were produced in the selected flex4 interaction design.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Centering for flex4
  # ---------------------------------------------------------------------------
  # flex4 allows each group to have its own FP powers.  Therefore each group
  # block needs centering constants computed with that group's selected powers.
  #
  # The fitting design must match create_z_variables()/transform_z_variables():
  #   * compute centers from the positive part only when zero_var = TRUE;
  #   * apply centering only to rows belonging to the corresponding group block;
  #   * restore structural-zero rows to 0 after centering;
  #   * fail on unexpected non-finite transformed values.

  group_levels <- sort(unique(as.vector(groupvar_vec)))
  group_idx    <- match(as.vector(groupvar_vec), group_levels)

  # Backscaled full cont_var for centering.  contvar_vec has already been
  # shifted and divided by scale_var upstream; multiplying by scale_var restores
  # x + shift, matching the basis used in z_linear_bs and transform_matrix().
  contvar_bs <- contvar_vec * scale_var
  zero_rows  <- if (isTRUE(zero_var)) as.vector(contvar_bs) <= 0 else rep(FALSE, nrow(contvar_bs))
  valid_rows <- !zero_rows

  if (!any(valid_rows)) {
    stop(
      "! No valid rows are available for flex4 FP transformation/centering.",
      call. = FALSE
    )
  }

  if (center) {
    center_vals <- numeric(ncol(transformed))
    names(center_vals) <- colnames(transformed)

    for (gi in seq_len(k)) {
      grp_name <- znames[gi]
      pw       <- bestfp_interaction[[grp_name]]
      n_terms  <- length(pw)
      cols_g   <- seq.int((gi - 1L) * n_terms + 1L, gi * n_terms)
      in_grp   <- group_idx == gi

      # FP-transform the full backscaled cont_var with this group's powers.
      # Passing zero = zero_var is essential: otherwise log/negative powers of
      # structural zeros can create non-finite centers, and zeros can influence
      # the mean even though they represent absence of the positive-part FP term.
      x_fp_full <- transform_vector_fp(
        x            = contvar_bs,
        power        = pw,
        shift        = 0,
        scale        = 1,
        zero         = isTRUE(zero_var),
        check_binary = FALSE
      )

      x_fp_full <- as.matrix(x_fp_full)
      storage.mode(x_fp_full) <- "double"

      if (ncol(x_fp_full) != n_terms) {
        stop(
          "! Internal flex4 centering error: transformed column count does not match selected powers.",
          call. = FALSE
        )
      }

      if (any(!is.finite(x_fp_full[valid_rows, , drop = FALSE]))) {
        stop(
          "! Non-finite values were produced among valid rows during flex4 centering.",
          call. = FALSE
        )
      }

      center_rows <- if (center_type == "grand") valid_rows else in_grp & valid_rows

      if (!any(center_rows)) {
        stop(
          paste0(
            "! Cannot compute flex4 centering constants for group ",
            group_levels[gi], ": no valid rows."
          ),
          call. = FALSE
        )
      }

      centers_g <- colMeans(x_fp_full[center_rows, , drop = FALSE])

      if (any(!is.finite(centers_g))) {
        stop(
          paste0(
            "! Could not compute finite flex4 centering constants for group ",
            group_levels[gi], "."
          ),
          call. = FALSE
        )
      }

      center_vals[cols_g] <- centers_g

      # Centre active rows in this group block only.  Out-of-group rows remain
      # zero because they encode inactive group blocks, not observed values.
      transformed[in_grp, cols_g] <- sweep(
        transformed[in_grp, cols_g, drop = FALSE],
        2L,
        centers_g,
        "-",
        check.margin = FALSE
      )

      # Structural zeros encode absence of the positive-part FP contribution.
      # They must stay at 0 after centering; otherwise they become -center.
      if (isTRUE(zero_var) && any(in_grp & zero_rows)) {
        transformed[in_grp & zero_rows, cols_g] <- 0
      }
    }
  } else {
    center_vals <- NULL
  }

  x_interaction <- cbind(group_dummies, transformed, xadj)

  # group names
  coefficient_groups <- stats::setNames(
    vector("list", k),
    as.character(group_levels)
  )

  col_pos <- 1L

  for (gi in seq_len(k)) {
    grp_name <- znames[gi]
    n_cols_g <- length(bestfp_interaction[[grp_name]])

    cols_g <- seq.int(col_pos, col_pos + n_cols_g - 1L)

    coefficient_groups[[gi]] <- colnames(transformed)[cols_g]

    col_pos <- col_pos + n_cols_g
  }

  if (col_pos != ncol(transformed) + 1L) {
    stop(
      "! Internal flex4 error: coefficient_groups do not cover the full interaction design.",
      call. = FALSE
    )
  }

  # Step 4: Main-effects design matrix via flex1 (pooled FP powers) -----------
  fit_main <- flex1(
    x             = x, y = y, cont_var = cont_var, group_var = group_var,
    group_dummies = group_dummies,
    xadj          = xadj, criterion = criterion, ties = ties, degree = degree,
    family        = family, family_string = family_string, fitter = fitter,
    fp_cand       = fp_cand, use_ftest = use_ftest,
    center        = center, center_type = center_type,
    scale_var     = scale_var,
    xorder        = xorder, weights = weights, offset = offset,
    strata        = strata, control = control, nocenter = nocenter,
    cycles        = cycles, zero_var = zero_var, spike_var = spike_var,
    min_saz_component_prop = min_saz_component_prop,
    has_offset    = has_offset,
    run_test = FALSE
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
      family_string      = family_string,
      fitter             = fitter,
      weights            = weights,
      offset             = offset,
      ties               = ties,
      strata             = strata,
      control            = control,
      nocenter           = nocenter,
      has_offset         = has_offset
    )
  }

  list(
    bestfp_main        = bestfp_main,
    bestfp_interaction = bestfp_interaction,
    center_vals        = center_vals,
    xmain              = x_main,
    xinteraction       = x_interaction,
    znames             = znames,
    coefficient_groups = coefficient_groups,
    test_results       = test_results
  )
}