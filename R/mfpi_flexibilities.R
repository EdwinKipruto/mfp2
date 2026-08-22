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
#' @param ties Character string; Cox tie-handling - `"breslow"` or `"efron"`.
#'   `"exact"` is rejected by the public MFP/MFPI interfaces before selection.
#'   Ignored for non-Cox families.
#' @param degree Integer; FP degree - `0` (linear), `1` (FP1), or `2` (FP2).
#'   Values below 1 are silently treated as `0` and routed to `flex0`.
#' @param family Character string; `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"negbin"`, or `"cox"`.
#' @param fp_cand Numeric vector of candidate FP powers for `cont_var`.
#'   Corresponds to one element of the `fp_powers` list in [mfp2::mfpi()]. For
#'   FP1, the conventional candidate set is retained, including power = 1.
#'   A prespecified linear interaction is still handled separately by `flex0`;
#'   thus `linear` means power 1 fixed with no power search, whereas an FP1
#'   search may select power 1 as its best member.
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
#'   as computed and applied in \code{mfpi.default()}. The flex functions
#'   multiply the pre-scaled covariate by \code{scale_var} before final FP
#'   transformation so interaction-model coefficients are on the
#'   \eqn{\phi(x + \text{shift})} scale used by the adjustment model and
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
  n_total <- length(vnames)

  df_vec <- stats::setNames(
    c(
      2L * degree,               # cont_var: FP of requested degree
      rep(1L, k - 1L),           # group dummies: linear
      if (n_adj > 0L) rep(1L, n_adj)
    ),
    vnames
  )

  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  pw_list <- setNames(replicate(length(vnames), default_powers, simplify = FALSE),
                      vnames)
  pw_list[[cont_var]] <- fp_cand

  select_vec <- stats::setNames(rep(1, n_total), vnames)
  alpha_vec  <- stats::setNames(rep(1, n_total), vnames)
  center_vec <- stats::setNames(
    c(center, rep(FALSE, n_total - 1L)),
    vnames
  )
  shift_vec <- stats::setNames(rep(0, n_total), vnames)
  scale_vec <- stats::setNames(rep(1, n_total), vnames)
  zero_vec <- stats::setNames(
    c(zero_var, rep(FALSE, n_total - 1L)),
    vnames
  )
  spike_vec <- stats::setNames(
    c(spike_var, rep(FALSE, n_total - 1L)),
    vnames
  )
  acd_vec <- stats::setNames(rep(FALSE, n_total), vnames)
  catzero <- stats::setNames(rep(FALSE, n_total), vnames)
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
    select        = select_vec,
    alpha         = alpha_vec,
    keep          = vnames,
    force_max_fp  = force_max_fp,
    # MFPI follows the conventional FP1 search: if power 1 is present in the
    # supplied FP1 candidate set, retain it as a possible selected power. This
    # internal switch is used only for this MFPI power-selection fit; ordinary
    # mfp2() closed-test behaviour is unchanged.
    retain_linear_fp1 = TRUE,
    method        = ties,
    family        = family,
    family_string = family_string,
    fitter        = fitter,
    powers        = pw_list,
    ftest         = use_ftest,
    center        = center_vec,
    shift         = shift_vec,
    scale         = scale_vec,
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

  # Defensive guard: the forced fixed-degree MFPI search must return exactly
  # `degree` selected powers. For FP1, power 1 is a valid selected result: it
  # means the FP1 search selected its linear member, not that the user requested
  # the separate prespecified `linear` interaction form.
  if (length(bestfp) != degree || anyNA(bestfp)) {
    stop(
      "! Internal MFPI error: FP power selection did not return the requested degree.",
      call. = FALSE
    )
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
#' @details
#' The FP candidate search uses the same compact-basis representation as the
#' optimized `mfp2()` selector. Instead of materializing one full
#' `n x (number of groups * FP degree)` interaction matrix for every candidate
#' power combination, the function computes each distinct FP basis column once.
#' `generate_transformations_fp_basis_cpp()` returns (i) that shared basis and
#' (ii) an integer `candidate_map` telling `flex2()` which basis columns belong
#' to each power combination.
#'
#' For every candidate, `fill_mfpi_fp_candidate_cpp()` scatters only the selected
#' basis columns into the active group block of one reusable design matrix. The
#' same helper also applies grand-mean or within-group centering and preserves
#' structural-zero rows. This avoids retaining all candidate interaction
#' matrices and avoids repeated `cbind()` allocation in the model-selection loop.
#'
#' After the candidate with the smallest deviance is identified, the reusable
#' interaction block is filled once more for the winning candidate. For ordinary
#' continuous covariates, the pooled main-effect FP columns are taken directly
#' from the same compact basis. Binary covariates use `transform_vector_fp()` for
#' the pooled main effect because the compact candidate map deliberately collapses
#' the interaction-side binary basis to a single passthrough column.
#'
#' Structural-zero rows are excluded from centering means and remain exact zeros
#' in the final group-specific design. Grand centering uses one mean per FP term
#' across all valid rows; group centering uses a separate mean for each group and
#' FP term.
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

  # Step 1: Build one compact FP basis shared by every candidate ---------------
  # Store each distinct FP basis term once together with a small integer map
  # from candidate power combinations to basis columns. Each candidate is then
  # scattered into one reusable group-specific fitting matrix.
  # `contvar_vec` has already been shifted and divided by `scale_var` by the
  # upstream MFPI preprocessing. The compact FP basis must be evaluated on the
  # required x + shift scale used by both MFPI model designs.
  # Multiplication is therefore a restoration step, not an additional rescaling.
  contvar_main <- contvar_vec

  if (!is.null(scale_var) && length(scale_var) == 1L &&
      is.finite(scale_var) && scale_var != 1) {
    contvar_main <- contvar_main * scale_var
  }

  # Convert the one-column matrices to simple vectors before crossing the R/C++
  # boundary. `group_idx` is deliberately converted to a dense 1..G integer
  # index because the C++ scatter kernel uses it for direct block addressing.
  contvar_main_vec <- as.numeric(contvar_main)
  group_vec        <- as.vector(groupvar_vec)
  group_levels     <- sort(unique(group_vec))
  group_idx        <- match(group_vec, group_levels)
  group_idx_int    <- as.integer(group_idx)
  n_groups         <- length(group_levels)
  n_groups_int     <- as.integer(n_groups)
  n                 <- length(contvar_main_vec)

  # Structural-zero semantics: when `zero_var` is active, x <= 0 does not
  # participate in the positive-part FP transform or in centering. Those rows
  # are later written back as exact zeros in every group-specific FP block.
  zero_rows  <- if (isTRUE(zero_var)) contvar_main_vec <= 0 else rep(FALSE, n)
  valid_rows <- !zero_rows

  if (!any(valid_rows)) {
    stop(
      "! No valid rows are available for flex2 FP transformation/centering.",
      call. = FALSE
    )
  }

  powers_matrix <- generate_powers_fp(
    degree = degree,
    powers = fp_cand
  )
  powers_matrix <- as.matrix(powers_matrix)
  storage.mode(powers_matrix) <- "double"

  # Compact representation
  # ----------------------
  # `basis`: n x B matrix containing each distinct FP basis term exactly once.
  # `candidate_map`: C x d integer matrix; row i gives the 1-based basis columns
  #                  needed for candidate i. B is typically much smaller than
  #                  C * d because candidates reuse the same powers.
  compact <- generate_transformations_fp_basis_cpp(
    x      = contvar_main_vec,
    powers = powers_matrix,
    zero   = isTRUE(zero_var)
  )

  basis         <- compact$basis
  candidate_map <- compact$candidate_map
  n_candidates  <- nrow(powers_matrix)
  n_terms       <- ncol(candidate_map)

  if (nrow(candidate_map) != n_candidates || n_terms < 1L) {
    stop(
      "! Internal flex2 error: compact FP candidate map is inconsistent with the power matrix.",
      call. = FALSE
    )
  }

  xname <- colnames(contvar_vec)
  if (is.null(xname) || length(xname) != 1L) {
    xname <- cont_var
  }

  # Column order is group-major: all FP terms for group 1, then all terms for
  # group 2, and so on. The C++ helper and `coefficient_groups` use this exact
  # layout, so keep the naming logic synchronized with `focal_width`.
  z_names <- paste0(
    xname,
    rep(as.character(group_levels), each = n_terms),
    rep(seq_len(n_terms), times = n_groups)
  )

  # `znames` denotes the conceptual group-specific continuous variables and is
  # therefore one name per group, irrespective of the FP degree.
  znames <- paste0(xname, as.character(group_levels), "1")

  coefficient_groups <- stats::setNames(
    lapply(seq_len(n_groups), function(g) {
      cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
      z_names[cols_g]
    }),
    as.character(group_levels)
  )

  # Step 2: Select powers with one reusable search design ----------------------
  fixed_x <- cbind(group_dummies, xadj)
  focal_width <- n_groups * n_terms

  # Allocate the model matrix once. The leading focal block starts at zero and
  # is overwritten candidate-by-candidate; the adjustment/group columns to its
  # right never change. Candidate-search memory therefore stays bounded by one
  # focal design matrix rather than growing with the number of FP candidates.
  search_x <- cbind(
    matrix(0, nrow = n, ncol = focal_width),
    fixed_x
  )
  colnames(search_x) <- c(z_names, colnames(fixed_x))

  deviance_vec <- numeric(n_candidates)

  for (i in seq_len(n_candidates)) {
    # The kernel mutates/rebuilds only the active group-specific focal columns.
    # It returns `target` explicitly because Rcpp's copy-on-write boundary should
    # not be relied on as an undocumented in-place side effect at the R level.
    candidate <- fill_mfpi_fp_candidate_cpp(
      target           = search_x,
      basis            = basis,
      source_cols      = as.integer(candidate_map[i, ]),
      group_idx        = group_idx_int,
      n_groups         = n_groups_int,
      target_start_col = 1L,
      center           = isTRUE(center),
      group_center     = isTRUE(center) && center_type == "group",
      valid_rows       = valid_rows
    )
    search_x <- candidate$target

    fit <- fit_model(
      x             = search_x,
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
    deviance_vec[i] <- -2 * fit$logl
  }

  # Deviance is -2 log-likelihood, so the smallest value identifies the
  # preferred common FP power vector across all groups.
  best_idx <- which.min(deviance_vec)
  bestfp   <- as.numeric(powers_matrix[best_idx, ])

  # Refill the reusable focal block for the selected candidate only. This gives
  # us the final interaction design and its centering constants without retaining
  # all candidate matrices or all candidate center vectors.
  selected <- fill_mfpi_fp_candidate_cpp(
    target           = search_x,
    basis            = basis,
    source_cols      = as.integer(candidate_map[best_idx, ]),
    group_idx        = group_idx_int,
    n_groups         = n_groups_int,
    target_start_col = 1L,
    center           = isTRUE(center),
    group_center     = isTRUE(center) && center_type == "group",
    valid_rows       = valid_rows
  )
  search_x <- selected$target

  center_vals <- if (center) {
    out <- as.numeric(selected$centers)
    names(out) <- z_names
    out
  } else {
    NULL
  }

  z_best <- search_x[, seq_len(focal_width), drop = FALSE]
  colnames(z_best) <- z_names

  if (any(!is.finite(z_best))) {
    stop(
      "! Non-finite values were produced in the selected flex2 interaction design.",
      call. = FALSE
    )
  }

  if (!identical(
    unname(unlist(coefficient_groups, use.names = FALSE)),
    colnames(z_best)
  )) {
    stop(
      "! Internal flex2 error: coefficient_groups do not match the selected interaction design.",
      call. = FALSE
    )
  }

  # Step 3: Reuse the winning columns for the pooled main-effects model --------
  # For ordinary continuous variables, each selected FP term maps directly to
  # one compact-basis column, so no second FP transformation is required after
  # model selection. For binary covariates, the compact interaction basis can
  # collapse to one passthrough column; the pooled main effect is therefore
  # transformed directly with check_binary = FALSE.
  if (ncol(candidate_map) == ncol(powers_matrix)) {
    contvar_transformed <- basis[
      , as.integer(candidate_map[best_idx, ]), drop = FALSE
    ]
  } else {
    contvar_transformed <- transform_vector_fp(
      x            = contvar_main,
      power        = bestfp,
      name         = cont_var,
      scale        = 1,
      shift        = 0,
      zero         = zero_var,
      check_binary = FALSE
    )
  }

  contvar_transformed <- as.matrix(contvar_transformed)
  storage.mode(contvar_transformed) <- "double"
  colnames(contvar_transformed) <- name_transformed_variables(
    cont_var,
    ncol(contvar_transformed)
  )

  if (center) {
    # Main-effect centering is pooled across groups. When structural zeros are
    # enabled, use only valid positive-part rows so structural zeros do not
    # contribute to the pooled centering constants.
    center_rows <- if (isTRUE(zero_var)) valid_rows else rep(TRUE, n)
    ct_means <- colMeans(
      contvar_transformed[center_rows, , drop = FALSE],
      na.rm = TRUE
    )

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

  if (isTRUE(zero_var) && any(zero_rows)) {
    contvar_transformed[zero_rows, ] <- 0
  }

  # Step 4: Assemble final design matrices -------------------------------------
  x_main        <- cbind(group_dummies, contvar_transformed, xadj)
  x_interaction <- cbind(group_dummies, z_best,              xadj)

  bestfp_interaction <- setNames(
    replicate(length(znames), bestfp, simplify = FALSE),
    znames
  )

  # Step 5: Interaction test ---------------------------------------------------
  test_results <- NULL
  if (run_test) {
    test_results <- test_interaction(
      y                  = y,
      cont_var           = contvar_vec,
      group_var          = groupvar_vec,
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
#' @details
#' The final interaction design is built once with the selected group-specific
#' powers. Centering is then applied without unnecessarily re-transforming the
#' full continuous variable for every group. With `center_type = "group"`, the
#' active rows of each already-transformed group block contain exactly the basis
#' needed to compute that group's centering constants, so their column means are
#' used directly. With `center_type = "grand"`, the selected FP basis must be
#' evaluated on all valid observations; these full-variable transformations are
#' cached by an exact power-vector key so groups selecting identical powers share
#' the same result. Structural-zero rows are excluded from centering and restored
#' to zero after subtraction.
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

  df_vec <- stats::setNames(
    c(
      rep(2L * degree, k),     # one FP term per group
      rep(1L, k - 1L),         # group dummies
      if (n_adj > 0L) rep(1L, n_adj)
    ),
    vnames
  )

  default_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  pw_list <- setNames(replicate(n_total, default_powers, simplify = FALSE), vnames)
  pw_list[znames] <- replicate(k, fp_cand, simplify = FALSE)

  select_vec <- stats::setNames(rep(1, n_total), vnames)
  alpha_vec  <- stats::setNames(rep(1, n_total), vnames)
  # Centre only the within-group cont_var columns; adjustment already centred
  center_vec <- stats::setNames(
    c(rep(center, k), rep(FALSE, n_total - k)),
    vnames
  )
  shift_vec <- stats::setNames(rep(0, n_total), vnames)
  scale_vec <- stats::setNames(rep(1, n_total), vnames)
  # Developer note: Non-positive values within groups are handled via zero = TRUE"),
  # this is unrelated to the user's zero_var flag.
  zero_vec <- stats::setNames(
    c(rep(TRUE, k), rep(FALSE, n_total - k)),
    vnames
  )
  spike_vec <- stats::setNames(
    c(rep(spike_var, k), rep(FALSE, n_total - k)),
    vnames
  )
  acd_vec     <- stats::setNames(rep(FALSE, n_total), vnames)
  catzero_vec <- stats::setNames(rep(FALSE, n_total), vnames)
  # force_max_fp: named logical vector, TRUE for the k per-group cont_var
  # columns, FALSE for group dummies and adjustment variables (df = 1,
  # nothing to force). Same rationale as flex1.
  force_max_fp <- setNames(c(rep(TRUE, k), rep(FALSE, n_total - k)), vnames)


  fit_int <- fit_mfp(
    x             = xnew,
    y             = y,
    df            = df_vec,
    criterion     = criterion,
    select        = select_vec,
    alpha         = alpha_vec,
    keep          = vnames,
    force_max_fp  = force_max_fp,
    # MFPI follows the conventional FP1 search: if power 1 is present in the
    # supplied FP1 candidate set, retain it as a possible selected power. This
    # internal switch is used only for this MFPI power-selection fit; ordinary
    # mfp2() closed-test behaviour is unchanged.
    retain_linear_fp1 = TRUE,
    method        = ties,
    family        = family,
    family_string = family_string,
    fitter        = fitter,
    powers        = pw_list,
    shift         = shift_vec,         # x already shifted upstream
    scale         = scale_vec,         # intentionally 1: power selection only;
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
  # The fitting design follows the same structural-zero and centering rules as
  # the other MFPI group-specific FP designs:
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

    # Grand centering needs each selected FP basis evaluated on the full
    # continuous variable. Cache by power vector so groups that selected the
    # same FP do not repeat the transformation. Group centering does not need a
    # full-variable transformation at all: the active rows of `transformed`
    # already contain exactly that group's selected FP basis.
    # An environment provides O(1)-style lookup without copying cached matrices.
    # Keys are deterministic text encodings of the selected power vectors.
    grand_basis_cache <- new.env(parent = emptyenv())

    for (gi in seq_len(k)) {
      grp_name <- znames[gi]
      pw       <- bestfp_interaction[[grp_name]]
      n_terms  <- length(pw)
      cols_g   <- seq.int((gi - 1L) * n_terms + 1L, gi * n_terms)
      in_grp   <- group_idx == gi
      group_valid_rows <- in_grp & valid_rows

      if (center_type == "group") {
        # No new FP transform is necessary here: `transformed` already stores
        # this group's selected basis in `cols_g`, and only `in_grp` rows are
        # active for that block. Restricting to `group_valid_rows` also excludes
        # structural zeros from the centering denominator.
        if (!any(group_valid_rows)) {
          stop(
            paste0(
              "! Cannot compute flex4 centering constants for group ",
              group_levels[gi], ": no valid rows."
            ),
            call. = FALSE
          )
        }

        centers_g <- colMeans(
          transformed[group_valid_rows, cols_g, drop = FALSE]
        )
      } else {
        # Grand centering requires the group's selected FP basis on the full
        # valid sample rather than only its active group rows. Cache that full
        # transform so identical selected power vectors are evaluated once.
        power_key <- paste(
          formatC(as.numeric(pw), digits = 17L, format = "fg", flag = "#"),
          collapse = "|"
        )

        if (!exists(power_key, envir = grand_basis_cache, inherits = FALSE)) {
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

          assign(power_key, x_fp_full, envir = grand_basis_cache)
        }

        x_fp_full <- get(power_key, envir = grand_basis_cache, inherits = FALSE)
        centers_g <- colMeans(x_fp_full[valid_rows, , drop = FALSE])
      }

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

      # Centre active rows in this group block only. Out-of-group rows remain
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