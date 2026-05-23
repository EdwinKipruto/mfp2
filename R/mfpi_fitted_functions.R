# Fitted-function computation for MFPI
#
# gen_fitted_values_per_group() computes group-specific FP fitted functions,
# their pointwise standard errors, and 95% confidence intervals from a fitted
# interaction model. It is called by flex0–flex4 whenever compute_fitted = TRUE
# and is intended to feed the plotting layer and, eventually, predict.mfpi().
#
# The `data_input` argument (original vs. equidistant grid) belongs in
# predict.mfpi() once that method is implemented; it is retained here for
# backward compatibility but marked clearly in the documentation.
#
# Naming conventions match the rest of the package:
#   cont_var          — the continuous variable matrix
#   group_var         — the grouping variable matrix
#   group_fp_powers   — named list of per-group FP powers (was `powers`)
#   interaction_model — fitted model object from test_interaction()


# -----------------------------------------------------------------------------
# gen_fitted_values_per_group() -----------------------------------------------
# -----------------------------------------------------------------------------

#' Compute Group-Specific Fitted FP Functions and Confidence Intervals
#'
#' Given a fitted interaction model and the estimated FP powers for each group,
#' computes the fitted linear predictor \eqn{\hat{f}_j(x)} for each level
#' \eqn{j} of `group_var`, together with pointwise standard errors and 95%
#' confidence intervals. The function-difference \eqn{\hat{f}_j - \hat{f}_0}
#' (relative to the reference group) and its standard error are also returned.
#'
#' Standard errors for the differences are computed via the delta method:
#' \deqn{\text{Var}(\hat{f}_j - \hat{f}_0) =
#'   \mathbf{d}^\top \, \text{Cov}(\hat{\boldsymbol{\beta}}) \, \mathbf{d},}
#' where \eqn{\mathbf{d}} is the vector of partial derivatives of the
#' difference with respect to the model coefficients.
#'
#' @note The `use_grid` argument (equidistant evaluation grid) is a temporary
#'   placeholder. Once `predict.mfpi()` is implemented, grid-based prediction
#'   should be handled there and this argument should be removed.
#'
#' @param cont_var A one-column numeric matrix of the continuous variable.
#'   Must have a column name and have been shifted if it contains
#'   non-positive values.
#' @param group_fp_powers Named list of FP power vectors, one per group of
#'   `group_var`. The names must follow the `<varname><group><power_index>`
#'   convention produced by \code{create_z_variables()}.
#' @param interaction_model The fitted interaction model object returned by
#'   `test_interaction()`. Must expose `$coefficients` and support [vcov()].
#' @param group_var A one-column numeric matrix of the grouping variable.
#'   Must have a column name.
#' @param family Character string; the regression family used to fit the
#'   interaction model — `"gaussian"`, `"binomial"`, `"poisson"`, or `"cox"`.
#'   For Cox models the intercept is set to zero.
#' @param transform Logical. Whether to FP-transform `cont_var` using
#'   `group_fp_powers` before computing fitted values. Default `TRUE`.
#' @param center Logical. Whether to mean-centre the transformed variables.
#'   Should match the `center` setting used in [mfp2::mfpi()]. Default `FALSE`.
#' @param use_grid Logical. If `TRUE`, replaces `cont_var` with a 200-point
#'   equidistant sequence spanning its observed range before computing fitted
#'   values. Useful for smooth plotting curves. Default `FALSE`.
#'   **This argument is provisional and will move to `predict.mfpi()` in a
#'   future release.**
#'
#' @return A numeric matrix with one row per observation (or per grid point
#'   when `use_grid = TRUE`) and the following columns, where \eqn{j} ranges
#'   over all group levels:
#' \describe{
#'   \item{`<varname>`}{The (possibly grid-replaced) values of `cont_var`.}
#'   \item{`f0`, `f1`, …}{Fitted linear predictor for each group.}
#'   \item{`se(f0)`, `se(f1)`, …}{Pointwise standard errors.}
#'   \item{`f0_lower`, `f0_upper`, …}{95% confidence interval bounds.}
#'   \item{`f1-f0`, `f2-f0`, …}{Pointwise difference relative to group 0.}
#'   \item{`se(f1-f0)`, …}{Standard errors of the differences.}
#'   \item{`(f1-f0)_lower`, `(f1-f0)_upper`, …}{95% CI bounds for differences.}
#' }
#'
#' @seealso \code{create_z_variables()}, \code{compute_std_errors_diff()},
#'   \code{var_group()}
#'
#' @keywords internal
#' @noRd
gen_fitted_values_per_group <- function(cont_var,
                                        group_fp_powers,
                                        interaction_model,
                                        group_var,
                                        family,
                                        family_string,
                                        transform  = TRUE,
                                        center     = FALSE,
                                        use_grid   = FALSE) {
  
  # Input validation -----------------------------------------------------------
  if (!is.matrix(cont_var) || ncol(cont_var) != 1L)
    stop("`cont_var` must be a one-column numeric matrix.", call. = FALSE)
  cont_name <- colnames(cont_var)
  if (is.null(cont_name))
    stop("`cont_var` must have a column name.", call. = FALSE)
  
  if (!is.matrix(group_var) || ncol(group_var) != 1L)
    stop("`group_var` must be a one-column numeric matrix.", call. = FALSE)
  group_name <- colnames(group_var)
  if (is.null(group_name))
    stop("`group_var` must have a column name.", call. = FALSE)
  
  if (!is.list(group_fp_powers) || length(group_fp_powers) == 0L)
    stop("`group_fp_powers` must be a non-empty named list.", call. = FALSE)
  
  # Extract model coefficients -------------------------------------------------
  coef_vec   <- interaction_model$coefficients
  coef_names <- names(coef_vec)
  
  # Cox models have no intercept; all others include one
  intercept <- if (family_string == "cox") 0 else coef_vec["(Intercept)"]
  
  # Identify which coefficient blocks belong to each group --------------------
  groups     <- var_group(cont_name, coef_names)
  n_groups   <- length(groups)
  grp_levels <- sort(unique(as.vector(group_var)))
  
  if (n_groups != length(grp_levels)) {
    stop(
      paste0("Number of coefficient groups (", n_groups, ") does not match ",
             "the number of levels in `group_var` (", length(grp_levels), ")."),
      call. = FALSE
    )
  }
  
  # Optionally replace cont_var with an equidistant grid ----------------------
  # TODO: move this block to predict.mfpi() once that method is implemented.
  if (use_grid) {
    nonzero     <- as.vector(cont_var)[as.vector(cont_var) != 0]
    grid_range  <- range(nonzero, na.rm = TRUE)
    cont_var    <- matrix(seq(grid_range[1L], grid_range[2L], length.out = 200L),
                          ncol = 1L)
    colnames(cont_var) <- cont_name
  }
  
  n_obs <- nrow(cont_var)
  
  # Replicate cont_var once per group (same x, different FP powers) -----------
  x_split           <- cont_var[, rep(1L, n_groups), drop = FALSE]
  colnames(x_split) <- sprintf("%s%d%d", cont_name,
                               rep(grp_levels, each = 1L),
                               rep(1L,         times = n_groups))
  znames <- colnames(x_split)
  
  # FP transformation ----------------------------------------------------------
  if (transform) {
    center_map <- stats::setNames(rep(FALSE, n_groups), znames)
    acd_map    <- stats::setNames(rep(FALSE, n_groups), znames)
    
    x_split <- mfp2::transform_matrix(
      x          = x_split,
      power_list = group_fp_powers,
      center     = center_map,
      acdx       = acd_map
    )$x_transformed
    
    # Align column names with the corresponding model coefficients
    x_split_names <- grep(paste0("^", cont_name, "\\d+"), coef_names, value = TRUE)
    colnames(x_split) <- x_split_names
  }
  
  # Optional mean-centring of transformed variables ---------------------------
  if (center) {
    col_means <- colMeans(x_split, na.rm = TRUE)
    x_split   <- sweep(x_split, 2L, col_means, "-", check.margin = FALSE)
  }
  
  # Dummy-variable names for the grouping variable (non-reference levels) -----
  group_dummy_names <- grep(paste0("^", group_name, "\\d+"), coef_names,
                            value = TRUE)
  
  # Variance-covariance matrix of the fitted model ----------------------------
  vcov_mat <- vcov(interaction_model$fit)
  
  # Allocate output matrices --------------------------------------------------
  make_mat <- function(col_fmt) {
    m <- matrix(0, nrow = n_obs, ncol = n_groups)
    colnames(m) <- sprintf(col_fmt, grp_levels)
    m
  }
  
  fitted_vals    <- make_mat("f%d")
  fitted_var     <- make_mat("f%d_var")
  fitted_lower   <- make_mat("f%d_lower")
  fitted_upper   <- make_mat("f%d_upper")
  
  crit_val <- qnorm(0.975)   # 1.96 for 95% CI
  
  # Compute fitted values and variances per group -----------------------------
  for (i in seq_len(n_groups)) {
    grp_coef_names <- groups[[i]]
    grp_coef       <- coef_vec[grp_coef_names]
    x_grp          <- x_split[, grp_coef_names, drop = FALSE]
    
    # Build the predictor matrix used for variance calculation.
    # It extends x_grp with a column of 1s for the intercept (non-Cox) and,
    # for non-reference groups, a column of 1s for the group dummy coefficient.
    x_var       <- x_grp
    x_var_names <- grp_coef_names
    
    if (!identical(intercept, 0)) {
      x_var       <- cbind(1, x_var)
      x_var_names <- c("(Intercept)", x_var_names)
    }
    
    group_coef_offset <- 0
    if (i > 1L) {
      group_coef_offset <- coef_vec[group_dummy_names[i - 1L]]
      x_var             <- cbind(x_var, 1)
      x_var_names       <- c(x_var_names, group_dummy_names[i - 1L])
    }
    
    # Linear predictor: intercept + group-specific FP terms + group offset
    fitted_vals[, i] <- intercept +
      as.vector(x_grp %*% grp_coef) +
      group_coef_offset
    
    # Pointwise variance via quadratic form: diag(X V X^T)
    vcov_sub <- vcov_mat[x_var_names, x_var_names, drop = FALSE]
    fitted_var[, i] <- vapply(
      seq_len(n_obs),
      function(k) {
        xk <- x_var[k, , drop = FALSE]
        as.numeric(xk %*% vcov_sub %*% t(xk))
      },
      numeric(1L)
    )
  }
  
  # Standard errors and confidence intervals for fitted values ----------------
  fitted_se    <- sqrt(fitted_var)
  colnames(fitted_se) <- sprintf("se(f%d)", grp_levels)
  
  fitted_lower <- fitted_vals - crit_val * fitted_se
  colnames(fitted_lower) <- sprintf("f%d_lower", grp_levels)
  
  fitted_upper <- fitted_vals + crit_val * fitted_se
  colnames(fitted_upper) <- sprintf("f%d_upper", grp_levels)
  
  # Function differences relative to the reference group (group 0) ------------
  ref_col   <- fitted_vals[, 1L, drop = FALSE]
  diff_cols <- fitted_vals[, -1L, drop = FALSE] - ref_col[, rep(1L, n_groups - 1L)]
  colnames(diff_cols) <- sprintf("f%d-f%d", grp_levels[-1L], grp_levels[1L])
  
  # Standard errors of differences via delta method ---------------------------
  diff_se <- compute_diff_standard_errors(
    coefx        = coef_vec,
    cov_betas    = vcov_mat,
    groups       = groups,
    group_name   = group_name,
    xtransformed = x_split,
    group_fp_powers = group_fp_powers
  )
  colnames(diff_se) <- sprintf("se(f%d-f%d)", grp_levels[-1L], grp_levels[1L])
  
  # 95% CI for differences ----------------------------------------------------
  diff_lower <- diff_cols - crit_val * diff_se
  colnames(diff_lower) <- sprintf("(f%d-f%d)_lower", grp_levels[-1L], grp_levels[1L])
  
  diff_upper <- diff_cols + crit_val * diff_se
  colnames(diff_upper) <- sprintf("(f%d-f%d)_upper", grp_levels[-1L], grp_levels[1L])
  
  # Assemble and return --------------------------------------------------------
  cbind(
    cont_var,
    fitted_vals,
    fitted_se,
    fitted_lower,
    fitted_upper,
    diff_cols,
    diff_se,
    diff_lower,
    diff_upper
  )
}