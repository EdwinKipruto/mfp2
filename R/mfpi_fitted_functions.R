# Fitted-function computation for MFPI
#
# gen_fitted_values_per_group() computes group-specific FP fitted functions,
# their pointwise standard errors, and 95% confidence intervals from a fitted
# interaction model. It is called by flex0-flex4 whenever compute_fitted = TRUE
# and is intended to feed the plotting layer and, eventually, predict.mfpi().
#
# The `data_input` argument (original vs. equidistant grid) belongs in
# predict.mfpi() once that method is implemented; it is retained here for
# backward compatibility but marked clearly in the documentation.
#
# Naming conventions match the rest of the package:
#   cont_var          - the continuous variable matrix
#   group_var         - the grouping variable matrix
#   group_fp_powers   - named list of per-group FP powers (was `powers`)
#   interaction_model - fitted model object from test_interaction()


# -----------------------------------------------------------------------------
# gen_fitted_values_per_group() -----------------------------------------------
# -----------------------------------------------------------------------------

#' Compute Group-Specific Fitted FP Functions and Confidence Intervals
#'
#' Given a fitted interaction model and the estimated FP powers for each group,
#' computes the fitted linear predictor \eqn{\hat{f}_j(x)} for each level
#' \eqn{j} of \code{group_var}, together with pointwise standard errors and
#' 95\% confidence intervals. The function-difference
#' \eqn{\hat{f}_j - \hat{f}_0} (relative to the reference group) and its
#' standard error are also returned.
#'
#' Standard errors for the differences are computed via the delta method:
#' \deqn{\text{Var}(\hat{f}_j - \hat{f}_0) =
#'   \mathbf{d}^\top \, \text{Cov}(\hat{\boldsymbol{\beta}}) \, \mathbf{d},}
#' where \eqn{\mathbf{d}} is the vector of partial derivatives of the
#' difference with respect to the model coefficients.
#'
#' @note The \code{use_grid} argument (equidistant evaluation grid) is a
#'   temporary placeholder. Once \code{predict.mfpi()} is implemented,
#'   grid-based prediction should be handled there and this argument will
#'   be removed.
#'
#' @param cont_var A one-column numeric matrix of the continuous variable.
#'   Must have a column name. Should already be shifted and scaled as it
#'   was at model fitting time (i.e. the same \code{x[, cont_var]} matrix
#'   passed to the flex functions).
#' @param group_fp_powers Named list of FP power vectors, one element per
#'   group. Names are the group-specific column names following the
#'   \code{<varname><group>1} convention (e.g. \code{cavol01}, \code{cavol11}
#'   for variable \code{cavol} with groups 0 and 1). Each element is a
#'   numeric vector of FP powers for that group. For flex0-flex3, all groups
#'   share the same powers; for flex4, each group may have its own powers.
#' @param interaction_model The fitted interaction model object returned by
#'   \code{test_interaction()}. Must expose \code{$coefficients} and support
#'   \code{vcov()}.
#' @param group_var A one-column numeric matrix of the grouping variable.
#'   Must have a column name.
#' @param family The regression family object used to fit the interaction
#'   model (a \code{glm} family object for non-Cox families, or the
#'   character string \code{"cox"} for Cox models). Used internally by the
#'   fitting routines.
#' @param family_string Character string identifying the family:
#'   \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, or
#'   \code{"cox"}. Used to determine whether to include an intercept
#'   (Cox models have no intercept; the intercept is set to zero).
#' @param transform Logical. Whether to FP-transform \code{cont_var} using
#'   \code{group_fp_powers} before computing fitted values. Default
#'   \code{TRUE}. Set to \code{FALSE} only if \code{cont_var} has already
#'   been transformed.
#' @param center Logical. Whether to mean-centre the FP-transformed
#'   variables before computing fitted values. Must match the \code{center}
#'   setting used when fitting the interaction model. Default \code{FALSE}.
#' @param center_vals Named numeric vector of centering constants, one per
#'   column of the FP-transformed evaluation matrix. When supplied (as
#'   returned by the flex functions via \code{create_z_variables()$center_vals}
#'   or the flex4 per-group centering block), these exact fit-time constants
#'   are subtracted from \code{x_split}. This is the recommended path: it
#'   guarantees that \eqn{\hat{f}_0(x)} and \eqn{\hat{f}_j(x)} are on the
#'   same scale as the fitted model coefficients regardless of whether a grid
#'   is used. If \code{NULL} and \code{center = TRUE}, falls back to
#'   recomputing from \code{x_split} using \code{center_type}; note that
#'   the fallback may diverge from fit-time centering when
#'   \code{use_grid = TRUE} because the grid replaces the original
#'   observations.
#' @param center_type Character string; \code{"grand"} (default) or
#'   \code{"group"}. Used only in the fallback path when
#'   \code{center_vals = NULL}. \code{"grand"} subtracts the grand mean
#'   of \code{x_split} (all observations); \code{"group"} subtracts the
#'   within-group mean using the group membership mask.
#' @param use_grid Logical. If \code{TRUE}, replaces \code{cont_var} with a
#'   200-point equidistant sequence spanning its observed range before
#'   computing fitted values, producing smooth curves for plotting. Default
#'   \code{FALSE}. \strong{This argument is provisional and will move to
#'   \code{predict.mfpi()} in a future release.} Note that when
#'   \code{use_grid = TRUE} and \code{center_vals = NULL}, the fallback
#'   centering is computed on the grid rather than the original sample,
#'   which may differ from fit-time centering; supply \code{center_vals}
#'   to avoid this.
#'
#' @return A numeric matrix with one row per observation (or per grid point
#'   when \code{use_grid = TRUE}) and the following columns, where \eqn{j}
#'   ranges over the non-reference group levels and group 0 is the reference:
#' \describe{
#'   \item{\code{<varname>}}{The (possibly grid-replaced) values of
#'     \code{cont_var} as supplied (shifted/scaled).}
#'   \item{\code{f0}, \code{f1}, \ldots}{Fitted linear predictor for each
#'     group level.}
#'   \item{\code{se(f0)}, \code{se(f1)}, \ldots}{Pointwise standard errors
#'     of the fitted values.}
#'   \item{\code{f0_lower}, \code{f0_upper}, \ldots}{Lower and upper bounds
#'     of the 95\% confidence interval for each group's fitted curve.}
#'   \item{\code{f1-f0}, \code{f2-f0}, \ldots}{Pointwise difference of each
#'     non-reference group's curve relative to group 0.}
#'   \item{\code{se(f1-f0)}, \ldots}{Standard errors of the differences,
#'     computed via the delta method.}
#'   \item{\code{(f1-f0)_lower}, \code{(f1-f0)_upper}, \ldots}{95\% CI
#'     bounds for each difference curve.}
#' }
#'   Three attributes are attached to the returned matrix for use in
#'   \code{predict.mfpi()}: \code{fp_centers} (the centering constants
#'   applied, or \code{NULL} if \code{center = FALSE}); \code{center_type}
#'   (the centering strategy used, or \code{NULL} if \code{center = FALSE});
#'   and \code{group_fp_powers} (the \code{group_fp_powers} list as passed in).
#'
#' @seealso \code{create_z_variables()}, \code{compute_diff_standard_errors()},
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
                                        transform   = TRUE,
                                        center      = FALSE,
                                        center_vals = NULL,
                                        center_type = c("grand", "group"),
                                        use_grid    = FALSE) {
  
  center_type <- match.arg(center_type)
  
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
  
  # Cox models have no intercept; all others include one.
  # Use an explicit flag rather than extracting the intercept value, so that
  # downstream checks are not affected by named-vector comparison issues.
  has_intercept <- family_string != "cox" && "(Intercept)" %in% names(coef_vec)
  intercept     <- if (has_intercept) unname(coef_vec["(Intercept)"]) else 0
  
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
  
  # Initialise col_means so it is always defined for the attribute assignment
  col_means <- NULL
  
  # ---------------------------------------------------------------------------
  # Centering of x_split
  # ---------------------------------------------------------------------------
  # When center_vals is supplied (from fit time via create_z_variables), reuse
  # the exact constants so that f0(x) and f1(x) are on the same scale as the
  # fitted model coefficients.
  #
  # Fallback (center_vals = NULL): recompute from x_split using center_type.
  # x_split has n_obs rows and n_groups * n_fp columns. The group membership
  # mask is used for "group" centering -- NOT a zero-value test -- to correctly
  # handle valid FP values of zero (e.g. log(1) = 0).
  if (center) {
    if (!is.null(center_vals)) {
      # Primary path: reuse exact fit-time constants
      col_means <- center_vals
    } else {
      n_fp      <- ncol(x_split) %/% n_groups
      col_means <- numeric(ncol(x_split))
      names(col_means) <- colnames(x_split)
      grp_vec   <- as.vector(group_var)
      
      for (gi in seq_len(n_groups)) {
        grp_cols <- seq.int((gi - 1L) * n_fp + 1L, gi * n_fp)
        if (center_type == "grand") {
          col_means[grp_cols] <- colMeans(x_split[, grp_cols, drop = FALSE],
                                          na.rm = TRUE)
        } else {
          in_grp              <- grp_vec == grp_levels[gi]
          col_means[grp_cols] <- colMeans(
            x_split[in_grp, grp_cols, drop = FALSE], na.rm = TRUE
          )
        }
      }
    }
    x_split <- sweep(x_split, 2L, col_means, "-", check.margin = FALSE)
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
    
    if (has_intercept) {
      x_var       <- cbind(1, x_var)
      x_var_names <- c("(Intercept)", x_var_names)
    }
    
    group_coef_offset <- 0
    if (i > 1L) {
      group_coef_offset <- unname(coef_vec[group_dummy_names[i - 1L]])
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
  fitted <- cbind(
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
  
  # Attach centering metadata as attributes for reuse in predict.mfpi()
  attr(fitted, "fp_centers")      <- if (center) col_means else NULL
  attr(fitted, "center_type")     <- if (center) center_type else NULL
  attr(fitted, "group_fp_powers") <- group_fp_powers
  
  fitted
}