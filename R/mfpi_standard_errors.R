# Delta-method standard errors for MFPI function differences
#
# These three functions work together to compute pointwise standard errors for
# the difference fj(x) - f0(x) between a non-reference group j and the
# reference group 0, using the delta method applied to the fitted interaction
# model coefficients.
#
# Call chain (external -> internal):
#   compute_fitted_standard_errors()      -- iterates over non-reference groups
#     +-- compute_fitted_standard_errors() -- applies delta method for one group pair
#           +-- compute_group_difference() -- evaluates fj(x) - f0(x) given coefs
#
# Naming conventions match the rest of the package:
#   group_fp_powers  -- named list of per-group FP powers (was `powers`)
#   group_name       -- name of the categorical grouping variable
#   xtransformed     -- FP-transformed predictor matrix from
#                       gen_fitted_values_per_group()


# -----------------------------------------------------------------------------
# compute_group_difference() --------------------------------------------------
# -----------------------------------------------------------------------------

#' Evaluate the Pointwise Group-Difference Function \eqn{f_j(x) - f_0(x)}
#'
#' Computes the pointwise difference between the fitted linear predictor for a
#' non-reference group \eqn{j} and the reference group \eqn{0} at each
#' observed value of the continuous variable. This function is passed to
#' [numDeriv::jacobian()] so that the partial derivatives needed for the delta
#' method can be obtained numerically.
#'
#' @section Interaction model structure:
#' The fitted interaction model has the form
#' \deqn{
#'   \eta(x, g) =
#'     \alpha
#'     + \boldsymbol{\beta}_{g}^{\top} \phi(x)
#'     + \gamma_{g} \, \mathbf{1}[g > 0],
#' }
#' where
#' \itemize{
#'   \item \eqn{\alpha} is the intercept (zero for Cox models),
#'   \item \eqn{\phi(x) = \bigl(\phi_1(x), \ldots, \phi_m(x)\bigr)^{\top}}
#'     collects the \eqn{m} FP basis functions evaluated at \eqn{x}
#'     (\eqn{m = 1} for FP1, \eqn{m = 2} for FP2),
#'   \item \eqn{\boldsymbol{\beta}_{g} \in \mathbb{R}^{m}} are the FP
#'     regression coefficients specific to group \eqn{g},
#'   \item \eqn{\gamma_g} is the main-effect offset for group \eqn{g}
#'     (\eqn{\gamma_0 = 0} by convention, i.e. group 0 is the reference).
#' }
#'
#' @section Group-difference function:
#' Because the intercept \eqn{\alpha} is shared across groups it cancels in
#' the difference, giving
#' \deqn{
#'   D_j(x)
#'   \;=\; f_j(x) - f_0(x)
#'   \;=\;
#'   \boldsymbol{\beta}_{j}^{\top} \phi(x) + \gamma_j
#'   \;-\;
#'   \boldsymbol{\beta}_{0}^{\top} \phi(x).
#' }
#' Note that \eqn{\phi(x)} is the same vector for both groups; the
#' group-specific FP-transformed columns \eqn{\mathbf{x}_0} and
#' \eqn{\mathbf{x}_j} in `xtransformed` contain \eqn{\phi(x)} evaluated once
#' per group so that [numDeriv::jacobian()] can differentiate with respect to
#' \eqn{\boldsymbol{\beta}_0} and \eqn{\boldsymbol{\beta}_j} independently.
#'
#' @section Coefficient layout of `coef_vec`:
#' `coef_vec` is a concatenation of three blocks:
#' \deqn{
#'   \mathbf{c} =
#'   \bigl(
#'     \underbrace{\beta_{0,1}, \ldots, \beta_{0,m}}_{\text{group 0 FP}},\;
#'     \underbrace{\beta_{j,1}, \ldots, \beta_{j,m}}_{\text{group } j \text{ FP}},\;
#'     \underbrace{\gamma_j}_{\text{group } j \text{ offset}}
#'   \bigr)^{\top},
#' }
#' so `length(coef_vec)` equals \eqn{2m + 1}.  The corresponding predictor
#' matrix `xtransformed` has \eqn{2m} columns (the last element of `coef_vec`
#' multiplies an implicit constant 1 and has no column in `xtransformed`).
#'
#' @param coef_vec Numeric vector of length \eqn{2m + 1} arranged as described
#'   above: \eqn{m} group-0 FP coefficients, \eqn{m} group-\eqn{j} FP
#'   coefficients, and the scalar group-\eqn{j} offset \eqn{\gamma_j}.
#' @param xtransformed Numeric matrix with \eqn{n} rows and \eqn{2m} columns.
#'   The first \eqn{m} columns hold \eqn{\phi(x_i)} for group 0; the next
#'   \eqn{m} columns hold \eqn{\phi(x_i)} for group \eqn{j}.
#' @param group_fp_powers Named list of FP power vectors, one per group. Only
#'   the length of the first element is used here to determine \eqn{m}.
#'
#' @return Numeric vector of length \eqn{n} giving
#'   \eqn{D_j(x_i) = f_j(x_i) - f_0(x_i)} for each observation \eqn{i}.
#'
#' @keywords internal
#' @noRd
compute_group_difference <- function(coef_vec, xtransformed, group_fp_powers) {
  
  n_fp_terms <- length(group_fp_powers[[1L]])   # m: 1 for FP1, 2 for FP2
  
  # Partition coef_vec into [beta_0 | beta_j | gamma_j]
  idx_base  <- seq_len(n_fp_terms)
  idx_grpj  <- seq_len(n_fp_terms) + n_fp_terms
  idx_dummy <- length(coef_vec)                 # last element
  
  beta_0  <- coef_vec[idx_base]
  beta_j  <- coef_vec[idx_grpj]
  gamma_j <- coef_vec[idx_dummy]
  
  x_base <- xtransformed[, idx_base,  drop = FALSE]
  x_grpj <- xtransformed[, idx_grpj,  drop = FALSE]
  
  # D_j(x) = beta_j' phi(x) + gamma_j - beta_0' phi(x)
  as.vector(x_grpj %*% beta_j) + gamma_j - as.vector(x_base %*% beta_0)
}


# -----------------------------------------------------------------------------
# compute_diff_standard_errors() ----------------------------------------------
# -----------------------------------------------------------------------------

#' Compute Delta-Method Standard Errors for All Group Differences
#'
#' Iterates over the \eqn{K - 1} non-reference groups, assembles the relevant
#' coefficient sub-vector and covariance sub-matrix for each group pair
#' \eqn{(0, j)}, and delegates to \code{compute_standard_errors()} to obtain
#' pointwise standard errors for \eqn{\hat{D}_j(x) = \hat{f}_j(x) -
#' \hat{f}_0(x)}.
#'
#' @section Delta-method overview:
#' Let \eqn{\hat{\boldsymbol{\theta}}_j} denote the sub-vector of estimated
#' coefficients that contribute to \eqn{D_j(x)}, and let
#' \eqn{\hat{V}_j = \widehat{\mathrm{Cov}}(\hat{\boldsymbol{\theta}}_j)} be
#' the corresponding covariance sub-matrix extracted from the full model
#' covariance. By the delta method, the pointwise variance of
#' \eqn{\hat{D}_j(x_i)} is
#' \deqn{
#'   \widehat{\mathrm{Var}}\!\bigl[\hat{D}_j(x_i)\bigr]
#'   \;=\;
#'   \mathbf{g}_i^{\top}\, \hat{V}_j\, \mathbf{g}_i,
#' }
#' where
#' \deqn{
#'   \mathbf{g}_i
#'   = \frac{\partial D_j(x_i)}{\partial \boldsymbol{\theta}_j}
#'   \bigg|_{\boldsymbol{\theta}_j = \hat{\boldsymbol{\theta}}_j}
#' }
#' is the gradient vector of \code{group_diff_function()} with respect to
#' \eqn{\boldsymbol{\theta}_j}, evaluated numerically by
#' [numDeriv::jacobian()].
#'
#' @section Coefficient sub-vector \eqn{\hat{\boldsymbol{\theta}}_j}:
#' For the \eqn{j}-th non-reference group, the sub-vector is
#' \deqn{
#'   \hat{\boldsymbol{\theta}}_j =
#'   \bigl(
#'     \hat{\boldsymbol{\beta}}_0^{\top},\;
#'     \hat{\boldsymbol{\beta}}_j^{\top},\;
#'     \hat{\gamma}_j
#'   \bigr)^{\top} \in \mathbb{R}^{2m+1},
#' }
#' where \eqn{m} is the number of FP terms per group. The intercept
#' \eqn{\alpha} is excluded because it cancels in \eqn{D_j}.
#'
#' @param coefx Named numeric vector of all regression coefficients from the
#'   fitted interaction model, including the intercept and group dummies.
#' @param cov_betas Named square matrix. Full estimated covariance matrix of
#'   `coefx`, as returned by [vcov()].
#' @param groups List of character vectors of length \eqn{K}, one per group.
#'   Each element gives the coefficient names of the FP terms belonging to
#'   that group, in ascending group order (reference group first). As returned
#'   by \code{var_group()}.
#' @param group_name Character string. Column name of the grouping variable in
#'   the original data. Used to identify the group dummy coefficients
#'   \eqn{\hat{\gamma}_j} in `coefx` via a regex pattern
#'   `"^<group_name>\d+"`.
#' @param xtransformed Numeric matrix of FP-transformed predictor columns for
#'   all groups, as produced by \code{gen_fitted_values_per_group()}. Columns for
#'   every group are present side by side in the same order as `groups`.
#' @param group_fp_powers Named list of FP power vectors, one per group.
#'   Passed through to \code{group_diff_function()} to determine \eqn{m}.
#'
#' @return Numeric matrix with \eqn{n} rows and \eqn{K - 1} columns. Column
#'   \eqn{j} contains
#'   \deqn{
#'     \widehat{\mathrm{SE}}\!\bigl[\hat{D}_j(x_i)\bigr]
#'     = \sqrt{\mathbf{g}_i^{\top}\, \hat{V}_j\, \mathbf{g}_i}
#'   }
#'   for \eqn{i = 1, \ldots, n}.
#'
#' @seealso \code{compute_standard_errors()}, \code{group_diff_function()},
#'   \code{gen_fitted_values_per_group()}
#'
#' @keywords internal
#' @noRd
compute_diff_standard_errors <- function(coefx, cov_betas, groups,
                                    group_name, xtransformed, group_fp_powers) {
  
  coef_names <- names(coefx)
  
  # Positional index of each group's FP coefficients in coefx
  coef_indices <- lapply(groups, function(grp) which(coef_names %in% grp))
  
  # Names of the group dummy coefficients gamma_j (one per non-reference group)
  dummy_names <- grep(paste0("^", group_name, "\\d+"), coef_names, value = TRUE)
  n_dummies   <- length(dummy_names)
  
  if (n_dummies != length(groups) - 1L) {
    stop(
      paste0(
        "Number of group dummy coefficients (", n_dummies, ") does not match ",
        "the expected number of non-reference groups (", length(groups) - 1L, ")."
      ),
      call. = FALSE
    )
  }
  
  # FP coefficients beta_0 for the reference group
  coef_base <- coefx[coef_indices[[1L]]]
  
  se_list <- vector("list", length = n_dummies)
  
  for (i in seq_len(n_dummies)) {
    
    # Assemble theta_j = (beta_0, beta_j, gamma_j)
    coef_grpj  <- coefx[coef_indices[[i + 1L]]]
    coef_dummy <- coefx[dummy_names[i]]
    coef_final <- c(coef_base, coef_grpj, coef_dummy)
    
    # Extract V_j: covariance sub-matrix for theta_j
    cov_final <- cov_betas[names(coef_final), names(coef_final), drop = FALSE]
    
    # Predictor columns for groups 0 and j only (2m columns)
    x_sub <- xtransformed[, c(names(coef_base), names(coef_grpj)), drop = FALSE]
    
    se_list[[i]] <- compute_fitted_standard_errors(
      coef_final      = coef_final,
      cov_final       = cov_final,
      xtransformed    = x_sub,
      group_fp_powers = group_fp_powers
    )
  }
  
  do.call(cbind, se_list)
}


# -----------------------------------------------------------------------------
# compute_fitted_standard_errors() ---------------------------------------------
# -----------------------------------------------------------------------------

#' Apply the Delta Method to Compute Pointwise Standard Errors for One Group Pair
#'
#' Given the coefficient sub-vector \eqn{\hat{\boldsymbol{\theta}}_j} and its
#' covariance \eqn{\hat{V}_j}, evaluates the Jacobian of
#' \code{group_diff_function()} numerically and returns the pointwise standard error
#' of \eqn{\hat{D}_j(x) = \hat{f}_j(x) - \hat{f}_0(x)}.
#'
#' @section Delta-method computation:
#' Define the \eqn{n \times (2m+1)} Jacobian matrix
#' \deqn{
#'   J \;=\; \frac{\partial \mathbf{D}_j}{\partial \boldsymbol{\theta}_j^{\top}}
#'   \bigg|_{\boldsymbol{\theta}_j = \hat{\boldsymbol{\theta}}_j},
#'   \qquad
#'   J_{ik} = \frac{\partial D_j(x_i)}{\partial \theta_{j,k}},
#' }
#' where \eqn{\mathbf{D}_j = \bigl(D_j(x_1), \ldots, D_j(x_n)\bigr)^{\top}}
#' is the vector of pointwise differences. Because \code{group_diff_function()} is
#' linear in \eqn{\boldsymbol{\theta}_j} the Jacobian rows have a closed form,
#' but it is evaluated numerically via [numDeriv::jacobian()] for generality.
#'
#' The pointwise variance and standard error are then
#' \deqn{
#'   \widehat{\mathrm{Var}}\!\bigl[\hat{D}_j(x_i)\bigr]
#'   \;=\; J_i\, \hat{V}_j\, J_i^{\top},
#'   \qquad
#'   \widehat{\mathrm{SE}}\!\bigl[\hat{D}_j(x_i)\bigr]
#'   \;=\; \sqrt{J_i\, \hat{V}_j\, J_i^{\top}},
#' }
#' where \eqn{J_i} is the \eqn{i}-th row of \eqn{J}. To avoid allocating the
#' full \eqn{n \times n} matrix \eqn{J \hat{V}_j J^{\top}}, the quadratic form
#' is evaluated row-by-row.
#'
#' @section Jacobian structure for FP models:
#' For an FP\eqn{m} model the group-difference function is linear in
#' \eqn{\boldsymbol{\theta}_j}, so the Jacobian has the closed-form rows
#' \deqn{
#'   J_i =
#'   \Bigl(
#'     -\phi_1(x_i),\, \ldots,\, -\phi_m(x_i),\;
#'      \phi_1(x_i),\, \ldots,\,  \phi_m(x_i),\;
#'      1
#'   \Bigr),
#' }
#' corresponding to derivatives with respect to
#' \eqn{(\boldsymbol{\beta}_0, \boldsymbol{\beta}_j, \gamma_j)}. The numerical
#' Jacobian from [numDeriv::jacobian()] recovers this exactly.
#'
#' @param coef_final Named numeric vector \eqn{\hat{\boldsymbol{\theta}}_j} of
#'   length \eqn{2m + 1}, arranged as group-0 FP coefficients
#'   \eqn{\hat{\boldsymbol{\beta}}_0}, then group-\eqn{j} FP coefficients
#'   \eqn{\hat{\boldsymbol{\beta}}_j}, then the group-\eqn{j} offset
#'   \eqn{\hat{\gamma}_j}. See \code{group_diff_function()} for the full layout.
#' @param cov_final Named square matrix \eqn{\hat{V}_j} of dimension
#'   \eqn{(2m+1) \times (2m+1)}. The estimated covariance sub-matrix for
#'   exactly the coefficients in `coef_final`, extracted from the full model
#'   covariance by \code{compute_std_errors_diff()}.
#' @param xtransformed Numeric matrix with \eqn{n} rows and \eqn{2m} columns.
#'   The first \eqn{m} columns contain \eqn{\phi(x_i)} for group 0; the next
#'   \eqn{m} columns contain \eqn{\phi(x_i)} for group \eqn{j}.
#' @param group_fp_powers Named list of FP power vectors passed to
#'   \code{group_diff_function()} to determine \eqn{m}.
#'
#' @return Numeric vector of length \eqn{n}:
#'   \deqn{
#'     \Bigl(
#'       \widehat{\mathrm{SE}}\!\bigl[\hat{D}_j(x_1)\bigr],\;
#'       \ldots,\;
#'       \widehat{\mathrm{SE}}\!\bigl[\hat{D}_j(x_n)\bigr]
#'     \Bigr)^{\top}.
#'   }
#'
#' @seealso [numDeriv::jacobian()], \code{compute_std_errors_diff()}
#'
#' @keywords internal
#' @noRd
compute_fitted_standard_errors <- function(coef_final, cov_final,
                                    xtransformed, group_fp_powers) {
  
  # Fix xtransformed and group_fp_powers; let numDeriv vary only coef_final
  diff_fn <- function(coef_vec) {
    compute_group_difference(
      coef_vec        = coef_vec,
      xtransformed    = xtransformed,
      group_fp_powers = group_fp_powers
    )
  }
  
  # J: n x (2m+1) Jacobian matrix  dD_j/d(theta_j)
  jac <- numDeriv::jacobian(func = diff_fn, x = coef_final)
  
  # SE_i = sqrt(J_i V_j J_i'), computed row-by-row to avoid the n x n matrix
  sqrt(vapply(
    seq_len(nrow(jac)),
    function(k) {
      Ji <- jac[k, , drop = FALSE]
      as.numeric(Ji %*% cov_final %*% t(Ji))
    },
    numeric(1L)
  ))
}