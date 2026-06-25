# Interaction test for MFPI
#
# test_interaction() fits the main-effects and interaction models, computes
# the likelihood-ratio statistic, and assembles the full table of evaluation
# metrics (deviance, df, p-value, AIC, BIC).  It is called by every flex
# function (flex0-flex4) after the design matrices have been built.
#
# print.best_model_metrics() provides a readable console representation
# of the tibble returned inside the list.
#
# Naming conventions match the rest of the package:
#   cont_var            - continuous variable matrix (was `contvar`)
#   group_var           - grouping variable matrix   (was `catvar`)
#   ties                - Cox tie-handling            (was `method`)
#   use_ftest           - F-test flag                 (was `ftest`)
#   bestfp_main         - FP powers for main model
#   bestfp_interaction  - per-group FP powers for interaction model


# -----------------------------------------------------------------------------
# test_interaction() ----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Test for an Interaction Effect Between a Continuous and a Grouping Variable
#'
#' Fits the main-effects model and the interaction model using pre-built design
#' matrices, then compares them via a likelihood-ratio test. AIC and BIC for
#' both models are also computed. The function is called by every flex
#' implementation (flex0-flex4) and is not exported.
#'
#' @section Model specification:
#' Both models share the same adjustment terms (included as columns of `xmain`
#' and `xinteraction`). Let \eqn{n} be the number of observations,
#' \eqn{K} the number of groups, and \eqn{m} the FP degree (`degree`).
#'
#' **Main-effects model** (no interaction):
#' \deqn{
#'   \eta_{\text{main}}
#'   = \alpha
#'   + \boldsymbol{\beta}^{\top} \phi(x)
#'   + \sum_{k=1}^{K-1} \gamma_k \, \mathbf{1}[g = k]
#'   + \boldsymbol{\delta}^{\top} \mathbf{z}_{\text{adj}},
#' }
#' where \eqn{\phi(x)} collects the \eqn{m} FP basis functions, \eqn{\gamma_k}
#' are group main effects, and \eqn{\mathbf{z}_{\text{adj}}} are the
#' pre-transformed adjustment covariates.
#'
#' **Interaction model** (group-specific FP slopes):
#' \deqn{
#'   \eta_{\text{int}}
#'   = \alpha
#'   + \sum_{k=0}^{K-1} \boldsymbol{\beta}_k^{\top} \phi_k(x)
#'   + \sum_{k=1}^{K-1} \gamma_k \, \mathbf{1}[g = k]
#'   + \boldsymbol{\delta}^{\top} \mathbf{z}_{\text{adj}},
#' }
#' where \eqn{\phi_k(x)} uses the FP powers specific to group \eqn{k} (which
#' may or may not differ across groups depending on `flex`).
#'
#' @section Degrees of freedom:
#' The degrees of freedom for the likelihood-ratio test are determined by
#' `interaction_model_df(n_groups, degree, flex)`:
#' \deqn{
#'   df_{\text{int}}
#'   = p_{\text{interaction}} - p_{\text{main}},
#' }
#' where \eqn{p} counts model parameters excluding the intercept and
#' adjustment terms. Specifically:
#' \itemize{
#'   \item \eqn{p_{\text{main}} = m + (K - 1)}: \eqn{m} shared FP terms
#'     plus \eqn{K-1} group dummies.
#'   \item \eqn{p_{\text{interaction}} = Km + (K - 1)}: \eqn{m} FP terms
#'     per group plus \eqn{K-1} group dummies.
#'   \item Hence \eqn{df_{\text{int}} = (K-1)m}.
#' }
#' Under `flex3` and `flex4` the main and interaction models may use different
#' FP families, so the models are non-nested and the degrees of freedom
#' calculation differs; see `interaction_model_df()` for details.
#'
#' @section Likelihood-ratio and F-test statistics:
#' Let \eqn{\ell_{\text{main}}} and \eqn{\ell_{\text{int}}} denote the
#' maximised log-likelihoods of the two models. The deviance difference is
#' \deqn{
#'   T = -2\ell_{\text{main}} - \bigl(-2\ell_{\text{int}}\bigr)
#'     = -2\bigl(\ell_{\text{main}} - \ell_{\text{int}}\bigr) \;\geq\; 0.
#' }
#' When \code{use_ftest = FALSE} (default, or for non-Gaussian families),
#' the p-value uses the chi-square approximation:
#' \deqn{p = \Pr\!\bigl[\chi^2(df_{\text{int}}) > T\bigr].}
#' When \code{use_ftest = TRUE} and \code{family = "gaussian"}, the F-test
#' of Royston and Sauerbrei (see \code{calculate_f_test()}) is used
#' instead:
#' \deqn{F = \frac{d_2}{d_1}
#'   \left(\exp\!\left(\frac{T}{n}\right) - 1\right),}
#' where \eqn{d_1 = df_{\text{resid,main}} - df_{\text{resid,int}}} is the
#' difference in residual degrees of freedom between the main-effects and
#' interaction models, and \eqn{d_2 = df_{\text{resid,int}}} is the residual
#' degrees of freedom of the interaction model. Both are taken directly from
#' the fitted model objects to correctly account for all parameters including
#' adjustment variables. The p-value is
#' \eqn{\Pr[F(d_1, d_2) > F_{\text{obs}}]}.
#'
#' @section Information criteria:
#' AIC and BIC penalise model complexity differently:
#' \deqn{
#'   \mathrm{AIC} = -2\ell + 2p, \qquad
#'   \mathrm{BIC} = -2\ell + p\log(n^*),
#' }
#' where \eqn{p} is the number of model parameters (excluding intercepts and
#' adjustment terms, as these are common across models and do not affect
#' comparisons), and \eqn{n^*} is the effective sample size
#' (\eqn{n^* =} number of events for Cox models; \eqn{n^* = n} otherwise).
#' Positive values of
#' \eqn{\mathrm{AIC}_{\text{main}} - \mathrm{AIC}_{\text{int}}} and
#' \eqn{\mathrm{BIC}_{\text{main}} - \mathrm{BIC}_{\text{int}}} indicate
#' that the interaction model fits better.
#'
#' @param y Response vector or [survival::Surv()] object. For
#'   `family = "binomial"`, a numeric vector with exactly two distinct values.
#'   For `family = "cox"`, a two-column right-censored Surv object.
#' @param cont_var A one-column numeric matrix of the continuous variable.
#'   Must have a column name.
#' @param group_var A one-column numeric matrix of the grouping variable.
#'   Must have a column name. Values must not be expanded into dummies before
#'   being passed here (that is done internally).
#' @param xmain Numeric design matrix for the main-effects model. Columns
#'   contain the FP-transformed continuous variable, group dummy variables,
#'   and pre-transformed adjustment covariates.
#' @param xinteraction Numeric design matrix for the interaction model.
#'   Columns contain the group-specific FP-transformed variables, group dummy
#'   variables, and adjustment covariates.
#' @param degree Non-negative integer. FP degree: `0` for linear, `1` for
#'   FP1, `2` for FP2.
#' @param bestfp_main Numeric vector of the FP powers used in the main-effects
#'   model for `cont_var`.
#' @param bestfp_interaction Named list of FP power vectors, one per group,
#'   used in the interaction model.
#' @param flex Character string; `"flex0"`, `"flex1"`, `"flex2"`, `"flex3"`,
#'   or `"flex4"`. Passed to `interaction_model_df()`.
#' @param use_ftest Logical. If \code{TRUE} and \code{family = "gaussian"},
#'   use an F-test rather than a chi-square likelihood-ratio test for the
#'   interaction p-value. The F-statistic uses the formula from the Stata FP
#'   manual (Royston and Sauerbrei), with residual df taken directly from the
#'   fitted interaction model to correctly account for adjustment variables.
#'   Ignored (chi-square used) for non-Gaussian families. Default \code{FALSE}.
#' @param family Character string; `"gaussian"`, `"binomial"`, `"poisson"`,
#'   or `"cox"`.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param ties Character string; Cox tie-handling method - `"breslow"`,
#'   `"efron"`, or `"exact"`. Ignored for non-Cox families.
#' @param strata Integer stratum vector for stratified Cox models, or `NULL`.
#' @param control Fitting control list from [stats::glm.control()] or
#'   [survival::coxph.control()].
#' @param nocenter Numeric vector for Cox centring suppression; see
#'   [survival::coxph()].
#' @param digits Positive integer. Number of decimal places to which metrics
#'   are rounded in the output tibble.
#'
#' @return A list with four components:
#' \describe{
#'   \item{`evaluation_metrics`}{A tibble of class
#'     `"best_model_metrics"` with one row per call, containing:
#'     `variable`, `fp_powers_main`, `fp_powers_int`,
#'     `deviance_int` (\eqn{-2\ell_{\text{int}}}),
#'     `deviance_diff` (\eqn{T}),
#'     `df_interaction` (\eqn{df_{\text{int}}}),
#'     `pvalue`,
#'     `df_total` (total df of the interaction model),
#'     `AIC_main`, `AIC_interaction`,
#'     `AIC_main_minus_int` (\eqn{\mathrm{AIC}_{\text{main}} - \mathrm{AIC}_{\text{int}}}),
#'     `BIC_main`, `BIC_interaction`,
#'     `BIC_main_minus_int` (\eqn{\mathrm{BIC}_{\text{main}} - \mathrm{BIC}_{\text{int}}}).}
#'   \item{`interaction_model`}{The fitted interaction model object returned by
#'     `fit_model()`, with `fast = FALSE` so that the full coefficient
#'     vector and covariance matrix are available.}
#'   \item{`deviance_models`}{Named list with elements `main`
#'     (\eqn{-2\ell_{\text{main}}}) and `interaction`
#'     (\eqn{-2\ell_{\text{int}}}).}
#'   \item{`df`}{Named list with elements `main` (\eqn{p_{\text{main}}}),
#'     `interaction` (\eqn{p_{\text{interaction}}}), and `interaction_only`
#'     (\eqn{df_{\text{int}}}).}
#' }
#'
#' @importFrom tibble tibble
#' @keywords internal
#' @noRd
test_interaction <- function(y, cont_var, group_var, xmain, xinteraction,
                             degree, bestfp_main, bestfp_interaction,
                             flex, use_ftest, family, family_string,
                             weights, offset, ties, strata, control,
                             nocenter, has_offset, digits) {
  
  cont_name  <- colnames(cont_var)
  n_groups   <- length(unique(as.vector(group_var)))
  
  # ---------------------------------------------------------------------------
  # Fit main-effects and interaction models
  # ---------------------------------------------------------------------------
  fit_main <- fit_model(
    x        = xmain,
    y        = y,
    family   = family,
    family_string = family_string,
    weights  = weights,
    offset   = offset,
    method   = ties,
    strata   = strata,
    control  = control,
    rownames = NULL,
    nocenter = nocenter,
    has_offset = has_offset,
    fast     = TRUE    # log-likelihood only; no vcov needed
  )
  
  fit_interaction <- fit_model(
    x        = xinteraction,
    y        = y,
    family   = family,
    family_string = family_string,
    weights  = weights,
    offset   = offset,
    method   = ties,
    strata   = strata,
    control  = control,
    rownames = NULL,
    nocenter = nocenter,
    has_offset = has_offset,
    fast     = FALSE   # full fit: coefficients and vcov required downstream
  )
  
  # ---------------------------------------------------------------------------
  # Deviance: -2 * log-likelihood
  # ---------------------------------------------------------------------------
  dev_main        <- -2 * fit_main$logl
  dev_interaction <- -2 * fit_interaction$logl
  deviance_diff        <- dev_main - dev_interaction   # T = -2(l_main - l_int) >= 0
  
  # ---------------------------------------------------------------------------
  # Degrees of freedom
  # df_int  = (K-1)*m  for flex1/flex2 (nested models)
  # differs for flex3/flex4; delegated to interaction_model_df()
  # ---------------------------------------------------------------------------
  deg_freedom <- interaction_model_df(n_groups, degree, flex)
  df_main     <- deg_freedom$dfmain       # p_main (excl. intercept & adj.)
  df_total    <- deg_freedom$total_df     # p_interaction (excl. intercept & adj.)
  df_int      <- deg_freedom$dfint        # df_total - df_main = (K-1)*m
  
  # ---------------------------------------------------------------------------
  # P-value: F-test (Gaussian only) or chi-square LRT
  # ---------------------------------------------------------------------------
  if (use_ftest && family_string == "gaussian") {
    # Use calculate_f_test with convention 1: pass both residual dfs
    # as a vector and let the function compute d1 = dfs_resid[1] - dfs_resid[2].
    # Residual dfs are taken from the fitted objects to correctly account for
    # all parameters including adjustment variables.
    df_resid_main <- fit_main$fit$df.residual
    df_resid_int  <- fit_interaction$fit$df.residual
    n_obs         <- nrow(xinteraction)
    ftest_result  <- calculate_f_test(
      deviances = c(dev_main, dev_interaction),
      dfs_resid = c(df_resid_main, df_resid_int),
      n_obs     = n_obs
    )
    pvalue <- ftest_result$pvalue
  } else {
    pvalue <- stats::pchisq(deviance_diff, df = df_int, lower.tail = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # AIC and BIC
  # AIC = -2l + 2p;  BIC = -2l + p * log(n*)
  # n* = number of events for Cox; n otherwise
  # ---------------------------------------------------------------------------
  n_eff              <- if (family_string == "cox") sum(y[, "status"]) else length(y)
  AIC_main           <- dev_main        + 2          * df_main
  AIC_interaction    <- dev_interaction + 2          * df_total
  BIC_main           <- dev_main        + log(n_eff) * df_main
  BIC_interaction    <- dev_interaction + log(n_eff) * df_total
  
  # ---------------------------------------------------------------------------
  # Assemble evaluation metrics tibble
  # ---------------------------------------------------------------------------
  fp_powers_main <- setNames(list(bestfp_main), cont_name)
  fp_powers_int  <- bestfp_interaction
  
  metrics <- tibble::tibble(
    variable           = cont_name,
    fp_powers_main     = fp_powers_main,
    fp_powers_int      = list(fp_powers_int),
    deviance_int       = round(dev_interaction,            digits),
    deviance_diff      = round(deviance_diff,              digits),
    df_interaction     = df_int,
    pvalue             = round(pvalue,                     digits),
    df_total           = df_total,
    AIC_main           = round(AIC_main,                   digits),
    AIC_interaction    = round(AIC_interaction,            digits),
    AIC_main_minus_int = round(AIC_main - AIC_interaction, digits),
    BIC_main           = round(BIC_main,                   digits),
    BIC_interaction    = round(BIC_interaction,            digits),
    BIC_main_minus_int = round(BIC_main - BIC_interaction, digits)
  )
  class(metrics) <- c("best_model_metrics", class(metrics))
  
  list(
    evaluation_metrics = metrics,
    interaction_model  = fit_interaction,
    deviance_models    = list(main = dev_main, interaction = dev_interaction),
    df                 = list(
    main               = df_main,
    interaction        = df_total,
    interaction_only   = df_int
    )
  )
}


# -----------------------------------------------------------------------------
# print.best_model_metrics() --------------------------------------------
# -----------------------------------------------------------------------------

#' Print Method for Model Evaluation Metrics
#'
#' Formats the list-columns `fp_powers_main` and `fp_powers_int` of a
#' `"best_model_metrics"` tibble as human-readable strings, then
#' delegates to [print.data.frame()].
#'
#' @section Output columns:
#' The printed table contains the following columns:
#' \describe{
#'   \item{`variable`}{Name of the continuous variable tested.}
#'   \item{`fp_powers_main`}{FP powers in the main-effects model, formatted as
#'     `"(p1)"` (FP1) or `"(p1, p2)"` (FP2).}
#'   \item{`fp_powers_int`}{FP powers in the interaction model, formatted as one
#'     parenthesised tuple per group, comma-separated.}
#'   \item{`deviance_int`}{\eqn{-2\ell_{\text{int}}}: deviance of the interaction
#'     model.}
#'   \item{`deviance_diff`}{\eqn{T = -2\ell_{\text{main}} - (-2\ell_{\text{int}})
#'     \geq 0}: likelihood-ratio test statistic.}
#'   \item{`df_interaction`}{\eqn{df_{\text{int}} = (K-1)m}: degrees of freedom for the
#'     LRT (may differ for `flex3`/`flex4`).}
#'   \item{`pvalue`}{\eqn{p = \Pr[\chi^2(df_{\text{int}}) > T]}.}
#'   \item{`df_total`}{Total parameters in the interaction model (excluding
#'     intercept and adjustment terms).}
#'   \item{`AIC_main`}{\eqn{\mathrm{AIC}_{\text{main}} =
#'     -2\ell_{\text{main}} + 2\,p_{\text{main}}}: AIC of the main-effects
#'     model (no interaction).}
#'   \item{`AIC_interaction`}{\eqn{\mathrm{AIC}_{\text{int}} =
#'     -2\ell_{\text{int}} + 2\,p_{\text{int}}}.}
#'   \item{`AIC_main_minus_int`}{\eqn{\mathrm{AIC}_{\text{main}} -
#'     \mathrm{AIC}_{\text{int}}}. Positive values favour the interaction
#'     model.}
#'   \item{`BIC_main`}{\eqn{\mathrm{BIC}_{\text{main}} =
#'     -2\ell_{\text{main}} + p_{\text{main}}\log(n^*)}: BIC of the main-effects
#'     model (no interaction).}
#'   \item{`BIC_interaction`}{\eqn{\mathrm{BIC}_{\text{int}} =
#'     -2\ell_{\text{int}} + p_{\text{int}}\log(n^*)}.}
#'   \item{`BIC_main_minus_int`}{\eqn{\mathrm{BIC}_{\text{main}} -
#'     \mathrm{BIC}_{\text{int}}}. Positive values favour the interaction
#'     model.}
#' }
#'
#' @param x An object of class `"best_model_metrics"`, as returned by
#'   \code{test_interaction()}.
#' @param ... Additional arguments passed to [print.data.frame()].
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.best_model_metrics <- function(x, ...) {
  # Guard: power columns may be absent if the user subset the data frame
  # (e.g. `x[, c("type", "pvalue")]` drops fp_powers_main / fp_powers_int).
  # Only format columns that exist and are non-empty.
  if (!is.null(x$fp_powers_main) && length(x$fp_powers_main) > 0L) {
    x$fp_powers_main <- vapply(x$fp_powers_main, function(p) {
      if (is.null(p) || length(p) == 0L) "()"
      else paste0("(", paste(p, collapse = ", "), ")")
    }, character(1L))
  }
  
  if (!is.null(x$fp_powers_int) && length(x$fp_powers_int) > 0L) {
    x$fp_powers_int <- vapply(x$fp_powers_int, function(p_list) {
      if (is.null(p_list) || length(p_list) == 0L) return("()")
      # Flex1-3: p_list is an atomic vector; flex4: list of per-group vectors
      if (is.list(p_list)) {
        paste(
          vapply(p_list, function(p) {
            if (is.null(p) || length(p) == 0L) "()"
            else paste0("(", paste(p, collapse = ", "), ")")
          }, character(1L)),
          collapse = ", "
        )
      } else {
        paste0("(", paste(p_list, collapse = ", "), ")")
      }
    }, character(1L))
  }
  
  print.data.frame(as.data.frame(x), row.names = FALSE, ...)
  invisible(x)
}