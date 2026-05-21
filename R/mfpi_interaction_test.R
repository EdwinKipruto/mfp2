# Interaction test for MFPI
#
# test_interaction() fits the main-effects and interaction models, computes
# the likelihood-ratio statistic, and assembles the full table of evaluation
# metrics (deviance, df, p-value, AIC, BIC).  It is called by every flex
# function (flex0–flex4) after the design matrices have been built.
#
# print.model_evaluation_metrics() provides a readable console representation
# of the tibble returned inside the list.
#
# Naming conventions match the rest of the package:
#   cont_var            — continuous variable matrix (was `contvar`)
#   group_var           — grouping variable matrix   (was `catvar`)
#   ties                — Cox tie-handling            (was `method`)
#   use_ftest           — F-test flag                 (was `ftest`)
#   bestfp_main         — FP powers for main model
#   bestfp_interaction  — per-group FP powers for interaction model


# -----------------------------------------------------------------------------
# test_interaction() ----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Test for an Interaction Effect Between a Continuous and a Grouping Variable
#'
#' Fits the main-effects model and the interaction model using pre-built design
#' matrices, then compares them via a likelihood-ratio test. AIC and BIC for
#' both models are also computed. The function is called by every flex
#' implementation (flex0–flex4) and is not exported.
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
#' @section Likelihood-ratio statistic:
#' Let \eqn{\ell_{\text{main}}} and \eqn{\ell_{\text{int}}} denote the
#' maximised log-likelihoods of the two models. The test statistic is
#' \deqn{
#'   T = -2\ell_{\text{main}} - \bigl(-2\ell_{\text{int}}\bigr)
#'     = -2\bigl(\ell_{\text{main}} - \ell_{\text{int}}\bigr) \;\geq\; 0,
#' }
#' which under \eqn{H_0} (no interaction) is approximately
#' \eqn{\chi^2(df_{\text{int}})}. The p-value is
#' \deqn{
#'   p = \Pr\!\bigl[\chi^2(df_{\text{int}}) > T\bigr].
#' }
#'
#' @section Information criteria:
#' AIC and BIC penalise model complexity differently:
#' \deqn{
#'   \mathrm{AIC} = -2\ell + 2p, \qquad
#'   \mathrm{BIC} = -2\ell + p\log(n^*),
#' }
#' where \eqn{p} is the total number of parameters in the model
#' (including intercept and adjustment terms) and \eqn{n^*} is the
#' effective sample size (\eqn{n^* =} number of events for Cox models,
#' \eqn{n^* = n} otherwise). Positive values of
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
#' @param use_ftest Logical. If `TRUE` and `family = "gaussian"`, use an
#'   F-test rather than a chi-square test for the interaction p-value.
#'   **Not yet implemented for the interaction test**; currently falls back to
#'   the chi-square test with a warning.
#' @param family Character string; `"gaussian"`, `"binomial"`, `"poisson"`,
#'   or `"cox"`.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param ties Character string; Cox tie-handling method — `"breslow"`,
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
#'     `"model_evaluation_metrics"` with one row per call, containing:
#'     `variable`, `pow_main`, `pow_int`, `dev_int` (\eqn{-2\ell_{\text{int}}}),
#'     `dev_diff` (\eqn{T}), `df` (\eqn{df_{\text{int}}}), `pvalue`,
#'     `tdf` (total df of the interaction model), `AIC_int`, `BIC_int`,
#'     `AIC_diff` (\eqn{\mathrm{AIC}_{\text{main}} - \mathrm{AIC}_{\text{int}}}),
#'     `BIC_diff` (\eqn{\mathrm{BIC}_{\text{main}} - \mathrm{BIC}_{\text{int}}}).}
#'   \item{`interaction_model`}{The fitted interaction model object returned by
#'     `mfp2:::fit_model()`, with `fast = FALSE` so that the full coefficient
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
                             flex, use_ftest, family, weights, offset,
                             ties, strata, control, nocenter, digits) {
  
  cont_name  <- colnames(cont_var)
  n_groups   <- length(unique(as.vector(group_var)))
  
  # ---------------------------------------------------------------------------
  # Fit main-effects and interaction models
  # ---------------------------------------------------------------------------
  fit_main <- mfp2:::fit_model(
    x        = xmain,
    y        = y,
    family   = family,
    weights  = weights,
    offset   = offset,
    method   = ties,
    strata   = strata,
    control  = control,
    rownames = NULL,
    nocenter = nocenter,
    fast     = TRUE    # log-likelihood only; no vcov needed
  )
  
  fit_interaction <- mfp2:::fit_model(
    x        = xinteraction,
    y        = y,
    family   = family,
    weights  = weights,
    offset   = offset,
    method   = ties,
    strata   = strata,
    control  = control,
    rownames = NULL,
    nocenter = nocenter,
    fast     = FALSE   # full fit: coefficients and vcov required downstream
  )
  
  # ---------------------------------------------------------------------------
  # Deviance: -2 * log-likelihood
  # ---------------------------------------------------------------------------
  dev_main        <- -2 * fit_main$logl
  dev_interaction <- -2 * fit_interaction$logl
  dev_diff        <- dev_main - dev_interaction   # T = -2(l_main - l_int) >= 0
  
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
  # P-value: chi-square LRT (F-test not yet implemented for interaction term)
  # ---------------------------------------------------------------------------
  if (use_ftest && family == "gaussian") {
    warning(
      "! F-test for the interaction term is not yet implemented; ",
      "reverting to chi-square test.",
      call. = FALSE
    )
  }
  pvalue <- stats::pchisq(dev_diff, df = df_int, lower.tail = FALSE)
  
  # ---------------------------------------------------------------------------
  # AIC and BIC
  # AIC = -2l + 2p;  BIC = -2l + p * log(n*)
  # n* = number of events for Cox; n otherwise
  # ---------------------------------------------------------------------------
  n_eff   <- if (family == "cox") sum(y[, "status"]) else length(y)
  AIC_main <- dev_main        + 2          * df_main
  AIC_int  <- dev_interaction + 2          * df_total
  BIC_main <- dev_main        + log(n_eff) * df_main
  BIC_int  <- dev_interaction + log(n_eff) * df_total
  
  # ---------------------------------------------------------------------------
  # Assemble evaluation metrics tibble
  # ---------------------------------------------------------------------------
  # Store powers as list-columns so the structure is losslessly preserved;
  # print.model_evaluation_metrics() formats them as readable strings.
  pow_main <- setNames(list(bestfp_main), cont_name)
  pow_int  <- bestfp_interaction
  
  metrics <- tibble::tibble(
    variable = cont_name,
    pow_main = pow_main,
    pow_int  = list(pow_int),
    dev_int  = round(dev_interaction, digits),
    dev_diff = round(dev_diff,        digits),
    df       = df_int,
    pvalue   = round(pvalue,          digits),
    tdf      = df_total,
    AIC_int  = round(AIC_int,         digits),
    BIC_int  = round(BIC_int,         digits),
    AIC_diff = round(AIC_main - AIC_int, digits),
    BIC_diff = round(BIC_main - BIC_int, digits)
  )
  class(metrics) <- c("model_evaluation_metrics", class(metrics))
  
  list(
    evaluation_metrics = metrics,
    interaction_model  = fit_interaction,
    deviance_models    = list(main = dev_main, interaction = dev_interaction),
    df                 = list(
      main             = df_main,
      interaction      = df_total,
      interaction_only = df_int
    )
  )
}


# -----------------------------------------------------------------------------
# print.model_evaluation_metrics() --------------------------------------------
# -----------------------------------------------------------------------------

#' Print Method for Model Evaluation Metrics
#'
#' Formats the list-columns `pow_main` and `pow_int` of a
#' `"model_evaluation_metrics"` tibble as human-readable strings, then
#' delegates to [print.data.frame()].
#'
#' @section Output columns:
#' The printed table contains the following columns:
#' \describe{
#'   \item{`variable`}{Name of the continuous variable tested.}
#'   \item{`pow_main`}{FP powers in the main-effects model, formatted as
#'     `"(p1)"` (FP1) or `"(p1, p2)"` (FP2).}
#'   \item{`pow_int`}{FP powers in the interaction model, formatted as one
#'     parenthesised tuple per group, comma-separated.}
#'   \item{`dev_int`}{\eqn{-2\ell_{\text{int}}}: deviance of the interaction
#'     model.}
#'   \item{`dev_diff`}{\eqn{T = -2\ell_{\text{main}} - (-2\ell_{\text{int}})
#'     \geq 0}: likelihood-ratio test statistic.}
#'   \item{`df`}{\eqn{df_{\text{int}} = (K-1)m}: degrees of freedom for the
#'     LRT (may differ for `flex3`/`flex4`).}
#'   \item{`pvalue`}{\eqn{p = \Pr[\chi^2(df_{\text{int}}) > T]}.}
#'   \item{`tdf`}{Total parameters in the interaction model (excluding
#'     intercept and adjustment terms).}
#'   \item{`AIC_int`}{\eqn{\mathrm{AIC}_{\text{int}} =
#'     -2\ell_{\text{int}} + 2\,p_{\text{int}}}.}
#'   \item{`BIC_int`}{\eqn{\mathrm{BIC}_{\text{int}} =
#'     -2\ell_{\text{int}} + p_{\text{int}}\log(n^*)}.}
#'   \item{`AIC_diff`}{\eqn{\mathrm{AIC}_{\text{main}} -
#'     \mathrm{AIC}_{\text{int}}}. Positive values favour the interaction
#'     model.}
#'   \item{`BIC_diff`}{\eqn{\mathrm{BIC}_{\text{main}} -
#'     \mathrm{BIC}_{\text{int}}}. Positive values favour the interaction
#'     model.}
#' }
#'
#' @param x An object of class `"model_evaluation_metrics"`, as returned by
#'   `test_interaction()`.
#' @param ... Additional arguments passed to [print.data.frame()].
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.model_evaluation_metrics <- function(x, ...) {
  x$pow_main <- vapply(x$pow_main, function(p) {
    paste0("(", paste(p, collapse = ", "), ")")
  }, character(1L))
  
  x$pow_int <- vapply(x$pow_int, function(p_list) {
    paste(
      vapply(p_list, function(p)
        paste0("(", paste(p, collapse = ", "), ")"),
        character(1L)
      ),
      collapse = ", "
    )
  }, character(1L))
  
  print.data.frame(as.data.frame(x), row.names = FALSE, ...)
  invisible(x)
}