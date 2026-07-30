#' Construct a best-model metrics row
#'
#' Creates a one-row base data frame containing the model comparison metrics
#' for a selected MFPI interaction candidate. Fractional polynomial powers are
#' stored as list-columns to preserve their numeric representation for later
#' printing, summarising, and prediction-related use.
#'
#' @param variable Character scalar. Name of the continuous variable tested for
#'   interaction.
#' @param type Character scalar. Functional form tested for `variable`, such as
#'   `"linear"`, `"fp1"`, or `"fp2"`.
#' @param fp_powers_main Numeric vector. Fractional polynomial powers selected
#'   for the main-effect function of `variable`.
#' @param fp_powers_int Numeric vector. Fractional polynomial powers selected
#'   for the interaction function of `variable`.
#' @param deviance_int Numeric scalar. Deviance of the interaction model.
#' @param deviance_diff Numeric scalar. Difference in deviance between the
#'   main-effects model and the interaction model.
#' @param df_int Numeric scalar. Degrees of freedom associated with the
#'   interaction test.
#' @param pvalue Numeric scalar. P-value for the interaction test.
#' @param df_total Numeric scalar. Total degrees of freedom used by the selected
#'   interaction model.
#' @param AIC_main Numeric scalar. Akaike information criterion of the
#'   main-effects model.
#' @param AIC_interaction Numeric scalar. Akaike information criterion of the
#'   interaction model.
#' @param BIC_main Numeric scalar. Bayesian information criterion of the
#'   main-effects model.
#' @param BIC_interaction Numeric scalar. Bayesian information criterion of the
#'   interaction model.
#'
#' @return A one-row data frame of class `"best_model_metrics"` with columns:
#'   `type`, `variable`, `fp_powers_main`, `fp_powers_int`, `deviance_int`,
#'   `deviance_diff`, `df_int`, `pvalue`, `df_total`, `AIC_main`,
#'   `AIC_interaction`, `AIC_main_minus_int`, `BIC_main`, `BIC_interaction`,
#'   and `BIC_main_minus_int`. The `fp_powers_main` and `fp_powers_int`
#'   columns are list-columns.
#'
#' @keywords internal
#' @noRd
make_best_model_metrics <- function(variable,
                                    type = NA_character_,
                                    fp_powers_main,
                                    fp_powers_int,
                                    deviance_int,
                                    deviance_diff,
                                    df_int,
                                    pvalue,
                                    df_total,
                                    AIC_main,
                                    AIC_interaction,
                                    BIC_main,
                                    BIC_interaction) {
  out <- data.frame(
    type               = type,
    variable           = variable,
    deviance_int       = deviance_int,
    deviance_diff      = deviance_diff,
    df_int             = df_int,
    pvalue             = pvalue,
    df_total           = df_total,
    AIC_main           = AIC_main,
    AIC_interaction    = AIC_interaction,
    AIC_main_minus_int = AIC_main - AIC_interaction,
    BIC_main           = BIC_main,
    BIC_interaction    = BIC_interaction,
    BIC_main_minus_int = BIC_main - BIC_interaction,
    stringsAsFactors   = FALSE,
    check.names        = FALSE
  )

  out$fp_powers_main <- I(list(fp_powers_main))
  out$fp_powers_int  <- I(list(fp_powers_int))

  out <- out[
    c(
      "type",
      "variable",
      "fp_powers_main",
      "fp_powers_int",
      "deviance_int",
      "deviance_diff",
      "df_int",
      "pvalue",
      "df_total",
      "AIC_main",
      "AIC_interaction",
      "AIC_main_minus_int",
      "BIC_main",
      "BIC_interaction",
      "BIC_main_minus_int"
    )
  ]

  class(out) <- c("best_model_metrics", "data.frame")
  out
}


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
#'   `"negbin"`, or `"cox"`.
#' @param weights Numeric vector of observation weights, length \eqn{n}.
#' @param offset Numeric vector of linear-predictor offsets, length \eqn{n}.
#' @param ties Character string; Cox tie-handling method - `"breslow"`,
#'   `"efron"`, or `"exact"`. Ignored for non-Cox families.
#' @param strata Integer stratum vector for stratified Cox models, or `NULL`.
#' @param control Fitting control list from [stats::glm.control()] or
#'   [survival::coxph.control()].
#' @param nocenter Numeric vector for Cox centring suppression; see
#'   [survival::coxph()].
#'
#' @return A list with five components:
#' \describe{
#'   \item{`evaluation_metrics`}{A data frame of class
#'     `"best_model_metrics"` with one row per call, containing:
#'     `variable`, `fp_powers_main`, `fp_powers_int`,
#'     `deviance_int` (\eqn{-2\ell_{\text{int}}}),
#'     `deviance_diff` (\eqn{T}),
#'     `df_int` (\eqn{df_{\text{int}}}),
#'     `pvalue`,
#'     `df_total` (total df of the interaction model),
#'     `AIC_main`, `AIC_interaction`,
#'     `AIC_main_minus_int` (\eqn{\mathrm{AIC}_{\text{main}} - \mathrm{AIC}_{\text{int}}}),
#'     `BIC_main`, `BIC_interaction`,
#'     `BIC_main_minus_int` (\eqn{\mathrm{BIC}_{\text{main}} - \mathrm{BIC}_{\text{int}}}).}
#'   \item{`interaction_model`}{The fitted interaction model object returned by
#'     `fit_model()`, with `fast = FALSE` so that the full coefficient
#'     vector and covariance matrix are available.}
#'   \item{`main_model`}{The fitted main-effects (no-interaction) model object
#'     returned by `fit_model()`. Preserved so that `plot.mfpi()` can draw the
#'     constant contrast \eqn{\hat\alpha^{(M)}} implied by the no-interaction
#'     model as a horizontal reference line on difference plots.}
#'   \item{`deviance_models`}{Named list with elements `main`
#'     (\eqn{-2\ell_{\text{main}}}) and `interaction`
#'     (\eqn{-2\ell_{\text{int}}}).}
#'   \item{`df`}{Named list with elements `main` (\eqn{p_{\text{main}}}),
#'     `interaction` (\eqn{p_{\text{interaction}}}), and `interaction_only`
#'     (\eqn{df_{\text{int}}}).}
#' }
#'
#' @keywords internal
#' @noRd
test_interaction <- function(y, cont_var, group_var, xmain, xinteraction,
                             degree, bestfp_main, bestfp_interaction,
                             flex, use_ftest, family, family_string,
                             weights, offset, ties, strata, control,
                             nocenter, has_offset, fitter = "base") {

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
    fitter   = fitter,
    weights  = weights,
    offset   = offset,
    method   = ties,
    strata   = strata,
    control  = control,
    rownames = NULL,
    nocenter = nocenter,
    has_offset = has_offset,
    fast     = TRUE,
    calculate_gaussian_deviance = isTRUE(
      use_ftest && identical(family_string, "gaussian")
    ),
    keep_fit = TRUE    # retained for downstream coefficient/vcov extraction
  )

  fit_interaction <- fit_model(
    x        = xinteraction,
    y        = y,
    family   = family,
    family_string = family_string,
    fitter   = fitter,
    weights  = weights,
    offset   = offset,
    method   = ties,
    strata   = strata,
    control  = control,
    rownames = NULL,
    nocenter = nocenter,
    has_offset = has_offset,
    calculate_gaussian_deviance = isTRUE(
      use_ftest && identical(family_string, "gaussian")
    ),
    fast     = FALSE   # full fit: coefficients and vcov required downstream
  )

  # ---------------------------------------------------------------------------
  # Deviance: -2 * log-likelihood
  # ---------------------------------------------------------------------------
  dev_main        <- -2 * fit_main$logl
  dev_interaction <- -2 * fit_interaction$logl
  deviance_diff   <- dev_main - dev_interaction   # T = -2(l_main - l_int) >= 0

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
  # ---------------------------------------------------------------------------
  # P-value: F-test (Gaussian only) or chi-square LRT
  # ---------------------------------------------------------------------------
  if (deviance_diff < -sqrt(.Machine$double.eps)) {

    warning(
      sprintf(
        paste0(
          "Interaction test for variable '%s' using %s was not assigned a p-value: ",
          "the fitted interaction model has a larger deviance than the main-effects model ",
          "(main deviance = %.6g, interaction deviance = %.6g, difference = %.6g). ",
          "This can occur with non-nested selected FP transformations, rank deficiency, ",
          "or numerical convergence issues. For non-nested MFPI comparisons, especially ",
          "when using flex3 or flex4, AIC or BIC can be more appropriate than a ",
          "likelihood-ratio or F-test p-value. The interaction p-value is set to NA; ",
          "AIC/BIC values are still reported."
        ),
        cont_name,
        flex,
        dev_main,
        dev_interaction,
        deviance_diff
      ),
      call. = FALSE
    )

    pvalue <- NA_real_

  } else {

    deviance_diff <- max(deviance_diff, 0)

    if (use_ftest && family_string == "gaussian") {

      # MFP2-style Gaussian deviance, not -2 * logLik. fit_model() computes
      # and stores only the scalar when F-test support is requested.
      dev_main_f <- fit_main$deviance_gaussian
      dev_int_f <- fit_interaction$deviance_gaussian

      if (!is.finite(dev_main_f) || !is.finite(dev_int_f)) {
        stop(
          "Internal error: Gaussian deviance is unavailable for MFPI F-test.",
          call. = FALSE
        )
      }

      power_df_interaction <- switch(
        flex,
        flex0 = 0L,
        flex1 = degree,
        flex2 = degree,
        flex3 = degree,
        flex4 = n_groups * degree,
        stop("Internal error: unknown flexibility level.", call. = FALSE)
      )

      df_resid_int <- fit_interaction$fit$df.residual - power_df_interaction

      ftest_result <- calculate_f_test(
        deviances = c(dev_main_f, dev_int_f),
        dfs_resid = df_resid_int,
        n_obs     = nrow(xinteraction),
        d1        = df_int
      )

      pvalue <- ftest_result$pvalue

    } else {

      pvalue <- stats::pchisq(
        deviance_diff,
        df = df_int,
        lower.tail = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # AIC and BIC
  # AIC = -2l + 2p;  BIC = -2l + p * log(n*)
  # n* = number of events for Cox; n otherwise
  # ---------------------------------------------------------------------------
  if (family_string == "cox") {
    # `mfpi.default()` validates that Cox responses are right-censored Surv
    # objects before interaction testing is reached. For those objects, the
    # event/status indicator is stored in the final matrix column, but that
    # column is not guaranteed to be named "status". Use the final column
    # positionally rather than y[, "status"].
    status <- y[, ncol(y)]
    n_eff <- sum(!is.na(status) & status > 0)

    if (!is.finite(n_eff) || n_eff <= 0L) {
      stop(
        "Cox BIC requires at least one observed event.",
        call. = FALSE
      )
    }
  } else {
    n_eff <- NROW(y)
  }

  AIC_main           <- dev_main        + 2          * df_main
  AIC_interaction    <- dev_interaction + 2          * df_total
  BIC_main           <- dev_main        + log(n_eff) * df_main
  BIC_interaction    <- dev_interaction + log(n_eff) * df_total

  # ---------------------------------------------------------------------------
  # Assemble evaluation metrics tibble
  # ---------------------------------------------------------------------------
  #fp_powers_main <- setNames(list(bestfp_main), cont_name)
  #fp_powers_int  <- bestfp_interaction

  metric_type <- switch(
    as.character(degree),
    `0` = "linear",
    `1` = "fp1",
    `2` = "fp2",
    paste0("fp", degree)
  )

  metrics <- make_best_model_metrics(
    variable        = cont_name,
    type            = metric_type,
    fp_powers_main  = bestfp_main,
    fp_powers_int   = bestfp_interaction,
    deviance_int    = dev_interaction,
    deviance_diff   = deviance_diff,
    df_int          = df_int,
    pvalue          = pvalue,
    df_total        = df_total,
    AIC_main        = AIC_main,
    AIC_interaction = AIC_interaction,
    BIC_main        = BIC_main,
    BIC_interaction = BIC_interaction
  )

  list(
    evaluation_metrics = metrics,
    interaction_model  = fit_interaction,
    main_model         = fit_main,
    deviance_models    = list(main = dev_main, interaction = dev_interaction),
    df                 = list(
      main               = df_main,
      interaction        = df_total,
      interaction_only   = df_int
    )
  )
}

#' Print MFPI model-evaluation metrics
#'
#' Internal S3 method used to format model-evaluation results displayed by
#' [mfpi()] and its print methods.
#'
#' @param x An object of class `"best_model_metrics"`.
#' @param ... Additional arguments passed to [print.data.frame()].
#'
#' @return Invisibly returns `x`.
#'
#' @method print best_model_metrics
#' @export
#' @noRd
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