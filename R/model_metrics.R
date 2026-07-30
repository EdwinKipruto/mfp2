#' Deviance computation for Gaussian models as used by Stata's `mfp`
#'
#' @description
#' Internal helper for computing the Gaussian-model deviance used in
#' fractional-polynomial model comparisons.
#'
#' @param residuals Numeric vector of model residuals.
#' @param weights Numeric vector of prior/case weights corresponding to
#'   `residuals`.
#'
#' @details
#' This is not the usual R Gaussian deviance. It follows the normal-error
#' deviance formula used by Stata's `fp`/`mfp` implementation:
#'
#' `D = n * (1 - l + log((2 * pi * rss) / n))`
#'
#' where `rss` is the weighted residual sum of squares:
#'
#' `rss = sum(weights * residuals^2)`
#'
#' and `l` is the mean of the log-normalized weights:
#'
#' `l = mean(log(weights / mean(weights)))`
#'
#' If all weights are equal, this term is zero because
#' `weights / mean(weights) = 1` and `log(1) = 0`.
#'
#' Observations with non-finite or missing residuals or weights are excluded
#' consistently from `rss`, `n`, and the weight-normalization term.
#'
#' See the Stata fractional-polynomial reference manual:
#' <https://www.stata.com/manuals/rfp.pdf>
#'
#' This deviance is intended for normal-error models only and should not be used
#' for other GLM families.
#'
#' @return A numeric value representing the Stata-style Gaussian deviance, or
#'   `NULL` if residuals or weights are unavailable.
#'
#' @references
#' StataCorp. `fp — Fractional polynomial regression`.
#' <https://www.stata.com/manuals/rfp.pdf>
#'
#' @keywords internal
#' @noRd
deviance_gaussian <- function(residuals, weights) {

  if (is.null(residuals) || is.null(weights)) {
    return(NULL)
  }

  if (length(residuals) != length(weights)) {
    stop(
      "`residuals` and `weights` must have the same length.",
      call. = FALSE
    )
  }

  ok <- !is.na(residuals) & !is.na(weights) &
    is.finite(residuals) & is.finite(weights)

  residuals <- residuals[ok]
  weights <- weights[ok]

  if (!length(residuals)) {
    return(NULL)
  }

  if (any(weights < 0)) {
    stop("`weights` must be non-negative.", call. = FALSE)
  }

  if (!any(weights > 0)) {
    return(NULL)
  }

  # Zero-weight observations do not contribute to the weighted RSS and cannot
  # enter log(weights / mean(weights)); exclude them consistently.
  positive <- weights > 0
  residuals <- residuals[positive]
  weights <- weights[positive]

  rss <- sum(weights * residuals^2)
  n <- length(residuals)

  if (!is.finite(rss) || rss <= 0) {
    return(NULL)
  }

  meanwts <- mean(log(weights / mean(weights)))
  k <- log((2 * pi * rss) / n)

  n * (1 - meanwts + k)
}



#' Function to compute model metrics to be used within `mfp2`
#'
#' Mostly used within an mfp step to compare between the different fp models
#' of a variable.
#'
#' @param obj a list returned by \code{fit_model()} representing a glm or Cox model
#' fit.
#' @param n_obs a numeric value indicating the number of observations for the
#' data used to fit `obj`.
#' @param df_additional a numeric value indicating the number of additional
#' degrees of freedom to be accounted for in the computations of AIC and BIC.
#' These may be necessary when a model uses FP terms, as these add another
#' degree of freedom per estimated power.
#'
#' @return
#' A numeric vector with the following entries:
#'
#' * `df`: number of degrees of freedom of model (i.e. coefficients plus
#' `df_additional`).
#' * `deviance_rs`: "deviance", i.e. minus twice the log likelihood.
#' This is not the usual definition of deviance used by R, which is defined as
#' twice the difference between the log likelihoods of the saturated model
#' (one parameter per observation) and the null (or reduced) model.
#' It is, however, the definition used in Royston and Sauerbrei (2008) and in
#' `mfp`. For selection of fps this does not really play a role, as the common
#' factor would be cancelled anyway when comparing models based on deviances.
#' * `deviance_gaussian`: deviance computed by \code{deviance_gaussian()},
#' applicable to Gaussian models and used for F-test computations.
#' * `aic`: Akaike information criterion, defined as
#' `-2logL + 2(df + df_additional)`.
#' * `bic`: Bayesian information criterion, defined as
#' `-2logL + log(n_obs)(df + df_additional)`.
#' * `df_resid`: residual degrees of freedom. Gaussian scale and negative-
#' binomial theta are excluded, while regression and FP transformation degrees
#' of freedom are subtracted from `n_obs`.
#'
#' @references
#' Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building:
#' A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials
#' for Modelling Continuous Variables. John Wiley & Sons.}\cr
#' @keywords internal
#' @noRd
calculate_model_metrics <- function(obj,
                                    n_obs,
                                    df_additional = 0) {
  # Collect the core model metrics returned by fit_model().
  #
  # obj$df is the model degrees of freedom reported by fit_model().
  # For Gaussian models, fit_model() includes the estimated scale parameter
  # in obj$df, and negative-binomial models include theta. Their nuisance-df
  # contribution is derived below from obj$df and the stored regression rank.
  #
  # df_additional is used to add extra degrees of freedom for FP/ACD
  # transformations when the fitted model has more transformation parameters
  # than ordinary linear terms.
  res <- c(
    logl = obj$logl,
    df = obj$df + df_additional,
    deviance_rs = -2 * obj$logl
  )

  # Residual df should subtract only regression and transformation
  # parameters. Gaussian scale and negative-binomial theta are nuisance
  # parameters included in obj$df but not in the regression residual df.
  # Derive their count from the stored regression rank rather than carrying
  # family-specific boolean flags in every fit_model() result.
  fit_rank <- obj$rank
  regression_df <- if (is.numeric(fit_rank) && length(fit_rank) == 1L &&
                       !is.na(fit_rank) && is.finite(fit_rank)) {
    unname(fit_rank)
  } else {
    sum(!is.na(obj$coefficients))
  }

  nuisance_df <- max(0, obj$df - regression_df)
  df_for_resid <- res[["df"]] - nuisance_df

  # fit_model() computes this scalar only when the caller requests Gaussian
  # F-test support. LR, AIC, and BIC paths deliberately leave it unavailable
  # so they avoid unnecessary O(n) work and residual/weight-vector copies.
  gaussian_deviance <- obj$deviance_gaussian
  if (is.null(gaussian_deviance) || length(gaussian_deviance) == 0L) {
    gaussian_deviance <- NA_real_
  } else if (!is.numeric(gaussian_deviance) ||
             length(gaussian_deviance) != 1L) {
    stop(
      "Internal error: `deviance_gaussian` must be a numeric scalar.",
      call. = FALSE
    )
  } else {
    gaussian_deviance <- unname(gaussian_deviance)
  }

  c(res,
    deviance_gaussian = gaussian_deviance,
    aic = res[["deviance_rs"]] + 2 * res[["df"]],
    bic = res[["deviance_rs"]] + log(n_obs) * res[["df"]],
    df_resid = n_obs - df_for_resid
  )
}