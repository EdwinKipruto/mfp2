# Variable-ordering helpers for the MFP backfitting algorithm
#
# The full linear model has two distinct roles in mfp2:
#
#   1. It provides the null- and full-linear-model deviances retained for later
#      reporting.
#   2. When significance-based ordering is requested, it is the reference model
#      against which each leave-one-predictor-out model is compared.
#
# Importantly, the full linear model itself does not depend on `xorder`.
# `xorder` affects only the order in which predictors are visited during the
# subsequent MFP backfitting cycles.


# -----------------------------------------------------------------------------
# order_variables() -----------------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the Full Linear Reference Model and Determine Predictor Visiting Order
#'
#' Fits the model containing every candidate predictor as an ordinary linear
#' term and determines the order in which predictors are visited by the MFP
#' backfitting algorithm.
#'
#' @details
#' The full linear reference model is fitted unconditionally because its null
#' and fitted-model deviances are retained by \code{fit_mfp()}.
#' This fit is invariant to \code{xorder}: changing the visiting order does not
#' change the predictors, likelihood, coefficients, or deviance of the full
#' linear model.
#'
#' Predictor ranking is a separate operation. It is performed only when all of
#' the following are true:
#' \itemize{
#'   \item more than one predictor is present; and
#'   \item \code{xorder} is \code{"ascending"} or \code{"descending"}.
#' }
#'
#' For significance-based ordering, one reduced model is fitted per predictor by
#' omitting that predictor from the full linear model. A likelihood-ratio test
#' compares the reduced model with the full model. Predictors are then ordered by
#' the resulting p-values:
#' \describe{
#'   \item{\code{"ascending"}}{Smallest p-value first; the predictor whose
#'     omission most strongly worsens model fit is visited first.}
#'   \item{\code{"descending"}}{Largest p-value first; the least significant
#'     predictor is visited first.}
#'   \item{\code{"original"}}{The original column order of \code{x} is
#'     retained and no reduced models are fitted.}
#' }
#'
#' With a single predictor, no ranking problem exists. The full linear reference
#' model is still fitted so that its deviance is available, and the sole
#' predictor is returned unchanged.
#'
#' @section Deviance convention:
#' Deviance is family-specific and is taken directly from the full reference
#' fit returned by \code{fit_model()}. For GLMs, \code{null_deviance} and
#' \code{linear_deviance} are the \code{null.deviance} and \code{deviance}
#' values computed by \code{stats::glm.fit()} or \code{stats::glm()}. For Cox
#' models, they are minus twice the null and fitted partial log-likelihoods.
#'
#' These reported deviances are separate from the log-likelihoods used for the
#' leave-one-predictor-out likelihood-ratio tests below.
#'
#' @param xorder Character scalar controlling the predictor visiting order.
#'   Supported values are \code{"ascending"}, \code{"descending"}, and
#'   \code{"original"}.
#' @param x Numeric design matrix with one column per candidate predictor and
#'   one row per observation. The matrix excludes the intercept for both GLM and
#'   Cox models. Column names identify the predictors.
#' @param y Response used to fit the models. For GLMs, this may be a numeric
#'   vector, a factor response accepted by \code{stats::glm()}, or a two-column
#'   matrix of grouped binomial counts. For Cox models, this must be a
#'   \code{survival::Surv()} object.
#' @param family GLM family object used by \code{fit_model()}, or character
#'   \code{"cox"} for a Cox proportional-hazards model.
#' @param family_string Normalized character family name, for example
#'   \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, or
#'   \code{"cox"}.
#' @param weights Optional observation weights passed to \code{fit_model()}.
#' @param offset Optional linear-predictor offset passed to \code{fit_model()}.
#' @param strata Optional Cox stratification object. Ignored for GLMs.
#' @param method Cox tie-handling method. Ignored for GLMs.
#' @param control Model-fitting control object passed to the GLM or Cox fitting
#'   path.
#' @param nocenter Cox centering-suppression argument. Ignored for GLMs.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{variables_ordered}}{Character vector containing every
#'     predictor name exactly once, in the requested visiting order.}
#'   \item{\code{null_deviance}}{Family-specific deviance of the null model
#'     associated with the full linear reference fit.}
#'   \item{\code{linear_deviance}}{Family-specific deviance of the full
#'     linear reference model.}
#' }
#'
#' @keywords internal
#' @noRd
order_variables <- function(xorder = "ascending",
                            x,
                            y,
                            family,
                            family_string,
                            weights = NULL,
                            offset = NULL,
                            strata = NULL,
                            method = NULL,
                            control = NULL,
                            nocenter = NULL) {
  predictor_names <- colnames(x)
  n_predictors <- ncol(x)
  
  # The full linear reference model is required independently of predictor
  # ordering. Fit it once and reuse its log-likelihood for all reduced-model
  # comparisons below.
  full_reference <- fit_full_linear_reference(
    x = x,
    y = y,
    family = family,
    family_string = family_string,
    weights = weights,
    offset = offset,
    strata = strata,
    method = method,
    control = control,
    nocenter = nocenter
  )
  
  # No reduced-model fits are needed when the user requests the original order
  # or when only one predictor is available. The full reference fit above is
  # nevertheless retained because its deviance is required for reporting.
  rank_predictors <- n_predictors > 1L && !identical(xorder, "original")
  
  variables_ordered <- if (rank_predictors) {
    order_variables_by_significance(
      xorder = xorder,
      x = x,
      y = y,
      family = family,
      family_string = family_string,
      weights = weights,
      offset = offset,
      strata = strata,
      method = method,
      control = control,
      nocenter = nocenter,
      full_reference = full_reference
    )
  } else {
    predictor_names
  }
  
  list(
    variables_ordered = variables_ordered,
    null_deviance = full_reference$null_deviance,
    linear_deviance = full_reference$model_deviance
  )
}


# -----------------------------------------------------------------------------
# fit_full_linear_reference() -------------------------------------------------
# -----------------------------------------------------------------------------

#' Fit the Full Linear Reference Model
#'
#' Fits a model containing all candidate predictors as ordinary linear terms.
#' The fitted model supplies both the null and full-linear-model deviances and,
#' when
#' significance-based ordering is requested, the full-model likelihood and
#' degrees of freedom used in leave-one-predictor-out likelihood-ratio tests.
#'
#' @inheritParams order_variables
#'
#' @return A model-fit wrapper returned by \code{fit_model()}, including
#'   \code{logl}, \code{df}, \code{null_deviance},
#'   \code{model_deviance}, coefficients, residual information, and the
#'   underlying fitted object.
#'
#' @keywords internal
#' @noRd
fit_full_linear_reference <- function(x,
                                      y,
                                      family,
                                      family_string,
                                      weights,
                                      offset,
                                      strata,
                                      method,
                                      control,
                                      nocenter) {
  fit_model(
    x = x,
    y = y,
    family = family,
    family_string = family_string,
    weights = weights,
    offset = offset,
    method = method,
    strata = strata,
    control = control,
    rownames = rownames(x),
    nocenter = nocenter,
    fast = TRUE
  )
}


# -----------------------------------------------------------------------------
# order_variables_by_significance() -------------------------------------------
# -----------------------------------------------------------------------------

#' Order Predictors by Leave-One-Out Significance
#'
#' Orders predictors using likelihood-ratio tests comparing the full linear
#' reference model with models obtained by removing one predictor at a time.
#'
#' @details
#' This helper does not fit the full linear model. The caller supplies
#' \code{full_reference}, ensuring that the invariant full-model fit is computed
#' only once and reused for every reduced-model comparison.
#'
#' A predictor may correspond to a single design-matrix column. The test degrees
#' of freedom are calculated as the difference between the parameter counts of
#' the full and reduced fits. If that difference is not positive, or if a valid
#' likelihood-ratio statistic cannot be formed, the predictor receives
#' \code{NA} as its ordering p-value and is placed after predictors with valid
#' p-values. No predictor is dropped from the returned order.
#'
#' Ties retain the original column order of \code{x}.
#'
#' @inheritParams order_variables
#' @param full_reference Model-fit wrapper returned by
#'   \code{fit_full_linear_reference()} for the model containing all predictors
#'   linearly.
#'
#' @return Character vector containing all predictor names in significance-based
#'   visiting order.
#'
#' @keywords internal
#' @noRd
order_variables_by_significance <- function(xorder,
                                            x,
                                            y,
                                            family,
                                            family_string,
                                            weights,
                                            offset,
                                            strata,
                                            method,
                                            control,
                                            nocenter,
                                            full_reference) {
  predictor_names <- colnames(x)
  n_predictors <- ncol(x)
  
  # Initialize with NA rather than zero. A failed or non-identifiable comparison
  # must not be interpreted as overwhelming evidence against the predictor.
  p_values <- stats::setNames(
    rep(NA_real_, n_predictors),
    predictor_names
  )
  
  for (predictor_index in seq_len(n_predictors)) {
    # Remove exactly one predictor while preserving matrix structure. This loop
    # is entered only when n_predictors > 1, so reduced_x remains non-empty.
    reduced_x <- x[, -predictor_index, drop = FALSE]
    
    reduced_fit <- fit_model(
      x = reduced_x,
      y = y,
      family = family,
      family_string = family_string,
      weights = weights,
      offset = offset,
      method = method,
      strata = strata,
      control = control,
      rownames = rownames(x),
      nocenter = nocenter,
      fast = TRUE
    )
    
    lrt_df <- full_reference$df - reduced_fit$df
    lrt_statistic <- 2 * (full_reference$logl - reduced_fit$logl)
    
    # A valid nested-model likelihood-ratio test requires a positive df
    # difference and finite likelihoods. Small negative statistics can occur
    # from numerical rounding, so truncate such values to zero.
    if (is.finite(lrt_df) && lrt_df > 0L && is.finite(lrt_statistic)) {
      p_values[predictor_index] <- stats::pchisq(
        q = max(0, lrt_statistic),
        df = lrt_df,
        lower.tail = FALSE
      )
    }
  }
  
  # Convert the requested p-value direction into a single ascending score.
  # A secondary index keeps the original column order when p-values are tied.
  ordering_score <- if (identical(xorder, "descending")) {
    -p_values
  } else {
    p_values
  }
  
  ordering_index <- order(
    ordering_score,
    seq_along(ordering_score),
    na.last = TRUE
  )
  
  predictor_names[ordering_index]
}
