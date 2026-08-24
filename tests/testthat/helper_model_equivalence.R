# Shared GLM and Cox equivalence oracles.

# 8.1 Prediction equivalence against stats::glm()
# -----------------------------------------------------------------------------
# These tests deliberately disable FP selection/transformation complexity
# so that mfp2() should reduce to the corresponding base glm() fit:
#   - df = 1: linear terms only
#   - select = 1 and alpha = 1: keep all ordinary predictors
#   - shift = 0, scale = 1, center = FALSE: no preprocessing-induced change
#   - no spike-at-zero, zero/catzero handling, or ACD transformation
#
# mfp2 stores ordinary transformed columns as x.1, x1.1, and so on. These
# tests compare numerical values and set xorder = "original", so coefficient
# and covariance order must agree with glm() without name canonicalization.


# Shared assertion helper for ordinary GLM equivalence tests.
#
# The helper verifies three increasingly independent layers:
#   1. mfp2() and stats::glm() fit the same statistical model;
#   2. both prediction methods return the same link values, responses, and
#      link-scale standard errors; and
#   3. those values agree with direct matrix algebra using X %*% beta,
#      the formula offset, the inverse-link function, and X V X'.
#
# Using a manual oracle is important because two prediction methods can agree
# while sharing the same reconstruction error. The X beta calculation checks
# the fitted coefficient order, factor expansion, offset handling, link
# inversion, and covariance propagation independently.
expect_mfp2_glm_predictions_equal <- function(dat,
                                              formula,
                                              family_name,
                                              newdata_cols,
                                              tolerance = 1e-8) {
  family_fun <- switch(
    family_name,
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    stop("Unsupported test family: ", family_name, call. = FALSE)
  )

  fit_mfp2 <- mfp2(
    formula,
    data = dat,
    family = family_name,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  fit_glm <- stats::glm(
    formula,
    data = dat,
    family = family_fun
  )

  nd <- dat[1:25, newdata_cols, drop = FALSE]

  # Request link-scale standard errors from both methods. Standard errors are
  # naturally calculated on the linear-predictor scale by predict.glm().
  pred_mfp2_link <- predict(
    fit_mfp2,
    newdata = nd,
    type = "link",
    se.fit = TRUE
  )
  pred_glm_link <- predict(
    fit_glm,
    newdata = nd,
    type = "link",
    se.fit = TRUE
  )

  pred_mfp2_response <- predict(
    fit_mfp2,
    newdata = nd,
    type = "response"
  )
  pred_glm_response <- predict(
    fit_glm,
    newdata = nd,
    type = "response"
  )

  # Build the reference model frame from raw newdata. This evaluates factor
  # contrasts and formula offsets with the same terms object used by glm().
  reference_terms <- stats::delete.response(stats::terms(fit_glm))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = nd,
    xlev = fit_glm$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = fit_glm$contrasts
  )

  beta <- stats::coef(fit_glm)
  beta_vcov <- stats::vcov(fit_glm)

  # xorder = "original" makes the model-matrix and coefficient order
  # positional, so no name-based reordering is required.
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  # model.offset() returns NULL when the formula has no offset. In that case the
  # additive offset contribution is exactly zero for every prediction row.
  manual_offset <- stats::model.offset(reference_frame)
  if (is.null(manual_offset)) {
    manual_offset <- rep(0, nrow(manual_x))
  }

  # Independent prediction calculations.
  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_response <- as.numeric(fit_glm$family$linkinv(manual_link))
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  # Fitted-model equivalence checks. These detect differences that may be hidden
  # when predictions happen to be evaluated at only a small set of rows.
  expect_mfp2_glm_parameters_equal(
    fit_mfp2,
    fit_glm,
    tolerance = tolerance
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_glm)),
    tolerance = tolerance
  )
  expect_equal(
    unname(stats::fitted(fit_mfp2)),
    unname(stats::fitted(fit_glm)),
    tolerance = tolerance
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_glm$linear.predictors),
    tolerance = tolerance
  )

  # mfp2() must agree with glm() for both prediction scales and link-scale SEs.
  expect_equal(
    unname(pred_mfp2_link$fit),
    unname(pred_glm_link$fit),
    tolerance = tolerance
  )
  expect_equal(
    unname(pred_mfp2_link$se.fit),
    unname(pred_glm_link$se.fit),
    tolerance = tolerance
  )
  expect_equal(
    unname(pred_mfp2_response),
    unname(pred_glm_response),
    tolerance = tolerance
  )

  # Both methods must also agree with the independently constructed oracle.
  expect_equal(unname(pred_mfp2_link$fit), manual_link, tolerance = tolerance)
  expect_equal(unname(pred_glm_link$fit), manual_link, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_link$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_glm_link$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_response), manual_response, tolerance = tolerance)
  expect_equal(unname(pred_glm_response), manual_response, tolerance = tolerance)
}

# Shared helper for tests that use non-default family/link objects or special
# formula constructions. It performs the same manual X beta checks on two
# already fitted GLM objects.
expect_glm_objects_and_manual_prediction_equal <- function(fit_mfp2,
                                                           fit_glm,
                                                           newdata,
                                                           tolerance = 1e-8) {
  pred_mfp2_link <- predict(
    fit_mfp2,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  pred_glm_link <- predict(
    fit_glm,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  pred_mfp2_response <- predict(fit_mfp2, newdata = newdata, type = "response")
  pred_glm_response <- predict(fit_glm, newdata = newdata, type = "response")

  reference_terms <- stats::delete.response(stats::terms(fit_glm))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = newdata,
    xlev = fit_glm$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = fit_glm$contrasts
  )

  beta <- stats::coef(fit_glm)
  beta_vcov <- stats::vcov(fit_glm)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_offset <- stats::model.offset(reference_frame)
  if (is.null(manual_offset)) {
    manual_offset <- rep(0, nrow(manual_x))
  }

  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_response <- as.numeric(fit_glm$family$linkinv(manual_link))
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_mfp2_glm_parameters_equal(
    fit_mfp2,
    fit_glm,
    tolerance = tolerance
  )
  expect_equal(as.numeric(logLik(fit_mfp2)), as.numeric(logLik(fit_glm)), tolerance = tolerance)
  expect_equal(unname(fitted(fit_mfp2)), unname(fitted(fit_glm)), tolerance = tolerance)
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_glm$linear.predictors),
    tolerance = tolerance
  )

  expect_equal(unname(pred_mfp2_link$fit), unname(pred_glm_link$fit), tolerance = tolerance)
  expect_equal(unname(pred_mfp2_link$se.fit), unname(pred_glm_link$se.fit), tolerance = tolerance)
  expect_equal(unname(pred_mfp2_response), unname(pred_glm_response), tolerance = tolerance)

  expect_equal(unname(pred_mfp2_link$fit), manual_link, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_link$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_response), manual_response, tolerance = tolerance)
}

# Test purpose: 8.1.1 Gaussian GLM equivalence without offset.

# 8.2 Prediction equivalence against survival::coxph()
# -----------------------------------------------------------------------------
# These tests use the simplest Cox configuration where mfp2() should reduce to
# the corresponding survival::coxph() fit:
#   - df = 1: linear terms only
#   - select = 1 and alpha = 1: keep all ordinary predictors
#   - shift = 0, scale = 1, center = FALSE: no preprocessing-induced change
#   - matched tie handling and reference scale
#   - formula-level strata reconstructed during prediction

# Compare a forced-linear mfp2 Cox fit with survival::coxph() by coefficient
# position. Every caller uses xorder = "original", so no name repair or
# coefficient reordering is required. The helper checks coefficients,
# coefficient standard errors, training linear predictors, newdata LP and risk
# predictions, newdata prediction standard errors, partial log-likelihood, and
# an independent X beta / X V X' calculation from the coxph model matrix.
expect_mfp2_cox_predictions_equal <- function(fit_mfp2,
                                              fit_coxph,
                                              newdata,
                                              mfp2_newdata = newdata,
                                              mfp2_predict_args = list(),
                                              tolerance = 1e-8) {
  coef_mfp2 <- stats::coef(fit_mfp2)
  coef_coxph <- stats::coef(fit_coxph)
  vcov_mfp2 <- stats::vcov(fit_mfp2)
  vcov_coxph <- stats::vcov(fit_coxph)

  expect_length(coef_mfp2, length(coef_coxph))
  expect_equal(dim(vcov_mfp2), dim(vcov_coxph))
  expect_equal(unname(coef_mfp2), unname(coef_coxph), tolerance = tolerance)
  expect_equal(unname(vcov_mfp2), unname(vcov_coxph), tolerance = tolerance)
  expect_equal(
    unname(sqrt(diag(vcov_mfp2))),
    unname(sqrt(diag(vcov_coxph))),
    tolerance = tolerance
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_coxph)),
    tolerance = tolerance
  )

  # Training predictions verify that fitted offsets and strata are stored and
  # applied consistently without requiring any newdata reconstruction.
  training_lp_mfp2 <- predict(fit_mfp2, type = "lp")
  training_lp_coxph <- predict(
    fit_coxph,
    type = "lp",
    reference = "zero"
  )
  expect_equal(
    unname(training_lp_mfp2),
    unname(training_lp_coxph),
    tolerance = tolerance
  )

  # Formula fits reconstruct strata and offsets from raw newdata. Matrix fits
  # supply them explicitly through mfp2_predict_args.
  pred_mfp2 <- do.call(
    stats::predict,
    c(
      list(
        object = fit_mfp2,
        newdata = mfp2_newdata,
        type = "lp",
        se.fit = TRUE
      ),
      mfp2_predict_args
    )
  )
  pred_coxph <- predict(
    fit_coxph,
    newdata = newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  risk_mfp2 <- do.call(
    stats::predict,
    c(
      list(
        object = fit_mfp2,
        newdata = mfp2_newdata,
        type = "risk"
      ),
      mfp2_predict_args
    )
  )
  risk_coxph <- predict(
    fit_coxph,
    newdata = newdata,
    type = "risk",
    reference = "zero"
  )

  # model.matrix.coxph() removes the non-estimated intercept and formula
  # specials such as strata() and offset(), leaving the coefficient design in
  # its fitted order. This is the independent X beta oracle for newdata.
  manual_x <- stats::model.matrix(fit_coxph, data = newdata)
  if ("(Intercept)" %in% colnames(manual_x)) {
    manual_x <- manual_x[, colnames(manual_x) != "(Intercept)", drop = FALSE]
  }
  expect_equal(ncol(manual_x), length(coef_coxph))
  expect_equal(dim(vcov_coxph), c(length(coef_coxph), length(coef_coxph)))

  reference_terms <- stats::delete.response(stats::terms(fit_coxph))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = newdata,
    xlev = fit_coxph$xlevels,
    na.action = stats::na.pass
  )
  new_offset <- stats::model.offset(reference_frame)
  training_frame <- stats::model.frame(fit_coxph)
  training_offset <- stats::model.offset(training_frame)
  expected_offset_reference <- if (is.null(training_offset)) {
    0
  } else {
    mean(training_offset)
  }

  # fit_mfp() stores the same offset origin used internally by
  # predict.coxph() once on the final mfp2 object. Package-owned manual Cox
  # prediction paths reuse this scalar rather than recomputing it in fit_model().
  expect_equal(
    fit_mfp2$cox_offset_reference,
    unname(expected_offset_reference),
    tolerance = tolerance
  )

  # predict.coxph() always recentres an offset by subtracting its mean in the
  # training model frame. This applies even when cox_reference = "zero"; that
  # option controls covariate centering, not offset centering. Reproduce that
  # convention explicitly in the independent oracle.
  if (is.null(new_offset)) {
    manual_offset <- rep(0, nrow(manual_x))
  } else {
    expect_false(is.null(training_offset))
    manual_offset <- as.numeric(new_offset - expected_offset_reference)
  }

  manual_lp <- as.numeric(manual_x %*% coef_coxph + manual_offset)
  manual_risk <- exp(manual_lp)
  manual_variance <- rowSums((manual_x %*% vcov_coxph) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_equal(
    unname(pred_mfp2$fit),
    unname(pred_coxph$fit),
    tolerance = tolerance
  )
  expect_equal(
    unname(pred_mfp2$se.fit),
    unname(pred_coxph$se.fit),
    tolerance = tolerance
  )
  expect_equal(unname(pred_mfp2$fit), manual_lp, tolerance = tolerance)
  expect_equal(unname(pred_coxph$fit), manual_lp, tolerance = tolerance)
  expect_equal(unname(pred_mfp2$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_coxph$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(risk_mfp2), unname(risk_coxph), tolerance = tolerance)
  expect_equal(unname(risk_mfp2), manual_risk, tolerance = tolerance)
  expect_equal(unname(risk_coxph), manual_risk, tolerance = tolerance)
}

# Exact partial likelihood is deliberately rejected because the low-level Cox
# candidate fitter used during MFP/MFPI selection does not implement it. These
# tests cover both matrix and formula entry points so the unsupported method
# cannot bypass public argument validation.

# Generate one stable Cox data set containing continuous and categorical
# predictors, two potential stratification factors, and an offset. Individual
# equivalence tests use different subsets of these model components.
make_cox_equivalence_data <- function(n = 480L, seed = 8200L) {
  set.seed(seed)

  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    stratum1 = factor(sample(c("S1", "S2", "S3"), n, replace = TRUE)),
    stratum2 = factor(sample(c("T1", "T2"), n, replace = TRUE)),
    off = stats::rnorm(n, mean = 0, sd = 0.25)
  )

  group_effect <- c(A = 0, B = 0.45, C = -0.35)[as.character(dat$group)]
  baseline_multiplier <-
    c(S1 = 0.75, S2 = 1.00, S3 = 1.35)[as.character(dat$stratum1)] *
    c(T1 = 0.85, T2 = 1.20)[as.character(dat$stratum2)]
  eta <- 0.38 * dat$x1 - 0.25 * dat$x2 + group_effect + dat$off
  event_time <- stats::rexp(
    n,
    rate = 0.020 * baseline_multiplier * exp(eta)
  )
  censor_time <- stats::rexp(n, rate = 0.010)
  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  dat
}
