# Complete forced-linear Cox equivalence --------------------------------------
#
# These tests make test_linear_equivalence_cox.R the authoritative ordinary
# Cox regression contract. They compare the final fit, likelihood metadata,
# residuals, relative and absolute predictions, prediction standard errors,
# reference scales, an independent matrix oracle, and mfp2's term/contrast
# paths. Formula, stratified, matrix, Breslow, and Efron routes are separate so
# a failure identifies the affected interface immediately.


expect_complete_cox_linear_equivalence <- function(
    fit_mfp2,
    fit_coxph,
    relative_newdata,
    absolute_newdata,
    mfp2_relative_newdata = relative_newdata,
    mfp2_absolute_newdata = absolute_newdata,
    mfp2_predict_args = list(),
    references = c("zero", "sample"),
    term_newdata = mfp2_relative_newdata,
    info,
    tolerance = 1e-8) {
  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "coxph")
  expect_identical(fit_mfp2$family_string, "cox", info = info)
  expect_identical(fit_mfp2$fitter, "base", info = info)
  expect_identical(fit_mfp2$method, fit_coxph$method, info = info)
  expect_true(isTRUE(fit_mfp2$convergence_mfp), info = info)

  beta_mfp2 <- stats::coef(fit_mfp2)
  beta_coxph <- stats::coef(fit_coxph)
  covariance_mfp2 <- stats::vcov(fit_mfp2)
  covariance_coxph <- stats::vcov(fit_coxph)
  expect_equal(length(beta_mfp2), length(beta_coxph), info = info)
  expect_equal(
    unname(beta_mfp2),
    unname(beta_coxph),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(covariance_mfp2),
    unname(covariance_coxph),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(sqrt(diag(covariance_mfp2))),
    unname(sqrt(diag(covariance_coxph))),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$loglik),
    unname(fit_coxph$loglik),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_coxph$linear.predictors),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$residuals),
    unname(fit_coxph$residuals),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$means),
    unname(fit_coxph$means),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$n, fit_coxph$n, info = info)
  expect_equal(fit_mfp2$nevent, fit_coxph$nevent, info = info)
  expect_equal(fit_mfp2$nevents, fit_coxph$nevent, info = info)
  expect_equal(fit_mfp2$iter, fit_coxph$iter, info = info)
  expect_equal(
    unname(fit_mfp2$weights),
    unname(fit_coxph$weights),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$concordance),
    unname(fit_coxph$concordance),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$score),
    unname(fit_coxph$score),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$wald.test),
    unname(fit_coxph$wald.test),
    tolerance = tolerance,
    info = info
  )

  for (residual_type in c(
    "martingale", "deviance", "score", "schoenfeld",
    "dfbeta", "dfbetas", "scaledsch"
  )) {
    expect_equal(
      unname(stats::residuals(fit_mfp2, type = residual_type)),
      unname(stats::residuals(fit_coxph, type = residual_type)),
      tolerance = tolerance,
      info = paste(info, residual_type, "residuals")
    )
  }

  loglik_mfp2 <- stats::logLik(fit_mfp2)
  loglik_coxph <- stats::logLik(fit_coxph)
  expect_equal(
    as.numeric(loglik_mfp2),
    as.numeric(loglik_coxph),
    tolerance = tolerance,
    info = info
  )
  expect_equal(attr(loglik_mfp2, "df"), attr(loglik_coxph, "df"), info = info)
  expect_equal(
    attr(loglik_mfp2, "nobs"),
    attr(loglik_coxph, "nobs"),
    info = info
  )
  expect_equal(stats::nobs(fit_mfp2), stats::nobs(fit_coxph), info = info)
  expect_equal(
    stats::AIC(fit_mfp2),
    stats::AIC(fit_coxph),
    tolerance = tolerance,
    info = info
  )

  expect_equal(
    fit_mfp2$null_logl,
    unname(fit_coxph$loglik[1L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_logl,
    unname(fit_coxph$loglik[2L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$mfp_df, attr(loglik_coxph, "df"), info = info)
  expect_equal(
    fit_mfp2$null_deviance,
    -2 * unname(fit_coxph$loglik[1L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_deviance,
    -2 * unname(fit_coxph$loglik[2L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$linear_logl,
    unname(fit_coxph$loglik[2L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$linear_df, attr(loglik_coxph, "df"), info = info)
  expect_equal(
    fit_mfp2$linear_deviance,
    -2 * unname(fit_coxph$loglik[2L]),
    tolerance = tolerance,
    info = info
  )

  mfp2_predict <- function(newdata = NULL, type, se.fit, reference = NULL) {
    arguments <- list(object = fit_mfp2, type = type, se.fit = se.fit)
    if (!is.null(newdata)) arguments$newdata <- newdata
    if (!is.null(reference)) arguments$cox_reference <- reference
    if (!is.null(newdata)) arguments <- c(arguments, mfp2_predict_args)
    do.call(stats::predict, arguments)
  }

  # Relative predictions are checked for every meaningful reference on both the
  # training rows and raw newdata. Risk SEs matter because predict.coxph() has a
  # distinct risk-scale SE convention; checking only exp(lp) would miss it.
  newdata_predictions <- list()
  for (reference in references) {
    for (type in c("lp", "risk")) {
      training_mfp2 <- mfp2_predict(
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      training_coxph <- stats::predict(
        fit_coxph,
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      expect_native_prediction_equal(
        training_mfp2,
        training_coxph,
        tolerance = tolerance,
        info = paste(info, "training", type, reference)
      )

      new_mfp2 <- mfp2_predict(
        newdata = mfp2_relative_newdata,
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      new_coxph <- stats::predict(
        fit_coxph,
        newdata = relative_newdata,
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      expect_native_prediction_equal(
        new_mfp2,
        new_coxph,
        tolerance = tolerance,
        info = paste(info, "newdata", type, reference)
      )
      newdata_predictions[[paste(type, reference, sep = ":")]] <- new_mfp2
    }
  }

  # Independent newdata oracle for the two global reference scales. For a
  # stratum-specific reference, predict.coxph() remains the oracle because each
  # stratum has a different training mean.
  manual_x <- stats::model.matrix(fit_coxph, data = relative_newdata)
  if ("(Intercept)" %in% colnames(manual_x)) {
    manual_x <- manual_x[, colnames(manual_x) != "(Intercept)", drop = FALSE]
  }
  manual_x <- manual_x[, names(beta_coxph), drop = FALSE]
  reference_terms <- stats::delete.response(stats::terms(fit_coxph))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = relative_newdata,
    xlev = fit_coxph$xlevels,
    na.action = stats::na.pass
  )
  new_offset <- stats::model.offset(reference_frame)
  training_offset <- stats::model.offset(stats::model.frame(fit_coxph))
  offset_origin <- if (is.null(training_offset)) 0 else mean(training_offset)
  manual_offset <- if (is.null(new_offset)) {
    rep(0, nrow(manual_x))
  } else {
    as.numeric(new_offset - offset_origin)
  }
  expect_equal(
    fit_mfp2$cox_offset_reference,
    unname(offset_origin),
    tolerance = tolerance,
    info = info
  )

  for (reference in intersect(references, c("zero", "sample"))) {
    center <- if (identical(reference, "zero")) {
      rep(0, ncol(manual_x))
    } else {
      unname(fit_coxph$means[names(beta_coxph)])
    }
    centered_x <- sweep(manual_x, 2L, center, FUN = "-")
    manual_lp <- as.numeric(centered_x %*% beta_coxph + manual_offset)
    manual_variance <- unname(rowSums(
      (centered_x %*% covariance_coxph) * centered_x
    ))
    manual_lp_se <- sqrt(pmax(manual_variance, 0))
    manual_risk <- exp(manual_lp)
    manual_risk_se <- sqrt(manual_risk) * manual_lp_se

    got_lp <- newdata_predictions[[paste("lp", reference, sep = ":")]]
    got_risk <- newdata_predictions[[paste("risk", reference, sep = ":")]]
    expect_equal(
      unname(got_lp$fit),
      manual_lp,
      tolerance = tolerance,
      info = paste(info, "manual lp", reference)
    )
    expect_equal(
      unname(got_lp$se.fit),
      manual_lp_se,
      tolerance = tolerance,
      info = paste(info, "manual lp SE", reference)
    )
    expect_equal(
      unname(got_risk$fit),
      manual_risk,
      tolerance = tolerance,
      info = paste(info, "manual risk", reference)
    )
    expect_equal(
      unname(got_risk$se.fit),
      manual_risk_se,
      tolerance = tolerance,
      info = paste(info, "manual risk SE", reference)
    )
  }

  zero_lp <- newdata_predictions[["lp:zero"]]$fit
  expect_equal(
    unname(mfp2_predict(
      newdata = mfp2_relative_newdata,
      type = "lp",
      se.fit = FALSE
    )),
    unname(zero_lp),
    tolerance = tolerance,
    info = paste(info, "default lp")
  )
  expect_equal(
    unname(mfp2_predict(
      newdata = mfp2_relative_newdata,
      type = "link",
      se.fit = FALSE
    )),
    unname(zero_lp),
    tolerance = tolerance,
    info = paste(info, "link alias")
  )

  # Absolute Cox predictions verify the baseline-hazard route independently of
  # relative reference centering. Both point predictions and native SEs are
  # compared, then survival = exp(-expected) and its delta-method SE are checked.
  for (prediction_data in c("training", "newdata")) {
    is_training <- identical(prediction_data, "training")
    mfp2_data <- if (is_training) NULL else mfp2_absolute_newdata
    coxph_data <- if (is_training) NULL else absolute_newdata

    expected_mfp2 <- mfp2_predict(
      newdata = mfp2_data,
      type = "expected",
      se.fit = TRUE
    )
    survival_mfp2 <- mfp2_predict(
      newdata = mfp2_data,
      type = "survival",
      se.fit = TRUE
    )
    expected_coxph <- if (is_training) {
      stats::predict(fit_coxph, type = "expected", se.fit = TRUE)
    } else {
      stats::predict(
        fit_coxph,
        newdata = coxph_data,
        type = "expected",
        se.fit = TRUE
      )
    }
    survival_coxph <- if (is_training) {
      stats::predict(fit_coxph, type = "survival", se.fit = TRUE)
    } else {
      stats::predict(
        fit_coxph,
        newdata = coxph_data,
        type = "survival",
        se.fit = TRUE
      )
    }

    expect_native_prediction_equal(
      expected_mfp2,
      expected_coxph,
      tolerance = tolerance,
      info = paste(info, prediction_data, "expected")
    )
    expect_native_prediction_equal(
      survival_mfp2,
      survival_coxph,
      tolerance = tolerance,
      info = paste(info, prediction_data, "survival")
    )
    expect_equal(
      unname(survival_mfp2$fit),
      exp(-unname(expected_mfp2$fit)),
      tolerance = tolerance,
      info = paste(info, prediction_data, "survival identity")
    )
    expect_equal(
      unname(survival_mfp2$se.fit),
      unname(survival_mfp2$fit * expected_mfp2$se.fit),
      tolerance = tolerance,
      info = paste(info, prediction_data, "survival SE identity")
    )
  }

  expect_linear_term_and_contrast_equal(
    object = fit_mfp2,
    newdata = term_newdata,
    term = "x1",
    ref = 0.25,
    has_intercept = FALSE,
    tolerance = tolerance
  )
}


test_that("df = 1 Cox formula fits match coxph completely for both tie methods", {
  dat <- make_cox_equivalence_data(n = 480L, seed = 4201L)
  dat$time <- pmax(round(dat$time, digits = 1L), 0.1)
  expect_true(anyDuplicated(dat$time[dat$status == 1L]) > 0L)

  relative_newdata <- dat[
    seq_len(36L),
    c("x1", "x2", "group", "off"),
    drop = FALSE
  ]
  absolute_newdata <- dat[
    seq_len(36L),
    c("time", "status", "x1", "x2", "group", "off"),
    drop = FALSE
  ]

  for (ties in c("breslow", "efron")) {
    fit_mfp2 <- mfp2(
      survival::Surv(time, status) ~ x1 + x2 + group + offset(off),
      data = dat,
      family = "cox",
      weights = case_weight,
      ties = ties,
      nocenter = NULL,
      cycles = 1,
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      xorder = "original",
      verbose = FALSE
    )
    fit_coxph <- survival::coxph(
      survival::Surv(time, status) ~ x1 + x2 + group + offset(off),
      data = dat,
      weights = case_weight,
      ties = ties,
      nocenter = NULL,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )

    expect_equal(
      sort(get_selected_variable_names(fit_mfp2)),
      sort(c("x1", "x2", "group")),
      info = ties
    )
    expect_complete_cox_linear_equivalence(
      fit_mfp2,
      fit_coxph,
      relative_newdata = relative_newdata,
      absolute_newdata = absolute_newdata,
      info = paste("formula", ties),
      tolerance = 1e-8
    )
  }
})


test_that("df = 1 Cox formula fit with formula strata and data-mask offset matches coxph completely", {
  dat <- make_cox_equivalence_data(n = 480L, seed = 4202L)
  rows <- c(2L, 11L, 25L, 48L, 79L, 116L, 175L, 244L, 331L, 408L)
  relative_newdata <- dat[
    rows,
    c("x1", "x2", "group", "stratum1", "off"),
    drop = FALSE
  ]
  absolute_newdata <- dat[
    rows,
    c("time", "status", "x1", "x2", "group", "stratum1", "off"),
    drop = FALSE
  ]

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 + group + strata(stratum1),
    data = dat,
    family = "cox",
    weights = case_weight,
    offset = off,
    ties = "breslow",
    nocenter = NULL,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~
      x1 + x2 + group + strata(stratum1) + offset(off),
    data = dat,
    weights = case_weight,
    ties = "breslow",
    nocenter = NULL,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  expect_setequal(
    get_selected_variable_names(fit_mfp2),
    c("x1", "x2", "group")
  )
  expect_complete_cox_linear_equivalence(
    fit_mfp2,
    fit_coxph,
    relative_newdata = relative_newdata,
    absolute_newdata = absolute_newdata,
    mfp2_predict_args = list(
      newoffset = relative_newdata$off
    ),
    references = c("zero", "sample", "strata"),
    info = "formula strata and data-mask offset",
    tolerance = 1e-8
  )
})


test_that("df = 1 Cox matrix fit with strata and offset matches coxph completely", {
  dat <- make_cox_equivalence_data(n = 480L, seed = 4203L)
  x <- as.matrix(dat[, c("x1", "x2"), drop = FALSE])
  y <- survival::Surv(dat$time, dat$status)
  rows <- c(3L, 17L, 36L, 71L, 122L, 189L, 278L, 367L, 451L)
  newx <- x[rows, , drop = FALSE]
  absolute_mfp2 <- as.data.frame(newx, check.names = FALSE)
  absolute_mfp2$prediction_response <- I(y[rows])
  reference_newdata <- dat[
    rows,
    c("time", "status", "x1", "x2", "stratum1", "off"),
    drop = FALSE
  ]
  prediction_args <- list(
    strata = dat$stratum1[rows],
    newoffset = dat$off[rows]
  )

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "cox",
    weights = dat$case_weight,
    strata = dat$stratum1,
    offset = dat$off,
    ties = "breslow",
    nocenter = NULL,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~
      x1 + x2 + strata(stratum1) + offset(off),
    data = dat,
    weights = case_weight,
    ties = "breslow",
    nocenter = NULL,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  expect_setequal(
    get_selected_variable_names(fit_mfp2),
    c("x1", "x2")
  )
  expect_complete_cox_linear_equivalence(
    fit_mfp2,
    fit_coxph,
    relative_newdata = reference_newdata,
    absolute_newdata = reference_newdata,
    mfp2_relative_newdata = newx,
    mfp2_absolute_newdata = absolute_mfp2,
    mfp2_predict_args = prediction_args,
    references = c("zero", "sample", "strata"),
    term_newdata = as.data.frame(newx, check.names = FALSE),
    info = "matrix strata offset",
    tolerance = 1e-8
  )
})
