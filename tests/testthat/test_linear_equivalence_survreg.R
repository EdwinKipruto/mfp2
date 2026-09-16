# Complete forced-linear parametric-survival equivalence ----------------------
#
# The first test loops over every distribution exported by survival. Additional
# tests isolate separate scale strata and the matrix interface with a fixed
# scale. Together they compare fitted quantities, likelihood metadata, all
# native complete-model survreg prediction types and aliases with SEs, a manual
# location-scale predictor, and mfp2's own term/contrast implementation.


expect_complete_survreg_linear_equivalence <- function(
    fit_mfp2,
    fit_survreg,
    reference_newdata,
    mfp2_newdata = reference_newdata,
    mfp2_predict_args = list(),
    term_newdata = mfp2_newdata,
    info,
    tolerance = 1e-7) {
  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "survreg")
  expect_identical(fit_mfp2$family_string, "survreg", info = info)
  expect_identical(fit_mfp2$fitter, "base", info = info)
  expect_null(fit_mfp2$family$prepared, info = info)
  expect_identical(fit_mfp2$dist, fit_survreg$dist, info = info)
  expect_true(isTRUE(fit_mfp2$convergence_mfp), info = info)

  beta_mfp2 <- stats::coef(fit_mfp2)
  beta_survreg <- stats::coef(fit_survreg)
  covariance_mfp2 <- stats::vcov(fit_mfp2)
  covariance_survreg <- stats::vcov(fit_survreg)
  expect_equal(length(beta_mfp2), length(beta_survreg), info = info)
  expect_equal(
    unname(beta_mfp2),
    unname(beta_survreg),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(covariance_mfp2),
    unname(covariance_survreg),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(sqrt(diag(covariance_mfp2))),
    unname(sqrt(diag(covariance_survreg))),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$scale),
    unname(fit_survreg$scale),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$loglik),
    unname(fit_survreg$loglik),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_survreg$linear.predictors),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$residuals),
    unname(fit_survreg$residuals),
    tolerance = tolerance,
    info = info
  )
  expect_equal(unname(fit_mfp2$df), unname(fit_survreg$df), info = info)
  expect_equal(unname(fit_mfp2$idf), unname(fit_survreg$idf), info = info)
  expect_equal(fit_mfp2$iter, fit_survreg$iter, info = info)
  expect_equal(
    unname(fit_mfp2$score),
    unname(fit_survreg$score),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$weights),
    unname(fit_survreg$weights),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$parms),
    unname(fit_survreg$parms),
    tolerance = tolerance,
    info = info
  )

  for (residual_type in c(
    "response", "deviance", "dfbeta", "dfbetas", "working",
    "ldcase", "ldresp", "ldshape", "matrix"
  )) {
    expect_equal(
      unname(stats::residuals(fit_mfp2, type = residual_type)),
      unname(stats::residuals(fit_survreg, type = residual_type)),
      tolerance = tolerance,
      info = paste(info, residual_type, "residuals")
    )
  }

  loglik_mfp2 <- stats::logLik(fit_mfp2)
  loglik_survreg <- stats::logLik(fit_survreg)
  expect_equal(
    as.numeric(loglik_mfp2),
    as.numeric(loglik_survreg),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    attr(loglik_mfp2, "df"),
    attr(loglik_survreg, "df"),
    info = info
  )
  expect_equal(
    attr(loglik_mfp2, "nobs"),
    attr(loglik_survreg, "nobs"),
    info = info
  )
  expect_equal(stats::nobs(fit_mfp2), stats::nobs(fit_survreg), info = info)
  expect_equal(
    stats::AIC(fit_mfp2),
    stats::AIC(fit_survreg),
    tolerance = tolerance,
    info = info
  )

  expect_equal(
    fit_mfp2$null_logl,
    unname(fit_survreg$loglik[1L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_logl,
    as.numeric(loglik_survreg),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$mfp_df, attr(loglik_survreg, "df"), info = info)
  expect_equal(
    fit_mfp2$null_deviance,
    -2 * unname(fit_survreg$loglik[1L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_deviance,
    -2 * as.numeric(loglik_survreg),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$linear_logl,
    as.numeric(loglik_survreg),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$linear_df, attr(loglik_survreg, "df"), info = info)
  expect_equal(
    fit_mfp2$linear_deviance,
    -2 * as.numeric(loglik_survreg),
    tolerance = tolerance,
    info = info
  )

  mfp2_predict <- function(newdata = NULL, type, se.fit, dots = list()) {
    arguments <- list(object = fit_mfp2, type = type, se.fit = se.fit)
    if (!is.null(newdata)) {
      arguments$newdata <- newdata
      arguments <- c(arguments, mfp2_predict_args)
    }
    do.call(stats::predict, c(arguments, dots))
  }

  survreg_predict <- function(newdata = NULL, type, se.fit, dots = list()) {
    arguments <- list(object = fit_survreg, type = type, se.fit = se.fit)
    if (!is.null(newdata)) arguments$newdata <- newdata
    do.call(stats::predict, c(arguments, dots))
  }

  prediction_specs <- list(
    link = list(),
    lp = list(),
    linear = list(),
    response = list(),
    quantile = list(p = c(0.25, 0.5, 0.75)),
    uquantile = list(p = c(0.25, 0.5, 0.75))
  )
  newdata_link <- NULL

  for (context in c("training", "newdata")) {
    is_training <- identical(context, "training")
    mfp2_data <- if (is_training) NULL else mfp2_newdata
    reference_data <- if (is_training) NULL else reference_newdata

    for (prediction_type in names(prediction_specs)) {
      dots <- prediction_specs[[prediction_type]]
      got <- mfp2_predict(
        newdata = mfp2_data,
        type = prediction_type,
        se.fit = TRUE,
        dots = dots
      )
      expected <- survreg_predict(
        newdata = reference_data,
        type = prediction_type,
        se.fit = TRUE,
        dots = dots
      )
      expect_native_prediction_equal(
        got,
        expected,
        tolerance = tolerance,
        info = paste(info, context, prediction_type)
      )
      if (!is_training && identical(prediction_type, "link")) {
        newdata_link <- got
      }
    }
  }

  # Construct the independent location-scale prediction from the reference
  # design. Scale parameters are nuisance parameters, so only the regression
  # coefficient covariance block contributes to the location SE.
  reference_terms <- stats::delete.response(stats::terms(fit_survreg))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = reference_newdata,
    xlev = fit_survreg$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = fit_survreg$contrasts
  )
  manual_x <- manual_x[, names(beta_survreg), drop = FALSE]
  coefficient_covariance <- covariance_survreg[
    names(beta_survreg),
    names(beta_survreg),
    drop = FALSE
  ]
  # Match predict.survreg() exactly: for supplied newdata its location
  # prediction is constructed from the model matrix and fitted coefficients.
  # The independent oracle therefore uses X beta; offset parity is still tested
  # directly above by comparing mfp2 with the native method on training and new
  # rows.
  manual_link <- as.numeric(manual_x %*% beta_survreg)
  manual_variance <- unname(rowSums(
    (manual_x %*% coefficient_covariance) * manual_x
  ))
  manual_link_se <- sqrt(pmax(manual_variance, 0))
  distribution <- survival::survreg.distributions[[fit_survreg$dist]]
  inverse_transform <- distribution$itrans
  if (is.null(inverse_transform)) inverse_transform <- identity
  manual_response <- as.numeric(inverse_transform(manual_link))

  expect_equal(
    unname(newdata_link$fit),
    manual_link,
    tolerance = tolerance,
    info = paste(info, "manual link")
  )
  expect_equal(
    unname(newdata_link$se.fit),
    manual_link_se,
    tolerance = tolerance,
    info = paste(info, "manual link SE")
  )
  expect_equal(
    unname(mfp2_predict(
      newdata = mfp2_newdata,
      type = "response",
      se.fit = FALSE
    )),
    manual_response,
    tolerance = tolerance,
    info = paste(info, "manual response")
  )
  expect_equal(
    unname(mfp2_predict(
      newdata = mfp2_newdata,
      type = "link",
      se.fit = FALSE
    )),
    unname(newdata_link$fit),
    tolerance = tolerance,
    info = paste(info, "default link target")
  )
  expect_equal(
    unname(mfp2_predict(
      newdata = mfp2_newdata,
      type = "lp",
      se.fit = FALSE
    )),
    unname(newdata_link$fit),
    tolerance = tolerance,
    info = paste(info, "lp alias")
  )
  default_arguments <- list(object = fit_mfp2, newdata = mfp2_newdata)
  default_arguments <- c(default_arguments, mfp2_predict_args)
  expect_equal(
    unname(do.call(stats::predict, default_arguments)),
    unname(newdata_link$fit),
    tolerance = tolerance,
    info = paste(info, "default prediction")
  )

  expect_linear_term_and_contrast_equal(
    object = fit_mfp2,
    newdata = term_newdata,
    term = "x1",
    ref = 0.25,
    has_intercept = TRUE,
    tolerance = tolerance
  )
}


make_survreg_linear_equivalence_data <- function(n = 360L, seed = 4102L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    scale_group = factor(rep(c("low", "high"), length.out = n)),
    off = stats::rnorm(n, mean = 0.05, sd = 0.12),
    case_weight = stats::runif(n, 0.8, 1.2)
  )
  group_effect <- c(A = 0, B = 0.25, C = -0.20)[as.character(dat$group)]
  log_event_time <- 1.8 + 0.30 * dat$x1 - 0.18 * dat$x2 +
    group_effect + dat$off + stats::rnorm(n, sd = 0.42)
  event_time <- exp(log_event_time)
  censor_time <- stats::rexp(n, rate = 0.025)
  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  dat
}


test_that("df = 1 survreg fits match every native survival distribution completely", {
  dat <- make_survreg_linear_equivalence_data()
  newdata <- dat[
    seq_len(32L),
    c("x1", "x2", "group", "off"),
    drop = FALSE
  ]

  for (dist_name in names(survival::survreg.distributions)) {
    fit_mfp2 <- mfp2(
      survival::Surv(time, status) ~ x1 + x2 + group + offset(off),
      data = dat,
      family = survreg_family(dist = dist_name),
      weights = case_weight,
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
    fit_survreg <- survival::survreg(
      survival::Surv(time, status) ~ x1 + x2 + group + offset(off),
      data = dat,
      weights = case_weight,
      dist = dist_name,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )

    expect_equal(
      sort(get_selected_variable_names(fit_mfp2)),
      sort(c("x1", "x2", "group")),
      info = dist_name
    )
    expect_complete_survreg_linear_equivalence(
      fit_mfp2,
      fit_survreg,
      reference_newdata = newdata,
      info = dist_name,
      tolerance = 1e-7
    )
  }
})


test_that("df = 1 survreg formula fit with formula scale strata and data-mask offset matches survreg completely", {
  dat <- make_survreg_linear_equivalence_data(seed = 4103L)
  rows <- c(2L, 9L, 21L, 44L, 87L, 133L, 201L, 288L, 344L)
  newdata <- dat[
    rows,
    c("x1", "x2", "group", "scale_group", "off"),
    drop = FALSE
  ]

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 + group + strata(scale_group),
    data = dat,
    family = survreg_family(dist = "weibull"),
    weights = case_weight,
    offset = off,
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
  fit_survreg <- survival::survreg(
    survival::Surv(time, status) ~
      x1 + x2 + group + strata(scale_group) + offset(off),
    data = dat,
    weights = case_weight,
    dist = "weibull",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  expect_setequal(
    get_selected_variable_names(fit_mfp2),
    c("x1", "x2", "group")
  )
  expect_complete_survreg_linear_equivalence(
    fit_mfp2,
    fit_survreg,
    reference_newdata = newdata,
    mfp2_predict_args = list(
      newoffset = newdata$off
    ),
    info = "formula scale strata and data-mask offset",
    tolerance = 1e-7
  )
})


test_that("df = 1 survreg matrix fit with fixed scale matches survreg completely", {
  dat <- make_survreg_linear_equivalence_data(seed = 4104L)
  x <- as.matrix(dat[, c("x1", "x2"), drop = FALSE])
  y <- survival::Surv(dat$time, dat$status)
  rows <- c(4L, 15L, 39L, 73L, 119L, 184L, 253L, 327L)
  newx <- x[rows, , drop = FALSE]
  reference_newdata <- dat[rows, c("x1", "x2", "off"), drop = FALSE]

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = survreg_family(dist = "weibull", scale = 0.75),
    weights = dat$case_weight,
    offset = dat$off,
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
  fit_survreg <- survival::survreg(
    survival::Surv(time, status) ~ x1 + x2 + offset(off),
    data = dat,
    weights = case_weight,
    dist = "weibull",
    scale = 0.75,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  expect_setequal(
    get_selected_variable_names(fit_mfp2),
    c("x1", "x2")
  )
  expect_complete_survreg_linear_equivalence(
    fit_mfp2,
    fit_survreg,
    reference_newdata = reference_newdata,
    mfp2_newdata = newx,
    mfp2_predict_args = list(newoffset = dat$off[rows]),
    term_newdata = as.data.frame(newx, check.names = FALSE),
    info = "fixed-scale matrix",
    tolerance = 1e-7
  )
})
