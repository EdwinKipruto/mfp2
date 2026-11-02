# Complete forced-linear Fine--Gray equivalence -------------------------------
#
# survival implements Fine--Gray regression as finegray() data expansion plus a
# weighted, subject-clustered Cox model. These tests construct that native
# composition independently and compare every relative prediction exposed by
# mfp2, including SEs and all reference scales. Absolute Cox predictions are not
# tested because mfp2 deliberately does not assign them a competing-risks
# interpretation.


make_finegray_linear_equivalence_data <- function(n = 300L, seed = 4301L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    off = stats::runif(n, -0.2, 0.2),
    case_weight = sample(c(1, 2), n, replace = TRUE),
    row_id = seq_len(n)
  )
  group_effect <- c(A = 0, B = 0.40, C = -0.30)[as.character(dat$group)]
  relapse_time <- stats::rexp(
    n,
    rate = 0.055 * exp(
      0.42 * dat$x1 - 0.25 * dat$x2 + group_effect + dat$off
    )
  )
  death_time <- stats::rexp(
    n,
    rate = 0.035 * exp(-0.18 * dat$x1 + 0.12 * dat$x2)
  )
  censor_time <- stats::rexp(n, rate = 0.030)
  raw_time <- pmin(relapse_time, death_time, censor_time)
  dat$time <- pmax(round(raw_time, digits = 1L), 0.1)
  dat$event <- factor(
    ifelse(
      relapse_time <= death_time & relapse_time <= censor_time,
      "relapse",
      ifelse(death_time <= censor_time, "death", "censor")
    ),
    levels = c("censor", "relapse", "death")
  )
  dat
}


expect_complete_finegray_linear_equivalence <- function(
    fit_mfp2,
    fit_reference,
    expanded,
    original_data,
    reference_newdata,
    mfp2_newdata,
    mfp2_predict_args = list(),
    term_newdata = mfp2_newdata,
    info,
    tolerance = 1e-7) {
  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "coxph")
  expect_identical(fit_mfp2$family_string, "finegray", info = info)
  expect_identical(fit_mfp2$fitter, "base", info = info)
  expect_identical(fit_mfp2$family$etype, "relapse", info = info)
  expect_null(fit_mfp2$family$prepared, info = info)
  expect_identical(fit_mfp2$method, fit_reference$method, info = info)
  expect_true(isTRUE(fit_mfp2$convergence_mfp), info = info)
  expect_false(is.null(fit_mfp2$naive.var), info = info)
  expect_false(is.null(fit_reference$naive.var), info = info)

  row_map <- as.integer(expanded$row_id)
  n_original <- nrow(original_data)
  first_pseudo_row <- !duplicated(row_map)
  collapse_reference <- function(value) {
    result <- rep(NA_real_, n_original)
    result[row_map[first_pseudo_row]] <- as.numeric(value[first_pseudo_row])
    result
  }
  collapse_prediction <- function(value) {
    if (is.list(value)) {
      value$fit <- collapse_reference(value$fit)
      value$se.fit <- collapse_reference(value$se.fit)
      return(value)
    }
    collapse_reference(value)
  }

  expect_equal(fit_mfp2$mfp2_finegray_row_map, row_map, info = info)
  expect_equal(fit_mfp2$mfp2_finegray_n_original, n_original, info = info)
  expect_equal(fit_mfp2$nobs, n_original, info = info)
  expect_equal(
    fit_mfp2$mfp2_finegray_event,
    attr(expanded, "event", exact = TRUE),
    info = info
  )

  beta_mfp2 <- stats::coef(fit_mfp2)
  beta_reference <- stats::coef(fit_reference)
  covariance_mfp2 <- stats::vcov(fit_mfp2)
  covariance_reference <- stats::vcov(fit_reference)
  expect_equal(length(beta_mfp2), length(beta_reference), info = info)
  expect_equal(
    unname(beta_mfp2),
    unname(beta_reference),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(covariance_mfp2),
    unname(covariance_reference),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$naive.var),
    unname(fit_reference$naive.var),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(sqrt(diag(covariance_mfp2))),
    unname(sqrt(diag(covariance_reference))),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$loglik),
    unname(fit_reference$loglik),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_reference$linear.predictors),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$residuals),
    unname(fit_reference$residuals),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$means),
    unname(fit_reference$means),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$weights),
    unname(fit_reference$weights),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$n, fit_reference$n, info = info)
  expect_equal(fit_mfp2$nevent, fit_reference$nevent, info = info)
  expect_equal(fit_mfp2$iter, fit_reference$iter, info = info)
  expect_equal(
    unname(fit_mfp2$concordance),
    unname(fit_reference$concordance),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$nevents,
    sum(original_data$event == "relapse"),
    info = info
  )
  expect_equal(
    unname(fit_mfp2$score),
    unname(fit_reference$score),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$wald.test),
    unname(fit_reference$wald.test),
    tolerance = tolerance,
    info = info
  )

  for (residual_type in c(
    "martingale", "deviance", "score", "schoenfeld",
    "dfbeta", "dfbetas", "scaledsch"
  )) {
    expect_equal(
      unname(stats::residuals(fit_mfp2, type = residual_type)),
      unname(stats::residuals(fit_reference, type = residual_type)),
      tolerance = tolerance,
      info = paste(info, residual_type, "residuals")
    )
  }

  loglik_mfp2 <- stats::logLik(fit_mfp2)
  loglik_reference <- stats::logLik(fit_reference)
  expect_equal(
    as.numeric(loglik_mfp2),
    as.numeric(loglik_reference),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    attr(loglik_mfp2, "df"),
    attr(loglik_reference, "df"),
    info = info
  )
  expect_equal(
    attr(loglik_mfp2, "nobs"),
    attr(loglik_reference, "nobs"),
    info = info
  )
  expect_equal(stats::nobs(fit_mfp2), stats::nobs(fit_reference), info = info)
  expect_equal(stats::AIC(fit_mfp2), stats::AIC(fit_reference), tolerance = tolerance, info = info)
  expect_equal(
    fit_mfp2$null_logl,
    unname(fit_reference$loglik[1L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_logl,
    as.numeric(loglik_reference),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$mfp_df, attr(loglik_reference, "df"), info = info)
  expect_equal(
    fit_mfp2$null_deviance,
    -2 * unname(fit_reference$loglik[1L]),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_deviance,
    -2 * as.numeric(loglik_reference),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$linear_logl,
    as.numeric(loglik_reference),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$linear_df, attr(loglik_reference, "df"), info = info)
  expect_equal(
    fit_mfp2$linear_deviance,
    -2 * as.numeric(loglik_reference),
    tolerance = tolerance,
    info = info
  )

  expect_equal(
    fit_mfp2$mfp2_original_offset,
    original_data$off,
    tolerance = tolerance,
    info = info
  )
  expanded_offset <- expanded$off
  expect_equal(
    fit_mfp2$cox_offset_reference,
    mean(expanded_offset),
    tolerance = tolerance,
    info = info
  )

  mfp2_predict <- function(newdata = NULL, type, se.fit, reference = NULL) {
    arguments <- list(object = fit_mfp2, type = type, se.fit = se.fit)
    if (!is.null(newdata)) {
      arguments$newdata <- newdata
      arguments <- c(arguments, mfp2_predict_args)
    }
    if (!is.null(reference)) arguments$cox_reference <- reference
    do.call(stats::predict, arguments)
  }

  newdata_predictions <- list()
  for (reference in c("zero", "sample", "strata")) {
    for (type in c("lp", "risk")) {
      got_training <- mfp2_predict(
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      expected_training <- collapse_prediction(stats::predict(
        fit_reference,
        type = type,
        se.fit = TRUE,
        reference = reference
      ))
      expect_native_prediction_equal(
        got_training,
        expected_training,
        tolerance = tolerance,
        info = paste(info, "training", type, reference)
      )

      got_newdata <- mfp2_predict(
        newdata = mfp2_newdata,
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      expected_newdata <- stats::predict(
        fit_reference,
        newdata = reference_newdata,
        type = type,
        se.fit = TRUE,
        reference = reference
      )
      expect_native_prediction_equal(
        got_newdata,
        expected_newdata,
        tolerance = tolerance,
        info = paste(info, "newdata", type, reference)
      )
      if (identical(type, "risk")) {
        expect_equal(
          unname(got_newdata$fit),
          exp(unname(newdata_predictions[[paste("lp", reference, sep = ":")]]$fit)),
          tolerance = tolerance,
          info = paste(info, "risk identity", reference)
        )
      }
      newdata_predictions[[paste(type, reference, sep = ":")]] <- got_newdata
    }
  }

  # Independent X beta and X V X' calculations on the original new rows.
  manual_x <- stats::model.matrix(fit_reference, data = reference_newdata)
  if ("(Intercept)" %in% colnames(manual_x)) {
    manual_x <- manual_x[, colnames(manual_x) != "(Intercept)", drop = FALSE]
  }
  manual_x <- manual_x[, names(beta_reference), drop = FALSE]
  reference_terms <- stats::delete.response(stats::terms(fit_reference))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = reference_newdata,
    xlev = fit_reference$xlevels,
    na.action = stats::na.pass
  )
  new_offset <- stats::model.offset(reference_frame)
  training_offset <- stats::model.offset(stats::model.frame(fit_reference))
  offset_origin <- if (is.null(training_offset)) 0 else mean(training_offset)
  manual_offset <- if (is.null(new_offset)) {
    rep(0, nrow(manual_x))
  } else {
    as.numeric(new_offset - offset_origin)
  }

  for (reference in c("zero", "sample", "strata")) {
    center <- if (identical(reference, "zero")) {
      rep(0, ncol(manual_x))
    } else {
      unname(fit_reference$means[names(beta_reference)])
    }
    centered_x <- sweep(manual_x, 2L, center, FUN = "-")
    manual_lp <- as.numeric(centered_x %*% beta_reference + manual_offset)
    manual_variance <- unname(rowSums(
      (centered_x %*% covariance_reference) * centered_x
    ))
    manual_lp_se <- sqrt(pmax(manual_variance, 0))
    manual_risk <- exp(manual_lp)
    manual_risk_se <- sqrt(manual_risk) * manual_lp_se
    got_lp <- newdata_predictions[[paste("lp", reference, sep = ":")]]
    got_risk <- newdata_predictions[[paste("risk", reference, sep = ":")]]

    expect_equal(unname(got_lp$fit), manual_lp, tolerance = tolerance, info = info)
    expect_equal(unname(got_lp$se.fit), manual_lp_se, tolerance = tolerance, info = info)
    expect_equal(unname(got_risk$fit), manual_risk, tolerance = tolerance, info = info)
    expect_equal(
      unname(got_risk$se.fit),
      manual_risk_se,
      tolerance = tolerance,
      info = info
    )
  }

  expect_equal(
    unname(mfp2_predict(newdata = mfp2_newdata, type = "lp", se.fit = FALSE)),
    unname(newdata_predictions[["lp:zero"]]$fit),
    tolerance = tolerance,
    info = paste(info, "default lp")
  )
  expect_equal(
    unname(mfp2_predict(newdata = mfp2_newdata, type = "link", se.fit = FALSE)),
    unname(newdata_predictions[["lp:zero"]]$fit),
    tolerance = tolerance,
    info = paste(info, "link alias")
  )
  default_arguments <- list(object = fit_mfp2, newdata = mfp2_newdata)
  default_arguments <- c(default_arguments, mfp2_predict_args)
  expect_equal(
    unname(do.call(stats::predict, default_arguments)),
    unname(newdata_predictions[["lp:zero"]]$fit),
    tolerance = tolerance,
    info = paste(info, "default prediction")
  )

  expect_linear_term_and_contrast_equal(
    object = fit_mfp2,
    newdata = term_newdata,
    term = "x1",
    ref = 0.25,
    has_intercept = FALSE,
    tolerance = tolerance
  )
}


test_that("df = 1 Fine--Gray formula fits match the native composition for both ties", {
  dat <- make_finegray_linear_equivalence_data()
  expect_true(anyDuplicated(dat$time[dat$event == "relapse"]) > 0L)
  rows <- seq_len(34L)
  mfp2_newdata <- dat[rows, c("x1", "x2", "group", "off"), drop = FALSE]

  for (ties in c("breslow", "efron")) {
    fit_mfp2 <- mfp2(
      survival::Surv(time, event) ~ x1 + x2 + group + offset(off),
      data = dat,
      family = finegray_family(etype = "relapse"),
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
    expanded <- survival::finegray(
      survival::Surv(time, event) ~ x1 + x2 + group + off + row_id,
      data = dat,
      weights = case_weight,
      etype = "relapse"
    )
    fit_reference <- survival::coxph(
      survival::Surv(fgstart, fgstop, fgstatus) ~
        x1 + x2 + group + offset(off) + cluster(row_id),
      data = expanded,
      weights = fgwt,
      ties = ties,
      nocenter = NULL,
      robust = TRUE,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )
    reference_newdata <- mfp2_newdata
    reference_newdata$row_id <- nrow(dat) + seq_len(nrow(reference_newdata))

    expect_equal(
      sort(get_selected_variable_names(fit_mfp2)),
      sort(c("x1", "x2", "group")),
      info = ties
    )
    expect_complete_finegray_linear_equivalence(
      fit_mfp2,
      fit_reference,
      expanded = expanded,
      original_data = dat,
      reference_newdata = reference_newdata,
      mfp2_newdata = mfp2_newdata,
      info = paste("formula", ties),
      tolerance = 1e-7
    )
  }
})


test_that("df = 1 Fine--Gray matrix fit matches the native composition completely", {
  dat <- make_finegray_linear_equivalence_data(seed = 4302L)
  x <- as.matrix(dat[, c("x1", "x2"), drop = FALSE])
  y <- survival::Surv(dat$time, dat$event)
  rows <- c(2L, 13L, 28L, 57L, 96L, 141L, 204L, 267L, 294L)
  newx <- x[rows, , drop = FALSE]

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = finegray_family(etype = "relapse"),
    weights = dat$case_weight,
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
  expanded <- survival::finegray(
    survival::Surv(time, event) ~ x1 + x2 + off + row_id,
    data = dat,
    weights = case_weight,
    etype = "relapse"
  )
  fit_reference <- survival::coxph(
    survival::Surv(fgstart, fgstop, fgstatus) ~
      x1 + x2 + offset(off) + cluster(row_id),
    data = expanded,
    weights = fgwt,
    ties = "breslow",
    nocenter = NULL,
    robust = TRUE,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  reference_newdata <- data.frame(
    x1 = newx[, "x1"],
    x2 = newx[, "x2"],
    off = dat$off[rows],
    row_id = nrow(dat) + seq_along(rows)
  )

  expect_setequal(
    get_selected_variable_names(fit_mfp2),
    c("x1", "x2")
  )
  expect_complete_finegray_linear_equivalence(
    fit_mfp2,
    fit_reference,
    expanded = expanded,
    original_data = dat,
    reference_newdata = reference_newdata,
    mfp2_newdata = newx,
    mfp2_predict_args = list(newoffset = dat$off[rows]),
    term_newdata = as.data.frame(newx, check.names = FALSE),
    info = "matrix",
    tolerance = 1e-7
  )
})
