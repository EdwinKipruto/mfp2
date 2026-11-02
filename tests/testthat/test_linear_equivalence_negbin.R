# Complete forced-linear negative-binomial equivalence ------------------------
#
# Negative-binomial fitting is intentionally separate from base glm: mfp2 uses
# fastglm to estimate theta, while MASS::glm.nb() supplies the independent
# reference. The looser tolerances reflect the different optimizers, not a
# weaker contract. Both ordinary and weighted-offset routes check the fitted
# model, likelihood accounting, complete predictions and SEs, a direct matrix
# oracle, and package-owned term/contrast predictions.


expect_complete_negbin_linear_equivalence <- function(
    fit_mfp2,
    fit_mass,
    reference_newdata,
    mfp2_newdata = reference_newdata,
    mfp2_predict_args = list(),
    info,
    coefficient_tolerance = 5e-4,
    likelihood_tolerance = 2e-3) {
  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "fastglm_nb")
  expect_identical(fit_mfp2$family_string, "negbin", info = info)
  expect_identical(fit_mfp2$fitter, "fastglm", info = info)
  expect_identical(fit_mfp2$family$link, "log", info = info)
  expect_true(isTRUE(fit_mfp2$convergence_mfp), info = info)
  expect_true(isTRUE(fit_mfp2$converged), info = info)

  expect_negbin_mfp2_mass_equal(
    fit_mfp2,
    fit_mass,
    coefficient_tolerance = coefficient_tolerance,
    likelihood_tolerance = likelihood_tolerance
  )
  expect_equal(fit_mfp2$rank, fit_mass$rank, info = info)
  expect_equal(fit_mfp2$df.residual, fit_mass$df.residual, info = info)
  expect_equal(
    unname(sqrt(diag(stats::vcov(fit_mfp2)))),
    unname(sqrt(diag(stats::vcov(fit_mass)))),
    tolerance = coefficient_tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_mass$linear.predictors),
    tolerance = coefficient_tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$prior.weights),
    unname(fit_mass$prior.weights),
    tolerance = coefficient_tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$null.deviance,
    fit_mass$null.deviance,
    tolerance = likelihood_tolerance,
    info = info
  )
  for (residual_type in c("deviance", "pearson", "response")) {
    expect_equal(
      unname(stats::residuals(fit_mfp2, type = residual_type)),
      unname(stats::residuals(fit_mass, type = residual_type)),
      tolerance = likelihood_tolerance,
      info = paste(info, residual_type, "residuals")
    )
  }

  loglik_mfp2 <- stats::logLik(fit_mfp2)
  loglik_mass <- stats::logLik(fit_mass)
  expect_equal(
    as.numeric(loglik_mfp2),
    as.numeric(loglik_mass),
    tolerance = likelihood_tolerance,
    info = info
  )
  expect_equal(attr(loglik_mfp2, "df"), attr(loglik_mass, "df"), info = info)
  expect_equal(
    attr(loglik_mfp2, "nobs"),
    attr(loglik_mass, "nobs"),
    info = info
  )
  expect_equal(
    stats::AIC(fit_mfp2),
    stats::AIC(fit_mass),
    tolerance = likelihood_tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$twologlik / 2,
    as.numeric(loglik_mass),
    tolerance = likelihood_tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_df,
    attr(loglik_mass, "df"),
    info = info
  )
  expect_equal(
    fit_mfp2$mfp_deviance,
    stats::deviance(fit_mass),
    tolerance = likelihood_tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$linear_logl,
    fit_mfp2$mfp_logl,
    tolerance = coefficient_tolerance,
    info = info
  )
  expect_equal(fit_mfp2$linear_df, fit_mfp2$mfp_df, info = info)
  expect_equal(
    fit_mfp2$linear_deviance,
    fit_mfp2$mfp_deviance,
    tolerance = likelihood_tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$null_deviance,
    fit_mass$null.deviance,
    tolerance = likelihood_tolerance,
    info = info
  )
  expect_true(is.na(fit_mfp2$null_logl), info = info)

  mfp2_predict <- function(newdata = NULL, type, se.fit) {
    arguments <- list(object = fit_mfp2, type = type, se.fit = se.fit)
    if (!is.null(newdata)) {
      arguments$newdata <- newdata
      arguments <- c(arguments, mfp2_predict_args)
    }
    do.call(stats::predict, arguments)
  }

  for (context in c("training", "newdata")) {
    is_training <- identical(context, "training")
    mfp2_data <- if (is_training) NULL else mfp2_newdata
    mass_data <- if (is_training) NULL else reference_newdata

    for (prediction_type in c("link", "response")) {
      got <- mfp2_predict(
        newdata = mfp2_data,
        type = prediction_type,
        se.fit = TRUE
      )
      expected <- if (is_training) {
        stats::predict(
          fit_mass,
          type = prediction_type,
          se.fit = TRUE
        )
      } else {
        stats::predict(
          fit_mass,
          newdata = mass_data,
          type = prediction_type,
          se.fit = TRUE
        )
      }
      expect_native_prediction_equal(
        got,
        expected,
        tolerance = coefficient_tolerance,
        info = paste(info, context, prediction_type)
      )
    }
  }

  reference_terms <- stats::delete.response(stats::terms(fit_mass))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = reference_newdata,
    xlev = fit_mass$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = fit_mass$contrasts
  )
  beta <- stats::coef(fit_mass)
  covariance <- stats::vcov(fit_mass)
  manual_x <- manual_x[, names(beta), drop = FALSE]
  manual_offset <- stats::model.offset(reference_frame)
  if (is.null(manual_offset)) manual_offset <- rep(0, nrow(manual_x))
  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_variance <- unname(rowSums((manual_x %*% covariance) * manual_x))
  manual_link_se <- sqrt(pmax(manual_variance, 0))
  manual_response <- exp(manual_link)
  manual_response_se <- manual_response * manual_link_se

  new_link <- mfp2_predict(
    newdata = mfp2_newdata,
    type = "link",
    se.fit = TRUE
  )
  new_response <- mfp2_predict(
    newdata = mfp2_newdata,
    type = "response",
    se.fit = TRUE
  )
  expect_equal(
    unname(new_link$fit),
    manual_link,
    tolerance = coefficient_tolerance,
    info = paste(info, "manual link")
  )
  expect_equal(
    unname(new_link$se.fit),
    manual_link_se,
    tolerance = coefficient_tolerance,
    info = paste(info, "manual link SE")
  )
  expect_equal(
    unname(new_response$fit),
    manual_response,
    tolerance = coefficient_tolerance,
    info = paste(info, "manual response")
  )
  expect_equal(
    unname(new_response$se.fit),
    manual_response_se,
    tolerance = coefficient_tolerance,
    info = paste(info, "manual response SE")
  )
  expect_equal(
    unname(mfp2_predict(newdata = mfp2_newdata, type = "lp", se.fit = FALSE)),
    unname(new_link$fit),
    tolerance = 1e-10,
    info = paste(info, "lp alias")
  )
  default_arguments <- c(
    list(object = fit_mfp2, newdata = mfp2_newdata),
    mfp2_predict_args
  )
  expect_equal(
    unname(do.call(stats::predict, default_arguments)),
    unname(new_link$fit),
    tolerance = 1e-10,
    info = paste(info, "default prediction")
  )

  expect_linear_term_and_contrast_equal(
    object = fit_mfp2,
    newdata = as.data.frame(mfp2_newdata, check.names = FALSE),
    term = "x1",
    ref = 1,
    has_intercept = TRUE,
    tolerance = 1e-8
  )
}


test_that("df = 1 negative-binomial fit matches MASS completely", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  fits <- get_negbin_reference_fits()
  newdata <- fits$data[seq_len(24L), c("x1", "x2", "x3"), drop = FALSE]

  expect_setequal(
    get_selected_variable_names(fits$mfp2),
    c("x1", "x2", "x3")
  )
  expect_complete_negbin_linear_equivalence(
    fits$mfp2,
    fits$mass,
    reference_newdata = newdata,
    info = "ordinary",
    coefficient_tolerance = 1e-4,
    likelihood_tolerance = 1e-3
  )
})


test_that("df = 1 weighted-offset negative-binomial fit matches MASS completely", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  set.seed(1201)
  n <- 280L
  dat <- data.frame(
    x1 = stats::runif(n, 0.5, 2.5),
    x2 = stats::rnorm(n),
    exposure = stats::runif(n, 0.6, 2.2),
    w = sample(c(1, 2), n, replace = TRUE)
  )
  dat$log_exposure <- log(dat$exposure)
  dat$y <- stats::rnbinom(
    n,
    mu = exp(0.15 + 0.4 * dat$x1 - 0.2 * dat$x2 + dat$log_exposure),
    size = 3.2
  )

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = "negbin",
    fitter = "fastglm",
    weights = w,
    offset = log_exposure,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    xorder = "original",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  fit_mass <- MASS::glm.nb(
    y ~ x1 + x2 + offset(log_exposure),
    data = dat,
    weights = w,
    link = log,
    control = stats::glm.control(maxit = 100L)
  )
  reference_newdata <- dat[
    seq_len(20L),
    c("x1", "x2", "log_exposure"),
    drop = FALSE
  ]
  mfp2_newdata <- reference_newdata[, c("x1", "x2"), drop = FALSE]

  expect_setequal(
    get_selected_variable_names(fit_mfp2),
    c("x1", "x2")
  )
  expect_complete_negbin_linear_equivalence(
    fit_mfp2,
    fit_mass,
    reference_newdata = reference_newdata,
    mfp2_newdata = mfp2_newdata,
    mfp2_predict_args = list(newoffset = reference_newdata$log_exposure),
    info = "weighted offset"
  )

  shifted <- predict(
    fit_mfp2,
    newdata = mfp2_newdata,
    newoffset = reference_newdata$log_exposure + 0.35,
    type = "link"
  )
  baseline <- predict(
    fit_mfp2,
    newdata = mfp2_newdata,
    newoffset = reference_newdata$log_exposure,
    type = "link"
  )
  expect_equal(
    unname(shifted - baseline),
    rep(0.35, nrow(reference_newdata)),
    tolerance = 1e-10
  )
})
