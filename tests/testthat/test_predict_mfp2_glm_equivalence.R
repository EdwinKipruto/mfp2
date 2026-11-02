# Distinct forced-linear GLM interface cases not covered by the canonical
# all-family formula and matrix contracts in test_linear_equivalence_glm.R.


# Test purpose: 8.1.5 Grouped-binomial cbind(successes, failures) models
# should reduce exactly to stats::glm() when all predictors are forced linear
# and preprocessing is disabled. Besides coefficients and predictions, this
# test verifies the covariance matrix, log-likelihood, fitted probabilities,
# and an independent manual calculation of eta = X beta and its standard error.
test_that("8.1.5 Binomial matrix response: mfp2 matches glm and manual calculation without offset", {
  set.seed(8015)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220),
    trials = sample(8:25, 220, replace = TRUE)
  )
  eta <- -0.7 + 0.30 * dat$x1 - 0.45 * dat$x2
  dat$successes <- stats::rbinom(
    nrow(dat),
    size = dat$trials,
    prob = stats::plogis(eta)
  )
  dat$failures <- dat$trials - dat$successes

  fit_mfp2 <- mfp2(
    cbind(successes, failures) ~ x1 + x2,
    data = dat,
    family = "binomial",
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
    cbind(successes, failures) ~ x1 + x2,
    data = dat,
    family = stats::binomial()
  )

  # Both fitters should estimate the same model, not merely produce similar
  # predictions on one selected set of rows.
  expect_mfp2_glm_parameters_equal(fit_mfp2, fit_glm, tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit_mfp2)), as.numeric(stats::logLik(fit_glm)), tolerance = 1e-8)
  expect_equal(unname(stats::fitted(fit_mfp2)), unname(stats::fitted(fit_glm)), tolerance = 1e-8)

  nd <- dat[1:25, c("x1", "x2"), drop = FALSE]

  pred_mfp2_link <- predict(fit_mfp2, newdata = nd, type = "link", se.fit = TRUE)
  pred_glm_link <- predict(fit_glm, newdata = nd, type = "link", se.fit = TRUE)
  pred_mfp2_response <- predict(fit_mfp2, newdata = nd, type = "response")
  pred_glm_response <- predict(fit_glm, newdata = nd, type = "response")

  expect_equal(as.numeric(pred_mfp2_link$fit), as.numeric(pred_glm_link$fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), as.numeric(pred_glm_link$se.fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_response), as.numeric(pred_glm_response), tolerance = 1e-8)

  # Independent manual calculation. model.matrix() creates the intercept and
  # linear predictor columns, but the multiplication below is performed
  # directly rather than by predict.glm().
  manual_x <- stats::model.matrix(~ x1 + x2, data = nd)
  beta <- stats::coef(fit_mfp2)
  beta_vcov <- stats::vcov(fit_mfp2)

  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_link <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  manual_response <- stats::plogis(manual_link)

  expect_equal(as.numeric(pred_mfp2_link$fit), manual_link, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), manual_se, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_response), manual_response, tolerance = 1e-10)
})


# Test purpose: 8.1.9 Poisson matrix-interface models with an explicit offset
# should match glm() in coefficients, covariance, likelihood, fitted values,
# link/response predictions, and link-scale standard errors. A separate manual
# calculation verifies eta = X beta + offset and mu = exp(eta).
test_that("8.1.9 Poisson matrix interface: mfp2 offset model matches glm and manual calculation", {
  set.seed(8019)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220),
    exposure = runif(220, 0.5, 3)
  )
  eta <- 0.15 + 0.12 * dat$x1 - 0.25 * dat$x2 + log(dat$exposure)
  dat$y <- stats::rpois(nrow(dat), lambda = exp(eta))

  x <- as.matrix(dat[, c("x1", "x2")])
  training_offset <- log(dat$exposure)

  fit_mfp2 <- mfp2(
    x = x,
    y = dat$y,
    family = "poisson",
    offset = training_offset,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  dat$log_exposure <- training_offset
  fit_glm <- stats::glm(
    y ~ x1 + x2 + offset(log_exposure),
    data = dat,
    family = stats::poisson()
  )

  expect_mfp2_glm_parameters_equal(fit_mfp2, fit_glm, tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit_mfp2)), as.numeric(stats::logLik(fit_glm)), tolerance = 1e-8)
  expect_equal(unname(stats::fitted(fit_mfp2)), unname(stats::fitted(fit_glm)), tolerance = 1e-8)

  nd <- dat[1:25, , drop = FALSE]
  newx <- as.matrix(nd[, c("x1", "x2")])
  newoffset <- log(nd$exposure)

  pred_mfp2_link <- predict(
    fit_mfp2,
    newdata = newx,
    newoffset = newoffset,
    type = "link",
    se.fit = TRUE
  )
  pred_glm_link <- predict(fit_glm, newdata = nd, type = "link", se.fit = TRUE)
  pred_mfp2_response <- predict(
    fit_mfp2,
    newdata = newx,
    newoffset = newoffset,
    type = "response"
  )
  pred_glm_response <- predict(fit_glm, newdata = nd, type = "response")

  expect_equal(as.numeric(pred_mfp2_link$fit), as.numeric(pred_glm_link$fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), as.numeric(pred_glm_link$se.fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_response), as.numeric(pred_glm_response), tolerance = 1e-8)

  # Manual offset calculation. The offset is added after X beta and has no
  # variance term because it is supplied as known data rather than estimated.
  manual_x <- cbind(`(Intercept)` = 1, newx)
  beta <- stats::coef(fit_mfp2)
  beta_vcov <- stats::vcov(fit_mfp2)

  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_link <- as.numeric(manual_x %*% beta + newoffset)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  manual_response <- exp(manual_link)

  expect_equal(as.numeric(pred_mfp2_link$fit), manual_link, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), manual_se, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_response), manual_response, tolerance = 1e-10)
})


# Test purpose: 8.1.10 The default/matrix interface must handle a grouped
# binomial cbind(successes, failures) response and an explicit offset exactly as
# glm(). This test compares model estimates and also reconstructs predictions
# manually from X beta + offset, including link-scale standard errors.
test_that("8.1.10 Binomial matrix response and offset: mfp2 matches glm and manual calculation", {
  set.seed(8020)

  dat <- data.frame(
    x1 = runif(240, 1, 5),
    x2 = rnorm(240),
    off = rnorm(240, mean = 0.1, sd = 0.15),
    trials = sample(8:25, 240, replace = TRUE)
  )
  eta <- -0.7 + 0.30 * dat$x1 - 0.45 * dat$x2 + dat$off
  dat$successes <- stats::rbinom(
    nrow(dat),
    size = dat$trials,
    prob = stats::plogis(eta)
  )
  dat$failures <- dat$trials - dat$successes

  x <- as.matrix(dat[, c("x1", "x2")])
  y <- cbind(dat$successes, dat$failures)

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "binomial",
    offset = dat$off,
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
    cbind(successes, failures) ~ x1 + x2 + offset(off),
    data = dat,
    family = stats::binomial()
  )

  expect_mfp2_glm_parameters_equal(fit_mfp2, fit_glm, tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit_mfp2)), as.numeric(stats::logLik(fit_glm)), tolerance = 1e-8)
  expect_equal(unname(stats::fitted(fit_mfp2)), unname(stats::fitted(fit_glm)), tolerance = 1e-8)

  nd <- dat[1:25, , drop = FALSE]
  newx <- as.matrix(nd[, c("x1", "x2")])
  newoffset <- nd$off

  pred_mfp2_link <- predict(
    fit_mfp2,
    newdata = newx,
    newoffset = newoffset,
    type = "link",
    se.fit = TRUE
  )
  pred_glm_link <- predict(fit_glm, newdata = nd, type = "link", se.fit = TRUE)
  pred_mfp2_response <- predict(
    fit_mfp2,
    newdata = newx,
    newoffset = newoffset,
    type = "response"
  )
  pred_glm_response <- predict(fit_glm, newdata = nd, type = "response")

  expect_equal(as.numeric(pred_mfp2_link$fit), as.numeric(pred_glm_link$fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), as.numeric(pred_glm_link$se.fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_response), as.numeric(pred_glm_response), tolerance = 1e-8)

  # Manual grouped-binomial prediction. Trial counts affect estimation but do
  # not enter the newdata linear predictor; response predictions are event
  # probabilities obtained by applying plogis() to X beta + offset.
  manual_x <- cbind(`(Intercept)` = 1, newx)
  beta <- stats::coef(fit_mfp2)
  beta_vcov <- stats::vcov(fit_mfp2)

  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_link <- as.numeric(manual_x %*% beta + newoffset)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  manual_response <- stats::plogis(manual_link)

  expect_equal(as.numeric(pred_mfp2_link$fit), manual_link, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), manual_se, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_response), manual_response, tolerance = 1e-10)
})


# Test purpose: 8.1.11 Gaussian GLM equivalence with a non-default log link.
test_that("8.1.11 Gaussian log link: mfp2 predictions match glm", {
  set.seed(8021)

  dat <- data.frame(
    x1 = runif(180, 1, 5),
    x2 = rnorm(180)
  )
  eta <- 0.2 + 0.10 * dat$x1 - 0.08 * dat$x2
  dat$y <- exp(eta + rnorm(nrow(dat), sd = 0.05))

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = stats::gaussian(link = "log"),
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
    y ~ x1 + x2,
    data = dat,
    family = stats::gaussian(link = "log")
  )

  nd <- dat[1:25, c("x1", "x2"), drop = FALSE]

  # The shared oracle checks coefficients, covariance, fitted values, both
  # prediction scales, link-scale SEs, and the manual inverse-link calculation.
  expect_glm_objects_and_manual_prediction_equal(fit_mfp2, fit_glm, nd)

})


# Test purpose: 8.1.12 Binomial GLM equivalence with a non-default probit link.
test_that("8.1.12 Binomial probit link: mfp2 predictions match glm", {
  set.seed(8022)

  dat <- data.frame(
    x1 = runif(240, 1, 5),
    x2 = rnorm(240)
  )
  eta <- -0.6 + 0.25 * dat$x1 - 0.35 * dat$x2
  dat$y <- stats::rbinom(nrow(dat), size = 1, prob = stats::pnorm(eta))

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = stats::binomial(link = "probit"),
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
    y ~ x1 + x2,
    data = dat,
    family = stats::binomial(link = "probit")
  )

  nd <- dat[1:25, c("x1", "x2"), drop = FALSE]

  # The shared oracle checks coefficients, covariance, fitted values, both
  # prediction scales, link-scale SEs, and the manual inverse-link calculation.
  expect_glm_objects_and_manual_prediction_equal(fit_mfp2, fit_glm, nd)

})


# Test purpose: 8.1.13 Poisson GLM equivalence with a non-default sqrt link.
test_that("8.1.13 Poisson sqrt link: mfp2 predictions match glm", {
  set.seed(8023)

  dat <- data.frame(
    x1 = runif(240, 1, 5),
    x2 = runif(240, 0, 2)
  )
  eta <- 1.5 + 0.12 * dat$x1 + 0.10 * dat$x2
  dat$y <- stats::rpois(nrow(dat), lambda = eta^2)

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = stats::poisson(link = "sqrt"),
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
    y ~ x1 + x2,
    data = dat,
    family = stats::poisson(link = "sqrt")
  )

  nd <- dat[1:25, c("x1", "x2"), drop = FALSE]

  # The shared oracle checks coefficients, covariance, fitted values, both
  # prediction scales, link-scale SEs, and the manual inverse-link calculation.
  expect_glm_objects_and_manual_prediction_equal(fit_mfp2, fit_glm, nd)

})


# Test purpose: 8.1.15 Namespace-qualified stats::offset() is normalized by
# mfp2() to true formula-offset semantics, matching glm() with bare offset().
test_that("8.1.15 stats::offset expression: mfp2 predictions match glm offset semantics", {
  set.seed(8015)

  n <- 220
  dat <- data.frame(
    x1 = runif(n, 1, 8),
    x2 = rnorm(n),
    exposure = runif(n, 0.5, 2.5)
  )

  eta <- 0.25 + 0.08 * dat$x1 - 0.27 * dat$x2 + log(dat$exposure)
  dat$y <- stats::rpois(n, lambda = exp(eta))

  fit_mfp2 <- mfp2(
    y ~ fp(x1, df = 1, center = FALSE) +
      fp(x2, df = 1, center = FALSE) +
      stats::offset(log(exposure)),
    data = dat,
    family = stats::poisson(),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    xorder = "original",
    center = FALSE,
    verbose = FALSE
  )

  fit_glm <- stats::glm(
    y ~ x1 + x2 + offset(log(exposure)),
    data = dat,
    family = stats::poisson()
  )

  nd <- dat[1:25, c("x1", "x2", "exposure"), drop = FALSE]

  # The manual oracle evaluates log(exposure) from raw newdata and verifies
  # eta = X beta + log(exposure), mu = exp(eta), and X V X' standard errors.
  expect_glm_objects_and_manual_prediction_equal(fit_mfp2, fit_glm, nd)
})


# Test purpose: 8.1.22 Grouped-binomial counts with a categorical predictor and
# formula offset should match stats::glm() in the same complete set of fitting
# and prediction quantities.
test_that("8.1.22 Grouped binomial categorical predictor with offset matches glm completely", {
  set.seed(8122)
  n <- 360

  dat <- data.frame(
    x = stats::rnorm(n),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    off = stats::rnorm(n, mean = 0, sd = 0.20),
    trials = sample(8:25, n, replace = TRUE)
  )
  group_effect <- c(A = 0, B = 0.50, C = -0.40)[as.character(dat$group)]
  eta <- -0.40 + 0.45 * dat$x + group_effect + dat$off
  dat$successes <- stats::rbinom(
    n,
    size = dat$trials,
    prob = stats::plogis(eta)
  )
  dat$failures <- dat$trials - dat$successes

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = cbind(successes, failures) ~ x + group + offset(off),
    family_name = "binomial",
    newdata_cols = c("x", "group", "off")
  )
})
