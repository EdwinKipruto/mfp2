# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: 8.1.1 Gaussian GLM equivalence without offset.
test_that("8.1.1 Gaussian: mfp2 predictions match glm without offset", {
  set.seed(8011)

  dat <- data.frame(
    y = rnorm(160),
    x1 = runif(160, 1, 5),
    x2 = rnorm(160)
  )

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x1 + x2,
    family_name = "gaussian",
    newdata_cols = c("x1", "x2")
  )
})


# Test purpose: 8.1.2 Gaussian GLM equivalence with a formula offset.
test_that("8.1.2 Gaussian: mfp2 predictions match glm with formula offset", {
  set.seed(8012)

  dat <- data.frame(
    x1 = runif(160, 1, 5),
    x2 = rnorm(160),
    off = rnorm(160, mean = 0.2, sd = 0.1)
  )
  dat$y <- 0.5 + 0.8 * dat$x1 - 0.4 * dat$x2 + dat$off + rnorm(160, sd = 0.5)

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x1 + x2 + offset(off),
    family_name = "gaussian",
    newdata_cols = c("x1", "x2", "off")
  )
})


# Test purpose: 8.1.3 Binomial GLM equivalence without offset.
test_that("8.1.3 Binomial: mfp2 predictions match glm without offset", {
  set.seed(8013)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220)
  )
  eta <- -0.6 + 0.35 * dat$x1 - 0.55 * dat$x2
  dat$y <- stats::rbinom(nrow(dat), size = 1, prob = stats::plogis(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x1 + x2,
    family_name = "binomial",
    newdata_cols = c("x1", "x2")
  )
})


# Test purpose: 8.1.4 Binomial GLM equivalence with a formula offset.
test_that("8.1.4 Binomial: mfp2 predictions match glm with formula offset", {
  set.seed(8014)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220),
    off = rnorm(220, mean = 0.1, sd = 0.2)
  )
  eta <- -0.6 + 0.35 * dat$x1 - 0.55 * dat$x2 + dat$off
  dat$y <- stats::rbinom(nrow(dat), size = 1, prob = stats::plogis(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x1 + x2 + offset(off),
    family_name = "binomial",
    newdata_cols = c("x1", "x2", "off")
  )
})


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


# Test purpose: 8.1.6 Grouped-binomial formula models with an offset should
# match stats::glm() in coefficients, covariance, likelihood, fitted values,
# predictions, and standard errors. The offset is also added manually to X beta
# so the test independently verifies the formula-offset prediction contract.
test_that("8.1.6 Binomial matrix response: mfp2 matches glm and manual calculation with formula offset", {
  set.seed(8016)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220),
    off = rnorm(220, mean = 0.1, sd = 0.15),
    trials = sample(8:25, 220, replace = TRUE)
  )
  eta <- -0.7 + 0.30 * dat$x1 - 0.45 * dat$x2 + dat$off
  dat$successes <- stats::rbinom(
    nrow(dat),
    size = dat$trials,
    prob = stats::plogis(eta)
  )
  dat$failures <- dat$trials - dat$successes

  fit_mfp2 <- mfp2(
    cbind(successes, failures) ~ x1 + x2 + offset(off),
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
    cbind(successes, failures) ~ x1 + x2 + offset(off),
    data = dat,
    family = stats::binomial()
  )

  expect_mfp2_glm_parameters_equal(fit_mfp2, fit_glm, tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit_mfp2)), as.numeric(stats::logLik(fit_glm)), tolerance = 1e-8)
  expect_equal(unname(stats::fitted(fit_mfp2)), unname(stats::fitted(fit_glm)), tolerance = 1e-8)

  nd <- dat[1:25, c("x1", "x2", "off"), drop = FALSE]

  pred_mfp2_link <- predict(fit_mfp2, newdata = nd, type = "link", se.fit = TRUE)
  pred_glm_link <- predict(fit_glm, newdata = nd, type = "link", se.fit = TRUE)
  pred_mfp2_response <- predict(fit_mfp2, newdata = nd, type = "response")
  pred_glm_response <- predict(fit_glm, newdata = nd, type = "response")

  expect_equal(as.numeric(pred_mfp2_link$fit), as.numeric(pred_glm_link$fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), as.numeric(pred_glm_link$se.fit), tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_response), as.numeric(pred_glm_response), tolerance = 1e-8)

  # The offset is fixed, so it changes the linear predictor but contributes no
  # coefficient uncertainty. Therefore the manual variance uses X V X' only.
  manual_x <- stats::model.matrix(~ x1 + x2, data = nd)
  beta <- stats::coef(fit_mfp2)
  beta_vcov <- stats::vcov(fit_mfp2)

  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_link <- as.numeric(manual_x %*% beta + nd$off)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  manual_response <- stats::plogis(manual_link)

  expect_equal(as.numeric(pred_mfp2_link$fit), manual_link, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_link$se.fit), manual_se, tolerance = 1e-10)
  expect_equal(as.numeric(pred_mfp2_response), manual_response, tolerance = 1e-10)
})


# Test purpose: 8.1.7 Poisson GLM equivalence without offset.
test_that("8.1.7 Poisson: mfp2 predictions match glm without offset", {
  set.seed(8017)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220)
  )
  eta <- 0.15 + 0.12 * dat$x1 - 0.25 * dat$x2
  dat$y <- stats::rpois(nrow(dat), lambda = exp(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x1 + x2,
    family_name = "poisson",
    newdata_cols = c("x1", "x2")
  )
})


# Test purpose: 8.1.8 Poisson GLM equivalence with a formula offset expression.
test_that("8.1.8 Poisson: mfp2 predictions match glm with formula offset", {
  set.seed(8018)

  dat <- data.frame(
    x1 = runif(220, 1, 5),
    x2 = rnorm(220),
    exposure = runif(220, 0.5, 3)
  )
  eta <- 0.15 + 0.12 * dat$x1 - 0.25 * dat$x2 + log(dat$exposure)
  dat$y <- stats::rpois(nrow(dat), lambda = exp(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x1 + x2 + offset(log(exposure)),
    family_name = "poisson",
    newdata_cols = c("x1", "x2", "exposure")
  )
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


# Test purpose: 8.1.14 Formula-interface factor expansion should match glm()
# when all terms are forced linear and retained.
test_that("8.1.14 Gaussian categorical predictor: coefficients, LPs, SEs, and predictions match glm", {
  set.seed(8024)
  n <- 180

  dat <- data.frame(
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  group_effect <- c(A = 0, B = 0.6, C = -0.4)[as.character(dat$group)]
  dat$y <- 0.5 + 0.3 * dat$x + group_effect + rnorm(n, sd = 0.3)

  fit_mfp2 <- mfp2(
    y ~ x + group,
    data = dat,
    family = "gaussian",
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
    y ~ x + group,
    data = dat,
    family = stats::gaussian()
  )

  nd <- dat[1:25, c("x", "group"), drop = FALSE]

  # The manual model matrix verifies the treatment-contrast dummy columns and
  # their coefficient ordering, not only the final delegated predictions.
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


# Test purpose: 8.1.16 Simple Gaussian coefficient and log-likelihood
# equivalence when mfp2() is forced to the same linear model as glm().
test_that("8.1.16 Gaussian: mfp2 coefficients and logLik match glm", {
  set.seed(8026)

  dat <- data.frame(
    y = rnorm(160),
    x1 = runif(160, 1, 5),
    x2 = rnorm(160)
  )

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = "gaussian",
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
    family = stats::gaussian()
  )

  # Use the complete training data as the prediction set so this test also
  # verifies the design matrix and direct X beta calculation.
  nd <- dat[, c("x1", "x2"), drop = FALSE]
  expect_glm_objects_and_manual_prediction_equal(fit_mfp2, fit_glm, nd)
})


# Test purpose: 8.1.17 A retained three-level categorical predictor should
# produce the same binomial coefficients, coefficient standard errors,
# training linear predictors, newdata link/response predictions, and
# prediction standard errors as stats::glm().
test_that("8.1.17 Binomial categorical predictor: coefficients, LPs, SEs, and predictions match glm", {
  set.seed(8017)
  n <- 420

  dat <- data.frame(
    x = stats::rnorm(n),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
  group_effect <- c(A = 0, B = 0.65, C = -0.45)[as.character(dat$group)]
  eta <- -0.35 + 0.55 * dat$x + group_effect
  dat$y <- stats::rbinom(n, size = 1, prob = stats::plogis(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x + group,
    family_name = "binomial",
    newdata_cols = c("x", "group")
  )
})


# Test purpose: 8.1.18 A retained three-level categorical predictor should
# produce the same Poisson coefficients, coefficient standard errors, training
# linear predictors, newdata link/response predictions, and prediction
# standard errors as stats::glm().
test_that("8.1.18 Poisson categorical predictor: coefficients, LPs, SEs, and predictions match glm", {
  set.seed(8018)
  n <- 360

  dat <- data.frame(
    x = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
  group_effect <- c(A = 0, B = 0.30, C = -0.25)[as.character(dat$group)]
  eta <- 0.40 + 0.35 * dat$x + group_effect
  dat$y <- stats::rpois(n, lambda = exp(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x + group,
    family_name = "poisson",
    newdata_cols = c("x", "group")
  )
})


# Test purpose: 8.1.19 A Gaussian model containing a retained categorical
# predictor and a formula offset should match stats::glm() in coefficients,
# coefficient standard errors, training linear predictors, newdata link and
# response predictions, and prediction standard errors.
test_that("8.1.19 Gaussian categorical predictor with offset matches glm completely", {
  set.seed(8119)
  n <- 260

  dat <- data.frame(
    x = stats::runif(n, -1, 2),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    off = stats::rnorm(n, mean = 0.15, sd = 0.20)
  )
  group_effect <- c(A = 0, B = 0.70, C = -0.45)[as.character(dat$group)]
  dat$y <- 0.6 + 0.35 * dat$x + group_effect + dat$off +
    stats::rnorm(n, sd = 0.35)

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x + group + offset(off),
    family_name = "gaussian",
    newdata_cols = c("x", "group", "off")
  )
})


# Test purpose: 8.1.20 A binomial model containing a retained categorical
# predictor and a formula offset should match stats::glm() for estimates,
# linear predictors, both prediction scales, and standard errors.
test_that("8.1.20 Binomial categorical predictor with offset matches glm completely", {
  set.seed(8120)
  n <- 440

  dat <- data.frame(
    x = stats::rnorm(n),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    off = stats::rnorm(n, mean = 0, sd = 0.25)
  )
  group_effect <- c(A = 0, B = 0.55, C = -0.50)[as.character(dat$group)]
  eta <- -0.35 + 0.50 * dat$x + group_effect + dat$off
  dat$y <- stats::rbinom(n, size = 1, prob = stats::plogis(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x + group + offset(off),
    family_name = "binomial",
    newdata_cols = c("x", "group", "off")
  )
})


# Test purpose: 8.1.21 A Poisson model containing a retained categorical
# predictor and a log-exposure offset should match stats::glm() completely.
test_that("8.1.21 Poisson categorical predictor with offset matches glm completely", {
  set.seed(8121)
  n <- 380

  dat <- data.frame(
    x = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    exposure = stats::runif(n, 0.4, 3.0)
  )
  group_effect <- c(A = 0, B = 0.30, C = -0.25)[as.character(dat$group)]
  eta <- 0.30 + 0.32 * dat$x + group_effect + log(dat$exposure)
  dat$y <- stats::rpois(n, lambda = exp(eta))

  expect_mfp2_glm_predictions_equal(
    dat = dat,
    formula = y ~ x + group + offset(log(exposure)),
    family_name = "poisson",
    newdata_cols = c("x", "group", "exposure")
  )
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
