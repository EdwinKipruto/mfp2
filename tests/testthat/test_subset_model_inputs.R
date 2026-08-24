# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Automatic continuous preprocessing must be estimated from the
# complete formula design even when extreme rows are not retained for fitting.
test_that("24.20 formula subset preserves full-data shift and scale", {
  set.seed(2420)
  n <- 120L
  x <- c(-10000, seq(-2, 2, length.out = n - 2L), 10000)
  dat <- data.frame(
    y = 1 + 0.35 * x + rnorm(n, sd = 0.3),
    x = x
  )
  rows <- 2:(n - 1L)

  fit <- mfp2(
    y ~ x,
    data = dat,
    subset = rows,
    keep = "x",
    df = 2,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = NULL,
    scale = NULL,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  full_shift <- find_shift_factor(dat$x)
  full_scale <- find_scale_factor(dat$x + full_shift)
  retained_shift <- find_shift_factor(dat$x[rows])
  retained_scale <- find_scale_factor(dat$x[rows] + retained_shift)

  expect_equal(
    as.numeric(fit$transformations["x", "shift"]),
    full_shift,
    tolerance = 0
  )
  expect_equal(
    as.numeric(fit$transformations["x", "scale"]),
    full_scale,
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(full_shift, retained_shift)))
  expect_false(isTRUE(all.equal(full_scale, retained_scale)))
})


# Test purpose: Externally supplied weights and offsets must follow the exact
# retained-row order in both formula methods.
test_that("24.21 formula subset aligns external weights and offsets", {
  set.seed(2421)
  n <- 150L
  dat <- data.frame(
    trt = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 5),
    z = rnorm(n),
    w = runif(n, 0.5, 2.5),
    off = rnorm(n, sd = 0.2)
  )
  dat$y <- 0.8 + 0.5 * dat$x - 0.35 * dat$z + dat$off +
    rnorm(n, sd = 0.3)
  rows <- c(121:150, 1:90)

  fit_mfp2 <- mfp2(
    y ~ x + z,
    data = dat,
    subset = rows,
    weights = dat$w,
    offset = dat$off,
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  reference <- stats::glm(
    y ~ x + z,
    data = dat[rows, , drop = FALSE],
    weights = w,
    offset = off
  )

  expect_mfp2_glm_parameters_equal(fit_mfp2, reference, tolerance = 1e-8)
  expect_equal(unname(fit_mfp2$prior.weights), unname(dat$w[rows]), tolerance = 0)
  expect_equal(unname(fit_mfp2$offset), unname(dat$off[rows]), tolerance = 0)
  expect_equal(unname(fitted(fit_mfp2)), unname(fitted(reference)), tolerance = 1e-8)

  fit_mfpi <- mfpi(
    y ~ trt + x + z,
    data = dat,
    subset = rows,
    weights = dat$w,
    offset = dat$off,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "z",
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_equal(
    unname(fit_mfpi$adjustment_model$prior.weights),
    unname(dat$w[rows]),
    tolerance = 0
  )
  expect_equal(
    unname(fit_mfpi$adjustment_model$offset),
    unname(dat$off[rows]),
    tolerance = 0
  )

  interaction_fit <-
    fit_mfpi$var_winners[["x"]]$fit$test_results$interaction_model$fit
  expect_equal(
    unname(interaction_fit$prior.weights),
    unname(dat$w[rows]),
    tolerance = 0
  )
  expect_equal(
    unname(interaction_fit$offset),
    unname(dat$off[rows]),
    tolerance = 0
  )
})


# Test purpose: Formula offsets and two-column grouped-binomial responses must
# be reconstructed from the same retained rows used by the predictor matrix.
test_that("24.22 subset aligns grouped-binomial response and formula offset", {
  set.seed(2422)
  n <- 180L
  dat <- data.frame(
    x = runif(n, 1, 5),
    z = rnorm(n),
    off = rnorm(n, sd = 0.2),
    trials = sample(5:12, n, replace = TRUE)
  )
  eta <- -0.4 + 0.35 * dat$x - 0.25 * dat$z + dat$off
  dat$successes <- stats::rbinom(n, size = dat$trials, prob = stats::plogis(eta))
  dat$failures <- dat$trials - dat$successes
  dat$keep_row <- rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)

  fit <- mfp2(
    cbind(successes, failures) ~ x + z + offset(off),
    data = dat,
    subset = keep_row,
    family = stats::binomial(),
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  retained <- dat[dat$keep_row, , drop = FALSE]
  reference <- stats::glm(
    cbind(successes, failures) ~ x + z + offset(off),
    data = retained,
    family = stats::binomial()
  )

  expect_mfp2_glm_parameters_equal(fit, reference, tolerance = 1e-8)
  expect_equal(unname(fitted(fit)), unname(fitted(reference)), tolerance = 1e-8)
  expect_equal(
    unname(fit$prior.weights),
    unname(reference$prior.weights),
    tolerance = 0
  )
  expect_equal(unname(fit$offset), unname(retained$off), tolerance = 0)
  expect_equal(
    unname(predict(fit, newdata = retained, type = "link")),
    unname(predict(reference, newdata = retained, type = "link")),
    tolerance = 1e-8
  )
})


# Test purpose: External Cox strata and formula strata() must remain aligned
# after an ordered numeric subset. The MFPI formula path must reconstruct its
# strata special from retained-level metadata during prediction.
test_that("24.23 subset aligns external and formula Cox strata", {
  set.seed(2423)
  n <- 180L
  dat <- data.frame(
    trt = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 5),
    z = rnorm(n),
    stratum = factor(rep(c("S1", "S2", "S3"), length.out = n))
  )
  lp <- 0.35 * dat$x - 0.2 * dat$z + 0.25 * (dat$trt == "treated")
  event_time <- stats::rexp(n, rate = 0.03 * exp(lp))
  censor_time <- stats::rexp(n, rate = 0.018)
  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  rows <- c(121:180, 1:90)
  retained <- dat[rows, , drop = FALSE]

  fit_external <- mfp2(
    survival::Surv(time, status) ~ x + z,
    data = dat,
    subset = rows,
    family = "cox",
    strata = dat$stratum,
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_formula <- mfp2(
    survival::Surv(time, status) ~ x + z + strata(stratum),
    data = dat,
    subset = rows,
    family = "cox",
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  reference <- survival::coxph(
    survival::Surv(time, status) ~ x + z + strata(stratum),
    data = retained,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )

  expect_equal(unname(coef(fit_external)), unname(coef(reference)), tolerance = 1e-8)
  expect_equal(unname(coef(fit_formula)), unname(coef(reference)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit_external)), as.numeric(logLik(reference)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit_formula)), as.numeric(logLik(reference)), tolerance = 1e-8)

  nd <- retained[1:15, c("x", "z", "stratum"), drop = FALSE]
  expect_equal(
    unname(predict(fit_formula, newdata = nd, type = "lp", cox_reference = "zero")),
    unname(predict(reference, newdata = nd, type = "lp", reference = "zero")),
    tolerance = 1e-8
  )

  fit_mfpi <- mfpi(
    survival::Surv(time, status) ~ trt + x + z + strata(stratum),
    data = dat,
    subset = rows,
    family = "cox",
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "z",
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )
  expect_s3_class(fit_mfpi, "mfpi")
  expect_false(is.null(fit_mfpi$formula_strata_terms))

  mfpi_nd <- retained[1:12, c("trt", "x", "z", "stratum"), drop = FALSE]
  mfpi_prediction <- predict(
    fit_mfpi,
    newdata = mfpi_nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expect_true(all(is.finite(mfpi_prediction$predictions$fit)))
})


# Test purpose: A supplied numeric categorical block must remain untouched when
# the retained rows preserve its estimable dimension in both matrix interfaces.
test_that("24.24 matrix grouped terms that retain rank keep supplied coding", {
  set.seed(2424)
  n <- 180L
  stage <- factor(rep(c("A", "B", "C"), each = n / 3L))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  x_cont <- runif(n, 1, 5)
  z <- rnorm(n)
  rows <- seq_len(n) %% 5L != 0L

  x_mfp2 <- cbind(x = x_cont, stage_mm, z = z)
  y <- 1 + 0.45 * x_cont + 0.6 * stage_mm[, 1L] -
    0.25 * stage_mm[, 2L] + 0.2 * z + rnorm(n, sd = 0.3)

  expect_equal(
    qr(cbind(1, stage_mm))$rank,
    qr(cbind(1, stage_mm[rows, , drop = FALSE]))$rank
  )

  fit_mfp2 <- mfp2(
    x_mfp2,
    y,
    subset = rows,
    term_groups = list(stage = colnames(stage_mm)),
    keep = c("x", "stage", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  expect_equal(
    unname(fit_mfp2$x_original),
    unname(x_mfp2[rows, , drop = FALSE]),
    tolerance = 0
  )

  trt <- rep(c(0, 1), length.out = n)
  x_mfpi <- cbind(trt = trt, x = x_cont, stage_mm, z = z)
  fit_mfpi <- mfpi(
    x_mfpi,
    y,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = colnames(stage_mm)),
    keep = c("stage", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_s3_class(fit_mfpi, "mfpi")
  expect_equal(
    unname(fit_mfpi$x_train_internal[, colnames(stage_mm), drop = FALSE]),
    unname(stage_mm[rows, , drop = FALSE]),
    tolerance = 0
  )
})
