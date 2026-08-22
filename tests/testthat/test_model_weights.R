# =============================================================================
# Observation-weight validation and model compatibility
# =============================================================================
#
# mfp2() and mfpi() deliberately use one strict weight contract for every
# supported family: supplied observation weights must be finite and strictly
# positive. This avoids family-specific edge cases in likelihood-based MFP
# selection, including the infinite Gaussian AIC/log-likelihood produced by
# zero prior weights in stats::glm().

library(testthat)
library(survival)
library(mfp2)


make_weight_test_cox_data <- function(seed = 2601L, n = 140L) {
  set.seed(seed)

  dat <- data.frame(
    group = factor(rep(c("A", "B"), length.out = n)),
    x1 = stats::runif(n, 0.5, 3.0),
    x2 = stats::rnorm(n)
  )

  eta <- 0.35 * dat$x1 - 0.25 * dat$x2 +
    0.30 * (dat$group == "B")
  event_time <- stats::rexp(n, rate = 0.025 * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.012)

  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  dat
}


# The shared helper enforces the same strictly-positive contract regardless of
# model family. This test calls the internal helper directly so later interface
# changes cannot silently weaken the underlying rule.
test_that("weight validator requires strictly positive weights", {
  expect_silent(
    mfp2:::validate_model_weights(
      weights = c(0.5, 1, 2),
      nobs = 3
    )
  )

  expect_error(
    mfp2:::validate_model_weights(
      weights = c(0, 1, 2),
      nobs = 3
    ),
    "strictly positive"
  )

  expect_error(
    mfp2:::validate_model_weights(
      weights = c(-1, 1, 2),
      nobs = 3
    ),
    "strictly positive"
  )
})


# Missing and infinite values are diagnosed before the positivity check so the
# reported error identifies the actual input defect.
test_that("weight validator rejects malformed observation weights", {
  expect_error(
    mfp2:::validate_model_weights(
      weights = c(1, NA_real_, 1),
      nobs = 3
    ),
    "missing"
  )

  expect_error(
    mfp2:::validate_model_weights(
      weights = c(1, Inf, 1),
      nobs = 3
    ),
    "finite"
  )

  expect_error(
    mfp2:::validate_model_weights(
      weights = c(1, 1),
      nobs = 3
    ),
    "one value per observation"
  )

  expect_error(
    mfp2:::validate_model_weights(
      weights = c("1", "2", "3"),
      nobs = 3
    ),
    "must be numeric"
  )
})


# Positive Cox weights should reach the same forced-linear model as coxph().
test_that("weighted Cox mfp2 matches coxph with positive weights", {
  dat <- make_weight_test_cox_data()
  dat$w <- stats::runif(nrow(dat), 0.4, 2.5)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2,
    data = dat,
    family = "cox",
    weights = dat$w,
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2,
    data = dat,
    weights = w,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )

  expect_equal(
    unname(stats::coef(fit_mfp2)),
    unname(stats::coef(fit_coxph)),
    tolerance = 1e-8
  )
  expect_equal(
    unname(stats::vcov(fit_mfp2)),
    unname(stats::vcov(fit_coxph)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_coxph)),
    tolerance = 1e-8
  )
})


# Positive Gaussian weights remain supported; the stricter rule only removes
# zero/negative cases that can invalidate likelihood-based model comparisons.
test_that("weighted Gaussian mfp2 matches glm with positive weights", {
  set.seed(2602L)
  n <- 120L
  dat <- data.frame(
    x1 = stats::runif(n, 0.5, 3),
    x2 = stats::rnorm(n)
  )
  dat$y <- 0.5 + 0.8 * dat$x1 - 0.3 * dat$x2 + stats::rnorm(n, sd = 0.25)
  dat$w <- stats::runif(n, 0.5, 2)

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = "gaussian",
    weights = dat$w,
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
  fit_glm <- stats::glm(
    y ~ x1 + x2,
    data = dat,
    family = stats::gaussian(),
    weights = w
  )

  expect_equal(
    unname(stats::coef(fit_mfp2)),
    unname(stats::coef(fit_glm)),
    tolerance = 1e-8
  )
})


# mfp2 rejects zero weights for both GLM and Cox models, in both public
# interfaces. The contract applies to the complete supplied vector, so a zero
# weight remains invalid even if that row would later be excluded by `subset`.
test_that("mfp2 rejects zero weights for every family and interface", {
  dat <- make_weight_test_cox_data(seed = 2603L)
  dat$y_gaussian <- 0.5 + 0.7 * dat$x1 - 0.2 * dat$x2 +
    stats::rnorm(nrow(dat), sd = 0.3)
  w <- rep(1, nrow(dat))
  w[1L] <- 0
  rows <- 2:nrow(dat)

  x <- as.matrix(dat[, c("x1", "x2")])
  y_cox <- survival::Surv(dat$time, dat$status)

  expect_error(
    mfp2(
      x,
      y_cox,
      family = "cox",
      weights = w,
      df = 1,
      select = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )

  expect_error(
    mfp2(
      survival::Surv(time, status) ~ x1 + x2,
      data = dat,
      family = "cox",
      weights = w,
      df = 1,
      select = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )

  expect_error(
    mfp2(
      y_gaussian ~ x1 + x2,
      data = dat,
      family = "gaussian",
      weights = w,
      df = 1,
      select = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )

  expect_error(
    mfp2(
      x,
      dat$y_gaussian,
      family = "gaussian",
      weights = w,
      df = 1,
      select = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )

  expect_error(
    mfp2(
      survival::Surv(time, status) ~ x1 + x2,
      data = dat,
      subset = rows,
      family = "cox",
      weights = w,
      df = 1,
      select = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )
})


# MFPI uses the same helper as mfp2. Check both Cox and Gaussian entry points so
# a future refactor cannot accidentally restore family-specific zero handling.
test_that("mfpi rejects zero weights for Cox and Gaussian models", {
  dat <- make_weight_test_cox_data(seed = 2604L)
  dat$y_gaussian <- 0.5 + 0.6 * dat$x1 - 0.15 * dat$x2 +
    stats::rnorm(nrow(dat), sd = 0.35)
  w <- rep(1, nrow(dat))
  w[8L] <- 0

  expect_error(
    mfpi(
      survival::Surv(time, status) ~ group + x1 + x2,
      data = dat,
      family = "cox",
      group_var = "group",
      cont_vars = "x1",
      cont_var_forms = c(x1 = "linear"),
      weights = w,
      flex = "flex1",
      df = 1,
      select = 1,
      alpha = 1,
      p_interact = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )

  expect_error(
    mfpi(
      y_gaussian ~ group + x1 + x2,
      data = dat,
      family = "gaussian",
      group_var = "group",
      cont_vars = "x1",
      cont_var_forms = c(x1 = "linear"),
      weights = w,
      flex = "flex1",
      df = 1,
      select = 1,
      alpha = 1,
      p_interact = 1,
      verbose = FALSE
    ),
    "strictly positive"
  )
})
