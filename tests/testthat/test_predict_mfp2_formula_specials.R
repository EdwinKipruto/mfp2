# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# -----------------------------------------------------------------------------
# 8.3 Formula-special prediction reconstruction and error paths
# -----------------------------------------------------------------------------
# These tests target formula specials that are removed from the model matrix and
# therefore must be reconstructed from raw newdata before calling predict.glm()
# or predict.coxph().

# Test purpose: 8.3.1 Formula-level Cox strata missing from newdata should fail
# with the dedicated reconstruction message.
test_that("8.3.1 Formula-level strata missing from newdata errors clearly", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]

  fit <- mfp2(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:5, c("age", "sex", "inst"), drop = FALSE]
  expect_true(all(is.finite(predict(fit, newdata = nd, type = "lp"))))

  expect_error(
    predict(
      fit,
      newdata = nd,
      strata = nd$inst,
      type = "lp"
    ),
    "formula-level strata.*original strata variable"
  )

  expect_error(
    predict(fit, newdata = dat[1:5, c("age", "sex"), drop = FALSE], type = "lp"),
    "formula-level Cox strata|strata term could not be reconstructed"
  )
})


# Test purpose: Formula fits created through the deprecated external-strata
# compatibility path have no formula strata term to reconstruct. They must
# therefore retain the separate prediction `strata` argument until that fitting
# compatibility path is removed.
test_that("8.3.1a Deprecated external formula strata retain prediction compatibility", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(
    dat[, c("time", "status", "age", "sex", "inst")]
  ), ]

  expect_warning(
    fit <- mfp2(
      survival::Surv(time, status) ~ age + sex,
      data = dat,
      family = "cox",
      strata = dat$inst,
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      ties = "breslow",
      verbose = FALSE
    ),
    "deprecated"
  )

  expect_false(prediction_has_formula_strata(fit))
  rows <- 1:5
  prediction <- predict(
    fit,
    newdata = dat[rows, c("age", "sex"), drop = FALSE],
    strata = dat$inst[rows],
    type = "lp"
  )

  expect_length(prediction, length(rows))
  expect_true(all(is.finite(prediction)))
})


# Test purpose: 8.3.2 Formula-level offset missing from newdata should fail with
# the dedicated reconstruction message.
test_that("8.3.2 Formula-level offset missing from newdata errors clearly", {
  set.seed(8032)

  dat <- data.frame(
    x1 = runif(180, 1, 5),
    x2 = rnorm(180),
    exposure = runif(180, 0.5, 3)
  )
  eta <- 0.15 + 0.12 * dat$x1 - 0.25 * dat$x2 + log(dat$exposure)
  dat$y <- stats::rpois(nrow(dat), lambda = exp(eta))

  fit <- mfp2(
    y ~ x1 + x2 + stats::offset(log(exposure)),
    data = dat,
    family = "poisson",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_error(
    predict(fit, newdata = dat[1:5, c("x1", "x2"), drop = FALSE], type = "link"),
    "formula-level offset|offset could not be reconstructed"
  )
})
