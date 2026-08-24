# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 23. Selection truth, SAZ decisions, weights, and serialization
# =============================================================================
# These tests target the remaining high-risk contracts: the model-selection
# decision itself, the three SAZ stage-2 representations, rejection of a false
# interaction, weighted-fit equivalence, and persistence of fitted objects.

# Test purpose: A strong linear signal should be retained as an ordinary linear
# effect, not replaced by a nonlinear FP1 candidate. The fitted model should be
# numerically equivalent to glm(y ~ x) when preprocessing is disabled.
test_that("23.1 strong linear signal is selected as linear and matches glm", {
  set.seed(2301)
  n <- 300
  x <- seq(0.5, 10, length.out = n)
  y <- 1.25 + 2.4 * x + rnorm(n, sd = 0.02)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))

  fit <- mfp2(
    xmat,
    y,
    powers = list(x = c(0, 1)),
    df = 2,
    select = 0.05,
    alpha = 0.05,
    criterion = "pvalue",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ x)

  expect_true(fit$fp_terms["x", "selected"])
  expect_equal(as.numeric(fit$fp_powers[["x"]]), 1)
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(reference)),
               tolerance = 1e-8)
})


# Test purpose: A strong logarithmic signal with candidate powers restricted to
# 0 and 1 should select power 0. This verifies that the FP1 search and closed
# testing procedure prefer the known nonlinear generating function.
test_that("23.2 strong logarithmic signal selects FP1 power zero", {
  set.seed(2302)
  n <- 350
  x <- exp(seq(log(0.4), log(20), length.out = n))
  y <- 0.8 + 3.1 * log(x) + rnorm(n, sd = 0.02)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))

  fit <- mfp2(
    xmat,
    y,
    powers = list(x = c(0, 1)),
    df = 2,
    select = 0.05,
    alpha = 0.05,
    criterion = "pvalue",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ log(x))

  expect_true(fit$fp_terms["x", "selected"])
  expect_equal(as.numeric(fit$fp_powers[["x"]]), 0)
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})


# Test purpose: With a single non-unity candidate and x named in
# force_max_fp_vars, the selected FP2 basis must be the repeated-power pair
# (2, 2), corresponding to
# x^2 and x^2 log(x). This exercises repeated-power use in the complete fitter.
test_that("23.3 forced repeated FP2 uses the expected power pair and basis", {
  set.seed(2303)
  n <- 280
  x <- seq(0.5, 6, length.out = n)
  y <- 1 + 1.8 * x^2 - 0.7 * x^2 * log(x) + rnorm(n, sd = 0.02)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))

  fit <- mfp2(
    xmat,
    y,
    powers = list(x = 2),
    df = 4,
    select = 1,
    criterion = "aic",
    force_max_fp_vars = "x",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ I(x^2) + I(x^2 * log(x)))

  expect_equal(as.numeric(fit$fp_powers[["x"]]), c(2, 2))
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})
