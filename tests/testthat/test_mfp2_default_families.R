# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 1. mfp2.default() — core fitting across families
# =============================================================================

# Test purpose: Fits the default Gaussian model and verifies core classes,
# convergence, and MFP metadata are present.
test_that("mfp2.default() returns an mfp2 object for Gaussian family", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_s3_class(fit, "glm")
  expect_equal(fit$family_string, "gaussian")
  expect_identical(fit$fitter, "base")
  expect_true(fit$convergence_mfp)
  expect_true(is.data.frame(fit$fp_terms))
  expect_true(is.list(fit$fp_powers))
  expect_true(is.data.frame(fit$transformations))
  expect_equal(nrow(fit$fp_terms), ncol(x_prostate))
  expect_equal(length(fit$fp_powers), ncol(x_prostate))
})


# Test purpose: Fits a binomial model with a numeric binary response and checks
#  the returned object metadata.
test_that("mfp2.default() returns correct object for binomial family", {
  data("pima", package = "mfp2")
  x_pima <- as.matrix(pima[, c("glucose", "bmi", "age", "pregnant")])
  y_pima <- pima$y

  fit <- mfp2(x_pima, y_pima, family = "binomial", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_s3_class(fit, "glm")
  expect_equal(fit$family_string, "binomial")
  expect_true(fit$convergence_mfp)
})


# Test purpose: Checks that a two-level factor response is accepted for binomial
#  models.
test_that("mfp2.default() works with two-level factor binomial response", {
  data("pima", package = "mfp2")
  x_pima <- as.matrix(pima[, c("glucose", "bmi", "age")])
  y_factor <- factor(pima$y, levels = c(0, 1), labels = c("no", "yes"))

  fit <- mfp2(x_pima, y_factor, family = "binomial", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_equal(fit$family_string, "binomial")
})


# Test purpose: Checks that grouped binomial counts supplied as successes/failures
# are accepted.
test_that("mfp2.default() works with grouped binomial response", {
  set.seed(42)
  n <- 100
  x <- cbind(x1 = rnorm(n, 10, 2), x2 = rnorm(n, 5, 1))
  trials <- sample(10:20, n, replace = TRUE)
  successes <- rbinom(n, trials, plogis(-2 + 0.1 * x[, 1]))
  y <- cbind(successes, trials - successes)

  fit <- mfp2(
    x,
    y,
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
  fit_glm <- stats::glm.fit(
    x = cbind(`(Intercept)` = 1, x),
    y = y,
    family = stats::binomial()
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(fit$family_string, "binomial")
  expect_length(stats::coef(fit), length(fit_glm$coefficients))
  expect_equal(
    unname(stats::coef(fit)),
    unname(fit_glm$coefficients),
    tolerance = 1e-8
  )
})


# Test purpose: Fits a Poisson model on simulated count data and checks the
# resolved family.
test_that("mfp2.default() works with Poisson family", {
  set.seed(1)
  n <- 200
  x <- cbind(x1 = runif(n, 1, 10), x2 = runif(n, 1, 5))
  y <- rpois(n, exp(0.5 + 0.1 * x[, 1]))

  fit <- mfp2(x, y, family = "poisson", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_equal(fit$family_string, "poisson")
})


# Test purpose: User-supplied glm.control() settings must reach both the fast
# candidate path and the final formula-based base GLM fit. The deliberately
# restrictive one-iteration limit should therefore trigger the package's
# fail-fast non-convergence error in both paths.
test_that("GLM control is propagated to candidate and final base fits", {
  set.seed(1200)
  n <- 100L
  x <- cbind(
    x1 = seq(-2.5, 2.5, length.out = n),
    x2 = stats::rnorm(n)
  )
  eta <- -0.3 + 1.8 * x[, "x1"] - 1.1 * x[, "x2"]
  y <- stats::rbinom(n, 1L, stats::plogis(eta))
  control <- stats::glm.control(epsilon = 1e-12, maxit = 1L)

  expect_error(
    fit_model(
      x = x,
      y = y,
      family = stats::binomial(),
      family_string = "binomial",
      control = control,
      fast = TRUE,
      keep_fit = TRUE,
      fitter = "base"
    ),
    "did not converge"
  )

  expect_error(
    fit_model(
      x = x,
      y = y,
      family = stats::binomial(),
      family_string = "binomial",
      control = control,
      fast = FALSE,
      keep_fit = TRUE,
      fitter = "base"
    ),
    "did not converge"
  )
})


# Test purpose: Control lists are normalized with glm.control() semantics so
# partial lists receive defaults and malformed controls fail with an mfp2 error.
test_that("GLM control normalization is explicit", {
  control <- normalize_glm_control(list(epsilon = 1e-7, maxit = 17L))

  expect_equal(control$epsilon, 1e-7)
  expect_identical(control$maxit, 17L)
  expect_false(control$trace)

  expect_error(
    normalize_glm_control("not-a-control-list"),
    "must be `NULL` or a list"
  )
  expect_error(
    normalize_glm_control(list(maxit = 0L)),
    "Invalid GLM `control`"
  )
})


# Test purpose: fastglm has no trace equivalent. Reject the setting before
# backend dispatch rather than silently dropping glm.control(trace = TRUE).
test_that("fastglm rejects unsupported GLM iteration tracing", {
  x <- cbind(x1 = seq_len(20L))
  y <- stats::rpois(20L, lambda = 2)

  expect_error(
    fit_model(
      x = x,
      y = y,
      family = stats::poisson(),
      family_string = "poisson",
      control = stats::glm.control(trace = TRUE),
      fast = TRUE,
      fitter = "fastglm"
    ),
    "not supported with `fitter = \"fastglm\"`"
  )
})


# Test purpose: Fits a negative-binomial model through the public matrix
# interface and checks the family-specific metadata returned by mfp2().
test_that("mfp2.default() works with negative-binomial family", {
  skip_if_not_installed("fastglm")

  dat <- make_negbin_test_data(n = 180L, seed = 1101L)
  x <- as.matrix(dat[, c("x1", "x2", "x3")])

  fit <- mfp2(
    x,
    dat$y,
    family = "negbin",
    fitter = "fastglm",
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

  expect_s3_class(fit, "mfp2")
  expect_s3_class(fit, "fastglm_nb")
  expect_identical(fit$family_string, "negbin")
  expect_identical(fit$fitter, "fastglm")
  expect_true(is.finite(fit$theta) && fit$theta > 0)
  expect_true(is.finite(fit$mfp_logl))
  expect_true(is.na(fit$null_logl))
  expect_true(fit$convergence_mfp)
  expect_true(all(fit$fp_terms[c("x1", "x2", "x3"), "selected"]))
  expect_equal(fit$mfp_df, fit$rank + 1L)
  expect_equal(fit$aic, -2 * fit$mfp_logl + 2 * fit$mfp_df,
               tolerance = 1e-10)
})


# Test purpose: Exercises the repeated FP candidate-search path for negative
# binomial rather than only the forced-linear special case.
test_that("negative-binomial fitting completes an FP candidate search", {
  skip_if_not_installed("fastglm")

  set.seed(1102)
  n <- 220L
  x <- cbind(
    x1 = stats::runif(n, 0.5, 3),
    x2 = stats::rbinom(n, 1L, 0.45)
  )
  mu <- exp(0.2 + 0.8 * log(x[, "x1"]) + 0.35 * x[, "x2"])
  y <- stats::rnbinom(n, mu = mu, size = 3.5)

  fit <- mfp2(
    x,
    y,
    family = "negbin",
    fitter = "fastglm",
    cycles = 5,
    df = c(x1 = 4, x2 = 1),
    select = 1,
    alpha = 0.05,
    xorder = "original",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_s3_class(fit, "fastglm_nb")
  expect_true(fit$convergence_mfp)
  expect_true(is.finite(fit$theta) && fit$theta > 0)
  expect_true(is.finite(fit$mfp_logl))
  expect_true(is.na(fit$null_logl))
  expect_true(fit$fp_terms["x1", "selected"])
  expect_true(length(fit$fp_powers[["x1"]]) %in% c(1L, 2L))
  expect_equal(fit$mfp_df, fit$rank + 1L)
})


# Test purpose: Fits a Cox proportional hazards model and checks Cox-specific
#  class and metadata.
test_that("mfp2.default() works with Cox family", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes", "pgr", "er")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_s3_class(fit, "coxph")
  expect_equal(fit$family_string, "cox")
  expect_true(fit$convergence_mfp)
})


# Test purpose: Checks public validation of the fitter argument.
test_that("mfp2 rejects invalid fitter values", {
  expect_error(
    mfp2(
      x_prostate,
      y_prostate,
      fitter = "unknown",
      verbose = FALSE
    ),
    "one of"
  )
})


# Test purpose: Negative binomial has no stats::glm.fit() backend capable of
# estimating theta, so requesting the base fitter must fail clearly.
test_that("negative binomial rejects the base fitter", {
  x <- cbind(x1 = seq_len(8L), x2 = rep(c(0, 1), 4L))
  y <- c(0, 1, 2, 3, 5, 2, 1, 4)

  expect_error(
    mfp2(
      x,
      y,
      family = "negbin",
      fitter = "base",
      verbose = FALSE
    ),
    "available only with.*fitter = \"fastglm\""
  )
})


# Test purpose: Confirms that NB count validation is enforced through the public
# mfp2 interface, not only by the internal validation helper.
test_that("negative binomial rejects non-integer responses publicly", {
  skip_if_not_installed("fastglm")

  x <- cbind(x1 = seq_len(8L), x2 = rep(c(0, 1), 4L))
  y <- c(0, 1, 2, 3.5, 5, 2, 1, 4)

  expect_error(
    mfp2(
      x,
      y,
      family = "negbin",
      fitter = "fastglm",
      verbose = FALSE
    ),
    "integer"
  )
})
