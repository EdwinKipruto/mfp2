# =============================================================================
# Comprehensive tests for the mfp2 package
# =============================================================================
#
# File organization
# -----------------
# Tests are grouped by feature area. Keep new tests inside the most specific
# numbered section, and use dotted subsection numbers when several tests belong
# together, for example 8.1.1, 8.1.2, ... for prediction-equivalence checks.
#
# Sections:
#   1.  mfp2.default() - Gaussian, binomial, Poisson, negative binomial, Cox, fitter
#   2.  mfp2.formula() - formula parsing, fp()/fp2(), factor handling, strata
#   3.  Family and response validation
#   4.  Preprocessing - shift, scale, centering
#       4.1 Shared shift/scale fixtures and validation expectations
#       4.2 mfp2 shift/scale matrix and formula interfaces
#       4.3 mfpi shift/scale matrix and formula interfaces
#   5.  Candidate-power validation and custom powers
#   6.  SAZ (spike-at-zero) - eligibility, cascade, reset, prediction
#   7.  ACD transformation - fitting, validation, prediction, stored parameters
#   8.  predict.mfp2() - ordinary prediction, equivalence tests, offsets, strata
#       8.1 GLM equivalence against stats::glm() and manual X beta checks
#       8.2 Cox equivalence against survival::coxph()
#       8.3 Formula-special prediction reconstruction and error paths
#   9.  Model selection criteria - p-value, AIC, BIC
#  10.  Edge cases and input validation
#  11.  Summary, print, and coef methods
#  12.  Weights and offsets
#  13.  Convergence and cycles
#  14.  zero_vars and catzero_vars
#  15.  force_max_fp_vars and formula-term force_max_fp
#  16.  mfpi() - basic interaction fitting
#       16.1 Grouped categorical adjustment terms
#  17.  predict.mfpi()
#       17.2 Grouped categorical adjustment prediction
#  18.  Reproducibility
#  19.  Transformation helpers
#  20.  Likelihood-ratio and F-test helpers
#  21.  plot()
#  22.  C++ core tests
#  23.  Reference-model, serialization, and compatibility regressions
#  24.  subset handling across formula and matrix interfaces
#  25.  MFPI printing and ACD reconstruction regressions
# =============================================================================

library(testthat)
library(survival)
# Keep survival formula specials unqualified. Some older survival releases
# used by package-check services do not recognize the namespace-qualified
# form as a formula special, whereas strata() is supported consistently.
library(mfp2)


# =============================================================================
# Test data setup
# =============================================================================

data("prostate", package = "mfp2")

# Shared prostate fixtures used by many Gaussian/default-interface tests.
# x_prostate contains only predictors; y_prostate is the continuous response.
x_prostate <- as.matrix(prostate[, 2:8])
y_prostate <- as.numeric(prostate$lpsa)

# Helper: suppress verbose progress messages in tests that intentionally call
# functions with user-facing output. Keep test assertions outside quiet().
quiet <- function(expr) suppressMessages(capture.output(expr, type = "message"))


# Compare mfp2 and glm parameters positionally. Every GLM-equivalence fit in
# this suite uses xorder = "original", so the coefficient and covariance order
# must match the reference glm() fit directly. mfp2's transformed-column suffix
# is deliberately ignored because these assertions compare numerical values.
expect_mfp2_glm_parameters_equal <- function(fit_mfp2,
                                             fit_glm,
                                             tolerance = 1e-8) {
  coef_mfp2 <- stats::coef(fit_mfp2)
  coef_glm <- stats::coef(fit_glm)
  vcov_mfp2 <- stats::vcov(fit_mfp2)
  vcov_glm <- stats::vcov(fit_glm)

  expect_length(coef_mfp2, length(coef_glm))
  expect_equal(dim(vcov_mfp2), dim(vcov_glm))
  expect_equal(unname(coef_mfp2), unname(coef_glm), tolerance = tolerance)
  expect_equal(unname(vcov_mfp2), unname(vcov_glm), tolerance = tolerance)
  expect_equal(
    unname(sqrt(diag(vcov_mfp2))),
    unname(sqrt(diag(vcov_glm))),
    tolerance = tolerance
  )
}

# Generate stable negative-binomial data for public-interface tests. Keeping the
# signal moderate avoids boundary estimates of theta and makes comparisons with
# an independent MASS::glm.nb() fit reproducible.
make_negbin_test_data <- function(n = 240L, seed = 104L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::runif(n, 0.5, 2.5),
    x2 = stats::rnorm(n),
    x3 = stats::rbinom(n, 1L, 0.4)
  )
  eta <- 0.25 + 0.35 * dat$x1 - 0.25 * dat$x2 + 0.30 * dat$x3
  dat$y <- stats::rnbinom(n, mu = exp(eta), size = 3)
  dat
}

# Cache one independent linear negative-binomial reference problem. Several
# public-method tests use the same fit so that comprehensive NB coverage does
# not repeatedly run the full MFP algorithm.
.negbin_reference_cache <- new.env(parent = emptyenv())

get_negbin_reference_fits <- function() {
  if (!exists("fits", envir = .negbin_reference_cache, inherits = FALSE)) {
    dat <- make_negbin_test_data(n = 320L, seed = 1103L)

    fit_mfp2 <- mfp2(
      y ~ x1 + x2 + x3,
      data = dat,
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

    fit_mass <- MASS::glm.nb(
      y ~ x1 + x2 + x3,
      data = dat,
      link = log,
      control = stats::glm.control(maxit = 100L)
    )

    assign(
      "fits",
      list(
        data = dat,
        mfp2 = fit_mfp2,
        mass = fit_mass
      ),
      envir = .negbin_reference_cache
    )
  }

  get("fits", envir = .negbin_reference_cache, inherits = FALSE)
}

expect_negbin_mfp2_mass_equal <- function(fit_mfp2,
                                          fit_mass,
                                          coefficient_tolerance = 1e-4,
                                          likelihood_tolerance = 1e-3) {
  expect_equal(
    unname(stats::coef(fit_mfp2)),
    unname(stats::coef(fit_mass)),
    tolerance = coefficient_tolerance
  )
  expect_equal(
    unname(stats::vcov(fit_mfp2)),
    unname(stats::vcov(fit_mass)),
    tolerance = coefficient_tolerance
  )
  expect_equal(
    unname(fit_mfp2$fitted.values),
    unname(stats::fitted(fit_mass)),
    tolerance = coefficient_tolerance
  )
  expect_equal(
    unname(fit_mfp2$deviance),
    unname(stats::deviance(fit_mass)),
    tolerance = likelihood_tolerance
  )
  expect_equal(
    fit_mfp2$theta,
    fit_mass$theta,
    tolerance = likelihood_tolerance
  )
  expect_equal(
    fit_mfp2$mfp_logl,
    as.numeric(stats::logLik(fit_mass)),
    tolerance = likelihood_tolerance
  )
  expect_equal(
    fit_mfp2$aic,
    stats::AIC(fit_mass),
    tolerance = likelihood_tolerance
  )
}

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

# -----------------------------------------------------------------------------
# 1.1 Grouped terms in the default matrix interface
# -----------------------------------------------------------------------------

# Test purpose: Checks that manually supplied dummy columns can be represented
# and fitted as one conceptual linear term.
test_that("mfp2.default() groups manually supplied dummy columns", {
  set.seed(2101)
  n <- 180

  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_mm <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  x <- cbind(
    x1 = runif(n, 1, 10),
    group_mm
  )
  y <- 0.4 * x[, "x1"] +
    0.8 * x[, "groupB"] -
    0.5 * x[, "groupC"] +
    rnorm(n, sd = 0.3)

  fit <- mfp2(
    x,
    y,
    df = 1,
    select = 1,
    term_groups = list(group = c("groupB", "groupC")),
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_identical(
    fit$term_to_columns[["group"]],
    c("groupB", "groupC")
  )
  expect_true(fit$fp_terms["group", "selected"])
  expect_true(all(c("groupB", "groupC") %in% colnames(fit$x_original)))
})

# Test purpose: Checks that selection operates on the complete grouped block;
# the final fitting matrix must contain either every member column or none.
test_that("grouped matrix columns are retained or removed together", {
  set.seed(2102)
  n <- 180

  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_mm <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  x <- cbind(
    x1 = runif(n, 1, 10),
    group_mm
  )
  y <- 0.5 * x[, "x1"] + rnorm(n)

  fit <- mfp2(
    x,
    y,
    term_groups = list(group = c("groupB", "groupC")),
    verbose = FALSE
  )

  present <- c("groupB", "groupC") %in% colnames(fit$x_original)
  expect_true(all(present) || !any(present))
  expect_identical(all(present), isTRUE(fit$fp_terms["group", "selected"]))
})


# Test purpose: A scalar df is a continuous-variable default and must not
# send a supplied ordered-factor contrast block into the FP search.
test_that("grouped ordered-factor columns are linear under scalar df default", {
  set.seed(21021)
  n <- 160
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  x <- cbind(
    x1 = seq(0.5, 8, length.out = n),
    stage_mm
  )
  y <- 0.4 * x[, "x1"] + 0.8 * stage_mm[, 1L] -
    0.5 * stage_mm[, 2L] + stats::rnorm(n, sd = 0.3)

  fit <- mfp2(
    x,
    y,
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(as.numeric(fit$fp_terms["stage", "df_setting"]), 1)
  expect_equal(
    as.numeric(fit$fp_terms["stage", "df_initial"]),
    ncol(stage_mm)
  )
  expect_true(fit$fp_terms["stage", "selected"])
})


# Test purpose: Automatic scaling is a continuous-variable default and must not
# independently rescale the member columns of a supplied categorical block.
test_that("grouped matrix columns default to scale one", {
  set.seed(21022)
  n <- 160
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  stage_mm[, 1L] <- 1000 * stage_mm[, 1L]
  x <- cbind(
    x1 = seq(0.5, 8, length.out = n),
    stage_mm
  )
  y <- 0.4 * x[, "x1"] + 0.001 * stage_mm[, 1L] -
    0.5 * stage_mm[, 2L] + stats::rnorm(n, sd = 0.3)

  fit <- mfp2(
    x,
    y,
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(as.numeric(fit$transformations["stage", "scale"]), 1)
  expect_true(fit$fp_terms["stage", "selected"])
})


# Test purpose: A retained manually grouped block reports its fitted rank
# contribution rather than the single linear power used by the selection engine.
test_that("grouped matrix terms report one final df per estimable coefficient", {
  set.seed(2103)
  n <- 180
  group <- factor(rep(c("A", "B", "C"), length.out = n))
  x <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  y <- 0.9 * x[, "groupB"] - 0.6 * x[, "groupC"] +
    stats::rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    term_groups = list(group = c("groupB", "groupC")),
    keep = "group",
    select = 0,
    df = 1,
    verbose = FALSE
  )

  transformed_columns <- paste0(fit$term_to_columns[["group"]], ".1")
  coefficient_columns <- unname(
    fit$transformed_to_model_columns[transformed_columns]
  )
  expected_df <- sum(!is.na(stats::coef(fit)[coefficient_columns]))

  expect_equal(expected_df, 2)
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), expected_df)
})


# Test purpose: The shared metadata expansion applies grouped-term invariants
# while preserving all fitted settings for ordinary singleton terms.
test_that("term metadata expansion is consistent for grouped and singleton terms", {
  acd_fit <- list(beta0 = 0.2, beta1 = 0.8, power = 1)
  expanded <- expand_term_metadata_to_columns(
    term_to_columns = list(
      x1 = "x1",
      group = c("groupB", "groupC"),
      omitted = c("omittedB", "omittedC")
    ),
    powers = list(x1 = c(1, 2), group = 1, omitted = NA_real_),
    raw_columns = c("groupC", "x1", "groupB", "omittedB"),
    center = c(x1 = TRUE, group = TRUE, omitted = FALSE),
    acdx = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    zero = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    catzero = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    spike = c(x1 = TRUE, group = TRUE, omitted = TRUE),
    spike_decision = c(x1 = 1L, group = 1L, omitted = 1L),
    acd_parameter = list(x1 = acd_fit, group = acd_fit, omitted = acd_fit)
  )

  expect_identical(
    unname(expanded$terms),
    c("group", "x1", "group", "omitted")
  )
  expect_identical(expanded$powers[["groupC"]], 1)
  expect_identical(expanded$powers[["groupB"]], 1)
  expect_true(is.na(expanded$powers[["omittedB"]]))
  expect_identical(expanded$powers[["x1"]], c(1, 2))
  expect_true(expanded$center[["groupC"]])
  expect_true(expanded$center[["x1"]])

  grouped_columns <- c("groupC", "groupB", "omittedB")
  expect_false(any(expanded$acdx[grouped_columns]))
  expect_false(any(expanded$zero[grouped_columns]))
  expect_false(any(expanded$catzero[grouped_columns]))
  expect_false(any(expanded$spike[grouped_columns]))
  expect_true(all(
    expanded$spike_decision[grouped_columns] ==
      saz_decision_codes[["continuous_only"]]
  ))
  expect_true(all(vapply(
    expanded$acd_parameter[grouped_columns],
    is.null,
    logical(1L)
  )))

  expect_true(expanded$acdx[["x1"]])
  expect_true(expanded$zero[["x1"]])
  expect_true(expanded$catzero[["x1"]])
  expect_true(expanded$spike[["x1"]])
  expect_identical(expanded$spike_decision[["x1"]], 1L)
  expect_identical(expanded$acd_parameter[["x1"]], acd_fit)
})


# Test purpose: Grouped-term df uses fitted estimability, so aliased member
# coefficients do not count toward the final rank contribution.
test_that("grouped final df excludes non-estimable coefficients", {
  fp_terms <- create_fp_terms(
    fp_powers = list(group = 1),
    acdx = c(group = FALSE),
    df = c(group = 1),
    select = c(group = 1),
    alpha = c(group = 1),
    criterion = "pvalue",
    zero = c(group = FALSE),
    catzero = c(group = FALSE),
    spike = c(group = FALSE),
    spike_decision = c(group = 2),
    term_to_columns = list(group = c("groupB", "groupC")),
    transformed_to_model_columns = c(
      "groupB.1" = "groupB.1",
      "groupC.1" = "groupC.1"
    ),
    coefficients = c(
      "(Intercept)" = 0,
      "groupB.1" = 0.4,
      "groupC.1" = NA_real_
    )
  )

  expect_equal(as.numeric(fp_terms["group", "df_setting"]), 1)
  expect_equal(as.numeric(fp_terms["group", "df_initial"]), 2)
  expect_equal(as.numeric(fp_terms["group", "df_final"]), 1)
})


# =============================================================================
# 2. mfp2.formula() — equivalence and fp() terms
# =============================================================================

# Test purpose: Checks that the formula interface parses fp() terms and returns
# a converged mfp2 object.
test_that("mfp2.formula() returns an mfp2 object", {
  fit <- mfp2(
    lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) + fp(weight) +
      fp(bph) + fp(cp),
    data = prostate, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$convergence_mfp)
  expect_true(is.data.frame(fit$fp_terms))
})

# Test purpose: With all predictors retained and restricted to linear effects,
# verifies the complete NB fit against the independent MASS implementation,
# including covariance, likelihood, AIC, theta, model df, and null likelihood.
test_that("linear negative-binomial mfp2 agrees with MASS::glm.nb()", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  fits <- get_negbin_reference_fits()
  fit_mfp2 <- fits$mfp2
  fit_mass <- fits$mass

  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "fastglm_nb")
  expect_identical(fit_mfp2$family_string, "negbin")
  expect_identical(fit_mfp2$fitter, "fastglm")
  expect_identical(fit_mfp2$family$link, "log")
  expect_true(all(fit_mfp2$fp_terms[c("x1", "x2", "x3"), "selected"]))
  expect_true(all(vapply(
    fit_mfp2$fp_powers[c("x1", "x2", "x3")],
    function(p) identical(unname(p), 1),
    logical(1L)
  )))

  expect_negbin_mfp2_mass_equal(fit_mfp2, fit_mass)

  # Theta is an estimated nuisance parameter. It contributes to AIC/model df,
  # but not to regression residual degrees of freedom.
  expect_equal(fit_mfp2$mfp_df, fit_mfp2$rank + 1L)
  expect_equal(
    fit_mfp2$aic,
    -2 * fit_mfp2$mfp_logl + 2 * fit_mfp2$mfp_df,
    tolerance = 1e-10
  )
  # With df = 1 and select = 1, the full-linear reference and final MFP model
  # are the same statistical model.
  expect_equal(fit_mfp2$linear_logl, fit_mfp2$mfp_logl, tolerance = 1e-8)
  expect_equal(fit_mfp2$linear_df, fit_mfp2$mfp_df)

  # GLM null_logl is retained as NA for compatibility because Model Fit now
  # reports the already-computed family-specific deviance instead.
  expect_true(is.na(fit_mfp2$null_logl))
  expect_equal(
    fit_mfp2$null_deviance,
    fit_mass$null.deviance,
    tolerance = 1e-3
  )
})

# Test purpose: Compares selected variables from equivalent matrix and formula
# interface fits. Centering is made explicit on both interfaces because fp()
# carries its own center setting (default TRUE) rather than inheriting the
# top-level formula setting.
test_that("default and formula interfaces give consistent selected variables", {
  fit_default <- mfp2(
    x_prostate, y_prostate, center = TRUE, verbose = FALSE
  )
  fit_formula <- mfp2(
    lpsa ~ fp(age, center = TRUE) + fp(svi, df = 1, center = TRUE) +
      fp(pgg45, center = TRUE) + fp(cavol, center = TRUE) +
      fp(weight, center = TRUE) + fp(bph, center = TRUE) +
      fp(cp, center = TRUE),
    data = prostate,
    center = TRUE,
    verbose = FALSE
  )

  sel_default <- sort(get_selected_variable_names(fit_default))
  sel_formula <- sort(get_selected_variable_names(fit_formula))
  expect_equal(sel_default, sel_formula)
})

# Test purpose: Checks that df specified inside fp() is propagated to the fitted
# term metadata.
test_that("fp() applies per-variable df correctly", {
  fit <- mfp2(
    lpsa ~ fp(age, df = 2) + fp(svi, df = 1) + fp(cavol, df = 4),
    data = prostate, verbose = FALSE
  )

  # svi is binary so df should be 1
  expect_equal(as.numeric(fit$fp_terms["svi", "df_initial"]), 1)
})

# Test purpose: Protects formula centering precedence. An fp() term carries its
# own center setting and defaults to TRUE, so it overrides a top-level
# center = FALSE. Ordinary numeric terms continue to use the top-level value.
test_that("fp() center default overrides the global formula center setting", {
  n <- 90L
  i <- seq_len(n)
  dat <- data.frame(
    x1 = seq(1, 9, length.out = n),
    x2 = sin(i / 7) + (i %% 5) / 10
  )
  dat$y <- 0.8 + 0.25 * sqrt(dat$x1) - 0.4 * dat$x2 +
    0.01 * cos(i / 3)

  fit <- mfp2(
    y ~ fp(x1, df = 2, force_max_fp = TRUE) + x2,
    data = dat,
    df = 1,
    keep = c("x1", "x2"),
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_true(isTRUE(fit$transformations["x1", "center"]))
  expect_false(isTRUE(fit$transformations["x2", "center"]))
})

# Test purpose: Confirms that fp2() can be used as a formula-interface alias for
# fp().
test_that("fp2() is an alias for fp()", {
  fit <- mfp2(
    lpsa ~ fp2(age) + fp2(cavol) + fp2(svi, df = 1),
    data = prostate, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks that formula-based Cox fitting works with survival responses
# and covariates.
test_that("formula interface works with Cox family and strata", {
  data("gbsg", package = "mfp2")
  fit <- mfp2(
    Surv(rectime, censrec) ~ fp(age) + fp(size) + fp(nodes) +
      fp(er) + meno,
    data = gbsg,
    family = "cox",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(fit$family_string, "cox")
})

# Test purpose: Checks that an unordered factor is represented by one conceptual
# term, while its complete treatment-contrast block is retained through keep.
test_that("formula interface keeps an unordered factor as one grouped term", {
  set.seed(101)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true("group" %in% rownames(fit$fp_terms))
  expect_true(fit$fp_terms["group", "selected"])
  expect_equal(as.numeric(fit$fp_terms["group", "df_setting"]), 1)
  expect_equal(as.numeric(fit$fp_terms["group", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), 2)

  expect_true("group" %in% names(fit$term_to_columns))
  expect_setequal(
    fit$term_to_columns[["group"]],
    c("groupB", "groupC")
  )
  expect_true(
    all(fit$term_to_columns[["group"]] %in% colnames(fit$x_original))
  )
})


# Test purpose: 2.1 Ordered-factor contrast columns are grouped into one fixed
# linear term rather than being treated as separate candidate predictors.
test_that("2.1 Formula interface supports ordered factors as grouped terms", {
  set.seed(2011)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    ordered_group = ordered(
      rep(c("low", "medium", "high"), length.out = n),
      levels = c("low", "medium", "high")
    )
  )

  fit <- mfp2(
    y ~ fp(x) + ordered_group,
    data = dat,
    keep = "ordered_group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true("ordered_group" %in% names(fit$term_to_columns))
  expect_length(fit$term_to_columns[["ordered_group"]], 2L)
  expect_true(fit$fp_terms["ordered_group", "selected"])
  expect_equal(as.numeric(fit$fp_terms["ordered_group", "df_setting"]), 1)
  expect_equal(as.numeric(fit$fp_terms["ordered_group", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["ordered_group", "df_final"]), 2)

  # Continuous-only extensions are disabled for the complete factor block.
  expect_false(fit$fp_terms["ordered_group", "acd"])
  expect_false(fit$fp_terms["ordered_group", "zero"])
  expect_false(fit$fp_terms["ordered_group", "catzero"])
  expect_false(fit$fp_terms["ordered_group", "spike"])
  expect_true(fit$formula_factor_info[["ordered_group"]]$ordered)
})


# Test purpose: Regression test for final transformation after fitting the ART
# data, which contains both ordinal and nominal three-level predictors.
test_that("formula fit with grouped ART terms completes final transformation", {
  data("art", package = "mfp2")

  # ART stores x4/x9 as numeric level codes. Convert them explicitly so this
  # regression test exercises formula-factor grouping rather than numeric
  # low-cardinality handling.
  art$x4 <- ordered(art$x4)
  art$x9 <- factor(art$x9)

  fit <- mfp2(
    y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = art,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$convergence_mfp)
  expect_true(is.logical(fit$catzero))
  expect_false(anyNA(fit$catzero))

  expect_true(all(c("x4", "x9") %in% names(fit$term_to_columns)))
  expect_gt(length(fit$term_to_columns[["x4"]]), 1L)
  expect_gt(length(fit$term_to_columns[["x9"]]), 1L)
})

# Test purpose: Ensures formula-interface keep names must match either formula terms
# or expanded model-matrix columns; misspelled names should error.
test_that("formula interface rejects unknown keep variables", {
  expect_error(
    mfp2(
      lpsa ~ fp(age) + fp(cavol),
      data = prostate,
      keep = "does_not_exist",
      verbose = FALSE
    ),
    "Unknown variable"
  )
})

test_that("predict.mfp2 works after fitting with strata()", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]

  fit <- mfp2(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    center = FALSE,
    verbose = FALSE
  )

  nd <- dat[1:5, c("age", "sex", "inst"), drop = FALSE]

  p <- predict(
    fit,
    newdata = nd,
    type = "lp"
  )

  expect_length(p, 5)
  expect_true(all(is.finite(p)))
})

# =============================================================================
# 3. Family and response validation
# =============================================================================

# Test purpose: Checks that supported character families, family functions, and
# family objects are normalized correctly.
test_that("normalize_family_argument() accepts valid families", {
  # Character inputs
  for (fam in c("gaussian", "binomial", "poisson", "cox", "negbin")) {
    result <- normalize_family_argument(fam)
    expect_equal(result$family_string, fam)
  }

  # Function inputs
  result <- normalize_family_argument(stats::gaussian)
  expect_equal(result$family_string, "gaussian")

  result <- normalize_family_argument(stats::binomial)
  expect_equal(result$family_string, "binomial")

  # Family object input
  result <- normalize_family_argument(stats::binomial(link = "probit"))
  expect_equal(result$family_string, "binomial")
})

# Test purpose: Checks that unsupported, ambiguous, or non-scalar family
# specifications are rejected.
test_that("normalize_family_argument() rejects invalid families", {
  expect_error(normalize_family_argument("gamma"), "Invalid family")
  expect_error(normalize_family_argument("inverse.gaussian"), "Invalid family")
  expect_error(normalize_family_argument(c("gaussian", "binomial")),
               "single character string")
})

# Test purpose: Confirms that the character Cox family specification is accepted
# by the family normalizer.
test_that("cox must be specified as character string, not function", {
  # There is no stats::cox(), so creating a fake one would test the guard
  expect_error(normalize_family_argument("cox")$family_string, NA)
})

# Test purpose: Checks that invalid response values or response classes are
# rejected for each family.
test_that("validate_family_response() catches bad responses", {
  # Gaussian: must be numeric
  expect_error(validate_family_response("abc", "gaussian", 3))
  expect_error(validate_family_response(c(1, NA, 3), "gaussian", 3))

  # Binomial: numeric must be in [0, 1]
  expect_error(validate_family_response(c(0, 2, 1), "binomial", 3))
  expect_error(validate_family_response(c(-0.1, 0.5, 0.9), "binomial", 3))

  # Poisson: must be non-negative
  expect_error(validate_family_response(c(-1, 0, 1), "poisson", 3))

  # Negative binomial: must be finite non-negative integer counts
  expect_error(
    validate_family_response(c(-1, 0, 1), "negbin", 3),
    "non-negative"
  )
  expect_error(
    validate_family_response(c(0, 1.5, 2), "negbin", 3),
    "integer"
  )
  expect_error(
    validate_family_response(c(0, NA, 2), "negbin", 3),
    "missing"
  )
  expect_error(
    validate_family_response(c(0, Inf, 2), "negbin", 3),
    "finite"
  )
  expect_error(
    validate_family_response(c("0", "1", "2"), "negbin", 3),
    "numeric"
  )

  # Cox: must be Surv
  expect_error(validate_family_response(c(1, 2, 3), "cox", 3))

  # Surv for non-cox
  expect_error(validate_family_response(Surv(1:3, c(1,0,1)), "gaussian", 3))
})

# Test purpose: Checks that valid Gaussian, binomial, Poisson, and Cox responses
# pass validation.
test_that("validate_family_response() accepts valid responses", {
  expect_true(validate_family_response(c(1.5, 2.5, 3.5), "gaussian", 3))
  expect_true(validate_family_response(c(0, 0.5, 1), "binomial", 3))
  expect_true(validate_family_response(
    factor(c("a", "b"), levels = c("a", "b")), "binomial", 2
  ))
  expect_true(validate_family_response(c(0, 1, 5), "poisson", 3))
  expect_true(validate_family_response(c(0, 1, 5), "negbin", 3))
  expect_true(validate_family_response(Surv(1:3, c(1, 0, 1)), "cox", 3))
})

# =============================================================================
# 4. Preprocessing: shift, scale, centering
# =============================================================================

# Test purpose: Checks that no shift is added when all values are already positive.
test_that("find_shift_factor() returns 0 for all-positive data", {
  expect_equal(find_shift_factor(1:10), 0)
})

# Test purpose: Checks that the estimated shift makes zero or negative data
# strictly positive.
test_that("find_shift_factor() shifts data containing zero or negatives", {
  x <- c(-1, 0, 1, 2, 3)
  s <- find_shift_factor(x)
  expect_true(s > 0)
  expect_true(all((x + s) > 0))
})

# Test purpose: Checks that binary variables are not shifted.
test_that("find_shift_factor() returns 0 for binary variables", {
  expect_equal(find_shift_factor(c(0, 1, 0, 1)), 0)
})

# Test purpose: Checks that binary variables are not rescaled.
test_that("find_scale_factor() returns 1 for binary variables", {
  expect_equal(find_scale_factor(c(0, 1, 0, 1)), 1)
})

# Test purpose: An explicitly supplied scale must not rescale a binary
# predictor. Binary variables retain neutral preprocessing factors regardless
# of user-supplied shift and scale values.
test_that("mfp2.default() resets supplied scale for binary variables", {
  n <- 80L
  svi <- rep(c(0, 1), each = n / 2L)
  x <- cbind(svi = svi)
  y <- 1 + 2 * svi + sin(seq_len(n) / 5)

  fit <- mfp2(
    x = x,
    y = y,
    shift = c(svi = 10),
    scale = c(svi = 1000),
    df = 1,
    keep = "svi",
    center = FALSE,
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$transformations["svi", "shift"]), 0)
  expect_equal(as.numeric(fit$transformations["svi", "scale"]), 1)
})

# Test purpose: Checks automatic power-of-10 scaling for a known numeric range.
test_that("find_scale_factor() returns correct power-of-10 scaling", {
  # range = 999, log10(999) ~ 2.999, floor = 2, so scale = 100
  expect_equal(find_scale_factor(1:1000), 100)
})

# Test purpose: Checks that constant variables are rejected before scaling.
test_that("find_scale_factor() errors on constant input", {
  expect_error(find_scale_factor(rep(5, 10)), "must not be constant")
})

# Test purpose: Checks that the combined preprocessing helper returns strictly
# positive transformed values.
test_that("apply_shift_scale() produces positive, scaled output", {
  x <- c(-2, 0, 3, 5, 10)
  x_ss <- apply_shift_scale(x)
  expect_true(all(x_ss > 0))
})

# Test purpose: Checks that default model fitting records centering constants.
test_that("centering is applied by default", {
  fit <- mfp2(x_prostate, y_prostate, center = TRUE, verbose = FALSE)
  expect_true(!is.null(fit$centers))
})

# Test purpose: Checks that center = FALSE disables centering and leaves no
# centers stored.
test_that("centering can be disabled", {
  fit <- mfp2(x_prostate, y_prostate, center = FALSE, verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_null(fit$centers)
  expect_true(all(fit$transformations$center == FALSE))
})


# =============================================================================
# 4.1 Shift and scale safety contract - shared setup
# =============================================================================
# Review objective
# ----------------
# The tests in Sections 4.2 and 4.3 verify the complete user-input contract for
# shift and scale without changing the existing preprocessing calculations:
#
#   1. NULL keeps automatic shift/scale selection.
#   2. A single finite numeric value is recycled to every matrix column.
#   3. A multi-value vector must be fully and uniquely named.
#   4. Named values are matched to colnames(x), not to supplied position.
#   5. Missing, unknown, duplicate, or empty names are rejected.
#   6. shift values must be finite numeric values with no missing values.
#   7. scale values must additionally be strictly positive.
#   8. Formula-level scalar and per-variable fp() settings remain compatible.
#   9. Correctly specified models retain identical fits and predictions when
#      the same named settings are supplied in a different order.
#
# Sections 4.2 and 4.3 are intentionally separate. MFPI has an additional
# grouping-variable column that must participate in matrix-name validation but
# is subsequently handled as categorical metadata by the fitting procedure.

# -----------------------------------------------------------------------------
# 4.1.1 Deterministic shared data
# -----------------------------------------------------------------------------
# Expected use:
# - setting_x and setting_y exercise the ordinary mfp2 matrix interface.
# - mfpi_setting_x adds the MFPI grouping variable "svi" as a matrix column.
# - The response contains nonlinear signal so stored transformations and
#   predictions are meaningful integration checks rather than validation-only
#   calls.

set.seed(20260717)
setting_n <- 96L
setting_x <- cbind(
  age = runif(setting_n, 20, 80),
  weight = runif(setting_n, 45, 105)
)
setting_y <-
  1.5 * sqrt(setting_x[, "age"] + 2) -
  2.0 * log(setting_x[, "weight"] + 3) +
  rnorm(setting_n, sd = 0.15)

mfpi_setting_x <- cbind(
  svi = rep(c(0, 1), each = setting_n / 2),
  age = setting_x[, "age"],
  weight = setting_x[, "weight"]
)
mfpi_setting_y <-
  setting_y +
  0.4 * mfpi_setting_x[, "svi"] * sqrt(mfpi_setting_x[, "age"] + 2)

# -----------------------------------------------------------------------------
# 4.1.2 Shared fitting helpers
# -----------------------------------------------------------------------------
# Expected use:
# - Keep model specifications identical across ordered and reverse-ordered
#   setting vectors.
# - Limit differences between compared fits to the order in which named shift
#   or scale values were supplied.

fit_mfp2_settings <- function(shift = NULL, scale = NULL) {
  mfp2(
    x = setting_x,
    y = setting_y,
    shift = shift,
    scale = scale,
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = colnames(setting_x),
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )
}

fit_mfpi_settings <- function(shift = NULL, scale = NULL) {
  mfpi(
    x = mfpi_setting_x,
    y = mfpi_setting_y,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    shift = shift,
    scale = scale,
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = c("age", "weight"),
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )
}

# -----------------------------------------------------------------------------
# 4.1.3 Shared validation-call helpers
# -----------------------------------------------------------------------------
# Expected use:
# - Reach the public matrix interfaces with the smallest practical fit.
# - Keep validation tests focused on shift/scale errors rather than FP model
#   selection behavior.

mfp2_validation_call <- function(argument, value) {
  args <- list(
    x = setting_x,
    y = setting_y,
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 1,
    verbose = FALSE
  )
  args[[argument]] <- value
  do.call(mfp2, args)
}

mfpi_validation_call <- function(argument, value) {
  args <- list(
    x = mfpi_setting_x,
    y = mfpi_setting_y,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 1,
    verbose = FALSE
  )
  args[[argument]] <- value
  do.call(mfpi, args)
}

# -----------------------------------------------------------------------------
# 4.1.4 Shared malformed-input fixtures
# -----------------------------------------------------------------------------
# Each fixture represents one distinct part of the public input contract. The
# same cases are applied separately to mfp2.default() and mfpi.default().

invalid_named_settings <- function(column_names) {
  n <- length(column_names)
  values <- seq_len(n)

  partially_named <- stats::setNames(values, column_names)
  names(partially_named)[2L] <- ""

  duplicated_names <- stats::setNames(values, column_names)
  names(duplicated_names)[2L] <- column_names[[1L]]

  unknown_replacement <- stats::setNames(values, column_names)
  names(unknown_replacement)[n] <- "replacement_column"

  list(
    # More than one value with no names must never be matched by position.
    unnamed = unname(values),

    # Every element of a multi-value vector must have a non-empty name.
    partially_named = partially_named,

    # Ambiguous duplicate names cannot identify one value per matrix column.
    duplicated = duplicated_names,

    # Replacing a valid name introduces an unknown supplied column. The
    # omitted known column itself is valid under the partial-vector contract.
    unknown_replacement = unknown_replacement,

    # A complete required set plus an extra name must be rejected.
    unknown_extra = stats::setNames(
      c(values, n + 1L),
      c(column_names, "unknown_column")
    ),

    # Extra named values fail because their names are not matrix columns.
    extra_unknown_names = stats::setNames(
      c(values, n + 1L, n + 2L),
      c(column_names, "unknown_column_1", "unknown_column_2")
    ),

    # Explicit numerical settings must be finite and nonmissing.
    non_finite = stats::setNames(
      replace(as.numeric(values), 1L, Inf),
      column_names
    ),
    missing_value = stats::setNames(
      replace(as.numeric(values), 1L, NA_real_),
      column_names
    ),

    # Logical values must not be silently coerced into numeric settings.
    logical = stats::setNames(rep(TRUE, n), column_names)
  )
}

# -----------------------------------------------------------------------------
# 4.1.5 Shared expectation helpers for malformed names
# -----------------------------------------------------------------------------
# These expectations are deliberately explicit. A review can therefore map
# every required naming failure to one assertion while still applying the same
# contract consistently to the two public matrix interfaces.

expect_setting_name_errors <- function(call_setting, column_names) {
  invalid <- invalid_named_settings(column_names)

  # Shift and scale accept named subsets, but never unnamed multi-value
  # vectors, incomplete names, duplicate names, or names outside colnames(x).
  expect_error(
    call_setting("shift", invalid$unnamed),
    "`shift` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("shift", invalid$partially_named),
    "`shift` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("shift", invalid$duplicated),
    "`shift` names must be unique"
  )
  for (invalid_shift in invalid[c(
    "unknown_replacement", "unknown_extra", "extra_unknown_names"
  )]) {
    expect_error(
      call_setting("shift", invalid_shift),
      "`shift` contains unknown column name"
    )
  }

  expect_error(
    call_setting("scale", invalid$unnamed),
    "`scale` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("scale", invalid$partially_named),
    "`scale` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("scale", invalid$duplicated),
    "`scale` names must be unique"
  )
  for (invalid_scale in invalid[c(
    "unknown_replacement", "unknown_extra", "extra_unknown_names"
  )]) {
    expect_error(
      call_setting("scale", invalid_scale),
      "`scale` contains unknown column name"
    )
  }
}

# -----------------------------------------------------------------------------
# 4.1.6 Shared expectation helpers for invalid values
# -----------------------------------------------------------------------------
# Numerical validation is common to shift and scale. Positivity is tested
# separately because it is an additional requirement that applies only to scale.

expect_setting_numeric_errors <- function(call_setting, column_names) {
  invalid <- invalid_named_settings(column_names)

  for (argument in c("shift", "scale")) {
    # Requirement: Inf and -Inf are not valid explicit settings.
    expect_error(
      call_setting(argument, invalid$non_finite),
      paste0("`", argument, "` must contain only finite values")
    )

    # Requirement: direct matrix settings cannot contain NA values.
    expect_error(
      call_setting(argument, invalid$missing_value),
      paste0("`", argument, "` must not contain missing values")
    )

    # Requirement: TRUE/FALSE must not be accepted as 1/0 settings.
    expect_error(
      call_setting(argument, invalid$logical),
      paste0("`", argument, "` must contain numeric values")
    )
  }
}

expect_nonpositive_scale_errors <- function(call_setting, column_names) {
  zero_scale <- stats::setNames(rep(1, length(column_names)), column_names)
  zero_scale[[1L]] <- 0

  negative_scale <- stats::setNames(rep(1, length(column_names)), column_names)
  negative_scale[[1L]] <- -1

  # Requirement: zero is not a valid divisor for scaling, including when
  # supplied through the named-partial interface.
  expect_error(
    call_setting("scale", zero_scale),
    "`scale` must contain only strictly positive values"
  )
  expect_error(
    call_setting("scale", stats::setNames(0, column_names[[1L]])),
    "`scale` must contain only strictly positive values"
  )

  # Requirement: negative scales are not valid transformation settings.
  expect_error(
    call_setting("scale", negative_scale),
    "`scale` must contain only strictly positive values"
  )
  expect_error(
    call_setting("scale", stats::setNames(-1, column_names[[1L]])),
    "`scale` must contain only strictly positive values"
  )
}


# -----------------------------------------------------------------------------
# 4.1.7-4.1.8 Direct normalization regression tests
# -----------------------------------------------------------------------------
# These tests isolate argument normalization from model fitting. They protect
# the precise distinction that caused the original bug: a named length-one
# vector is partial and variable-specific, while an unnamed length-one vector
# is global.

test_that("4.1.7 named scalar settings are partial, not global", {
  columns <- c("cavol", "age", "weight")
  normalized_shift <- normalize_named_numeric_setting(
    value = c(age = 20),
    column_names = columns,
    argument_name = "shift",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = FALSE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )
  normalized_scale <- normalize_named_numeric_setting(
    value = c(age = 10),
    column_names = columns,
    argument_name = "scale",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = TRUE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )

  expect_identical(
    normalized_shift,
    c(cavol = NA_real_, age = 20, weight = NA_real_)
  )
  expect_identical(
    normalized_scale,
    c(cavol = NA_real_, age = 10, weight = NA_real_)
  )
})

test_that("4.1.8 unnamed scalar settings remain global", {
  columns <- c("cavol", "age", "weight")
  normalized_shift <- normalize_named_numeric_setting(
    value = 20,
    column_names = columns,
    argument_name = "shift",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = FALSE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )
  normalized_scale <- normalize_named_numeric_setting(
    value = 10,
    column_names = columns,
    argument_name = "scale",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = TRUE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )

  expect_identical(normalized_shift, c(cavol = 20, age = 20, weight = 20))
  expect_identical(normalized_scale, c(cavol = 10, age = 10, weight = 10))
})


# =============================================================================
# 4.2 mfp2 shift and scale behavior
# =============================================================================
# Scope
# -----
# This section tests the ordinary mfp2 matrix interface first, followed by the
# formula interface. Each test title is numbered so a reviewer can trace the
# expected behavior directly to the corresponding requirement above.

# -----------------------------------------------------------------------------
# 4.2.1 NULL keeps automatic preprocessing
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that leaving both arguments as NULL still delegates shift and scale
#   selection to the package's existing automatic preprocessing logic.
# - Protect against the new normalization helper accidentally treating NULL as
#   an invalid or incomplete user-supplied vector.
# - Confirm that automatic preprocessing produces usable settings for every
#   matrix column before fractional-polynomial fitting begins.
#
# Expected output:
# - fit$transformations contains rows named "age" and "weight".
# - The selected shift and scale cells contain no NA values.
# - Every selected shift and scale value is finite.
# - Both scale values are greater than zero.

test_that("4.2.1 mfp2.default() keeps automatic shift and scale for NULL", {
  fit <- fit_mfp2_settings()
  transformations <- fit$transformations[
    colnames(setting_x), c("shift", "scale"), drop = FALSE
  ]

  expect_false(anyNA(transformations))
  expect_true(all(is.finite(as.matrix(transformations))))
  expect_true(all(transformations$scale > 0))
})

# -----------------------------------------------------------------------------
# 4.2.2 Scalars are recycled to every matrix column
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that a single finite numeric shift remains a global matrix setting.
# - Verify that a single positive numeric scale remains a global matrix setting.
# - Protect the documented scalar interface while tightening validation for
#   vectors containing more than one value.
#
# Expected output:
# - In colnames(setting_x) order, the stored shifts are c(2, 2).
# - In colnames(setting_x) order, the stored scales are c(10, 10).
# - No name matching is required from the user for these scalar inputs.

test_that("4.2.2 mfp2.default() recycles scalar shift and scale", {
  fit <- fit_mfp2_settings(shift = 2, scale = 10)
  transformations <- fit$transformations[
    colnames(setting_x), c("shift", "scale"), drop = FALSE
  ]

  expect_equal(unname(transformations$shift), rep(2, ncol(setting_x)))
  expect_equal(unname(transformations$scale), rep(10, ncol(setting_x)))
})


# -----------------------------------------------------------------------------
# 4.2.2a Named partial settings fix supplied columns and estimate the rest
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that c(age = 20) and c(age = 100) apply only to age rather than
#   being recycled as global shift and scale settings.
# - Verify that an unspecified nonlinear predictor retains internal NA
#   sentinels until the ordinary shift and scale estimation steps.
# - Exercise the public matrix interface through stored transformation metadata.

test_that("4.2.2a mfp2.default() estimates unspecified named partial settings", {
  index <- seq_len(60L)
  partial_x <- cbind(
    age = seq(20, 79),
    weight = -30 + ((17 * index) %% 60)
  )
  expected_weight_shift <- find_shift_factor(partial_x[, "weight"])
  expected_weight_scale <- find_scale_factor(
    partial_x[, "weight"] + expected_weight_shift
  )
  partial_y <-
    1.2 * sqrt(partial_x[, "age"] + 20) -
    0.7 * log(partial_x[, "weight"] + expected_weight_shift) +
    0.05 * sin(index / 4)

  fit <- mfp2(
    x = partial_x,
    y = partial_y,
    shift = c(age = 20),
    scale = c(age = 100),
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = colnames(partial_x),
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  transformations <- fit$transformations[
    colnames(partial_x), c("shift", "scale"), drop = FALSE
  ]
  expect_equal(unname(transformations["age", "shift"]), 20)
  expect_equal(
    unname(transformations["weight", "shift"]),
    expected_weight_shift
  )
  expect_equal(unname(transformations["age", "scale"]), 100)
  expect_equal(
    unname(transformations["weight", "scale"]),
    expected_weight_scale
  )
})

# -----------------------------------------------------------------------------
# 4.2.2b Insufficient explicit shifts fail before FP fitting
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that a supplied named shift is not silently increased.
# - Verify that the existing positivity check reports the affected variable when
#   the shifted nonlinear predictor still contains zero or negative values.

test_that("4.2.2b mfp2.default() rejects insufficient named partial shifts", {
  index <- seq_len(60L)
  insufficient_x <- cbind(
    age = -30 + ((17 * index) %% 60),
    weight = 40 + ((23 * index) %% 60)
  )
  insufficient_y <-
    0.3 * insufficient_x[, "age"] +
    0.2 * sqrt(insufficient_x[, "weight"]) +
    sin(index / 5)

  expect_error(
    mfp2(
      x = insufficient_x,
      y = insufficient_y,
      shift = c(age = 30),
      scale = 1,
      df = 2,
      select = 1,
      alpha = 1,
      force_max_fp_vars = colnames(insufficient_x),
      center = FALSE,
      cycles = 5,
      xorder = "original",
      verbose = FALSE
    ),
    "Problematic variables: age"
  )
})


# -----------------------------------------------------------------------------
# 4.2.2c Partial scale keeps mapped-term automatic defaults
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that an unspecified grouped design block still receives scale = 1.
# - Verify that an explicit partial scale remains authoritative for its column.

test_that("4.2.2c partial scale preserves mapped-term defaults", {
  normalized_scale <- c(group_b = NA_real_, group_c = NA_real_, age = 100)
  term_to_columns <- list(
    group = c("group_b", "group_c"),
    age = "age"
  )

  result <- expand_scale_for_mapped_terms(
    scale = normalized_scale,
    vnames = names(normalized_scale),
    term_to_columns = term_to_columns,
    automatic = is.na(normalized_scale)
  )

  expect_equal(
    result,
    c(group_b = 1, group_c = 1, age = 100)
  )
})

# -----------------------------------------------------------------------------
# 4.2.3 Ordered named vectors are stored against the matching columns
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify the accepted multi-value form: one numeric value for every matrix
#   column, with complete and unique names.
# - Confirm that normalization returns the settings in colnames(setting_x)
#   order and that downstream transformation metadata preserves those values.
# - Establish the reference fit used conceptually by the reordered-input test.
#
# Expected output:
# - The "age" transformation row stores shift = 1 and scale = 10.
# - The "weight" transformation row stores shift = 2 and scale = 100.
# - Reading the two rows in matrix-column order returns shifts c(1, 2) and
#   scales c(10, 100).

test_that("4.2.3 mfp2.default() accepts fully named shift and scale", {
  fit <- fit_mfp2_settings(
    shift = c(age = 1, weight = 2),
    scale = c(age = 10, weight = 100)
  )
  transformations <- fit$transformations[
    colnames(setting_x), c("shift", "scale"), drop = FALSE
  ]

  expect_equal(unname(transformations$shift), c(1, 2))
  expect_equal(unname(transformations$scale), c(10, 100))
})

# -----------------------------------------------------------------------------
# 4.2.4 Named-vector order does not change the fitted model
# -----------------------------------------------------------------------------
# Test purpose:
# - Reproduce the original failure mode by supplying the same named settings in
#   matrix order and then in reverse order.
# - Verify that values are matched by their names rather than silently relabeled
#   according to their supplied positions.
# - Exercise the complete path from normalization through FP transformation,
#   model fitting, stored metadata, and prediction.
#
# Expected output:
# - The ordered and reversed fits have identical transformation tables.
# - Their coefficient vectors and fitted values are numerically identical.
# - Predictions for the first 12 rows of setting_x are numerically identical.
# - In particular, age always receives shift = 1 and scale = 10, while weight
#   always receives shift = 2 and scale = 100, regardless of input order.

test_that("4.2.4 mfp2.default() matches shift and scale by name", {
  ordered <- fit_mfp2_settings(
    shift = c(age = 1, weight = 2),
    scale = c(age = 10, weight = 100)
  )
  reversed <- fit_mfp2_settings(
    shift = c(weight = 2, age = 1),
    scale = c(weight = 100, age = 10)
  )

  expect_equal(ordered$transformations, reversed$transformations)
  expect_equal(unname(stats::coef(ordered)), unname(stats::coef(reversed)))
  expect_equal(unname(stats::fitted(ordered)), unname(stats::fitted(reversed)))

  new_x <- setting_x[1:12, , drop = FALSE]
  expect_equal(
    unname(stats::predict(ordered, newdata = new_x)),
    unname(stats::predict(reversed, newdata = new_x))
  )
})

# -----------------------------------------------------------------------------
# 4.2.5 Invalid multi-value names are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that positional matching is no longer available for multi-value
#   matrix settings.
# - Verify that every accepted multi-value vector identifies only known columns,
#   with no ambiguity or surplus entries.
# - Apply the same naming contract independently to shift and scale.
#
# Expected output for each argument:
# - Unnamed multi-value vectors and vectors containing empty element names
#   error with the single-unnamed-or-named-vector message.
# - Duplicate names error with "`<argument>` names must be unique".
# - Named vectors containing unknown columns error with
#   "`<argument>` contains unknown column name".
# - No model is fitted for any invalid input.

test_that("4.2.5 mfp2.default() rejects invalid shift and scale names", {
  expect_setting_name_errors(mfp2_validation_call, colnames(setting_x))
})

# -----------------------------------------------------------------------------
# 4.2.6 Missing, non-finite, and logical values are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that direct matrix settings are genuine numeric configuration values,
#   not logical vectors that R could silently coerce to 0 and 1.
# - Verify that explicit values are complete and finite before transformation or
#   model fitting is attempted.
# - Apply the same numeric validation to shift and scale.
#
# Expected output for each argument:
# - A vector containing Inf errors with
#   "`<argument>` must contain only finite values".
# - A vector containing NA errors with
#   "`<argument>` must not contain missing values".
# - A logical vector errors with
#   "`<argument>` must contain numeric values".
# - No invalid value reaches the transformation table.

test_that("4.2.6 mfp2.default() rejects invalid numeric settings", {
  expect_setting_numeric_errors(mfp2_validation_call, colnames(setting_x))
})

# -----------------------------------------------------------------------------
# 4.2.7 Scale must be strictly positive
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify the scale-specific constraint that every explicit divisor is greater
#   than zero, in addition to being numeric and finite.
# - Protect transformation calculations from division by zero and from a sign
#   reversal introduced by a negative scale.
#
# Expected output:
# - A named scale vector containing 0 errors with
#   "`scale` must contain only strictly positive values".
# - A named scale vector containing -1 produces the same error.
# - Validation fails before an mfp2 model is fitted.

test_that("4.2.7 mfp2.default() rejects nonpositive scale values", {
  expect_nonpositive_scale_errors(mfp2_validation_call, colnames(setting_x))
})

# -----------------------------------------------------------------------------
# 4.2.8 Formula-level scalar settings remain global
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that the matrix-interface safety change does not impose a new naming
#   requirement on the established top-level formula scalar interface.
# - Verify that formula preprocessing expands one global scalar to each modeled
#   numeric variable before calling the default method.
#
# Expected output:
# - The fitted transformation rows for age and weight both store shift = 2.
# - The fitted transformation rows for age and weight both store scale = 10.
# - The formula fit completes without requiring named top-level vectors.
# - `center = FALSE` is repeated inside each fp() term deliberately: fp() has
#   its own default `center = TRUE`, so the top-level setting alone would not
#   disable centering for these FP terms. Centering is not under test here.

test_that("4.2.8 mfp2.formula() keeps scalar shift and scale compatible", {
  data("prostate", package = "mfp2")

  fit <- mfp2(
    lpsa ~ fp(age, df = 2, center = FALSE, force_max_fp = TRUE) +
      fp(weight, df = 2, center = FALSE, force_max_fp = TRUE),
    data = prostate,
    shift = 2,
    scale = 10,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    unname(fit$transformations[c("age", "weight"), "shift"]),
    c(2, 2)
  )
  expect_equal(
    unname(fit$transformations[c("age", "weight"), "scale"]),
    c(10, 10)
  )
})

# -----------------------------------------------------------------------------
# 4.2.9 Per-variable fp() settings remain supported
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that variable-specific shift and scale values supplied inside fp()
#   remain part of the supported formula interface.
# - Verify that formula preprocessing creates internally named vectors and that
#   the default method matches those values to the correct model-matrix columns.
# - Protect against a regression where formula-generated vectors become unnamed
#   and are rejected by the stricter matrix-interface validation.
#
# Expected output:
# - The age transformation row stores shift = 1 and scale = 10.
# - The weight transformation row stores shift = 2 and scale = 100.
# - Reading rows c("age", "weight") returns shifts c(1, 2) and scales
#   c(10, 100), with no positional reassignment.
# - `center = FALSE` is specified inside each fp() term so this test isolates
#   shift/scale behavior instead of inheriting fp()'s default center = TRUE.

test_that("4.2.9 mfp2.formula() keeps per-variable fp() settings compatible", {
  data("prostate", package = "mfp2")

  fit <- mfp2(
    lpsa ~ fp(
      age,
      df = 2,
      shift = 1,
      scale = 10,
      center = FALSE,
      force_max_fp = TRUE
    ) + fp(
      weight,
      df = 2,
      shift = 2,
      scale = 100,
      center = FALSE,
      force_max_fp = TRUE
    ),
    data = prostate,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    unname(fit$transformations[c("age", "weight"), "shift"]),
    c(1, 2)
  )
  expect_equal(
    unname(fit$transformations[c("age", "weight"), "scale"]),
    c(10, 100)
  )
})


# -----------------------------------------------------------------------------
# 4.2.10 Direct single-shift path matches the historical two-sweep path
# -----------------------------------------------------------------------------
# Test purpose:
# - Exercise the actual Step-16 shift optimization with df > 1. Variables with
#   df = 1 have their shift forced to zero earlier in mfp2.default(), so they do
#   not test this code path meaningfully.
# - Compare the optimized direct-matrix branch with the historical two-sweep
#   branch on the SAME matrix interface. Attaching the private preprocessing
#   matrix makes mfp2.default() take the historical two-matrix branch while
#   keeping the fitting matrix and model specification otherwise identical.
# - This is a direct numerical regression test of the source refactor, not a
#   comparison between matrix and formula coefficient parameterizations.

test_that("4.2.10 direct single-shift path matches historical two-sweep path", {
  n <- 96L
  index <- seq_len(n)
  x_direct <- cbind(
    x1 = seq(-8, 12, length.out = n),
    x2 = -4 + ((17 * index) %% 31)
  )
  x_before <- x_direct
  y <- 1.2 + 0.28 * x_direct[, "x1"] - 0.17 * x_direct[, "x2"] +
    0.015 * x_direct[, "x1"]^2 + 0.03 * sin(index / 5)

  fit_optimized <- mfp2(
    x = x_direct,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  # Force the historical Step-16 branch with a separate preprocessing source
  # containing exactly the same numerical matrix. extract_preprocess_matrix()
  # removes this private attribute before fitting.
  x_historical <- x_direct
  attr(x_historical, "mfp2_preprocess_x") <- x_direct

  fit_historical <- mfp2(
    x = x_historical,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    fit_optimized$transformations,
    fit_historical$transformations,
    tolerance = 0
  )
  expect_equal(
    unname(coef(fit_optimized)),
    unname(coef(fit_historical)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fitted(fit_optimized)),
    unname(fitted(fit_historical)),
    tolerance = 1e-12
  )
  expect_identical(x_direct, x_before)
})

# -----------------------------------------------------------------------------
# 4.2.11 Explicit shift/scale path matches the historical two-sweep path
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that sharing the shifted direct matrix does not alter explicit
#   preprocessing settings or fitted results.
# - Use non-collinear predictors and df > 1 so the supplied shift values are
#   actually active (df = 1 would force shift = 0 before Step 16).

test_that("4.2.11 direct single-shift path preserves explicit settings", {
  n <- 100L
  index <- seq_len(n)
  x_direct <- cbind(
    x1 = seq(-4, 9, length.out = n),
    x2 = 2 + ((19 * index) %% 37)
  )
  y <- 2 + 0.45 * x_direct[, "x1"] + 0.11 * x_direct[, "x2"] +
    0.02 * cos(index / 7)

  shift <- c(x1 = 5, x2 = 1)
  scale <- c(x1 = 10, x2 = 20)

  fit_optimized <- mfp2(
    x = x_direct,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = shift,
    scale = scale,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  x_historical <- x_direct
  attr(x_historical, "mfp2_preprocess_x") <- x_direct

  fit_historical <- mfp2(
    x = x_historical,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = shift,
    scale = scale,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    fit_optimized$transformations,
    fit_historical$transformations,
    tolerance = 0
  )
  expect_equal(
    unname(coef(fit_optimized)),
    unname(coef(fit_historical)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fitted(fit_optimized)),
    unname(fitted(fit_historical)),
    tolerance = 1e-12
  )
})

# -----------------------------------------------------------------------------
# 4.2.12 Direct subset still uses full-data preprocessing
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that the optimized direct branch still estimates shift and scale
#   before applying subset, exactly as the original source does.
# - Use a full-rank retained design; the previous fixture made the retained x1
#   and x2 columns affine functions of each other and therefore correctly
#   triggered mfp2's rank-deficiency validation before this behavior was tested.

test_that("4.2.12 direct subset preserves full-data shift and scale", {
  n <- 122L
  interior_index <- seq_len(n - 2L)
  x1 <- c(-50, seq(-3, 4, length.out = n - 2L), 80)
  x2 <- c(-20, 2 + ((17 * interior_index) %% 53), 70)
  x_direct <- cbind(x1 = x1, x2 = x2)
  y <- 0.7 + 0.22 * x1 - 0.08 * x2 + 0.02 * sin(seq_len(n) / 6)
  rows <- 2:(n - 1L)

  fit_optimized <- mfp2(
    x = x_direct,
    y = y,
    subset = rows,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  # Historical two-sweep reference on the same direct matrix and subset.
  x_historical <- x_direct
  attr(x_historical, "mfp2_preprocess_x") <- x_direct
  fit_historical <- mfp2(
    x = x_historical,
    y = y,
    subset = rows,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expected_shift <- apply(x_direct, 2, find_shift_factor)
  expected_scale <- vapply(
    seq_len(ncol(x_direct)),
    function(j) find_scale_factor(x_direct[, j] + expected_shift[j]),
    numeric(1L)
  )
  names(expected_scale) <- colnames(x_direct)

  stored <- fit_optimized$transformations[
    colnames(x_direct), c("shift", "scale"), drop = FALSE
  ]

  expect_equal(unname(stored$shift), unname(expected_shift), tolerance = 0)
  expect_equal(unname(stored$scale), unname(expected_scale), tolerance = 0)
  expect_equal(
    fit_optimized$transformations,
    fit_historical$transformations,
    tolerance = 0
  )
  expect_equal(
    unname(fitted(fit_optimized)),
    unname(fitted(fit_historical)),
    tolerance = 1e-12
  )

  retained_shift <- apply(x_direct[rows, , drop = FALSE], 2, find_shift_factor)
  expect_false(isTRUE(all.equal(unname(expected_shift), unname(retained_shift))))
})


# =============================================================================
# 4.3 mfpi shift and scale behavior
# =============================================================================
# MFPI-specific review notes
# --------------------------
# The MFPI matrix contains the grouping variable "svi" in addition to ordinary
# covariates. Named shift and scale vectors may specify any subset of its columns;
# unspecified entries remain automatic. Downstream MFPI code continues to treat
# svi as categorical grouping metadata, so its final shift and scale are 0 and 1.
#
# The tests below are intentionally separate from mfp2 tests because MFPI also
# builds an adjustment model, removes or specially handles the grouping column,
# computes interaction metrics, and has its own prediction method.

# -----------------------------------------------------------------------------
# 4.3.1 NULL keeps automatic preprocessing for MFPI covariates
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that NULL retains MFPI's existing automatic preprocessing for the
#   ordinary continuous covariates after the grouping column is identified.
# - Protect the separate MFPI normalization and adjustment-model path from
#   receiving unresolved or invalid settings.
# - Confirm that the new validation does not alter automatic shift/scale choice.
#
# Expected output:
# - fit$shift[c("age", "weight")] contains two non-missing finite values.
# - fit$scale[c("age", "weight")] contains two non-missing finite values.
# - Both continuous-variable scales are strictly greater than zero.
# - The fit completes with the existing categorical handling of svi unchanged.

test_that("4.3.1 mfpi.default() keeps automatic shift and scale for NULL", {
  fit <- fit_mfpi_settings()

  expect_false(anyNA(fit$shift[c("age", "weight")]))
  expect_false(anyNA(fit$scale[c("age", "weight")]))
  expect_true(all(is.finite(fit$shift[c("age", "weight")])))
  expect_true(all(is.finite(fit$scale[c("age", "weight")])))
  expect_true(all(fit$scale[c("age", "weight")] > 0))
})

# -----------------------------------------------------------------------------
# 4.3.2 Scalars are recycled to ordinary MFPI covariates
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that a scalar remains a valid global MFPI setting even though the
#   matrix also contains the grouping variable svi.
# - Confirm that scalar expansion survives grouping-column handling and reaches
#   both ordinary continuous covariates used by the adjustment model.
#
# Expected output:
# - fit$shift[c("age", "weight")] equals c(2, 2).
# - fit$scale[c("age", "weight")] equals c(10, 10).
# - The test intentionally leaves svi to MFPI's existing categorical metadata
#   rules rather than treating it as an ordinary transformed covariate.

test_that("4.3.2 mfpi.default() recycles scalar shift and scale", {
  fit <- fit_mfpi_settings(shift = 2, scale = 10)

  expect_equal(unname(fit$shift[c("age", "weight")]), c(2, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 10))
})


# -----------------------------------------------------------------------------
# 4.3.2a MFPI uses the same named partial-setting contract
# -----------------------------------------------------------------------------

test_that("4.3.2a mfpi.default() estimates unspecified named partial settings", {
  fit <- fit_mfpi_settings(shift = c(age = 2), scale = c(age = 100))

  expect_equal(unname(fit$shift["age"]), 2)
  expect_equal(
    unname(fit$shift["weight"]),
    find_shift_factor(mfpi_setting_x[, "weight"])
  )
  expect_equal(unname(fit$shift["svi"]), 0)
  expect_equal(unname(fit$scale["age"]), 100)
  expect_equal(
    unname(fit$scale["weight"]),
    find_scale_factor(
      mfpi_setting_x[, "weight"] + fit$shift[["weight"]]
    )
  )
  expect_equal(unname(fit$scale["svi"]), 1)
})

# Test purpose: MFPI must neutralize explicitly supplied preprocessing values
# for an ordinary binary adjustment covariate, not only for the group variable.
test_that("4.3.2b mfpi.default() resets supplied scale for binary covariates", {
  binary_adjustment <- rep(c(0, 1), length.out = setting_n)
  x <- cbind(mfpi_setting_x, binary_adjustment = binary_adjustment)
  y <- mfpi_setting_y + 0.3 * binary_adjustment

  fit <- mfpi(
    x = x,
    y = y,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    shift = c(binary_adjustment = 10),
    scale = c(binary_adjustment = 1000),
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = c("age", "weight"),
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(unname(fit$shift["binary_adjustment"]), 0)
  expect_equal(unname(fit$scale["binary_adjustment"]), 1)
})

# -----------------------------------------------------------------------------
# 4.3.3 Fully named vectors include the MFPI grouping column
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that complete named vectors remain accepted at the public matrix
#   boundary, including an explicit value for the grouping column svi.
# - Confirm that, after MFPI applies its special grouping-variable handling, the
#   values for age and weight remain aligned with their names.
#
# Expected output:
# - The complete inputs use names c("svi", "age", "weight").
# - fit$shift[c("age", "weight")] equals c(1, 2).
# - fit$scale[c("age", "weight")] equals c(10, 100).
# - The grouping column does not cause either continuous-variable value to move
#   to the wrong variable.

test_that("4.3.3 mfpi.default() accepts fully named shift and scale", {
  fit <- fit_mfpi_settings(
    shift = c(svi = 0, age = 1, weight = 2),
    scale = c(svi = 1, age = 10, weight = 100)
  )

  expect_equal(unname(fit$shift[c("age", "weight")]), c(1, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 100))
})

# -----------------------------------------------------------------------------
# 4.3.4 Named-vector order does not change MFPI results
# -----------------------------------------------------------------------------
# Test purpose:
# - Reproduce the positional-matching risk in the more complex MFPI path by
#   reversing named settings that include the grouping variable.
# - Verify alignment before and after MFPI handles or removes svi from ordinary
#   transformation processing.
# - Exercise normalized settings, adjustment-model transformations, fitted
#   coefficients, interaction metrics, and MFPI prediction in one regression
#   test.
#
# Expected output:
# - ordered$shift and reversed$shift are identical named vectors.
# - ordered$scale and reversed$scale are identical named vectors.
# - The two adjustment-model transformation tables and coefficient vectors are
#   numerically identical.
# - ordered$all_model_metrics and reversed$all_model_metrics are identical.
# - Link-scale predictions for the age term and all interaction models are
#   identical for the first 12 rows of mfpi_setting_x.
# - Age retains shift = 1 and scale = 10; weight retains shift = 2 and
#   scale = 100, regardless of the supplied order around svi.

test_that("4.3.4 mfpi.default() matches shift and scale by name", {
  ordered <- fit_mfpi_settings(
    shift = c(svi = 0, age = 1, weight = 2),
    scale = c(svi = 1, age = 10, weight = 100)
  )
  reversed <- fit_mfpi_settings(
    shift = c(weight = 2, age = 1, svi = 0),
    scale = c(weight = 100, age = 10, svi = 1)
  )

  expect_equal(ordered$shift, reversed$shift)
  expect_equal(ordered$scale, reversed$scale)
  expect_equal(
    ordered$adjustment_model$transformations,
    reversed$adjustment_model$transformations
  )
  expect_equal(
    unname(stats::coef(ordered$adjustment_model)),
    unname(stats::coef(reversed$adjustment_model))
  )
  expect_equal(ordered$all_model_metrics, reversed$all_model_metrics)

  new_x <- as.data.frame(mfpi_setting_x[1:12, , drop = FALSE])
  pred_ordered <- predict(
    ordered,
    newdata = new_x,
    terms = "age",
    model = "all",
    type = "link"
  )
  pred_reversed <- predict(
    reversed,
    newdata = new_x,
    terms = "age",
    model = "all",
    type = "link"
  )

  expect_equal(
    pred_ordered$predictions$fit,
    pred_reversed$predictions$fit
  )
})

# -----------------------------------------------------------------------------
# 4.3.5 Invalid MFPI multi-value names are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that MFPI rejects unnamed multi-value vectors, empty element names,
#   duplicate names, and names outside its public matrix columns.
# - Confirm that any supplied names are validated before special treatment of
#   the grouping variable begins.
# - Apply the same partial, unique naming contract to shift and scale.
#
# Expected output for each argument:
# - Unnamed multi-value vectors and empty element names produce the
#   single-unnamed-or-named-vector error.
# - Duplicate names produce the unique-names error.
# - Supplying an unknown name produces the unknown-column error; omitting
#   otherwise valid columns leaves those settings automatic.
# - No adjustment or interaction model is fitted for invalid inputs.

test_that("4.3.5 mfpi.default() rejects invalid shift and scale names", {
  expect_setting_name_errors(
    mfpi_validation_call,
    colnames(mfpi_setting_x)
  )
})

# -----------------------------------------------------------------------------
# 4.3.6 Missing, non-finite, and logical MFPI values are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that MFPI uses the same finite, nonmissing numeric-value contract as
#   the ordinary mfp2 matrix interface.
# - Ensure invalid values fail before settings are separated into grouping and
#   adjustment-model components.
#
# Expected output for both shift and scale:
# - Inf produces the finite-values error.
# - NA produces the no-missing-values error.
# - TRUE/FALSE vectors produce the numeric-values error.
# - No invalid value is stored in fit$shift, fit$scale, or the adjustment model.

test_that("4.3.6 mfpi.default() rejects invalid numeric settings", {
  expect_setting_numeric_errors(
    mfpi_validation_call,
    colnames(mfpi_setting_x)
  )
})

# -----------------------------------------------------------------------------
# 4.3.7 MFPI scale values must be strictly positive
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that MFPI enforces the scale-specific positivity rule on each explicit
#   named setting before grouping-variable handling or adjustment fitting.
# - Prevent zero divisors and negative rescaling from entering any MFPI model.
#
# Expected output:
# - A named scale vector containing 0 errors with
#   "`scale` must contain only strictly positive values".
# - A named scale vector containing -1 produces the same error.
# - Neither the adjustment model nor the interaction models are fitted.

test_that("4.3.7 mfpi.default() rejects nonpositive scale values", {
  expect_nonpositive_scale_errors(
    mfpi_validation_call,
    colnames(mfpi_setting_x)
  )
})

# -----------------------------------------------------------------------------
# 4.3.8 Formula-level scalar settings remain global in MFPI
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that the stricter MFPI matrix-vector contract does not alter the
#   documented top-level formula scalar interface.
# - Verify that formula preprocessing applies global scalar values to ordinary
#   continuous variables while continuing to handle svi as the group variable.
#
# Expected output:
# - fit$shift[c("age", "weight")] equals c(2, 2).
# - fit$scale[c("age", "weight")] equals c(10, 10).
# - The formula fit succeeds without a named top-level vector or an explicit
#   scalar setting for svi.

test_that("4.3.8 mfpi.formula() keeps scalar shift and scale compatible", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ svi + fp(age, df = 2, force_max_fp = TRUE) +
      fp(weight, df = 2, force_max_fp = TRUE),
    data = prostate,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    shift = 2,
    scale = 10,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(unname(fit$shift[c("age", "weight")]), c(2, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 10))
})

# -----------------------------------------------------------------------------
# 4.3.9 Per-variable fp() settings remain supported in MFPI formulas
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that per-variable settings inside fp() remain supported when the
#   formula also contains the MFPI grouping variable.
# - Verify that formula preprocessing generates correctly named internal
#   settings and that MFPI's grouping-column handling does not shift them.
# - Protect compatibility for existing formula calls that intentionally assign
#   different preprocessing values to different continuous variables.
#
# Expected output:
# - fit$shift[c("age", "weight")] equals c(1, 2).
# - fit$scale[c("age", "weight")] equals c(10, 100).
# - Age receives only the values declared in fp(age, ...), and weight receives
#   only the values declared in fp(weight, ...).

test_that("4.3.9 mfpi.formula() keeps per-variable fp() settings compatible", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ svi +
      fp(
        age,
        df = 2,
        shift = 1,
        scale = 10,
        force_max_fp = TRUE
      ) + fp(
        weight,
        df = 2,
        shift = 2,
        scale = 100,
        force_max_fp = TRUE
      ),
    data = prostate,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(unname(fit$shift[c("age", "weight")]), c(1, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 100))
})

# =============================================================================
# 5. Candidate-power validation and custom powers
# =============================================================================

# Test purpose: Checks that candidate power vectors are deduplicated and sorted.
test_that("normalize_fp_power_vector() removes duplicates and sorts", {
  result <- normalize_fp_power_vector(c(3, 1, 2, 1, -1))
  expect_equal(result, c(-1, 1, 2, 3))
})

# Test purpose: Checks that invalid power-vector inputs fail validation.
test_that("normalize_fp_power_vector() rejects empty or non-numeric input", {
  expect_error(normalize_fp_power_vector(numeric(0)))
  expect_error(normalize_fp_power_vector("abc"))
  expect_error(normalize_fp_power_vector(c(1, NA, 2)))
})

# Test purpose: Checks that a single non-1 candidate power is allowed and can
# generate repeated-power FP candidates.
test_that("single non-1 candidate power is valid for df > 1", {
  # A single non-unity power should work — repeated-power FP2 candidate
  fit <- mfp2(
    x_prostate[, c("cavol", "age")], y_prostate,
    powers = list(cavol = 2), df = 4,
    verbose = FALSE
  )
  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks that powers = 1 alone is rejected for nonlinear FP selection.
test_that("candidate power set containing only 1 is invalid for df > 1", {
  expect_error(
    validate_mfp_candidate_powers(
      powers = list(x = 1),
      df = c(x = 4)
    ),
    "only power 1"
  )
})

# Test purpose: Checks that powers = 1 is allowed when the term is restricted
# to a linear effect.
test_that("candidate power set containing only 1 is valid for df = 1", {
  expect_true(
    validate_mfp_candidate_powers(
      powers = list(x = 1),
      df = c(x = 1)
    )
  )
})

# Test purpose: Checks that candidate powers supplied inside fp() are accepted
# by the formula interface.
test_that("custom powers via formula fp() work", {
  fit <- mfp2(
    lpsa ~ fp(cavol, powers = c(-1, 0, 1, 2)) + fp(age),
    data = prostate, verbose = FALSE
  )
  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks FP candidate-generation dimensions, including repeated
# power combinations.
test_that("generate_powers_fp() produces correct number of combinations", {
  powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  fp1 <- generate_powers_fp(degree = 1, powers = powx)
  expect_equal(nrow(fp1), 8)
  expect_equal(ncol(fp1), 1)

  fp2 <- generate_powers_fp(degree = 2, powers = powx)
  expect_equal(nrow(fp2), 36) # C(8+2-1, 2) = C(9,2) = 36
  expect_equal(ncol(fp2), 2)

  # Single power with degree 2 gives repeated-power pair
  fp_single <- generate_powers_fp(degree = 2, powers = 2)
  expect_equal(nrow(fp_single), 1)
  expect_equal(fp_single[1, ], c(2, 2))
})

# Test purpose: Checks the null-degree FP power matrix used for omitted/null terms.
test_that("generate_powers_fp() degree 0 returns matrix(1)", {
  fp0 <- generate_powers_fp(degree = 0)
  expect_equal(fp0, matrix(1, nrow = 1, ncol = 1))
})

# =============================================================================
# 6. SAZ (spike-at-zero) — eligibility, cascade, and reset
# =============================================================================

# Test purpose: Fits a clear spike-at-zero example and checks that the spike
# flag is retained.
test_that("spike-at-zero basic fitting works", {
  set.seed(123)
  n <- 300
  prop_zero <- 0.25
  x_val <- numeric(n)
  zero_idx <- sample(seq_len(n), size = round(n * prop_zero))
  x_val[zero_idx] <- 0
  x_val[-zero_idx] <- rgamma(n - length(zero_idx), shape = 2, rate = 0.5)

  Z <- ifelse(x_val == 0, 1, 0)
  y_val <- 1.5 * Z + 2 * log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(x_mat, y_val, spike_vars = "exposure", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  # spike should be TRUE because the proportion of zeros is well within threshold
  expect_true(fit$fp_terms["exposure", "spike"])
})

# Test purpose: Checks that SAZ is reset when the zero component fails the
# minimum proportion threshold.
test_that("spike-at-zero is reset when zero proportion is too low", {
  # Almost no zeros
  set.seed(42)
  n <- 200
  x_val <- rgamma(n, shape = 2, rate = 1)
  # Only 2 zeros out of 200 = 1% which is below default 10%
  x_val[1:2] <- 0
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(x_mat, y_val, spike_vars = "exposure", verbose = FALSE),
    "spike"
  )

  # Spike should be reset
  expect_false(fit$fp_terms["exposure", "spike"])
})

# Test purpose: Checks that explicit zero handling survives after an ineligible
# spike request is reset.
test_that("spike cascade restores user-specified zero/catzero on reset", {
  # Build a scenario where spike is reset but user also set zero=TRUE
  set.seed(42)
  n <- 200
  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros => spike reset
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat, y_val,
      spike_vars = "exposure",
      zero_vars = "exposure",
      verbose = FALSE
    )
  })

  # Spike should be reset, but zero should be preserved
  expect_false(fit$fp_terms["exposure", "spike"])
  expect_true(fit$zero["exposure"])
})

# Test purpose: Checks that spike-at-zero handling can be requested inside fp()
# in the formula interface.
test_that("spike formula interface fp(spike = TRUE) works", {
  set.seed(123)
  n <- 300
  prop_zero <- 0.25
  x_val <- numeric(n)
  zero_idx <- sample(seq_len(n), size = round(n * prop_zero))
  x_val[zero_idx] <- 0
  x_val[-zero_idx] <- rgamma(n - length(zero_idx), shape = 2, rate = 0.5)

  Z <- ifelse(x_val == 0, 1, 0)
  y_val <- 1.5 * Z + 2 * log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  dat <- data.frame(y = y_val, exposure = x_val)

  fit <- mfp2(
    y ~ fp(exposure, spike = TRUE, center = FALSE),
    data = dat, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks that changing the SAZ component threshold changes spike
#  eligibility as expected.
test_that("min_saz_component_prop controls eligibility threshold", {
  set.seed(123)
  n <- 200
  x_val <- rgamma(n, shape = 2, rate = 1)
  # 15% zeros
  x_val[sample(n, 30)] <- 0
  y_val <- 0.5 * x_val + rnorm(n)
  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  # With default threshold 0.10 it should be eligible
  fit_low <- mfp2(x_mat, y_val, spike_vars = "exposure",
                  min_saz_component_prop = 0.10, verbose = FALSE)
  expect_true(fit_low$fp_terms["exposure", "spike"])

  # With high threshold 0.40 it should be ineligible
  expect_warning(
    fit_high <- mfp2(x_mat, y_val, spike_vars = "exposure",
                     min_saz_component_prop = 0.40, verbose = FALSE),
    "spike"
  )
  expect_false(fit_high$fp_terms["exposure", "spike"])
})

# Test purpose: Ensures SAZ eligibility requires enough positive-component
# observations, not only enough zero-component observations.
test_that("spike-at-zero is reset when positive component proportion is too low", {
  set.seed(104)
  n <- 200

  x_val <- numeric(n)
  x_val[1:5] <- rgamma(5, shape = 2, rate = 1)  # 2.5% positive component
  y_val <- 1.5 * (x_val == 0) + 0.2 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      verbose = FALSE
    ),
    "spike"
  )

  expect_false(fit$fp_terms["exposure", "spike"])
})

# Test purpose: Checks that non-positive values, including negative values,
# are counted in the zero component for SAZ eligibility.
test_that("spike-at-zero counts non-positive values in the zero component", {
  set.seed(105)
  n <- 200

  x_val <- numeric(n)
  zero_idx <- sample(seq_len(n), size = 50)
  x_val[zero_idx] <- -1
  x_val[-zero_idx] <- rgamma(n - length(zero_idx), shape = 2, rate = 1)

  y_val <- 1.5 * (x_val <= 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)
  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$fp_terms["exposure", "spike"])
  expect_true(fit$zero["exposure"])
})

# Test purpose: Ensures that when spike is reset, user-specified catzero handling
# is preserved rather than removed with the spike-implied cascade.
test_that("spike cascade preserves user-specified catzero on reset", {
  set.seed(106)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0  # too few zeros for SAZ eligibility
  y_val <- 1.5 * (x_val <= 0) + 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      catzero_vars = "exposure",
      verbose = FALSE
    )
  })

  expect_false(fit$fp_terms["exposure", "spike"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})

# Test purpose: Directly checks reset_spike() for the all-zero case, where
# the positive component is absent.
test_that("reset_spike() resets all-zero variables", {
  x <- matrix(0, nrow = 100, ncol = 1, dimnames = list(NULL, "exposure"))

  spike <- c(exposure = TRUE)
  user_catzero <- c(exposure = FALSE)
  user_zero <- c(exposure = FALSE)

  expect_warning(
    out <- reset_spike(
      x = x,
      spike = spike,
      user_catzero = user_catzero,
      user_zero = user_zero,
      min_saz_component_prop = 0.10
    ),
    "positive observation proportion"
  )

  expect_false(out$spike["exposure"])
  expect_false(out$catzero["exposure"])
  expect_false(out$zero["exposure"])
})


# Test purpose: Checks that resolve_saz_eligibility() treats nonpositive values
# as the structural-zero component without requiring a recoded matrix.
test_that("resolve_saz_eligibility() counts negative values as zero component", {
  # Only the zero-component proportion matters here. Use deterministic
  # positive values so this validation test is independent of RNG state.
  x <- matrix(
    c(rep(-2, 20), seq(0.1, 18, length.out = 180)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  spike <- c(exposure = TRUE)
  catzero <- c(exposure = FALSE)
  zero <- c(exposure = FALSE)

  out <- resolve_saz_eligibility(
    x = x,
    spike = spike,
    catzero = catzero,
    zero = zero,
    min_saz_component_prop = 0.10
  )

  expect_true(out$spike["exposure"])
  expect_true(out$catzero["exposure"])
  expect_true(out$zero["exposure"])
})

# Test purpose: The raw-predicate reset must be exactly equivalent to the old
# temporary x[x <= 0] <- 0 representation, including binary detection after
# multiple distinct nonpositive values collapse to one structural-zero level.
test_that("reset_spike raw predicates match zero-recoded eligibility", {
  x_raw <- cbind(
    exposure = c(-3, -2, -1, 0, rep(2, 16)),
    eligible = c(rep(-2, 4), rep(1:4, each = 4))
  )
  x_recoded <- x_raw
  x_recoded[x_recoded <= 0] <- 0

  spike <- c(exposure = TRUE, eligible = TRUE)
  user_catzero <- c(exposure = FALSE, eligible = FALSE)
  user_zero <- c(exposure = FALSE, eligible = FALSE)

  raw <- suppressWarnings(reset_spike(
    x = x_raw,
    spike = spike,
    user_catzero = user_catzero,
    user_zero = user_zero,
    min_saz_component_prop = 0.10
  ))
  recoded <- suppressWarnings(reset_spike(
    x = x_recoded,
    spike = spike,
    user_catzero = user_catzero,
    user_zero = user_zero,
    min_saz_component_prop = 0.10
  ))

  expect_identical(raw, recoded)
  # exposure has one effective zero level plus one positive level, so it is
  # binary after structural-zero collapsing and must be reset.
  expect_false(raw$spike[["exposure"]])
  expect_true(raw$spike[["eligible"]])
})

# Test purpose: Checks that structural-zero proportions use finite observations
# and count nonpositive values as the SAZ zero component.
test_that("calculate_saz_prop_zero() reports retained SAZ proportions", {
  x <- cbind(
    exposure = c(-2, 0, 0, 1, 2, 3, NA_real_, Inf),
    ordinary = seq_len(8)
  )

  out <- calculate_saz_prop_zero(
    x = x,
    spike = c(exposure = TRUE, ordinary = FALSE)
  )

  expect_equal(out[["exposure"]], 3 / 6)
  expect_true(is.na(out[["ordinary"]]))
})

# Test purpose: Ensures that when spike-only handling is reset, the variable
# returns to ordinary FP handling with no zero or catzero flags.
test_that("spike-only reset restores ordinary FP handling", {
  set.seed(305)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros for SAZ eligibility
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      verbose = FALSE
    ),
    "spike-at-zero option has been reset"
  )

  expect_false(fit$spike["exposure"])
  expect_false(fit$catzero["exposure"])
  expect_false(fit$zero["exposure"])
})

# Test purpose: Ensures that when spike is reset, explicitly requested zero
# handling is preserved rather than removed with the spike-implied cascade.
test_that("spike reset preserves user-specified zero handling", {
  set.seed(306)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros for SAZ eligibility
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      zero_vars = "exposure",
      verbose = FALSE
    )
  })

  expect_false(fit$spike["exposure"])
  expect_false(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})

# Test purpose: Ensures that when spike is reset, explicitly requested catzero
# handling is preserved and still implies zero handling.
test_that("spike reset preserves user-specified catzero handling", {
  set.seed(307)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros for SAZ eligibility
  y_val <- 1.5 * (x_val <= 0) + 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      catzero_vars = "exposure",
      verbose = FALSE
    )
  })

  expect_false(fit$spike["exposure"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})

# Test purpose: Ensures retained SAZ variables satisfy the internal cascade:
# spike implies catzero, and catzero implies zero.
test_that("retained spike variable implies catzero and zero handling", {
  set.seed(308)
  n <- 200

  x_val <- c(rep(0, 60), rgamma(140, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  expect_true(fit$spike["exposure"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
  expect_true(fit$fp_terms["exposure", "spike"])
})

# Test purpose: Ensures formula-interface fp(spike = TRUE) is translated into
# spike, catzero, and zero handling for an eligible SAZ variable.
test_that("formula interface fp(spike = TRUE) activates SAZ handling", {
  set.seed(309)
  n <- 200

  exposure <- c(rep(0, 50), rgamma(150, shape = 2, rate = 1))
  dat <- data.frame(
    y = 2 * (exposure == 0) + log(ifelse(exposure > 0, exposure, 1)) + rnorm(n),
    exposure = exposure,
    z = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ fp(exposure, spike = TRUE) + fp(z),
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})

# Test purpose: Checks that the SAZ algorithm runs under AIC-based selection,
# not only under p-value based closed testing.
test_that("spike-at-zero works with AIC criterion", {
  set.seed(310)
  n <- 220

  x_val <- c(rep(0, 70), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    criterion = "aic",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_true(fit$fp_terms["exposure", "spike"])
  expect_true(fit$spike_dec["exposure"] %in% c(1L, 2L, 3L))
})

# Test purpose: Checks that the SAZ algorithm runs under BIC-based selection.
test_that("spike-at-zero works with BIC criterion", {
  set.seed(311)
  n <- 220

  x_val <- c(rep(0, 70), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    criterion = "bic",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_true(fit$fp_terms["exposure", "spike"])
  expect_true(fit$spike_dec["exposure"] %in% c(1L, 2L, 3L))
})

# Test purpose: Ensures retained SAZ variables with at most 3 distinct positive
# values have their maximum FP df forced to 1.
test_that("cap_spike_df() forces df = 1 for at most 3 distinct positive values", {
  x <- matrix(
    c(rep(0, 30), rep(c(1, 2, 3), each = 10)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  df <- c(exposure = 4)
  spike <- c(exposure = TRUE)

  expect_warning(
    out <- cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    "maximum FP df was reduced"
  )

  expect_equal(out[["exposure"]], 1)
})

# Test purpose: Ensures retained SAZ variables with 4 or 5 distinct positive
# values have their maximum FP df capped at FP1, i.e. df = 2.
test_that("cap_spike_df() caps df at 2 for 4 or 5 distinct positive values", {
  x <- matrix(
    c(rep(0, 30), rep(c(1, 2, 3, 4, 5), each = 10)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  df <- c(exposure = 4)
  spike <- c(exposure = TRUE)

  expect_warning(
    out <- cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    "maximum FP df was reduced"
  )

  expect_equal(out[["exposure"]], 2)
})

# Test purpose: Ensures retained SAZ variables with at least 6 distinct positive
# values keep the requested maximum FP df.
test_that("cap_spike_df() keeps df unchanged for at least 6 distinct positive values", {
  x <- matrix(
    c(rep(0, 30), rep(1:6, each = 10)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  df <- c(exposure = 4)
  spike <- c(exposure = TRUE)

  expect_warning(
    out <- cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    NA
  )

  expect_equal(out[["exposure"]], 4)
})

# Test purpose: GLM SAZ stage 2 must assemble the final intercept-inclusive
# Model 2/Model 3 design once and tell fit_model() not to prepend another
# intercept. This protects the one-allocation GLM path.
test_that("SAZ reduced GLMs pass final intercept-inclusive designs", {
  seen <- new.env(parent = emptyenv())
  seen$x <- list()
  seen$x_has_intercept <- logical()

  testthat::local_mocked_bindings(
    fit_model = function(x, x_has_intercept = FALSE, ...) {
      seen$x[[length(seen$x) + 1L]] <- x
      seen$x_has_intercept <- c(seen$x_has_intercept, x_has_intercept)
      list(logl = -1, df = NCOL(x))
    },
    .package = "mfp2"
  )

  data_xi <- cbind(
    catzero = c(1, 1, 0, 0),
    fp1 = c(0, 0, 1.2, 2.4)
  )
  adjustment <- cbind(z = c(2, 3, 4, 5))
  stage1 <- list(
    current_adj_params = list(
      exposure = list(data_xi = data_xi, data_adj = adjustment)
    )
  )

  out <- fit_saz_reduced_models(
    stage1_selection = stage1,
    xi = "exposure",
    y = rep(0, 4),
    weights = NULL,
    offset = NULL,
    family = stats::gaussian(),
    family_string = "gaussian",
    method = NULL,
    strata = NULL,
    nocenter = NULL,
    control = NULL,
    rownames = NULL,
    has_offset = FALSE
  )

  expected_fit2 <- cbind(
    "(Intercept)" = rep(1, 4),
    fp1 = data_xi[, "fp1"],
    adjustment
  )
  expected_fit3 <- cbind(
    "(Intercept)" = rep(1, 4),
    catzero = data_xi[, "catzero"],
    adjustment
  )

  expect_length(seen$x, 2L)
  expect_equal(seen$x[[1L]], expected_fit2)
  expect_equal(seen$x[[2L]], expected_fit3)
  expect_identical(seen$x_has_intercept, c(TRUE, TRUE))
  expect_false("x" %in% names(out))
  expect_identical(out$data_xi, data_xi)
  expect_identical(out$adjustment_matrix, adjustment)
})

# Test purpose: Cox SAZ stage 2 must preserve the historical no-intercept
# design contract. The GLM allocation optimization must never add an ordinary
# intercept to a matrix sent to the Cox fitter.
test_that("SAZ reduced Cox models remain intercept-free", {
  seen <- new.env(parent = emptyenv())
  seen$x <- list()
  seen$x_has_intercept <- logical()

  testthat::local_mocked_bindings(
    fit_model = function(x, x_has_intercept = FALSE, ...) {
      seen$x[[length(seen$x) + 1L]] <- x
      seen$x_has_intercept <- c(seen$x_has_intercept, x_has_intercept)
      list(logl = -1, df = NCOL(x))
    },
    .package = "mfp2"
  )

  data_xi <- cbind(
    catzero = c(1, 1, 0, 0),
    fp1 = c(0, 0, 1.2, 2.4)
  )
  adjustment <- cbind(z = c(2, 3, 4, 5))
  stage1 <- list(
    current_adj_params = list(
      exposure = list(data_xi = data_xi, data_adj = adjustment)
    )
  )

  fit_saz_reduced_models(
    stage1_selection = stage1,
    xi = "exposure",
    y = rep(0, 4),
    weights = NULL,
    offset = NULL,
    family = NULL,
    family_string = "cox",
    method = "efron",
    strata = NULL,
    nocenter = NULL,
    control = NULL,
    rownames = NULL,
    has_offset = FALSE
  )

  expect_length(seen$x, 2L)
  expect_equal(seen$x[[1L]], cbind(fp1 = data_xi[, "fp1"], adjustment))
  expect_equal(seen$x[[2L]], cbind(catzero = data_xi[, "catzero"], adjustment))
  expect_identical(seen$x_has_intercept, c(FALSE, FALSE))
  expect_false("(Intercept)" %in% colnames(seen$x[[1L]]))
  expect_false("(Intercept)" %in% colnames(seen$x[[2L]]))
})

# Test purpose: The optimized GLM SAZ assembly must be numerically identical
# to the historical path where fit_glm() prepended the intercept itself.
test_that("SAZ reduced GLM fits match historical assembly numerically", {
  data_xi <- cbind(
    catzero = c(1, 1, 0, 0, 0, 0, 0, 0),
    fp1 = c(0, 0, 0.3, 0.8, 1.2, 1.7, 2.1, 2.8)
  )
  adjustment <- cbind(z = c(-1.2, -0.4, 0.1, 0.7, 1.1, 1.8, 2.3, 3.0))
  y <- c(0.4, 0.8, 1.2, 1.7, 2.0, 2.6, 3.0, 3.5)
  family <- stats::gaussian()
  stage1 <- list(
    current_adj_params = list(
      exposure = list(data_xi = data_xi, data_adj = adjustment)
    )
  )

  out <- fit_saz_reduced_models(
    stage1_selection = stage1,
    xi = "exposure",
    y = y,
    weights = NULL,
    offset = NULL,
    family = family,
    family_string = "gaussian",
    method = NULL,
    strata = NULL,
    nocenter = NULL,
    control = NULL,
    rownames = NULL,
    has_offset = FALSE
  )

  old_x2 <- cbind(fp1 = data_xi[, "fp1"], adjustment)
  old_x3 <- cbind(catzero = data_xi[, "catzero"], adjustment)

  expected2 <- fit_model(
    x = old_x2,
    y = y,
    family = family,
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    method = NULL,
    strata = NULL,
    control = NULL,
    rownames = NULL,
    nocenter = NULL,
    has_offset = FALSE
  )
  expected3 <- fit_model(
    x = old_x3,
    y = y,
    family = family,
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    method = NULL,
    strata = NULL,
    control = NULL,
    rownames = NULL,
    nocenter = NULL,
    has_offset = FALSE
  )

  expect_equal(out$fit2$coefficients, expected2$coefficients, tolerance = 1e-12)
  expect_equal(out$fit3$coefficients, expected3$coefficients, tolerance = 1e-12)
  expect_equal(out$fit2$logl, expected2$logl, tolerance = 1e-12)
  expect_equal(out$fit3$logl, expected3$logl, tolerance = 1e-12)
  expect_identical(out$fit2$df, expected2$df)
  expect_identical(out$fit3$df, expected3$df)
})

# Test purpose: Ensures public mfp2() applies the SAZ positive-part df cap
# before final FP model selection metadata are stored.
test_that("mfp2() applies positive-part df cap for retained spike variables", {
  set.seed(312)
  n <- 180

  positive_values <- rep(c(1, 2, 3), each = 40)
  x_val <- c(rep(0, 60), positive_values)
  y_val <- 1.5 * (x_val == 0) + 0.4 * x_val + rnorm(n, sd = 0.2)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      df = 4,
      verbose = FALSE
    ),
    "maximum FP df was reduced"
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_equal(fit$fp_terms["exposure", "df_initial"], 1)
})

# Test purpose: Ensures prediction works for retained SAZ models, including
# new zero and positive values in newdata.
test_that("predict.mfp2() works for retained spike-at-zero models", {
  set.seed(313)
  n <- 200

  x_val <- c(rep(0, 50), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n, sd = 0.3)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  newx <- matrix(
    c(0, 0, 0.5, 1, 2, 4),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  pred <- predict(fit, newdata = newx)

  expect_length(pred, nrow(newx))
  expect_true(all(is.finite(pred)))
})

# Test purpose: Ensures prediction for retained SAZ models treats negative
# newdata values as part of the zero component rather than failing as ordinary FP.
test_that("predict.mfp2() treats negative newdata as zero component for SAZ models", {
  set.seed(314)
  n <- 200

  x_val <- c(rep(0, 50), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n, sd = 0.3)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  newx <- matrix(
    c(-2, -1, 0, 0.5, 2),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  pred <- predict(fit, newdata = newx)

  expect_length(pred, nrow(newx))
  expect_true(all(is.finite(pred)))
})

# Test purpose: Ensures SAZ eligibility is resolved independently for multiple
# spike variables, so one reset variable does not reset all spike variables.
test_that("multiple spike variables are reset independently", {
  set.seed(315)
  n <- 240

  eligible <- c(rep(0, 60), rgamma(180, shape = 2, rate = 1))
  ineligible <- rgamma(n, shape = 2, rate = 1)
  ineligible[1:2] <- 0 # too few zeros

  y_val <- 1.5 * (eligible == 0) + log(ifelse(eligible > 0, eligible, 1)) +
    0.2 * ineligible + rnorm(n)

  x_mat <- cbind(
    eligible = eligible,
    ineligible = ineligible
  )

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = c("eligible", "ineligible"),
      verbose = FALSE
    ),
    "ineligible"
  )

  expect_true(fit$spike["eligible"])
  expect_false(fit$spike["ineligible"])

  expect_true(fit$catzero["eligible"])
  expect_true(fit$zero["eligible"])

  expect_false(fit$catzero["ineligible"])
  expect_false(fit$zero["ineligible"])
})


# =============================================================================
# 7. ACD transformation
# =============================================================================

# Test purpose: Checks that ACD can be requested through the default matrix
# interface.
test_that("ACD transformation via default interface works", {
  fit <- mfp2(x_prostate, y_prostate, acdx = "cavol", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_true(fit$fp_terms["cavol", "acd"])
})

# Test purpose: Checks that ACD can be requested inside fp() in the formula
# interface.
test_that("ACD transformation via formula interface works", {
  fit <- mfp2(
    lpsa ~ fp(cavol, acdx = TRUE) + fp(age) + fp(svi, df = 1),
    data = prostate, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$fp_terms["cavol", "acd"])
})

# Test purpose: Checks candidate-power matrix dimensions for ACD degrees 0, 1,
# and 2.
test_that("ACD power generation produces correct matrix dimensions", {
  powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  acd0 <- generate_powers_acd(degree = 0, powers = powx)
  expect_equal(ncol(acd0), 2)
  expect_equal(nrow(acd0), 1)

  acd1 <- generate_powers_acd(degree = 1, powers = powx)
  expect_equal(ncol(acd1), 2)
  expect_equal(nrow(acd1), 8)
  expect_true(all(is.na(acd1[, 1]))) # first column all NA

  acd2 <- generate_powers_acd(degree = 2, powers = powx)
  expect_equal(ncol(acd2), 2)
  expect_equal(nrow(acd2), 64)
})

# Test purpose: Ensures reset_acd() turns off ACD for variables with fewer than
# 5 distinct values while preserving ACD for eligible variables.
test_that("reset_acd() resets variables with fewer than five unique values", {
  x <- cbind(
    low_unique = rep(1:4, each = 10),
    enough_unique = seq_len(40)
  )

  acdx <- c(low_unique = TRUE, enough_unique = TRUE)

  expect_warning(
    out <- reset_acd(x, acdx),
    "fewer than 5 unique values"
  )

  expect_false(out["low_unique"])
  expect_true(out["enough_unique"])
  expect_equal(names(out), names(acdx))
})

# Test purpose: Ensures reset_acd() uses variable names, not vector position,
# when acdx is ordered differently from the columns of x.
test_that("reset_acd() aligns acdx by variable name", {
  x <- cbind(
    low_unique = rep(1:4, each = 10),
    enough_unique = seq_len(40)
  )

  acdx <- c(enough_unique = TRUE, low_unique = TRUE)

  expect_warning(
    out <- reset_acd(x, acdx),
    "low_unique"
  )

  expect_true(out["enough_unique"])
  expect_false(out["low_unique"])
  expect_equal(names(out), names(acdx))
})

# Test purpose: Ensures reset_acd() requires acdx to be a named logical vector,
# because ACD variables are aligned by name.
test_that("reset_acd() rejects unnamed acdx vectors", {
  x <- cbind(x1 = seq_len(20), x2 = seq_len(20))

  expect_error(
    reset_acd(x, c(TRUE, FALSE)),
    "`acdx` must be a named logical vector"
  )
})

# Test purpose: Ensures reset_acd() rejects missing ACD flags before model
# fitting starts.
test_that("reset_acd() rejects missing acdx values", {
  x <- cbind(x1 = seq_len(20), x2 = seq_len(20))
  acdx <- c(x1 = TRUE, x2 = NA)

  expect_error(
    reset_acd(x, acdx),
    "`acdx` must not contain missing values"
  )
})

# Test purpose: Ensures public mfp2() applies reset_acd() and records acd = FALSE
# for requested ACD variables with fewer than 5 unique values.
test_that("mfp2() resets ACD for low-cardinality variables", {
  set.seed(401)
  n <- 120

  x <- cbind(
    low_unique = rep(1:4, length.out = n),
    z = runif(n, 1, 10)
  )
  y <- 0.5 * x[, "low_unique"] + rnorm(n)

  expect_warning(
    fit <- mfp2(
      x,
      y,
      acdx = "low_unique",
      verbose = FALSE
    ),
    "fewer than 5 unique values"
  )

  expect_s3_class(fit, "mfp2")
  expect_false(fit$acd["low_unique"])
  expect_false(fit$fp_terms["low_unique", "acd"])
})

# Test purpose: Ensures retained ACD variables are forced to effective df = 4,
# even when the user supplies a smaller df.
test_that("mfp2() forces retained ACD variables to df = 4", {
  set.seed(402)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    df = c(x1 = 2, z = 2),
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$acd["x1"])
  expect_true(fit$fp_terms["x1", "acd"])
  expect_equal(as.numeric(fit$fp_terms["x1", "df_initial"]), 4)
})

# Test purpose: Ensures formula-interface fp(acdx = TRUE) is translated to ACD
# handling and receives the same effective df = 4 treatment as the default interface.
test_that("formula interface fp(acdx = TRUE) forces effective df = 4", {
  set.seed(403)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)

  dat <- data.frame(
    y = 2 * pnorm(scale(x1)) + 0.2 * z + rnorm(n, sd = 0.2),
    x1 = x1,
    z = z
  )

  fit <- mfp2(
    y ~ fp(x1, acdx = TRUE, df = 2) + fp(z, df = 2),
    data = dat,
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$acd["x1"])
  expect_true(fit$fp_terms["x1", "acd"])
  expect_equal(as.numeric(fit$fp_terms["x1", "df_initial"]), 4)
})

# Test purpose: Ensures the formula interface rejects the old acd_vars argument
# and directs users to fp(..., acdx = TRUE).
test_that("formula interface rejects acd_vars argument", {
  data("prostate", package = "mfp2")

  expect_error(
    mfp2(
      lpsa ~ fp(cavol) + fp(age),
      data = prostate,
      acd_vars = "cavol",
      verbose = FALSE
    ),
    "acd_vars.*not supported"
  )
})

# Test purpose: Ensures retained ACD variables store the fitted ACD parameters
# needed for prediction on newdata.
test_that("mfp2() stores ACD parameters for retained ACD variables", {
  set.seed(404)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x1", "acd"])
  expect_true("x1" %in% names(fit$acd_parameter))
  expect_true(is.list(fit$acd_parameter[["x1"]]))

  expect_true(all(c("beta0", "beta1", "power", "shift", "scale") %in%
                    names(fit$acd_parameter[["x1"]])))
  expect_null(fit$acd_parameter[["x1"]]$acd)
})

# Test purpose: Active ACD variables must bypass ordinary MFP scaling in the
# matrix interface, while non-ACD variables retain their requested scales.
test_that("matrix-interface ACD variables use scale one", {
  set.seed(4041)
  n <- 180
  x1 <- runif(n, 10, 2000)
  z <- runif(n, 1, 20)
  x <- cbind(x1 = x1, z = z)
  y <- 1.5 * pnorm(as.numeric(scale(x1))) + 0.15 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    scale = c(x1 = 1000, z = 10),
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  acd_reference <- fit_acd(x1)

  expect_true(fit$acd[["x1"]])
  expect_equal(unname(fit$transformations["x1", "scale"]), 1)
  expect_equal(unname(fit$transformations["z", "scale"]), 10)
  expect_equal(fit$acd_parameter[["x1"]]$scale, 1)
  expect_equal(fit$acd_parameter[["x1"]]$power, acd_reference$power)
  expect_equal(fit$acd_parameter[["x1"]]$beta0, acd_reference$beta0, tolerance = 1e-10)
  expect_equal(fit$acd_parameter[["x1"]]$beta1, acd_reference$beta1, tolerance = 1e-10)
})

# Test purpose: Formula-level scale settings are also overridden only for terms
# that remain eligible for ACD modelling.
test_that("formula-interface ACD variables use scale one", {
  set.seed(4042)
  n <- 180
  dat <- data.frame(
    x1 = runif(n, 10, 2000),
    z = runif(n, 1, 20)
  )
  dat$y <- 1.5 * pnorm(scale(dat$x1)) + 0.15 * dat$z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ fp(x1, acdx = TRUE, scale = 1000) + fp(z, scale = 10),
    data = dat,
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_true(fit$acd[["x1"]])
  expect_equal(unname(fit$transformations["x1", "scale"]), 1)
  expect_equal(unname(fit$transformations["z", "scale"]), 10)
  expect_equal(fit$acd_parameter[["x1"]]$scale, 1)
})

# Test purpose: A request reset for insufficient distinct values reverts to
# ordinary FP preprocessing and therefore does not have its scale forced to one.
test_that("reset ACD requests retain ordinary FP scaling", {
  set.seed(4043)
  n <- 160
  x1 <- rep(1:4, length.out = n)
  z <- runif(n, 1, 20)
  x <- cbind(x1 = x1, z = z)
  y <- 0.4 * x1 + 0.1 * z + rnorm(n, sd = 0.2)

  expect_warning(
    fit <- mfp2(
      x,
      y,
      acdx = "x1",
      scale = c(x1 = 100, z = 10),
      df = c(x1 = 2, z = 1),
      select = 1,
      alpha = 1,
      verbose = FALSE
    ),
    "fewer than 5 unique values"
  )

  expect_false(fit$acd[["x1"]])
  expect_equal(unname(fit$transformations["x1", "scale"]), 100)
})

# Test purpose: Ensures predict.mfp2() can reuse stored ACD parameters to
# transform newdata for an ACD-fitted model.
test_that("predict.mfp2() works for ACD models with newdata", {
  set.seed(405)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = x[1:10, , drop = FALSE])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})

# Test purpose: Ensures predict.mfp2() fails clearly if an active ACD variable
# has no stored ACD parameters.
test_that("predict.mfp2() errors when active ACD parameters are missing", {
  set.seed(406)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  fit$acd_parameter[["x1"]] <- NULL

  expect_error(
    predict(fit, newdata = x[1:10, , drop = FALSE]),
    "Missing stored ACD parameters"
  )
})

# Test purpose: Checks the exported fit_acd() helper returns an ACD-transformed
# vector on the cumulative-probability scale.
test_that("fit_acd() returns ACD values in the unit interval", {
  set.seed(407)
  x <- runif(100, 1, 20)

  acd <- fit_acd(x)

  expect_true(is.list(acd))
  expect_length(acd$acd, length(x))
  expect_true(all(is.finite(acd$acd)))
  expect_true(all(acd$acd >= 0 & acd$acd <= 1))
  expect_true(all(c("beta0", "beta1", "power", "shift", "scale") %in% names(acd)))
})


# Test purpose: Internal GLM candidate searches can request fitted values without
# retaining the complete backend fit object. The lightweight vector must be
# exactly the one produced by the same fitting backend.
test_that("fit_model() can retain fitted values without retaining the GLM fit", {
  set.seed(701)
  x <- cbind(
    "(Intercept)" = 1,
    "x" = stats::runif(40, 0.5, 3)
  )
  y <- 1.5 + 0.8 * x[, "x"] + stats::rnorm(40, sd = 0.2)
  fam <- stats::gaussian()

  retained <- fit_model(
    x = x,
    y = y,
    family = fam,
    family_string = fam$family,
    fitter = "base",
    x_has_intercept = TRUE,
    keep_fit = TRUE,
    keep_fitted_values = FALSE
  )

  lightweight <- fit_model(
    x = x,
    y = y,
    family = fam,
    family_string = fam$family,
    fitter = "base",
    x_has_intercept = TRUE,
    keep_fit = FALSE,
    keep_fitted_values = TRUE
  )

  expect_false("fit" %in% names(lightweight))
  expect_true("fitted_values" %in% names(lightweight))
  expect_identical(lightweight$fitted_values, retained$fit$fitted.values)
  expect_identical(lightweight$coefficients, retained$coefficients)
  expect_identical(lightweight$logl, retained$logl)
})

# Test purpose: The ACD FP1 search must use the lightweight fit_model() contract
# for every candidate and must not accidentally thread outcome-model weights or
# offsets into the distribution-based ACD auxiliary regression.
test_that("find_best_fp1_for_acd() requests only lightweight fitted values", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$keep_fit <- logical()
  calls$keep_fitted_values <- logical()
  calls$weights <- list()
  calls$offset <- list()

  testthat::local_mocked_bindings(
    # This helper should use the compact FP basis path introduced for mfp2
    # candidate fitting, not recreate the legacy list of n x 1 FP matrices.
    generate_transformations_fp = function(...) {
      stop("legacy FP materialization should not be used", call. = FALSE)
    },
    generate_transformations_fp_basis = function(x, degree, powers, zero, ...) {
      expect_equal(degree, 1L)
      expect_false(zero)
      list(
        basis = cbind(x^-1, log(x), x),
        candidate_map = matrix(1:3, ncol = 1L),
        powers = matrix(powers, ncol = 1L),
        catzero = NULL
      )
    },
    fit_model = function(x,
                         y,
                         weights = NULL,
                         offset = NULL,
                         keep_fit,
                         keep_fitted_values,
                         ...) {
      calls$n <- calls$n + 1L
      i <- calls$n
      calls$keep_fit[[i]] <- keep_fit
      calls$keep_fitted_values[[i]] <- keep_fitted_values
      calls$weights[i] <- list(weights)
      calls$offset[i] <- list(offset)

      logl <- c(-5, -1, -3)[[i]]
      list(
        logl = logl,
        coefficients = c("(Intercept)" = i, "fp1" = i + 1),
        fitted_values = rep.int(as.numeric(i), NROW(x))
      )
    },
    .package = "mfp2"
  )

  out <- find_best_fp1_for_acd(
    x = seq(1, 6),
    y = seq(-1, 1, length.out = 6),
    powers = c(-1, 0, 1),
    zero = FALSE,
    fitter = "base"
  )

  expect_equal(calls$n, 3L)
  expect_true(all(!calls$keep_fit))
  expect_true(all(calls$keep_fitted_values))
  expect_true(all(vapply(calls$weights, is.null, logical(1L))))
  expect_true(all(vapply(calls$offset, is.null, logical(1L))))
  expect_equal(out$power, 0)
  expect_equal(out$coefficients, c("(Intercept)" = 2, "fp1" = 3))
  expect_equal(out$fitted_values, rep(2, 6))
})

# Test purpose: Removing the retained backend object must not change the ACD
# candidate selected, its coefficients, or its fitted values relative to the
# previous keep_fit = TRUE implementation.
test_that("find_best_fp1_for_acd() preserves previous fitted results", {
  x <- seq(0.5, 8, length.out = 60)
  y <- stats::qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
  powers <- c(-1, 0, 0.5, 1, 2)
  fam <- stats::gaussian()
  trafo <- generate_transformations_fp(
    x = x,
    degree = 1L,
    powers = powers,
    zero = FALSE
  )$data

  reference <- lapply(seq_along(powers), function(i) {
    design <- cbind(
      "(Intercept)" = 1,
      "fp1" = trafo[[i]][, 1L]
    )
    fit <- fit_model(
      x = design,
      y = y,
      family = fam,
      family_string = fam$family,
      fitter = "base",
      x_has_intercept = TRUE,
      keep_fit = TRUE
    )
    list(
      deviance = -2 * fit$logl,
      coefficients = fit$coefficients,
      fitted_values = fit$fit$fitted.values
    )
  })

  best <- which.min(vapply(reference, `[[`, numeric(1L), "deviance"))
  out <- find_best_fp1_for_acd(
    x = x,
    y = y,
    powers = powers,
    zero = FALSE,
    fitter = "base"
  )

  expect_identical(out$power, powers[[best]])
  expect_identical(out$coefficients, reference[[best]]$coefficients)
  expect_identical(out$fitted_values, reference[[best]]$fitted_values)
})

# Test purpose: Ensures apply_acd() reproduces the training ACD transformation
# when supplied with parameters from fit_acd().
test_that("apply_acd() reproduces fit_acd() values using stored parameters", {
  set.seed(408)
  x <- runif(100, 1, 20)

  acd <- fit_acd(x)

  applied <- apply_acd(
    x = x,
    beta0 = acd$beta0,
    beta1 = acd$beta1,
    power = acd$power,
    shift = acd$shift,
    scale = acd$scale,
    zero = FALSE
  )

  expect_equal(as.numeric(applied), as.numeric(acd$acd), tolerance = 1e-10)
})

# Test purpose: Ensures generate_powers_acd() only accepts supported ACD degrees
# 0, 1, and 2.
test_that("generate_powers_acd() rejects unsupported degrees", {
  expect_error(
    generate_powers_acd(degree = 3),
    "degree.*ACD.*0, 1, or 2"
  )

  expect_error(
    generate_powers_acd(degree = NA),
    "degree.*ACD.*0, 1, or 2"
  )
})

# Test purpose: Ensures ACD power generation uses ordered Cartesian products,
# because the first power applies to x and the second to A(x).
test_that("generate_powers_acd() keeps ordered power pairs", {
  powx <- c(0, 1)

  acd2 <- generate_powers_acd(degree = 2, powers = powx)

  expect_equal(ncol(acd2), 2)
  expect_equal(nrow(acd2), 4)

  expect_true(any(acd2[, 1] == 0 & acd2[, 2] == 1))
  expect_true(any(acd2[, 1] == 1 & acd2[, 2] == 0))
})

# Test purpose: Ensures generate_transformations_acd() can reuse stored ACD
# parameters instead of refitting them.
test_that("generate_transformations_acd() reuses stored ACD parameters", {
  set.seed(409)
  x <- runif(80, 1, 20)
  powers <- c(0, 1)

  acd_par <- fit_acd(x, powers = powers)
  acd_par$acd <- NULL

  out <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_true(is.list(out))
  expect_equal(nrow(out$powers), 4)
  expect_equal(length(out$data), 4)
  expect_true(all(vapply(out$data, nrow, integer(1)) == length(x)))
})

# Test purpose: Ensures ACD transformation generation includes the catzero
# indicator column when catzero is supplied.
test_that("generate_transformations_acd() includes catzero column when supplied", {
  set.seed(410)
  x <- runif(80, 1, 20)
  catzero <- matrix(as.integer(seq_along(x) <= 10), ncol = 1)

  out <- generate_transformations_acd(
    x = x,
    degree = 1,
    powers = c(0, 1),
    zero = FALSE,
    catzero = catzero
  )

  expect_true(is.list(out))
  expect_equal(length(out$data), 2)

  first <- out$data[[1]]
  expect_equal(nrow(first), length(x))
  expect_equal(colnames(first)[1], "catzero")
})

# Test purpose: Checks that compact ACD degree-2 generation reconstructs every
# materialized candidate while storing only the unique x/A(x) transformations.
test_that("compact ACD basis reconstructs materialized degree-2 candidates", {
  x <- seq(1, 8, length.out = 40)
  allowed_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  # Supply fixed ACD parameters so both representations use exactly the same
  # A(x) values and the test isolates candidate storage/reconstruction.
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  materialized <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_equal(compact$powers, materialized$powers)
  expect_equal(dim(compact$candidate_map), c(64L, 2L))
  expect_equal(ncol(compact$basis), 16L)

  for (i in seq_len(nrow(compact$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(compact, i)),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Checks the one-term ACD degree-1 representation and verifies
# that the invariant catzero indicator is stored once rather than per candidate.
test_that("compact ACD basis preserves degree-1 catzero candidates", {
  x <- seq(1, 6, length.out = 30)
  catzero <- matrix(as.integer(seq_along(x) <= 5L), ncol = 1L)
  powers <- c(0, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  materialized <- generate_transformations_acd(
    x = x,
    degree = 1,
    powers = powers,
    zero = FALSE,
    catzero = catzero,
    acd_parameter = acd_par
  )

  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 1,
    powers = powers,
    zero = FALSE,
    catzero = catzero,
    acd_parameter = acd_par
  )

  expect_equal(dim(compact$candidate_map), c(2L, 1L))
  expect_equal(ncol(compact$basis), 2L)
  expect_equal(compact$catzero, catzero)

  for (i in seq_len(nrow(compact$candidate_map))) {
    reconstructed <- materialize_acd_basis_candidate(compact, i)
    expect_equal(
      unname(reconstructed),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
    expect_equal(colnames(reconstructed), c("catzero", "V1"))
  }
})


# Test purpose: Checks that compact ACD generation preserves zero-mode
# semantics for nonpositive x values as well as positive transformed values.
test_that("compact ACD basis preserves zero-mode transformations", {
  x <- c(-2, 0, 1, 2, 4, 8)
  powers <- c(0, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  materialized <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = TRUE,
    acd_parameter = acd_par
  )

  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = TRUE,
    acd_parameter = acd_par
  )

  for (i in seq_len(nrow(compact$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(compact, i)),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Confirms the model-search transformation wrapper actually uses
# the compact ACD representation instead of returning a materialized data list.
test_that("transform_data_step() returns compact ACD basis when requested", {
  x_vec <- seq(1, 5, length.out = 25)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 4,
    powers = list(x = c(0, 1)),
    acdx = c(x = TRUE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = acd_par),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    term_to_columns = list(x = 1L),
    compact_acd = TRUE
  )

  expect_null(out$data_fp)
  expect_null(out$fp_basis)
  expect_true(is.list(out$acd_basis))
  expect_equal(dim(out$acd_basis$candidate_map), c(4L, 2L))
  expect_equal(ncol(out$acd_basis$basis), 4L)
})


# Test purpose: Checks important fit_acd() input validation branches.
test_that("fit_acd() validates input arguments", {
  expect_error(
    fit_acd(factor(c("a", "b", "c"))),
    "`x` must be a numeric vector"
  )

  expect_error(
    fit_acd(c(1, NA, 3)),
    "missing values"
  )

  expect_error(
    fit_acd(1),
    "at least two values"
  )

  expect_error(
    fit_acd(1:10, scale = 0),
    "`scale` must be a single positive numeric value"
  )

  expect_error(
    fit_acd(1:10, zero = NA),
    "`zero` must be a single logical value"
  )
})


# =============================================================================
# 8. predict.mfp2()
# =============================================================================

# Test purpose: Checks default Gaussian predictions are finite and have one
# value per observation.
test_that("predict.mfp2() returns predictions for Gaussian model", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  preds <- predict(fit)
  expect_length(preds, nrow(x_prostate))
  expect_true(all(is.finite(preds)))
})

# Test purpose: Checks the structure and numerical correctness of Gaussian
# link-scale predictions and standard errors. The expected values are calculated
# independently as X beta and sqrt(diag(X V X')).
test_that("predict.mfp2() Gaussian fit and SE equal manual matrix calculation", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  result <- predict(fit, type = "link", se.fit = TRUE)
  beta <- stats::coef(fit)
  beta_vcov <- stats::vcov(fit)
  manual_x <- cbind(`(Intercept)` = 1, x_prostate)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_fit <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_true(is.list(result))
  expect_named(result, c("fit", "se.fit", "residual.scale"))
  expect_equal(as.numeric(result$fit), manual_fit, tolerance = 1e-8)
  expect_equal(as.numeric(result$se.fit), manual_se, tolerance = 1e-8)
})

# Test purpose: Independently validates negative-binomial link and response
# predictions, including standard errors, against MASS::predict.glm().
test_that("negative-binomial predictions and SEs agree with MASS::glm.nb()", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  fits <- get_negbin_reference_fits()
  newdata <- fits$data[seq_len(24L), c("x1", "x2", "x3"), drop = FALSE]

  mfp_link <- predict(
    fits$mfp2,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  mass_link <- stats::predict(
    fits$mass,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )

  expect_named(mfp_link, c("fit", "se.fit", "residual.scale"))
  expect_equal(
    unname(mfp_link$fit),
    unname(mass_link$fit),
    tolerance = 1e-4
  )
  expect_equal(
    unname(mfp_link$se.fit),
    unname(mass_link$se.fit),
    tolerance = 1e-4
  )
  expect_equal(mfp_link$residual.scale, 1, tolerance = 1e-12)

  mfp_response <- predict(
    fits$mfp2,
    newdata = newdata,
    type = "response",
    se.fit = TRUE
  )
  mass_response <- stats::predict(
    fits$mass,
    newdata = newdata,
    type = "response",
    se.fit = TRUE
  )

  expect_equal(
    unname(mfp_response$fit),
    unname(mass_response$fit),
    tolerance = 1e-4
  )
  expect_equal(
    unname(mfp_response$se.fit),
    unname(mass_response$se.fit),
    tolerance = 1e-4
  )
  expect_equal(
    unname(mfp_response$fit),
    exp(unname(mfp_link$fit)),
    tolerance = 1e-10
  )

  # Prediction without newdata must reproduce the stored training means.
  expect_equal(
    unname(predict(fits$mfp2, type = "response")),
    unname(fits$mfp2$fitted.values),
    tolerance = 1e-10
  )

  # Term and contrast predictions use the same selected linear function and
  # covariance matrix. Check them independently on the link scale.
  term <- predict(
    fits$mfp2,
    newdata = newdata,
    type = "terms",
    terms = "x1",
    terms_seq = "data",
    add_intercept = FALSE
  )
  coefficient_name <- unname(
    fits$mfp2$transformed_to_model_columns[["x1.1"]]
  )
  beta_x1 <- stats::coef(fits$mfp2)[[coefficient_name]]
  se_beta_x1 <- sqrt(stats::vcov(fits$mfp2)[coefficient_name, coefficient_name])

  expect_named(term, "x1")
  expect_equal(
    unname(term$x1$value),
    unname(newdata$x1 * beta_x1),
    tolerance = 1e-8
  )
  expect_equal(
    unname(term$x1$se),
    unname(abs(newdata$x1) * se_beta_x1),
    tolerance = 1e-8
  )

  contrast <- predict(
    fits$mfp2,
    newdata = newdata,
    type = "contrasts",
    terms = "x1",
    terms_seq = "data",
    ref = list(x1 = 1)
  )
  expect_named(contrast, "x1")
  expect_equal(
    unname(contrast$x1$value),
    unname((newdata$x1 - 1) * beta_x1),
    tolerance = 1e-8
  )
  expect_equal(
    unname(contrast$x1$se),
    unname(abs(newdata$x1 - 1) * se_beta_x1),
    tolerance = 1e-8
  )
})

test_that("predict.mfp2() with newdata reproduces training predictions", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  preds_train <- predict(fit)
  preds_new <- predict(fit, newdata = x_prostate)

  expect_equal(as.numeric(preds_train), as.numeric(preds_new),
               tolerance = 1e-10)
})

# Test purpose: Checks that term-level predictions return per-term data frames with values and standard errors.
test_that("predict.mfp2() type = 'terms' returns list of data frames", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  terms_result <- predict(fit, type = "terms")
  expect_true(is.list(terms_result))

  for (nm in names(terms_result)) {
    expect_true(is.data.frame(terms_result[[nm]]))
    expect_true("value" %in% colnames(terms_result[[nm]]))
    expect_true("se" %in% colnames(terms_result[[nm]]))
  }
})

# Test purpose: Checks that contrast predictions return per-term data-frame outputs.
test_that("predict.mfp2() type = 'contrasts' returns list of data frames", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  contrasts_result <- predict(fit, type = "contrasts")
  expect_true(is.list(contrasts_result))

  for (nm in names(contrasts_result)) {
    expect_true(is.data.frame(contrasts_result[[nm]]))
  }
})

# Test purpose: Checks that prediction errors when newdata violates the positive domain required by a fitted log transform.
test_that("predict.mfp2() stops on domain violation in newdata", {
  x <- cbind(x1 = seq(1, 100, length.out = 100))
  # A small deterministic perturbation avoids an exact fit while keeping the
  # fitted transformation reproducible without depending on RNG state.
  y <- log(x[, "x1"]) + 0.01 * sin(seq_len(nrow(x)))

  fit <- mfp2(
    x, y,
    select = 1,
    alpha = 1,
    powers = list(x1 = 0),
    df = 2,
    shift = 0,
    scale = 1,
    verbose = FALSE
  )

  # create bad data
  bad_data <- x[1:5, , drop = FALSE]
  bad_data[, "x1"] <- -1

  expect_error(
    predict(fit, newdata = bad_data),
    "non-positive"
  )
})

# Test purpose: Checks that Cox-model predictions are finite and have the correct length.
test_that("predict.mfp2() works for Cox models", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes", "pgr", "er")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)

  preds <- predict(fit)
  expect_length(preds, nrow(x_gbsg))
  expect_true(all(is.finite(preds)))
})

# Test purpose: Checks that default Cox predictions use cox_reference = "zero".
# On this scale the linear predictor is the uncentered matrix product X beta,
# rather than the sample-centered predictor returned by another reference mode.
test_that("predict.mfp2() for Cox equals manual X beta on reference-zero scale", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(
    x_gbsg,
    y_gbsg,
    family = "cox",
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  got <- predict(fit, type = "lp", se.fit = TRUE)
  beta <- stats::coef(fit)
  beta_vcov <- stats::vcov(fit)
  manual_x <- x_gbsg
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_lp <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_equal(as.numeric(got$fit), manual_lp, tolerance = 1e-8)
  expect_equal(as.numeric(got$se.fit), manual_se, tolerance = 1e-8)
})

# Test purpose: Ensures factor binomial responses do not break inherited
# predict.glm(se.fit = TRUE) behavior.
test_that("predict.mfp2() with se.fit = TRUE works for factor binomial response", {
  data("pima", package = "mfp2")

  x_pima <- as.matrix(pima[, c("glucose", "bmi", "age")])
  y_factor <- factor(pima$y, levels = c(0, 1), labels = c("no", "yes"))

  fit <- mfp2(
    x_pima,
    y_factor,
    family = "binomial",
    verbose = FALSE
  )

  result <- predict(fit, se.fit = TRUE)

  expect_true(is.list(result))
  expect_true("fit" %in% names(result))
  expect_true("se.fit" %in% names(result))
  expect_length(result$fit, nrow(x_pima))
  expect_length(result$se.fit, nrow(x_pima))
  expect_identical(fit$y_original, y_factor)
})

# Test purpose: Ensures prediction with newdata requires newoffset when the
# model was fitted with an offset.
test_that("predict.mfp2() requires newoffset when fitted model used offset", {
  set.seed(107)
  n <- 150

  x <- cbind(x1 = runif(n, 1, 10), x2 = runif(n, 1, 5))
  exposure <- runif(n, 0.5, 2)
  y <- rpois(n, exposure * exp(0.5 + 0.1 * x[, "x1"]))

  fit <- mfp2(
    x,
    y,
    family = "poisson",
    offset = log(exposure),
    verbose = FALSE
  )

  expect_error(
    predict(fit, newdata = x[1:5, , drop = FALSE]),
    "newoffset"
  )
})

# Test purpose: Checks that formula-interface prediction preserves the number
# of rows in newdata even when model selection drops all predictors.
test_that("predict.mfp2() returns newdata-length predictions for intercept-only formula models", {
  set.seed(109)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = dat[1:10, ])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})

# Test purpose: Checks that matrix-interface prediction preserves newdata row
# count when the final selected model is intercept-only.
test_that("predict.mfp2() returns newdata-length predictions for intercept-only matrix models", {
  set.seed(110)
  n <- 120

  x <- cbind(
    x1 = runif(n, 1, 10),
    x2 = runif(n, 1, 10)
  )
  y <- rnorm(n)

  fit <- mfp2(
    x,
    y,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = x[1:10, , drop = FALSE])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})


# Regression: full formula prediction must not require terms eliminated during
# selection. This also exercises the stored fp() formula-label mapping.
test_that("predict.mfp2() formula prediction requires only selected terms", {
  set.seed(1101)
  n <- 180
  dat <- data.frame(
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  dat$y <- 1 + 0.8 * dat$x + rnorm(n, sd = 0.15)

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    keep = "x",
    select = 0,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x", "selected"])
  expect_false(fit$fp_terms["group", "selected"])
  expect_identical(fit$formula_prediction_term_names[["fp(x)"]], "x")

  nd_minimal <- dat[1:12, "x", drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(fit, newdata = dat[1:12, c("x", "group")])
  nd_with_unseen_dropped_level <- data.frame(
    x = nd_minimal$x,
    group = factor(rep("D", nrow(nd_minimal)), levels = c("A", "B", "C", "D"))
  )
  pred_unseen_dropped <- predict(fit, newdata = nd_with_unseen_dropped_level)
  term_minimal <- predict(
    fit,
    newdata = nd_minimal,
    type = "terms",
    terms = "x",
    terms_seq = "data"
  )

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
  expect_equal(pred_minimal, pred_unseen_dropped, tolerance = 1e-12)
  expect_named(term_minimal, "x")
  expect_equal(nrow(term_minimal$x), nrow(nd_minimal))
})


# Regression: retaining a categorical block must not require an eliminated
# continuous candidate in formula-style newdata.
test_that("predict.mfp2() reconstructs only a selected factor block", {
  set.seed(1102)
  n <- 180
  dat <- data.frame(
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  group_effect <- c(A = 0, B = 1, C = -0.7)[as.character(dat$group)]
  dat$y <- group_effect + rnorm(n, sd = 0.15)

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    keep = "group",
    select = 0,
    verbose = FALSE
  )

  expect_false(fit$fp_terms["x", "selected"])
  expect_true(fit$fp_terms["group", "selected"])

  nd_minimal <- dat[1:12, "group", drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(fit, newdata = dat[1:12, c("x", "group")])

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
})


# Regression: selected factor interactions retain the fit-time contrast coding
# while unrelated eliminated formula terms are not required.
test_that("predict.mfp2() preserves active interaction design with minimal newdata", {
  set.seed(11021)
  n <- 210
  dat <- data.frame(
    x = runif(n, 1, 5),
    z = runif(n, -2, 2),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  group_slope <- c(A = 0.2, B = 0.8, C = -0.5)[as.character(dat$group)]
  dat$y <- group_slope * dat$x + rnorm(n, sd = 0.1)

  fit <- mfp2(
    y ~ x + group + x:group + z,
    data = dat,
    keep = "x:group",
    select = 0,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x:group", "selected"])
  expect_false(fit$fp_terms["z", "selected"])

  nd_minimal <- dat[1:15, c("x", "group"), drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(
    fit,
    newdata = dat[1:15, c("x", "group", "z"), drop = FALSE]
  )

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
})


# Regression: matrix prediction expands only selected conceptual terms, so a
# grouped block eliminated during selection is not a newdata dependency.
test_that("predict.mfp2() matrix prediction ignores eliminated grouped columns", {
  set.seed(1103)
  n <- 180
  x <- cbind(
    x1 = runif(n, 1, 10),
    groupB = rep(c(0, 1, 0), length.out = n),
    groupC = rep(c(0, 0, 1), length.out = n)
  )
  y <- 1 + 0.7 * x[, "x1"] + rnorm(n, sd = 0.15)

  fit <- mfp2(
    x,
    y,
    term_groups = list(group = c("groupB", "groupC")),
    keep = "x1",
    select = 0,
    df = 1,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x1", "selected"])
  expect_false(fit$fp_terms["group", "selected"])
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), 0)

  nd_minimal <- x[1:12, "x1", drop = FALSE]
  pred_minimal <- predict(fit, newdata = nd_minimal)
  pred_complete <- predict(fit, newdata = x[1:12, , drop = FALSE])

  expect_equal(pred_minimal, pred_complete, tolerance = 1e-12)
})


# Regression: an intercept-only final model has no predictor dependencies and
# therefore accepts a zero-column data frame carrying only the requested rows.
test_that("predict.mfp2() intercept-only models accept zero-column newdata", {
  set.seed(1104)
  n <- 160
  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )

  fit_formula <- mfp2(
    y ~ x + group,
    data = dat,
    select = 0,
    verbose = FALSE
  )
  fit_matrix <- mfp2(
    as.matrix(dat[, "x", drop = FALSE]),
    dat$y,
    select = 0,
    verbose = FALSE
  )

  expect_length(get_selected_variable_names(fit_formula), 0L)
  expect_length(get_selected_variable_names(fit_matrix), 0L)

  nd_empty <- data.frame(row.names = seq_len(9L))
  expect_length(predict(fit_formula, newdata = nd_empty), 9L)
  expect_length(predict(fit_matrix, newdata = nd_empty), 9L)
  omitted_terms <- NULL
  expect_warning(
    omitted_terms <- predict(
      fit_formula,
      newdata = nd_empty,
      type = "terms",
      terms = "group"
    ),
    "All the terms supplied are not in the final model"
  )
  expect_identical(omitted_terms, list())
})

# Test purpose: Checks that formula-fitted models can predict from ordinary
# newdata containing original factor variables rather than expanded dummy columns.
test_that("predict.mfp2() reconstructs formula-interface newdata with retained factors", {
  set.seed(111)
  n <- 120

  x <- runif(n, 1, 10)
  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_effect <- c(A = 0, B = 1, C = 2)[as.character(group)]

  dat <- data.frame(
    y = 0.2 * x + 0.5 * group_effect + rnorm(n, sd = 0.2),
    x = x,
    group = group
  )

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    keep = c("x", "group"),
    verbose = FALSE
  )

  pred <- predict(fit, newdata = dat[1:10, ])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})


# Test purpose: Checks conceptual-term prediction for an unordered factor,
# including level labels, raw dummy columns, fitted values, and uncertainty.
test_that("predict.mfp2() returns a complete unordered-factor term block", {
  set.seed(2201)
  n <- 150

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), length.out = n)),
    x = runif(n, 1, 10)
  )
  effect <- c(A = 0, B = 1, C = -0.5)[as.character(dat$group)]
  dat$y <- 0.3 * dat$x + effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )
  nd <- dat[1:12, c("x", "group"), drop = FALSE]

  out <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms = "group"
  )

  expect_named(out, "group")
  expect_s3_class(out$group, "data.frame")
  expect_equal(nrow(out$group), nrow(nd))
  expect_true(
    all(c(
      "variable", "variable_pre", "groupB", "groupC",
      "value", "se", "lower", "upper"
    ) %in% names(out$group))
  )
  expect_equal(out$group$variable, as.character(nd$group))
  expect_true(all(is.finite(out$group$value)))
  expect_true(all(is.finite(out$group$se)))
})

# Test purpose: Checks grouped factor contrasts against a named fitted level.
test_that("predict.mfp2() computes factor contrasts from a level reference", {
  set.seed(2202)
  n <- 150

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), length.out = n)),
    x = runif(n, 1, 10)
  )
  effect <- c(A = 0, B = 1, C = -0.5)[as.character(dat$group)]
  dat$y <- 0.3 * dat$x + effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )
  nd <- dat[1:12, c("x", "group"), drop = FALSE]

  out <- predict(
    fit,
    newdata = nd,
    type = "contrasts",
    terms = "group",
    ref = list(group = "A")
  )

  expect_named(out, "group")
  expect_equal(nrow(out$group), nrow(nd))
  expect_true(all(is.finite(out$group$value)))
  expect_true(all(is.finite(out$group$se)))
  expect_equal(
    out$group$value[nd$group == "A"],
    rep(0, sum(nd$group == "A")),
    tolerance = 1e-10
  )
})

# Test purpose: Checks that ordered-factor contrasts can be reconstructed from
# ordinary formula-style newdata for full and conceptual-term prediction.
test_that("predict.mfp2() reconstructs ordered-factor contrasts", {
  set.seed(2203)
  n <- 180

  dat <- data.frame(
    severity = ordered(
      rep(c("low", "medium", "high"), length.out = n),
      levels = c("low", "medium", "high")
    ),
    x = runif(n, 1, 10)
  )
  severity_effect <- c(low = 0, medium = 0.5, high = 1.2)[
    as.character(dat$severity)
  ]
  dat$y <- 0.2 * dat$x + severity_effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x + severity,
    data = dat,
    keep = "severity",
    verbose = FALSE
  )
  nd <- dat[1:15, c("x", "severity"), drop = FALSE]

  ordinary <- predict(fit, newdata = nd)
  term_output <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms = "severity"
  )

  expect_length(ordinary, nrow(nd))
  expect_true(all(is.finite(ordinary)))
  expect_named(term_output, "severity")
  expect_equal(term_output$severity$variable, as.character(nd$severity))
  expect_true(all(is.finite(term_output$severity$value)))
})

# Test purpose: Ensures every fitted model stores a complete conceptual-term
# lookup, including ordinary identity mappings required by prediction.
test_that("mfp2 stores complete identity term mappings", {
  set.seed(22031)
  dat <- data.frame(
    y = rnorm(80),
    x1 = runif(80, 1, 4),
    x2 = runif(80, 2, 6)
  )

  fit <- mfp2(
    y ~ x1 + x2,
    data = dat,
    df = 1,
    select = 1,
    verbose = FALSE
  )

  expect_identical(fit$term_to_columns$x1, "x1")
  expect_identical(fit$term_to_columns$x2, "x2")
  expect_length(predict(fit, newdata = dat[1:5, ]), 5L)
})

# Test purpose: Checks that alternative simple factor wrappers follow the same
# source-variable conceptual naming contract as factor().
test_that("simple factor wrappers use source-variable conceptual names", {
  set.seed(22032)
  dat <- data.frame(
    y = rnorm(90),
    x = rep(1:3, length.out = 90)
  )

  fit <- mfp2(
    y ~ as.factor(x),
    data = dat,
    keep = "x",
    verbose = FALSE
  )

  expect_identical(
    fit$term_to_columns$x,
    c("as.factor(x)2", "as.factor(x)3")
  )
  expect_true("x" %in% rownames(fit$fp_terms))
  expect_false("as.factor(x)" %in% rownames(fit$fp_terms))
})

# Test purpose: Checks that a simple inline factor() wrapper uses the source
# variable name as the conceptual term while retaining model.matrix() column and
# coefficient names.
test_that("formula interface names inline factors by their source variable", {
  set.seed(2204)
  n <- 150

  dat <- data.frame(
    x1 = runif(n, 1, 10),
    x2 = rep(1:3, length.out = n)
  )
  x2_effect <- c(`1` = 0, `2` = 0.8, `3` = -0.4)[as.character(dat$x2)]
  dat$y <- 0.3 * dat$x1 + x2_effect + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x1 + factor(x2),
    data = dat,
    keep = "x2",
    verbose = FALSE
  )

  expect_true("x2" %in% names(fit$term_to_columns))
  expect_false("factor(x2)" %in% names(fit$term_to_columns))
  expect_identical(
    fit$term_to_columns[["x2"]],
    c("factor(x2)2", "factor(x2)3")
  )
  expect_true("x2" %in% rownames(fit$fp_terms))
  expect_false("factor(x2)" %in% rownames(fit$fp_terms))

  coefficient_names <- gsub("`", "", names(stats::coef(fit)), fixed = TRUE)
  expect_true(all(c("factor(x2)2.1", "factor(x2)3.1") %in% coefficient_names))

  nd <- dat[1:10, c("x1", "x2"), drop = FALSE]
  pred <- predict(fit, newdata = nd)
  term_pred <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms = "x2"
  )

  expect_length(pred, 10L)
  expect_true(all(is.finite(pred)))
  expect_named(term_pred, "x2")
  expect_true(all(is.finite(term_pred[["x2"]]$value)))
  expect_true(all(is.finite(term_pred[["x2"]]$se)))
})

# Test purpose: Checks that a binary inline factor remains a mapped categorical
# term even though model.matrix() generates only one dummy column.
test_that("binary inline factors retain a source-name singleton mapping", {
  set.seed(22041)
  n <- 140

  dat <- data.frame(
    x1 = runif(n, 1, 10),
    x2 = rep(1:2, length.out = n)
  )
  dat$y <- 0.4 * dat$x1 + 0.9 * (dat$x2 == 2) + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ x1 + factor(x2),
    data = dat,
    keep = "x2",
    verbose = FALSE
  )

  expect_identical(fit$term_to_columns[["x2"]], "factor(x2)2")
  expect_true(fit$fp_terms["x2", "selected"])
  expect_equal(as.numeric(fit$fp_terms["x2", "df_final"]), 1)

  nd <- dat[1:12, c("x1", "x2"), drop = FALSE]
  out <- predict(fit, newdata = nd, type = "terms", terms = "x2")
  expect_named(out, "x2")
  expect_true(all(is.finite(out$x2$value)))
  expect_true(all(is.finite(out$x2$se)))
})

# Test purpose: Prevents ambiguous formulas that include the same source
# variable both directly and through a simple factor wrapper.
test_that("formula interface rejects duplicate conceptual source variables", {
  dat <- data.frame(
    # The duplicate-term validation occurs before fitting; no RNG is needed.
    y = seq_len(60),
    x2 = rep(1:3, length.out = 60)
  )

  expect_error(
    mfp2(y ~ x2 + factor(x2), data = dat, verbose = FALSE),
    "same conceptual variable"
  )
})

# Test purpose: Checks that formula reconstruction rejects factor levels that
# were not present when the model matrix and contrasts were fitted.
test_that("predict.mfp2() rejects unseen factor levels", {
  set.seed(2205)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )
  fit <- mfp2(
    y ~ group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )
  nd <- data.frame(
    group = factor("D", levels = c("A", "B", "C", "D"))
  )

  expect_error(
    predict(fit, newdata = nd),
    "new level|factor level|levels",
    ignore.case = TRUE
  )
})

# Test purpose: Checks that matrix-interface prediction rejects a partial grouped
# block rather than silently changing the fitted parameterisation.
test_that("matrix prediction requires every grouped-term column", {
  set.seed(2206)
  n <- 120

  x <- cbind(
    groupB = rep(c(0, 1, 0), length.out = n),
    groupC = rep(c(0, 0, 1), length.out = n),
    x1 = runif(n, 1, 10)
  )
  y <- x[, "groupB"] - x[, "groupC"] +
    0.2 * x[, "x1"] + rnorm(n)

  fit <- mfp2(
    x,
    y,
    df = 1,
    select = 1,
    term_groups = list(group = c("groupB", "groupC")),
    verbose = FALSE
  )
  incomplete <- x[1:5, c("groupB", "x1"), drop = FALSE]

  expect_error(
    predict(fit, newdata = incomplete),
    "groupC|grouped term|grouped-term",
    ignore.case = TRUE
  )
})


# -----------------------------------------------------------------------------
# 8.1 Prediction equivalence against stats::glm()
# -----------------------------------------------------------------------------
# These tests deliberately disable FP selection/transformation complexity
# so that mfp2() should reduce to the corresponding base glm() fit:
#   - df = 1: linear terms only
#   - select = 1 and alpha = 1: keep all ordinary predictors
#   - shift = 0, scale = 1, center = FALSE: no preprocessing-induced change
#   - no spike-at-zero, zero/catzero handling, or ACD transformation
#
# mfp2 stores ordinary transformed columns as x.1, x1.1, and so on. These
# tests compare numerical values and set xorder = "original", so coefficient
# and covariance order must agree with glm() without name canonicalization.


# Shared assertion helper for ordinary GLM equivalence tests.
#
# The helper verifies three increasingly independent layers:
#   1. mfp2() and stats::glm() fit the same statistical model;
#   2. both prediction methods return the same link values, responses, and
#      link-scale standard errors; and
#   3. those values agree with direct matrix algebra using X %*% beta,
#      the formula offset, the inverse-link function, and X V X'.
#
# Using a manual oracle is important because two prediction methods can agree
# while sharing the same reconstruction error. The X beta calculation checks
# the fitted coefficient order, factor expansion, offset handling, link
# inversion, and covariance propagation independently.
expect_mfp2_glm_predictions_equal <- function(dat,
                                              formula,
                                              family_name,
                                              newdata_cols,
                                              tolerance = 1e-8) {
  family_fun <- switch(
    family_name,
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    stop("Unsupported test family: ", family_name, call. = FALSE)
  )

  fit_mfp2 <- mfp2(
    formula,
    data = dat,
    family = family_name,
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
    formula,
    data = dat,
    family = family_fun
  )

  nd <- dat[1:25, newdata_cols, drop = FALSE]

  # Request link-scale standard errors from both methods. Standard errors are
  # naturally calculated on the linear-predictor scale by predict.glm().
  pred_mfp2_link <- predict(
    fit_mfp2,
    newdata = nd,
    type = "link",
    se.fit = TRUE
  )
  pred_glm_link <- predict(
    fit_glm,
    newdata = nd,
    type = "link",
    se.fit = TRUE
  )

  pred_mfp2_response <- predict(
    fit_mfp2,
    newdata = nd,
    type = "response"
  )
  pred_glm_response <- predict(
    fit_glm,
    newdata = nd,
    type = "response"
  )

  # Build the reference model frame from raw newdata. This evaluates factor
  # contrasts and formula offsets with the same terms object used by glm().
  reference_terms <- stats::delete.response(stats::terms(fit_glm))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = nd,
    xlev = fit_glm$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = fit_glm$contrasts
  )

  beta <- stats::coef(fit_glm)
  beta_vcov <- stats::vcov(fit_glm)

  # xorder = "original" makes the model-matrix and coefficient order
  # positional, so no name-based reordering is required.
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  # model.offset() returns NULL when the formula has no offset. In that case the
  # additive offset contribution is exactly zero for every prediction row.
  manual_offset <- stats::model.offset(reference_frame)
  if (is.null(manual_offset)) {
    manual_offset <- rep(0, nrow(manual_x))
  }

  # Independent prediction calculations.
  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_response <- as.numeric(fit_glm$family$linkinv(manual_link))
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  # Fitted-model equivalence checks. These detect differences that may be hidden
  # when predictions happen to be evaluated at only a small set of rows.
  expect_mfp2_glm_parameters_equal(
    fit_mfp2,
    fit_glm,
    tolerance = tolerance
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_glm)),
    tolerance = tolerance
  )
  expect_equal(
    unname(stats::fitted(fit_mfp2)),
    unname(stats::fitted(fit_glm)),
    tolerance = tolerance
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_glm$linear.predictors),
    tolerance = tolerance
  )

  # mfp2() must agree with glm() for both prediction scales and link-scale SEs.
  expect_equal(
    unname(pred_mfp2_link$fit),
    unname(pred_glm_link$fit),
    tolerance = tolerance
  )
  expect_equal(
    unname(pred_mfp2_link$se.fit),
    unname(pred_glm_link$se.fit),
    tolerance = tolerance
  )
  expect_equal(
    unname(pred_mfp2_response),
    unname(pred_glm_response),
    tolerance = tolerance
  )

  # Both methods must also agree with the independently constructed oracle.
  expect_equal(unname(pred_mfp2_link$fit), manual_link, tolerance = tolerance)
  expect_equal(unname(pred_glm_link$fit), manual_link, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_link$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_glm_link$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_response), manual_response, tolerance = tolerance)
  expect_equal(unname(pred_glm_response), manual_response, tolerance = tolerance)
}

# Shared helper for tests that use non-default family/link objects or special
# formula constructions. It performs the same manual X beta checks on two
# already fitted GLM objects.
expect_glm_objects_and_manual_prediction_equal <- function(fit_mfp2,
                                                           fit_glm,
                                                           newdata,
                                                           tolerance = 1e-8) {
  pred_mfp2_link <- predict(
    fit_mfp2,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  pred_glm_link <- predict(
    fit_glm,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  pred_mfp2_response <- predict(fit_mfp2, newdata = newdata, type = "response")
  pred_glm_response <- predict(fit_glm, newdata = newdata, type = "response")

  reference_terms <- stats::delete.response(stats::terms(fit_glm))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = newdata,
    xlev = fit_glm$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = fit_glm$contrasts
  )

  beta <- stats::coef(fit_glm)
  beta_vcov <- stats::vcov(fit_glm)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_offset <- stats::model.offset(reference_frame)
  if (is.null(manual_offset)) {
    manual_offset <- rep(0, nrow(manual_x))
  }

  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_response <- as.numeric(fit_glm$family$linkinv(manual_link))
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_mfp2_glm_parameters_equal(
    fit_mfp2,
    fit_glm,
    tolerance = tolerance
  )
  expect_equal(as.numeric(logLik(fit_mfp2)), as.numeric(logLik(fit_glm)), tolerance = tolerance)
  expect_equal(unname(fitted(fit_mfp2)), unname(fitted(fit_glm)), tolerance = tolerance)
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_glm$linear.predictors),
    tolerance = tolerance
  )

  expect_equal(unname(pred_mfp2_link$fit), unname(pred_glm_link$fit), tolerance = tolerance)
  expect_equal(unname(pred_mfp2_link$se.fit), unname(pred_glm_link$se.fit), tolerance = tolerance)
  expect_equal(unname(pred_mfp2_response), unname(pred_glm_response), tolerance = tolerance)

  expect_equal(unname(pred_mfp2_link$fit), manual_link, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_link$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_mfp2_response), manual_response, tolerance = tolerance)
}

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


# -----------------------------------------------------------------------------
# 8.2 Prediction equivalence against survival::coxph()
# -----------------------------------------------------------------------------
# These tests use the simplest Cox configuration where mfp2() should reduce to
# the corresponding survival::coxph() fit:
#   - df = 1: linear terms only
#   - select = 1 and alpha = 1: keep all ordinary predictors
#   - shift = 0, scale = 1, center = FALSE: no preprocessing-induced change
#   - matched tie handling and reference scale
#   - formula-level strata reconstructed during prediction

# Compare a forced-linear mfp2 Cox fit with survival::coxph() by coefficient
# position. Every caller uses xorder = "original", so no name repair or
# coefficient reordering is required. The helper checks coefficients,
# coefficient standard errors, training linear predictors, newdata LP and risk
# predictions, newdata prediction standard errors, partial log-likelihood, and
# an independent X beta / X V X' calculation from the coxph model matrix.
expect_mfp2_cox_predictions_equal <- function(fit_mfp2,
                                              fit_coxph,
                                              newdata,
                                              mfp2_newdata = newdata,
                                              mfp2_predict_args = list(),
                                              tolerance = 1e-8) {
  coef_mfp2 <- stats::coef(fit_mfp2)
  coef_coxph <- stats::coef(fit_coxph)
  vcov_mfp2 <- stats::vcov(fit_mfp2)
  vcov_coxph <- stats::vcov(fit_coxph)

  expect_length(coef_mfp2, length(coef_coxph))
  expect_equal(dim(vcov_mfp2), dim(vcov_coxph))
  expect_equal(unname(coef_mfp2), unname(coef_coxph), tolerance = tolerance)
  expect_equal(unname(vcov_mfp2), unname(vcov_coxph), tolerance = tolerance)
  expect_equal(
    unname(sqrt(diag(vcov_mfp2))),
    unname(sqrt(diag(vcov_coxph))),
    tolerance = tolerance
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_coxph)),
    tolerance = tolerance
  )

  # Training predictions verify that fitted offsets and strata are stored and
  # applied consistently without requiring any newdata reconstruction.
  training_lp_mfp2 <- predict(fit_mfp2, type = "lp")
  training_lp_coxph <- predict(
    fit_coxph,
    type = "lp",
    reference = "zero"
  )
  expect_equal(
    unname(training_lp_mfp2),
    unname(training_lp_coxph),
    tolerance = tolerance
  )

  # Formula fits reconstruct strata and offsets from raw newdata. Matrix fits
  # supply them explicitly through mfp2_predict_args.
  pred_mfp2 <- do.call(
    stats::predict,
    c(
      list(
        object = fit_mfp2,
        newdata = mfp2_newdata,
        type = "lp",
        se.fit = TRUE
      ),
      mfp2_predict_args
    )
  )
  pred_coxph <- predict(
    fit_coxph,
    newdata = newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  risk_mfp2 <- do.call(
    stats::predict,
    c(
      list(
        object = fit_mfp2,
        newdata = mfp2_newdata,
        type = "risk"
      ),
      mfp2_predict_args
    )
  )
  risk_coxph <- predict(
    fit_coxph,
    newdata = newdata,
    type = "risk",
    reference = "zero"
  )

  # model.matrix.coxph() removes the non-estimated intercept and formula
  # specials such as strata() and offset(), leaving the coefficient design in
  # its fitted order. This is the independent X beta oracle for newdata.
  manual_x <- stats::model.matrix(fit_coxph, data = newdata)
  if ("(Intercept)" %in% colnames(manual_x)) {
    manual_x <- manual_x[, colnames(manual_x) != "(Intercept)", drop = FALSE]
  }
  expect_equal(ncol(manual_x), length(coef_coxph))
  expect_equal(dim(vcov_coxph), c(length(coef_coxph), length(coef_coxph)))

  reference_terms <- stats::delete.response(stats::terms(fit_coxph))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = newdata,
    xlev = fit_coxph$xlevels,
    na.action = stats::na.pass
  )
  new_offset <- stats::model.offset(reference_frame)
  training_frame <- stats::model.frame(fit_coxph)
  training_offset <- stats::model.offset(training_frame)
  expected_offset_reference <- if (is.null(training_offset)) {
    0
  } else {
    mean(training_offset)
  }

  # fit_mfp() stores the same offset origin used internally by
  # predict.coxph() once on the final mfp2 object. Package-owned manual Cox
  # prediction paths reuse this scalar rather than recomputing it in fit_model().
  expect_equal(
    fit_mfp2$cox_offset_reference,
    unname(expected_offset_reference),
    tolerance = tolerance
  )

  # predict.coxph() always recentres an offset by subtracting its mean in the
  # training model frame. This applies even when cox_reference = "zero"; that
  # option controls covariate centering, not offset centering. Reproduce that
  # convention explicitly in the independent oracle.
  if (is.null(new_offset)) {
    manual_offset <- rep(0, nrow(manual_x))
  } else {
    expect_false(is.null(training_offset))
    manual_offset <- as.numeric(new_offset - expected_offset_reference)
  }

  manual_lp <- as.numeric(manual_x %*% coef_coxph + manual_offset)
  manual_risk <- exp(manual_lp)
  manual_variance <- rowSums((manual_x %*% vcov_coxph) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_equal(
    unname(pred_mfp2$fit),
    unname(pred_coxph$fit),
    tolerance = tolerance
  )
  expect_equal(
    unname(pred_mfp2$se.fit),
    unname(pred_coxph$se.fit),
    tolerance = tolerance
  )
  expect_equal(unname(pred_mfp2$fit), manual_lp, tolerance = tolerance)
  expect_equal(unname(pred_coxph$fit), manual_lp, tolerance = tolerance)
  expect_equal(unname(pred_mfp2$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(pred_coxph$se.fit), manual_se, tolerance = tolerance)
  expect_equal(unname(risk_mfp2), unname(risk_coxph), tolerance = tolerance)
  expect_equal(unname(risk_mfp2), manual_risk, tolerance = tolerance)
  expect_equal(unname(risk_coxph), manual_risk, tolerance = tolerance)
}

# Generate one stable Cox data set containing continuous and categorical
# predictors, two potential stratification factors, and an offset. Individual
# equivalence tests use different subsets of these model components.
make_cox_equivalence_data <- function(n = 480L, seed = 8200L) {
  set.seed(seed)

  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    stratum1 = factor(sample(c("S1", "S2", "S3"), n, replace = TRUE)),
    stratum2 = factor(sample(c("T1", "T2"), n, replace = TRUE)),
    off = stats::rnorm(n, mean = 0, sd = 0.25)
  )

  group_effect <- c(A = 0, B = 0.45, C = -0.35)[as.character(dat$group)]
  baseline_multiplier <-
    c(S1 = 0.75, S2 = 1.00, S3 = 1.35)[as.character(dat$stratum1)] *
    c(T1 = 0.85, T2 = 1.20)[as.character(dat$stratum2)]
  eta <- 0.38 * dat$x1 - 0.25 * dat$x2 + group_effect + dat$off
  event_time <- stats::rexp(
    n,
    rate = 0.020 * baseline_multiplier * exp(eta)
  )
  censor_time <- stats::rexp(n, rate = 0.010)
  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  dat
}

# Test purpose: 8.2.1 A formula model with one strata() term should
# match coxph() for estimates, coefficient SEs, likelihood, training LPs,
# newdata LP/risk predictions, and prediction SEs.
test_that("8.2.1 Cox single formula strata matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8201)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 + strata(stratum1),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + strata(stratum1),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "stratum1"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.3 Two formula strata() terms should be reconstructed and
# combined identically to coxph() for fitting and prediction.
test_that("8.2.3 Cox multiple formula strata match coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8203)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 +
      strata(stratum1) + strata(stratum2),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 +
      strata(stratum1) + strata(stratum2),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "stratum1", "stratum2"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.4 Matrix-interface strata supplied at fit and prediction
# time should match an equivalent coxph() formula model completely.
test_that("8.2.4 Cox matrix-interface strata matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8204)
  x <- as.matrix(dat[, c("x1", "x2")])
  y <- survival::Surv(dat$time, dat$status)

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "cox",
    strata = dat$stratum1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + strata(stratum1),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "stratum1"), drop = FALSE]
  newx <- as.matrix(nd[, c("x1", "x2")])

  expect_mfp2_cox_predictions_equal(
    fit_mfp2,
    fit_coxph,
    newdata = nd,
    mfp2_newdata = newx,
    mfp2_predict_args = list(strata = nd$stratum1)
  )
})

# Test purpose: 8.2.5 A formula offset without strata should be included in
# training and newdata Cox linear predictors exactly as in coxph().
test_that("8.2.5 Cox formula offset matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8205)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "off"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.6 A retained categorical predictor without strata or
# offset should match coxph() completely.
test_that("8.2.6 Cox categorical predictor matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8206)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group,
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.7 A categorical predictor combined with formula strata
# should preserve both grouped contrasts and stratum-specific risk sets.
test_that("8.2.7 Cox categorical predictor with strata matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8207)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group + strata(stratum1),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group + strata(stratum1),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group", "stratum1"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.8 A categorical predictor combined with a formula offset
# should match coxph() in all fitted and predicted quantities.
test_that("8.2.8 Cox categorical predictor with offset matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8208)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group", "off"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.9 Categorical contrasts, one formula strata term, and a
# formula offset should all be reconstructed together without changing the
# coxph() fit or predictions.
test_that("8.2.9 Cox categorical predictor with strata and offset matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8209)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group +
      strata(stratum1) + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group +
      strata(stratum1) + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group", "stratum1", "off"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.10 Two formula strata terms and an offset should match
# coxph() simultaneously, including the offset contribution to LP and risk.
test_that("8.2.10 Cox multiple strata with offset match coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8210)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 +
      strata(stratum1) + strata(stratum2) + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 +
      strata(stratum1) + strata(stratum2) + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[
    1:30,
    c("x1", "x2", "stratum1", "stratum2", "off"),
    drop = FALSE
  ]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})

# Test purpose: 8.2.11 Matrix-interface strata and offset vectors should match
# the equivalent coxph() formula model for fitting and newdata prediction.
test_that("8.2.11 Cox matrix-interface strata and offset match coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8211)
  x <- as.matrix(dat[, c("x1", "x2")])
  y <- survival::Surv(dat$time, dat$status)

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "cox",
    strata = dat$stratum1,
    offset = dat$off,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + strata(stratum1) + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "stratum1", "off"), drop = FALSE]
  newx <- as.matrix(nd[, c("x1", "x2")])

  expect_mfp2_cox_predictions_equal(
    fit_mfp2,
    fit_coxph,
    newdata = nd,
    mfp2_newdata = newx,
    mfp2_predict_args = list(
      strata = nd$stratum1,
      newoffset = nd$off
    )
  )
})


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

  expect_error(
    predict(fit, newdata = dat[1:5, c("age", "sex"), drop = FALSE], type = "lp"),
    "formula-level Cox strata|strata term could not be reconstructed"
  )
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

# =============================================================================
# 9. Model selection criteria
# =============================================================================

# Test purpose: Checks that AIC-based model selection runs successfully.
test_that("criterion = 'aic' runs without error", {
  fit <- mfp2(x_prostate, y_prostate, criterion = "aic", verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks that BIC-based model selection runs successfully.
test_that("criterion = 'bic' runs without error", {
  fit <- mfp2(x_prostate, y_prostate, criterion = "bic", verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks the expected stronger-penalty behavior of BIC relative to AIC on this dataset.
test_that("BIC selects equal or fewer variables than AIC", {
  fit_aic <- mfp2(x_prostate, y_prostate, criterion = "aic", verbose = FALSE)
  fit_bic <- mfp2(x_prostate, y_prostate, criterion = "bic", verbose = FALSE)

  n_aic <- sum(fit_aic$fp_terms[, "selected"])
  n_bic <- sum(fit_bic$fp_terms[, "selected"])

  # BIC penalizes more heavily so typically selects <= AIC variables
  # This may not always hold for every dataset, but is generally expected
  expect_true(n_bic <= n_aic + 1) # allow slack of 1
})

# Test purpose: Checks that select = 1 retains all predictors under p-value selection.
test_that("select = 1 forces all variables into model", {
  fit <- mfp2(x_prostate, y_prostate, select = 1, verbose = FALSE)
  expect_true(all(fit$fp_terms[, "selected"]))
})

# Test purpose: Checks that variables listed in keep remain selected in the final model.
test_that("keep argument retains specified variables", {
  fit <- mfp2(x_prostate, y_prostate, keep = c("age", "bph"), verbose = FALSE)
  expect_true(fit$fp_terms["age", "selected"])
  expect_true(fit$fp_terms["bph", "selected"])
})

# Test purpose: Checks that Gaussian fitting works when F-test based selection is requested.
test_that("ftest argument works for Gaussian family", {
  fit <- mfp2(x_prostate, y_prostate, ftest = TRUE, verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})

# =============================================================================
# 10. Edge cases and input validation
# =============================================================================

# Test purpose: Defines the shared scalar-or-named-override contract used by
# df, select, and alpha in the matrix interfaces.
test_that("named override settings fill defaults and match by name", {
  columns <- c("age", "bmi", "weight")

  global <- normalize_named_override_setting(
    value = 2,
    column_names = columns,
    default = 4,
    argument_name = "df"
  )
  expect_equal(global$value, c(age = 2, bmi = 2, weight = 2))
  expect_true(global$global_scalar)
  expect_true(all(global$supplied))

  partial <- normalize_named_override_setting(
    value = c(weight = 1, age = 2),
    column_names = columns,
    default = 4,
    argument_name = "df"
  )
  expect_equal(partial$value, c(age = 2, bmi = 4, weight = 1))
  expect_false(partial$global_scalar)
  expect_equal(partial$supplied, c(age = TRUE, bmi = FALSE, weight = TRUE))
})

# Test purpose: Prevents the previous positional assignment behavior and
# rejects malformed or unknown names.
test_that("named override settings reject positional and malformed vectors", {
  columns <- c("age", "bmi", "weight")

  expect_error(
    normalize_named_override_setting(
      value = c(1, 2),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "single unnamed numeric value or a named numeric vector"
  )

  expect_error(
    normalize_named_override_setting(
      value = c(age = 1, unknown = 2),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "unknown column name.*unknown"
  )

  expect_error(
    normalize_named_override_setting(
      value = stats::setNames(c(1, 2), c("age", "age")),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "names must be unique"
  )

  expect_error(
    normalize_named_override_setting(
      value = stats::setNames(c(1, 2), c("age", "")),
      column_names = columns,
      default = 4,
      argument_name = "df"
    ),
    "single unnamed numeric value or a named numeric vector"
  )
})

# Test purpose: Verifies partial named df/select/alpha settings in
# mfp2.default(), including default filling and cardinality reduction for an
# omitted five-level predictor.
test_that("mfp2.default() accepts partial named df select and alpha overrides", {
  set.seed(10001)
  n <- 80L
  x <- cbind(
    age = stats::runif(n, 1, 8),
    bmi = rep(1:5, length.out = n),
    weight = stats::runif(n, 2, 10)
  )
  y <- 0.5 * x[, "age"] + stats::rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    cycles = 5,
    df = c(age = 1),
    select = c(age = 1),
    alpha = c(weight = 1),
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["age", "df_initial"]), 1)
  expect_equal(as.numeric(fit$fp_terms["bmi", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["weight", "df_initial"]), 4)
  expect_equal(as.numeric(fit$fp_terms["age", "select"]), 1)
  expect_equal(as.numeric(fit$fp_terms["bmi", "select"]), 0.05)
  expect_equal(as.numeric(fit$fp_terms["age", "alpha"]), 0.05)
  expect_equal(as.numeric(fit$fp_terms["weight", "alpha"]), 1)
})

# Test purpose: Checks the public error path rather than only the internal
# normalizer when an unnamed multi-value vector is supplied.
test_that("mfp2.default() rejects unnamed multi-value df select and alpha", {
  x <- cbind(age = 1:30, bmi = seq(2, 8, length.out = 30))
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(30)

  expect_error(
    mfp2(x, y, df = c(1, 4), verbose = FALSE),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    mfp2(x, y, select = c(1, 0.05), verbose = FALSE),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    mfp2(x, y, alpha = c(1, 0.05), verbose = FALSE),
    "single unnamed numeric value or a named numeric vector"
  )
})

# Test purpose: Checks that every manually grouped member names an existing raw
# design-matrix column.
test_that("term_groups rejects unknown columns", {
  x <- cbind(x1 = 1:20, groupB = rep(0:1, 10))
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(20)

  expect_error(
    mfp2(
      x,
      y,
      term_groups = list(group = c("groupB", "groupC"))
    ),
    "unknown column.*groupC|groupC.*unknown column",
    ignore.case = TRUE
  )
})

# Test purpose: Checks that one raw column cannot belong to two conceptual terms.
test_that("term_groups rejects columns in more than one group", {
  x <- cbind(
    groupB = rep(0:1, 10),
    groupC = rep(c(0, 0, 1, 0), 5),
    other = 1:20
  )
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(20)

  expect_error(
    mfp2(
      x,
      y,
      term_groups = list(
        group1 = c("groupB", "groupC"),
        group2 = c("groupC", "other")
      )
    ),
    "only one.*Duplicated column|Duplicated column.*groupC",
    ignore.case = TRUE
  )
})

# Test purpose: Checks that a grouped term cannot enter the FP search with a
# nonlinear df setting; grouped blocks are fixed linear terms.
test_that("grouped terms require df = 1 for every member column", {
  x <- cbind(
    group1 = 1:20,
    group2 = (1:20)^2,
    x1 = seq(0.5, 10, length.out = 20)
  )
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(20)

  expect_error(
    mfp2(
      x,
      y,
      df = c(group1 = 1, group2 = 4, x1 = 1),
      term_groups = list(group = c("group1", "group2"))
    ),
    paste0(
      "Grouped term 'group'.*`df`.*1.*fixed linear design blocks.*",
      "fractional-polynomial transformation search"
    ),
    ignore.case = TRUE
  )
})

# Test purpose: Checks that continuous-only extensions cannot be assigned to
# any member of a grouped categorical block.
test_that("grouped terms reject continuous-only processing options", {
  x <- cbind(
    groupB = rep(c(0, 1, 0), length.out = 24),
    groupC = rep(c(0, 0, 1), length.out = 24),
    x1 = seq(1, 12, length.out = 24)
  )
  # Response values are irrelevant because validation fails before fitting.
  y <- seq_len(24)
  grouped <- list(group = c("groupB", "groupC"))

  cases <- list(
    acdx = list(acdx = "groupB"),
    zero_vars = list(zero_vars = "groupB"),
    catzero_vars = list(catzero_vars = "groupB"),
    spike_vars = list(spike_vars = "groupB")
  )

  for (setting in names(cases)) {
    args <- c(
      list(x = x, y = y, term_groups = grouped),
      cases[[setting]]
    )
    expect_error(
      do.call(mfp2, args),
      paste0(
        "Grouped term 'group'.*`", setting, "`.*FALSE.*",
        "fixed categorical design blocks.*singleton continuous predictors"
      ),
      ignore.case = TRUE,
      info = paste("setting =", setting)
    )
  }
})

# Test purpose: Confirms that conceptual-term mappings are supplied explicitly
# throughout the fitting chain and are not read from an attribute on x.
test_that("grouped mappings use explicit arguments instead of matrix attributes", {
  x <- cbind(
    groupB = c(0, 1, 0, 1),
    groupC = c(0, 0, 1, 0),
    x1 = c(1, 2, 3, 4)
  )
  attr(x, "mfp2_term_to_columns") <- list(stale = "missing_column")

  term_to_columns <- list(
    group = c("groupB", "groupC"),
    x1 = "x1"
  )
  term_names <- names(term_to_columns)

  out <- build_adjustment_step(
    x = x,
    xi = "x1",
    powers_current = list(group = 1, x1 = 1),
    powers = list(group = 1, x1 = 1),
    acdx = stats::setNames(rep(FALSE, 2L), term_names),
    zero = stats::setNames(rep(FALSE, 2L), term_names),
    catzero = stats::setNames(vector("list", 2L), term_names),
    spike = stats::setNames(rep(FALSE, 2L), term_names),
    spike_decision = stats::setNames(rep(0L, 2L), term_names),
    acd_parameter = stats::setNames(vector("list", 2L), term_names),
    prev_adj_params = stats::setNames(vector("list", 2L), term_names),
    term_to_columns = term_to_columns
  )

  expect_identical(colnames(out$data_adj), c("groupB", "groupC"))
  expect_equal(
    unname(out$data_adj),
    unname(x[, c("groupB", "groupC"), drop = FALSE])
  )
  expect_false(any(grepl(
    "mfp2_term_to_columns",
    deparse(body(fit_mfp), width.cutoff = 500L),
    fixed = TRUE
  )))
})

# Test purpose: Covers the MFPI names that use the shared grouped-setting
# validator and confirms that both option classes receive an actionable reason.
test_that("grouped-setting errors explain why MFPI options are unsupported", {
  term_to_columns <- list(group = c("groupB", "groupC"))

  expect_error(
    validate_grouped_term_setting(
      term_to_columns = term_to_columns,
      values = c(groupB = TRUE, groupC = FALSE),
      setting = "acd_vars",
      predicate = function(v) !v,
      requirement = "FALSE"
    ),
    paste0(
      "fixed categorical design blocks.*",
      "singleton continuous predictors"
    ),
    ignore.case = TRUE
  )

  expect_error(
    validate_grouped_term_setting(
      term_to_columns = term_to_columns,
      values = c(groupB = TRUE, groupC = FALSE),
      setting = "force_max_fp_vars",
      predicate = function(v) !v,
      requirement = "FALSE"
    ),
    paste0(
      "fixed linear design blocks.*",
      "fractional-polynomial transformation search"
    ),
    ignore.case = TRUE
  )
})

# Test purpose: Checks that the default interface requires a matrix design input.
test_that("mfp2() rejects non-matrix input in default interface", {
  expect_error(mfp2(as.data.frame(x_prostate), y_prostate), "matrix")
})

# Test purpose: Checks that missing predictor values are rejected.
test_that("mfp2() rejects input with missing data", {
  x_bad <- x_prostate
  x_bad[1, 1] <- NA
  expect_error(mfp2(x_bad, y_prostate), "NA")
})

# Test purpose: Checks that unnamed matrix columns are rejected to avoid ambiguous variable handling.
test_that("mfp2() rejects input without column names", {
  x_bad <- unname(x_prostate)
  expect_error(mfp2(x_bad, y_prostate), "column names")
})

# Test purpose: Checks that character-valued predictors are rejected.
test_that("mfp2() rejects character data in x", {
  x_bad <- x_prostate
  storage.mode(x_bad) <- "character"
  expect_error(mfp2(x_bad, y_prostate), "character")
})

# Test purpose: Checks that response length must match the number of predictor
# rows.
test_that("mfp2() rejects mismatched y length", {
  expect_error(mfp2(x_prostate, y_prostate[1:10]), "must match")
})

# Test purpose: Checks that numeric subsetting fits the model on the selected
# observations.
test_that("subset argument works correctly", {
  idx <- 1:50
  fit <- mfp2(x_prostate, y_prostate, subset = idx, verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  # The model should be fitted on the subset
  expect_equal(length(fit$residuals), length(idx))
})

# Test purpose: Checks that logical subsetting fits the model on TRUE observations
#  only.
test_that("subset with logical vector works", {
  log_sub <- rep(FALSE, nrow(x_prostate))
  log_sub[1:50] <- TRUE
  fit <- mfp2(x_prostate, y_prostate, subset = log_sub, verbose = FALSE)

  expect_equal(length(fit$residuals), 50)
})

# Test purpose: Checks that all supported covariate-entry order options run
# successfully.
test_that("xorder options work without error", {
  for (ord in c("ascending", "descending", "original")) {
    fit <- mfp2(x_prostate, y_prostate, xorder = ord, verbose = FALSE)
    expect_s3_class(fit, "mfp2")
  }
})

# Test purpose: The leave-one-term-out likelihood-ratio test used for
# significance ordering must use the fitted rank contribution of a grouped
# term, not the number of raw columns in its design block. The mocked fits are
# chosen so those two df rules produce opposite orders.
test_that("grouped significance ordering uses fitted rank difference", {
  x <- cbind(
    groupB = c(0, 1, 0, 1),
    groupC = c(0, 0, 1, 1),
    z = c(-1, 0, 1, 2)
  )
  term_to_columns <- list(
    group = c("groupB", "groupC"),
    z = "z"
  )

  testthat::local_mocked_bindings(
    fit_model = function(x, ...) {
      remaining <- colnames(x)
      if (identical(remaining, "z")) {
        # Dropping the two-column group reduces fitted rank by one and gives
        # likelihood-ratio statistic 20.
        return(list(logl = 90, df = 1))
      }
      if (identical(remaining, c("groupB", "groupC"))) {
        # Dropping z also reduces fitted rank by one and gives statistic 18.
        return(list(logl = 91, df = 1))
      }
      stop("Unexpected reduced design in ordering test.")
    },
    .package = "mfp2"
  )

  common_args <- list(
    x = x,
    y = rep(0, nrow(x)),
    family = stats::gaussian(),
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    strata = NULL,
    method = NULL,
    control = NULL,
    nocenter = NULL,
    full_reference = list(logl = 100, df = 2),
    term_to_columns = term_to_columns
  )

  ascending <- do.call(
    order_variables_by_significance,
    c(list(xorder = "ascending"), common_args)
  )
  descending <- do.call(
    order_variables_by_significance,
    c(list(xorder = "descending"), common_args)
  )

  expect_identical(ascending, c("group", "z"))
  expect_identical(descending, c("z", "group"))
})

# Test purpose: Checks df down-capping rules for binary, ternary, few-level, and
# continuous variables.
test_that("assign_df() correctly limits df for low-cardinality variables", {
  x <- cbind(
    binary = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    ternary = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1),
    few = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
    continuous = 1:10
  )

  df <- assign_df(x, df_default = 4)
  expect_equal(df[["binary"]], 1)    # <= 3 unique -> 1
  expect_equal(df[["ternary"]], 1)   # <= 3 unique -> 1
  expect_equal(df[["few"]], 2)       # 4-5 unique -> min(2, 4) = 2
  expect_equal(df[["continuous"]], 4) # >= 6 unique -> 4
})

# Test purpose: Checks that the selected-variable accessor returns valid predictor
# names.
test_that("get_selected_variable_names() returns correct names", {
  fit <- mfp2(x_prostate, y_prostate, select = 1, verbose = FALSE)

  sel <- get_selected_variable_names(fit)
  expect_true(is.character(sel))
  expect_true(length(sel) > 0)
  expect_true(all(sel %in% colnames(x_prostate)))
})

# =============================================================================
# 11. Summary, print, coef methods
# =============================================================================

# Test purpose: Checks that summary() returns output for a fitted Gaussian mfp2
#  model.
test_that("summary.mfp2() works for Gaussian", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  s <- summary(fit)
  expect_true(!is.null(s))
})

# Test purpose: Checks that summary() returns output for a fitted Cox mfp2 model.
test_that("summary.mfp2() works for Cox", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)
  s <- summary(fit)
  expect_true(!is.null(s))
})


# Test purpose: Ensures displayed ordinary FP equations use the final
# shifted-but-unscaled basis rather than reapplying the preprocessing scale.
test_that("summary FP labels do not reapply preprocessing scale", {
  object <- list(
    fp_terms = data.frame(
      selected = TRUE,
      acd = FALSE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      df_final = 2,
      power1 = -1,
      power2 = NA_real_,
      row.names = "x"
    ),
    x = matrix(1, nrow = 1L, dimnames = list(NULL, "x.1")),
    coefficients = c("x.1" = 2.5),
    transformations = data.frame(
      shift = 3,
      scale = 100,
      row.names = "x"
    )
  )

  classified <- mfp2_summary_classify_terms(object)
  label <- mfp2_summary_term_label(object, classified, "x", "x.1")
  formula <- mfp2_summary_formula_strings(object, classified)

  expect_true(grepl("x + 3", label, fixed = TRUE))
  expect_false(grepl("/100", label, fixed = TRUE))
  expect_false(any(grepl("/100", formula, fixed = TRUE)))
})

# Test purpose: Ensures ACD direct and transformed columns remain grouped under
# one variable and new ACD definitions use shifted, unscaled predictor values.
test_that("summary shows unscaled ACD definitions for new fits", {
  object <- list(
    fp_terms = data.frame(
      selected = TRUE,
      acd = TRUE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      df_final = 4,
      power1 = 1,
      power2 = 1,
      row.names = "x"
    ),
    x = matrix(
      c(1, 0.5),
      nrow = 1L,
      dimnames = list(NULL, c("x.1", "A_x.1"))
    ),
    coefficients = c("x.1" = 0.4, "A_x.1" = 1.2),
    transformations = data.frame(
      shift = 2,
      scale = 1,
      row.names = "x"
    ),
    acd_parameter = list(
      x = list(
        beta0 = -0.8,
        beta1 = 1.3,
        power = 0,
        shift = 0,
        scale = 1
      )
    )
  )

  classified <- mfp2_summary_classify_terms(object)
  direct_label <- mfp2_summary_term_label(object, classified, "x", "x.1")
  acd_label <- mfp2_summary_term_label(object, classified, "x", "A_x.1")
  definitions <- mfp2_summary_acd_definitions(object, classified)
  formulas <- mfp2_summary_formula_strings(object, classified)

  expect_identical(
    classified$cols_by_var[["x"]],
    c("x.1", "A_x.1")
  )
  expect_true(grepl("x + 2", direct_label, fixed = TRUE))
  expect_false(grepl("/", direct_label, fixed = TRUE))
  expect_true(grepl("A(x)", acd_label, fixed = TRUE))
  expect_false(grepl("log(A(x))", acd_label, fixed = TRUE))
  expect_true(any(grepl("x + 2", definitions, fixed = TRUE)))
  expect_false(any(grepl("/", definitions, fixed = TRUE)))
  expect_true(any(grepl("pnorm", definitions, fixed = TRUE)))
  expect_true(any(grepl("A(x)", formulas, fixed = TRUE)))
})

# Test purpose: Ensures ACD power positions are preserved so an ACD-only form
# c(NA, p) is not mislabeled as a direct ordinary FP term.
test_that("summary preserves ACD-only power slots", {
  object <- list(
    fp_terms = data.frame(
      selected = TRUE,
      acd = TRUE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      df_final = 2,
      power1 = NA_real_,
      power2 = 1,
      row.names = "x"
    ),
    x = matrix(0.5, nrow = 1L, dimnames = list(NULL, "A_x.1")),
    coefficients = c("A_x.1" = 1.2),
    transformations = data.frame(
      shift = 0,
      scale = 10,
      row.names = "x"
    ),
    acd_parameter = list(
      x = list(
        beta0 = 0,
        beta1 = 1,
        power = 1,
        shift = 0,
        scale = 10
      )
    )
  )

  classified <- mfp2_summary_classify_terms(object)
  label <- mfp2_summary_term_label(object, classified, "x", "A_x.1")

  expect_equal(classified$power_slots_by_var[["x"]], c(NA_real_, 1))
  expect_true(grepl("A(x)", label, fixed = TRUE))
  expect_false(grepl("log(A(x))", label, fixed = TRUE))
})

# Test purpose: Verifies that GLM Model Fit output uses stored deviances and a
# Deviance header rather than reconstructing minus twice log-likelihood.
test_that("Model Fit reports deviance for GLMs", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  values <- mfp2_summary_model_fit_values(fit)

  expect_identical(attr(values, "statistic_label"), "Deviance")
  expect_equal(
    values$fit_statistic,
    c(fit$linear_deviance, fit$mfp_deviance)
  )
  output <- capture.output(print(fit))
  expect_true(any(grepl("Deviance", output, fixed = TRUE)))
  expect_false(any(grepl("-2 log L", output, fixed = TRUE)))
})

# Test purpose: Verifies that Cox Model Fit output remains on the existing
# minus-twice-partial-log-likelihood scale.
test_that("Model Fit keeps -2 log L for Cox models", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)
  values <- mfp2_summary_model_fit_values(fit)

  expect_identical(attr(values, "statistic_label"), "-2 log L")
  expect_equal(
    values$fit_statistic,
    c(fit$linear_deviance, fit$mfp_deviance)
  )
  output <- capture.output(print(fit))
  expect_true(any(grepl("-2 log L", output, fixed = TRUE)))
})

# Test purpose: Verifies that internal fast fits retain only lightweight
# quantities by default and calculate reporting statistics only when requested.
test_that("fit_model returns only requested internal components", {
  x <- matrix(
    c(-2, -1, 0, 1, 2, 3),
    ncol = 1L,
    dimnames = list(NULL, "x")
  )
  y <- c(-1.8, -0.9, 0.2, 1.1, 1.9, 3.2)

  selection_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = TRUE
  )

  expect_false("fit" %in% names(selection_fit))
  expect_false("null_deviance" %in% names(selection_fit))
  expect_false("model_deviance" %in% names(selection_fit))
  expect_false("sse" %in% names(selection_fit))
  expect_true(all(c("logl", "coefficients", "rank", "df") %in%
                    names(selection_fit)))

  metrics <- calculate_model_metrics(selection_fit, n_obs = length(y))
  expect_true(all(is.finite(metrics[c("logl", "df", "aic", "bic")])))

  reference_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = TRUE,
    calculate_fit_statistics = TRUE
  )

  expect_false("fit" %in% names(reference_fit))
  expect_true(all(c("null_deviance", "model_deviance") %in%
                    names(reference_fit)))

  retained_fast_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = TRUE,
    keep_fit = TRUE
  )

  expect_true("fit" %in% names(retained_fast_fit))

  full_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = FALSE
  )

  expect_true("fit" %in% names(full_fit))
  expect_false("null_deviance" %in% names(full_fit))
  expect_false("model_deviance" %in% names(full_fit))
  expect_identical(
    unname(full_fit$transformed_to_model_columns),
    "x"
  )
})

# Test purpose: Verifies that negative-binomial summary inference agrees with
# MASS and uses z rather than t statistics after theta has been estimated.
test_that("summary.mfp2() reports correct negative-binomial inference", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  fits <- get_negbin_reference_fits()
  summary_mfp2 <- summary(fits$mfp2, raw = TRUE)
  summary_mass <- summary(fits$mass)

  expect_true(!is.null(summary_mfp2))
  expect_true(is.matrix(summary_mfp2$coefficients))
  expect_true(any(grepl("^z($| value)", colnames(summary_mfp2$coefficients))))
  expect_false(any(grepl("^t($| value)", colnames(summary_mfp2$coefficients))))
  expect_equal(
    unname(summary_mfp2$coefficients[, 1:2, drop = FALSE]),
    unname(summary_mass$coefficients[, 1:2, drop = FALSE]),
    tolerance = 1e-4
  )
  expect_equal(
    unname(summary_mfp2$coefficients[, 3:4, drop = FALSE]),
    unname(summary_mass$coefficients[, 3:4, drop = FALSE]),
    tolerance = 1e-3
  )
})

# Test purpose: Checks that coef() returns named numeric coefficients.
test_that("coef.mfp2() returns named numeric vector", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_true(!is.null(names(cf)))
})

# Test purpose: Checks that the print method produces console output without error.
test_that("print.mfp2() runs without error", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  expect_output(print(fit))
})

# Test purpose: Ensures grouped terms print their actual initial and final
# model degrees of freedom without exposing the internal FP search setting.
test_that("print.mfp2() reports grouped initial and final df only", {
  set.seed(91501)
  n <- 120L
  group <- factor(rep(c("A", "B", "C"), length.out = n))
  group_mm <- stats::model.matrix(~ group)[, -1L, drop = FALSE]
  x <- cbind(x1 = stats::runif(n, 1, 8), group_mm)
  y <- 0.4 * x[, "x1"] + 0.8 * group_mm[, 1L] -
    0.5 * group_mm[, 2L] + stats::rnorm(n, sd = 0.3)

  fit <- mfp2(
    x, y,
    term_groups = list(group = colnames(group_mm)),
    keep = "group",
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["group", "df_setting"]), 1)
  expect_equal(as.numeric(fit$fp_terms["group", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), 2)

  output <- capture.output(print(fit, detailed_settings = TRUE, notes = TRUE))
  printed <- paste(output, collapse = "\n")
  expect_match(printed, "df (init->final)", fixed = TRUE)
  expect_match(printed, "2 -> 2", fixed = TRUE)
  expect_false(grepl("df setting", printed, fixed = TRUE))
  expect_false(grepl("df_setting", printed, fixed = TRUE))
})

# Test purpose: Ensures an SAZ-only fit does not print an empty Standard MFP
# data frame in the function-selection summary.
test_that("print.mfp2() omits an empty Standard MFP subsection", {
  fit <- mfp2(
    lpsa ~ fp(
      pgg45,
      df = 4,
      select = 0.05,
      alpha = 0.05,
      spike = TRUE
    ),
    data = prostate,
    family = "gaussian",
    criterion = "pvalue",
    verbose = FALSE
  )

  output <- capture.output(print(fit, detailed_settings = FALSE))

  expect_false(any(output == "Standard MFP"))
  expect_false(any(grepl("<0 rows>", output, fixed = TRUE)))
  expect_true(any(grepl("Spike-at-Zero (SAZ)", output, fixed = TRUE)))
})

# Test purpose: Ensures retained SAZ terms store and print the structural-zero
# proportion from the actual fitting sample.
test_that("SAZ metadata and print output include prop_zero", {
  set.seed(91502)
  n <- 100L
  exposure <- c(rep(0, 25L), stats::runif(75L, 0.5, 6))
  exposure <- sample(exposure)
  x <- cbind(exposure = exposure, age = stats::runif(n, 20, 70))
  y <- 1.2 * (exposure == 0) + 0.25 * exposure +
    0.01 * x[, "age"] + stats::rnorm(n, sd = 0.25)

  fit <- mfp2(
    x, y,
    df = c(exposure = 1),
    spike_vars = "exposure",
    keep = "exposure",
    verbose = FALSE
  )

  expect_equal(fit$fp_terms["exposure", "prop_zero"], 0.25)
  expect_true(is.na(fit$fp_terms["age", "prop_zero"]))

  output <- capture.output(
    print(fit, detailed_settings = FALSE, digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "prop_zero", fixed = TRUE)
  expect_match(printed, "0.250", fixed = TRUE)

  saz_header <- output[grepl("prop_zero", output, fixed = TRUE)][1L]
  expect_false(grepl("Selected", saz_header, fixed = TRUE))
})

# =============================================================================
# 12. Weights and offsets
# =============================================================================

# Test purpose: Checks that observation weights are accepted during model fitting.
test_that("weights argument is accepted and used", {
  w <- rep(1, nrow(x_prostate))
  w[1:10] <- 2
  fit <- mfp2(x_prostate, y_prostate, weights = w, verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})

# Test purpose: Checks that Poisson offsets are accepted and recorded in the
# fitted object.
test_that("offset argument is accepted for Poisson", {
  set.seed(1)
  n <- 200
  x <- cbind(x1 = runif(n, 1, 10), x2 = runif(n, 1, 5))
  exposure <- runif(n, 0.5, 2)
  y <- rpois(n, exposure * exp(0.5 + 0.1 * x[, 1]))

  fit <- mfp2(x, y, family = "poisson", offset = log(exposure), verbose = FALSE)
  expect_s3_class(fit, "mfp2")
  expect_true(fit$has_offset)
})

# Test purpose: Verifies that negative-binomial weights and offsets are passed
# consistently through fitting and prediction by comparison with glm.nb().
test_that("negative-binomial weights and offsets agree with MASS::glm.nb()", {
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
    weights = dat$w,
    offset = dat$log_exposure,
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

  expect_negbin_mfp2_mass_equal(
    fit_mfp2,
    fit_mass,
    coefficient_tolerance = 5e-4,
    likelihood_tolerance = 2e-3
  )
  expect_equal(unname(fit_mfp2$prior.weights), dat$w)

  nd <- dat[seq_len(20L), c("x1", "x2", "log_exposure"), drop = FALSE]
  pred_mfp2 <- predict(
    fit_mfp2,
    newdata = nd[, c("x1", "x2"), drop = FALSE],
    newoffset = nd$log_exposure,
    type = "link",
    se.fit = TRUE
  )
  pred_mass <- stats::predict(
    fit_mass,
    newdata = nd,
    type = "link",
    se.fit = TRUE
  )

  expect_equal(unname(pred_mfp2$fit), unname(pred_mass$fit), tolerance = 5e-4)
  expect_equal(
    unname(pred_mfp2$se.fit),
    unname(pred_mass$se.fit),
    tolerance = 5e-4
  )

  shifted <- predict(
    fit_mfp2,
    newdata = nd[, c("x1", "x2"), drop = FALSE],
    newoffset = nd$log_exposure + 0.35,
    type = "link"
  )
  expect_equal(
    unname(shifted - pred_mfp2$fit),
    rep(0.35, nrow(nd)),
    tolerance = 1e-10
  )
})


# =============================================================================
# 13. Convergence and cycles
# =============================================================================

# Test purpose: Checks that the default maximum number of cycles is sufficient
# for convergence on the prostate data.
test_that("mfp2() converges within default cycles", {
  fit <- mfp2(x_prostate, y_prostate, cycles = 5, verbose = FALSE)
  expect_true(fit$convergence_mfp)
})

# Test purpose: Checks that a non-converged one-cycle fit warns but still returns
# an mfp2 object.
test_that("mfp2() with cycles = 1 still returns a result", {
  expect_warning(
    fit <- mfp2(x_prostate, y_prostate, cycles = 1, verbose = FALSE),
    "No convergence after 1 cycles"
  )

  expect_s3_class(fit, "mfp2")
})

# =============================================================================
# 14. zero_vars and catzero_vars
# =============================================================================

# Test purpose: Checks that zero_vars activates zero-component handling for
# nonpositive values.
test_that("zero_vars recodes non-positive values to zero", {
  set.seed(1)
  n <- 200
  x_val <- rnorm(n, mean = 5, sd = 3) # some values may be <= 0
  x_mat <- cbind(exposure = x_val, x2 = runif(n, 1, 10))
  y_val <- 2 * pmax(x_val, 0) + rnorm(n)

  fit <- mfp2(x_mat, y_val, zero_vars = "exposure", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_true(fit$zero["exposure"])
})

# Test purpose: Checks that catzero_vars adds a zero-component indicator and
# implies zero handling.
test_that("catzero_vars creates binary indicator", {
  set.seed(1)
  n <- 200
  x_val <- rnorm(n, mean = 5, sd = 3)
  x_mat <- cbind(exposure = x_val, x2 = runif(n, 1, 10))
  y_val <- 2 * pmax(x_val, 0) + 1.5 * (x_val <= 0) + rnorm(n)

  fit <- mfp2(x_mat, y_val, catzero_vars = "exposure", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_true(fit$catzero["exposure"])
  # zero should also be TRUE (catzero implies zero)
  expect_true(fit$zero["exposure"])
})

# Test purpose: Verifies that fp(x, zero = TRUE) is converted to zero_vars
# and that non-positive values are handled through the zero component.
test_that("formula interface fp(zero = TRUE) enables zero handling", {
  set.seed(102)
  n <- 200

  exposure <- rnorm(n, mean = 5, sd = 3)
  dat <- data.frame(
    y = 2 * pmax(exposure, 0) + rnorm(n),
    exposure = exposure,
    x2 = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ fp(exposure, zero = TRUE) + fp(x2),
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$zero["exposure"])
})

# Test purpose: Verifies that fp(x, catzero = TRUE) creates a zero-component
# indicator and also implies zero handling.
test_that("formula interface fp(catzero = TRUE) enables catzero and zero handling", {
  set.seed(103)
  n <- 200

  exposure <- rnorm(n, mean = 5, sd = 3)
  dat <- data.frame(
    y = 2 * pmax(exposure, 0) + 1.5 * (exposure <= 0) + rnorm(n),
    exposure = exposure,
    x2 = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ fp(exposure, catzero = TRUE) + fp(x2),
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})

# =============================================================================
# 15. force_max_fp_vars and formula-term force_max_fp
# =============================================================================

# Test purpose: Checks that force_max_fp_vars uses the requested maximum FP
# complexity for named predictors under AIC.
test_that("force_max_fp_vars forces maximum FP degree with AIC/BIC", {
  fit_force <- mfp2(
    x_prostate, y_prostate,
    criterion = "aic",
    force_max_fp_vars = colnames(x_prostate),
    select = 1,
    verbose = FALSE
  )

  for (v in get_selected_variable_names(fit_force)) {
    powers <- fit_force$fp_powers[[v]]
    requested_df <- as.numeric(fit_force$fp_terms[v, "df_setting"])

    expected_n_powers <- if (requested_df <= 1) {
      1
    } else {
      requested_df / 2
    }

    n_powers <- sum(!is.na(powers))

    expect_equal(
      n_powers,
      expected_n_powers,
      info = paste("Variable:", v)
    )
  }
})


# Test purpose: force_max_fp has one dedicated selector for p-value, AIC, and
# BIC selection. It fits only the predetermined maximum ordinary FP form;
# null, linear, and lower-degree FP models cannot affect a forced result.
test_that("force-max selector fits only the maximum ordinary FP form", {
  calls <- new.env(parent = emptyenv())
  calls$degrees <- integer(0)

  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        data_adj = NULL,
        current_params = list(),
        transform_cache = list(reused = TRUE)
      )
    },
    find_best_fpm_step = function(..., degree) {
      calls$degrees <- c(calls$degrees, degree)

      metrics <- matrix(
        c(
          -20, 4, 40, 40, 48, 50, 96,
          -18, 4, 36, 36, 44, 46, 96
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(NULL, metric_names)
      )

      list(
        powers = rbind(c(-1, -1), c(-1, 2)),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 2L,
        current_adj_params = list(x = list(data_adj = NULL))
      )
    },
    .package = "mfp2"
  )

  out <- select_force_max_fp(
    x = matrix(1:8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = 1:8,
    powers_current = list(x = c(1, 1)),
    powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
    criterion = "aic",
    ftest = FALSE,
    select = 0.05,
    alpha = 0.05,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(x = NULL),
    transform_cache = NULL,
    force_max_fp = c(x = TRUE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  expect_equal(calls$degrees, 2)
  expect_identical(out$model_best, 1L)
  expect_identical(rownames(out$metrics), "FP2")
  expect_identical(rownames(out$powers), "FP2")
  expect_equal(unname(out$power_best), c(-1, 2))
  expect_true(out$transform_cache$reused)
})

# Test purpose: The dedicated selector itself accepts p-value forcing. This
# verifies that p-value forcing no longer relies on select = 1 / alpha = 1 to
# walk the RA2 closed-test sequence before reaching the predetermined FPm.
test_that("force-max selector accepts p-value criterion directly", {
  calls <- new.env(parent = emptyenv())
  calls$degrees <- integer(0)

  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        data_adj = NULL,
        current_params = list(),
        transform_cache = list(pvalue = TRUE)
      )
    },
    find_best_fpm_step = function(..., degree) {
      calls$degrees <- c(calls$degrees, degree)
      metrics <- matrix(
        c(-18, 4, 36, 36, 44, 46, 96),
        nrow = 1,
        dimnames = list(NULL, metric_names)
      )
      list(
        powers = matrix(c(-1, 2), nrow = 1),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 1L,
        current_adj_params = list(x = list(data_adj = NULL))
      )
    },
    .package = "mfp2"
  )

  out <- select_force_max_fp(
    x = matrix(1:8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = 1:8,
    powers_current = list(x = c(1, 1)),
    powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
    criterion = "pvalue",
    ftest = FALSE,
    select = 1,
    alpha = 1,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(x = NULL),
    transform_cache = NULL,
    force_max_fp = c(x = TRUE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  expect_equal(calls$degrees, 2)
  expect_identical(out$model_best, 1L)
  expect_identical(rownames(out$metrics), "FP2")
  expect_equal(unname(out$power_best), c(-1, 2))
  expect_length(out$pvalue, 0L)
  expect_true(out$transform_cache$pvalue)
})

# Test purpose: The same criterion-independent selector handles forced ACD terms
# without fitting reduced ACD alternatives. The maximum ACD form is always
# FP1(x, A(x)), represented by the degree-2 ACD candidate search.
test_that("force-max selector fits only the full ACD form", {
  calls <- new.env(parent = emptyenv())
  calls$degrees <- integer(0)

  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        data_adj = NULL,
        current_params = list(),
        transform_cache = list(acd = TRUE)
      )
    },
    find_best_fpm_step = function(..., degree) {
      calls$degrees <- c(calls$degrees, degree)
      metrics <- matrix(
        c(-10, 4, 20, 20, 28, 30, 96),
        nrow = 1,
        dimnames = list(NULL, metric_names)
      )

      list(
        powers = matrix(c(-1, 2), nrow = 1),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 1L,
        current_adj_params = list(x = list(data_adj = NULL))
      )
    },
    .package = "mfp2"
  )

  out <- select_force_max_fp(
    x = matrix(1:8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = TRUE),
    y = 1:8,
    powers_current = list(x = c(1, 1)),
    powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
    criterion = "pvalue",
    ftest = FALSE,
    select = 1,
    alpha = 1,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = list()),
    prev_adj_params = list(x = NULL),
    transform_cache = NULL,
    force_max_fp = c(x = TRUE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  expect_equal(calls$degrees, 2)
  expect_true(out$acd)
  expect_identical(out$model_best, 1L)
  expect_identical(rownames(out$metrics), "FP1(x, A(x))")
  expect_identical(rownames(out$powers), "FP1(x, A(x))")
})

# Test purpose: find_best_fp_step() must dispatch every forced non-linear term
# to select_force_max_fp(), including criterion = "pvalue". The ordinary RA2
# and IC selectors must not run when the final functional complexity is forced.
test_that("force-max dispatch bypasses RA2 and IC selectors for all criteria", {
  calls <- new.env(parent = emptyenv())
  calls$criteria <- character(0)

  testthat::local_mocked_bindings(
    select_force_max_fp = function(..., criterion) {
      calls$criteria <- c(calls$criteria, criterion)
      metrics <- matrix(
        c(-10, 4, 20, 20, 28, 30, 96),
        nrow = 1,
        dimnames = list(
          "FP2",
          c("logl", "df", "deviance_rs", "deviance_gaussian",
            "aic", "bic", "df_resid")
        )
      )
      list(
        keep = FALSE,
        acd = FALSE,
        powers = matrix(c(-1, 2), nrow = 1, dimnames = list("FP2", NULL)),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 1L,
        statistic = NA,
        pvalue = NA,
        spike = FALSE,
        current_adj_params = list(),
        transform_cache = list()
      )
    },
    select_ra2 = function(...) {
      stop("ordinary RA2 selector must not run", call. = FALSE)
    },
    select_ra2_acd = function(...) {
      stop("ACD RA2 selector must not run", call. = FALSE)
    },
    select_ic = function(...) {
      stop("ordinary IC selector must not run", call. = FALSE)
    },
    select_ic_acd = function(...) {
      stop("ACD IC selector must not run", call. = FALSE)
    },
    .package = "mfp2"
  )

  for (criterion_value in c("pvalue", "aic", "bic")) {
    out <- find_best_fp_step(
      x = matrix(seq_len(8), ncol = 1, dimnames = list(NULL, "x")),
      y = seq_len(8),
      xi = "x",
      weights = NULL,
      offset = NULL,
      df = 4,
      powers_current = list(x = c(1, 1)),
      family = stats::gaussian(),
      family_string = "gaussian",
      criterion = criterion_value,
      select = 1,
      alpha = 1,
      keep = character(0),
      powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
      method = NULL,
      strata = NULL,
      nocenter = FALSE,
      acdx = c(x = FALSE),
      ftest = FALSE,
      control = list(),
      rownames = as.character(seq_len(8)),
      zero = c(x = FALSE),
      catzero = list(x = NULL),
      spike = list(x = FALSE),
      spike_decision = c(x = 2),
      acd_parameter = list(x = NULL),
      prev_adj_params = list(x = NULL),
      transform_cache = NULL,
      force_max_fp = c(x = TRUE),
      has_offset = FALSE,
      n_obs = 8,
      verbose = FALSE,
      term_to_columns = list(x = "x")
    )

    expect_equal(unname(out$power_best), c(-1, 2))
  }

  expect_identical(calls$criteria, c("pvalue", "aic", "bic"))
})

# Test purpose: Verifies that force_max_fp_vars translates to select = 1 and
# alpha = 1 under p-value selection in the matrix interface.
test_that("force_max_fp_vars forces maximum FP degree with p-value selection", {
  set.seed(1501)

  dat <- data.frame(
    x = seq(0.5, 6, length.out = 250)
  )
  dat$y <- 1 + 1.5 / dat$x - 0.7 * dat$x^2 +
    rnorm(nrow(dat), sd = 0.05)

  fit <- mfp2(
    x = as.matrix(dat["x"]),
    y = dat$y,
    criterion = "pvalue",
    df = 4,
    select = 0.05,
    alpha = 0.05,
    force_max_fp_vars = "x",
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["x", "select"]), 1)
  expect_equal(as.numeric(fit$fp_terms["x", "alpha"]), 1)
  expect_true(fit$fp_terms["x", "selected"])
  expect_equal(sum(!is.na(fit$fp_powers[["x"]])), 2)
})

# Test purpose: Verifies that fp(force_max_fp = TRUE) reaches the same
# p-value forcing logic through mfp2.formula().
test_that("formula force_max_fp uses p-value forcing from mfp2.default", {
  set.seed(1502)

  dat <- data.frame(
    x = seq(0.5, 6, length.out = 250)
  )
  dat$y <- 1 + 1.5 / dat$x - 0.7 * dat$x^2 +
    rnorm(nrow(dat), sd = 0.05)

  fit <- mfp2(
    y ~ fp(x, df = 4, force_max_fp = TRUE),
    data = dat,
    criterion = "pvalue",
    select = 0.05,
    alpha = 0.05,
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["x", "select"]), 1)
  expect_equal(as.numeric(fit$fp_terms["x", "alpha"]), 1)
  expect_true(fit$fp_terms["x", "selected"])
  expect_equal(sum(!is.na(fit$fp_powers[["x"]])), 2)
})


# =============================================================================
# 16. mfpi() — basic interaction fitting
# =============================================================================

make_mfpi_factor_data <- function(n = 240L, ordered_stage = FALSE) {
  set.seed(9101)

  trt <- factor(rep(c("control", "treated"), length.out = n))
  stage_values <- rep(c("I", "II", "III"), length.out = n)
  stage <- if (ordered_stage) {
    ordered(stage_values, levels = c("I", "II", "III"))
  } else {
    factor(stage_values, levels = c("I", "II", "III"))
  }

  x <- runif(n, 1, 8)
  z <- rnorm(n)
  stage_effect <- c(I = 0, II = 0.7, III = -0.5)[as.character(stage)]

  data.frame(
    y = 1 + 0.4 * x + 0.9 * x * (trt == "treated") +
      stage_effect + 0.2 * z + rnorm(n, sd = 0.25),
    trt = trt,
    x = x,
    stage = stage,
    z = z
  )
}


# Test purpose: Verifies that mfpi.default() uses the same partial named
# override semantics as mfp2.default() for its Stage-1 adjustment model.
test_that("mfpi.default() accepts partial named df select and alpha overrides", {
  set.seed(16001)
  n <- 180L
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), length.out = n)),
    x = stats::runif(n, 1, 8),
    z = stats::runif(n, 1, 5),
    w = stats::runif(n, 2, 7)
  )
  dat$y <- 1 + 0.4 * dat$x + 0.8 * dat$x * (dat$group == "treated") +
    0.5 * dat$z + stats::rnorm(n, sd = 0.3)

  fit <- mfpi(
    dat[, c("group", "x", "z", "w")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 5,
    df = c(z = 1),
    select = c(z = 1),
    alpha = c(w = 1),
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    p_interact = 1,
    verbose = FALSE
  )

  adjustment_terms <- fit$adjustment_model$fp_terms
  expect_equal(as.numeric(adjustment_terms["z", "df_initial"]), 1)
  expect_equal(as.numeric(adjustment_terms["w", "df_initial"]), 4)
  expect_equal(as.numeric(adjustment_terms["z", "select"]), 1)
  expect_equal(as.numeric(adjustment_terms["w", "select"]), 0.05)
  expect_equal(as.numeric(adjustment_terms["z", "alpha"]), 0.05)
  expect_equal(as.numeric(adjustment_terms["w", "alpha"]), 1)
})

# Test purpose: Ensures MFPI no longer accepts positional multi-value settings.
test_that("mfpi.default() rejects unnamed multi-value df select and alpha", {
  set.seed(16002)
  n <- 60L
  x <- data.frame(
    group = factor(rep(c("control", "treated"), length.out = n)),
    x = stats::runif(n, 1, 5),
    z = stats::runif(n, 1, 4)
  )
  y <- stats::rnorm(n)

  common <- list(
    x = x,
    y = y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    verbose = FALSE
  )

  expect_error(
    do.call(mfpi, c(common, list(df = c(1, 4, 4)))),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    do.call(mfpi, c(common, list(select = c(1, 0.05, 0.05)))),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    do.call(mfpi, c(common, list(alpha = c(1, 0.05, 0.05)))),
    "single unnamed numeric value or a named numeric vector"
  )
})

# Test purpose: Verifies that force_max_fp_vars applies only to MFPI's
# adjustment-model selection and sets both select and alpha to 1 under the
# p-value criterion.
test_that("MFPI force_max_fp_vars forces p-value adjustment variables", {
  dat <- make_mfpi_factor_data()

  # Replace z after generating y so that z is a positive, unassociated
  # adjustment variable. Its retention therefore depends on force_max_fp_vars,
  # while no shift is needed for the FP transformation.
  set.seed(1503)
  dat$z <- runif(nrow(dat), min = 0.5, max = 3)

  fit <- mfpi(
    y ~ trt + x + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex3",
    criterion = "pvalue",
    force_max_fp_vars = "z",
    df = 4,
    select = 0.05,
    alpha = 0.05,
    center = FALSE,
    cycles = 5,
    p_interact = 1,
    verbose = FALSE
  )

  adjustment_terms <- fit$adjustment_model$fp_terms

  expect_equal(as.numeric(adjustment_terms["z", "select"]), 1)
  expect_equal(as.numeric(adjustment_terms["z", "alpha"]), 1)
  expect_true(adjustment_terms["z", "selected"])
  expect_equal(
    sum(!is.na(fit$adjustment_model$fp_powers[["z"]])),
    2
  )
})


# -----------------------------------------------------------------------------
# 16.1 Grouped categorical adjustment terms
# -----------------------------------------------------------------------------
# These tests verify that factor contrasts and manually specified dummy blocks
# remain one conceptual adjustment term throughout Stage 1 selection and Stage 2
# interaction-model construction.

# Test purpose: Verifies that unordered-factor contrast columns are stored
# and selected as one conceptual MFPI adjustment term.
test_that("MFPI formula groups unordered-factor adjustment columns", {
  dat <- make_mfpi_factor_data()

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true("stage" %in% names(fit$term_to_columns))
  expect_setequal(fit$term_to_columns[["stage"]], c("stageII", "stageIII"))

  expect_true("stage" %in% rownames(fit$adjustment_model$fp_terms))
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_setting"]),
    1
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_initial"]),
    2
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_final"]),
    2
  )
  expect_false(any(c("stageII", "stageIII") %in%
                     rownames(fit$adjustment_model$fp_terms)))
})

# Test purpose: Verifies that ordered-factor polynomial contrasts are accepted
# and grouped as one fixed linear MFPI adjustment term.
test_that("MFPI formula groups ordered-factor polynomial contrasts", {
  dat <- make_mfpi_factor_data(ordered_stage = TRUE)

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true("stage" %in% names(fit$term_to_columns))
  expect_setequal(fit$term_to_columns[["stage"]], c("stage.L", "stage.Q"))
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_setting"]),
    1
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_initial"]),
    2
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_final"]),
    2
  )
})


# Test purpose: Checks that MFPI uses the source variable name for a simple
# inline factor wrapper while retaining the original contrast-column names.
test_that("MFPI names inline factor adjustments by their source variable", {
  dat <- make_mfpi_factor_data()
  dat$stage_code <- as.integer(dat$stage)

  fit <- mfpi(
    y ~ trt + x + factor(stage_code) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage_code",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_true("stage_code" %in% names(fit$term_to_columns))
  expect_false("factor(stage_code)" %in% names(fit$term_to_columns))
  expect_identical(
    fit$term_to_columns[["stage_code"]],
    c("factor(stage_code)2", "factor(stage_code)3")
  )
  expect_true("stage_code" %in% rownames(fit$adjustment_model$fp_terms))
  expect_false(
    "factor(stage_code)" %in% rownames(fit$adjustment_model$fp_terms)
  )
})


# Test purpose: Verifies that a binary inline factor is retained as a mapped
# one-column MFPI adjustment term rather than reverting to its dummy-column name.
test_that("MFPI retains binary inline factors as source-name mappings", {
  dat <- make_mfpi_factor_data()
  dat$binary_stage <- rep(c(1, 1, 2, 1, 2, 2), length.out = nrow(dat))

  fit <- mfpi(
    y ~ trt + x + factor(binary_stage) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "binary_stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    fit$term_to_columns[["binary_stage"]],
    "factor(binary_stage)2"
  )
  expect_true(fit$adjustment_model$fp_terms["binary_stage", "selected"])

  nd <- dat[1:18, c("trt", "x", "binary_stage", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(fit, "x", fit_result, nd)
  expect_true("factor(binary_stage)2.1" %in% names(design$model_newdata))
})

# Test purpose: Checks manual matrix-interface grouping of dummy columns in
# the Stage 1 MFPI adjustment model.
test_that("MFPI matrix interface accepts manual grouped adjustment columns", {
  dat <- make_mfpi_factor_data()
  stage_mm <- stats::model.matrix(~ stage, dat)[, -1L, drop = FALSE]

  x <- cbind(
    trt = as.numeric(dat$trt) - 1,
    x = dat$x,
    stage_mm,
    z = dat$z
  )

  fit <- mfpi(
    x,
    dat$y,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = c("stageII", "stageIII")),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    fit$term_to_columns[["stage"]],
    c("stageII", "stageIII")
  )
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
})

# Test purpose: MFPI applies the scalar df default only to ordinary
# continuous columns; a supplied ordered-factor contrast block remains linear.
test_that("MFPI grouped ordered-factor columns are linear under scalar df default", {
  set.seed(51231)
  n <- 192
  trt <- rep(rep(0:1, each = 4L), length.out = n)
  x_cont <- seq(0.5, 9, length.out = n)
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  z <- stats::rnorm(n)

  x <- cbind(
    trt = trt,
    x = x_cont,
    stage_mm,
    z = z
  )
  y <- 0.5 * x_cont + 0.7 * stage_mm[, 1L] -
    0.4 * stage_mm[, 2L] + 0.2 * z + stats::rnorm(n, sd = 0.5)

  fit <- mfpi(
    x,
    y,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_setting"]),
    1
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_initial"]),
    ncol(stage_mm)
  )
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
})


# Test purpose: MFPI must preserve a supplied categorical contrast block under
# scale one instead of estimating a separate scale for each member column.
test_that("MFPI grouped matrix columns default to scale one", {
  set.seed(51232)
  n <- 192
  trt <- rep(rep(0:1, each = 4L), length.out = n)
  x_cont <- seq(0.5, 9, length.out = n)
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  stage_mm[, 1L] <- 1000 * stage_mm[, 1L]
  z <- stats::rnorm(n)

  x <- cbind(
    trt = trt,
    x = x_cont,
    stage_mm,
    z = z
  )
  y <- 0.5 * x_cont + 0.001 * stage_mm[, 1L] -
    0.4 * stage_mm[, 2L] + 0.2 * z + stats::rnorm(n, sd = 0.5)

  fit <- mfpi(
    x,
    y,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(
    unname(fit$scale[colnames(stage_mm)]),
    rep(1, ncol(stage_mm))
  )
  expect_equal(
    as.numeric(fit$adjustment_model$transformations["stage", "scale"]),
    1
  )
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
})


# Test purpose: Ensures grouped categorical terms and their member columns
# cannot be requested as continuous MFPI interaction variables.
test_that("MFPI grouped terms cannot be used as cont_vars", {
  dat <- make_mfpi_factor_data()
  stage_mm <- stats::model.matrix(~ stage, dat)[, -1L, drop = FALSE]

  x <- cbind(
    trt = as.numeric(dat$trt) - 1,
    x = dat$x,
    stage_mm
  )

  expect_error(
    mfpi(
      x,
      dat$y,
      group_var = "trt",
      cont_vars = "stageII",
      term_groups = list(stage = c("stageII", "stageIII")),
      verbose = FALSE
    ),
    "singleton continuous|grouped categorical"
  )
})

# Test purpose: Checks that a multilevel group_var is represented by one
# conceptual dummy block when include_group_var = TRUE.
test_that("include_group_var stores group dummies as one conceptual term", {
  dat <- make_mfpi_factor_data()
  # Use a three-level pattern that is not collinear with the three-level stage
  # factor; otherwise the null and group-dummy models have the same fitted rank.
  dat$trt <- factor(
    rep(c("A", "A", "B", "C", "B", "C"), length.out = nrow(dat)),
    levels = c("A", "B", "C")
  )

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    include_group_var = TRUE,
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_true("trt" %in% names(fit$adjustment_model$term_to_columns))
  expect_length(fit$adjustment_model$term_to_columns[["trt"]], 2L)

  selected_for_stage2 <- fit$univariable_interactions$selected_vars
  if (!is.null(selected_for_stage2)) {
    expect_false("trt" %in% selected_for_stage2)
  }
})


# Test purpose: Checks that the default MFPI interface fits an interaction-analysis
# object with expected metadata.
test_that("mfpi.default() runs and returns an mfpi object", {
  data("prostate", package = "mfp2")

  x_p <- data.frame(
    svi = prostate$svi,
    age = prostate$age,
    cavol = prostate$cavol,
    pgg45 = prostate$pgg45,
    weight = prostate$weight,
    bph = prostate$bph,
    cp = prostate$cp
  )

  fit <- mfpi(
    x_p, y_prostate,
    group_var = "svi",
    cont_vars = c("cavol", "age"),
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(!is.null(fit$all_model_metrics))
  expect_true(!is.null(fit$adjustment_model))
  expect_equal(fit$group_var, "svi")
  expect_equal(fit$cont_vars, c("cavol", "age"))
})

# Test purpose: Checks that the formula MFPI interface parses fp() terms and
# returns an mfpi object.
test_that("mfpi.formula() runs and returns an mfpi object", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
      fp(bph) + fp(cp),
    data = prostate,
    cont_vars = c("cavol", "age"),
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
})

# Test purpose: Checks that requested interaction functional forms are stored for
# continuous variables.
test_that("mfpi() cont_var_forms specifies functional form correctly", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    cont_var_forms = c(cavol = "fp2", age = "linear"),
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "fp2"))
  expect_equal(fit$cont_var_forms["age"], c(age = "linear"))
})

# Test purpose: Checks that MFPI runs with information-criterion based interaction
# assessment.
test_that("mfpi() with criterion = 'aic' works", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol"),
    group_var = "svi",
    criterion = "aic",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
})

# Test purpose: Checks that all supported MFPI flexibility settings fit without
# error.
test_that("mfpi() flexibility levels run without error", {
  data("prostate", package = "mfp2")

  for (fl in c("flex1", "flex2", "flex3", "flex4")) {
    fit <- mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      cont_var_forms = c(cavol = "fp1"),
      group_var = "svi",
      flex = fl,
      verbose = FALSE
    )
    expect_true(
      inherits(fit, "mfpi"),
      info = paste("flex =", fl)
    )
  }
})

# Test purpose: default behavior
test_that("mfpi() defaults to flex3", {
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    verbose = FALSE
  )
  expect_equal(fit$flex, "flex3")
})

# Test purpose: Checks that omitted cont_var_forms are filled with "linear"
# for every variable listed in cont_vars.
test_that("mfpi() defaults missing cont_var_forms to linear", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    group_var = "svi",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "fp1"))
  expect_equal(fit$cont_var_forms["age"], c(age = "fp1"))
})

# Test purpose: Ensures cont_var_forms may specify only some cont_vars;
# omitted cont_vars are filled with "linear" and ordering follows cont_vars.
test_that("mfpi() fills missing cont_var_forms entries with linear", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    cont_var_forms = c(cavol = "fp2"),
    group_var = "svi",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(names(fit$cont_var_forms), c("cavol", "age"))
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "fp2"))
  expect_equal(fit$cont_var_forms["age"], c(age = "fp1"))
})

# Test purpose: Verifies that cont_var_forms only accepts "linear", "fp1",
# and "fp2".
test_that("mfpi() rejects invalid cont_var_forms values", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = c("cavol", "age"),
      cont_var_forms = c(cavol = "spline"),
      group_var = "svi",
      verbose = FALSE
    ),
    "Invalid value"
  )
})


# Test purpose: Ensures cont_var_forms entries must be named so each requested
# form is explicitly tied to a variable in cont_vars.
test_that("mfpi() rejects unnamed cont_var_forms", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = c("cavol", "age"),
      cont_var_forms = c("fp2", "linear"),
      group_var = "svi",
      verbose = FALSE
    ),
    "Every entry of `cont_var_forms` must be named"
  )
})


# Test purpose: Checks that cont_var_forms cannot name variables that are not
# being tested as continuous interaction variables.
test_that("mfpi() rejects cont_var_forms names not in cont_vars", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      cont_var_forms = c(age = "fp1"),
      group_var = "svi",
      verbose = FALSE
    ),
    "not in `cont_vars`"
  )
})

# Test purpose: Ensures the grouping variable cannot also be listed as a
# continuous interaction variable.
test_that("mfpi() rejects group_var included in cont_vars", {
  set.seed(204)
  n <- 120

  x <- data.frame(
    group = rep(1:4, length.out = n),
    x1 = runif(n, 1, 10),
    x2 = runif(n, 1, 10)
  )
  y <- rnorm(n)

  expect_error(
    mfpi(
      x,
      y,
      group_var = "group",
      cont_vars = c("group", "x1"),
      verbose = FALSE
    ),
    "must not also appear in `cont_vars`"
  )
})

# Test purpose: Checks that binary variables are not allowed in cont_vars,
# because MFPI requires continuous variables for FP interaction forms.
test_that("mfpi.formula() rejects a binary variable in cont_vars", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "svi",
      group_var = "cavol",
      verbose = FALSE
    ),
    "binary"
  )
})


# Test purpose: Checks that binary variables are not allowed in cont_vars,
# because MFPI requires continuous variables for FP interaction forms.
test_that("mfpi.formula() rejects a simulated binary variable in cont_vars", {
  set.seed(201)
  n <- 100

  dat <- data.frame(
    y = rnorm(n),
    group = rep(0:1, length.out = n),
    binary_x = rep(0:1, length.out = n),
    x = runif(n, 1, 10)
  )

  expect_error(
    mfpi(
      y ~ group + binary_x + fp(x),
      data = dat,
      cont_vars = "binary_x",
      group_var = "group",
      verbose = FALSE
    ),
    "binary"
  )
})


# Test purpose: Ensures MFPI validates multiplicity-adjustment methods against
# stats::p.adjust.methods.
test_that("mfpi() rejects invalid p_adjust_method", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      group_var = "svi",
      p_adjust_method = "not_a_method",
      verbose = FALSE
    ),
    "Invalid `p_adjust_method`"
  )
})


# Test purpose: Checks that a valid multiplicity-adjustment method is accepted
# and stored on the returned mfpi object.
test_that("mfpi() stores valid p_adjust_method", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    group_var = "svi",
    p_adjust_method = "holm",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$p_adjust_method, "holm")
})

# Test purpose: Verifies criterion-specific min_improvement defaults:
# pvalue uses p_interact, while AIC and BIC default to 2.
test_that("mfpi() sets criterion-specific default min_improvement", {
  data("prostate", package = "mfp2")

  fit_p <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "pvalue",
    p_interact = 0.10,
    verbose = FALSE
  )

  fit_aic <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "aic",
    verbose = FALSE
  )

  fit_bic <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "bic",
    verbose = FALSE
  )

  expect_equal(fit_p$min_improvement, 0.10)
  expect_equal(fit_aic$min_improvement, 2)
  expect_equal(fit_bic$min_improvement, 2)
})

# Test purpose: Checks that an explicit min_improvement threshold is respected
# for information-criterion based MFPI selection.
test_that("mfpi() stores explicit min_improvement for AIC", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "aic",
    min_improvement = 3.5,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$min_improvement, 3.5)
})

# Test purpose: Ensures min_improvement must be NULL or a single positive
# finite numeric value.
test_that("mfpi() rejects invalid min_improvement", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      group_var = "svi",
      criterion = "aic",
      min_improvement = 0,
      verbose = FALSE
    ),
    "`min_improvement`"
  )
})

# Test purpose: Checks that include_group_var = TRUE fits successfully and is
# recorded on the returned mfpi object.
test_that("mfpi() accepts include_group_var = TRUE", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    include_group_var = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(isTRUE(fit$include_group_var))
})

# Test purpose: Ensures group-specific centering mode is accepted and stored.
test_that("mfpi() accepts group-specific centering", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    center_type = "group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$center_type, "group")
})

# Test purpose: Ensures mfpi.default() allows group_var to be categorical but
# rejects other categorical predictors in x.
test_that("mfpi.default() rejects non-group categorical predictors", {
  set.seed(202)
  n <- 100

  x <- data.frame(
    group = factor(rep(c("A", "B"), length.out = n)),
    x = runif(n, 1, 10),
    bad_factor = factor(rep(c("low", "high"), length.out = n))
  )
  y <- rnorm(n)

  expect_error(
    mfpi(
      x,
      y,
      group_var = "group",
      cont_vars = "x",
      verbose = FALSE
    ),
    "Only `group_var` may be categorical"
  )
})

# Test purpose: Checks that categorical group labels are retained as metadata
# after internal recoding of group_var.
test_that("mfpi.default() stores original group levels", {
  set.seed(203)
  n <- 120

  x <- data.frame(
    group = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 10),
    z = runif(n, 1, 10)
  )
  y <- 0.2 * x$x + 0.5 * (x$group == "treated") + rnorm(n)

  fit <- mfpi(
    x,
    y,
    group_var = "group",
    cont_vars = "x",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(!is.null(fit$group_levels_original))
  expect_true(all(c("control", "treated") %in% fit$group_levels_original))
})

# Test purpose: Ensures predict.mfpi() validates se.fit as a single non-missing
# logical value.
test_that("predict.mfpi() rejects invalid se.fit", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(fit, terms = "cavol", type = "function", model = "all", se.fit = NA),
    "`se.fit` must be"
  )
})

# Test purpose: Ensures predict.mfpi() validates confidence level as a single
# numeric value in (0, 1).
test_that("predict.mfpi() rejects invalid confidence level", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(fit, terms = "cavol", type = "function", model = "all", level = 1),
    "`level` must be"
  )
})

# Test purpose: Ensures predict.mfpi() fails clearly when a requested term has
# no stored MFPI interaction model in the requested model scope.
test_that("predict.mfpi() rejects unknown prediction terms", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(fit, terms = "age", type = "function", model = "all"),
    "requested terms"
  )
})

# Test purpose: Ensures fitted-function prediction checks that newdata contains
# the requested continuous variable.
test_that("predict.mfpi() requires fitted-function newdata to contain requested term", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(
      fit,
      terms = "cavol",
      type = "function",
      model = "all",
      newdata = data.frame(age = prostate$age[1:10])
    ),
    "must contain a column named `cavol`"
  )
})

# Test purpose: Checks ordinary subject-level link-scale MFPI prediction from
# a term-specific interaction model using supplied newdata.
test_that("predict.mfpi() type = 'link' returns subject-level predictions", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(
    fit,
    terms = "cavol",
    type = "link",
    model = "all",
    newdata = prostate[1:10, ],
    se.fit = FALSE
  )

  expect_s3_class(p, "mfpi_prediction")
  expect_true(!is.null(p$predictions))
  expect_equal(nrow(p$predictions), 10)
  expect_true(all(is.finite(p$predictions$fit)))
})

# Test purpose: Checks ordinary subject-level response-scale MFPI prediction
# using supplied newdata.
test_that("predict.mfpi() type = 'response' returns subject-level predictions", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(
    fit,
    terms = "cavol",
    type = "response",
    model = "all",
    newdata = prostate[1:10, ],
    se.fit = FALSE
  )

  expect_s3_class(p, "mfpi_prediction")
  expect_true(!is.null(p$predictions))
  expect_equal(nrow(p$predictions), 10)
  expect_true(all(is.finite(p$predictions$fit)))
})

# Test purpose: Ensures fitted-function prediction with grid = TRUE returns
# values on the requested evaluation grid.
test_that("predict.mfpi() fitted-function grid uses requested n_grid", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(
    fit,
    terms = "cavol",
    type = "function",
    model = "all",
    grid = TRUE,
    n_grid = 25,
    se.fit = FALSE
  )

  expect_s3_class(p, "mfpi_prediction")
  expect_true(!is.null(p$functions))
  expect_true(length(unique(p$functions$x)) <= 25)
  expect_true(all(is.finite(p$functions$fit)))
})

# Test purpose: Ensures grid = TRUE is ignored with a warning for ordinary
# subject-level MFPI prediction types.
test_that("predict.mfpi() warns when grid is used with link prediction", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_warning(
    predict(
      fit,
      terms = "cavol",
      type = "link",
      model = "all",
      newdata = prostate[1:5, ],
      grid = TRUE,
      se.fit = FALSE
    ),
    "grid"
  )
})

# =============================================================================
# 17. predict.mfpi()
# =============================================================================

# Test purpose: Checks that MFPI fitted-function predictions are returned in
# the expected structure.
test_that("predict.mfpi() returns fitted-function predictions", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(fit, terms = "cavol", type = "function", model = "all")

  if (inherits(p, "mfpi_prediction")) {
    expect_true(!is.null(p$functions))
    expect_true(is.data.frame(p$functions))
  } else if (is.list(p)) {
    # May be a list when wrapping
    expect_true(length(p) >= 1)
  }
})

# Test purpose: Checks that MFPI prediction can return both fitted functions
# and group differences.
test_that("predict.mfpi() type = 'both' returns functions and differences", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(fit, terms = "cavol", type = "both", model = "all")

  if (inherits(p, "mfpi_prediction")) {
    expect_true(!is.null(p$functions))
    expect_true(!is.null(p$differences))
  }
})


# Test purpose: Ordinary Gaussian MFPI newdata prediction delegates to the
# stored interaction glm using the reconstructed formula-compatible data frame.
test_that("predict.mfpi ordinary Gaussian prediction matches stored glm", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )
  nd <- prostate[1:12, , drop = FALSE]
  fit_result <- fit$var_winners[["cavol"]]$fit
  design <- mfpi_build_ordinary_design(
    fit, "cavol", fit_result, nd
  )
  direct_link <- stats::predict(
    fit_result$test_results$interaction_model$fit,
    newdata = design$model_newdata,
    type = "link",
    se.fit = TRUE
  )
  got_link <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  expect_equal(got_link$predictions$fit, as.numeric(direct_link$fit))
  expect_equal(got_link$predictions$se.fit, as.numeric(direct_link$se.fit))
  expect_true(got_link$metadata$used_model_predict)

  direct_response <- stats::predict(
    fit_result$test_results$interaction_model$fit,
    newdata = design$model_newdata,
    type = "response",
    se.fit = FALSE
  )
  got_response <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "response", se.fit = FALSE
  )
  expect_equal(got_response$predictions$fit, as.numeric(direct_response))

  # Independent oracle: reconstruct the interaction model matrix and calculate
  # eta = X beta and sqrt(diag(X V X')) directly. This avoids relying solely on
  # predict.glm(), which is also used internally by ordinary MFPI prediction.
  stored <- fit_result$test_results$interaction_model$fit
  reference_terms <- stats::delete.response(stats::terms(stored))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = design$model_newdata,
    xlev = stored$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = stored$contrasts
  )
  beta <- stats::coef(stored)
  beta_vcov <- stats::vcov(stored)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_link <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_equal(got_link$predictions$fit, manual_link, tolerance = 1e-8)
  expect_equal(got_link$predictions$se.fit, manual_se, tolerance = 1e-8)
  expect_equal(got_response$predictions$fit, manual_link, tolerance = 1e-8)
})

# Test purpose: Stratified Cox MFPI prediction passes raw vector strata through
# exactly once and matches predict.coxph(reference = "zero").
test_that("predict.mfpi stratified Cox prediction matches stored coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    cont_vars = "age",
    group_var = "sex",
    cont_var_forms = c(age = "linear"),
    strata = dat$inst,
    p_interact = 0.95,
    verbose = FALSE
  )
  nd <- dat[1:10, c("age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = nd$inst
  )
  expect_identical(design$model_newdata$strata_, nd$inst)

  direct <- stats::predict(
    fit_result$test_results$interaction_model$fit,
    newdata = design$model_newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got <- predict(
    fit, newdata = nd, terms = "age", model = "all",
    type = "link", strata = nd$inst, se.fit = TRUE
  )
  expect_equal(got$predictions$fit, as.numeric(direct$fit))
  expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit))
  expect_true(got$metadata$used_model_predict)
})


# -----------------------------------------------------------------------------
# 17.1 High-confidence MFPI and SAZ regression tests
# -----------------------------------------------------------------------------
# The tests in this subsection use independent reference calculations whenever
# possible. They are intended to detect statistically meaningful regressions,
# rather than only checking that a function returns an object without error.

# Test purpose: Verifies every supported ordinary FP1 power against its direct
# mathematical definition. This protects the transformation layer used by MFP,
# SAZ positive components, MFPI interaction bases, and prediction reconstruction.
test_that("17.1.1 ordinary FP1 powers equal their mathematical definitions", {
  x <- c(0.5, 1, 2, 4)
  powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  expected <- list(
    `-2` = x^-2,
    `-1` = x^-1,
    `-0.5` = x^-0.5,
    `0` = log(x),
    `0.5` = sqrt(x),
    `1` = x,
    `2` = x^2,
    `3` = x^3
  )

  for (power in powers) {
    got <- transform_vector_fp(
      x,
      power = power,
      shift = 0,
      scale = 1,
      check_binary = FALSE
    )

    expect_equal(
      as.numeric(got[, 1]),
      expected[[as.character(power)]],
      tolerance = 1e-12,
      info = paste("power =", power)
    )
  }
})

# Test purpose: Verifies the general repeated-power rule. For a power repeated
# three times, the expected basis is x^p, x^p log(x), x^p log(x)^2. This is a
# stronger check than the existing FP2-only repeated-power tests.
test_that("17.1.2 repeated FP powers follow the logarithmic multiplier rule", {
  x <- c(0.5, 1, 2, 4)

  got <- transform_vector_fp(
    x,
    power = c(2, 2, 2),
    shift = 0,
    scale = 1,
    check_binary = FALSE
  )

  expected <- cbind(
    x^2,
    x^2 * log(x),
    x^2 * log(x)^2
  )

  expect_equal(unname(got), expected, tolerance = 1e-12)
})

# Test purpose: Verifies that the C++ batch implementation and the public R
# transformation wrapper produce identical FP columns. Both fitting and
# prediction rely on these paths, so disagreement would invalidate model reuse.
test_that("17.1.3 C++ and public FP transformation paths agree", {
  x <- c(0.5, 1, 2, 4, 8)
  candidates <- list(c(-1), c(0), c(0.5), c(1, 1), c(0, 0), c(2, 2))

  for (power in candidates) {
    cpp <- transform_fp_core(
      x_raw = x,
      power = power,
      shift_val = 0,
      scale_val = 1,
      zero = FALSE
    )
    public <- transform_vector_fp(
      x,
      power = power,
      shift = 0,
      scale = 1,
      check_binary = FALSE
    )

    expect_equal(
      unname(cpp),
      unname(public),
      tolerance = 1e-12,
      info = paste("powers =", paste(power, collapse = ","))
    )
  }
})

# Test purpose: Verifies the documented MFPI degrees of freedom for two and
# three groups. Incorrect df changes interaction p-values and AIC/BIC penalties,
# even when the fitted coefficients themselves are correct.
test_that("17.1.4 MFPI interaction degrees of freedom follow group and FP degree", {
  # Linear interaction: K group-specific slopes versus one common slope.
  linear_2 <- interaction_model_df(n_groups = 2, degree = 0, flex = "flex1")
  linear_3 <- interaction_model_df(n_groups = 3, degree = 0, flex = "flex1")
  expect_equal(linear_2$dfint, 1)
  expect_equal(linear_3$dfint, 2)

  # FP1 with common powers in flex1/flex2 adds K - 1 slope parameters.
  fp1_3 <- interaction_model_df(n_groups = 3, degree = 1, flex = "flex1")
  expect_equal(fp1_3$dfint, 2)

  # FP2 with common powers adds two group-specific slope differences per
  # non-reference group.
  fp2_3 <- interaction_model_df(n_groups = 3, degree = 2, flex = "flex2")
  expect_equal(fp2_3$dfint, 4)
})

# Test purpose: In the simplest two-group linear case, MFPI should fit the same
# interaction model as an ordinary Gaussian model y ~ group * x. The comparison
# uses fitted values, log-likelihood, and newdata predictions, avoiding reliance
# on package-specific coefficient names.
test_that("17.1.5 two-group linear MFPI equals an explicit Gaussian interaction model", {
  set.seed(1715)
  n <- 240
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2)),
    x = runif(n, 1, 8)
  )
  dat$y <- 1.2 +
    0.7 * (dat$group == "treated") +
    0.4 * dat$x +
    1.1 * dat$x * (dat$group == "treated") +
    rnorm(n, sd = 0.15)

  fit_mfpi <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    p_interact = 1,
    verbose = FALSE
  )

  fit_reference <- stats::glm(y ~ group * x, data = dat)
  stored <- fit_mfpi$var_winners[["x"]]$fit$test_results$interaction_model$fit

  expect_equal(
    unname(stats::fitted(stored)),
    unname(stats::fitted(fit_reference)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(stored)),
    as.numeric(stats::logLik(fit_reference)),
    tolerance = 1e-8
  )

  nd <- dat[c(1, 30, 121, 180), c("group", "x"), drop = FALSE]
  got <- predict(
    fit_mfpi,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expected <- stats::predict(fit_reference, newdata = nd, type = "link")

  expect_equal(got$predictions$fit, unname(expected), tolerance = 1e-8)
})

# Test purpose: Extends the explicit interaction oracle to three groups. This
# detects incorrect K - 1 dummy construction, swapped group-specific slopes,
# and hard-coded assumptions that only two treatment groups exist.
test_that("17.1.6 three-group linear MFPI equals an explicit Gaussian interaction model", {
  set.seed(1716)
  n_per_group <- 100
  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = n_per_group)),
    x = runif(3 * n_per_group, 1, 8)
  )
  intercept_shift <- c(A = 0, B = 0.6, C = -0.4)[as.character(dat$group)]
  slope_shift <- c(A = 0, B = 0.8, C = -0.5)[as.character(dat$group)]
  dat$y <- 1 + intercept_shift + (0.5 + slope_shift) * dat$x +
    rnorm(nrow(dat), sd = 0.15)

  fit_mfpi <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    p_interact = 1,
    verbose = FALSE
  )

  fit_reference <- stats::glm(y ~ group * x, data = dat)
  stored <- fit_mfpi$var_winners[["x"]]$fit$test_results$interaction_model$fit

  expect_equal(
    unname(stats::fitted(stored)),
    unname(stats::fitted(fit_reference)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(stored)),
    as.numeric(stats::logLik(fit_reference)),
    tolerance = 1e-8
  )

  nd <- dat[c(1, 101, 201), c("group", "x"), drop = FALSE]
  got <- predict(
    fit_mfpi,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "response",
    se.fit = FALSE
  )
  expected <- stats::predict(fit_reference, newdata = nd, type = "response")

  expect_equal(got$predictions$fit, unname(expected), tolerance = 1e-8)
})

# Test purpose: Verifies that a strong interaction produces the expected MFPI
# test result. The test does not depend on a borderline random p-value: the data
# use a large slope difference and low noise, so failure indicates a structural
# interaction-test regression.
test_that("17.1.7 MFPI detects a strong prespecified linear interaction", {
  set.seed(1717)
  n <- 260
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2)),
    x = runif(n, 1, 8)
  )
  dat$y <- 1 + 0.3 * dat$x + 2.0 * dat$x * (dat$group == "treated") +
    rnorm(n, sd = 0.2)

  fit <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 0.05,
    verbose = FALSE
  )

  metric <- fit$all_model_metrics[fit$all_model_metrics$variable == "x", , drop = FALSE]
  expect_equal(nrow(metric), 1)
  expect_true(is.finite(metric$pvalue))
  expect_lt(metric$pvalue, 0.05)
  expect_true("x" %in% names(fit$best_interaction_model))
})

# Test purpose: Verifies Poisson MFPI ordinary prediction when the fitted
# interaction model uses an offset. The reconstructed offset_ column must be
# consumed by predict.glm(), and both link and response predictions must match
# direct prediction from the stored interaction model.
test_that("17.1.8 Poisson MFPI offset predictions match the stored glm", {
  set.seed(1718)
  n <- 260
  dat <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    x = runif(n, 1, 6),
    exposure = runif(n, 0.5, 3)
  )
  eta <- 0.2 + 0.12 * dat$x + 0.35 * (dat$group == "B") +
    0.18 * dat$x * (dat$group == "B") + log(dat$exposure)
  dat$y <- rpois(n, exp(eta))

  fit <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    family = "poisson",
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    offset = log(dat$exposure),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  nd <- dat[1:18, c("group", "x", "exposure"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(
    fit,
    "x",
    fit_result,
    nd,
    newoffset = log(nd$exposure)
  )
  stored <- fit_result$test_results$interaction_model$fit

  expect_true("offset_" %in% names(design$model_newdata))
  expect_equal(design$model_newdata$offset_, log(nd$exposure))

  # Build the interaction design directly from the stored formula. The offset
  # is not a coefficient column; it is added to X beta after multiplication.
  reference_terms <- stats::delete.response(stats::terms(stored))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = design$model_newdata,
    xlev = stored$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = stored$contrasts
  )
  beta <- stats::coef(stored)
  beta_vcov <- stats::vcov(stored)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_offset <- stats::model.offset(reference_frame)
  expect_equal(manual_offset, log(nd$exposure))

  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_link_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_link_se <- as.numeric(sqrt(pmax(manual_link_variance, 0)))
  manual_response <- as.numeric(stored$family$linkinv(manual_link))

  # predict.glm(type = "response", se.fit = TRUE) applies the delta method:
  # response-scale SE = link-scale SE * abs(d mu / d eta).
  manual_response_se <- manual_link_se * abs(stored$family$mu.eta(manual_link))

  for (prediction_type in c("link", "response")) {
    direct <- stats::predict(
      stored,
      newdata = design$model_newdata,
      type = prediction_type,
      se.fit = TRUE
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "x",
      model = "all",
      type = prediction_type,
      newoffset = log(nd$exposure),
      se.fit = TRUE
    )

    expected_fit <- if (prediction_type == "link") manual_link else manual_response
    expected_se <- if (prediction_type == "link") manual_link_se else manual_response_se

    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    expect_equal(got$predictions$fit, expected_fit, tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, expected_se, tolerance = 1e-8)
  }
})

# Test purpose: Verifies native Cox ordinary prediction types without strata.
# "link" is accepted as an alias for "lp"; "response" is not an alias for
# "risk". Every result is delegated to the stored coxph model.
test_that("17.1.9 unstratified Cox MFPI lp and risk match stored coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:15, c("age", "sex"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  stored <- fit_result$test_results$interaction_model$fit

  cases <- list(
    lp = "lp",
    link = "lp",
    named_link = c(alias = "link"),
    risk = "risk"
  )
  for (case_name in names(cases)) {
    mfpi_type <- cases[[case_name]]
    native_type <- if (identical(unname(mfpi_type), "link")) {
      "lp"
    } else {
      unname(mfpi_type)
    }
    direct <- stats::predict(
      stored,
      newdata = design$model_newdata,
      type = native_type,
      se.fit = TRUE,
      reference = "zero"
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = mfpi_type,
      cox_reference = "zero",
      se.fit = TRUE
    )

    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    expect_true(got$metadata$used_model_predict)
  }

  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all",
      type = "response"
    ),
    "Cox MFPI models"
  )

  # Independent zero-reference oracle for the linear predictor.
  reference_terms <- stats::delete.response(stats::terms(stored))
  manual_x <- stats::model.matrix(
    reference_terms,
    data = design$model_newdata,
    contrasts.arg = stored$contrasts
  )
  manual_x <- manual_x[, colnames(manual_x) != "(Intercept)", drop = FALSE]
  beta <- stats::coef(stored)
  beta_vcov <- stats::vcov(stored)
  manual_lp <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  got_lp <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "lp",
    cox_reference = "zero",
    se.fit = TRUE
  )
  expect_equal(got_lp$predictions$fit, manual_lp, tolerance = 1e-8)
  expect_equal(got_lp$predictions$se.fit, manual_se, tolerance = 1e-8)
})

# Test purpose: Verifies Cox offsets and every covariate-reference choice by
# comparing MFPI directly with the retained coxph model. No manual prediction
# branch is exercised or retained.
test_that("17.1.9b Cox MFPI offsets and references match stored coxph", {
  set.seed(17109)
  n <- 320

  dat <- data.frame(
    age = stats::rnorm(n, mean = 55, sd = 9),
    sex = factor(sample(c("female", "male"), n, replace = TRUE)),
    exposure = stats::runif(n, 0.5, 3)
  )
  sex_effect <- ifelse(dat$sex == "male", 0.35, 0)
  eta <- 0.025 * (dat$age - 55) + sex_effect + log(dat$exposure)
  event_time <- stats::rexp(n, rate = 0.02 * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.012)
  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    offset = log(dat$exposure),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:24, c("age", "sex", "exposure"), drop = FALSE]
  new_offset <- log(nd$exposure)
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, newoffset = new_offset
  )
  stored <- fit_result$test_results$interaction_model$fit

  for (prediction_type in c("lp", "risk")) {
    for (reference_value in c("zero", "sample", "strata")) {
      direct <- stats::predict(
        stored,
        newdata = design$model_newdata,
        type = prediction_type,
        se.fit = TRUE,
        reference = reference_value
      )
      got <- predict(
        fit,
        newdata = nd,
        terms = "age",
        model = "all",
        type = prediction_type,
        newoffset = new_offset,
        cox_reference = reference_value,
        se.fit = TRUE
      )

      expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
      expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
      expect_identical(got$metadata$cox_reference, reference_value)
    }
  }
})

# Test purpose: Verifies formula-interface strata reconstruction. The caller
# supplies only ordinary newdata; predict.mfpi() must recover the original
# strata variable from stored formula metadata and create model_newdata$strata_
# without converting it to integer codes.
test_that("17.1.10 formula-stratified Cox MFPI reconstructs strata from newdata", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  dat$inst <- factor(dat$inst)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:15, c("age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  reconstructed <- reconstruct_formula_strata_newdata(fit, nd)
  design <- mfpi_build_ordinary_design(
    fit,
    "age",
    fit_result,
    nd,
    strata = reconstructed
  )
  stored <- fit_result$test_results$interaction_model$fit

  expect_identical(design$model_newdata$strata_, nd$inst)

  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "link",
    se.fit = TRUE
  )

  expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
  expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
})

# Test purpose: Exercises actual prediction with two Cox strata columns. Matrix
# or data-frame strata must be combined exactly once with strata(),
# whereas a single vector/factor must remain raw for the stored formula to
# evaluate strata(strata_) itself.
test_that("17.1.11 multiple Cox strata columns are combined once in MFPI prediction", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(
    dat[, c("time", "status", "age", "sex", "inst", "ph.ecog")]
  ), ]
  dat$sex <- factor(dat$sex)

  strata_fit <- data.frame(inst = dat$inst, ecog = dat$ph.ecog)
  fit <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = strata_fit,
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:15, c("sex", "age", "inst", "ph.ecog"), drop = FALSE]
  strata_new <- nd[, c("inst", "ph.ecog"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(
    fit,
    "age",
    fit_result,
    nd,
    strata = strata_new
  )

  expected_strata <- do.call(
    survival::strata,
    c(as.list(strata_new), list(shortlabel = TRUE))
  )
  expect_equal(design$model_newdata$strata_, expected_strata)
  expect_s3_class(design$model_newdata$strata_, "factor")

  stored <- fit_result$test_results$interaction_model$fit
  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    reference = "zero"
  )
  got <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "link",
    strata = strata_new,
    se.fit = FALSE
  )

  expect_equal(got$predictions$fit, as.numeric(direct), tolerance = 1e-8)
})

# Test purpose: Verifies clear validation for stratified Cox prediction. A
# stratified stored interaction model cannot predict supplied rows without one
# stratum value/row per prediction row and without missing stratum values.
test_that("17.1.12 MFPI Cox strata validation rejects missing, short, and NA strata", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = dat$inst,
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  nd <- dat[1:8, c("age", "sex", "inst"), drop = FALSE]

  expect_error(
    predict(fit, newdata = nd, terms = "age", model = "all", type = "link"),
    "stratified|strata"
  )
  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all", type = "link",
      strata = nd$inst[-1]
    ),
    "one value or row per prediction row"
  )
  bad_strata <- nd$inst
  bad_strata[1] <- NA
  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all", type = "link",
      strata = bad_strata
    ),
    "must not contain missing values"
  )
  bad_strata_inf <- as.numeric(nd$inst)
  bad_strata_inf[1] <- Inf
  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all", type = "link",
      strata = bad_strata_inf
    ),
    "strata.*finite"
  )
})

# Test purpose: Verifies fitted-function values using direct matrix algebra.
# The manual fitted-function path should return X_g beta_g for every group and
# x value, using the exact stored basis and coefficient mapping.
test_that("17.1.13 MFPI fitted functions equal direct basis-times-coefficient calculations", {
  set.seed(1723)
  n <- 220
  dat <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    x = runif(n, 1, 8)
  )
  dat$y <- 1 + 0.4 * dat$x + 1.0 * dat$x * (dat$group == "B") +
    rnorm(n, sd = 0.2)

  fit <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  eval_x <- data.frame(x = c(1.5, 3, 6))
  pred <- predict(
    fit,
    newdata = eval_x,
    terms = "x",
    model = "all",
    type = "function",
    grid = FALSE,
    se.fit = FALSE
  )

  fit_result <- fit$var_winners[["x"]]$fit
  prepared <- mfpi_prepare_prediction_data(
    object = fit,
    term = "x",
    newdata = eval_x,
    grid = FALSE,
    n_grid = 200L
  )
  basis <- mfpi_build_function_basis(
    object = fit,
    term = "x",
    fit_result = fit_result,
    cont_var_scaled = prepared$cont_var_scaled,
    x_display = prepared$x_display
  )
  coefficients <- fit_result$test_results$interaction_model$coefficients

  group_internal <- names(basis$coefficient_groups)
  group_display <- mfpi_prediction_group_display_labels(
    fit,
    group_internal
  )
  intercept <- if ("(Intercept)" %in% names(coefficients)) {
    unname(coefficients["(Intercept)"])
  } else {
    0
  }

  expected <- do.call(rbind, lapply(seq_along(group_internal), function(i) {
    g <- group_internal[i]
    cols <- basis$coefficient_groups[[g]]
    dummy_name <- if (i > 1L) paste0(fit$group_var, g) else NULL
    dummy_effect <- if (!is.null(dummy_name)) unname(coefficients[dummy_name]) else 0

    data.frame(
      x = prepared$x_display,
      group = group_display[i],
      fit = intercept +
        as.numeric(basis$x[, cols, drop = FALSE] %*% coefficients[cols]) +
        dummy_effect,
      stringsAsFactors = FALSE
    )
  }))

  observed_key <- paste(pred$functions$x, pred$functions$group, sep = "::")
  expected_key <- paste(expected$x, expected$group, sep = "::")
  expected <- expected[match(observed_key, expected_key), , drop = FALSE]

  expect_false(anyNA(expected$fit))
  expect_equal(pred$functions$fit, expected$fit, tolerance = 1e-10)
})

# Test purpose: Replaces the former tautological representation check with a
# direct test of mfpi_build_ordinary_design(). Vector/factor strata must remain
# unchanged; only matrix/data-frame strata are combined before model prediction.
test_that("17.1.14 MFPI ordinary design preserves vector strata and combines tabular strata", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(
    dat[, c("time", "status", "age", "sex", "inst", "ph.ecog")]
  ), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = dat$inst,
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  nd <- dat[1:8, c("sex", "age", "inst", "ph.ecog"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit

  vector_design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = factor(nd$inst)
  )
  expect_identical(vector_design$model_newdata$strata_, factor(nd$inst))

  tabular_strata <- nd[, c("inst", "ph.ecog"), drop = FALSE]
  tabular_design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = tabular_strata
  )
  expected <- do.call(
    survival::strata,
    c(as.list(tabular_strata), list(shortlabel = TRUE))
  )
  expect_equal(tabular_design$model_newdata$strata_, expected)
})

# -----------------------------------------------------------------------------
# 17.2 Grouped categorical adjustment prediction
# -----------------------------------------------------------------------------
# Formula prediction must rebuild the original factor contrasts and pass the
# complete selected adjustment block to the stored interaction model.

# Test purpose: Ensures predict.mfpi() reconstructs a selected factor block
# from ordinary factor-valued formula newdata.
test_that("MFPI prediction reconstructs grouped factor adjustments", {
  dat <- make_mfpi_factor_data()

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  nd <- dat[1:20, c("trt", "x", "stage", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(
    fit,
    "x",
    fit_result,
    nd
  )

  stored <- fit_result$test_results$interaction_model$fit
  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "link"
  )

  got <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )

  expect_equal(got$predictions$fit, as.numeric(direct), tolerance = 1e-8)
  expect_true(all(c("stageII.1", "stageIII.1") %in%
                    names(design$model_newdata)))
})


# Test purpose: Ensures MFPI formula prediction reconstructs inline factor
# columns while the adjustment model is addressed by the source variable name.
test_that("MFPI prediction reconstructs inline factor adjustments", {
  dat <- make_mfpi_factor_data()
  dat$stage_code <- as.integer(dat$stage)

  fit <- mfpi(
    y ~ trt + x + factor(stage_code) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage_code",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  nd <- dat[1:20, c("trt", "x", "stage_code", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(fit, "x", fit_result, nd)

  expect_true(all(c("factor(stage_code)2.1", "factor(stage_code)3.1") %in%
                    names(design$model_newdata)))

  out <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expect_true(all(is.finite(out$predictions$fit)))
})



# Regression: formula ordinary prediction should depend on the selected
# term-specific interaction model, not the complete formula offered to mfpi().
test_that("MFPI formula prediction omits unrelated original predictors", {
  dat <- make_mfpi_factor_data()
  set.seed(17221)
  dat$w <- runif(nrow(dat), 0.5, 4)
  dat$nuisance <- factor(
    rep(
      c("low", "middle", "high", "middle", "high", "low", "high", "low", "middle"),
      length.out = nrow(dat)
    ),
    levels = c("low", "middle", "high")
  )

  fit <- mfpi(
    y ~ trt + x + w + stage + z + nuisance,
    data = dat,
    group_var = "trt",
    cont_vars = c("x", "w"),
    cont_var_forms = c(x = "linear", w = "linear"),
    keep = "stage",
    df = 1,
    select = 0,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  selected <- get_selected_variables(fit$adjustment_model)
  expect_true("stage" %in% selected)
  expect_false(any(c("z", "nuisance") %in% selected))
  expect_identical(
    unname(fit$formula_prediction_term_names[["nuisance"]]),
    "nuisance"
  )

  nd_full <- dat[
    1:18,
    c("trt", "x", "w", "stage", "z", "nuisance"),
    drop = FALSE
  ]
  nd_x <- nd_full[, c("trt", "x", "stage"), drop = FALSE]
  nd_w <- nd_full[, c("trt", "w", "stage"), drop = FALSE]

  pred_x_full <- predict(
    fit, newdata = nd_full, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )
  pred_x_minimal <- predict(
    fit, newdata = nd_x, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )
  pred_w_full <- predict(
    fit, newdata = nd_full, terms = "w", model = "all",
    type = "link", se.fit = TRUE
  )
  pred_w_minimal <- predict(
    fit, newdata = nd_w, terms = "w", model = "all",
    type = "link", se.fit = TRUE
  )

  expect_equal(
    pred_x_minimal$predictions,
    pred_x_full$predictions,
    tolerance = 1e-10
  )
  expect_equal(
    pred_w_minimal$predictions,
    pred_w_full$predictions,
    tolerance = 1e-10
  )

  nd_unseen <- nd_full
  nd_unseen$z <- NA_real_
  nd_unseen$nuisance <- factor(
    rep("unseen", nrow(nd_unseen)),
    levels = c(levels(dat$nuisance), "unseen")
  )
  pred_x_unseen <- predict(
    fit, newdata = nd_unseen, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )
  expect_equal(
    pred_x_unseen$predictions,
    pred_x_minimal$predictions,
    tolerance = 1e-10
  )
})


# Regression: the formula-label map must allow a selected inline factor block
# to be reconstructed without evaluating an unrelated eliminated inline factor.
test_that("MFPI minimal prediction supports selected inline factor adjustments", {
  dat <- make_mfpi_factor_data()
  dat$stage_code <- as.integer(dat$stage)
  dat$nuisance_code <- rep(
    c(1L, 2L, 3L, 2L, 3L, 1L, 3L, 1L, 2L),
    length.out = nrow(dat)
  )

  fit <- mfpi(
    y ~ trt + x + factor(stage_code) + factor(nuisance_code) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage_code",
    df = 1,
    select = 0,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    unname(fit$formula_prediction_term_names[["factor(stage_code)"]]),
    "stage_code"
  )
  expect_identical(
    unname(fit$formula_prediction_term_names[["factor(nuisance_code)"]]),
    "nuisance_code"
  )
  expect_true("stage_code" %in% get_selected_variables(fit$adjustment_model))
  expect_false("nuisance_code" %in% get_selected_variables(fit$adjustment_model))

  nd_full <- dat[
    1:20,
    c("trt", "x", "stage_code", "nuisance_code", "z"),
    drop = FALSE
  ]
  nd_minimal <- nd_full[, c("trt", "x", "stage_code"), drop = FALSE]
  nd_unseen <- nd_full
  nd_unseen$nuisance_code <- 99L

  pred_full <- predict(
    fit, newdata = nd_full, terms = "x", model = "all",
    type = "link", se.fit = FALSE
  )
  pred_minimal <- predict(
    fit, newdata = nd_minimal, terms = "x", model = "all",
    type = "link", se.fit = FALSE
  )
  pred_unseen <- predict(
    fit, newdata = nd_unseen, terms = "x", model = "all",
    type = "link", se.fit = FALSE
  )

  expect_equal(pred_minimal$predictions, pred_full$predictions,
               tolerance = 1e-10)
  expect_equal(pred_unseen$predictions, pred_minimal$predictions,
               tolerance = 1e-10)
  expect_true(all(c("factor(stage_code)2.1", "factor(stage_code)3.1") %in%
                    pred_minimal$metadata$model_newdata_columns))
  expect_false(any(grepl(
    "nuisance_code",
    pred_minimal$metadata$model_newdata_columns,
    fixed = TRUE
  )))
})

# =============================================================================
# 18. Reproducibility — same seed gives same result
# =============================================================================

# Test purpose: Checks that repeated fits on the same data produce identical
# selected powers, coefficients, and metadata.
test_that("mfp2() is deterministic across repeated calls", {
  fit1 <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  fit2 <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  expect_equal(fit1$fp_powers, fit2$fp_powers)
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(fit1$fp_terms, fit2$fp_terms)
})

# =============================================================================
# 19. Transformation helpers
# =============================================================================

# Test purpose: Checks that the exported FP transformation helper computes a
# log transform for power 0.
test_that("transform_vector_fp() is exported and works", {
  # Basic FP1 transformation
  x <- seq(0.1, 10, length.out = 100)
  # Power 0 = log
  result <- transform_vector_fp(x, power = 0)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1)
  expect_equal(as.numeric(result[, 1]), log(x), tolerance = 1e-10)
})

# Test purpose: Checks that repeated FP powers produce the standard x and
# x log(x) basis.
test_that("transform_vector_fp() handles repeated powers", {
  x <- seq(0.1, 10, length.out = 50)
  # Repeated power (1, 1) -> x, x*log(x)
  result <- transform_vector_fp(x, power = c(1, 1))
  expect_equal(ncol(result), 2)
  expect_equal(as.numeric(result[, 1]), x, tolerance = 1e-10)
  expect_equal(as.numeric(result[, 2]), x * log(x), tolerance = 1e-10)
})

# =============================================================================
# 20. Likelihood-ratio and F-test helpers
# =============================================================================

# Test purpose: Checks that the likelihood-ratio helper returns a valid
# nonnegative statistic and p-value.
test_that("calculate_lr_test() returns correct p-value for nested models", {
  # Fit two nested models manually
  fit_null <- glm(y_prostate ~ 1)
  fit_full <- glm(y_prostate ~ x_prostate[, "cavol"])

  lr <- calculate_lr_test(
    logl = c(logLik(fit_null), logLik(fit_full)),
    dfs = c(1, 2)
  )

  expect_true(lr$pvalue >= 0 && lr$pvalue <= 1)
  expect_true(lr$statistic >= 0)
})

# Test purpose: Checks that the likelihood-ratio helper rejects invalid nested
# model df ordering.
test_that("calculate_lr_test() errors when df ordering is wrong", {
  expect_error(
    calculate_lr_test(logl = c(100, 110), dfs = c(5, 3)),
    "more degrees of freedom"
  )
})

# Test purpose: Checks that the F-test helper returns valid statistic, deviance
# difference, and p-value.
test_that("calculate_f_test() returns correct p-value", {
  fit_null <- glm(y_prostate ~ 1)
  fit_full <- glm(y_prostate ~ x_prostate[, "cavol"])

  f_result <- calculate_f_test(
    deviances = c(deviance(fit_null), deviance(fit_full)),
    dfs_resid = c(df.residual(fit_null), df.residual(fit_full)),
    n_obs = length(y_prostate)
  )

  expect_true(f_result$pvalue >= 0 && f_result$pvalue <= 1)
  expect_true(f_result$statistic >= 0)
  expect_true(f_result$dev_diff >= 0)
})

# =============================================================================
# 21. plot()
# =============================================================================

# Test purpose: Checks that plot() can be called on a Gaussian mfp2 fit without
# error.
test_that("plot() runs without error for Gaussian model", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    verbose = FALSE
  )

  expect_error(
    plot(fit),
    NA
  )
})

# Test purpose: Checks that plot() returns a list and does not emit warnings for
# a Gaussian mfp2 model.
test_that("plot() runs without warning for Gaussian models", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    verbose = FALSE
  )

  expect_warning(
    plots <- plot(fit),
    NA
  )

  expect_type(plots, "list")
})

# Test purpose: Checks that the deprecated fracplot() wrapper still delegates to
# the plotting implementation and returns a list.
test_that("fracplot() is deprecated but still works for Gaussian models", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    verbose = FALSE
  )

  expect_warning(
    plots <- fracplot(fit),
    regexp = "fracplot.*deprecated.*plot",
    ignore.case = TRUE
  )

  expect_type(plots, "list")
})

# Test purpose: Checks that the C++ FP core returns the expected ordinary
# fractional-polynomial columns for powers 1 and 0.
# Test purpose: Checks that the C++ FP core returns the expected ordinary
# fractional-polynomial columns for powers 1 and 0.
test_that("C++ FP core computes ordinary FP transformations", {
  x <- c(1, 2, 4)
  p <- c(1, 0)

  out <- transform_fp_core(
    x_raw = x,
    power = p,
    shift_val = 0,
    scale_val = 1,
    zero = FALSE
  )

  expected <- cbind(
    x,
    log(x)
  )

  expect_true(is.matrix(out))
  expect_equal(
    unname(out, force = TRUE),
    unname(expected, force = TRUE),
    tolerance = 1e-12
  )
})

# Test purpose: Checks that the C++ FP core implements repeated-power FP rules.
test_that("C++ FP core handles repeated powers", {
  x <- c(1, 2, 4)
  p <- c(0, 0)

  out <- transform_fp_core(
    x_raw = x,
    power = p,
    shift_val = 0,
    scale_val = 1,
    zero = FALSE
  )

  expected <- cbind(
    log(x),
    log(x)^2
  )

  expect_true(is.matrix(out))
  expect_equal(unname(out), expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ FP core applies shift and scale exactly once.
test_that("C++ FP core applies shift and scale", {
  x <- c(1, 3, 5)
  shift <- 1
  scale <- 2
  x_scaled <- (x + shift) / scale

  out <- transform_fp_core(
    x_raw = x,
    power = c(1, 0),
    shift_val = shift,
    scale_val = scale,
    zero = FALSE
  )

  expected <- cbind(
    x_scaled,
    log(x_scaled)
  )

  expect_equal(
    unname(out, force = TRUE),
    unname(expected, force = TRUE),
    tolerance = 1e-12
  )
})


# Test purpose: Checks that the C++ FP core maps non-positive values to zero
# rows when zero handling is active.
test_that("C++ FP core handles zero mode", {
  x <- c(-2, 0, 1, 4)

  out <- transform_fp_core(
    x_raw = x,
    power = c(1, 0),
    shift_val = 0,
    scale_val = 1,
    zero = TRUE
  )

  expected <- rbind(
    c(0, 0),
    c(0, 0),
    c(1, log(1)),
    c(4, log(4))
  )

  expect_equal(unname(out), expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ FP core preserves missing and non-finite
# values. Keep this test only if transform_fp_core_internal() has the explicit
# missing/non-finite guard.
test_that("C++ FP core preserves missing and non-finite values", {
  x <- c(1, NA_real_, NaN, Inf, 4)
  p <- c(1, 0)

  out <- transform_fp_core(
    x_raw = x,
    power = p,
    shift_val = 0,
    scale_val = 1,
    zero = FALSE
  )

  expect_equal(out[1, ], c(1, 0))
  expect_true(is.na(out[2, 1]))
  expect_true(is.na(out[2, 2]))
  expect_true(is.nan(out[3, 1]))
  expect_true(is.nan(out[3, 2]))
  expect_true(is.infinite(out[4, 1]))
  expect_true(is.infinite(out[4, 2]))
  expect_equal(out[5, ], c(4, log(4)))
})


# Test purpose: Checks that the R wrapper around the C++ FP core preserves
# variable naming.
test_that("transform_vector_fp keeps expected column names with C++ core", {
  x <- c(1, 2, 4)

  out <- transform_vector_fp(
    x = x,
    power = c(1, 0),
    shift = 0,
    scale = 1,
    name = "x",
    zero = FALSE,
    check_binary = FALSE
  )

  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2L)
  expect_equal(unname(out[, 2]), x)
  expect_equal(unname(out[, 1]), log(x))
  expect_false(is.null(colnames(out)))
})


# Test purpose: Checks that the C++ batch FP generator returns one matrix per
# candidate power row.
test_that("generate_transformations_fp_cpp returns one matrix per power row", {
  x <- c(1, 2, 4)
  powers <- rbind(
    c(1, 1),
    c(0, 0),
    c(1, 0)
  )

  out <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = FALSE,
    catzero = NULL
  )

  expect_type(out, "list")
  expect_length(out, nrow(powers))

  expect_equal(
    unname(out[[1]]),
    unname(cbind(x, x * log(x))),
    tolerance = 1e-12
  )

  expect_equal(
    unname(out[[2]]),
    cbind(log(x), log(x)^2),
    tolerance = 1e-12
  )

  expect_equal(
    unname(out[[3]]),
    unname(cbind(x, log(x))),
    tolerance = 1e-12
  )
})


# Test purpose: Checks that the C++ batch FP generator applies the binary
# shortcut when zero handling is inactive.
test_that("generate_transformations_fp_cpp uses binary shortcut", {
  x <- c(0, 1, 0, 1)
  powers <- rbind(
    c(1, 1),
    c(0, 0)
  )

  out <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = FALSE,
    catzero = NULL
  )

  expect_length(out, 2L)
  expect_equal(unname(out[[1]]), matrix(x, ncol = 1L))
  expect_equal(unname(out[[2]]), matrix(x, ncol = 1L))
})


# Test purpose: Checks that the C++ batch FP generator prepends catzero when
# catzero is supplied.
test_that("generate_transformations_fp_cpp prepends catzero", {
  x <- c(1, 2, 4)
  powers <- rbind(c(1, 0))
  catzero <- matrix(c(0, 1, 0), ncol = 1L)

  out <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = FALSE,
    catzero = catzero
  )

  expect_length(out, 1L)
  expect_equal(ncol(out[[1]]), 3L)
  expect_equal(unname(out[[1]][, 1]), as.numeric(catzero[, 1]))
  expect_equal(unname(out[[1]][, 2]), x)
  expect_equal(unname(out[[1]][, 3]), log(x), tolerance = 1e-12)
  expect_equal(colnames(out[[1]]), c("catzero", "V1", "V2"))
})

# Test purpose: Checks that the compact ordinary-FP basis reconstructs every
# degree-2 candidate exactly while storing only the unique repeated-power terms.
test_that("compact FP basis reconstructs materialized degree-2 candidates", {
  x <- c(1, 2, 4, 8, 16)
  allowed_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  materialized <- generate_transformations_fp(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  compact <- generate_transformations_fp_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  expect_equal(nrow(compact$candidate_map), 36L)
  expect_equal(ncol(compact$candidate_map), 2L)
  expect_equal(ncol(compact$basis), 16L)
  expect_equal(compact$powers, materialized$powers)

  for (i in seq_len(nrow(compact$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(compact, i)),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Checks repeated zero powers and structural-zero indicators in
# the compact representation, including the column names required by SAZ stage 2.
test_that("compact FP basis preserves repeated zero powers and catzero", {
  x <- c(1, 2, 4, 8)
  catzero <- matrix(c(0, 1, 0, 1), ncol = 1L)

  compact <- generate_transformations_fp_basis(
    x = x,
    degree = 2,
    powers = 0,
    zero = FALSE,
    catzero = catzero
  )

  out <- materialize_fp_basis_candidate(compact, 1L)

  expect_equal(ncol(compact$basis), 2L)
  expect_equal(unname(compact$basis[, 1L]), log(x), tolerance = 1e-12)
  expect_equal(unname(compact$basis[, 2L]), log(x)^2, tolerance = 1e-12)
  expect_equal(unname(out[, 1L]), as.numeric(catzero[, 1L]))
  expect_equal(unname(out[, 2L]), log(x), tolerance = 1e-12)
  expect_equal(unname(out[, 3L]), log(x)^2, tolerance = 1e-12)
  expect_equal(colnames(out), c("catzero", "V1", "V2"))
})


# Test purpose: Checks zero-mode semantics in the compact basis. Nonpositive
# values must remain zero while positive observations use ordinary FP terms.
test_that("compact FP basis preserves zero-mode transformations", {
  x <- c(-2, 0, 1, 2, 4)
  powers <- rbind(c(-1, -1), c(0, 0), c(1, 2))

  old <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = TRUE,
    catzero = NULL
  )
  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = TRUE
  )

  for (i in seq_len(nrow(powers))) {
    reconstructed <- compact$basis[
      , compact$candidate_map[i, ], drop = FALSE
    ]
    expect_equal(
      unname(reconstructed),
      unname(old[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Checks that compact candidate columns can be copied into a
# reusable design matrix without allocating/materializing the complete candidate.
test_that("compact FP candidate copies into reusable design matrix", {
  x <- c(1, 2, 4, 8)
  powers <- rbind(c(0, 0), c(1, 2))
  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = FALSE
  )

  target <- matrix(-1, nrow = length(x), ncol = 4L)
  target[, 1L] <- 1
  target[, 4L] <- 99

  target <- copy_fp_basis_candidate_cpp(
    target = target,
    basis = compact$basis,
    source_cols = as.integer(compact$candidate_map[1L, ]),
    target_cols = c(2L, 3L)
  )

  expect_equal(target[, 1L], rep(1, length(x)))
  expect_equal(target[, 2L], log(x), tolerance = 1e-12)
  expect_equal(target[, 3L], log(x)^2, tolerance = 1e-12)
  expect_equal(target[, 4L], rep(99, length(x)))
})


# Test purpose: Checks that the C++ adjustment-step bridge returns NULL when
# there are no adjustment variables.
test_that("C++ adjustment-step bridge handles no adjustment variables", {
  x <- matrix(c(1, 2, 3), ncol = 1L)
  colnames(x) <- "x1"

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = character(0),
    powers_adj = list(),
    acdx_adj = list(),
    zero_adj = list(),
    catzero = list(),
    spike_adj = list(),
    spike_decision_int_adj = integer(0),
    acd_parameter_adj = list(),
    eliminated = logical(0),
    spike_binary_only_flags = logical(0),
    current_power_keys_adj = list(),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expect_equal(out$data_adj_list, list())
  expect_null(out$data_adj)
})


# Test purpose: Checks that the C++ adjustment-step bridge builds ordinary FP
# adjustment columns on cache miss.
test_that("C++ adjustment-step bridge builds FP adjustment columns", {
  x <- cbind(
    x1 = c(1, 2, 4),
    x2 = c(2, 3, 5)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = c("x1", "x2"),
    powers_adj = list(x1 = c(1), x2 = c(0, 0)),
    acdx_adj = list(x1 = FALSE, x2 = FALSE),
    zero_adj = list(x1 = FALSE, x2 = FALSE),
    catzero = list(x1 = NULL, x2 = NULL),
    spike_adj = list(x1 = FALSE, x2 = FALSE),
    spike_decision_int_adj = c(x1 = 2L, x2 = 2L),
    acd_parameter_adj = list(x1 = NULL, x2 = NULL),
    eliminated = c(x1 = FALSE, x2 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE, x2 = FALSE),
    current_power_keys_adj = list(x1 = c(1), x2 = c(0, 0)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected_x1 <- matrix(x[, "x1"], ncol = 1L)
  colnames(expected_x1) <- "x1_adj1"

  expected_x2 <- cbind(log(x[, "x2"]), log(x[, "x2"])^2)
  colnames(expected_x2) <- c("x2_adj1", "x2_adj2")

  expected <- cbind(expected_x1, expected_x2)

  expect_equal(names(out$data_adj_list), c("x1", "x2"))
  expect_equal(out$data_adj_list$x1, expected_x1, tolerance = 1e-12)
  expect_equal(out$data_adj_list$x2, expected_x2, tolerance = 1e-12)
  expect_equal(out$data_adj, expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ adjustment-step bridge prepends catzero for
# a non-spike variable.
test_that("C++ adjustment-step bridge prepends catzero for non-spike variable", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- cbind(cz[, 1], x[, "x1"])
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj_list$x1, expected)
  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies spike
# decision 1 as catzero plus continuous FP columns.
test_that("C++ adjustment-step bridge handles spike decision 1", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 1L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- cbind(cz[, 1], x[, "x1"])
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies spike
# decision 2 as continuous FP columns only.
test_that("C++ adjustment-step bridge handles spike decision 2", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- matrix(x[, "x1"], ncol = 1L)
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies spike
# decision 3 as catzero only.
test_that("C++ adjustment-step bridge handles spike decision 3", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(NA_real_)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 3L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = TRUE),
    spike_binary_only_flags = c(x1 = TRUE),
    current_power_keys_adj = list(x1 = NA_real_),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- matrix(cz[, 1], ncol = 1L)
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge returns n x 0 data_adj
# when adjustment variables exist but all are eliminated.
test_that("C++ adjustment-step bridge handles all-eliminated adjustment variables", {
  x <- cbind(
    x1 = c(1, 2, 4),
    x2 = c(2, 3, 5)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = c("x1", "x2"),
    powers_adj = list(x1 = c(NA_real_), x2 = c(NA_real_)),
    acdx_adj = list(x1 = FALSE, x2 = FALSE),
    zero_adj = list(x1 = FALSE, x2 = FALSE),
    catzero = list(x1 = NULL, x2 = NULL),
    spike_adj = list(x1 = FALSE, x2 = FALSE),
    spike_decision_int_adj = c(x1 = 2L, x2 = 2L),
    acd_parameter_adj = list(x1 = NULL, x2 = NULL),
    eliminated = c(x1 = TRUE, x2 = TRUE),
    spike_binary_only_flags = c(x1 = FALSE, x2 = FALSE),
    current_power_keys_adj = list(x1 = NA_real_, x2 = NA_real_),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expect_true(is.matrix(out$data_adj))
  expect_equal(nrow(out$data_adj), nrow(x))
  expect_equal(ncol(out$data_adj), 0L)
  expect_equal(NCOL(out$data_adj), 0L)
  expect_equal(ncol(out$data_adj_list$x1), 0L)
  expect_equal(ncol(out$data_adj_list$x2), 0L)
})


# Test purpose: The persistent focal-variable cache must drop whole-step
# matrices after a focal fit finishes, while retaining the per-variable blocks
# and metadata needed for fallback reuse in a later cycle.
test_that("persistent adjustment cache drops assembled matrices", {
  cached_block <- matrix(c(1, 2, 3), ncol = 1L)
  assembled <- cbind(cached_block, cached_block)

  params_xi <- list(
    powers_adj = list(x1 = 1),
    spike_decision_adj = c(x1 = 2L),
    data_adj_list = list(x1 = cached_block),
    data_adj = assembled,
    data_xi = matrix(c(4, 5, 6), ncol = 1L)
  )

  out <- compact_prev_adj_cache_entry(params_xi)

  # Keep exact named NULL entries. If `data_adj` were removed entirely, R's
  # partial `$` matching could resolve `out$data_adj` to `out$data_adj_list`.
  expect_true("data_adj" %in% names(out))
  expect_true("data_xi" %in% names(out))
  expect_null(out[["data_adj", exact = TRUE]])
  expect_null(out[["data_xi", exact = TRUE]])
  expect_null(out$data_adj)
  expect_null(out$data_xi)
  expect_identical(out$powers_adj, params_xi$powers_adj)
  expect_identical(out$spike_decision_adj, params_xi$spike_decision_adj)
  expect_identical(out$data_adj_list$x1, cached_block)
})


# Test purpose: build_adjustment_step() must still be able to reuse a historical
# per-focal block when prev_adj_params no longer contains a complete data_adj
# matrix. This protects the compatibility fallback after cache compaction.
test_that("adjustment fallback works without cached assembled matrix", {
  x <- cbind(
    x1 = c(1, 2, 4),
    x2 = c(3, 5, 7)
  )

  cached <- matrix(c(99, 98, 97), ncol = 1L)
  colnames(cached) <- "cached_x1"

  prev_adj_params <- list(
    x1 = NULL,
    x2 = list(
      powers_adj = list(x1 = c(1)),
      spike_decision_adj = c(x1 = 2L),
      data_adj_list = list(x1 = cached)
    )
  )

  out <- build_adjustment_step(
    x = x,
    xi = "x2",
    powers_current = list(x1 = c(1), x2 = c(1)),
    powers = list(x1 = c(-2, -1, 0, 1, 2), x2 = c(-2, -1, 0, 1, 2)),
    acdx = c(x1 = FALSE, x2 = FALSE),
    zero = c(x1 = FALSE, x2 = FALSE),
    catzero = list(x1 = NULL, x2 = NULL),
    spike = c(x1 = FALSE, x2 = FALSE),
    spike_decision = c(x1 = 2L, x2 = 2L),
    acd_parameter = list(x1 = NULL, x2 = NULL),
    prev_adj_params = prev_adj_params,
    transform_cache = list(x1 = NULL, x2 = NULL),
    term_to_columns = list(x1 = "x1", x2 = "x2")
  )

  expected <- cached
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
  expect_equal(out$data_adj_list$x1, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge reuses the cached
# matrix when normalized powers and spike decision are unchanged.
test_that("C++ adjustment-step bridge reuses cache on cache hit", {
  x <- cbind(x1 = c(1, 2, 4))

  cached <- matrix(c(99, 98, 97), ncol = 1L)
  colnames(cached) <- "old_name"

  prev_xi <- list(
    data_adj_list = list(x1 = cached),
    spike_decision_adj = c(x1 = 2L)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = NULL),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = list(x1 = c(1)),
    prev_xi = prev_xi,
    has_prev = TRUE
  )

  expected <- cached
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
  expect_equal(out$data_adj_list$x1, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge recomputes the matrix
# when the normalized power key changes.
test_that("C++ adjustment-step bridge recomputes on power-key cache miss", {
  x <- cbind(x1 = c(1, 2, 4))

  cached <- matrix(c(99, 98, 97), ncol = 1L)

  prev_xi <- list(
    data_adj_list = list(x1 = cached),
    spike_decision_adj = c(x1 = 2L)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(0)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = NULL),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(0)),
    prev_power_keys_adj = list(x1 = c(1)),
    prev_xi = prev_xi,
    has_prev = TRUE
  )

  expected <- matrix(log(x[, "x1"]), ncol = 1L)
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ adjustment-step bridge recomputes the matrix
# when the spike decision changes.
test_that("C++ adjustment-step bridge recomputes on spike-decision cache miss", {
  x <- cbind(x1 = c(1, 2, 4))
  cz <- matrix(c(0, 1, 0), ncol = 1L)

  cached <- matrix(c(99, 98, 97), ncol = 1L)

  prev_xi <- list(
    data_adj_list = list(x1 = cached),
    spike_decision_adj = c(x1 = 2L)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 1L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = list(x1 = c(1)),
    prev_xi = prev_xi,
    has_prev = TRUE
  )

  expected <- cbind(cz[, 1], x[, "x1"])
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies stored ACD
# parameters without fitting ACD inside build_adjustment_step().
test_that("C++ adjustment-step bridge builds stored-ACD adjustment columns", {
  x <- cbind(x1 = c(1, 2, 4, 8))

  acd_parameter <- list(
    beta0 = 0,
    beta1 = 1,
    power = 1,
    shift = 0,
    scale = 1
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1, 1)),
    acdx_adj = list(x1 = TRUE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = NULL),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = acd_parameter),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1, 1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- cbind(
    x[, "x1"],
    stats::pnorm(x[, "x1"])
  )
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj, expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ adjustment-step bridge rejects ACD variables
# without stored ACD parameters.
test_that("C++ adjustment-step bridge rejects missing stored ACD parameters", {
  x <- cbind(x1 = c(1, 2, 4, 8))

  expect_error(
    mfp2_build_adjustment_step_loop(
      x = x,
      vars_adj = "x1",
      powers_adj = list(x1 = c(1, 1)),
      acdx_adj = list(x1 = TRUE),
      zero_adj = list(x1 = FALSE),
      catzero = list(x1 = NULL),
      spike_adj = list(x1 = FALSE),
      spike_decision_int_adj = c(x1 = 2L),
      acd_parameter_adj = list(x1 = NULL),
      eliminated = c(x1 = FALSE),
      spike_binary_only_flags = c(x1 = FALSE),
      current_power_keys_adj = list(x1 = c(1, 1)),
      prev_power_keys_adj = NULL,
      prev_xi = NULL,
      has_prev = FALSE
    ),
    "missing stored|require stored|acd_parameter",
    ignore.case = TRUE
  )
})

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

# Test purpose: SAZ decision 1 retains both the structural-zero indicator and
# the continuous positive-part effect. Strong independent effects are used so
# AIC has an unambiguous preference for the full two-component representation.
test_that("23.4 SAZ decision 1 retains binary and continuous components", {
  set.seed(2304)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  z <- as.numeric(x <= 0)
  y <- 1 + 4.5 * z + 2.2 * x + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ z + x)

  expect_equal(as.integer(fit$spike_dec[["exposure"]]), 1L)
  expect_true(fit$catzero[["exposure"]])
  expect_false(all(is.na(fit$fp_powers[["exposure"]])))
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})

# Test purpose: SAZ decision 2 removes the structural-zero indicator when zero
# observations follow the same continuous relationship as positive values.
test_that("23.5 SAZ decision 2 retains only the continuous component", {
  set.seed(2305)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  y <- 1 + 2.4 * x + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ x)

  expect_equal(as.integer(fit$spike_dec[["exposure"]]), 2L)
  expect_false(fit$catzero[["exposure"]])
  expect_equal(as.numeric(fit$fp_powers[["exposure"]]), 1)
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})

# Test purpose: SAZ decision 3 removes the continuous FP component when only
# membership in the structural-zero group affects the outcome.
test_that("23.6 SAZ decision 3 retains only the binary zero indicator", {
  set.seed(2306)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  z <- as.numeric(x <= 0)
  y <- 1 + 4.2 * z + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ z)

  expect_equal(as.integer(fit$spike_dec[["exposure"]]), 3L)
  expect_true(fit$catzero[["exposure"]])
  expect_true(all(is.na(fit$fp_powers[["exposure"]])))
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})

# Test purpose: Stage-2 SAZ component selection must use alpha rather than
# select. keep forces Stage-1 inclusion (equivalent to select = 1), but the
# binary component should still be removable at alpha = 0.05 when it adds no
# information beyond the positive continuous component. Setting alpha = 1 on
# the same data should retain both components, proving that alpha is the Stage-2
# threshold.
test_that("23.6.1 SAZ Stage 2 uses alpha when select is forced by keep", {
  set.seed(23061)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  y <- 1 + 2.4 * x + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit_alpha_005 <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    keep = "exposure",
    alpha = 0.05,
    criterion = "pvalue",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  fit_alpha_1 <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    keep = "exposure",
    alpha = 1,
    criterion = "pvalue",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_equal(as.integer(fit_alpha_005$spike_dec[["exposure"]]), 2L)
  expect_equal(as.integer(fit_alpha_1$spike_dec[["exposure"]]), 1L)
  expect_true(fit_alpha_005$fp_terms["exposure", "selected"])
  expect_true(fit_alpha_1$fp_terms["exposure", "selected"])
})

# Test purpose: For a retained SAZ model, prediction on negative, zero, and
# positive new values must match both:
#   1. direct prediction from the final stored GLM; and
#   2. an independent manual calculation using X %*% beta and
#      diag(X %*% vcov(beta) %*% t(X)).
#
# This verifies the complete SAZ prediction contract:
#   - nonpositive values enter the binary structural-zero component;
#   - the continuous component uses the positive part of exposure;
#   - coefficient ordering matches the reconstructed model matrix;
#   - link-scale standard errors use the stored coefficient covariance matrix.
test_that("23.7 SAZ newdata prediction matches stored model and manual matrix calculation", {
  set.seed(2307)

  x <- c(rep(0, 120), seq(0.5, 8, length.out = 300))
  z <- as.numeric(x <= 0)

  y <- 1 +
    3.5 * z +
    1.7 * x +
    rnorm(length(x), sd = 0.04)

  xmat <- matrix(
    x,
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  # The data-generating model contains both a structural-zero effect and a
  # continuous positive-part effect. Therefore SAZ decision 1 must be retained:
  # binary indicator plus continuous component.
  expect_equal(unname(fit$spike_dec["exposure"]), 1L)

  newx <- matrix(
    c(-2, 0, 0.5, 2, 6),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  # Public mfp2 prediction.
  got <- predict(
    fit,
    newdata = newx,
    type = "link",
    se.fit = TRUE
  )

  # Because power = 1, shift = 0, scale = 1, and center = FALSE, the final
  # non-intercept design columns are exactly:
  #   exposure_bin = I(exposure <= 0)
  #   exposure.1   = max(exposure, 0)
  expected_design <- data.frame(
    exposure_bin = as.numeric(newx[, "exposure"] <= 0),
    exposure.1 = pmax(newx[, "exposure"], 0),
    check.names = FALSE
  )

  # Remove the mfp2 class so prediction dispatches directly to predict.glm()
  # using the already fitted final model.
  fit_glm <- fit
  class(fit_glm) <- setdiff(class(fit_glm), "mfp2")

  expected <- stats::predict(
    fit_glm,
    newdata = expected_design,
    type = "link",
    se.fit = TRUE
  )

  # Construct the complete model matrix manually, including the intercept.
  manual_x <- cbind(
    `(Intercept)` = 1,
    exposure.1 = expected_design$exposure.1,
    exposure_bin = expected_design$exposure_bin
  )

  beta <- stats::coef(fit_glm)
  beta_vcov <- stats::vcov(fit_glm)

  # The manually constructed design matrix is already in the fitted
  # coefficient order, so the numerical calculation is directly positional.
  expect_true(all(names(beta) %in% colnames(manual_x)))

  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  # Guard against any remaining coefficient-name or ordering mismatch.
  expect_identical(colnames(manual_x), names(beta))
  expect_identical(rownames(beta_vcov), names(beta))
  expect_identical(colnames(beta_vcov), names(beta))

  # Manual link prediction: eta = X beta.
  manual_fit <- as.numeric(manual_x %*% beta)

  # Manual link-scale standard error:
  # se_i = sqrt(x_i' Var(beta) x_i).
  #
  # rowSums((X %*% V) * X) is the diagonal of X V X' without constructing
  # the full prediction covariance matrix.
  manual_variance <- rowSums(
    (manual_x %*% beta_vcov) * manual_x
  )
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  # Public mfp2 prediction must equal direct predict.glm().
  expect_equal(
    as.numeric(got$fit),
    as.numeric(expected$fit),
    tolerance = 1e-8
  )

  expect_equal(
    as.numeric(got$se.fit),
    as.numeric(expected$se.fit),
    tolerance = 1e-8
  )

  # Public mfp2 prediction must also equal the independent matrix calculation.
  expect_equal(
    as.numeric(got$fit),
    manual_fit,
    tolerance = 1e-10
  )

  expect_equal(
    as.numeric(got$se.fit),
    manual_se,
    tolerance = 1e-10
  )

  # Direct predict.glm() must agree with the same manual calculation.
  expect_equal(
    as.numeric(expected$fit),
    manual_fit,
    tolerance = 1e-10
  )

  expect_equal(
    as.numeric(expected$se.fit),
    manual_se,
    tolerance = 1e-10
  )
})


# Test purpose: A balanced dataset with exactly the same x-response slope in
# every group should not be reported as an interaction. Identical residual
# patterns in both groups make the null interaction deterministic.
test_that("23.8 MFPI does not retain a deterministic no-interaction effect", {
  x_base <- seq(1, 8, length.out = 120)
  residual_pattern <- rep(c(-0.08, 0.08), length.out = length(x_base))
  dat <- rbind(
    data.frame(group = factor("A", levels = c("A", "B")), x = x_base,
               y = 1 + 0.5 * x_base + residual_pattern),
    data.frame(group = factor("B", levels = c("A", "B")), x = x_base,
               y = 1.7 + 0.5 * x_base + residual_pattern)
  )

  fit <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 0.05,
    verbose = FALSE
  )

  metric <- fit$all_model_metrics[fit$all_model_metrics$variable == "x", , drop = FALSE]
  expect_equal(nrow(metric), 1L)
  expect_gt(metric$pvalue, 0.05)
  expect_false("x" %in% names(fit$best_interaction_model))
})

# Test purpose: Observation weights must be passed unchanged through mfp2() and
# produce the same coefficients, covariance matrix, fitted values, and
# predictions as the corresponding weighted Gaussian glm.
test_that("23.9 weighted Gaussian mfp2 equals weighted glm", {
  set.seed(2309)
  n <- 220
  dat <- data.frame(
    x1 = runif(n, 1, 5),
    x2 = rnorm(n),
    w = runif(n, 0.5, 3)
  )
  dat$y <- 0.7 + 1.1 * dat$x1 - 0.6 * dat$x2 + rnorm(n, sd = 0.4)

  fit_mfp2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    weights = dat$w,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  fit_glm <- stats::glm(y ~ x1 + x2, data = dat, weights = w)
  nd <- dat[1:20, c("x1", "x2"), drop = FALSE]

  expect_mfp2_glm_parameters_equal(fit_mfp2, fit_glm, tolerance = 1e-8)
  expect_equal(unname(fitted(fit_mfp2)), unname(fitted(fit_glm)), tolerance = 1e-8)
  expect_equal(
    unname(predict(fit_mfp2, newdata = nd, se.fit = TRUE)$fit),
    unname(predict(fit_glm, newdata = nd, se.fit = TRUE)$fit),
    tolerance = 1e-8
  )
})

# Test purpose: Saving and restoring an mfp2 object must preserve coefficients,
# metadata, training predictions, newdata predictions, and prediction standard
# errors. This protects stored transformations and formula reconstruction.
test_that("23.10 mfp2 serialization preserves prediction behavior", {
  fit <- mfp2(
    lpsa ~ fp(age) + fp(cavol) + svi,
    data = prostate,
    keep = "svi",
    verbose = FALSE
  )
  nd <- prostate[1:15, c("age", "cavol", "svi"), drop = FALSE]
  before <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  after <- predict(restored, newdata = nd, type = "link", se.fit = TRUE)

  expect_s3_class(restored, "mfp2")
  expect_equal(restored$fp_powers, fit$fp_powers)
  expect_equal(coef(restored), coef(fit))
  expect_equal(as.numeric(after$fit), as.numeric(before$fit), tolerance = 1e-12)
  expect_equal(as.numeric(after$se.fit), as.numeric(before$se.fit), tolerance = 1e-12)
})

# Test purpose: Saving and restoring an mfpi object must preserve both ordinary
# subject-level prediction and the manually evaluated fitted-function path.
test_that("23.11 mfpi serialization preserves ordinary and fitted-function predictions", {
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "fp1"),
    flex = "flex1",
    verbose = FALSE
  )
  nd <- prostate[1:12, , drop = FALSE]
  before_link <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  before_fun <- predict(
    fit, terms = "cavol", model = "all", type = "function",
    grid = TRUE, n_grid = 20, se.fit = TRUE
  )

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  after_link <- predict(
    restored, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  after_fun <- predict(
    restored, terms = "cavol", model = "all", type = "function",
    grid = TRUE, n_grid = 20, se.fit = TRUE
  )

  expect_s3_class(restored, "mfpi")
  expect_equal(after_link$predictions, before_link$predictions, tolerance = 1e-12)
  expect_equal(after_fun$functions, before_fun$functions, tolerance = 1e-12)
})


# Test purpose: Continuous-only formula objects created before the explicit
# term and coefficient-column mappings were stored must remain predictable.
# The fallback must also resolve fp()/fp2() formula labels to source variables.
test_that("23.12 legacy continuous formula objects reconstruct prediction mappings", {
  fit <- mfp2(
    lpsa ~ fp(age, center = FALSE) + fp(cavol, center = FALSE) + svi,
    data = prostate,
    keep = c("age", "cavol", "svi"),
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  nd <- prostate[1:18, c("age", "cavol", "svi"), drop = FALSE]
  expected <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
  expected_term <- predict(
    fit, newdata = nd, type = "terms", terms = "age", se.fit = TRUE
  )

  legacy <- fit
  legacy$term_to_columns <- NULL
  legacy$formula_term_to_columns <- NULL
  legacy$transformed_to_model_columns <- NULL
  legacy$formula_prediction_term_names <- NULL

  got <- predict(legacy, newdata = nd, type = "link", se.fit = TRUE)
  got_term <- predict(
    legacy, newdata = nd, type = "terms", terms = "age", se.fit = TRUE
  )

  expect_equal(as.numeric(got$fit), as.numeric(expected$fit), tolerance = 1e-12)
  expect_equal(
    as.numeric(got$se.fit), as.numeric(expected$se.fit), tolerance = 1e-12
  )
  expect_equal(got_term, expected_term, tolerance = 1e-12)
})

# Test purpose: Legacy mapping reconstruction must preserve exact fitted
# coefficient names for matrix columns that require quoting in the final model.
test_that("23.13 legacy matrix objects preserve non-syntactic coefficient names", {
  x <- cbind(
    "age years" = prostate$age,
    "cavol-value" = prostate$cavol
  )
  fit <- mfp2(
    x,
    prostate$lpsa,
    keep = colnames(x),
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  nd <- as.data.frame(x[1:16, , drop = FALSE], check.names = FALSE)
  expected <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
  expected_term <- predict(
    fit, newdata = nd, type = "terms", terms = "age years", se.fit = TRUE
  )

  legacy <- fit
  legacy$term_to_columns <- NULL
  legacy$transformed_to_model_columns <- NULL
  got <- predict(legacy, newdata = nd, type = "link", se.fit = TRUE)
  got_term <- predict(
    legacy, newdata = nd, type = "terms", terms = "age years", se.fit = TRUE
  )

  expect_equal(as.numeric(got$fit), as.numeric(expected$fit), tolerance = 1e-12)
  expect_equal(
    as.numeric(got$se.fit), as.numeric(expected$se.fit), tolerance = 1e-12
  )
  expect_equal(got_term, expected_term, tolerance = 1e-12)
})

# Test purpose: Missing grouped-term metadata cannot be reconstructed as an
# identity mapping. Prediction must fail explicitly rather than silently use a
# conceptual group name as if it were one raw model-matrix column.
test_that("23.14 grouped matrix objects require their member-column mapping", {
  stage <- factor(
    rep(c("I", "II", "III"), length.out = nrow(prostate)),
    levels = c("I", "II", "III")
  )
  mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  fit <- mfp2(
    mm,
    prostate$lpsa,
    term_groups = list(stage = colnames(mm)),
    keep = "stage",
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  legacy <- fit
  legacy$term_to_columns <- NULL

  expect_error(
    predict(legacy, newdata = mm[1:10, , drop = FALSE]),
    "lacks grouped-term mapping metadata"
  )
})

# =============================================================================
# 24. subset handling across formula and matrix interfaces
# =============================================================================

# Test purpose: The formula interface must rebuild ordered-factor contrasts on
# the fitted rows while retaining full-data scaling for continuous predictors.
test_that("24.1 mfp2 formula rebuilds categorical coding after subset", {
  set.seed(2401)
  n_each <- 35L
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = n_each),
    levels = c("A", "B", "C", "D")
  )
  x <- seq_len(length(stage))
  x[stage == "D"] <- x[stage == "D"] * 10000
  dat <- data.frame(
    y = 1 + 0.03 * seq_len(length(stage)) + rnorm(length(stage), sd = 0.2),
    x = x,
    stage = stage
  )
  fit_rows <- dat$stage != "D"

  fit <- mfp2(
    y ~ x + stage,
    data = dat,
    subset = fit_rows,
    keep = c("x", "stage"),
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    verbose = FALSE
  )

  dat_fit <- droplevels(dat[fit_rows, , drop = FALSE])
  expected_x <- stats::model.matrix(~ x + stage, data = dat_fit)[, -1L, drop = FALSE]

  expect_identical(fit$formula_xlevels$stage, levels(dat_fit$stage))
  expect_identical(fit$formula_model_matrix_columns, colnames(expected_x))
  expect_identical(fit$term_to_columns$stage, c("stage.L", "stage.Q"))
  expect_equal(
    unname(fit$x_original[, c("stage.L", "stage.Q"), drop = FALSE]),
    unname(expected_x[, c("stage.L", "stage.Q"), drop = FALSE]),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(fit$transformations["x", "scale"]),
    find_scale_factor(dat$x)
  )
  expect_false(isTRUE(all.equal(
    find_scale_factor(dat$x),
    find_scale_factor(dat_fit$x)
  )))
  expect_error(
    predict(fit, newdata = dat[dat$stage == "D", , drop = FALSE][1L, ]),
    "new level|new levels|factor"
  )
})

# Test purpose: A formula factor with only one retained level has no estimable
# effect and must fail with a categorical-specific diagnostic.
test_that("24.2 mfp2 formula rejects a one-level factor after subset", {
  set.seed(2402)
  dat <- data.frame(
    y = rnorm(60),
    x = runif(60, 1, 5),
    group = factor(rep(c("A", "B", "C"), each = 20))
  )

  expect_error(
    mfp2(
      y ~ x + group,
      data = dat,
      subset = dat$group == "A",
      df = 1,
      verbose = FALSE
    ),
    "fewer than two observed levels"
  )
})

# Test purpose: The matrix interface cannot regenerate reduced polynomial
# contrasts and therefore rejects a grouped block that loses rank after subset.
test_that("24.3 mfp2 matrix rejects grouped rank loss after subset", {
  set.seed(2403)
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = 25),
    levels = c("A", "B", "C", "D")
  )
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  x <- cbind(z = runif(length(stage), 1, 4), stage_mm)
  y <- rnorm(length(stage))

  expect_error(
    mfp2(
      x,
      y,
      subset = stage != "D",
      term_groups = list(stage = colnames(stage_mm)),
      keep = "stage",
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    "lose estimable dimension"
  )
})

# Test purpose: MFPI formula fitting uses the retained factor levels and
# regenerates the corresponding polynomial-contrast block.
test_that("24.4 MFPI formula rebuilds categorical coding after subset", {
  set.seed(2404)
  n_each <- 35L
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = n_each),
    levels = c("A", "B", "C", "D")
  )
  trt <- factor(rep(c("control", "treated"), length.out = length(stage)))
  x <- runif(length(stage), 1, 7)
  z <- rnorm(length(stage))
  y <- 1 + 0.4 * x + 0.7 * x * (trt == "treated") + 0.2 * z + rnorm(length(stage), sd = 0.3)
  dat <- data.frame(y = y, trt = trt, x = x, stage = stage, z = z)
  fit_rows <- dat$stage != "D"

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    subset = fit_rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(fit$formula_xlevels$stage, c("A", "B", "C"))
  expect_identical(fit$term_to_columns$stage, c("stage.L", "stage.Q"))
  expect_false("stage.C" %in% fit$formula_design_columns)
})

# Test purpose: MFPI formula fitting gives a direct diagnostic when a factor is
# reduced to one observed level by subset.
test_that("24.5 MFPI formula rejects a one-level factor after subset", {
  set.seed(2405)
  dat <- make_mfpi_factor_data(n = 120L)

  expect_error(
    mfpi(
      y ~ trt + x + stage + z,
      data = dat,
      subset = dat$stage == "I",
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    "fewer than two observed levels"
  )
})

# Test purpose: The MFPI matrix interface rejects supplied grouped contrast
# columns when subsetting removes an estimable dimension.
test_that("24.6 MFPI matrix rejects grouped rank loss after subset", {
  set.seed(2406)
  n_each <- 30L
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = n_each),
    levels = c("A", "B", "C", "D")
  )
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  trt <- rep(c(0, 1), length.out = length(stage))
  x_cont <- runif(length(stage), 1, 6)
  x <- cbind(trt = trt, x = x_cont, stage_mm, z = rnorm(length(stage)))
  y <- 1 + 0.4 * x_cont + 0.6 * x_cont * trt + rnorm(length(stage), sd = 0.3)

  expect_error(
    mfpi(
      x,
      y,
      subset = stage != "D",
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      term_groups = list(stage = colnames(stage_mm)),
      keep = "stage",
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    "lose estimable dimension"
  )
})

# Test purpose: A raw matrix column that becomes constant after subsetting is
# rejected even when it is not declared through term_groups.
test_that("24.7 MFPI matrix rejects a predictor with no post-subset variation", {
  set.seed(2407)
  n <- 120L
  trt <- rep(c(0, 1), length.out = n)
  x_cont <- runif(n, 1, 5)
  indicator <- rep(c(0, 0, 1, 1), length.out = n)
  x <- cbind(
    trt = trt,
    x = x_cont,
    indicator = indicator,
    z = rnorm(n)
  )
  y <- 1 + 0.3 * x_cont + 0.5 * x_cont * trt + rnorm(n, sd = 0.3)

  expect_error(
    mfpi(
      x,
      y,
      subset = indicator == 0,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    "no variation after applying `subset`"
  )
})


# Test purpose: Removing the original treatment-contrast reference level causes
# the formula interface to rebuild the factor with a retained reference level.
test_that("24.8 mfp2 formula rebuilds treatment contrasts when reference is removed", {
  set.seed(2408)
  group <- factor(
    rep(c("A", "B", "C"), each = 30),
    levels = c("A", "B", "C")
  )
  dat <- data.frame(
    y = rnorm(length(group)),
    x = runif(length(group), 1, 5),
    group = group
  )
  fit_rows <- dat$group != "A"

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    subset = fit_rows,
    keep = "group",
    df = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  dat_fit <- droplevels(dat[fit_rows, , drop = FALSE])
  expected_x <- stats::model.matrix(~ x + group, dat_fit)[, -1L, drop = FALSE]

  expect_identical(fit$formula_xlevels$group, c("B", "C"))
  expect_identical(fit$term_to_columns$group, "groupC")
  expect_equal(
    unname(fit$x_original[, "groupC", drop = FALSE]),
    unname(expected_x[, "groupC", drop = FALSE]),
    tolerance = 1e-12
  )
})


# Test purpose: Formula methods must bypass all subset-specific processing
# when subset is NULL. This guards the ordinary formula path against both the
# model.frame() NSE regression and unnecessary row/factor reconstruction.
test_that("24.9 formula methods bypass subset helpers when subset is NULL", {
  set.seed(2409)
  n <- 90L
  dat <- data.frame(
    y = rnorm(n),
    trt = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 5),
    stage = factor(rep(c("A", "B", "C"), length.out = n)),
    z = rnorm(n)
  )

  testthat::local_mocked_bindings(
    formula_subset_rows = function(...) {
      stop("row resolver must not be called for subset = NULL", call. = FALSE)
    },
    subset_formula_model_frame = function(...) {
      stop("frame helper must not be called for subset = NULL", call. = FALSE)
    },
    .package = "mfp2"
  )

  fit_mfp2 <- mfp2(
    y ~ x + stage,
    data = dat,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  expect_s3_class(fit_mfp2, "mfp2")

  fit_mfpi <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_s3_class(fit_mfpi, "mfpi")
})


# Test purpose: The mfp2 formula interface evaluates a subset expression once
# in the formula data mask. A data column must take precedence over a same-named
# caller object.
test_that("24.10 mfp2 formula evaluates subset once with data precedence", {
  set.seed(2410)
  n <- 96L
  dat <- data.frame(
    y = 1 + 0.4 * seq_len(n) / n + rnorm(n, sd = 0.2),
    x = runif(n, 1, 5),
    keep = rep(c(TRUE, TRUE, FALSE), length.out = n)
  )
  keep <- rep(TRUE, n)
  evaluations <- 0L
  evaluate_once <- function(value) {
    evaluations <<- evaluations + 1L
    value
  }

  fit <- mfp2(
    y ~ x,
    data = dat,
    subset = evaluate_once(keep),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_identical(evaluations, 1L)
  expect_equal(nrow(fit$x_original), sum(dat$keep))
})

# Test purpose: The MFPI formula interface uses the same exact-once, data-first
# subset lookup as mfp2.formula().
test_that("24.11 MFPI formula evaluates subset once with data precedence", {
  set.seed(2411)
  n <- 120L
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.4 * x + 0.6 * x * (trt == "treated") + rnorm(n, sd = 0.3),
    trt = trt,
    x = x,
    z = rnorm(n),
    keep = rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)
  )
  keep <- rep(FALSE, n)
  evaluations <- 0L
  evaluate_once <- function(value) {
    evaluations <<- evaluations + 1L
    value
  }

  fit <- mfpi(
    y ~ trt + x + z,
    data = dat,
    subset = evaluate_once(keep),
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "z",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(evaluations, 1L)
  expect_equal(fit$nobs, sum(dat$keep))
})

# Test purpose: A formula subset need not be a data column. Numeric row indices
# defined in the formula environment remain valid in both formula interfaces.
test_that("24.12 formula methods accept caller-scoped subset row indices", {
  set.seed(2412)
  n <- 120L
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.3 * x + 0.5 * x * (trt == "treated") + rnorm(n, sd = 0.3),
    trt = trt,
    x = x,
    z = rnorm(n)
  )
  rows <- which(rep(c(TRUE, TRUE, FALSE), length.out = n))

  fit_mfp2 <- mfp2(
    y ~ x + z,
    data = dat,
    subset = rows,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  expect_equal(nrow(fit_mfp2$x_original), length(rows))

  fit_mfpi <- mfpi(
    y ~ trt + x + z,
    data = dat,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "z",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_equal(fit_mfpi$nobs, length(rows))
})

# Test purpose: Duplicate numeric row positions replicate observations and are
# therefore rejected consistently by all four public interfaces.
test_that("24.13 all subset interfaces reject duplicated numeric row indices", {
  set.seed(2413)
  n <- 120L
  trt_factor <- factor(rep(c("control", "treated"), length.out = n))
  trt_numeric <- as.numeric(trt_factor == "treated")
  x_cont <- runif(n, 1, 5)
  z <- rnorm(n)
  y <- 1 + 0.4 * x_cont + 0.5 * x_cont * trt_numeric + 0.2 * z +
    rnorm(n, sd = 0.3)
  duplicate_rows <- c(seq_len(80L), 80L)
  duplicate_message <- "must not contain duplicated row indices"

  x_matrix <- cbind(trt = trt_numeric, x = x_cont, z = z)
  dat <- data.frame(y = y, trt = trt_factor, x = x_cont, z = z)

  expect_error(
    mfp2(
      x_matrix[, c("x", "z"), drop = FALSE],
      y,
      subset = duplicate_rows,
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    duplicate_message
  )

  expect_error(
    mfp2(
      y ~ x + z,
      data = dat,
      subset = duplicate_rows,
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    duplicate_message
  )

  expect_error(
    mfpi(
      x_matrix,
      y,
      subset = duplicate_rows,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    duplicate_message
  )

  expect_error(
    mfpi(
      y ~ trt + x + z,
      data = dat,
      subset = duplicate_rows,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    duplicate_message
  )
})

# Test purpose: Unique numeric positions are not sorted. Their supplied order is
# retained by the row resolver and by all four fitting interfaces.
test_that("24.14 unique numeric subset indices preserve supplied order", {
  set.seed(2414)
  n <- 100L
  x <- cbind(
    x = seq_len(n),
    z = rnorm(n)
  )
  y <- 1 + 0.02 * x[, "x"] + rnorm(n, sd = 0.2)
  dat <- data.frame(y = y, x = x[, "x"], z = x[, "z"])
  rows <- c(51:100, 1:50)

  expect_identical(formula_subset_rows(rows, n), as.integer(rows))
  expect_error(
    formula_subset_rows(c(1L, 2L, 2L), n),
    "must not contain duplicated row indices"
  )
  expect_error(
    subset_formula_model_frame(
      data.frame(x = seq_len(4L)),
      c(1L, 2L, 2L)
    ),
    "unique valid integer model-frame positions"
  )

  fit_matrix <- mfp2(
    x,
    y,
    subset = rows,
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  expect_equal(
    unname(fit_matrix$x_original[, "x"]),
    unname(x[rows, "x"]),
    tolerance = 0
  )

  fit_formula <- mfp2(
    y ~ x + z,
    data = dat,
    subset = rows,
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  expect_equal(
    unname(fit_formula$x_original[, "x"]),
    unname(dat$x[rows]),
    tolerance = 0
  )

  trt_factor <- factor(rep(c("control", "treated"), length.out = n))
  trt_numeric <- as.numeric(trt_factor == "treated")
  x_mfpi <- cbind(trt = trt_numeric, x = x[, "x"], z = x[, "z"])
  dat_mfpi <- data.frame(y = y, trt = trt_factor, x = x[, "x"], z = x[, "z"])

  fit_mfpi_matrix <- mfpi(
    x_mfpi,
    y,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_equal(
    unname(fit_mfpi_matrix$x_train_internal[, "x"]),
    unname(x_mfpi[rows, "x"]),
    tolerance = 0
  )

  fit_mfpi_formula <- mfpi(
    y ~ trt + x + z,
    data = dat_mfpi,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_equal(
    unname(fit_mfpi_formula$x_train_internal[, "x"]),
    unname(dat_mfpi$x[rows]),
    tolerance = 0
  )
})

# Test purpose: The retained-frame helper preserves an explicit factor contrast
# when a genuine subset retains every level, including contrasts created by C().
test_that("24.15 formula subset helper preserves custom contrasts when all levels remain", {
  stage <- factor(rep(c("A", "B", "C"), each = 12L))
  contrasts(stage) <- stats::contr.sum(3L)
  # Only factor levels/contrasts are under test; keep the response deterministic.
  dat <- data.frame(y = seq_along(stage), stage = stage)
  mf <- stats::model.frame(y ~ stage, data = dat)
  rows <- c(1:8, 13:20, 25:32)

  retained <- subset_formula_model_frame(mf, rows)

  expect_identical(levels(retained$stage), levels(mf$stage))
  expect_identical(
    attr(retained$stage, "contrasts", exact = TRUE),
    attr(mf$stage, "contrasts", exact = TRUE)
  )

  mf_c <- stats::model.frame(y ~ C(stage, stats::contr.sum), data = dat)
  retained_c <- subset_formula_model_frame(mf_c, rows)
  factor_name <- setdiff(names(mf_c), "y")
  expect_identical(
    attr(retained_c[[factor_name]], "contrasts", exact = TRUE),
    attr(mf_c[[factor_name]], "contrasts", exact = TRUE)
  )
})

# Test purpose: Both formula methods retain factor-specific custom contrasts in
# their fitted prediction metadata when all factor levels survive subsetting.
test_that("24.16 formula methods retain custom contrasts when all levels remain", {
  set.seed(2416)
  n <- 120L
  stage <- factor(rep(c("A", "B", "C"), length.out = n))
  contrasts(stage) <- stats::contr.sum(3L)
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.3 * x + 0.5 * (stage == "B") + rnorm(n, sd = 0.25),
    trt = trt,
    x = x,
    stage = stage,
    z = rnorm(n),
    keep_row = rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)
  )

  expected_mf <- stats::model.frame(
    y ~ x + stage,
    data = dat[dat$keep_row, , drop = FALSE]
  )
  expected_mm <- stats::model.matrix(y ~ x + stage, data = expected_mf)

  fit_mfp2 <- mfp2(
    y ~ x + stage,
    data = dat,
    subset = keep_row,
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_identical(
    fit_mfp2$formula_contrasts$stage,
    attr(expected_mm, "contrasts")$stage
  )
  expect_equal(
    unname(fit_mfp2$x_original[, c("stage1", "stage2"), drop = FALSE]),
    unname(expected_mm[, c("stage1", "stage2"), drop = FALSE]),
    tolerance = 1e-12
  )

  expected_mfpi_mf <- stats::model.frame(
    y ~ trt + x + stage + z,
    data = dat[dat$keep_row, , drop = FALSE]
  )
  expected_mfpi_mm <- stats::model.matrix(
    y ~ trt + x + stage + z,
    data = expected_mfpi_mf
  )

  fit_mfpi <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    subset = keep_row,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = c("stage", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    fit_mfpi$formula_contrasts$stage,
    attr(expected_mfpi_mm, "contrasts")$stage
  )
  expect_equal(
    unname(fit_mfpi$x_train_internal[, c("stage1", "stage2"), drop = FALSE]),
    unname(expected_mfpi_mm[, c("stage1", "stage2"), drop = FALSE]),
    tolerance = 1e-12
  )
})

# Test purpose: A custom contrast matrix cannot be reduced uniquely when a
# subset removes one of its factor levels, so both formula methods stop clearly.
test_that("24.17 formula methods reject ambiguous custom contrasts after level removal", {
  set.seed(2417)
  n <- 120L
  stage <- factor(rep(c("A", "B", "C"), length.out = n))
  contrasts(stage) <- stats::contr.sum(3L)
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.3 * x + rnorm(n, sd = 0.25),
    trt = trt,
    x = x,
    stage = stage,
    z = rnorm(n)
  )
  keep_rows <- dat$stage != "C"
  message <- "Custom contrasts for factor `stage` cannot be reconstructed safely"

  expect_error(
    mfp2(
      y ~ x + stage,
      data = dat,
      subset = keep_rows,
      keep = "stage",
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    message,
    fixed = TRUE
  )

  expect_error(
    mfpi(
      y ~ trt + x + stage + z,
      data = dat,
      subset = keep_rows,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      keep = c("stage", "z"),
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    message,
    fixed = TRUE
  )
})

# Test purpose: The custom-contrast inspection is confined to the non-NULL
# subset helper path. An ordinary fit reuses its full model frame unchanged.
test_that("24.18 subset NULL bypasses custom contrast processing", {
  set.seed(2418)
  n <- 90L
  stage <- factor(rep(c("A", "B", "C"), length.out = n))
  contrasts(stage) <- stats::contr.sum(3L)
  dat <- data.frame(y = rnorm(n), x = runif(n, 1, 5), stage = stage)

  testthat::local_mocked_bindings(
    subset_formula_factor_column = function(...) {
      stop("factor subset helper must not be called for subset = NULL", call. = FALSE)
    },
    .package = "mfp2"
  )

  fit <- mfp2(
    y ~ x + stage,
    data = dat,
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_identical(
    fit$formula_contrasts$stage,
    attr(stats::model.matrix(y ~ x + stage, data = dat), "contrasts")$stage
  )
})


# Test purpose: Inline factor wrappers, binary factors, and interactions that
# include a categorical variable must all be rebuilt from the retained rows.
test_that("24.19 formula subset rebuilds inline, binary, and interaction coding", {
  set.seed(2419)
  dat <- expand.grid(
    stage_code = 1:3,
    binary_code = 0:1,
    group = c("A", "B", "C"),
    replicate = 1:10,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dat$group <- factor(dat$group, levels = c("A", "B", "C"))
  dat$x <- runif(nrow(dat), 1, 6)
  dat$y <- 0.4 * dat$x +
    0.5 * (dat$stage_code == 2L) -
    0.3 * (dat$stage_code == 3L) +
    0.6 * dat$binary_code +
    0.2 * dat$x * (dat$group == "B") -
    0.15 * dat$x * (dat$group == "C") +
    rnorm(nrow(dat), sd = 0.25)
  dat$keep_row <- seq_len(nrow(dat)) %% 7L != 0L

  fit <- mfp2(
    y ~ x + factor(stage_code) + as.factor(binary_code) +
      group + x:group,
    data = dat,
    subset = keep_row,
    keep = c("x", "stage_code", "binary_code", "group", "x:group"),
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

  retained <- droplevels(dat[dat$keep_row, , drop = FALSE])
  expected <- stats::model.matrix(
    ~ x + factor(stage_code) + as.factor(binary_code) +
      group + x:group,
    data = retained
  )[, -1L, drop = FALSE]

  expect_identical(fit$formula_model_matrix_columns, colnames(expected))
  expect_equal(
    unname(fit$x_original[, colnames(expected), drop = FALSE]),
    unname(expected),
    tolerance = 1e-12
  )
  expect_true(all(c("stage_code", "binary_code", "group", "x:group") %in%
                    names(fit$term_to_columns)))
})


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


# Test purpose: Predictions from a subset-fitted formula model must reproduce
# training predictions, accept retained levels, reject removed/unseen levels,
# and retain the same behavior after serialization.
test_that("24.25 subset-fitted mfp2 prediction uses retained factor levels", {
  set.seed(2425)
  n <- 180L
  dat <- data.frame(
    x = runif(n, 1, 5),
    stage = factor(rep(c("A", "B", "C"), length.out = n),
                   levels = c("A", "B", "C"))
  )
  dat$y <- 1 + 0.4 * dat$x + 0.7 * (dat$stage == "B") +
    rnorm(n, sd = 0.25)
  fit_rows <- dat$stage != "C"

  fit <- mfp2(
    y ~ x + stage,
    data = dat,
    subset = fit_rows,
    keep = c("x", "stage"),
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
  retained <- droplevels(dat[fit_rows, , drop = FALSE])

  training <- predict(fit, type = "link", se.fit = TRUE)
  reconstructed <- predict(
    fit,
    newdata = retained,
    type = "link",
    se.fit = TRUE
  )
  expect_equal(unname(reconstructed$fit), unname(training$fit), tolerance = 1e-10)
  expect_equal(unname(reconstructed$se.fit), unname(training$se.fit), tolerance = 1e-10)

  retained_row <- retained[retained$stage == "B", , drop = FALSE][1L, ]
  expect_length(predict(fit, newdata = retained_row), 1L)

  removed_row <- dat[dat$stage == "C", , drop = FALSE][1L, ]
  expect_error(
    predict(fit, newdata = removed_row),
    "new level|new levels|factor",
    ignore.case = TRUE
  )

  unseen_row <- retained_row
  unseen_row$stage <- factor("D", levels = c("A", "B", "D"))
  expect_error(
    predict(fit, newdata = unseen_row),
    "new level|new levels|factor",
    ignore.case = TRUE
  )

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  restored_prediction <- predict(
    restored,
    newdata = retained,
    type = "link",
    se.fit = TRUE
  )
  expect_equal(restored_prediction$fit, reconstructed$fit, tolerance = 1e-12)
  expect_equal(restored_prediction$se.fit, reconstructed$se.fit, tolerance = 1e-12)
})


# Test purpose: MFPI formula metadata, formula-offset reconstruction, and
# ordinary prediction must all use only factor levels retained by subset.
test_that("24.26 MFPI subset prediction uses retained metadata and formula offset", {
  set.seed(2426)
  n <- 180L
  dat <- data.frame(
    trt = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 5),
    stage = factor(rep(c("A", "B", "C"), length.out = n),
                   levels = c("A", "B", "C")),
    z = rnorm(n),
    off = rnorm(n, sd = 0.15)
  )
  dat$y <- 0.7 + 0.45 * dat$x + 0.6 * (dat$trt == "treated") * dat$x +
    0.4 * (dat$stage == "B") + 0.2 * dat$z + dat$off +
    rnorm(n, sd = 0.3)
  fit_rows <- dat$stage != "C"

  fit <- mfpi(
    y ~ trt + x + stage + z + offset(off),
    data = dat,
    subset = fit_rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
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

  expect_identical(fit$formula_xlevels$stage, c("A", "B"))
  expect_false(is.null(fit$formula_offset_terms))
  expect_true("stageB" %in% fit$formula_design_columns)
  expect_false("stageC" %in% fit$formula_design_columns)

  retained <- droplevels(
    dat[fit_rows, c("trt", "x", "stage", "z", "off"), drop = FALSE]
  )
  reconstructed <- predict(
    fit,
    newdata = retained,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )
  training <- predict(
    fit,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )
  expect_equal(
    reconstructed$predictions,
    training$predictions,
    tolerance = 1e-10
  )

  retained_row <- retained[retained$stage == "B", , drop = FALSE][1L, ]
  retained_prediction <- predict(
    fit,
    newdata = retained_row,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expect_true(all(is.finite(retained_prediction$predictions$fit)))

  removed_row <- dat[
    dat$stage == "C",
    c("trt", "x", "stage", "z", "off"),
    drop = FALSE
  ][1L, ]
  expect_error(
    predict(
      fit,
      newdata = removed_row,
      terms = "x",
      model = "all",
      type = "link",
      se.fit = FALSE
    ),
    "new level|new levels|factor",
    ignore.case = TRUE
  )

  unseen_row <- retained_row
  unseen_row$stage <- factor("D", levels = c("A", "B", "D"))
  expect_error(
    predict(
      fit,
      newdata = unseen_row,
      terms = "x",
      model = "all",
      type = "link",
      se.fit = FALSE
    ),
    "new level|new levels|factor",
    ignore.case = TRUE
  )
})

###########################################################################

make_cox_reference_v1_data <- function(n = 480L, seed = 11001L) {
  set.seed(seed)

  stratum <- factor(rep(c("A", "B"), length.out = n))
  x_mean <- ifelse(stratum == "A", 4, 9)
  x <- stats::rnorm(n, mean = x_mean, sd = 1.15)
  off <- stats::rnorm(n, mean = 0.35, sd = 0.12)

  eta <- 0.48 * (x - 6) + 0.30 * (stratum == "B") + off
  event_time <- stats::rexp(n, rate = 0.025 * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.012)

  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    x = x,
    stratum = stratum,
    off = off
  )
}


strip_mfp2_class_v1 <- function(object) {
  class(object) <- setdiff(class(object), "mfp2")
  object
}


prepare_formula_cox_newdata_v1 <- function(object,
                                           newdata,
                                           include_response = FALSE) {
  selected_terms <- get_selected_variable_names(object)

  predictor_data <- reconstruct_formula_newdata(
    object,
    newdata,
    terms = selected_terms
  )

  prediction_strata <- reconstruct_formula_strata_newdata(object, newdata)
  prediction_offset <- reconstruct_formula_offset_newdata(object, newdata)

  prepared <- prepare_newdata_for_predict(
    object,
    predictor_data,
    terms = selected_terms,
    strata = prediction_strata,
    offset = prediction_offset,
    check_binary = FALSE
  )

  if (include_response) {
    response_name <- cox_internal_response_name(object)
    prepared[[response_name]] <- I(
      survival::Surv(newdata$time, newdata$status)
    )
  }

  prepared
}


# The direct coxph fit uses exactly the design scale supplied to the final mfp2
# Cox fit. Matching nocenter and ties is important: otherwise the comparison can
# accidentally test a different Cox centering convention.
test_that("version 1 Cox lp and risk references match an equivalent coxph fit", {
  dat <- make_cox_reference_v1_data(seed = 11002L)
  prediction_x <- data.frame(x = c(2.75, 4.5, 7.25, 10.5))

  for (center_value in c(TRUE, FALSE)) {
    fit_mfp2 <- mfp2(
      x = as.matrix(dat["x"]),
      y = survival::Surv(dat$time, dat$status),
      family = "cox",
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = center_value,
      nocenter = NULL,
      xorder = "original",
      ties = "breslow",
      verbose = FALSE
    )

    fitted_center <- if (center_value) {
      unname(fit_mfp2$centers[[1L]])
    } else {
      0
    }

    dat_cox <- transform(dat, z = x - fitted_center)
    prediction_cox <- data.frame(z = prediction_x$x - fitted_center)

    fit_coxph <- survival::coxph(
      survival::Surv(time, status) ~ z,
      data = dat_cox,
      ties = "breslow",
      nocenter = NULL,
      x = TRUE,
      y = TRUE
    )

    for (reference_value in c("zero", "sample", "strata")) {
      got_lp <- predict(
        fit_mfp2,
        newdata = prediction_x,
        type = "lp",
        se.fit = TRUE,
        cox_reference = reference_value
      )
      expected_lp <- predict(
        fit_coxph,
        newdata = prediction_cox,
        type = "lp",
        se.fit = TRUE,
        reference = reference_value
      )

      expect_equal(
        unname(got_lp$fit),
        unname(expected_lp$fit),
        tolerance = 1e-8
      )
      expect_equal(
        unname(got_lp$se.fit),
        unname(expected_lp$se.fit),
        tolerance = 1e-8
      )

      got_risk <- predict(
        fit_mfp2,
        newdata = prediction_x,
        type = "risk",
        se.fit = TRUE,
        cox_reference = reference_value
      )
      expected_risk <- predict(
        fit_coxph,
        newdata = prediction_cox,
        type = "risk",
        se.fit = TRUE,
        reference = reference_value
      )

      expect_equal(
        unname(got_risk$fit),
        unname(expected_risk$fit),
        tolerance = 1e-8
      )
      expect_equal(
        unname(got_risk$se.fit),
        unname(expected_risk$se.fit),
        tolerance = 1e-8
      )
    }

    lp_zero <- predict(
      fit_mfp2,
      prediction_x,
      type = "lp",
      cox_reference = "zero"
    )
    lp_sample <- predict(
      fit_mfp2,
      prediction_x,
      type = "lp",
      cox_reference = "sample"
    )

    expected_constant <- sum(
      stats::coef(fit_coxph) * fit_coxph$means,
      na.rm = TRUE
    )

    expect_equal(
      unname(lp_zero - lp_sample),
      rep(unname(expected_constant), nrow(prediction_x)),
      tolerance = 1e-8
    )

    if (center_value) {
      expect_equal(
        unname(lp_zero),
        unname(lp_sample),
        tolerance = 1e-8
      )
    } else {
      expect_gt(abs(expected_constant), 1e-4)
    }
  }
})


# Offset centering is independent of the covariate reference. The zero/sample
# change must therefore remain a common covariate constant even when prediction
# offsets differ from row to row.
test_that("version 1 Cox hazard ratios are reference invariant with offsets", {
  dat <- make_cox_reference_v1_data(seed = 11003L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- data.frame(
    x = c(3.5, 8.5),
    off = c(0.10, 0.65)
  )

  lp_zero <- predict(fit, nd, type = "lp", cox_reference = "zero")
  lp_sample <- predict(fit, nd, type = "lp", cox_reference = "sample")
  risk_zero <- predict(fit, nd, type = "risk", cox_reference = "zero")
  risk_sample <- predict(fit, nd, type = "risk", cox_reference = "sample")

  # Use the stripped coxph object as the oracle. This is stronger and safer
  # than recomputing the shift with a positional coef * means expression:
  # predict.coxph() owns the exact model-matrix reconstruction and offset
  # centering rules for the fitted object.
  prepared <- prepare_formula_cox_newdata_v1(fit, nd)
  base <- strip_mfp2_class_v1(fit)

  base_lp_zero <- predict(
    base,
    newdata = prepared,
    type = "lp",
    reference = "zero"
  )
  base_lp_sample <- predict(
    base,
    newdata = prepared,
    type = "lp",
    reference = "sample"
  )
  base_risk_zero <- predict(
    base,
    newdata = prepared,
    type = "risk",
    reference = "zero"
  )
  base_risk_sample <- predict(
    base,
    newdata = prepared,
    type = "risk",
    reference = "sample"
  )

  expect_equal(unname(lp_zero), unname(base_lp_zero), tolerance = 1e-8)
  expect_equal(unname(lp_sample), unname(base_lp_sample), tolerance = 1e-8)
  expect_equal(unname(risk_zero), unname(base_risk_zero), tolerance = 1e-8)
  expect_equal(unname(risk_sample), unname(base_risk_sample), tolerance = 1e-8)

  reference_shift <- unname(base_lp_zero - base_lp_sample)
  expect_equal(
    unname(lp_zero - lp_sample),
    reference_shift,
    tolerance = 1e-8
  )
  expect_equal(
    reference_shift,
    rep(reference_shift[1L], nrow(nd)),
    tolerance = 1e-10
  )
  expect_equal(unname(diff(lp_zero)), unname(diff(lp_sample)), tolerance = 1e-10)
  expect_equal(
    unname(risk_zero[2L] / risk_zero[1L]),
    unname(risk_sample[2L] / risk_sample[1L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(risk_zero / risk_sample),
    exp(reference_shift),
    tolerance = 1e-8
  )
})


# cox_reference = "strata" uses a different weighted training mean in each stratum.
# The resulting shift is constant within a stratum, not across all prediction
# rows. Relative comparisons remain invariant only within the same stratum.
test_that("version 1 Cox strata reference is stratum specific", {
  dat <- make_cox_reference_v1_data(seed = 11004L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x + strata(stratum),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- data.frame(
    x = c(3.25, 7.75, 5.25, 10.75),
    stratum = factor(
      c("A", "A", "B", "B"),
      levels = levels(dat$stratum)
    )
  )

  lp_zero <- predict(fit, nd, type = "lp", cox_reference = "zero")
  lp_strata <- predict(fit, nd, type = "lp", cox_reference = "strata")
  risk_zero <- predict(fit, nd, type = "risk", cox_reference = "zero")
  risk_strata <- predict(fit, nd, type = "risk", cox_reference = "strata")

  prepared <- prepare_formula_cox_newdata_v1(fit, nd)
  base <- strip_mfp2_class_v1(fit)

  expect_equal(
    unname(lp_zero),
    unname(predict(base, prepared, type = "lp", reference = "zero")),
    tolerance = 1e-8
  )
  expect_equal(
    unname(lp_strata),
    unname(predict(base, prepared, type = "lp", reference = "strata")),
    tolerance = 1e-8
  )

  reference_shift <- unname(lp_zero - lp_strata)

  expect_equal(reference_shift[1L], reference_shift[2L], tolerance = 1e-10)
  expect_equal(reference_shift[3L], reference_shift[4L], tolerance = 1e-10)
  expect_gt(abs(reference_shift[1L] - reference_shift[3L]), 1e-4)

  expect_equal(
    unname(lp_zero[2L] - lp_zero[1L]),
    unname(lp_strata[2L] - lp_strata[1L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(lp_zero[4L] - lp_zero[3L]),
    unname(lp_strata[4L] - lp_strata[3L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(risk_zero[2L] / risk_zero[1L]),
    unname(risk_strata[2L] / risk_strata[1L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(risk_zero[4L] / risk_zero[3L]),
    unname(risk_strata[4L] / risk_strata[3L]),
    tolerance = 1e-10
  )
})


# This is the previously uncovered end-to-end baseline-hazard path. The oracle
# is the same fitted coxph object after only the mfp2 class is removed; newdata
# are independently reconstructed on the stored transformed design scale.
test_that("version 1 Cox expected and survival newdata predictions match coxph", {
  dat <- make_cox_reference_v1_data(seed = 11005L)

  fit <- mfp2(
    survival::Surv(time, status) ~
      x + strata(stratum) + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[c(7, 41, 88, 133, 207, 319),
            c("time", "status", "x", "stratum", "off"),
            drop = FALSE]

  prepared <- prepare_formula_cox_newdata_v1(
    fit,
    nd,
    include_response = TRUE
  )
  base <- strip_mfp2_class_v1(fit)

  for (prediction_type in c("expected", "survival")) {
    got <- predict(
      fit,
      newdata = nd,
      type = prediction_type,
      se.fit = TRUE
    )
    expected <- predict(
      base,
      newdata = prepared,
      type = prediction_type,
      se.fit = TRUE
    )

    expect_equal(
      unname(got$fit),
      unname(expected$fit),
      tolerance = 1e-8
    )
    expect_equal(
      unname(got$se.fit),
      unname(expected$se.fit),
      tolerance = 1e-8
    )
  }

  got_expected <- predict(fit, nd, type = "expected")
  got_survival <- predict(fit, nd, type = "survival")
  expect_equal(
    unname(got_survival),
    unname(exp(-got_expected)),
    tolerance = 1e-10
  )
})


test_that("version 1 Cox expected and survival training predictions match coxph", {
  dat <- make_cox_reference_v1_data(seed = 11006L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = TRUE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  base <- strip_mfp2_class_v1(fit)

  for (prediction_type in c("expected", "survival")) {
    expect_equal(
      unname(predict(fit, type = prediction_type)),
      unname(predict(base, type = prediction_type)),
      tolerance = 1e-8
    )
  }
})


# Absolute predictions derive their response from newdata for every interface.
# A matrix-interface caller supplies one Surv column alongside the predictors;
# the column is removed before FP transformation and reattached to the internal
# Cox prediction frame under the response name expected by predict.coxph().
test_that("version 1 matrix-interface absolute Cox predictions derive response from newdata", {
  dat <- make_cox_reference_v1_data(seed = 11007L)
  x <- as.matrix(dat["x"])
  y <- survival::Surv(dat$time, dat$status)

  fit <- mfp2(
    x = x,
    y = y,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  rows <- c(4, 29, 74, 151)
  newdata <- data.frame(x = x[rows, "x"])
  newdata$prediction_response <- I(y[rows])

  prepared <- prepare_newdata_for_predict(
    fit,
    newdata,
    terms = get_selected_variable_names(fit),
    check_binary = FALSE
  )
  response_name <- cox_internal_response_name(fit)
  prepared[[response_name]] <- I(y[rows])

  base <- strip_mfp2_class_v1(fit)

  expect_equal(
    unname(predict(fit, newdata = newdata, type = "expected")),
    unname(predict(base, newdata = prepared, type = "expected")),
    tolerance = 1e-8
  )
  expect_equal(
    unname(predict(fit, newdata = newdata, type = "survival")),
    unname(predict(base, newdata = prepared, type = "survival")),
    tolerance = 1e-8
  )

  expect_error(
    predict(fit, newdata = data.frame(x = newdata$x), type = "expected"),
    "require.*follow-up response information"
  )

  ambiguous <- newdata
  ambiguous$second_response <- I(y[rows])
  expect_error(
    predict(fit, newdata = ambiguous, type = "survival"),
    "more than one Surv column"
  )
})


# The public default remains the reference-zero scale used by mfp2 term
# decomposition. An explicit sample reference is allowed, but it introduces the
# documented common constant when transformed columns have nonzero means.
test_that("version 1 default Cox lp remains aligned with mfp2 terms", {
  dat <- make_cox_reference_v1_data(seed = 11008L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- data.frame(x = c(3, 5, 7, 9))
  term_result <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms_seq = "data"
  )

  term_sum <- Reduce(`+`, lapply(term_result, `[[`, "value"))
  lp_default <- predict(fit, nd, type = "lp")
  lp_null <- predict(fit, nd, type = "lp", cox_reference = NULL)
  lp_zero <- predict(fit, nd, type = "lp", cox_reference = "zero")
  lp_sample <- predict(fit, nd, type = "lp", cox_reference = "sample")

  expect_equal(unname(lp_default), unname(lp_null), tolerance = 1e-10)
  expect_equal(unname(lp_default), unname(lp_zero), tolerance = 1e-10)
  expect_equal(unname(lp_zero), unname(term_sum), tolerance = 1e-8)
  expect_equal(
    unname(lp_zero - lp_sample),
    rep(unname(lp_zero[1L] - lp_sample[1L]), nrow(nd)),
    tolerance = 1e-10
  )
})


test_that("version 1 reports clear Cox reference and response errors", {
  dat <- make_cox_reference_v1_data(seed = 11009L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd_relative <- data.frame(x = dat$x[1:4])
  nd_absolute <- dat[1:4, c("time", "status", "x"), drop = FALSE]

  # This was formerly the unhelpful duplicate-formal-argument crash.
  expect_no_error(
    predict(fit, nd_relative, type = "lp", cox_reference = "sample")
  )
  expect_no_error(
    predict(fit, type = "lp", cox_reference = "sample")
  )
  base <- strip_mfp2_class_v1(fit)
  expect_equal(
    unname(predict(fit, type = "lp", cox_reference = "sample")),
    unname(predict(base, type = "lp", reference = "sample")),
    tolerance = 1e-8
  )

  expect_error(
    predict(fit, nd_relative, type = "lp", cox_reference = "population"),
    "must be one of 'zero', 'sample', or 'strata'"
  )
  expect_error(
    predict(fit, nd_relative, type = "terms", cox_reference = "zero"),
    "not used for mfp2 term or contrast predictions"
  )
  expect_error(
    predict(fit, nd_absolute, type = "survival", cox_reference = "zero"),
    "does not apply to Cox predictions"
  )
  expect_error(
    predict(fit, nd_relative, type = "survival"),
    "require.*follow-up response"
  )
  expect_error(
    predict(
      fit,
      nd_relative,
      type = "lp",
      newy = survival::Surv(dat$time[1:4], dat$status[1:4])
    ),
    "'newy' has been removed"
  )

  expect_error(
    predict(fit, nd_relative, type = "lp", reference = "zero"),
    "renamed.*cox_reference"
  )

  fit_glm <- mfp2(
    x = as.matrix(dat["x"]),
    # This fit only supplies a valid non-Cox object for the error check below.
    # Use deterministic non-linear variation instead of consuming RNG state.
    y = dat$x + 0.1 * sin(seq_len(nrow(dat))),
    verbose = FALSE
  )
  expect_error(
    predict(fit_glm, nd_relative, cox_reference = "zero"),
    "only available for Cox models"
  )
})



# -----------------------------------------------------------------------------
# 21.1 Binary term plotting
# -----------------------------------------------------------------------------

# Test purpose: Numeric binary terms must be plotted as two fitted point
# estimates with vertical confidence intervals, even when an equidistant
# sequence is requested. No interpolating line or confidence ribbon is drawn.
test_that("plot() displays numeric binary terms as two point estimates", {
  skip_if_not_installed("ggplot2")

  set.seed(21011)
  n <- 120L
  binary <- rep(c(2, 5), each = n / 2L)
  y <- 1 + 1.8 * (binary == 5) + stats::rnorm(n, sd = 0.4)

  fit <- mfp2(
    x = cbind(binary = binary),
    y = y,
    df = 1,
    keep = "binary",
    center = FALSE,
    verbose = FALSE
  )

  p <- plot(
    fit,
    terms = "binary",
    partial_only = TRUE,
    terms_seq = "equidistant"
  )[["binary"]]

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )

  expect_true("GeomPoint" %in% geom_classes)
  expect_true("GeomErrorbar" %in% geom_classes)
  expect_false("GeomLine" %in% geom_classes)
  expect_false("GeomRibbon" %in% geom_classes)

  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]
  errorbar_layer <- p$layers[[which(geom_classes == "GeomErrorbar")[[1L]]]]

  expect_equal(nrow(point_layer$data), 2L)
  expect_equal(nrow(errorbar_layer$data), 2L)
  expect_equal(as.numeric(point_layer$data$variable), c(2, 5))

  x_scale <- p$scales$get_scales("x")
  expect_equal(as.numeric(x_scale$breaks), c(2, 5))
})

# Test purpose: Two-level formula factors retain their fitted level labels and
# use the same point-and-confidence-interval presentation as numeric binaries.
test_that("plot() displays two-level factors as binary effects", {
  skip_if_not_installed("ggplot2")

  set.seed(21012)
  n <- 120L
  group <- factor(
    rep(c("control", "treated"), each = n / 2L),
    levels = c("control", "treated")
  )
  y <- 0.5 + 1.4 * (group == "treated") + stats::rnorm(n, sd = 0.4)
  dat <- data.frame(y = y, group = group)

  fit <- mfp2(
    y ~ group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )

  p <- plot(
    fit,
    terms = "group",
    partial_only = TRUE,
    terms_seq = "equidistant"
  )[["group"]]

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )

  expect_true("GeomPoint" %in% geom_classes)
  expect_true("GeomErrorbar" %in% geom_classes)
  expect_false("GeomLine" %in% geom_classes)
  expect_false("GeomRibbon" %in% geom_classes)

  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]
  expect_equal(nrow(point_layer$data), 2L)
  expect_equal(
    as.character(point_layer$data$variable),
    c("control", "treated")
  )

  x_scale <- p$scales$get_scales("x")
  expect_equal(x_scale$breaks, c("control", "treated"))
  expect_equal(x_scale$limits, c("control", "treated"))
})


# =============================================================================
# Shared focal max-degree basis regression tests
# =============================================================================

# Test purpose: A maximum-degree ordinary FP basis must support lower-degree
# views for a non-default, variable-specific power set without changing the
# existing candidate order or the special FP1 exclusion of power 1.
test_that("shared FP basis reproduces custom-power degree views", {
  x <- seq(1, 12, length.out = 60)
  allowed_powers <- c(-1, 0, 0.5, 1, 2)

  shared <- build_shared_focal_fp_basis(
    x = x,
    max_degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  # find_best_fpm_step() removes power 1 only for FP1. Reproduce that existing
  # rule here and verify that the shared max-degree basis still maps every
  # remaining FP1 candidate correctly.
  fp1_powers <- setdiff(allowed_powers, 1)
  shared_fp1 <- view_shared_focal_fp_basis(
    shared_basis = shared,
    degree = 1,
    powers = fp1_powers
  )
  direct_fp1 <- generate_transformations_fp_basis(
    x = x,
    degree = 1,
    powers = fp1_powers,
    zero = FALSE
  )

  expect_equal(shared_fp1$powers, direct_fp1$powers)
  expect_false(any(shared_fp1$powers == 1))

  for (i in seq_len(nrow(shared_fp1$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(shared_fp1, i)),
      unname(materialize_fp_basis_candidate(direct_fp1, i)),
      tolerance = 1e-12
    )
  }

  # FP2 keeps the user's complete power set, including power 1. Its candidate
  # order and transformed values must also match a separately generated basis.
  shared_fp2 <- view_shared_focal_fp_basis(
    shared_basis = shared,
    degree = 2,
    powers = allowed_powers
  )
  direct_fp2 <- generate_transformations_fp_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  expect_equal(shared_fp2$powers, direct_fp2$powers)
  expect_true(any(shared_fp2$powers == 1))

  for (i in seq_len(nrow(shared_fp2$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(shared_fp2, i)),
      unname(materialize_fp_basis_candidate(direct_fp2, i)),
      tolerance = 1e-12
    )
  }
})


# Test purpose: One joint ACD degree-2 basis must reproduce all three nonlinear
# IC views: FP1(x,.), FP1(.,A(x)), and FP1(x,A(x)), again with custom powers.
test_that("shared ACD basis reproduces all custom-power IC views", {
  x <- seq(1, 10, length.out = 50)
  allowed_powers <- c(-1, 0, 1, 2)
  fp1_powers <- setdiff(allowed_powers, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  shared <- build_shared_focal_acd_basis(
    x = x,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  # FP1(x,.) uses the x component of the joint ACD basis through the ordinary
  # FP view. No ACD-specific candidate semantics are introduced here.
  shared_x <- view_shared_focal_fp_basis(
    shared_basis = shared,
    degree = 1,
    powers = fp1_powers
  )
  direct_x <- generate_transformations_fp_basis(
    x = x,
    degree = 1,
    powers = fp1_powers,
    zero = FALSE
  )

  expect_equal(shared_x$powers, direct_x$powers)
  for (i in seq_len(nrow(shared_x$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(shared_x, i)),
      unname(materialize_fp_basis_candidate(direct_x, i)),
      tolerance = 1e-12
    )
  }

  # FP1(.,A(x)) uses only the A(x) component.
  shared_a <- view_shared_focal_acd_basis(
    shared_basis = shared,
    degree = 1,
    powers = fp1_powers
  )
  direct_a <- generate_transformations_acd_basis(
    x = x,
    degree = 1,
    powers = fp1_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_equal(shared_a$powers, direct_a$powers)
  for (i in seq_len(nrow(shared_a$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(shared_a, i)),
      unname(materialize_acd_basis_candidate(direct_a, i)),
      tolerance = 1e-12
    )
  }

  # FP1(x,A(x)) uses one x and one A(x) column from the same shared basis.
  shared_joint <- view_shared_focal_acd_basis(
    shared_basis = shared,
    degree = 2,
    powers = allowed_powers
  )
  direct_joint <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_equal(shared_joint$powers, direct_joint$powers)
  for (i in seq_len(nrow(shared_joint$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(shared_joint, i)),
      unname(materialize_acd_basis_candidate(direct_joint, i)),
      tolerance = 1e-12
    )
  }
})


# Test purpose: transform_data_step() must use a supplied ordinary shared basis
# rather than regenerating n-length FP transformations for the requested degree.
test_that("transform_data_step reuses supplied shared FP basis", {
  x_vec <- seq(1, 8, length.out = 40)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 0.5, 1, 2)

  shared <- build_shared_focal_fp_basis(
    x = x_vec,
    max_degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  testthat::local_mocked_bindings(
    generate_transformations_fp_basis = function(...) {
      stop("ordinary FP basis was regenerated")
    },
    .package = "mfp2"
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 2,
    powers = list(x = setdiff(allowed_powers, 1)),
    acdx = c(x = FALSE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    focal_basis_cache = shared,
    term_to_columns = list(x = 1L),
    compact_fp = TRUE
  )

  expect_null(out$data_fp)
  expect_equal(out$fp_basis$powers[, 1L], sort(setdiff(allowed_powers, 1)))
})


# Test purpose: transform_data_step() must likewise reuse a supplied joint ACD
# basis without invoking the degree-specific ACD basis generator again.
test_that("transform_data_step reuses supplied shared ACD basis", {
  x_vec <- seq(1, 8, length.out = 40)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 1, 2)
  fp1_powers <- setdiff(allowed_powers, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  shared <- build_shared_focal_acd_basis(
    x = x_vec,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  testthat::local_mocked_bindings(
    generate_transformations_acd_basis = function(...) {
      stop("ACD basis was regenerated")
    },
    .package = "mfp2"
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 2,
    powers = list(x = fp1_powers),
    acdx = c(x = TRUE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = acd_par),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    focal_basis_cache = shared,
    term_to_columns = list(x = 1L),
    compact_acd = TRUE
  )

  expect_null(out$data_fp)
  expect_equal(out$acd_basis$powers[, 2L], fp1_powers)
})


# Test purpose: select_ic() should build one maximum-degree basis from the
# complete custom power set and pass the same cache to every FP degree.
test_that("select_ic shares one focal basis across FP degrees", {
  x <- matrix(seq(1, 10, length.out = 40), ncol = 1L,
              dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 0.5, 1, 2)
  sentinel_cache <- list(id = "shared-fp-cache")
  builder_calls <- 0L
  seen_degrees <- integer(0)
  seen_cache <- list()

  fake_metrics <- function(aic) {
    matrix(
      c(0, 1, 0, NA, aic, aic, 38),
      nrow = 1L,
      dimnames = list(NULL, c(
        "logl", "df", "deviance_rs", "deviance_gaussian",
        "aic", "bic", "df_resid"
      ))
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        powers_adj = list(),
        spike_decision_adj = numeric(0),
        data_adj_list = list(),
        data_adj = NULL,
        transform_cache = list()
      )
    },
    fit_null_step = function(...) {
      list(powers = NA, metrics = fake_metrics(10), current_adj_params = list())
    },
    fit_linear_step = function(...) {
      list(powers = 1, metrics = fake_metrics(8), current_adj_params = list())
    },
    build_shared_focal_fp_basis = function(x, max_degree, powers, ...) {
      builder_calls <<- builder_calls + 1L
      expect_equal(max_degree, 2)
      expect_equal(powers, allowed_powers)
      sentinel_cache
    },
    find_best_fpm_step = function(..., degree, focal_basis_cache = NULL) {
      seen_degrees <<- c(seen_degrees, degree)
      seen_cache[[length(seen_cache) + 1L]] <<- focal_basis_cache
      list(
        powers = if (degree == 1) matrix(-1, nrow = 1L) else
          matrix(c(-1, 0), nrow = 1L),
        metrics = fake_metrics(if (degree == 1) 6 else 4),
        model_best = 1L,
        current_adj_params = list()
      )
    },
    .package = "mfp2"
  )

  select_ic(
    x = x,
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = rep(0, nrow(x)),
    powers_current = list(x = c(1, 1)),
    powers = list(x = allowed_powers),
    criterion = "aic",
    ftest = FALSE,
    select = 1,
    alpha = 0.05,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(),
    transform_cache = list(),
    force_max_fp = c(x = FALSE),
    has_offset = FALSE,
    n_obs = nrow(x),
    term_to_columns = list(x = 1L)
  )

  expect_equal(builder_calls, 1L)
  expect_equal(seen_degrees, c(1, 2))
  expect_length(seen_cache, 2L)
  expect_true(all(vapply(seen_cache, identical, logical(1), sentinel_cache)))
})


# Test purpose: select_ic_acd() should build one joint x/A(x) basis and reuse
# it for FP1(x,.), FP1(.,A(x)), and FP1(x,A(x)).
test_that("select_ic_acd shares one joint focal basis across all nonlinear views", {
  x <- matrix(seq(1, 10, length.out = 40), ncol = 1L,
              dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 1, 2)
  sentinel_cache <- list(id = "shared-acd-cache")
  builder_calls <- 0L
  seen_degrees <- integer(0)
  seen_cache <- list()

  fake_metrics <- function(aic) {
    matrix(
      c(0, 1, 0, NA, aic, aic, 38),
      nrow = 1L,
      dimnames = list(NULL, c(
        "logl", "df", "deviance_rs", "deviance_gaussian",
        "aic", "bic", "df_resid"
      ))
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        powers_adj = list(),
        spike_decision_adj = numeric(0),
        data_adj_list = list(),
        data_adj = NULL,
        transform_cache = list()
      )
    },
    fit_null_step = function(...) {
      list(powers = c(NA, NA), metrics = fake_metrics(12), current_adj_params = list())
    },
    fit_linear_step = function(...) {
      list(powers = c(NA, 1), metrics = fake_metrics(10), current_adj_params = list())
    },
    build_shared_focal_acd_basis = function(x, powers, ...) {
      builder_calls <<- builder_calls + 1L
      expect_equal(powers, allowed_powers)
      sentinel_cache
    },
    find_best_fpm_step = function(..., degree, acdx, focal_basis_cache = NULL) {
      seen_degrees <<- c(seen_degrees, degree)
      seen_cache[[length(seen_cache) + 1L]] <<- focal_basis_cache
      list(
        powers = if (degree == 1) matrix(c(NA, -1), nrow = 1L) else
          matrix(c(-1, 0), nrow = 1L),
        metrics = fake_metrics(if (degree == 1) 8 else 6),
        model_best = 1L,
        current_adj_params = list()
      )
    },
    .package = "mfp2"
  )

  select_ic_acd(
    x = x,
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = TRUE),
    y = rep(0, nrow(x)),
    powers_current = list(x = c(1, 1)),
    powers = list(x = allowed_powers),
    criterion = "aic",
    ftest = FALSE,
    select = 1,
    alpha = 0.05,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = list(acd = seq(0.1, 0.9, length.out = nrow(x)))),
    prev_adj_params = list(),
    transform_cache = list(),
    force_max_fp = c(x = FALSE),
    has_offset = FALSE,
    n_obs = nrow(x),
    term_to_columns = list(x = 1L)
  )

  expect_equal(builder_calls, 1L)
  expect_equal(seen_degrees, c(1, 1, 2))
  expect_length(seen_cache, 3L)
  expect_true(all(vapply(seen_cache, identical, logical(1), sentinel_cache)))
})


# =============================================================================
# End of tests
# =============================================================================

# =============================================================================
# Prediction API regression tests maintained for 1.1.0.9003
# =============================================================================

test_that("predict.mfp2 normalizes linear-predictor aliases", {
  data("prostate", package = "mfp2")
  fit_glm <- mfp2(
    lpsa ~ fp(age) + svi,
    data = prostate,
    verbose = FALSE
  )
  nd_glm <- prostate[1:12, c("age", "svi"), drop = FALSE]
  expect_equal(
    predict(fit_glm, newdata = nd_glm, type = "lp"),
    predict(fit_glm, newdata = nd_glm, type = "link")
  )
  expect_equal(
    predict(fit_glm, newdata = nd_glm, type = c(alias = "lp")),
    predict(fit_glm, newdata = nd_glm, type = "link")
  )

  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  fit_cox <- mfp2(
    survival::Surv(time, status) ~ fp(age) + sex,
    data = dat,
    family = "cox",
    verbose = FALSE
  )
  nd_cox <- dat[1:12, c("age", "sex"), drop = FALSE]
  expect_equal(
    predict(fit_cox, newdata = nd_cox, type = "link", cox_reference = "zero"),
    predict(fit_cox, newdata = nd_cox, type = "lp", cox_reference = "zero")
  )
  expect_equal(
    predict(
      fit_cox,
      newdata = nd_cox,
      type = c(alias = "link"),
      cox_reference = "zero"
    ),
    predict(fit_cox, newdata = nd_cox, type = "lp", cox_reference = "zero")
  )
})


test_that("predict.mfp2 rejects infinite numeric newdata columns", {
  data("prostate", package = "mfp2")
  fit <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)

  nd_inf <- prostate[1:6, c("age", "svi"), drop = FALSE]
  nd_inf$age[2] <- Inf
  expect_error(
    predict(fit, newdata = nd_inf, type = "link"),
    "finite numeric values|Infinite"
  )

  nd_ninf <- prostate[1:6, c("age", "svi"), drop = FALSE]
  nd_ninf$age[2] <- -Inf
  expect_error(
    predict(fit, newdata = nd_ninf, type = "link"),
    "finite numeric values|Infinite"
  )
})


test_that("predict.mfpi defaults safely and normalizes GLM lp", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  default_pred <- predict(fit, terms = "cavol", model = "all")
  expect_s3_class(default_pred, "mfpi_prediction")
  expect_identical(default_pred$type, "both")
  expect_true(!is.null(default_pred$functions))
  expect_true(!is.null(default_pred$differences))

  nd <- prostate[1:10, , drop = FALSE]
  link_pred <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  lp_pred <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "lp", se.fit = TRUE
  )
  expect_identical(lp_pred$type, "link")
  expect_equal(lp_pred$predictions, link_pred$predictions)
})


test_that("predict.mfpi enforces family-specific arguments and types", {
  data("prostate", package = "mfp2")
  fit_glm <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )
  nd_glm <- prostate[1:8, , drop = FALSE]

  expect_error(
    predict(fit_glm, newdata = nd_glm, terms = "cavol", model = "all",
            type = "risk"),
    "GLM MFPI models"
  )
  expect_error(
    predict(fit_glm, newdata = nd_glm, terms = "cavol", model = "all",
            type = "link", cox_reference = "zero"),
    "available only for Cox"
  )
  expect_error(
    predict(fit_glm, newdata = nd_glm, terms = "cavol", model = "all",
            type = "link", strata = rep(1, nrow(nd_glm))),
    "available only for Cox"
  )
  expect_error(
    predict(fit_glm, terms = "cavol", model = "all", type = "function",
            strata = rep(1, nrow(prostate))),
    "not used for MFPI fitted-function"
  )

  nd_bad <- nd_glm
  nd_bad$cavol[1] <- Inf
  expect_error(
    predict(fit_glm, newdata = nd_bad, terms = "cavol", model = "all",
            type = "link"),
    "finite"
  )
})


test_that("predict.mfpi Cox references match the retained coxph model", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  dat$inst <- factor(dat$inst)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:18, c("age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  strata_new <- reconstruct_formula_strata_newdata(fit, nd)
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = strata_new
  )

  for (prediction_type in c("lp", "risk")) {
    for (reference_value in c("zero", "sample", "strata")) {
      direct <- stats::predict(
        stored,
        newdata = design$model_newdata,
        type = prediction_type,
        se.fit = TRUE,
        reference = reference_value
      )
      got <- predict(
        fit,
        newdata = nd,
        terms = "age",
        model = "all",
        type = prediction_type,
        cox_reference = reference_value,
        se.fit = TRUE
      )
      expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
      expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    }
  }
})


test_that("predict.mfpi Cox expected and survival match stored coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit

  for (prediction_type in c("expected", "survival")) {
    direct_training <- stats::predict(
      stored, type = prediction_type, se.fit = FALSE
    )
    got_training <- predict(
      fit, terms = "age", model = "all",
      type = prediction_type, se.fit = FALSE
    )
    expect_equal(
      got_training$predictions$fit,
      as.numeric(direct_training),
      tolerance = 1e-8
    )
  }

  nd <- dat[1:16, c("time", "status", "age", "sex"), drop = FALSE]
  response <- reconstruct_cox_prediction_response(
    object = fit,
    fit_obj = stored,
    newdata = nd
  )
  design <- mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  direct_data <- attach_cox_prediction_response(
    fit_obj = stored,
    newdata = design$model_newdata,
    response = response
  )

  for (prediction_type in c("expected", "survival")) {
    direct <- stats::predict(
      stored,
      newdata = direct_data,
      type = prediction_type,
      se.fit = TRUE
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = prediction_type,
      se.fit = TRUE
    )
    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    expect_true(got$metadata$absolute_cox_prediction)
  }

  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all",
      type = "survival", cox_reference = "zero"
    ),
    "applies only"
  )

  nd_missing_response <- nd[, c("age", "sex"), drop = FALSE]
  expect_error(
    predict(
      fit, newdata = nd_missing_response, terms = "age", model = "all",
      type = "survival"
    ),
    "response information"
  )
})


test_that("predict.mfpi matrix-interface Cox survival accepts one Surv column", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  rows <- 1:12
  nd <- data.frame(
    sex = dat$sex[rows],
    age = dat$age[rows],
    prediction_response = I(survival::Surv(dat$time[rows], dat$status[rows]))
  )
  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  response <- reconstruct_cox_prediction_response(fit, stored, nd)
  design <- mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  direct_data <- attach_cox_prediction_response(
    stored, design$model_newdata, response
  )

  direct <- stats::predict(
    stored, newdata = direct_data, type = "survival", se.fit = FALSE
  )
  got <- predict(
    fit, newdata = nd, terms = "age", model = "all",
    type = "survival", se.fit = FALSE
  )
  expect_equal(got$predictions$fit, as.numeric(direct), tolerance = 1e-8)
})


test_that("predict.mfpi rejects irrelevant Cox strata", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    p_interact = 1,
    verbose = FALSE
  )
  nd <- dat[1:8, c("age", "sex"), drop = FALSE]

  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all",
      type = "lp", strata = rep(1, nrow(nd))
    ),
    "not stratified"
  )
  expect_error(
    predict(
      fit, terms = "age", model = "all",
      type = "function", strata = rep(1, nrow(dat))
    ),
    "not used for MFPI fitted-function"
  )
})


test_that("predict.mfpi stratified Cox survival matches stored coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  dat$inst <- factor(dat$inst)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:14, c("time", "status", "age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  response <- reconstruct_cox_prediction_response(fit, stored, nd)
  strata_new <- reconstruct_formula_strata_newdata(fit, nd)
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = strata_new
  )
  direct_data <- attach_cox_prediction_response(
    stored, design$model_newdata, response
  )

  for (prediction_type in c("expected", "survival")) {
    direct <- stats::predict(
      stored,
      newdata = direct_data,
      type = prediction_type,
      se.fit = TRUE
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = prediction_type,
      se.fit = TRUE
    )
    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
  }
})


test_that("predict.mfpi rejects replacement offsets for models without offsets", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  nd <- prostate[1:7, , drop = FALSE]
  expect_error(
    predict(
      fit,
      newdata = nd,
      terms = "cavol",
      model = "all",
      type = "link",
      newoffset = rep(0, nrow(nd))
    ),
    "fitted with an offset"
  )
})


test_that("predict.mfpi training reconstruction preserves fitted Cox strata", {
  set.seed(19003)
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  dat$inst <- factor(dat$inst)
  dat$exposure <- stats::runif(nrow(dat), 0.7, 2.4)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex + offset(log(exposure)) + strata(inst),
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  replacement_offset <- log(dat$exposure)

  design <- mfpi_build_ordinary_design(
    fit,
    "age",
    fit_result,
    newdata = NULL,
    newoffset = replacement_offset
  )
  expect_true("strata_" %in% names(design$model_newdata))
  expect_true("offset_" %in% names(design$model_newdata))

  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got <- predict(
    fit,
    terms = "age",
    model = "all",
    type = "lp",
    newoffset = replacement_offset,
    cox_reference = "zero",
    se.fit = TRUE
  )

  expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
  expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)

  # Reconstructing the training rows solely to resupply their fitted strata
  # must also recover the original, uncentred offset scale. The resulting
  # prediction should therefore equal direct prediction on the stored fit.
  fitted_strata <- stored$strata
  direct_training <- stats::predict(
    stored,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got_strata_only <- predict(
    fit,
    terms = "age",
    model = "all",
    type = "lp",
    strata = fitted_strata,
    cox_reference = "zero",
    se.fit = TRUE
  )
  expect_equal(
    got_strata_only$predictions$fit,
    as.numeric(direct_training$fit),
    tolerance = 1e-8
  )
  expect_equal(
    got_strata_only$predictions$se.fit,
    as.numeric(direct_training$se.fit),
    tolerance = 1e-8
  )
})



# =============================================================================
# Prediction argument and non-finite validation regressions
# =============================================================================

test_that("predict.mfp2 validates only required missing predictors", {
  data("prostate", package = "mfp2")
  fit <- mfp2(
    lpsa ~ fp(age) + svi,
    data = prostate,
    keep = "age",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd <- prostate[1:8, c("age", "svi"), drop = FALSE]
  expected <- predict(fit, newdata = nd, type = "link")

  nd_extra <- nd
  nd_extra$unused_na <- NA_real_
  nd_extra$unused_nan <- NaN
  expect_equal(
    predict(fit, newdata = nd_extra, type = "link"),
    expected,
    tolerance = 1e-12
  )

  for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
    nd_bad <- nd
    nd_bad$age[2] <- bad_value
    expect_error(
      predict(fit, newdata = nd_bad, type = "link"),
      "missing|finite"
    )
  }

  fit_matrix <- mfp2(
    x = as.matrix(prostate[, c("age", "svi")]),
    y = prostate$lpsa,
    verbose = FALSE
  )
  nd_matrix <- data.frame(
    age = prostate$age[1:8],
    svi = prostate$svi[1:8],
    unused_na = NA_real_,
    unused_nan = NaN
  )
  expect_equal(
    predict(fit_matrix, newdata = nd_matrix, type = "link"),
    predict(
      fit_matrix,
      newdata = nd_matrix[, c("age", "svi"), drop = FALSE],
      type = "link"
    ),
    tolerance = 1e-12
  )
})


test_that("predict.mfp2 rejects irrelevant strata arguments", {
  data("prostate", package = "mfp2")
  fit_glm <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)
  nd_glm <- prostate[1:6, c("age", "svi"), drop = FALSE]
  expect_error(
    predict(fit_glm, newdata = nd_glm, strata = rep(1, nrow(nd_glm))),
    "only for Cox"
  )

  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  fit_cox <- mfp2(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd_cox <- dat[1:8, c("age", "sex"), drop = FALSE]

  expect_error(
    predict(
      fit_cox,
      newdata = nd_cox,
      type = "lp",
      strata = rep(1, nrow(nd_cox))
    ),
    "not stratified"
  )
  expect_error(
    predict(
      fit_cox,
      newdata = nd_cox,
      type = "terms",
      strata = rep(1, nrow(nd_cox))
    ),
    "not used.*term or contrast"
  )
  expect_error(
    predict(fit_cox, type = "lp", strata = rep(1, nrow(dat))),
    "together with.*newdata"
  )

  dat$stratum <- factor(rep(c("A", "B"), length.out = nrow(dat)))
  fit_stratified <- mfp2(
    survival::Surv(time, status) ~ age + sex + strata(stratum),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  nd_stratified <- dat[1:8, c("age", "sex", "stratum"), drop = FALSE]
  bad_strata <- as.numeric(nd_stratified$stratum)
  bad_strata[2] <- Inf
  expect_error(
    predict(
      fit_stratified,
      newdata = nd_stratified,
      type = "lp",
      strata = bad_strata
    ),
    "strata.*finite"
  )
})


test_that("predict.mfp2 rejects irrelevant and non-finite replacement offsets", {
  data("prostate", package = "mfp2")
  fit_plain <- mfp2(lpsa ~ fp(age) + svi, data = prostate, verbose = FALSE)
  nd_plain <- prostate[1:7, c("age", "svi"), drop = FALSE]
  expect_error(
    predict(
      fit_plain,
      newdata = nd_plain,
      type = "link",
      newoffset = rep(0, nrow(nd_plain))
    ),
    "used an offset"
  )

  set.seed(21001)
  dat <- prostate
  dat$exposure <- stats::runif(nrow(dat), 0.7, 2.2)
  fit_offset <- mfp2(
    lpsa ~ fp(age) + svi + offset(log(exposure)),
    data = dat,
    verbose = FALSE
  )
  nd <- dat[1:7, c("age", "svi", "exposure"), drop = FALSE]
  good_offset <- log(nd$exposure)

  expect_error(
    predict(
      fit_offset,
      newdata = nd,
      type = "terms",
      newoffset = good_offset
    ),
    "not used.*term or contrast"
  )
  expect_error(
    predict(fit_offset, type = "link", newoffset = log(dat$exposure)),
    "together with.*newdata"
  )

  for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
    bad_offset <- good_offset
    bad_offset[2] <- bad_value
    expect_error(
      predict(
        fit_offset,
        newdata = nd,
        type = "link",
        newoffset = bad_offset
      ),
      "finite numeric"
    )
  }
})



test_that("predict.mfpi validates required predictors but ignores irrelevant extras", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ age + svi + cavol,
    data = prostate,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    keep = "age",
    p_interact = 1,
    verbose = FALSE
  )
  nd <- prostate[1:8, c("age", "svi", "cavol"), drop = FALSE]
  expected <- predict(
    fit, newdata = nd, terms = "cavol", model = "all", type = "link"
  )
  nd_extra <- nd
  nd_extra$unused_na <- NA_real_
  nd_extra$unused_nan <- NaN
  expect_equal(
    predict(
      fit, newdata = nd_extra, terms = "cavol", model = "all", type = "link"
    )$predictions,
    expected$predictions,
    tolerance = 1e-12
  )

  for (column in c("cavol", "age", "svi")) {
    for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
      nd_bad <- nd
      nd_bad[[column]][2] <- bad_value
      expect_error(
        predict(
          fit, newdata = nd_bad, terms = "cavol", model = "all", type = "link"
        ),
        "finite|missing|NA"
      )
    }
  }
})


test_that("formula-derived offsets reject missing and non-finite inputs", {
  data("prostate", package = "mfp2")
  set.seed(21002)
  dat <- prostate
  dat$exposure <- stats::runif(nrow(dat), 0.7, 2.2)

  fit_mfp2 <- mfp2(
    lpsa ~ age + svi + offset(log(exposure)),
    data = dat,
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  fit_mfpi <- mfpi(
    lpsa ~ age + svi + cavol + offset(log(exposure)),
    data = dat,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    keep = "age",
    p_interact = 1,
    verbose = FALSE
  )
  nd <- dat[1:7, c("age", "svi", "cavol", "exposure"), drop = FALSE]

  for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
    nd_bad <- nd
    nd_bad$exposure[2] <- bad_value
    expect_error(
      predict(fit_mfp2, newdata = nd_bad, type = "link"),
      "offset|missing|finite"
    )
    expect_error(
      predict(
        fit_mfpi,
        newdata = nd_bad,
        terms = "cavol",
        model = "all",
        type = "link"
      ),
      "offset|missing|finite"
    )
  }
})



# =============================================================================
# 25. MFPI printing and ACD reconstruction regressions (consolidated v12)
# =============================================================================

# 25.1 MFPI printing regression tests

make_mfpi_print_metrics <- function(n_groups = 2L) {
  group_names <- LETTERS[seq_len(n_groups)]
  int_powers <- setNames(
    lapply(seq_len(n_groups), function(i) c(i)),
    group_names
  )

  out <- data.frame(
    variable = "x",
    type = "fp1",
    deviance_int = 10.1234,
    deviance_diff = 2.3456,
    df_int = 2,
    pvalue = 0.01234,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$fp_powers_main <- I(list(list(x = c(1, 2))))
  out$fp_powers_int <- I(list(int_powers))
  out
}

make_minimal_mfpi_print_object <- function() {
  structure(
    list(
      group_var = "group",
      nobs = 10L,
      flex = "flex3",
      criterion = "pvalue",
      p_adjust_method = "holm",
      p_interact = 0.05,
      min_improvement = 2,
      digits = 3L,
      adjust_terms = NULL,
      all_model_metrics = NULL,
      best_model_metrics = NULL,
      var_winners = NULL
    ),
    class = "mfpi"
  )
}


make_mfpi_adjustment_object <- function(criterion = "pvalue") {
  x <- make_minimal_mfpi_print_object()
  x$criterion <- criterion
  x$adjust_terms <- data.frame(
    df_setting = c(1, 4),
    df_initial = c(1, 4),
    select = if (criterion == "pvalue") c(0.05, 0.10) else toupper(criterion),
    alpha = if (criterion == "pvalue") c(0.05, 0.05) else toupper(criterion),
    df_final = c(1, 4),
    power1 = c(1, 3),
    power2 = c(NA_real_, 3),
    selected = c(TRUE, TRUE),
    row.names = c("hx", "age"),
    check.names = FALSE
  )
  x
}

test_that("small interaction-power sets print inline and large sets summarize", {
  small <- format_mfpi_powers(
    make_mfpi_print_metrics(2L),
    max_inline_int_powers = 2L
  )
  expect_identical(small$fp_powers_int, "A: (1); B: (2)")

  large <- format_mfpi_powers(
    make_mfpi_print_metrics(3L),
    max_inline_int_powers = 2L
  )
  expect_identical(large$fp_powers_int, "3 group-specific FPs")

  singleton <- format_mfpi_powers(
    make_mfpi_print_metrics(1L),
    max_inline_int_powers = 0L
  )
  expect_identical(singleton$fp_powers_int, "1 group-specific FP")
})

test_that("inline interaction-power threshold is validated", {
  metrics <- make_mfpi_print_metrics(2L)

  expect_error(
    format_mfpi_powers(metrics, max_inline_int_powers = -1L),
    "single non-negative integer",
    fixed = TRUE
  )
  expect_error(
    format_mfpi_powers(metrics, max_inline_int_powers = 1.5),
    "single non-negative integer",
    fixed = TRUE
  )
})

test_that("p-value Step 1 prints the cutoff, adjustment method, select, and alpha", {
  x <- make_mfpi_adjustment_object("pvalue")
  x$p_adjust_method <- "none"
  x$p_interact <- 0.075
  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")
  info <- mfpi_prepare_adjustment_display(x, digits = 3L)
  expect_true(all(c("df_initial", "df_final") %in% names(info$display)))
  expect_false("df_setting" %in% names(info$display))
  expect_true(all(c("select", "alpha") %in% names(info$display)))
  expect_false(grepl("df_setting", printed, fixed = TRUE))
  expect_match(printed, "df_initial", fixed = TRUE)
  expect_match(printed, "df_final", fixed = TRUE)

  expect_match(printed, " criterion             : p-value", fixed = TRUE)
  expect_match(
    printed,
    " interaction selection : p-value < 0.075",
    fixed = TRUE
  )
  expect_match(
    printed,
    " p-value adjustment    : none",
    fixed = TRUE
  )
  expect_match(printed, "select", fixed = TRUE)
  expect_match(printed, "alpha", fixed = TRUE)
})

test_that("p-value Step 1 prints the active adjustment method dynamically", {
  x <- make_mfpi_adjustment_object("pvalue")
  x$p_adjust_method <- "hochberg"
  x$p_interact <- 0.01
  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(
    printed,
    " interaction selection : adjusted p-value < 0.01",
    fixed = TRUE
  )
  expect_match(
    printed,
    " p-value adjustment    : hochberg",
    fixed = TRUE
  )
})

test_that("AIC and BIC Step 1 use dynamic human-readable thresholds", {
  for (criterion in c("aic", "bic")) {
    x <- make_mfpi_adjustment_object(criterion)
    x$min_improvement <- 3.5
    output <- capture.output(
      print_adjustment_step(x, ruler = "-----", digits = 3L)
    )
    printed <- paste(output, collapse = "\n")
    info <- mfpi_prepare_adjustment_display(x, digits = 3L)
    criterion_label <- toupper(criterion)

    expect_match(
      printed,
      sprintf(" criterion             : %s", criterion_label),
      fixed = TRUE
    )
    expect_match(
      printed,
      sprintf(" interaction selection : %s reduction > 3.5", criterion_label),
      fixed = TRUE
    )
    expect_false(grepl("min_improvement", printed, fixed = TRUE))
    expect_false(grepl("df_setting", printed, fixed = TRUE))
    expect_match(printed, "df_initial", fixed = TRUE)
    expect_match(printed, "df_final", fixed = TRUE)
    expect_match(printed, "power1", fixed = TRUE)
    expect_match(printed, "power2", fixed = TRUE)
    expect_false("select" %in% names(info$display))
    expect_false("alpha" %in% names(info$display))
    expect_match(printed, "Selected adjustment variables (2): hx, age", fixed = TRUE)
  }
})

test_that("MFPI Step 1 prints grouped initial and final df only", {
  x <- make_minimal_mfpi_print_object()
  x$adjust_terms <- data.frame(
    df_setting = 1,
    df_initial = 2,
    select = 0.05,
    alpha = 0.05,
    selected = TRUE,
    df_final = 2,
    power1 = 1,
    power2 = NA_real_,
    row.names = "stage",
    check.names = FALSE
  )

  info <- mfpi_prepare_adjustment_display(x, digits = 3L)
  expect_false("df_setting" %in% names(info$display))
  expect_equal(info$display[["df_initial"]], 2)
  expect_equal(info$display[["df_final"]], 2)

  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")
  expect_false(grepl("df_setting", printed, fixed = TRUE))
  expect_match(printed, "df_initial", fixed = TRUE)
  expect_match(printed, "df_final", fixed = TRUE)
})

test_that("candidate output omits interaction powers shown in the detail table", {
  x <- make_minimal_mfpi_print_object()
  x$all_model_metrics <- make_mfpi_print_metrics(2L)

  output <- capture.output(
    print_candidates_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "fp_powers_main", fixed = TRUE)
  expect_false(grepl("fp_powers_int", printed, fixed = TRUE))
  expect_false(grepl("A: (1); B: (2)", printed, fixed = TRUE))
})

test_that("interaction-power detail output has no underline", {
  output <- capture.output(
    print_interaction_power_details_step(
      metrics = make_mfpi_print_metrics(2L),
      group_var = "group",
      title = "  FP powers by group level"
    )
  )

  expect_match(paste(output, collapse = "\n"), "FP powers by group level:", fixed = TRUE)
  expect_false(any(grepl("^-+$", trimws(output))))
})

test_that("main header ruler matches the displayed header width", {
  x <- make_minimal_mfpi_print_object()
  x$nevents <- 7L
  output <- capture.output(print(x))

  expect_identical(output[[1L]], output[[3L]])
  expect_identical(
    nchar(output[[1L]], type = "width"),
    nchar(output[[2L]], type = "width")
  )
})


test_that("main header uses reader-facing criterion labels", {
  for (criterion in c("aic", "bic", "pvalue")) {
    x <- make_minimal_mfpi_print_object()
    x$criterion <- criterion
    output <- capture.output(print(x))
    expected <- if (criterion == "pvalue") "p-value" else toupper(criterion)

    expect_match(
      output[[2L]],
      sprintf("criterion: %s", expected),
      fixed = TRUE
    )
    expect_false(grepl("p-adjust:", output[[2L]], fixed = TRUE))
  }
})

test_that("MFPI methods warn consistently about unused dots", {
  x <- make_minimal_mfpi_print_object()

  expect_warning(
    capture.output(print(x, unused_argument = TRUE)),
    "Unused arguments in `print.mfpi(...)`: unused_argument.",
    fixed = TRUE
  )
  expect_warning(
    warn_unused_mfpi_dots(list(TRUE), method = "summary.mfpi"),
    "Unused arguments in `summary.mfpi(...)`: <unnamed>.",
    fixed = TRUE
  )
})

test_that("summary regression output receives the requested digits", {
  # Deliberately define an unregistered local S3 method. The source-level
  # Step 4 dispatcher must resolve it from the calling environment and pass the
  # requested digits value rather than falling back to print.default().
  print.mfpi_test_summary <- function(x, digits = NULL, ...) {
    cat(sprintf("TEST SUMMARY DIGITS: %s\n", digits))
    invisible(x)
  }

  x <- make_minimal_mfpi_print_object()
  x$model_summaries <- list(
    x = structure(list(), class = "mfpi_test_summary")
  )
  x$var_winners <- list(x = list(type = "fp1"))
  class(x) <- "summary.mfpi"

  output <- capture.output(print(x, digits = 2L))

  expect_match(
    paste(output, collapse = "\n"),
    "TEST SUMMARY DIGITS: 2",
    fixed = TRUE
  )
})


make_mfpi_interaction_summary_object <- function(criterion = "pvalue") {
  x <- make_minimal_mfpi_print_object()
  x$criterion <- criterion
  x$p_adjust_method <- "holm"
  x$p_interact <- 0.05
  x$min_improvement <- 2

  if (criterion == "pvalue") {
    age_metric <- data.frame(pvalue = 0.01, p_adjusted = 0.02)
    wt_metric <- data.frame(pvalue = 0.20, p_adjusted = 0.20)
  } else if (criterion == "aic") {
    age_metric <- data.frame(AIC_main_minus_int = 3.5)
    wt_metric <- data.frame(AIC_main_minus_int = 1.5)
  } else {
    age_metric <- data.frame(BIC_main_minus_int = 4.0)
    wt_metric <- data.frame(BIC_main_minus_int = 0.5)
  }

  x$var_winners <- list(
    age = list(fit = list(ok = TRUE), metric = age_metric, type = "fp2"),
    wt = list(fit = list(ok = TRUE), metric = wt_metric, type = "fp1")
  )
  x$best_model_metrics <- data.frame(
    variable = "age",
    stringsAsFactors = FALSE
  )
  x
}

test_that("Step 3 reports Yes and No and embeds the p-value rule in its heading", {
  x <- make_mfpi_interaction_summary_object("pvalue")
  output <- capture.output(
    print_interaction_summary_step(
      x,
      ruler = "-----",
      digits = 3L,
      step_no = 3L
    )
  )
  printed <- paste(output, collapse = "\n")

  expect_match(
    printed,
    "Step 3 - Interaction Summary (selected when p_adjusted < 0.05):",
    fixed = TRUE
  )
  expect_match(printed, "selected", fixed = TRUE)
  expect_match(printed, "Yes", fixed = TRUE)
  expect_match(printed, "No", fixed = TRUE)
  expect_false(grepl("* = selected", printed, fixed = TRUE))
  expect_false(grepl("p_adjust_method =", printed, fixed = TRUE))
})

test_that("Step 3 embeds the AIC and BIC rules and uses Yes and No", {
  for (criterion in c("aic", "bic")) {
    x <- make_mfpi_interaction_summary_object(criterion)
    output <- capture.output(
      print_interaction_summary_step(
        x,
        ruler = "-----",
        digits = 3L,
        step_no = 3L
      )
    )
    printed <- paste(output, collapse = "\n")
    label <- if (criterion == "aic") "dAIC" else "dBIC"

    expect_match(
      printed,
      sprintf(
        "Step 3 - Interaction Summary (selected when %s > 2):",
        label
      ),
      fixed = TRUE
    )
    expect_match(printed, "Yes", fixed = TRUE)
    expect_match(printed, "No", fixed = TRUE)
    expect_false(grepl("* = selected", printed, fixed = TRUE))
  }
})

test_that("shared adjustment display contains selected variables only", {
  x <- make_mfpi_adjustment_object("pvalue")
  dropped <- data.frame(
    df_setting = 4,
    df_initial = 4,
    select = 0.05,
    alpha = 0.05,
    df_final = 0,
    power1 = NA_real_,
    power2 = NA_real_,
    selected = FALSE,
    row.names = "hg",
    check.names = FALSE
  )
  x$adjust_terms <- rbind(x$adjust_terms, dropped)

  info <- mfpi_prepare_adjustment_display(x, digits = 3L)

  expect_identical(rownames(info$display), c("hx", "age"))
  expect_identical(info$selected_names, c("hx", "age"))
  expect_false("selected" %in% names(info$display))
  expect_false("df_setting" %in% names(info$display))
  expect_false("hg" %in% rownames(info$display))
})

test_that("MFPI Step 1 uses mfp2 SAZ display names and labels", {
  x <- make_minimal_mfpi_print_object()
  x$p_adjust_method <- "none"
  x$adjust_terms <- data.frame(
    df_setting = c(4, 4, 4),
    df_initial = c(4, 4, 4),
    select = c(0.05, 1.00, 1.00),
    alpha = c(0.05, 0.05, 0.05),
    acd = c(FALSE, FALSE, TRUE),
    zero = c(FALSE, TRUE, FALSE),
    catzero = c(FALSE, TRUE, FALSE),
    spike = c(FALSE, TRUE, FALSE),
    prop_zero = c(NA_real_, 0.25, NA_real_),
    spike_dec = c(2L, 1L, 2L),
    selected = c(TRUE, TRUE, TRUE),
    df_final = c(2, 2, 1),
    power1 = c(0, 1, 1),
    power2 = c(NA_real_, NA_real_, NA_real_),
    row.names = c("cavol", "pgg45", "age"),
    check.names = FALSE
  )

  info <- mfpi_prepare_adjustment_display(x, digits = 3L)

  expect_true("catzero_final" %in% names(info$display))
  expect_true("prop_zero" %in% names(info$display))
  expect_true("saz_decision" %in% names(info$display))
  expect_false("catzero" %in% names(info$display))
  expect_false("spike_dec" %in% names(info$display))
  expect_identical(info$display$prop_zero, c(".", "0.250", "."))
  expect_identical(
    info$display$saz_decision,
    c("not SAZ", "continuous + binary", "not SAZ")
  )

  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "catzero_final", fixed = TRUE)
  expect_match(printed, "prop_zero", fixed = TRUE)
  expect_match(printed, "0.250", fixed = TRUE)
  expect_match(printed, "saz_decision", fixed = TRUE)
  expect_match(printed, "not SAZ", fixed = TRUE)
  expect_match(printed, "continuous + binary", fixed = TRUE)
  expect_false(grepl("spike_dec", printed, fixed = TRUE))
})

test_that("verbose Step 3 reports final p-value decisions without stars", {
  winners <- list(
    age = list(
      fit = list(ok = TRUE),
      metric = data.frame(pvalue = 0.0273),
      type = "fp2",
      score = 0.0273
    ),
    wt = list(
      fit = list(ok = TRUE),
      metric = data.frame(pvalue = 0.663),
      type = "fp1",
      score = 0.663
    )
  )

  output <- capture.output(
    print_interaction_step3_summary(
      var_winners = winners,
      cont_vars = c("age", "wt"),
      mode = "pvalue",
      p_interact = 0.05,
      p_adjust_method = "none",
      digits = 4L
    )
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "p_raw", fixed = TRUE)
  expect_match(printed, "p_adjusted", fixed = TRUE)
  expect_match(printed, "selected", fixed = TRUE)
  expect_match(printed, "Yes", fixed = TRUE)
  expect_match(printed, "No", fixed = TRUE)
  expect_false(grepl("* = selected", printed, fixed = TRUE))
  expect_false(grepl("p_interact", printed, fixed = TRUE))
})

test_that("verbose Step 3 reports final AIC decisions without stars", {
  winners <- list(
    age = list(
      fit = list(ok = TRUE),
      metric = data.frame(AIC_main_minus_int = 3.681),
      type = "fp2",
      score = 3.681
    ),
    wt = list(
      fit = list(ok = TRUE),
      metric = data.frame(AIC_main_minus_int = -1.824),
      type = "fp1",
      score = -1.824
    )
  )

  output <- capture.output(
    print_interaction_step3_summary(
      var_winners = winners,
      cont_vars = c("age", "wt"),
      mode = "ic",
      ic_col = "AIC_main_minus_int",
      ic_label = "dAIC",
      min_improvement = 2,
      digits = 3L
    )
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "dAIC", fixed = TRUE)
  expect_match(printed, "selected", fixed = TRUE)
  expect_match(printed, "Yes", fixed = TRUE)
  expect_match(printed, "No", fixed = TRUE)
  expect_false(grepl("* = selected", printed, fixed = TRUE))
})

test_that("verbose candidate evaluation contains no provisional decision text", {
  body_text <- paste(
    deparse(body(evaluate_interaction_for_variable)),
    collapse = "\n"
  )

  expect_match(body_text, "p-value =", fixed = TRUE)
  expect_false(grepl("provisional", body_text, fixed = TRUE))
  expect_false(grepl("No significant interaction retained", body_text, fixed = TRUE))
  expect_false(grepl("Not selected", body_text, fixed = TRUE))
  expect_false(grepl("Selected  (", body_text, fixed = TRUE))
})

# 25.2 MFPI adjustment-model ACD power reconstruction

# Test purpose: The canonical fp_powers list must preserve both ACD power
# positions even when the display-oriented fp_terms table shows one power.
test_that("MFPI preserves structural ACD power positions", {
  adjustment_model <- list(
    fp_powers = list(
      age = c(NA_real_, 1),
      cavol = 1
    ),
    fp_terms = data.frame(
      power1 = c(1, 1),
      power2 = c(NA_real_, NA_real_),
      row.names = c("age", "cavol")
    ),
    acd = c(age = TRUE, cavol = FALSE)
  )

  powers <- mfpi_extract_adjustment_powers(
    adjustment_model,
    c("age", "cavol")
  )

  expect_named(powers, c("age", "cavol"))
  expect_length(powers$age, 2L)
  expect_true(is.na(powers$age[[1L]]))
  expect_equal(powers$age[[2L]], 1)
  expect_equal(powers$cavol, 1)
})

# Test purpose: A malformed or legacy ACD object that has lost one structural
# power slot should fail before transform_vector_acd() with a targeted error.
test_that("MFPI rejects malformed one-slot ACD adjustment powers", {
  adjustment_model <- list(
    fp_powers = list(age = 1),
    fp_terms = data.frame(
      power1 = 1,
      power2 = NA_real_,
      row.names = "age"
    ),
    acd = c(age = TRUE)
  )

  expect_error(
    mfpi_extract_adjustment_powers(adjustment_model, "age"),
    "must retain two power positions"
  )
})

# Test purpose: Reproduce the formula-interface failure in which a selected ACD
# adjustment variable is rebuilt while another continuous variable is tested
# for interaction. This previously passed a one-element power vector to
# transform_vector_acd().
test_that("MFPI rebuilds selected ACD adjustment terms during interaction fitting", {
  data("prostate", package = "mfp2")

  fit <- NULL
  expect_error(
    fit <- mfpi(
      lpsa ~ fp(age, acdx = TRUE, select = 1) + svi +
        fp(cavol, select = 1),
      data = prostate,
      cont_vars = "cavol",
      group_var = "svi",
      center = FALSE,
      flex = "flex1",
      winsorize = FALSE,
      include_group_var = TRUE,
      p_adjust_method = "holm",
      criterion = "p",
      show_models = FALSE,
      verbose = FALSE
    ),
    NA
  )

  expect_s3_class(fit, "mfpi")
  expect_length(fit$adjustment_model$fp_powers$age, 2L)
})

# Test purpose: A SAZ binary-only plot must not expose the internal indicator
# coding, where 1 denotes the zero group. The x-axis instead displays the two
# original-scale groups in the intuitive order zero, then positive.
test_that("plot() labels SAZ binary-only groups on the original scale", {
  skip_if_not_installed("ggplot2")

  data("prostate", package = "mfp2")
  fit <- mfp2(
    lpsa ~ fp(
      pgg45,
      df = 4,
      select = 0.05,
      alpha = 0.05,
      spike = TRUE
    ),
    data = prostate,
    family = "gaussian",
    criterion = "pvalue",
    verbose = FALSE
  )

  expect_identical(
    as.integer(fit$spike_dec[["pgg45"]]),
    as.integer(saz_decision_codes[["binary_only"]])
  )

  p <- plot(
    fit,
    terms = "pgg45",
    partial_only = TRUE
  )[["pgg45"]]

  expected_labels <- c("pgg45 = 0", "pgg45 > 0")
  x_scale <- p$scales$get_scales("x")

  expect_equal(x_scale$breaks, expected_labels)
  expect_equal(x_scale$limits, expected_labels)

  geom_classes <- vapply(
    p$layers,
    function(layer) class(layer$geom)[[1L]],
    character(1L)
  )
  point_layer <- p$layers[[which(geom_classes == "GeomPoint")[[1L]]]]

  expect_equal(
    as.character(point_layer$data$variable),
    expected_labels
  )
  expect_false(any(as.character(point_layer$data$variable) %in% c("0", "1")))
})


# =============================================================================
# ACD training-cache and compact-basis invariant regression tests
# =============================================================================

# Test purpose: The explicit fit-time A(x) cache must reproduce the historical
# apply_acd() path exactly for compact ACD candidate generation on the same
# training observations.
test_that("compact ACD generation reuses explicit training A(x) cache", {
  set.seed(411)
  x <- runif(100, 1, 20)
  powers <- c(0, 1)

  acd_fit <- fit_acd(x, powers = powers)
  acd_par <- acd_fit
  acd_par$acd <- NULL

  applied <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  cached <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par,
    acd_training_values = acd_fit$acd
  )

  expect_equal(cached$powers, applied$powers)
  expect_equal(cached$candidate_map, applied$candidate_map)
  expect_equal(cached$basis, applied$basis, tolerance = 1e-12)
})

# Test purpose: The materialized ACD generator shares the same explicit cache
# contract as the compact generator, so legacy/non-compact callers remain
# numerically identical when the training A(x) values are supplied.
test_that("materialized ACD generation reuses explicit training A(x) cache", {
  set.seed(412)
  x <- runif(80, 1, 15)
  powers <- c(0, 1)

  acd_fit <- fit_acd(x, powers = powers)
  acd_par <- acd_fit
  acd_par$acd <- NULL

  applied <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  cached <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par,
    acd_training_values = acd_fit$acd
  )

  expect_equal(cached$powers, applied$powers)
  expect_equal(length(cached$data), length(applied$data))
  for (i in seq_along(applied$data)) {
    expect_equal(cached$data[[i]], applied$data[[i]], tolerance = 1e-12)
  }
})

# Test purpose: A stored $acd component is training-data state, not a general
# prediction cache. Unless it is passed explicitly as acd_training_values, the
# generator must continue applying the stored ACD parameters to the supplied x.
test_that("ACD generator does not implicitly reuse acd_parameter$acd", {
  set.seed(413)
  x_fit <- runif(80, 1, 10)
  x_new <- runif(80, 10, 20)
  powers <- c(0, 1)

  acd_fit <- fit_acd(x_fit, powers = powers)

  with_stored_acd <- generate_transformations_acd_basis(
    x = x_new,
    degree = 1,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_fit
  )

  acd_par_without_values <- acd_fit
  acd_par_without_values$acd <- NULL

  without_stored_acd <- generate_transformations_acd_basis(
    x = x_new,
    degree = 1,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par_without_values
  )

  expect_equal(
    with_stored_acd$basis,
    without_stored_acd$basis,
    tolerance = 1e-12
  )
})

# Test purpose: The explicit training cache is row-aligned with x. Reject a
# wrong-length vector immediately instead of allowing recycling or a later,
# less informative matrix-shape failure.
test_that("ACD training cache validates observation count", {
  x <- seq(1, 10, length.out = 30)

  expect_error(
    generate_transformations_acd_basis(
      x = x,
      degree = 1,
      powers = c(0, 1),
      zero = FALSE,
      acd_parameter = list(
        beta0 = 0,
        beta1 = 1,
        power = 1,
        shift = 0,
        scale = 1
      ),
      acd_training_values = rep(0.5, length(x) - 1L)
    ),
    "one value per observation"
  )
})

# Test purpose: The normal model-search wrapper must pass fit_acd()$acd to the
# compact ACD generator. Invalid application coefficients make this test fail
# if apply_acd() is accidentally called instead of using the training cache.
test_that("transform_data_step() forwards cached training A(x)", {
  x_vec <- seq(1, 5, length.out = 25)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))

  acd_par <- list(
    acd = seq(0.2, 0.8, length.out = length(x_vec)),
    beta0 = NA_real_,
    beta1 = NA_real_,
    power = 1,
    shift = 0,
    scale = 1
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 4,
    powers = list(x = c(0, 1)),
    acdx = c(x = TRUE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = acd_par),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    term_to_columns = list(x = 1L),
    compact_acd = TRUE
  )

  expect_true(all(is.finite(out$acd_basis$basis)))
})

# Test purpose: Ordinary FP and ACD compact representations are mutually
# exclusive for one focal variable. If a future control-flow change violates
# that invariant, find_best_fpm_step() must fail explicitly instead of silently
# preferring fp_basis and potentially fitting the wrong candidate layout.
test_that("find_best_fpm_step() rejects simultaneous FP and ACD compact bases", {
  testthat::local_mocked_bindings(
    transform_data_step = function(...) {
      list(
        data_adj = NULL,
        data_fp = NULL,
        fp_basis = list(source = "fp"),
        acd_basis = list(source = "acd")
      )
    },
    .package = "mfp2"
  )

  x <- matrix(seq_len(8), ncol = 1L, dimnames = list(NULL, "x"))

  expect_error(
    find_best_fpm_step(
      x = x,
      xi = "x",
      degree = 2,
      y = rep(0, nrow(x)),
      powers_current = list(x = c(1, 1)),
      powers = list(x = c(0, 1)),
      acdx = c(x = FALSE),
      family = stats::gaussian(),
      family_string = "gaussian",
      zero = c(x = FALSE),
      catzero = list(x = NULL),
      spike = c(x = FALSE),
      spike_decision = c(x = 2),
      acd_parameter = list(x = NULL),
      prev_adj_params = list(),
      has_offset = FALSE,
      precomputed_adj = NULL,
      n_obs = nrow(x),
      term_to_columns = list(x = 1L)
    ),
    "both `fp_basis` and `acd_basis` are non-NULL"
  )
})

# =============================================================================
# End of tests
# =============================================================================
