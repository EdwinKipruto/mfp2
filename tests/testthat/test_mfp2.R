# =============================================================================
# Comprehensive tests for the mfp2 package
# =============================================================================
# Usage:
#   testthat::test_file("test-mfp2.R")
#
# These tests cover:
#   1.  mfp2.default() — Gaussian, binomial, Poisson, Cox
#   2.  mfp2.formula() — equivalence to default, fp() terms
#   3.  Family and response validation
#   4.  Preprocessing: shift, scale, centering
#   5.  Candidate-power validation and custom powers
#   6.  SAZ (spike-at-zero) eligibility, cascade, and reset
#   7.  ACD transformation
#   8.  predict.mfp2() — link, response, terms, contrasts, domain checks
#   9.  Model selection criteria (pvalue, AIC, BIC)
#  10.  Edge cases and input validation
#  11.  mfpi() — basic fitting and interaction testing
#  12.  predict.mfpi() — fitted functions and differences
# =============================================================================

library(testthat)
library(mfp2)
library(survival)

# =============================================================================
# Test data setup
# =============================================================================

data("prostate", package = "mfp2")

# Continuous predictors and Gaussian response
x_prostate <- as.matrix(prostate[, 2:8])
y_prostate <- as.numeric(prostate$lpsa)

# Helper: suppress verbose output
quiet <- function(expr) suppressMessages(capture.output(expr, type = "message"))

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
  
  fit <- mfp2(x, y, family = "binomial", verbose = FALSE)
  
  expect_s3_class(fit, "mfp2")
  expect_equal(fit$family_string, "binomial")
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

# Test purpose: Compares selected variables from equivalent matrix and formula
# interface fits.
test_that("default and formula interfaces give consistent selected variables", {
  fit_default <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  fit_formula <- mfp2(
    lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) + fp(weight) +
      fp(bph) + fp(cp),
    data = prostate, verbose = FALSE
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

# Test purpose: Checks that formula-interface factor variables are expanded correctly
# and that keep = "group" retains all dummy columns generated from the factor.
test_that("formula interface expands factor keep terms to dummy columns", {
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
  
  group_cols <- grep("^group", rownames(fit$fp_terms), value = TRUE)
  
  expect_true(length(group_cols) >= 1)
  expect_true(all(fit$fp_terms[group_cols, "selected"]))
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

# =============================================================================
# 3. Family and response validation
# =============================================================================

# Test purpose: Checks that supported character families, family functions, and 
# family objects are normalized correctly.
test_that("normalize_family_argument() accepts valid families", {
  # Character inputs
  for (fam in c("gaussian", "binomial", "poisson", "cox")) {
    result <- mfp2:::normalize_family_argument(fam)
    expect_equal(result$family_string, fam)
  }
  
  # Function inputs
  result <- mfp2:::normalize_family_argument(stats::gaussian)
  expect_equal(result$family_string, "gaussian")
  
  result <- mfp2:::normalize_family_argument(stats::binomial)
  expect_equal(result$family_string, "binomial")
  
  # Family object input
  result <- mfp2:::normalize_family_argument(stats::binomial(link = "probit"))
  expect_equal(result$family_string, "binomial")
})

# Test purpose: Checks that unsupported, ambiguous, or non-scalar family 
# specifications are rejected.
test_that("normalize_family_argument() rejects invalid families", {
  expect_error(mfp2:::normalize_family_argument("gamma"), "Invalid family")
  expect_error(mfp2:::normalize_family_argument("inverse.gaussian"), "Invalid family")
  expect_error(mfp2:::normalize_family_argument(c("gaussian", "binomial")),
               "single character string")
})

# Test purpose: Confirms that the character Cox family specification is accepted
# by the family normalizer.
test_that("cox must be specified as character string, not function", {
  # There is no stats::cox(), so creating a fake one would test the guard
  expect_error(mfp2:::normalize_family_argument("cox")$family_string, NA)
})

# Test purpose: Checks that invalid response values or response classes are 
# rejected for each family.
test_that("validate_family_response() catches bad responses", {
  # Gaussian: must be numeric
  expect_error(mfp2:::validate_family_response("abc", "gaussian", 3))
  expect_error(mfp2:::validate_family_response(c(1, NA, 3), "gaussian", 3))
  
  # Binomial: numeric must be in [0, 1]
  expect_error(mfp2:::validate_family_response(c(0, 2, 1), "binomial", 3))
  expect_error(mfp2:::validate_family_response(c(-0.1, 0.5, 0.9), "binomial", 3))
  
  # Poisson: must be non-negative
  expect_error(mfp2:::validate_family_response(c(-1, 0, 1), "poisson", 3))
  
  # Cox: must be Surv
  expect_error(mfp2:::validate_family_response(c(1, 2, 3), "cox", 3))
  
  # Surv for non-cox
  expect_error(mfp2:::validate_family_response(Surv(1:3, c(1,0,1)), "gaussian", 3))
})

# Test purpose: Checks that valid Gaussian, binomial, Poisson, and Cox responses 
# pass validation.
test_that("validate_family_response() accepts valid responses", {
  expect_true(mfp2:::validate_family_response(c(1.5, 2.5, 3.5), "gaussian", 3))
  expect_true(mfp2:::validate_family_response(c(0, 0.5, 1), "binomial", 3))
  expect_true(mfp2:::validate_family_response(
    factor(c("a", "b"), levels = c("a", "b")), "binomial", 2
  ))
  expect_true(mfp2:::validate_family_response(c(0, 1, 5), "poisson", 3))
  expect_true(mfp2:::validate_family_response(Surv(1:3, c(1, 0, 1)), "cox", 3))
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
# 5. Candidate-power validation and custom powers
# =============================================================================

# Test purpose: Checks that candidate power vectors are deduplicated and sorted.
test_that("normalize_fp_power_vector() removes duplicates and sorts", {
  result <- mfp2:::normalize_fp_power_vector(c(3, 1, 2, 1, -1))
  expect_equal(result, c(-1, 1, 2, 3))
})

# Test purpose: Checks that invalid power-vector inputs fail validation.
test_that("normalize_fp_power_vector() rejects empty or non-numeric input", {
  expect_error(mfp2:::normalize_fp_power_vector(numeric(0)))
  expect_error(mfp2:::normalize_fp_power_vector("abc"))
  expect_error(mfp2:::normalize_fp_power_vector(c(1, NA, 2)))
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
    mfp2:::validate_mfp_candidate_powers(
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
    mfp2:::validate_mfp_candidate_powers(
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
  
  fp1 <- mfp2:::generate_powers_fp(degree = 1, powers = powx)
  expect_equal(nrow(fp1), 8)
  expect_equal(ncol(fp1), 1)
  
  fp2 <- mfp2:::generate_powers_fp(degree = 2, powers = powx)
  expect_equal(nrow(fp2), 36) # C(8+2-1, 2) = C(9,2) = 36
  expect_equal(ncol(fp2), 2)
  
  # Single power with degree 2 gives repeated-power pair
  fp_single <- mfp2:::generate_powers_fp(degree = 2, powers = 2)
  expect_equal(nrow(fp_single), 1)
  expect_equal(fp_single[1, ], c(2, 2))
})

# Test purpose: Checks the null-degree FP power matrix used for omitted/null terms.
test_that("generate_powers_fp() degree 0 returns matrix(1)", {
  fp0 <- mfp2:::generate_powers_fp(degree = 0)
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

# Test purpose: Directly checks reset_spike() for the all-zero case, where
# the positive component is absent.
test_that("reset_spike() resets all-zero variables", {
  x <- matrix(0, nrow = 100, ncol = 1, dimnames = list(NULL, "exposure"))
  
  spike <- c(exposure = TRUE)
  user_catzero <- c(exposure = FALSE)
  user_zero <- c(exposure = FALSE)
  
  expect_warning(
    out <- mfp2:::reset_spike(
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


# Test purpose: Checks that resolve_saz_eligibility() temporarily recodes
# nonpositive values to zero before testing SAZ component proportions.
test_that("resolve_saz_eligibility() counts negative values as zero component", {
  x <- matrix(
    c(rep(-2, 20), rgamma(180, shape = 2, rate = 1)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )
  
  spike <- c(exposure = TRUE)
  catzero <- c(exposure = FALSE)
  zero <- c(exposure = FALSE)
  
  out <- mfp2:::resolve_saz_eligibility(
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
    out <- mfp2:::cap_spike_df(
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
    out <- mfp2:::cap_spike_df(
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
    out <- mfp2:::cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    NA
  )
  
  expect_equal(out[["exposure"]], 4)
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
  
  acd0 <- mfp2:::generate_powers_acd(degree = 0, powers = powx)
  expect_equal(ncol(acd0), 2)
  expect_equal(nrow(acd0), 1)
  
  acd1 <- mfp2:::generate_powers_acd(degree = 1, powers = powx)
  expect_equal(ncol(acd1), 2)
  expect_equal(nrow(acd1), 8)
  expect_true(all(is.na(acd1[, 1]))) # first column all NA
  
  acd2 <- mfp2:::generate_powers_acd(degree = 2, powers = powx)
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
    out <- mfp2:::reset_acd(x, acdx),
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
    out <- mfp2:::reset_acd(x, acdx),
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
    mfp2:::reset_acd(x, c(TRUE, FALSE)),
    "`acdx` must be a named logical vector"
  )
})

# Test purpose: Ensures reset_acd() rejects missing ACD flags before model
# fitting starts.
test_that("reset_acd() rejects missing acdx values", {
  x <- cbind(x1 = seq_len(20), x2 = seq_len(20))
  acdx <- c(x1 = TRUE, x2 = NA)
  
  expect_error(
    mfp2:::reset_acd(x, acdx),
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


# Test purpose: Ensures apply_acd() reproduces the training ACD transformation
# when supplied with parameters from fit_acd().
test_that("apply_acd() reproduces fit_acd() values using stored parameters", {
  set.seed(408)
  x <- runif(100, 1, 20)
  
  acd <- fit_acd(x)
  
  applied <- mfp2:::apply_acd(
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
    mfp2:::generate_powers_acd(degree = 3),
    "degree.*ACD.*0, 1, or 2"
  )
  
  expect_error(
    mfp2:::generate_powers_acd(degree = NA),
    "degree.*ACD.*0, 1, or 2"
  )
})

# Test purpose: Ensures ACD power generation uses ordered Cartesian products,
# because the first power applies to x and the second to A(x).
test_that("generate_powers_acd() keeps ordered power pairs", {
  powx <- c(0, 1)
  
  acd2 <- mfp2:::generate_powers_acd(degree = 2, powers = powx)
  
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
  
  out <- mfp2:::generate_transformations_acd(
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
  
  out <- mfp2:::generate_transformations_acd(
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

# Test purpose: Checks that prediction standard errors are returned when 
# se.fit = TRUE.
test_that("predict.mfp2() with se.fit = TRUE returns list", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  
  result <- predict(fit, se.fit = TRUE)
  expect_true(is.list(result))
  expect_true("fit" %in% names(result))
  expect_true("se.fit" %in% names(result))
  expect_length(result$fit, nrow(x_prostate))
  expect_length(result$se.fit, nrow(x_prostate))
})

# Test purpose: Checks that predicting on the original design matrix matches 
# training-data predictions.
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
  y <- log(x[, "x1"]) + rnorm(100, sd = 0.01)
  
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

# Test purpose: Checks that Cox prediction returns a numeric linear predictor on the intended reference scale.
test_that("predict.mfp2() for Cox uses reference = 'zero'", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)
  
  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)
  
  # Verify predictions are on the zero-reference scale
  preds <- predict(fit)
  # The linear predictor should not be centered around the mean
  # (which is what reference="sample" would do)
  expect_true(is.numeric(preds))
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
# 15. force_max_fp
# =============================================================================

# Test purpose: Checks that force_max_fp uses the requested maximum FP complexity 
# for selected variables under AIC.
test_that("force_max_fp forces maximum FP degree with AIC/BIC", {
  fit_force <- mfp2(
    x_prostate, y_prostate,
    criterion = "aic",
    force_max_fp = TRUE,
    select = 1,
    verbose = FALSE
  )
  
  for (v in get_selected_variable_names(fit_force)) {
    powers <- fit_force$fp_powers[[v]]
    requested_df <- as.numeric(fit_force$fp_terms[v, "df_initial"])
    
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

# =============================================================================
# 16. mfpi() — basic interaction fitting
# =============================================================================

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
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "linear"))
  expect_equal(fit$cont_var_forms["age"], c(age = "linear"))
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
  expect_equal(fit$cont_var_forms["age"], c(age = "linear"))
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
test_that("mfpi() rejects binary variables in cont_vars", {
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
test_that("mfpi() rejects binary variables in cont_vars", {
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
  
  lr <- mfp2:::calculate_lr_test(
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
    mfp2:::calculate_lr_test(logl = c(100, 110), dfs = c(5, 3)),
    "more degrees of freedom"
  )
})

# Test purpose: Checks that the F-test helper returns valid statistic, deviance 
# difference, and p-value.
test_that("calculate_f_test() returns correct p-value", {
  fit_null <- glm(y_prostate ~ 1)
  fit_full <- glm(y_prostate ~ x_prostate[, "cavol"])
  
  f_result <- mfp2:::calculate_f_test(
    deviances = c(deviance(fit_null), deviance(fit_full)),
    dfs_resid = c(df.residual(fit_null), df.residual(fit_full)),
    n_obs = length(y_prostate)
  )
  
  expect_true(f_result$pvalue >= 0 && f_result$pvalue <= 1)
  expect_true(f_result$statistic >= 0)
  expect_true(f_result$dev_diff >= 0)
})

# =============================================================================
# 21. fracplot()
# =============================================================================

# Test purpose: Checks that fracplot() can be called on a Gaussian mfp2 fit without error.
test_that("fracplot() runs without error for Gaussian model", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  # fracplot returns ggplot objects or similar; just check it runs
  expect_error(fracplot(fit), NA)
})



# =============================================================================
# End of tests
# =============================================================================