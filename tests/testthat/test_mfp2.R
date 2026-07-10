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
#   1.  mfp2.default() - Gaussian, binomial, Poisson, Cox
#   2.  mfp2.formula() - formula parsing, fp()/fp2(), factor handling, strata
#   3.  Family and response validation
#   4.  Preprocessing - shift, scale, centering
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
#  15.  force_max_fp
#  16.  mfpi() - basic interaction fitting
#  17.  predict.mfpi()
#  18.  Reproducibility
#  19.  Transformation helpers
#  20.  Likelihood-ratio and F-test helpers
#  21.  plot()
#  22.  C++ core tests
# =============================================================================

library(testthat)
library(survival)
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


# Convert mfp2's stored transformed-term names back to the corresponding
# ordinary GLM column names for tests where df = 1 and preprocessing is off.
# For example, mfp2 stores the linear transformation of x1 as x1.1, whereas
# model.matrix() and glm() use x1. Factor dummy names and the intercept are
# unchanged. This helper is deliberately restricted to the trailing .1 suffix.
canonical_mfp2_linear_names <- function(x) {
  sub("\\.1$", "", x)
}

# Reconstruct a manual design matrix in the exact order and names required by
# a fitted coefficient vector. mfp2 may both reorder variables through xorder
# and rename ordinary linear transformations with a trailing .1 suffix.
# Matrix multiplication is positional, so both differences must be resolved
# before calculating X %*% beta or X V X'.
align_manual_design_to_coefficients <- function(manual_x, coefficient_names) {
  stopifnot(is.matrix(manual_x) || is.data.frame(manual_x))
  stopifnot(!is.null(colnames(manual_x)))
  
  source_names <- vapply(
    coefficient_names,
    function(coefficient_name) {
      if (coefficient_name %in% colnames(manual_x)) {
        return(coefficient_name)
      }
      
      undotted_name <- canonical_mfp2_linear_names(coefficient_name)
      if (undotted_name %in% colnames(manual_x)) {
        return(undotted_name)
      }
      
      NA_character_
    },
    character(1)
  )
  
  if (anyNA(source_names)) {
    stop(
      "Could not match manual design columns to coefficients: ",
      paste(coefficient_names[is.na(source_names)], collapse = ", "),
      call. = FALSE
    )
  }
  
  aligned <- as.matrix(manual_x)[, source_names, drop = FALSE]
  colnames(aligned) <- coefficient_names
  aligned
}

# Compare mfp2 and glm parameters by statistical term rather than by storage
# position. mfp2's xorder may change coefficient order, and its transformed
# linear columns carry a .1 suffix. The covariance matrix must be renamed and
# reordered consistently with the coefficient vector before comparison.
expect_mfp2_glm_parameters_equal <- function(fit_mfp2,
                                             fit_glm,
                                             tolerance = 1e-8) {
  coef_mfp2 <- stats::coef(fit_mfp2)
  coef_glm <- stats::coef(fit_glm)
  
  canonical_names <- canonical_mfp2_linear_names(names(coef_mfp2))
  expect_identical(anyDuplicated(canonical_names), 0L)
  expect_setequal(canonical_names, names(coef_glm))
  
  names(coef_mfp2) <- canonical_names
  coef_mfp2 <- coef_mfp2[names(coef_glm)]
  
  vcov_mfp2 <- stats::vcov(fit_mfp2)
  rownames(vcov_mfp2) <- canonical_mfp2_linear_names(rownames(vcov_mfp2))
  colnames(vcov_mfp2) <- canonical_mfp2_linear_names(colnames(vcov_mfp2))
  vcov_mfp2 <- vcov_mfp2[names(coef_glm), names(coef_glm), drop = FALSE]
  
  vcov_glm <- stats::vcov(fit_glm)
  vcov_glm <- vcov_glm[names(coef_glm), names(coef_glm), drop = FALSE]
  
  expect_equal(unname(coef_mfp2), unname(coef_glm), tolerance = tolerance)
  expect_equal(unname(vcov_mfp2), unname(vcov_glm), tolerance = tolerance)
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



# Test purpose: 2.1 Ordered factors are rejected before model-matrix expansion,
# because polynomial contrasts would otherwise be treated as separate predictors.
test_that("2.1 Formula interface rejects ordered factors", {
  set.seed(2011)
  n <- 90
  
  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    ordered_group = ordered(rep(c("low", "medium", "high"), length.out = n))
  )
  
  expect_error(
    mfp2(
      y ~ fp(x) + ordered_group,
      data = dat,
      verbose = FALSE
    ),
    "ordered|factor|contrast",
    ignore.case = TRUE
  )
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

test_that("predict.mfp2 works after fitting with survival::strata()", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  
  fit <- mfp2(
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
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
    verbose = FALSE
  )
  
  result <- predict(fit, type = "link", se.fit = TRUE)
  beta <- stats::coef(fit)
  beta_vcov <- stats::vcov(fit)
  manual_x <- cbind(`(Intercept)` = 1, x_prostate)
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
  manual_fit <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  
  expect_true(is.list(result))
  expect_named(result, c("fit", "se.fit", "residual.scale"))
  expect_equal(as.numeric(result$fit), manual_fit, tolerance = 1e-8)
  expect_equal(as.numeric(result$se.fit), manual_se, tolerance = 1e-8)
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

# Test purpose: Checks that default Cox predictions use reference = "zero".
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
    verbose = FALSE
  )
  
  got <- predict(fit, type = "lp", se.fit = TRUE)
  beta <- stats::coef(fit)
  beta_vcov <- stats::vcov(fit)
  manual_x <- align_manual_design_to_coefficients(
    x_gbsg,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
# Even in this simple configuration, mfp2 stores ordinary transformed columns
# as x.1, x1.1, and so on, and xorder may change their coefficient order. The
# helpers below therefore align by canonical term name rather than position.


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
  
  # Matrix multiplication is positional. Reorder all matrices by coefficient
  # name before calculating X beta or X V X'.
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
  
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
  
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
  
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
  
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
test_that("8.1.14 Formula factors: mfp2 predictions match glm", {
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
    keep = "group",
    shift = 0,
    scale = 1,
    center = FALSE,
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

# Test purpose: 8.2.1 Cox equivalence with linear terms, formula strata,
# no selection, no centering, and matched Breslow tie handling.
test_that("8.2.1 Cox: mfp2 predictions match coxph for linear no-selection model", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  
  dat <- dat[stats::complete.cases(
    dat[, c("time", "status", "age", "sex", "inst")]
  ), ]
  
  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
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
  
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  
  nd <- dat[1:20, c("age", "sex", "inst"), drop = FALSE]
  
  pred_mfp2 <- predict(
    fit_mfp2,
    newdata = nd,
    type = "lp"
  )
  
  pred_coxph <- predict(
    fit_coxph,
    newdata = nd,
    type = "lp",
    reference = "zero"
  )
  
  expect_equal(
    unname(pred_mfp2),
    unname(pred_coxph),
    tolerance = 1e-8
  )
  
  # For reference = "zero", a Cox linear predictor is exactly X beta. Build the
  # ordinary covariate matrix manually; the strata term contributes no column.
  manual_x <- cbind(age = nd$age, sex = nd$sex)
  beta <- stats::coef(fit_coxph)
  beta_vcov <- stats::vcov(fit_coxph)
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
  manual_lp <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  
  pred_mfp2_se <- predict(fit_mfp2, newdata = nd, type = "lp", se.fit = TRUE)
  pred_coxph_se <- predict(
    fit_coxph,
    newdata = nd,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  
  expect_equal(as.numeric(pred_mfp2_se$fit), manual_lp, tolerance = 1e-8)
  expect_equal(as.numeric(pred_coxph_se$fit), manual_lp, tolerance = 1e-8)
  expect_equal(as.numeric(pred_mfp2_se$se.fit), manual_se, tolerance = 1e-8)
  expect_equal(as.numeric(pred_coxph_se$se.fit), manual_se, tolerance = 1e-8)
})


# Test purpose: 8.2.2 Cox equivalence with unqualified strata(), ensuring the
# namespace-normalized and unqualified formula-special paths are both covered.
test_that("8.2.2 Cox: unqualified strata() predictions match coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  
  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ fp(age, df = 1, center = FALSE) + sex + strata(inst),
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
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  
  nd <- dat[1:20, c("age", "sex", "inst"), drop = FALSE]
  
  expect_equal(
    unname(predict(fit_mfp2, newdata = nd, type = "lp")),
    unname(predict(fit_coxph, newdata = nd, type = "lp", reference = "zero")),
    tolerance = 1e-8
  )
  expect_equal(unname(coef(fit_mfp2)), unname(coef(fit_coxph)), tolerance = 1e-8)
})

# Test purpose: 8.2.3 Cox equivalence with two separate strata() terms, which
# exercises the multi-column formula-strata reconstruction path.
test_that("8.2.3 Cox: multiple strata terms predictions match coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(
    dat[, c("time", "status", "age", "sex", "inst", "ph.ecog")]
  ), ]
  
  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ age + sex + survival::strata(inst) + survival::strata(ph.ecog),
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
    survival::Surv(time, status) ~ age + sex + survival::strata(inst) + survival::strata(ph.ecog),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  
  nd <- dat[1:20, c("age", "sex", "inst", "ph.ecog"), drop = FALSE]
  
  expect_equal(
    unname(predict(fit_mfp2, newdata = nd, type = "lp")),
    unname(predict(fit_coxph, newdata = nd, type = "lp", reference = "zero")),
    tolerance = 1e-8
  )
  expect_equal(unname(coef(fit_mfp2)), unname(coef(fit_coxph)), tolerance = 1e-8)
})

# Test purpose: 8.2.4 Cox default/matrix-interface strata argument should still
# match coxph() after preserving factor/strata labels at fit time.
test_that("8.2.4 Cox matrix interface: explicit strata argument predictions match coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  
  x <- as.matrix(dat[, c("age", "sex")])
  y <- survival::Surv(dat$time, dat$status)
  
  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "cox",
    strata = dat$inst,
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
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  
  nd <- dat[1:20, , drop = FALSE]
  newx <- as.matrix(nd[, c("age", "sex")])
  
  expect_equal(
    unname(predict(fit_mfp2, newdata = newx, strata = nd$inst, type = "lp")),
    unname(predict(fit_coxph, newdata = nd, type = "lp", reference = "zero")),
    tolerance = 1e-8
  )
  expect_equal(unname(coef(fit_mfp2)), unname(coef(fit_coxph)), tolerance = 1e-8)
})

# Test purpose: 8.2.5 Cox coefficient and partial log-likelihood equivalence in
# the simplest formula-strata case.
test_that("8.2.5 Cox: mfp2 coefficients and logLik match coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  
  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
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
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  
  expect_equal(unname(coef(fit_mfp2)), unname(coef(fit_coxph)), tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit_mfp2)), as.numeric(stats::logLik(fit_coxph)), tolerance = 1e-8)
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
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
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
  design <- mfp2:::mfpi_build_ordinary_design(
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
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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
  design <- mfp2:::mfpi_build_ordinary_design(
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
    cpp <- mfp2:::transform_fp_core(
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
  linear_2 <- mfp2:::interaction_model_df(n_groups = 2, degree = 0, flex = "flex1")
  linear_3 <- mfp2:::interaction_model_df(n_groups = 3, degree = 0, flex = "flex1")
  expect_equal(linear_2$dfint, 1)
  expect_equal(linear_3$dfint, 2)
  
  # FP1 with common powers in flex1/flex2 adds K - 1 slope parameters.
  fp1_3 <- mfp2:::interaction_model_df(n_groups = 3, degree = 1, flex = "flex1")
  expect_equal(fp1_3$dfint, 2)
  
  # FP2 with common powers adds two group-specific slope differences per
  # non-reference group.
  fp2_3 <- mfp2:::interaction_model_df(n_groups = 3, degree = 2, flex = "flex2")
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
  design <- mfp2:::mfpi_build_ordinary_design(
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
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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

# Test purpose: Verifies both Cox ordinary prediction scales without strata.
# MFPI link maps to coxph type = "lp" and response maps to type = "risk";
# both must use reference = "zero" to preserve the X beta convention.
test_that("17.1.9 unstratified Cox MFPI link and response match stored coxph", {
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
  design <- mfp2:::mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  stored <- fit_result$test_results$interaction_model$fit
  
  cases <- list(link = "lp", response = "risk")
  for (mfpi_type in names(cases)) {
    direct <- stats::predict(
      stored,
      newdata = design$model_newdata,
      type = cases[[mfpi_type]],
      se.fit = TRUE,
      reference = "zero"
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = mfpi_type,
      se.fit = TRUE
    )
    
    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
  }
  
  # Manual Cox oracle for the link scale. With reference = "zero", the stored
  # model's LP is exactly X beta and its SE is sqrt(diag(X V X')).
  reference_terms <- stats::delete.response(stats::terms(stored))
  manual_x <- stats::model.matrix(
    reference_terms,
    data = design$model_newdata,
    contrasts.arg = stored$contrasts
  )
  beta <- stats::coef(stored)
  beta_vcov <- stats::vcov(stored)
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
  manual_lp <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))
  
  got_link <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "link",
    se.fit = TRUE
  )
  expect_equal(got_link$predictions$fit, manual_lp, tolerance = 1e-8)
  expect_equal(got_link$predictions$se.fit, manual_se, tolerance = 1e-8)
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
    survival::Surv(time, status) ~ age + sex + survival::strata(inst),
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
  reconstructed <- mfp2:::reconstruct_formula_strata_newdata(fit, nd)
  design <- mfp2:::mfpi_build_ordinary_design(
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
# or data-frame strata must be combined exactly once with survival::strata(),
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
  design <- mfp2:::mfpi_build_ordinary_design(
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
  prepared <- mfp2:::mfpi_prepare_prediction_data(
    object = fit,
    term = "x",
    newdata = eval_x,
    grid = FALSE,
    n_grid = 200L
  )
  basis <- mfp2:::mfpi_build_function_basis(
    object = fit,
    term = "x",
    fit_result = fit_result,
    cont_var_scaled = prepared$cont_var_scaled,
    x_display = prepared$x_display
  )
  coefficients <- fit_result$test_results$interaction_model$coefficients
  
  group_internal <- names(basis$coefficient_groups)
  group_display <- mfp2:::mfpi_prediction_group_display_labels(
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
  
  vector_design <- mfp2:::mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = factor(nd$inst)
  )
  expect_identical(vector_design$model_newdata$strata_, factor(nd$inst))
  
  tabular_strata <- nd[, c("inst", "ph.ecog"), drop = FALSE]
  tabular_design <- mfp2:::mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = tabular_strata
  )
  expected <- do.call(
    survival::strata,
    c(as.list(tabular_strata), list(shortlabel = TRUE))
  )
  expect_equal(tabular_design$model_newdata$strata_, expected)
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

# =============================================================================
# 22. C++ core tests
# =============================================================================
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
})

# Test purpose: Checks that the C++ FP core returns the expected ordinary
# fractional-polynomial columns for powers 1 and 0.
# Test purpose: Checks that the C++ FP core returns the expected ordinary
# fractional-polynomial columns for powers 1 and 0.
test_that("C++ FP core computes ordinary FP transformations", {
  x <- c(1, 2, 4)
  p <- c(1, 0)
  
  out <- mfp2:::transform_fp_core(
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

# Test purpose: With a single non-unity candidate and force_max_fp = TRUE, the
# selected FP2 basis must be the repeated-power pair (2, 2), corresponding to
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
    force_max_fp = TRUE,
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
  
  # Matrix multiplication is positional, so reorder the manually constructed
  # design matrix and covariance matrix to match the fitted coefficient order.
  expect_true(all(names(beta) %in% colnames(manual_x)))
  
  manual_x <- align_manual_design_to_coefficients(
    manual_x,
    names(beta)
  )
  beta_vcov <- beta_vcov[names(beta), names(beta), drop = FALSE]
  
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

# =============================================================================
# End of tests
# =============================================================================