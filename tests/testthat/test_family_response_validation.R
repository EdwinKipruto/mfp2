# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

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
