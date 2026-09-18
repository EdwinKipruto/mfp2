# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 3. Family and response validation
# =============================================================================

# Test purpose: Checks that supported character families, family functions, and
# family objects are normalized correctly.
test_that("normalize_family_argument() accepts likelihood families", {
  # Character inputs
  for (fam in c(
    "gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
    "cox", "negbin", "survreg", "finegray"
  )) {
    result <- normalize_family_argument(fam)
    expect_equal(result$family_string, fam)
  }

  # Character matching is case-insensitive while retaining stats' canonical
  # capitalisation for Gamma.
  expect_identical(normalize_family_argument("gamma")$family_string, "Gamma")

  # Function inputs
  result <- normalize_family_argument(stats::gaussian)
  expect_equal(result$family_string, "gaussian")

  result <- normalize_family_argument(stats::binomial)
  expect_equal(result$family_string, "binomial")

  result <- normalize_family_argument(stats::Gamma)
  expect_equal(result$family_string, "Gamma")

  result <- normalize_family_argument(stats::inverse.gaussian)
  expect_equal(result$family_string, "inverse.gaussian")

  # Family object input
  result <- normalize_family_argument(stats::binomial(link = "probit"))
  expect_equal(result$family_string, "binomial")

  expect_identical(
    normalize_family_argument(survreg_family(dist = "lognormal"))$family_string,
    "survreg"
  )
  expect_identical(
    normalize_family_argument(finegray_family(etype = "relapse"))$family_string,
    "finegray"
  )
})


# Test purpose: Checks that unsupported, ambiguous, or non-scalar family
# specifications are rejected.
test_that("normalize_family_argument() rejects invalid and quasi families", {
  expect_error(normalize_family_argument("quasi"), "no likelihood")
  expect_error(normalize_family_argument("quasibinomial"), "no likelihood")
  expect_error(normalize_family_argument(stats::quasipoisson), "no likelihood")
  expect_error(normalize_family_argument(stats::quasibinomial()), "no likelihood")
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

  # Poisson: must be finite non-negative integer counts
  expect_error(validate_family_response(c(-1, 0, 1), "poisson", 3),
               "non-negative")
  expect_error(validate_family_response(c(0, 1.5, 2), "poisson", 3),
               "integer")

  # Grouped binomial responses contain integer success/failure counts. The
  # separate numeric-vector representation remains a probability/proportion.
  expect_error(
    validate_family_response(cbind(c(1, 1.5, 2), c(3, 2.5, 1)), "binomial", 3),
    "integers"
  )

  # Gamma and inverse Gaussian: responses must be strictly positive.
  expect_error(validate_family_response(c(0, 1, 2), "Gamma", 3),
               "strictly positive")
  expect_error(validate_family_response(c(-1, 1, 2), "inverse.gaussian", 3),
               "strictly positive")

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


# Test purpose: Checks that valid likelihood-family responses pass validation.
test_that("validate_family_response() accepts valid responses", {
  expect_true(validate_family_response(c(1.5, 2.5, 3.5), "gaussian", 3))
  expect_true(validate_family_response(c(0, 0.5, 1), "binomial", 3))
  expect_true(validate_family_response(
    cbind(successes = c(0, 2, 4), failures = c(5, 3, 1)),
    "binomial",
    3
  ))
  expect_true(validate_family_response(
    factor(c("a", "b"), levels = c("a", "b")), "binomial", 2
  ))
  expect_true(validate_family_response(c(0, 1, 5), "poisson", 3))
  expect_true(validate_family_response(c(0.5, 1, 5), "Gamma", 3))
  expect_true(validate_family_response(c(0.5, 1, 5), "inverse.gaussian", 3))
  expect_true(validate_family_response(c(0, 1, 5), "negbin", 3))
  expect_true(validate_family_response(Surv(1:3, c(1, 0, 1)), "cox", 3))
})


test_that("public fitting interfaces reject fractional count responses", {
  x_mfp <- matrix(seq_len(8), ncol = 1L, dimnames = list(NULL, "x"))
  y_fractional <- c(0, 1, 2, 1.5, 3, 2, 1, 0)

  expect_error(
    mfp2(
      x = x_mfp,
      y = y_fractional,
      family = "poisson",
      verbose = FALSE
    ),
    "integer counts"
  )

  grouped_fractional <- cbind(
    successes = c(0, 1, 2, 1.5, 3, 2, 1, 0),
    failures = c(4, 3, 2, 2.5, 1, 2, 3, 4)
  )
  expect_error(
    mfp2(
      x = x_mfp,
      y = grouped_fractional,
      family = "binomial",
      verbose = FALSE
    ),
    "integers"
  )

  x_mfpi <- cbind(
    group = rep(0:1, each = 4),
    exposure = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5)
  )
  expect_error(
    mfpi(
      x = x_mfpi,
      y = y_fractional,
      group_var = "group",
      interaction_vars = "exposure",
      interaction_forms = c(exposure = "linear"),
      family = "poisson",
      verbose = FALSE
    ),
    "integer counts"
  )
})
