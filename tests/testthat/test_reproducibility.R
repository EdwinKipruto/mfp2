# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

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
