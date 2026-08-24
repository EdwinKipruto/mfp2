# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 13. Convergence and cycles
# =============================================================================

# Test purpose: Checks that the default maximum number of cycles is sufficient
# for convergence on the prostate data.
test_that("mfp2() converges within default cycles", {
  fit <- mfp2(x_prostate, y_prostate, cycles = 5, verbose = FALSE, warn_low_information = FALSE)
  expect_true(fit$convergence_mfp)
})


# Test purpose: Checks that a non-converged one-cycle fit warns but still returns
# an mfp2 object.
test_that("mfp2() with cycles = 1 still returns a result", {
  expect_warning(
    fit <- mfp2(x_prostate, y_prostate, cycles = 1, verbose = FALSE, warn_low_information = FALSE),
    "No convergence after 1 cycles"
  )

  expect_s3_class(fit, "mfp2")
})
