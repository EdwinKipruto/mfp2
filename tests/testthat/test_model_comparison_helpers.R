# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

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


# Test purpose: Stata compatibility must refer its displayed negative-Wald
# selection difference to the chi-square distribution.
test_that("calculate_deviance_test() uses the deviance difference directly", {
  result <- calculate_deviance_test(
    deviances = c(reduced = -9, full = -15),
    dfs = c(reduced = 2, full = 4)
  )

  expect_equal(result$statistic, 6)
  expect_equal(result$dev_diff, 6)
  expect_equal(
    result$pvalue,
    stats::pchisq(6, df = 2, lower.tail = FALSE)
  )
})


test_that("calculate_deviance_test() rejects invalid comparisons", {
  expect_error(
    calculate_deviance_test(c(-9, -15), c(4, 2)),
    "more degrees of freedom"
  )
  expect_error(
    calculate_deviance_test(c(-9, NA_real_), c(2, 4)),
    "two finite deviances"
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
