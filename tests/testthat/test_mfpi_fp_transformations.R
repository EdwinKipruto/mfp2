# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

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
