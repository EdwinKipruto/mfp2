# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

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
  fit <- mfp2(x_prostate, y_prostate, center = TRUE, verbose = FALSE, warn_low_information = FALSE)
  expect_true(!is.null(fit$centers))
})


# Test purpose: Checks that center = FALSE disables centering and leaves no
# centers stored.
test_that("centering can be disabled", {
  fit <- mfp2(x_prostate, y_prostate, center = FALSE, verbose = FALSE, warn_low_information = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_null(fit$centers)
  expect_true(all(fit$transformations$center == FALSE))
})
