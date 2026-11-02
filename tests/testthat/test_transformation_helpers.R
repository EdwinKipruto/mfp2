# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

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
