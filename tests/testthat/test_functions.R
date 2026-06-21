# Regression tests for spike-at-zero design matrix construction.
# These tests are intentionally focused on transform_matrix() rather than
# full model fitting, so they stay cheap and directly exercise the bug fixed
# in fit_mfp.R and predict.mfp2.R.

make_spike_design <- function(decision) {
  x <- matrix(c(0, 1, 2, 0, 4), ncol = 1L)
  colnames(x) <- "z"
  
  transform_matrix(
    x = x,
    power_list = list(z = c(1, 2)),
    center = c(z = FALSE),
    acdx = c(z = FALSE),
    zero = c(z = TRUE),
    catzero = c(z = TRUE),
    spike = c(z = TRUE),
    spike_decision = c(z = decision),
    check_binary = FALSE,
    reset_zero = FALSE
  )$x_transformed
}

test_that("spike decision 1 keeps FP columns plus z_bin", {
  m1 <- make_spike_design(1)
  
  expect_identical(colnames(m1), c("z.1", "z.2", "z_bin"))
  expect_identical(as.integer(m1[, "z_bin"]), c(1L, 0L, 0L, 1L, 0L))
})

test_that("spike decision 2 keeps FP columns only", {
  m2 <- make_spike_design(2)
  
  expect_identical(colnames(m2), c("z.1", "z.2"))
})

test_that("spike decision 3 keeps z_bin only", {
  m3 <- make_spike_design(3)
  
  expect_identical(colnames(m3), "z_bin")
  expect_identical(as.integer(m3[, "z_bin"]), c(1L, 0L, 0L, 1L, 0L))
})