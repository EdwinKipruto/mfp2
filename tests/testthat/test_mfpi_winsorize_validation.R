library(testthat)
library(mfp2)


test_that("Winsorisation rejects an all-missing continuous variable", {
  x <- cbind(age = rep(NA_real_, 5L))

  expect_error(
    mfp2:::winsorize_cont_vars(x, cont_vars = "age"),
    "variable `age`: it contains no non-missing values",
    fixed = TRUE
  )
})


test_that("Winsorisation leaves a single observed value unchanged", {
  x <- cbind(age = c(NA_real_, 5, NA_real_))

  result <- expect_silent(
    mfp2:::winsorize_cont_vars(x, cont_vars = "age")
  )

  expect_equal(result$x, x)
  expect_true(all(is.na(result$limits[, "age"])))
})


test_that("Winsorisation rejects coincident quantile limits", {
  # Both requested empirical quantiles equal zero. Applying them would replace
  # the sole upper-tail value with zero and collapse the variable completely.
  x <- cbind(age = c(rep(0, 9L), 1))

  expect_error(
    mfp2:::winsorize_cont_vars(
      x,
      cont_vars = "age",
      probs = c(0.25, 0.75)
    ),
    "lower and upper quantile limits are identical",
    fixed = TRUE
  )
})


test_that("Winsorisation retains ordinary distinct quantile limits", {
  x <- cbind(age = seq_len(100L))
  probs <- c(0.1, 0.9)
  expected_limits <- as.numeric(stats::quantile(x[, "age"], probs = probs))

  result <- expect_silent(
    mfp2:::winsorize_cont_vars(
      x,
      cont_vars = "age",
      probs = probs
    )
  )

  expect_equal(unname(result$limits[, "age"]), expected_limits)
  expect_equal(min(result$x[, "age"]), expected_limits[1L])
  expect_equal(max(result$x[, "age"]), expected_limits[2L])
  expect_gt(length(unique(result$x[, "age"])), 1L)
})
