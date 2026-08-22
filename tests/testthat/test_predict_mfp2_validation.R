test_that("predict.mfp2 validates scalar prediction controls", {
  dummy <- structure(list(), class = "mfp2")

  bad_alpha <- list(
    NA_real_, Inf, -Inf, 0, 1, -0.1, 1.1,
    c(0.05, 0.10), numeric(0), "0.05"
  )
  for (value in bad_alpha) {
    expect_error(
      stats::predict(dummy, alpha = value),
      "alpha"
    )
  }

  bad_nseq <- list(
    NA_real_, Inf, -Inf, 0, -1, 2.5,
    c(10, 20), numeric(0), 1e20, "100"
  )
  for (value in bad_nseq) {
    expect_error(
      stats::predict(dummy, nseq = value),
      "nseq"
    )
  }

  bad_logical <- list(NA, c(TRUE, FALSE), logical(0), 1L, "TRUE", NULL)
  for (value in bad_logical) {
    expect_error(
      stats::predict(dummy, se.fit = value),
      "se.fit"
    )
    expect_error(
      stats::predict(dummy, add_intercept = value),
      "add_intercept"
    )
  }
})


test_that("shared scalar prediction validators accept valid boundary-safe values", {
  expect_invisible(mfp2:::validate_open_probability_scalar(0.05, "alpha"))
  expect_invisible(mfp2:::validate_open_probability_scalar(.Machine$double.eps, "alpha"))
  expect_invisible(mfp2:::validate_open_probability_scalar(1 - .Machine$double.eps, "alpha"))

  expect_invisible(mfp2:::validate_positive_integer_scalar(1, "nseq"))
  expect_invisible(mfp2:::validate_positive_integer_scalar(100L, "nseq"))

  expect_invisible(
    mfp2:::validate_logical_vector(TRUE, "se.fit", allowed_lengths = 1L)
  )
  expect_invisible(
    mfp2:::validate_logical_vector(FALSE, "add_intercept", allowed_lengths = 1L)
  )
})
