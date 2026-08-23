# fit_mfp() is reused by mfp2 and MFPI. These tests call it directly with only
# `criterion`: validation must run before the fitting engine accesses any other
# required argument, providing a stable error instead of a later switch/indexing
# failure.

test_that("fit_mfp validates criterion at the engine boundary", {
  expect_error(
    mfp2:::fit_mfp(criterion = "invalid"),
    "must be one of 'pvalue', 'aic', or 'bic'"
  )

  expect_error(
    mfp2:::fit_mfp(criterion = c("aic", "bic")),
    "must be one character value"
  )

  expect_error(
    mfp2:::fit_mfp(criterion = NA_character_),
    "must be one character value"
  )

  expect_error(
    mfp2:::fit_mfp(criterion = 1),
    "must be one character value"
  )

  expect_error(
    mfp2:::fit_mfp(criterion = ""),
    "must be one character value"
  )
})
