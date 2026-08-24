# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Checks that the default interface requires a matrix design input.
test_that("mfp2() rejects non-matrix input in default interface", {
  expect_error(mfp2(as.data.frame(x_prostate), y_prostate), "matrix")
})


# Test purpose: Checks that missing predictor values are rejected.
test_that("mfp2() rejects input with missing data", {
  x_bad <- x_prostate
  x_bad[1, 1] <- NA
  expect_error(mfp2(x_bad, y_prostate), "NA")
})


# Test purpose: Checks that unnamed matrix columns are rejected to avoid ambiguous variable handling.
test_that("mfp2() rejects input without column names", {
  x_bad <- unname(x_prostate)
  expect_error(mfp2(x_bad, y_prostate), "column names")
})


# Test purpose: Checks that character-valued predictors are rejected.
test_that("mfp2() rejects character data in x", {
  x_bad <- x_prostate
  storage.mode(x_bad) <- "character"
  expect_error(mfp2(x_bad, y_prostate), "character")
})


# Test purpose: Checks that response length must match the number of predictor
# rows.
test_that("mfp2() rejects mismatched y length", {
  expect_error(mfp2(x_prostate, y_prostate[1:10]), "must match")
})


# Test purpose: Checks that numeric subsetting fits the model on the selected
# observations.
test_that("subset argument works correctly", {
  idx <- 1:50
  fit <- mfp2(x_prostate, y_prostate, subset = idx, verbose = FALSE, warn_low_information = FALSE)

  expect_s3_class(fit, "mfp2")
  # The model should be fitted on the subset
  expect_equal(length(fit$residuals), length(idx))
})


# Test purpose: Checks that logical subsetting fits the model on TRUE observations
#  only.
test_that("subset with logical vector works", {
  log_sub <- rep(FALSE, nrow(x_prostate))
  log_sub[1:50] <- TRUE
  fit <- mfp2(x_prostate, y_prostate, subset = log_sub, verbose = FALSE, warn_low_information = FALSE)

  expect_equal(length(fit$residuals), 50)
})


# Test purpose: Checks that all supported covariate-entry order options run
# successfully.
test_that("xorder options work without error", {
  for (ord in c("ascending", "descending", "original")) {
    fit <- mfp2(x_prostate, y_prostate, xorder = ord, verbose = FALSE, warn_low_information = FALSE)
    expect_s3_class(fit, "mfp2")
  }
})
