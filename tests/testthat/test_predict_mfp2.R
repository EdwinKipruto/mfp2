library(testthat)
library(mfp2)


make_shifted_log_prediction_fit <- function() {
  x <- cbind(x1 = seq(1, 100, length.out = 120))
  y <- log(x[, "x1"] + 5) + 0.01 * sin(seq_len(nrow(x)))

  mfp2(
    x,
    y,
    select = 1,
    alpha = 1,
    powers = list(x1 = 0),
    df = 2,
    shift = 5,
    scale = 1,
    center = FALSE,
    verbose = FALSE,
    warn_low_information = FALSE
  )
}


test_that("legacy fits without spike_dec support term and contrast predictions", {
  fit <- make_shifted_log_prediction_fit()

  # Simulate an object serialized before per-term SAZ decisions were stored.
  # Missing metadata means there is no evidence that this ordinary FP term is
  # binary-only; prediction should follow its stored continuous transformation.
  fit$spike_dec <- NULL
  newdata <- cbind(x1 = c(1, 2, 3))

  expect_false(mfp2:::prediction_term_is_binary_only(fit, "x1"))

  term_result <- expect_silent(
    predict(
      fit,
      newdata = newdata,
      type = "terms",
      terms = "x1",
      terms_seq = "data"
    )
  )
  expect_true(all(is.finite(term_result$x1$value)))

  contrast_result <- expect_silent(
    predict(
      fit,
      newdata = newdata,
      type = "contrasts",
      terms = "x1",
      terms_seq = "data",
      ref = list(x1 = 1)
    )
  )
  expect_true(all(is.finite(contrast_result$x1$value)))
})


test_that("binary-only prediction detection tolerates incomplete legacy metadata", {
  object <- list(spike_dec = NULL)
  expect_false(mfp2:::prediction_term_is_binary_only(object, "x"))

  object$spike_dec <- c(other = saz_decision_codes[["binary_only"]])
  expect_false(mfp2:::prediction_term_is_binary_only(object, "x"))

  object$spike_dec <- c(x = NA_integer_)
  expect_false(mfp2:::prediction_term_is_binary_only(object, "x"))

  object$spike_dec <- c(x = "malformed")
  expect_false(mfp2:::prediction_term_is_binary_only(object, "x"))

  object$spike_dec <- c(x = saz_decision_codes[["binary_only"]])
  expect_true(mfp2:::prediction_term_is_binary_only(object, "x"))
})


test_that("prediction rejects zero-row newdata before ranges or references", {
  fit <- make_shifted_log_prediction_fit()
  empty_df <- data.frame(x1 = numeric(0))
  empty_matrix <- matrix(
    numeric(0),
    nrow = 0L,
    ncol = 1L,
    dimnames = list(NULL, "x1")
  )

  # Full-model prediction reaches the same raw-newdata validator as the custom
  # term paths, so both documented data-frame and matrix inputs get one stable
  # diagnostic instead of backend-specific zero-length behavior.
  expect_error(
    predict(fit, newdata = empty_df, type = "link"),
    "must contain at least one row"
  )
  expect_error(
    predict(fit, newdata = empty_matrix, type = "response"),
    "must contain at least one row"
  )

  # Equidistant terms previously evaluated range(numeric(0)) and then seq() on
  # infinite endpoints. Data-sequence terms are covered separately because they
  # bypass range construction.
  expect_error(
    predict(
      fit,
      newdata = empty_df,
      type = "terms",
      terms = "x1",
      terms_seq = "equidistant"
    ),
    "must contain at least one row"
  )
  expect_error(
    predict(
      fit,
      newdata = empty_df,
      type = "terms",
      terms = "x1",
      terms_seq = "data"
    ),
    "must contain at least one row"
  )

  # With no explicit ref, contrasts formerly calculated mean(numeric(0)) and
  # silently obtained NaN. The row check must precede that default-reference
  # branch.
  expect_error(
    predict(
      fit,
      newdata = empty_df,
      type = "contrasts",
      terms = "x1",
      terms_seq = "data"
    ),
    "must contain at least one row"
  )
})


test_that("term and contrast predictions enforce the shifted FP domain", {
  fit <- make_shifted_log_prediction_fit()

  # Confirm the test premise: power zero represents log(x + shift).
  expect_equal(unname(fit$fp_powers[["x1"]]), 0)

  # These raw values become -1, 0, and 1 after applying the fitted shift.
  # Full-model, term, and contrast predictions must reject the same invalid
  # shifted domain instead of allowing NaN/Inf values into their results.
  # Keep newdata as a matrix to cover the documented matrix prediction
  # interface as well as the shifted-domain validation fixed by this test.
  bad_data <- cbind(x1 = c(-6, -5, -4))

  expect_error(
    predict(fit, newdata = bad_data, type = "link"),
    "non-positive"
  )

  expect_error(
    predict(
      fit,
      newdata = bad_data,
      type = "terms",
      terms = "x1",
      terms_seq = "data"
    ),
    "non-positive"
  )

  expect_error(
    predict(
      fit,
      newdata = bad_data,
      type = "terms",
      terms = "x1",
      terms_seq = "equidistant"
    ),
    "non-positive"
  )

  expect_error(
    predict(
      fit,
      type = "contrasts",
      terms = "x1",
      ref = list(x1 = -6)
    ),
    "non-positive"
  )

  # Positivity is checked after shifting, not on the raw values. Although these
  # raw values are negative, the fitted shift maps them to 1, 2, and 3, so term
  # predictions and a contrast against -4 must remain valid and finite.
  valid_data <- cbind(x1 = c(-4, -3, -2))

  term_result <- expect_silent(
    predict(
      fit,
      newdata = valid_data,
      type = "terms",
      terms = "x1",
      terms_seq = "data"
    )
  )
  expect_true(all(is.finite(
    as.matrix(term_result$x1[, c("value", "se", "lower", "upper")])
  )))

  contrast_result <- expect_silent(
    predict(
      fit,
      newdata = valid_data,
      type = "contrasts",
      terms = "x1",
      terms_seq = "data",
      ref = list(x1 = -4)
    )
  )
  expect_true(all(is.finite(
    as.matrix(contrast_result$x1[, c("value", "se", "lower", "upper")])
  )))
})


test_that("standard-error calculation clamps only numerical negative variance", {
  model <- list(
    family_string = "gaussian",
    transformed_to_model_columns = c(x = "x")
  )
  X <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(NULL, "x"))

  # A variance one machine epsilon below zero is compatible with floating-point
  # cancellation. It should yield an exact zero SE without sqrt() producing a
  # NaN or warning.
  testthat::local_mocked_bindings(
    vcov = function(object) {
      matrix(-.Machine$double.eps, 1L, 1L, dimnames = list("x", "x"))
    },
    .package = "mfp2"
  )

  expect_silent(
    expect_identical(
      mfp2:::calculate_standard_error(
        model,
        X,
        include_intercept = FALSE
      ),
      0
    )
  )
})


test_that("standard-error calculation rejects materially negative variance", {
  model <- list(
    family_string = "gaussian",
    transformed_to_model_columns = c(x = "x")
  )
  X <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(NULL, "x"))

  # This covariance block produces x' V x = -1, far beyond the numerical
  # tolerance. Returning NaN confidence limits would falsely suggest that a
  # standard error had been computed, so the calculation must stop clearly.
  testthat::local_mocked_bindings(
    vcov = function(object) {
      matrix(-1, 1L, 1L, dimnames = list("x", "x"))
    },
    .package = "mfp2"
  )

  expect_error(
    mfp2:::calculate_standard_error(
      model,
      X,
      include_intercept = FALSE
    ),
    "materially negative"
  )
})


test_that("standard-error calculation validates the required covariance block", {
  model <- list(
    family_string = "gaussian",
    transformed_to_model_columns = c(x = "x")
  )
  X <- matrix(2, nrow = 1L, ncol = 1L, dimnames = list(NULL, "x"))
  # diag() does not consistently preserve names from its input vector as
  # matrix dimnames. Supply the coefficient names explicitly so this fixture
  # has the same named-square-matrix contract as stats::vcov().
  selected_vcov <- matrix(
    c(4, 0, 0, NaN),
    nrow = 2L,
    dimnames = list(c("x", "unused"), c("x", "unused"))
  )

  testthat::local_mocked_bindings(
    vcov = function(object) selected_vcov,
    .package = "mfp2"
  )

  # The unused NaN belongs to another coefficient and must not invalidate the
  # requested x prediction. Its selected variance is 2^2 * 4 = 16.
  expect_silent(
    expect_equal(
      mfp2:::calculate_standard_error(
        model,
        X,
        include_intercept = FALSE
      ),
      4
    )
  )

  selected_vcov["x", "x"] <- NA_real_
  expect_error(
    mfp2:::calculate_standard_error(
      model,
      X,
      include_intercept = FALSE
    ),
    "non-finite entries for the required coefficients"
  )
})
