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
    verbose = FALSE
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


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 8. predict.mfp2()
# =============================================================================

# Test purpose: Checks default Gaussian predictions are finite and have one
# value per observation.
test_that("predict.mfp2() returns predictions for Gaussian model", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  preds <- predict(fit)
  expect_length(preds, nrow(x_prostate))
  expect_true(all(is.finite(preds)))
})


# Test purpose: Checks the structure and numerical correctness of Gaussian
# link-scale predictions and standard errors. The expected values are calculated
# independently as X beta and sqrt(diag(X V X')).
test_that("predict.mfp2() Gaussian fit and SE equal manual matrix calculation", {
  fit <- mfp2(
    x_prostate,
    y_prostate,
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  result <- predict(fit, type = "link", se.fit = TRUE)
  beta <- stats::coef(fit)
  beta_vcov <- stats::vcov(fit)
  manual_x <- cbind(`(Intercept)` = 1, x_prostate)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_fit <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_true(is.list(result))
  expect_named(result, c("fit", "se.fit", "residual.scale"))
  expect_equal(as.numeric(result$fit), manual_fit, tolerance = 1e-8)
  expect_equal(as.numeric(result$se.fit), manual_se, tolerance = 1e-8)
})


test_that("predict.mfp2() with newdata reproduces training predictions", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  preds_train <- predict(fit)
  preds_new <- predict(fit, newdata = x_prostate)

  expect_equal(as.numeric(preds_train), as.numeric(preds_new),
               tolerance = 1e-10)
})


# Test purpose: Checks that term-level predictions return per-term data frames with values and standard errors.
test_that("predict.mfp2() type = 'terms' returns list of data frames", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  terms_result <- predict(fit, type = "terms")
  expect_true(is.list(terms_result))

  for (nm in names(terms_result)) {
    expect_true(is.data.frame(terms_result[[nm]]))
    expect_true("value" %in% colnames(terms_result[[nm]]))
    expect_true("se" %in% colnames(terms_result[[nm]]))
  }
})


# Test purpose: Checks that contrast predictions return per-term data-frame outputs.
test_that("predict.mfp2() type = 'contrasts' returns list of data frames", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  contrasts_result <- predict(fit, type = "contrasts")
  expect_true(is.list(contrasts_result))

  for (nm in names(contrasts_result)) {
    expect_true(is.data.frame(contrasts_result[[nm]]))
  }
})


# Test purpose: Checks that prediction errors when newdata violates the positive domain required by a fitted log transform.
test_that("predict.mfp2() stops on domain violation in newdata", {
  x <- cbind(x1 = seq(1, 100, length.out = 100))
  # A small deterministic perturbation avoids an exact fit while keeping the
  # fitted transformation reproducible without depending on RNG state.
  y <- log(x[, "x1"]) + 0.01 * sin(seq_len(nrow(x)))

  fit <- mfp2(
    x, y,
    select = 1,
    alpha = 1,
    powers = list(x1 = 0),
    df = 2,
    shift = 0,
    scale = 1,
    verbose = FALSE
  )

  # create bad data
  bad_data <- x[1:5, , drop = FALSE]
  bad_data[, "x1"] <- -1

  expect_error(
    predict(fit, newdata = bad_data),
    "non-positive"
  )
})


# Test purpose: Checks that Cox-model predictions are finite and have the correct length.
test_that("predict.mfp2() works for Cox models", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes", "pgr", "er")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)

  preds <- predict(fit)
  expect_length(preds, nrow(x_gbsg))
  expect_true(all(is.finite(preds)))
})


# Test purpose: Ensures factor binomial responses do not break inherited
# predict.glm(se.fit = TRUE) behavior.
test_that("predict.mfp2() with se.fit = TRUE works for factor binomial response", {
  data("pima", package = "mfp2")

  x_pima <- as.matrix(pima[, c("glucose", "bmi", "age")])
  y_factor <- factor(pima$y, levels = c(0, 1), labels = c("no", "yes"))

  fit <- mfp2(
    x_pima,
    y_factor,
    family = "binomial",
    verbose = FALSE
  )

  result <- predict(fit, se.fit = TRUE)

  expect_true(is.list(result))
  expect_true("fit" %in% names(result))
  expect_true("se.fit" %in% names(result))
  expect_length(result$fit, nrow(x_pima))
  expect_length(result$se.fit, nrow(x_pima))
  expect_identical(fit$y_original, y_factor)
})


# Test purpose: Ensures prediction with newdata requires newoffset when the
# model was fitted with an offset.
test_that("predict.mfp2() requires newoffset when fitted model used offset", {
  set.seed(107)
  n <- 150

  x <- cbind(x1 = runif(n, 1, 10), x2 = runif(n, 1, 5))
  exposure <- runif(n, 0.5, 2)
  y <- rpois(n, exposure * exp(0.5 + 0.1 * x[, "x1"]))

  fit <- mfp2(
    x,
    y,
    family = "poisson",
    offset = log(exposure),
    verbose = FALSE
  )

  expect_error(
    predict(fit, newdata = x[1:5, , drop = FALSE]),
    "newoffset"
  )
})
