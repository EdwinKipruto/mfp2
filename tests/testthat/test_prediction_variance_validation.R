# The covariance fixtures below deliberately include invalid estimates. These
# tests exercise prediction-time validation, independently of fitting behavior.

test_that("shared prediction SEs preserve valid quadratic forms and row names", {
  X <- matrix(c(1, 2, -1, 1, 0, 3), ncol = 2L, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), c("x", "z")))
  V <- matrix(c(4, 1, 1, 9), 2L)
  # Independent per-row calculation includes the off-diagonal covariance.
  expected <- vapply(seq_len(nrow(X)), function(i) {
    sqrt(drop(t(X[i, ]) %*% V %*% X[i, ]))
  }, numeric(1L))
  names(expected) <- rownames(X)
  expect_equal(mfp2:::mfp2_prediction_standard_error(X, V), expected)
})

test_that("shared prediction SEs distinguish roundoff from invalid variance", {
  X <- matrix(1, 1L, 1L)
  expect_silent(expect_identical(
    mfp2:::mfp2_prediction_standard_error(X, matrix(-.Machine$double.eps)),
    0
  ))
  expect_error(
    mfp2:::mfp2_prediction_standard_error(X, matrix(-1)),
    "materially negative for 1 row"
  )

  # Cancellation tolerance must scale with the quadratic-form contributions.
  # Here the negative residual is tiny relative to the opposing 1e16 terms.
  X <- matrix(c(1, 1), 1L)
  V <- matrix(c(1e16, 0, 0, -1e16 - 2), 2L)
  expect_silent(expect_identical(
    mfp2:::mfp2_prediction_standard_error(X, V), 0
  ))
})

test_that("shared prediction SEs reject non-finite covariance and overflow", {
  X <- matrix(1, 1L, 1L)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(
      mfp2:::mfp2_prediction_standard_error(X, matrix(bad)),
      "non-finite entries for the required coefficients"
    )
  }
  # Finite inputs can still overflow during the quadratic-form calculation.
  expect_error(
    mfp2:::mfp2_prediction_standard_error(matrix(1e200), matrix(1)),
    "prediction variance contains non-finite values"
  )
})

test_that("GEE predictions validate covariance before the response delta method", {
  X <- matrix(c(1, 2), ncol = 1L, dimnames = list(NULL, "x"))
  object <- list(
    coefficients = c(x = 0.2),
    family = list(response_family = stats::poisson()),
    geese = list(X = X, gamma = 1),
    linear.predictors = c(0.2, 0.4)
  )
  V <- matrix(4, 1L, 1L, dimnames = list("x", "x"))
  testthat::local_mocked_bindings(
    vcov = function(object, ...) V,
    .package = "stats"
  )
  for (newx in list(NULL, X)) {
    for (type in c("link", "response")) {
      prediction <- function(se.fit = TRUE) {
        mfp2:::mfp2_predict_gee_matrix(
          object, newx = newx, type = type, se.fit = se.fit
        )
      }
      V[1L, 1L] <- 4
      result <- prediction()
      multiplier <- if (type == "response") exp(c(0.2, 0.4)) else 1
      expect_equal(unname(result$se.fit), c(2, 4) * multiplier)
      expect_equal(unname(result$fit), unname(prediction(FALSE)))

      V[1L, 1L] <- -.Machine$double.eps
      expect_silent(expect_equal(unname(prediction()$se.fit), c(0, 0)))
      V[1L, 1L] <- -1
      expect_error(prediction(), "materially negative")
      for (bad in c(NA_real_, NaN, Inf, -Inf)) {
        V[1L, 1L] <- bad
        expect_error(prediction(), "non-finite entries")
      }
      # Covariance validation applies only when standard errors are requested.
      expect_equal(unname(prediction(FALSE)), unname(result$fit))
    }
  }
})

test_that("fastglm predictions validate only the estimable covariance block", {
  X <- matrix(c(1, 0, 2, 0), ncol = 2L, byrow = TRUE,
              dimnames = list(NULL, c("x", "aliased")))
  object <- structure(list(
    coefficients = c(x = 0.2, aliased = NA_real_),
    x = X,
    linear.predictors = c(0.2, 0.4),
    family = stats::poisson(),
    dispersion = 2,
    cov.unscaled = matrix(c(2, NA_real_, NA_real_, NA_real_), 2L)
  ), class = "fastglm")
  for (newx in list(NULL, X)) {
    for (type in c("link", "response")) {
      prediction <- function(se.fit = TRUE) {
        mfp2:::mfp2_predict_fastglm_matrix(
          object, newx = newx, type = type, se.fit = se.fit
        )
      }
      object$cov.unscaled[1L, 1L] <- 2
      result <- prediction()
      multiplier <- if (type == "response") exp(c(0.2, 0.4)) else 1
      expect_equal(unname(result$se.fit), c(2, 4) * multiplier)
      expect_equal(unname(result$fit), unname(prediction(FALSE)))

      object$cov.unscaled[1L, 1L] <- -.Machine$double.eps
      expect_silent(expect_equal(unname(prediction()$se.fit), c(0, 0)))
      object$cov.unscaled[1L, 1L] <- -1
      expect_error(prediction(), "materially negative")
      for (bad in c(NA_real_, NaN, Inf, -Inf)) {
        object$cov.unscaled[1L, 1L] <- bad
        expect_error(prediction(), "non-finite entries")
      }
      expect_equal(unname(prediction(FALSE)), unname(result$fit))
    }
  }
})

test_that("multinomial terms and contrasts validate their required covariance", {
  # A balanced response avoids separation; fixed linear settings keep this
  # regression test focused on prediction rather than FP selection.
  dat <- data.frame(x = seq(1, 3, length.out = 90L),
                    y = factor(rep(c("a", "b", "c"), 30L)))
  fit <- mfp2(y ~ x, data = dat, family = "multinomial", df = 1,
              select = 1, alpha = 1, shift = 0, scale = 1,
              center = FALSE, cycles = 1, verbose = FALSE)
  original_covariance <- stats::vcov(fit)
  V <- original_covariance
  testthat::local_mocked_bindings(
    vcov = function(object, ...) V,
    .package = "stats"
  )
  for (type in c("terms", "contrasts")) {
    prediction <- function() {
      predict(fit, newdata = dat[1:3, ], type = type, terms = "x",
              terms_seq = "data", add_intercept = FALSE,
              ref = list(x = 0))$x
    }
    # Only the per-logit slope blocks are required here. Undefined intercept
    # covariance must not prevent an otherwise valid term/contrast prediction.
    slope <- !grepl(":(Intercept)", colnames(V), fixed = TRUE)
    V[,] <- 0
    diag(V) <- ifelse(slope, 4, NA_real_)
    result <- prediction()
    expect_equal(result$se, 2 * result$variable)
    diag(V) <- ifelse(slope, -.Machine$double.eps, NA_real_)
    expect_silent(expect_equal(prediction()$se, rep(0, nrow(result))))
    diag(V) <- ifelse(slope, -1, NA_real_)
    expect_error(prediction(), "materially negative")
    for (bad in c(NA_real_, NaN, Inf, -Inf)) {
      diag(V) <- ifelse(slope, bad, NA_real_)
      expect_error(prediction(), "non-finite entries")
    }
  }

  # The intercept becomes required when it is included in the reported term.
  V <- original_covariance
  intercept <- grepl(":(Intercept)", colnames(V), fixed = TRUE)
  diag(V)[intercept] <- NA_real_
  expect_error(
    predict(fit, type = "terms", terms = "x", add_intercept = TRUE),
    "non-finite entries"
  )
})
