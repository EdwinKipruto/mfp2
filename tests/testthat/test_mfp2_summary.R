library(testthat)
library(mfp2)


test_that("Cox summary LRT preserves weights, ties, and strata", {
  n <- 8L
  object <- list(
    x = cbind(
      "x.1" = seq_len(n),
      "x.2" = seq_len(n)^2,
      "z.1" = rep(c(0, 1), length.out = n)
    ),
    y = cbind(
      time = seq_len(n),
      status = rep(c(1, 0), length.out = n)
    ),
    family = "cox",
    family_string = "cox",
    fitter = "base",
    weights = c(0.5, 1, 1.5, 2, 0.75, 1.25, 1.75, 2.25),
    # Cox fits store case weights in `weights`, not `prior.weights`.
    prior.weights = rep(99, n),
    offset = rep(0, n),
    method = "breslow",
    strata = factor(rep(c("A", "B"), each = n / 2L))
  )

  classified <- list(
    variable_names = "x",
    df_final = 2,
    cols_by_var = list(x = c("x.1", "x.2"))
  )

  refit_calls <- list()
  testthat::local_mocked_bindings(
    fit_model = function(x,
                         weights = NULL,
                         method = NULL,
                         strata = NULL,
                         ...) {
      refit_calls[[length(refit_calls) + 1L]] <<- list(
        x = x,
        weights = weights,
        method = method,
        strata = strata
      )

      list(logl = if (ncol(x) == 3L) -100 else -108)
    },
    .package = "mfp2"
  )

  result <- mfp2:::mfp2_summary_lrt_drop_variable(
    object = object,
    classified = classified,
    v = "x"
  )

  expect_length(refit_calls, 2L)
  for (refit_call in refit_calls) {
    expect_identical(refit_call$weights, object$weights)
    expect_identical(refit_call$method, object$method)
    expect_identical(refit_call$strata, object$strata)
  }

  expect_identical(colnames(refit_calls[[1L]]$x), c("x.1", "x.2", "z.1"))
  expect_identical(colnames(refit_calls[[2L]]$x), "z.1")
  expect_equal(result$lr, 16)
  expect_equal(result$df, 2)
  expect_equal(result$p, stats::pchisq(16, df = 2, lower.tail = FALSE))
})


test_that("coefficient columns use exact GLM labels", {
  object <- list(coefficients = c("(Intercept)" = 1, x = 2))
  cmat <- cbind(
    "exp(coef)" = c(99, 98),
    "Estimate" = c(1, 2),
    "Std. Error" = c(0.1, 0.2),
    "t value" = c(10, 10),
    "Pr(>|t|)" = c(0.01, 0.02)
  )
  rownames(cmat) <- names(object$coefficients)

  result <- mfp2:::mfp2_summary_coef_matrix(
    object,
    raw_summary = list(coefficients = cmat)
  )

  # Exact matching must not confuse the deliberately preceding `exp(coef)`
  # column with the estimate column.
  expect_equal(unname(result[, "estimate"]), c(1, 2))
  expect_equal(unname(result[, "se"]), c(0.1, 0.2))
  expect_equal(unname(result[, "statistic"]), c(10, 10))
  expect_equal(unname(result[, "pvalue"]), c(0.01, 0.02))
})


test_that("robust Cox standard errors take exact-name precedence", {
  object <- list(coefficients = c(x = 0.5))
  cmat <- matrix(
    c(0.5, exp(0.5), 0.4, 0.2, 2.5, 0.012),
    nrow = 1L,
    dimnames = list(
      "x",
      c("coef", "exp(coef)", "se(coef)", "robust se", "z", "Pr(>|z|)")
    )
  )

  result <- mfp2:::mfp2_summary_coef_matrix(
    object,
    raw_summary = list(coefficients = cmat)
  )

  # summary.coxph() bases the displayed z and p-value on `robust se` when that
  # column exists, so the standardized table must select the same SE.
  expect_equal(unname(result[, "estimate"]), 0.5)
  expect_equal(unname(result[, "se"]), 0.2)
  expect_equal(unname(result[, "statistic"]), 2.5)
  expect_equal(unname(result[, "pvalue"]), 0.012)
})


test_that("unknown or ambiguous coefficient headers use the vcov fallback", {
  fit <- stats::lm(mpg ~ wt, data = mtcars)
  raw_summary <- summary(fit)
  colnames(raw_summary$coefficients) <- c(
    "estimate changed",
    "standard error changed",
    "statistic changed",
    "probability changed"
  )

  result <- mfp2:::mfp2_summary_coef_matrix(fit, raw_summary)
  expected_se <- sqrt(diag(stats::vcov(fit)))
  expected_statistic <- stats::coef(fit) / expected_se

  # Unknown future headers must not be used as NA indices. The existing
  # covariance path reconstructs finite statistics instead.
  expect_equal(result[, "estimate"], stats::coef(fit))
  expect_equal(result[, "se"], expected_se)
  expect_equal(result[, "statistic"], expected_statistic)
  expect_equal(
    result[, "pvalue"],
    2 * stats::pnorm(abs(expected_statistic), lower.tail = FALSE)
  )
  expect_false(anyNA(result))

  # Duplicate exact labels are ambiguous and must take the same safe path.
  ambiguous <- summary(fit)
  colnames(ambiguous$coefficients)[1:2] <- c("Estimate", "Estimate")
  ambiguous_result <- mfp2:::mfp2_summary_coef_matrix(fit, ambiguous)
  expect_equal(ambiguous_result, result)
})
