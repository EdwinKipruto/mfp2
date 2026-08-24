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


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 11. Summary, print, coef methods
# =============================================================================

# Test purpose: Checks that summary() returns output for a fitted Gaussian mfp2
#  model.
test_that("summary.mfp2() works for Gaussian", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  s <- summary(fit)
  expect_true(!is.null(s))
})


# Test purpose: Checks that summary() returns output for a fitted Cox mfp2 model.
test_that("summary.mfp2() works for Cox", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)
  s <- summary(fit)
  expect_true(!is.null(s))
})


# Test purpose: Ensures displayed ordinary FP equations use the final
# shifted-but-unscaled basis rather than reapplying the preprocessing scale.
test_that("summary FP labels do not reapply preprocessing scale", {
  object <- list(
    fp_terms = data.frame(
      selected = TRUE,
      acd = FALSE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      df_final = 2,
      power1 = -1,
      power2 = NA_real_,
      row.names = "x"
    ),
    x = matrix(1, nrow = 1L, dimnames = list(NULL, "x.1")),
    coefficients = c("x.1" = 2.5),
    transformations = data.frame(
      shift = 3,
      scale = 100,
      row.names = "x"
    )
  )

  classified <- mfp2_summary_classify_terms(object)
  label <- mfp2_summary_term_label(object, classified, "x", "x.1")
  formula <- mfp2_summary_formula_strings(object, classified)

  expect_true(grepl("x + 3", label, fixed = TRUE))
  expect_false(grepl("/100", label, fixed = TRUE))
  expect_false(any(grepl("/100", formula, fixed = TRUE)))
})


# Test purpose: Ensures ACD direct and transformed columns remain grouped under
# one variable and new ACD definitions use shifted, unscaled predictor values.
test_that("summary shows unscaled ACD definitions for new fits", {
  object <- list(
    fp_terms = data.frame(
      selected = TRUE,
      acd = TRUE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      df_final = 4,
      power1 = 1,
      power2 = 1,
      row.names = "x"
    ),
    x = matrix(
      c(1, 0.5),
      nrow = 1L,
      dimnames = list(NULL, c("x.1", "A_x.1"))
    ),
    coefficients = c("x.1" = 0.4, "A_x.1" = 1.2),
    transformations = data.frame(
      shift = 2,
      scale = 1,
      row.names = "x"
    ),
    acd_parameter = list(
      x = list(
        beta0 = -0.8,
        beta1 = 1.3,
        power = 0,
        shift = 0,
        scale = 1
      )
    )
  )

  classified <- mfp2_summary_classify_terms(object)
  direct_label <- mfp2_summary_term_label(object, classified, "x", "x.1")
  acd_label <- mfp2_summary_term_label(object, classified, "x", "A_x.1")
  definitions <- mfp2_summary_acd_definitions(object, classified)
  formulas <- mfp2_summary_formula_strings(object, classified)

  expect_identical(
    classified$cols_by_var[["x"]],
    c("x.1", "A_x.1")
  )
  expect_true(grepl("x + 2", direct_label, fixed = TRUE))
  expect_false(grepl("/", direct_label, fixed = TRUE))
  expect_true(grepl("A(x)", acd_label, fixed = TRUE))
  expect_false(grepl("log(A(x))", acd_label, fixed = TRUE))
  expect_true(any(grepl("x + 2", definitions, fixed = TRUE)))
  expect_false(any(grepl("/", definitions, fixed = TRUE)))
  expect_true(any(grepl("pnorm", definitions, fixed = TRUE)))
  expect_true(any(grepl("A(x)", formulas, fixed = TRUE)))
})


# Test purpose: Ensures ACD power positions are preserved so an ACD-only form
# c(NA, p) is not mislabeled as a direct ordinary FP term.
test_that("summary preserves ACD-only power slots", {
  object <- list(
    fp_terms = data.frame(
      selected = TRUE,
      acd = TRUE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      df_final = 2,
      power1 = NA_real_,
      power2 = 1,
      row.names = "x"
    ),
    x = matrix(0.5, nrow = 1L, dimnames = list(NULL, "A_x.1")),
    coefficients = c("A_x.1" = 1.2),
    transformations = data.frame(
      shift = 0,
      scale = 10,
      row.names = "x"
    ),
    acd_parameter = list(
      x = list(
        beta0 = 0,
        beta1 = 1,
        power = 1,
        shift = 0,
        scale = 10
      )
    )
  )

  classified <- mfp2_summary_classify_terms(object)
  label <- mfp2_summary_term_label(object, classified, "x", "A_x.1")

  expect_equal(classified$power_slots_by_var[["x"]], c(NA_real_, 1))
  expect_true(grepl("A(x)", label, fixed = TRUE))
  expect_false(grepl("log(A(x))", label, fixed = TRUE))
})


# Test purpose: Verifies that GLM Model Fit output uses stored deviances and a
# Deviance header rather than reconstructing minus twice log-likelihood.
test_that("Model Fit reports deviance for GLMs", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)
  values <- mfp2_summary_model_fit_values(fit)

  expect_identical(attr(values, "statistic_label"), "Deviance")
  expect_equal(
    values$fit_statistic,
    c(fit$linear_deviance, fit$mfp_deviance)
  )
  output <- capture.output(print(fit))
  expect_true(any(grepl("Deviance", output, fixed = TRUE)))
  expect_false(any(grepl("-2 log L", output, fixed = TRUE)))
})


# Test purpose: Verifies that Cox Model Fit output remains on the existing
# minus-twice-partial-log-likelihood scale.
test_that("Model Fit keeps -2 log L for Cox models", {
  data("gbsg", package = "mfp2")
  x_gbsg <- as.matrix(gbsg[, c("age", "size", "nodes")])
  y_gbsg <- Surv(gbsg$rectime, gbsg$censrec)

  fit <- mfp2(x_gbsg, y_gbsg, family = "cox", verbose = FALSE)
  values <- mfp2_summary_model_fit_values(fit)

  expect_identical(attr(values, "statistic_label"), "-2 log L")
  expect_equal(
    values$fit_statistic,
    c(fit$linear_deviance, fit$mfp_deviance)
  )
  output <- capture.output(print(fit))
  expect_true(any(grepl("-2 log L", output, fixed = TRUE)))
})


# Test purpose: Verifies that internal fast fits retain only lightweight
# quantities by default and calculate reporting statistics only when requested.
test_that("fit_model returns only requested internal components", {
  x <- matrix(
    c(-2, -1, 0, 1, 2, 3),
    ncol = 1L,
    dimnames = list(NULL, "x")
  )
  y <- c(-1.8, -0.9, 0.2, 1.1, 1.9, 3.2)

  selection_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = TRUE
  )

  expect_false("fit" %in% names(selection_fit))
  expect_false("null_deviance" %in% names(selection_fit))
  expect_false("model_deviance" %in% names(selection_fit))
  expect_false("sse" %in% names(selection_fit))
  expect_true(all(c("logl", "coefficients", "rank", "df") %in%
                    names(selection_fit)))

  metrics <- calculate_model_metrics(selection_fit, n_obs = length(y))
  expect_true(all(is.finite(metrics[c("logl", "df", "aic", "bic")])))

  reference_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = TRUE,
    calculate_fit_statistics = TRUE
  )

  expect_false("fit" %in% names(reference_fit))
  expect_true(all(c("null_deviance", "model_deviance") %in%
                    names(reference_fit)))

  retained_fast_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = TRUE,
    keep_fit = TRUE
  )

  expect_true("fit" %in% names(retained_fast_fit))

  full_fit <- fit_model(
    x = x,
    y = y,
    family = stats::gaussian(),
    family_string = "gaussian",
    fast = FALSE
  )

  expect_true("fit" %in% names(full_fit))
  expect_false("null_deviance" %in% names(full_fit))
  expect_false("model_deviance" %in% names(full_fit))
  expect_identical(
    unname(full_fit$transformed_to_model_columns),
    "x"
  )
})


# Test purpose: Verifies that negative-binomial summary inference agrees with
# MASS and uses z rather than t statistics after theta has been estimated.
test_that("summary.mfp2() reports correct negative-binomial inference", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  fits <- get_negbin_reference_fits()
  summary_mfp2 <- summary(fits$mfp2, raw = TRUE)
  summary_mass <- summary(fits$mass)

  expect_true(!is.null(summary_mfp2))
  expect_true(is.matrix(summary_mfp2$coefficients))
  expect_true(any(grepl("^z($| value)", colnames(summary_mfp2$coefficients))))
  expect_false(any(grepl("^t($| value)", colnames(summary_mfp2$coefficients))))
  expect_equal(
    unname(summary_mfp2$coefficients[, 1:2, drop = FALSE]),
    unname(summary_mass$coefficients[, 1:2, drop = FALSE]),
    tolerance = 1e-4
  )
  expect_equal(
    unname(summary_mfp2$coefficients[, 3:4, drop = FALSE]),
    unname(summary_mass$coefficients[, 3:4, drop = FALSE]),
    tolerance = 1e-3
  )
})


test_that("print and summary show the expanded df note unless notes is FALSE", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  direct_shown <- capture.output(
    print(fit, detailed_settings = FALSE, notes = TRUE)
  )
  direct_hidden <- capture.output(
    print(fit, detailed_settings = FALSE, notes = FALSE)
  )

  expect_true(any(grepl(
    "A retained catzero",
    direct_shown,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "binary-only SAZ uses 1 df",
    direct_shown,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "df counts fitted regression coefficients",
    direct_hidden,
    fixed = TRUE
  )))
  expect_true(any(grepl("Model Fit", direct_hidden, fixed = TRUE)))

  summary_shown <- summary(fit, notes = TRUE)
  summary_hidden <- summary(fit, notes = FALSE)

  expect_identical(summary_shown$notes, TRUE)
  expect_identical(summary_hidden$notes, FALSE)

  printed_summary_shown <- capture.output(print(summary_shown))
  printed_summary_hidden <- capture.output(print(summary_hidden))
  printed_summary_override <- capture.output(
    print(summary_shown, notes = FALSE)
  )

  expect_true(any(grepl(
    "A retained catzero",
    printed_summary_shown,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "df counts fitted regression coefficients",
    printed_summary_hidden,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "df counts fitted regression coefficients",
    printed_summary_override,
    fixed = TRUE
  )))
  expect_true(any(grepl("Model Fit", printed_summary_hidden, fixed = TRUE)))
})


test_that("notes must be one non-missing logical value", {
  fit <- mfp2(x_prostate, y_prostate, verbose = FALSE)

  expect_error(print(fit, notes = NA), "`notes`")
  expect_error(summary(fit, notes = NA), "`notes`")
  expect_error(print(summary(fit), notes = NA), "`notes`")
})
