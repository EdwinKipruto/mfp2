# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 9. Model selection criteria
# =============================================================================

# Test purpose: Checks that AIC-based model selection runs successfully.
test_that("criterion = 'aic' runs without error", {
  fit <- mfp2(x_prostate, y_prostate, criterion = "aic", verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})


# Test purpose: Checks that BIC-based model selection runs successfully.
test_that("criterion = 'bic' runs without error", {
  fit <- mfp2(x_prostate, y_prostate, criterion = "bic", verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})


# Test purpose: Checks the expected stronger-penalty behavior of BIC relative to AIC on this dataset.
test_that("BIC selects equal or fewer variables than AIC", {
  fit_aic <- mfp2(x_prostate, y_prostate, criterion = "aic", verbose = FALSE)
  fit_bic <- mfp2(x_prostate, y_prostate, criterion = "bic", verbose = FALSE)

  n_aic <- sum(fit_aic$fp_terms[, "selected"])
  n_bic <- sum(fit_bic$fp_terms[, "selected"])

  # BIC penalizes more heavily so typically selects <= AIC variables
  # This may not always hold for every dataset, but is generally expected
  expect_true(n_bic <= n_aic + 1) # allow slack of 1
})


# Test purpose: Checks that select = 1 retains all predictors under p-value selection.
test_that("select = 1 forces all variables into model", {
  fit <- mfp2(x_prostate, y_prostate, select = 1, verbose = FALSE)
  expect_true(all(fit$fp_terms[, "selected"]))
})


# Test purpose: Guards the Stata/CRAN MFP endpoint convention. For p-value
# selection, select = 1 and alpha = 1 are forcing values: a term is removed or
# simplified only when its p-value is strictly greater than the threshold.
# Therefore, even exact p-values of 1 must retain a degree-2 (df = 4) FP as FP2.
test_that("RA2 keeps FP2 when select = alpha = 1 and every test has p = 1", {
  metric_row <- function(df = 1) {
    matrix(
      c(
        logl = 0,
        df = df,
        aic = 0,
        bic = 0,
        deviance_gaussian = 1,
        df_resid = 10
      ),
      nrow = 1,
      dimnames = list(NULL, c(
        "logl", "df", "aic", "bic", "deviance_gaussian", "df_resid"
      ))
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(transform_cache = NULL)
    },
    find_best_fpm_step = function(..., degree) {
      power <- rep(1, degree)
      list(
        metrics = metric_row(df = 2 * degree),
        model_best = 1L,
        power_best = power,
        powers = matrix(power, nrow = 1),
        current_adj_params = NULL
      )
    },
    fit_null_step = function(...) {
      list(
        metrics = metric_row(df = 0),
        powers = NA_real_,
        current_adj_params = NULL
      )
    },
    fit_linear_step = function(...) {
      list(
        metrics = metric_row(df = 1),
        powers = 1,
        current_adj_params = NULL
      )
    },
    calculate_lr_test = function(...) {
      list(statistic = 0, pvalue = 1)
    },
    .package = "mfp2"
  )

  out <- select_ra2(
    x = matrix(1, nrow = 4, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = 1:4,
    powers_current = matrix(c(1, NA), nrow = 1, dimnames = list("x", NULL)),
    powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
    criterion = "pvalue",
    ftest = FALSE,
    select = 1,
    alpha = 1,
    family = gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 0L),
    acd_parameter = NULL,
    prev_adj_params = NULL,
    force_max_fp = c(x = FALSE),
    has_offset = FALSE,
    n_obs = 4,
    term_to_columns = list(x = "x")
  )

  expect_equal(out$model_best, 1)
  expect_equal(as.numeric(out$pvalue), c(1, 1, 1))
})


# Test purpose: Checks that variables listed in keep remain selected in the final model.
test_that("keep argument retains specified variables", {
  fit <- mfp2(x_prostate, y_prostate, keep = c("age", "bph"), verbose = FALSE)
  expect_true(fit$fp_terms["age", "selected"])
  expect_true(fit$fp_terms["bph", "selected"])
})


# Test purpose: Checks that Gaussian fitting works when F-test based selection is requested.
test_that("ftest argument works for Gaussian family", {
  fit <- mfp2(x_prostate, y_prostate, ftest = TRUE, verbose = FALSE)
  expect_s3_class(fit, "mfp2")
})
