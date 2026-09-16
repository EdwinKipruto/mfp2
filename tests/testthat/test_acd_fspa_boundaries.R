# Regression tests for the ACD function-selection procedure (FSPA).


acd_fspa_mock_fit <- function(pvalues, alpha = 0.05) {
  metric_row <- function(df = 1) {
    matrix(
      c(
        logl = 0,
        df = df,
        aic = 0,
        bic = 0,
        deviance_gaussian = 1,
        df_resid = 20
      ),
      nrow = 1,
      dimnames = list(NULL, c(
        "logl", "df", "aic", "bic", "deviance_gaussian", "df_resid"
      ))
    )
  }

  state <- new.env(parent = emptyenv())
  state$test <- 0L

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(transform_cache = NULL)
    },
    find_best_fpm_step = function(..., degree, acdx) {
      is_acd <- isTRUE(unname(acdx[["x"]]))
      selected_power <- if (degree == 2L) {
        c(1, 1)
      } else if (is_acd) {
        c(NA_real_, 1)
      } else {
        1
      }

      list(
        metrics = metric_row(df = 2 * degree),
        model_best = 1L,
        power_best = selected_power,
        powers = matrix(selected_power, nrow = 1L),
        current_adj_params = NULL
      )
    },
    fit_null_step = function(...) {
      list(
        metrics = metric_row(df = 0),
        powers = c(NA_real_, NA_real_),
        current_adj_params = NULL
      )
    },
    fit_linear_step = function(..., acdx) {
      is_acd <- isTRUE(unname(acdx[["x"]]))
      selected_power <- if (is_acd) c(NA_real_, 1) else 1
      list(
        metrics = metric_row(df = 1),
        powers = selected_power,
        current_adj_params = NULL
      )
    },
    calculate_lr_test = function(...) {
      state$test <- state$test + 1L
      list(statistic = 0, pvalue = pvalues[[state$test]])
    },
    .package = "mfp2"
  )

  result <- select_ra2_acd(
    x = matrix(1:8, nrow = 8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = TRUE),
    y = seq_len(8),
    powers_current = matrix(c(1, 1), nrow = 1, dimnames = list("x", NULL)),
    powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
    criterion = "pvalue",
    ftest = FALSE,
    select = 0.05,
    alpha = alpha,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 0L),
    acd_parameter = NULL,
    prev_adj_params = NULL,
    force_max_fp = c(x = FALSE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  list(result = result, tests_run = state$test)
}


# Equality is significant for retention: simplification is permitted only for
# a valid p-value strictly greater than the threshold. This also makes
# alpha = 1 a true forcing value, consistently with ordinary RA2.
test_that("ACD FSPA Test 4 retains M1 when p equals alpha", {
  out <- acd_fspa_mock_fit(c(0, 0, 0, 1), alpha = 1)

  expect_equal(out$result$model_best, 1L)
  expect_equal(out$tests_run, 4L)
  expect_equal(unname(out$result$pvalue[[4L]]), 1)
})


test_that("ACD FSPA Test 5 retains M3 when p equals alpha", {
  out <- acd_fspa_mock_fit(c(0, 0, 0, 0.50, 0.05), alpha = 0.05)

  expect_equal(out$result$model_best, 5L)
  expect_equal(out$tests_run, 5L)
  expect_equal(unname(out$result$pvalue[[5L]]), 0.05)
})


# A missing, non-finite, or out-of-domain p-value is not affirmative evidence
# for moving to a simpler model. Both late FSPA branches must fail safely.
test_that("ACD FSPA retains complexity for indeterminate late p-values", {
  for (bad_p in list(NA_real_, NaN, Inf, -Inf)) {
    at_test4 <- acd_fspa_mock_fit(c(0, 0, 0, bad_p), alpha = 0.05)
    expect_equal(at_test4$result$model_best, 1L)
    expect_equal(at_test4$tests_run, 4L)

    at_test5 <- acd_fspa_mock_fit(c(0, 0, 0, 0.50, bad_p), alpha = 0.05)
    expect_equal(at_test5$result$model_best, 5L)
    expect_equal(at_test5$tests_run, 5L)
  }
})
