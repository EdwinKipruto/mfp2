test_that("FP1 reuses a fixed linear fit at power 1 with FP selection df", {
  calls <- 0L
  powers <- matrix(c(-1, 1, 2), ncol = 1L)
  linear_metrics <- c(logl = 10, df = 2, deviance_rs = -20,
                      deviance_gaussian = NA_real_, aic = -16,
                      bic = -20 + 2 * log(20), df_resid = 18)

  testthat::local_mocked_bindings(
    transform_data_step = function(...) list(
      data_adj = NULL,
      data_fp = lapply(powers[, 1L], function(p) matrix(rep(p, 20), ncol = 1L)),
      fp_basis = NULL, acd_basis = NULL, powers_fp = powers,
      current_params = list(x = list())
    ),
    fit_model = function(...) {
      calls <<- calls + 1L
      list()
    },
    calculate_model_metrics = function(...) {
      c(logl = 0, df = 3, deviance_rs = 0,
        deviance_gaussian = NA_real_, aic = 6,
        bic = 3 * log(20), df_resid = 17)
    },
    .package = "mfp2"
  )

  fit <- find_best_fpm_step(
    x = matrix(seq_len(20), ncol = 1L, dimnames = list(NULL, "x")),
    xi = "x", degree = 1, y = seq_len(20),
    powers_current = list(x = 1), powers = list(x = powers[, 1L]),
    acdx = c(x = FALSE), family = NULL, family_string = "gaussian",
    zero = c(x = FALSE), catzero = list(x = NULL), spike = c(x = FALSE),
    spike_decision = c(x = 2), acd_parameter = list(x = NULL),
    prev_adj_params = list(), has_offset = FALSE,
    precomputed_adj = list(), n_obs = 20,
    term_to_columns = list(x = 1L),
    linear_fit = list(metrics = rbind(linear_metrics))
  )

  expect_equal(calls, 2L)
  expect_equal(fit$powers[, 1L], c(-1, 1, 2))
  expect_equal(fit$model_best, 2L)
  expect_equal(unname(fit$metrics[2L, "logl"]), 10)
  expect_equal(unname(fit$metrics[2L, "df"]), 3)
  expect_equal(unname(fit$metrics[2L, "aic"]), -14)
  expect_equal(unname(fit$metrics[2L, "df_resid"]), 17)
})

test_that("linear rung reuses power 1 when FP1 was searched first", {
  powers <- matrix(c(-1, 1, 2), ncol = 1L)
  fp_metrics <- matrix(0, nrow = 3L, ncol = 7L,
                       dimnames = list(NULL, c("logl", "df", "deviance_rs",
                                               "deviance_gaussian", "aic", "bic",
                                               "df_resid")))
  fp_metrics[2L, ] <- c(10, 3, -20, NA, -14,
                        -20 + 3 * log(20), 17)

  testthat::local_mocked_bindings(
    transform_data_step = function(...) list(
      data_fp = list(matrix(seq_len(20), ncol = 1L)),
      powers_fp = matrix(1, ncol = 1L),
      current_params = list(x = list())
    ),
    fit_model = function(...) stop("power 1 was fitted twice"),
    .package = "mfp2"
  )

  fit <- fit_linear_step(
    x = matrix(seq_len(20), ncol = 1L, dimnames = list(NULL, "x")),
    xi = "x", y = seq_len(20), powers_current = list(x = 1),
    powers = list(x = powers[, 1L]), acdx = c(x = FALSE),
    family = NULL, family_string = "gaussian", zero = c(x = FALSE),
    catzero = list(x = NULL), spike = c(x = FALSE),
    spike_decision = c(x = 2), acd_parameter = list(x = NULL),
    prev_adj_params = list(), has_offset = FALSE, n_obs = 20,
    precomputed_adj = list(), term_to_columns = list(x = 1L),
    fp1_fit = list(powers = powers, metrics = fp_metrics)
  )

  expect_equal(unname(fit$metrics[1L, "logl"]), 10)
  expect_equal(unname(fit$metrics[1L, "df"]), 2)
  expect_equal(unname(fit$metrics[1L, "aic"]), -16)
  expect_equal(unname(fit$metrics[1L, "df_resid"]), 18)
})
