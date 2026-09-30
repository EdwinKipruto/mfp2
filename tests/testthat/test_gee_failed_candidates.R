test_that("GEE FP search excludes failed powers without changing successful scores", {
  powers <- matrix(c(-1, 0, 1), ncol = 1L)
  fitted <- 0L
  scored <- 0L

  testthat::local_mocked_bindings(
    transform_data_step = function(...) list(
      data_adj = NULL,
      data_fp = lapply(powers[, 1L], function(p) matrix(rep(p, 20L), ncol = 1L)),
      fp_basis = NULL, acd_basis = NULL, powers_fp = powers,
      current_params = list(x = list())
    ),
    fit_model = function(x, gee_allow_failed_candidate, ...) {
      expect_true(gee_allow_failed_candidate)
      fitted <<- fitted + 1L
      list(converged = fitted != 1L, logl = -fitted,
           selection_deviance = -fitted)
    },
    calculate_model_metrics = function(obj, ...) {
      expect_true(obj$converged)
      scored <<- scored + 1L
      c(logl = obj$logl, df = 3, deviance_rs = obj$selection_deviance,
        deviance_gaussian = NA_real_, aic = 8, bic = 10, df_resid = 17)
    },
    .package = "mfp2"
  )

  result <- find_best_fpm_step(
    x = matrix(seq_len(20L), ncol = 1L, dimnames = list(NULL, "x")),
    xi = "x", degree = 1, y = seq_len(20L),
    powers_current = list(x = 1), powers = list(x = powers[, 1L]),
    acdx = c(x = FALSE), family = NULL, family_string = "gee",
    zero = c(x = FALSE), catzero = list(x = NULL), spike = c(x = FALSE),
    spike_decision = c(x = 2), acd_parameter = list(x = NULL),
    prev_adj_params = list(), has_offset = FALSE,
    precomputed_adj = list(), n_obs = 20,
    term_to_columns = list(x = 1L), criterion = "pvalue"
  )

  expect_equal(fitted, 3L)
  expect_equal(scored, 2L)
  expect_true(all(is.na(result$metrics[1L, ])))
  expect_equal(result$model_best, 3L)
  expect_equal(unname(result$power_best), 1)
})

test_that("GEE FP search reports an error if every power fails", {
  powers <- matrix(c(-1, 0), ncol = 1L)

  testthat::local_mocked_bindings(
    transform_data_step = function(...) list(
      data_adj = NULL,
      data_fp = lapply(powers[, 1L], function(p) matrix(rep(p, 20L), ncol = 1L)),
      fp_basis = NULL, acd_basis = NULL, powers_fp = powers,
      current_params = list(x = list())
    ),
    fit_model = function(...) list(converged = FALSE),
    calculate_model_metrics = function(...) stop("Failed fit must not be scored"),
    .package = "mfp2"
  )

  expect_error(
    find_best_fpm_step(
      x = matrix(seq_len(20L), ncol = 1L, dimnames = list(NULL, "x")),
      xi = "x", degree = 1, y = seq_len(20L),
      powers_current = list(x = 1), powers = list(x = powers[, 1L]),
      acdx = c(x = FALSE), family = NULL, family_string = "gee",
      zero = c(x = FALSE), catzero = list(x = NULL), spike = c(x = FALSE),
      spike_decision = c(x = 2), acd_parameter = list(x = NULL),
      prev_adj_params = list(), has_offset = FALSE,
      precomputed_adj = list(), n_obs = 20,
      term_to_columns = list(x = 1L), criterion = "pvalue"
    ),
    "No converged gee FP candidates"
  )
})

test_that("only FP power searches may return a failed GEE fit", {
  skip_if_not_installed("geepack")
  id <- rep(seq_len(12L), each = 3L)
  x <- seq(0.2, 1.8, length.out = length(id))
  y <- 2 + x + rep(seq_len(12L) / 12, each = 3L)
  family <- prepare_gee_family(
    gee_family(stats::gaussian(), corstr = "independence"),
    y = y, id = id
  )
  design <- cbind("(Intercept)" = 1, x = x)

  # Force an invalid score after an otherwise valid backend fit. This is
  # deterministic and tests the same eligibility decision as non-convergence.
  testthat::local_mocked_bindings(
    gee_quasi_likelihood = function(...) NA_real_,
    .package = "mfp2"
  )

  skipped <- fit_gee(
    x = design, y = y, family = family,
    x_has_intercept = TRUE, allow_failed_candidate = TRUE,
    selection_criterion = "pvalue"
  )
  expect_false(skipped$converged)
  expect_true(is.na(skipped$logl))
  expect_true(is.na(skipped$selection_deviance))

  expect_error(
    fit_gee(x = design, y = y, family = family,
            x_has_intercept = TRUE, selection_criterion = "pvalue"),
    "did not converge"
  )
})
