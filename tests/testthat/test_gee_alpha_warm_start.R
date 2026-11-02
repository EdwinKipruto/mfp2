test_that("FP candidates pass successful GEE alpha to the next candidate", {
  starts <- list()
  powers <- matrix(c(-1, 0, 1, 2), ncol = 1L)

  testthat::local_mocked_bindings(
    transform_data_step = function(...) list(
      data_adj = NULL,
      data_fp = lapply(powers[, 1L], function(p) matrix(rep(p, 20), ncol = 1L)),
      fp_basis = NULL, acd_basis = NULL, powers_fp = powers,
      current_params = list(x = list())
    ),
    fit_model = function(..., gee_alpha_start) {
      starts <<- c(starts, list(gee_alpha_start))
      # A failed candidate must not replace the previous good alpha.
      list(gee_alpha = if (length(starts) == 2L) NULL else 0.1 * length(starts))
    },
    calculate_model_metrics = function(...) {
      c(logl = -1, df = 3, deviance_rs = -1,
        deviance_gaussian = NA_real_, aic = 8, bic = 10, df_resid = 17)
    },
    .package = "mfp2"
  )

  find_best_fpm_step(
    x = matrix(seq_len(20), ncol = 1L, dimnames = list(NULL, "x")),
    xi = "x", degree = 1, y = seq_len(20),
    powers_current = list(x = 1), powers = list(x = powers[, 1L]),
    acdx = c(x = FALSE), family = NULL, family_string = "gee",
    zero = c(x = FALSE), catzero = list(x = NULL), spike = c(x = FALSE),
    spike_decision = c(x = 2), acd_parameter = list(x = NULL),
    prev_adj_params = list(), has_offset = FALSE,
    precomputed_adj = list(), n_obs = 20,
    term_to_columns = list(x = 1L)
  )

  expect_length(starts, 4L)
  expect_null(starts[[1L]])
  expect_equal(starts[[2L]], 0.1)
  expect_equal(starts[[3L]], 0.1)
  expect_equal(starts[[4L]], 0.3)
})

test_that("GEE alpha warm start still estimates the candidate correlation", {
  skip_if_not_installed("geepack")
  set.seed(842)
  id <- rep(seq_len(28), each = 4L)
  x <- seq(0.2, 2, length.out = length(id))
  y <- 1 + 0.6 * x + stats::rnorm(28)[id] +
    stats::rnorm(length(id), sd = 0.25)
  family <- prepare_gee_family(
    gee_family(stats::gaussian(), corstr = "exchangeable"),
    y = y, id = id
  )
  design <- cbind("(Intercept)" = 1, x = x)

  cold <- fit_gee(x = design, y = y, family = family,
                  x_has_intercept = TRUE)
  warm <- fit_gee(x = design, y = y, family = family,
                  x_has_intercept = TRUE, alpha_start = cold$gee_alpha)
  invalid <- fit_gee(x = design, y = y, family = family,
                     x_has_intercept = TRUE, alpha_start = Inf)

  expect_length(cold$gee_alpha, 1L)
  expect_true(is.finite(cold$gee_alpha))
  expect_equal(warm$logl, cold$logl, tolerance = 1e-5)
  expect_equal(warm$gee_alpha, cold$gee_alpha, tolerance = 1e-4)
  expect_equal(invalid$logl, cold$logl, tolerance = 1e-8)
})

test_that("GEE information-criterion candidates skip unused Wald and jackknife work", {
  skip_if_not_installed("geepack")
  set.seed(843)
  id <- rep(seq_len(18), each = 3L)
  x <- seq(0.3, 1.7, length.out = length(id))
  y <- 2 + x + stats::rnorm(18)[id] + stats::rnorm(length(id), sd = 0.3)
  family <- prepare_gee_family(
    gee_family(stats::gaussian(), corstr = "exchangeable", std.err = "jack"),
    y = y, id = id
  )
  design <- cbind("(Intercept)" = 1, x = x)

  ic <- fit_gee(x = design, y = y, family = family,
                x_has_intercept = TRUE, selection_criterion = "aic")
  pvalue <- fit_gee(x = design, y = y, family = family,
                    x_has_intercept = TRUE, selection_criterion = "pvalue")

  expect_equal(ic$logl, pvalue$logl, tolerance = 1e-6)
  expect_true(is.na(ic$selection_deviance))
  expect_true(all(is.na(ic$robust_vcov)))
  expect_true(is.finite(pvalue$selection_deviance))
  expect_true(all(is.finite(pvalue$robust_vcov)))
  expect_true(is.na(pvalue$family_deviance))
})
