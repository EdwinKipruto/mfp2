# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 15. force_max_fp_vars and formula-term force_max_fp
# =============================================================================

# Test purpose: Checks that force_max_fp_vars uses the requested maximum FP
# complexity for named predictors under AIC.
test_that("force_max_fp_vars forces maximum FP degree with AIC/BIC", {
  fit_force <- mfp2(
    x_prostate, y_prostate,
    criterion = "aic",
    force_max_fp_vars = colnames(x_prostate),
    select = 1,
    verbose = FALSE, warn_low_information = FALSE
  )

  for (v in get_selected_variable_names(fit_force)) {
    powers <- fit_force$fp_powers[[v]]
    requested_df <- as.numeric(fit_force$fp_terms[v, "df_setting"])

    expected_n_powers <- if (requested_df <= 1) {
      1
    } else {
      requested_df / 2
    }

    n_powers <- sum(!is.na(powers))

    expect_equal(
      n_powers,
      expected_n_powers,
      info = paste("Variable:", v)
    )
  }
})


# Test purpose: force_max_fp has one dedicated selector for p-value, AIC, and
# BIC selection. It fits only the predetermined maximum ordinary FP form;
# null, linear, and lower-degree FP models cannot affect a forced result.
test_that("force-max selector fits only the maximum ordinary FP form", {
  calls <- new.env(parent = emptyenv())
  calls$degrees <- integer(0)

  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        data_adj = NULL,
        current_params = list(),
        transform_cache = list(reused = TRUE)
      )
    },
    find_best_fpm_step = function(..., degree) {
      calls$degrees <- c(calls$degrees, degree)

      metrics <- matrix(
        c(
          -20, 4, 40, 40, 48, 50, 96,
          -18, 4, 36, 36, 44, 46, 96
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(NULL, metric_names)
      )

      list(
        powers = rbind(c(-1, -1), c(-1, 2)),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 2L,
        current_adj_params = list(x = list(data_adj = NULL))
      )
    },
    .package = "mfp2"
  )

  out <- select_force_max_fp(
    x = matrix(1:8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = 1:8,
    powers_current = list(x = c(1, 1)),
    powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
    criterion = "aic",
    ftest = FALSE,
    select = 0.05,
    alpha = 0.05,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(x = NULL),
    transform_cache = NULL,
    force_max_fp = c(x = TRUE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  expect_equal(calls$degrees, 2)
  expect_identical(out$model_best, 1L)
  expect_identical(rownames(out$metrics), "FP2")
  expect_identical(rownames(out$powers), "FP2")
  expect_equal(unname(out$power_best), c(-1, 2))
  expect_true(out$transform_cache$reused)
})


# Test purpose: The dedicated selector itself accepts p-value forcing. This
# verifies that p-value forcing no longer relies on select = 1 / alpha = 1 to
# walk the RA2 closed-test sequence before reaching the predetermined FPm.
test_that("force-max selector accepts p-value criterion directly", {
  calls <- new.env(parent = emptyenv())
  calls$degrees <- integer(0)

  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        data_adj = NULL,
        current_params = list(),
        transform_cache = list(pvalue = TRUE)
      )
    },
    find_best_fpm_step = function(..., degree) {
      calls$degrees <- c(calls$degrees, degree)
      metrics <- matrix(
        c(-18, 4, 36, 36, 44, 46, 96),
        nrow = 1,
        dimnames = list(NULL, metric_names)
      )
      list(
        powers = matrix(c(-1, 2), nrow = 1),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 1L,
        current_adj_params = list(x = list(data_adj = NULL))
      )
    },
    .package = "mfp2"
  )

  out <- select_force_max_fp(
    x = matrix(1:8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = 1:8,
    powers_current = list(x = c(1, 1)),
    powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
    criterion = "pvalue",
    ftest = FALSE,
    select = 1,
    alpha = 1,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(x = NULL),
    transform_cache = NULL,
    force_max_fp = c(x = TRUE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  expect_equal(calls$degrees, 2)
  expect_identical(out$model_best, 1L)
  expect_identical(rownames(out$metrics), "FP2")
  expect_equal(unname(out$power_best), c(-1, 2))
  expect_length(out$pvalue, 0L)
  expect_true(out$transform_cache$pvalue)
})


# Test purpose: The same criterion-independent selector handles forced ACD terms
# without fitting reduced ACD alternatives. The maximum ACD form is always
# FP1(x, A(x)), represented by the degree-2 ACD candidate search.
test_that("force-max selector fits only the full ACD form", {
  calls <- new.env(parent = emptyenv())
  calls$degrees <- integer(0)

  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        data_adj = NULL,
        current_params = list(),
        transform_cache = list(acd = TRUE)
      )
    },
    find_best_fpm_step = function(..., degree) {
      calls$degrees <- c(calls$degrees, degree)
      metrics <- matrix(
        c(-10, 4, 20, 20, 28, 30, 96),
        nrow = 1,
        dimnames = list(NULL, metric_names)
      )

      list(
        powers = matrix(c(-1, 2), nrow = 1),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 1L,
        current_adj_params = list(x = list(data_adj = NULL))
      )
    },
    .package = "mfp2"
  )

  out <- select_force_max_fp(
    x = matrix(1:8, ncol = 1, dimnames = list(NULL, "x")),
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = TRUE),
    y = 1:8,
    powers_current = list(x = c(1, 1)),
    powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
    criterion = "pvalue",
    ftest = FALSE,
    select = 1,
    alpha = 1,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = list(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = list()),
    prev_adj_params = list(x = NULL),
    transform_cache = NULL,
    force_max_fp = c(x = TRUE),
    has_offset = FALSE,
    n_obs = 8,
    term_to_columns = list(x = "x")
  )

  expect_equal(calls$degrees, 2)
  expect_true(out$acd)
  expect_identical(out$model_best, 1L)
  expect_identical(rownames(out$metrics), "FP1(x, A(x))")
  expect_identical(rownames(out$powers), "FP1(x, A(x))")
})


# Test purpose: find_best_fp_step() must dispatch every forced non-linear term
# to select_force_max_fp(), including criterion = "pvalue". The ordinary RA2
# and IC selectors must not run when the final functional complexity is forced.
test_that("force-max dispatch bypasses RA2 and IC selectors for all criteria", {
  calls <- new.env(parent = emptyenv())
  calls$criteria <- character(0)

  testthat::local_mocked_bindings(
    select_force_max_fp = function(..., criterion) {
      calls$criteria <- c(calls$criteria, criterion)
      metrics <- matrix(
        c(-10, 4, 20, 20, 28, 30, 96),
        nrow = 1,
        dimnames = list(
          "FP2",
          c("logl", "df", "deviance_rs", "deviance_gaussian",
            "aic", "bic", "df_resid")
        )
      )
      list(
        keep = FALSE,
        acd = FALSE,
        powers = matrix(c(-1, 2), nrow = 1, dimnames = list("FP2", NULL)),
        power_best = c(-1, 2),
        metrics = metrics,
        model_best = 1L,
        statistic = NA,
        pvalue = NA,
        spike = FALSE,
        current_adj_params = list(),
        transform_cache = list()
      )
    },
    select_ra2 = function(...) {
      stop("ordinary RA2 selector must not run", call. = FALSE)
    },
    select_ra2_acd = function(...) {
      stop("ACD RA2 selector must not run", call. = FALSE)
    },
    select_ic = function(...) {
      stop("ordinary IC selector must not run", call. = FALSE)
    },
    select_ic_acd = function(...) {
      stop("ACD IC selector must not run", call. = FALSE)
    },
    .package = "mfp2"
  )

  for (criterion_value in c("pvalue", "aic", "bic")) {
    out <- find_best_fp_step(
      x = matrix(seq_len(8), ncol = 1, dimnames = list(NULL, "x")),
      y = seq_len(8),
      xi = "x",
      weights = NULL,
      offset = NULL,
      df = 4,
      powers_current = list(x = c(1, 1)),
      family = stats::gaussian(),
      family_string = "gaussian",
      criterion = criterion_value,
      select = 1,
      alpha = 1,
      keep = character(0),
      powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
      method = NULL,
      strata = NULL,
      nocenter = FALSE,
      acdx = c(x = FALSE),
      ftest = FALSE,
      control = list(),
      rownames = as.character(seq_len(8)),
      zero = c(x = FALSE),
      catzero = list(x = NULL),
      spike = list(x = FALSE),
      spike_decision = c(x = 2),
      acd_parameter = list(x = NULL),
      prev_adj_params = list(x = NULL),
      transform_cache = NULL,
      force_max_fp = c(x = TRUE),
      has_offset = FALSE,
      n_obs = 8,
      verbose = FALSE,
      term_to_columns = list(x = "x")
    )

    expect_equal(unname(out$power_best), c(-1, 2))
  }

  expect_identical(calls$criteria, c("pvalue", "aic", "bic"))
})


# Test purpose: A forced eligible spike-at-zero term already has its final
# maximum representation after select_force_max_fp(): maximum continuous FP
# form plus binary zero indicator. SAZ Stage 2 must therefore be bypassed for
# p-value, AIC, and BIC so reduced component models cannot undo the force.
test_that("force-max spike terms bypass SAZ stage 2 for all criteria", {
  calls <- new.env(parent = emptyenv())
  calls$stage2 <- 0L
  calls$criteria <- character(0)

  testthat::local_mocked_bindings(
    select_force_max_fp = function(..., criterion) {
      calls$criteria <- c(calls$criteria, criterion)

      list(
        keep = FALSE,
        acd = FALSE,
        powers = matrix(c(-1, 2), nrow = 1,
                        dimnames = list("FP2 + Binary", NULL)),
        power_best = c(-1, 2),
        metrics = matrix(
          c(-10, 5, 20, 20, 30, 33, 95),
          nrow = 1,
          dimnames = list(
            "FP2 + Binary",
            c("logl", "df", "deviance_rs", "deviance_gaussian",
              "aic", "bic", "df_resid")
          )
        ),
        model_best = 1L,
        statistic = numeric(0),
        pvalue = numeric(0),
        spike = TRUE,
        current_adj_params = list(forced = TRUE),
        transform_cache = list(forced = TRUE)
      )
    },
    evaluate_saz_stage2 = function(...) {
      calls$stage2 <- calls$stage2 + 1L
      stop("SAZ stage 2 must not run for force_max_fp", call. = FALSE)
    },
    .package = "mfp2"
  )

  for (criterion_value in c("pvalue", "aic", "bic")) {
    out <- find_best_fp_step(
      x = matrix(seq_len(8), ncol = 1, dimnames = list(NULL, "x")),
      y = seq_len(8),
      xi = "x",
      weights = NULL,
      offset = NULL,
      df = 4,
      powers_current = list(x = c(1, 1)),
      family = stats::gaussian(),
      family_string = "gaussian",
      criterion = criterion_value,
      select = 1,
      alpha = 1,
      keep = character(0),
      powers = list(x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)),
      method = NULL,
      strata = NULL,
      nocenter = FALSE,
      acdx = c(x = FALSE),
      ftest = FALSE,
      control = list(),
      rownames = as.character(seq_len(8)),
      zero = c(x = TRUE),
      catzero = list(x = "x_binary"),
      spike = list(x = TRUE),
      spike_decision = c(x = saz_decision_codes[["binary_only"]]),
      acd_parameter = list(x = NULL),
      prev_adj_params = list(x = NULL),
      transform_cache = NULL,
      force_max_fp = c(x = TRUE),
      has_offset = FALSE,
      n_obs = 8,
      verbose = FALSE,
      term_to_columns = list(x = "x")
    )

    expect_equal(unname(out$power_best), c(-1, 2))
    expect_identical(
      unname(out$spike_decision[["x"]]),
      saz_decision_codes[["cont_binary"]]
    )
    expect_true(out$current_adj_params$forced)
    expect_true(out$transform_cache$forced)
  }

  expect_identical(calls$criteria, c("pvalue", "aic", "bic"))
  expect_identical(calls$stage2, 0L)
})


# Test purpose: Verifies that force_max_fp_vars translates to select = 1 and
# alpha = 1 under p-value selection in the matrix interface.
test_that("force_max_fp_vars forces maximum FP degree with p-value selection", {
  set.seed(1501)

  dat <- data.frame(
    x = seq(0.5, 6, length.out = 250)
  )
  dat$y <- 1 + 1.5 / dat$x - 0.7 * dat$x^2 +
    rnorm(nrow(dat), sd = 0.05)

  fit <- mfp2(
    x = as.matrix(dat["x"]),
    y = dat$y,
    criterion = "pvalue",
    df = 4,
    select = 0.05,
    alpha = 0.05,
    force_max_fp_vars = "x",
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["x", "select"]), 1)
  expect_equal(as.numeric(fit$fp_terms["x", "alpha"]), 1)
  expect_true(fit$fp_terms["x", "selected"])
  expect_equal(sum(!is.na(fit$fp_powers[["x"]])), 2)
})


# Test purpose: Verifies that fp(force_max_fp = TRUE) reaches the same
# p-value forcing logic through mfp2.formula().
test_that("formula force_max_fp uses p-value forcing from mfp2.default", {
  set.seed(1502)

  dat <- data.frame(
    x = seq(0.5, 6, length.out = 250)
  )
  dat$y <- 1 + 1.5 / dat$x - 0.7 * dat$x^2 +
    rnorm(nrow(dat), sd = 0.05)

  fit <- mfp2(
    y ~ fp(x, df = 4, force_max_fp = TRUE),
    data = dat,
    criterion = "pvalue",
    select = 0.05,
    alpha = 0.05,
    verbose = FALSE
  )

  expect_equal(as.numeric(fit$fp_terms["x", "select"]), 1)
  expect_equal(as.numeric(fit$fp_terms["x", "alpha"]), 1)
  expect_true(fit$fp_terms["x", "selected"])
  expect_equal(sum(!is.na(fit$fp_powers[["x"]])), 2)
})
