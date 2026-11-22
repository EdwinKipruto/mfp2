# Deterministic failures isolate selection policy from backend convergence.
spline_failure_args <- function(criterion = "aic", spike = FALSE,
                                 family_string = "gaussian", keep = character()) {
  columns <- paste0("x", seq_len(4L + as.integer(spike)))
  mapping <- list(x = columns)
  attr(mapping, "spline_specs") <- list(
    x = list(type = "rcs", df = 4L, indicator = spike, spike = spike)
  )
  list(
    x = matrix(seq_len(30L * length(columns)), nrow = 30L,
               dimnames = list(NULL, columns)),
    xi = "x", y = seq_len(30L), term_to_columns = mapping,
    criterion = criterion, ftest = FALSE, family_string = family_string,
    select = 1, alpha = 1, keep = keep, degree = 1,
    force_max_fp = c(x = FALSE), powers_current = list(x = 1),
    powers = list(x = c(0, 1)), acdx = c(x = FALSE), zero = c(x = FALSE),
    catzero = list(x = NULL), spike = c(x = FALSE),
    spike_decision = c(x = 2L), acd_parameter = list(x = NULL),
    prev_adj_params = list(x = NULL), transform_cache = NULL, n_obs = 30L
  )
}

spline_failure_fit <- function(width) {
  list(
    metrics = rbind(c(logl = -50 + 3 * width, df = width + 1,
                       deviance_rs = 100 - 6 * width,
                       deviance_gaussian = NA_real_,
                       aic = 102 - 4 * width, bic = 110 - 2 * width,
                       df_resid = 29 - width)),
    current_adj_params = list(x = list(width = width))
  )
}

spline_failure_mocks <- function(fail = character(), ordinary_error = FALSE,
                                null_failure = FALSE) {
  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) list(transform_cache = list()),
    fit_null_step = function(...) {
      if (null_failure) mfp2:::mfp2_stop_fit("required null failure")
      spline_failure_fit(0L)
    },
    fit_linear_step = function(...) {
      args <- list(...)
      columns <- args$term_to_columns$x
      indicator <- "x5" %in% columns
      df <- length(columns) - as.integer(indicator)
      key <- if (df == 0L) "binary" else
        paste0(df, if (indicator) "B" else "C")
      if (key %in% fail || "all" %in% fail) {
        if (ordinary_error) stop("unrelated programming error")
        mfp2:::mfp2_stop_fit(paste("injected failure", key))
      }
      spline_failure_fit(length(columns))
    },
    .package = "mfp2", .env = parent.frame()
  )
}

test_that("spline null and candidates receive the GEE selection criterion", {
  for (criterion in c("aic", "bic", "pvalue")) {
    seen <- character()
    testthat::local_mocked_bindings(
      build_adjustment_step = function(...) list(transform_cache = list()),
      fit_null_step = function(...) {
        seen <<- c(seen, list(...)$gee_selection_criterion)
        spline_failure_fit(0L)
      },
      fit_linear_step = function(...) {
        args <- list(...)
        seen <<- c(seen, args$gee_selection_criterion)
        spline_failure_fit(length(args$term_to_columns$x))
      }, .package = "mfp2"
    )
    result <- mfp2:::select_spline_step(
      spline_failure_args(criterion, family_string = "gee"), 4L
    )
    expect_identical(seen, rep(criterion, 5L))
    expect_identical(result$model_best, "4C")
  }
})

test_that("IC selection skips an optional failed spline and preserves scores", {
  for (criterion in c("aic", "bic")) {
    spline_failure_mocks(fail = "2C")
    expect_warning(
      result <- mfp2:::select_spline_step(spline_failure_args(criterion), 4L),
      "skipped \\(2C\\)", class = "mfp2_skipped_candidates"
    )
    expect_identical(result$candidate_keys, c("null", "4C", "1C", "3C"))
    expect_identical(result$model_best, "4C")
    expected <- do.call(rbind, lapply(c(0L, 4L, 1L, 3L), function(df) {
      spline_failure_fit(df)$metrics
    }))
    expect_equal(unname(result$metrics), unname(expected))
    expect_identical(result$current_adj_params$x$width, 4L)
  }
})

test_that("IC selection can select a smaller spline when the maximum fails", {
  spline_failure_mocks(fail = "4C")
  expect_warning(
    result <- mfp2:::select_spline_step(spline_failure_args("bic"), 4L),
    "skipped \\(4C\\)"
  )
  expect_identical(result$model_best, "3C")
  expect_identical(result$power_best, 3L)
})

test_that("IC selection drops a term when every optional form fails", {
  spline_failure_mocks(fail = "all")
  expect_warning(
    result <- mfp2:::select_spline_step(spline_failure_args(), 4L),
    "4 spline candidate"
  )
  expect_identical(result$model_best, "null")
  expect_true(is.na(result$power_best))
  expect_identical(result$candidate_keys, "null")
  expect_identical(result$current_adj_params$x$width, 0L)
})

test_that("a kept spline fails clearly if every eligible candidate fails", {
  spline_failure_mocks(fail = "all")
  expect_error(
    suppressWarnings(mfp2:::select_spline_step(
      spline_failure_args(keep = "x"), 4L
    )),
    "No finite AIC spline candidates.*keep constraint"
  )
})

test_that("spline null failures and unrelated errors remain strict", {
  spline_failure_mocks(null_failure = TRUE)
  expect_error(mfp2:::select_spline_step(spline_failure_args(), 4L),
                "required null failure", class = "mfp2_candidate_failure")
  spline_failure_mocks(fail = "2C", ordinary_error = TRUE)
  expect_error(mfp2:::select_spline_step(spline_failure_args(), 4L),
                "unrelated programming error")
})

test_that("closed testing never skips a failed reference or reduced fit", {
  for (failed in c("4C", "1C", "2C")) {
    spline_failure_mocks(fail = failed)
    expect_error(
      mfp2:::select_spline_step(spline_failure_args("pvalue"), 4L),
      paste("injected failure", failed), class = "mfp2_candidate_failure"
    )
  }
})

test_that("spike IC candidates retain the successful indicator-only form", {
  spline_failure_mocks(fail = c("4B", "1B", "2B", "3B", "4C", "1C", "2C", "3C"))
  expect_warning(
    result <- mfp2:::select_spline_step(spline_failure_args(spike = TRUE), 4L),
    "8 spline candidate"
  )
  expect_identical(result$candidate_keys, c("null", "binary"))
  expect_identical(result$model_best, "binary")
  expect_identical(result$power_best, 0L)
  expect_equal(result$spike_decision[["x"]],
                mfp2:::saz_decision_codes[["binary_only"]])
})

test_that("GEE spline IC fits avoid candidate Wald covariance work", {
  skip_if_not_installed("geepack")
  set.seed(624)
  n <- 180L
  d <- data.frame(id = rep(seq_len(45L), each = 4L),
                  x = runif(n, .2, 4), z = runif(n, .5, 2),
                  w = runif(n, .5, 1.5), off = runif(n, -.1, .1))
  d$y <- sin(2 * d$x) + .2 * d$z + d$off +
    rep(rnorm(45L, sd = .2), each = 4L) + rnorm(n, sd = .4)
  active <- FALSE
  original_step <- mfp2:::select_spline_step
  original_wald <- mfp2:::gee_selection_deviance
  original_vcov <- mfp2:::gee_selected_vcov
  testthat::local_mocked_bindings(
    select_spline_step = function(...) {
      active <<- TRUE
      on.exit({ active <<- FALSE }, add = TRUE)
      original_step(...)
    },
    gee_selection_deviance = function(...) {
      if (active) stop("unexpected candidate Wald calculation")
      original_wald(...)
    },
    gee_selected_vcov = function(...) {
      if (active) stop("unexpected candidate covariance calculation")
      original_vcov(...)
    }, .package = "mfp2"
  )
  for (criterion in c("aic", "bic")) {
    for (form in list(y ~ rcs(x, df = 3) + fp(z, df = 1) + offset(off),
                      y ~ bspline(x) + fp(z, df = 1) + offset(off))) {
      fit <- mfp2(form, data = d, weights = d$w, id = id,
                    family = gee_family(gaussian(), corstr = "exchangeable"),
                    criterion = criterion, xorder = "original", verbose = FALSE)
      expect_s3_class(fit, "mfp2")
      expect_true(all(is.finite(predict(fit))))
      expect_true(all(is.finite(vcov(fit))))
    }
  }
})
