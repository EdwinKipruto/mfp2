# Weighted ordinal fitting must use the same likelihood throughout selection,
# final fitting, and MFPI comparisons.

# rms reports a limitation of its validation/bootstrap methods when weights
# are supplied. Allow only that exact notice in weighted fitting oracles;
# convergence and all other warnings must remain visible to testthat.
with_orm_weight_notice <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (identical(
        conditionMessage(w),
        "currently weights are ignored in model validation and bootstrapping orm fits"
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

make_ordinal_weight_data <- function(n = 300L, seed = 8041L) {
  set.seed(seed)
  dat <- data.frame(
    x = stats::runif(n, 0.5, 5),
    z = stats::rnorm(n),
    group = factor(rep(c("control", "treated"), length.out = n)),
    off = stats::runif(n, -0.3, 0.3),
    w = rep(c(1, 2, 3), length.out = n)
  )
  eta <- with(dat, 0.6 * log(x) - 0.3 * z + off +
                0.4 * (group == "treated") +
                0.25 * x * (group == "treated"))
  dat$y <- ordered(cut(
    eta + stats::rlogis(n), breaks = c(-Inf, 0.5, 1.8, Inf),
    labels = c("low", "middle", "high")
  ))
  dat
}

test_that("ordinal preparation validates and retains case weights", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  prepare <- function(w) {
    mfp2:::prepare_ordinal_family(ordinal_family(), dat$y, weights = w)
  }
  expect_equal(prepare(dat$w)$prepared$weights, dat$w)
  expect_equal(prepare(NULL)$prepared$weights, rep(1, nrow(dat)))
  expect_error(prepare(rep(0, nrow(dat))), "strictly positive")
  expect_error(prepare(rep(-1, nrow(dat))), "strictly positive")
  expect_error(prepare(rep(NA_real_, nrow(dat))), "missing")
  expect_error(prepare(rep(Inf, nrow(dat))), "finite")
  expect_error(prepare(dat$w[-1L]), "one value per observation")
})

test_that("weighted ordinal mfp2 matches native rms with offsets and subsets", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  dat$w <- dat$w / 2
  rows <- seq_len(nrow(dat)) %% 4L != 0L

  for (link in c("logistic", "probit")) {
    fit <- with_orm_weight_notice(mfp2(
      y ~ x + z + offset(off), data = dat,
      family = ordinal_family(link = link), weights = w, subset = rows,
      df = 1, select = 1, center = FALSE, shift = 0, scale = 1,
      xorder = "original", control = list(maxit = 100, eps = 1e-8),
      verbose = FALSE
    ))
    reference <- with_orm_weight_notice(rms::orm(
      y ~ x + z + offset(off), data = dat[rows, ], family = link,
      weights = w, normwt = FALSE, x = TRUE, y = TRUE,
      maxit = 100, eps = 1e-8
    ))
    null_reference <- with_orm_weight_notice(rms::orm(
      y ~ offset(off), data = dat[rows, ], family = link,
      weights = w, normwt = FALSE, maxit = 100, eps = 1e-8
    ))
    expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
                 tolerance = 1e-5)
    expect_equal(unname(stats::vcov(fit)), unname(stats::vcov(reference)),
                 tolerance = 1e-5)
    expect_equal(as.numeric(stats::logLik(fit)),
                 as.numeric(stats::logLik(reference)), tolerance = 1e-5)
    expect_equal(fit$mfp_logl, as.numeric(stats::logLik(reference)),
                 tolerance = 1e-5)
    expect_equal(fit$null_deviance, -2 * as.numeric(stats::logLik(null_reference)),
                 tolerance = 1e-5)
    expect_equal(stats::AIC(fit), stats::AIC(reference), tolerance = 1e-5)
    expect_equal(fit$mfp2_family$prepared$weights, dat$w[rows])
  }
})

test_that("integer ordinal weights reproduce row replication during FP selection", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  repeated <- dat[rep(seq_len(nrow(dat)), times = dat$w), ]
  fit_data <- function(d, weights = NULL) {
    mfp2(
      as.matrix(d[, c("x", "z")]), d$y, weights = weights,
      family = ordinal_family(), offset = d$off,
      df = c(x = 2, z = 1), powers = list(x = c(0, 1)),
      force_max_fp_vars = "x", keep = "z", criterion = "aic",
      center = FALSE, shift = 0, scale = 1, xorder = "original",
      control = list(maxit = 100, eps = 1e-8), verbose = FALSE
    )
  }
  weighted <- with_orm_weight_notice(fit_data(dat, dat$w))
  expanded <- fit_data(repeated)

  expect_equal(weighted$fp_powers, expanded$fp_powers)
  expect_equal(stats::coef(weighted), stats::coef(expanded), tolerance = 1e-5)
  expect_equal(stats::vcov(weighted), stats::vcov(expanded), tolerance = 1e-5)
  expect_equal(weighted$mfp_logl, expanded$mfp_logl, tolerance = 1e-5)
  expect_equal(stats::AIC(weighted), stats::AIC(expanded), tolerance = 1e-5)
  nd <- as.matrix(dat[1:8, c("x", "z")])
  expect_equal(
    predict(weighted, newdata = nd, newoffset = dat$off[1:8], type = "response"),
    predict(expanded, newdata = nd, newoffset = dat$off[1:8], type = "response"),
    tolerance = 1e-5
  )
})

test_that("omitted and unit ordinal weights give identical fits", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  fit_data <- function(weights = NULL) {
    mfp2(
      y ~ x + z, data = dat, family = ordinal_family(), weights = weights,
      df = 1, select = 1, center = FALSE, shift = 0, scale = 1,
      xorder = "original", verbose = FALSE
    )
  }
  expect_no_warning(omitted <- fit_data())
  expect_no_warning(unit <- fit_data(rep(1, nrow(dat))))
  expect_equal(stats::coef(unit), stats::coef(omitted), tolerance = 0)
  expect_equal(stats::vcov(unit), stats::vcov(omitted), tolerance = 0)
  expect_equal(unit$mfp_logl, omitted$mfp_logl, tolerance = 0)
})

test_that("weighted intercept-only ordinal fits retain their offset and likelihood", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  family <- mfp2:::prepare_ordinal_family(
    ordinal_family(), dat$y, weights = dat$w, offset = dat$off
  )
  fit <- with_orm_weight_notice(mfp2:::fit_ordinal(
    x = NULL, family = family,
    control = mfp2:::normalize_ordinal_control(list(maxit = 100, eps = 1e-8)),
    fast = FALSE, calculate_fit_statistics = TRUE
  ))
  reference <- with_orm_weight_notice(rms::orm(
    y ~ offset(off), data = dat, weights = w, normwt = FALSE,
    maxit = 100, eps = 1e-8
  ))
  expect_equal(unname(stats::coef(fit$fit)), unname(stats::coef(reference)),
               tolerance = 1e-5)
  expect_equal(fit$logl, as.numeric(stats::logLik(reference)), tolerance = 1e-5)
  expect_equal(fit$null_logl, fit$logl, tolerance = 1e-5)
  expect_equal(fit$weights, dat$w)
})

test_that("a predictor named weights cannot mask ordinal case weights", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  dat$weights <- dat$x
  fit <- with_orm_weight_notice(mfp2(
    y ~ weights + z, data = dat, family = ordinal_family(), weights = w,
    df = 1, select = 1, center = FALSE, shift = 0, scale = 1,
    xorder = "original", verbose = FALSE
  ))
  reference <- with_orm_weight_notice(rms::orm(
    y ~ weights + z, data = dat, weights = w, normwt = FALSE
  ))
  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)),
               tolerance = 1e-5)
  expect_equal(fit$mfp_logl, as.numeric(stats::logLik(reference)),
               tolerance = 1e-5)
})

test_that("weighted ordinal MFPI comparisons match native interaction models", {
  skip_if_not_installed("rms", minimum_version = "7.0-0")
  dat <- make_ordinal_weight_data()
  fit <- with_orm_weight_notice(mfpi(
    y ~ group + x + z + offset(off), data = dat,
    family = ordinal_family(), weights = w,
    group_var = "group", interaction_vars = "x",
    interaction_forms = c(x = "linear"), include_group_var = TRUE,
    df = 1, select = 1, p_interact = 1, center = FALSE,
    shift = 0, scale = 1, xorder = "original",
    control = list(maxit = 100, eps = 1e-8), verbose = FALSE
  ))
  main <- with_orm_weight_notice(rms::orm(
    y ~ group + x + z + offset(off), data = dat, weights = w,
    normwt = FALSE, maxit = 100, eps = 1e-8
  ))
  interaction <- with_orm_weight_notice(rms::orm(
    y ~ group * x + z + offset(off), data = dat, weights = w,
    normwt = FALSE, maxit = 100, eps = 1e-8
  ))
  main_logl <- as.numeric(stats::logLik(main))
  int_logl <- as.numeric(stats::logLik(interaction))
  expected_diff <- 2 * (int_logl - main_logl)
  metrics <- fit$all_model_metrics

  expect_equal(as.numeric(stats::logLik(fit$adjustment_model)), main_logl,
               tolerance = 1e-5)
  expect_equal(fit$all_interaction_models$x$logl, int_logl, tolerance = 1e-5)
  expect_equal(metrics$deviance_diff, expected_diff, tolerance = 1e-5)
  expect_equal(metrics$pvalue, stats::pchisq(expected_diff, 1, lower.tail = FALSE),
               tolerance = 1e-5)
  expect_equal(metrics$AIC_main_minus_int,
               stats::AIC(main) - stats::AIC(interaction), tolerance = 1e-5)
  # BIC retains the package convention of counting fitted rows.
  expect_equal(metrics$BIC_main_minus_int,
               expected_diff - log(nrow(dat)), tolerance = 1e-5)
})
