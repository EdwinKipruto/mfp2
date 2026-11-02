# Tests for Cox nocenter semantics and defaults.
#
# mfp2/mfpi deliberately mirror survival::coxph(): nocenter is a set of
# design-matrix values, not predictor indices. Columns whose observed values
# are all contained in c(-1, 0, 1) are left uncentered by default.

test_that("Cox nocenter defaults match coxph", {
  expected <- c(-1, 0, 1)

  expect_identical(eval(formals(mfp2.default)$nocenter), expected)
  expect_identical(eval(formals(mfp2.formula)$nocenter), expected)
  expect_identical(eval(formals(mfpi.default)$nocenter), expected)
  expect_identical(eval(formals(mfpi.formula)$nocenter), expected)
  expect_identical(eval(formals(fit_model)$nocenter), expected)
})

make_nocenter_cox_data <- function(n = 320L, seed = 15101L) {
  set.seed(seed)

  dat <- data.frame(
    x_cont = stats::rnorm(n),
    x_bin = stats::rbinom(n, 1L, 0.45),
    x_signed = sample(c(-1, 0, 1), n, replace = TRUE)
  )

  eta <- 0.35 * dat$x_cont + 0.45 * dat$x_bin - 0.25 * dat$x_signed
  event_time <- stats::rexp(n, rate = 0.08 * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.04)

  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  dat
}

fit_nocenter_mfp2 <- function(dat, nocenter_missing = TRUE) {
  args <- list(
    x = as.matrix(dat[c("x_cont", "x_bin", "x_signed")]),
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  if (!nocenter_missing) {
    # Assign through single-bracket indexing so NULL is retained as an explicit
    # argument instead of deleting the list element and falling back to the default.
    args["nocenter"] <- list(NULL)
  }

  do.call(mfp2, args)
}

test_that("default Cox nocenter matches coxph fit and prediction semantics", {
  dat <- make_nocenter_cox_data()

  fit_mfp2 <- fit_nocenter_mfp2(dat, nocenter_missing = TRUE)
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x_cont + x_bin + x_signed,
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )

  expect_equal(
    unname(stats::coef(fit_mfp2)),
    unname(stats::coef(fit_coxph)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_coxph)),
    tolerance = 1e-8
  )
  expect_equal(unname(fit_mfp2$means), unname(fit_coxph$means), tolerance = 1e-12)

  # Under coxph's default nocenter set, indicator-like columns are not
  # recentered, so their stored means are exactly zero. mfp2 may retain
  # transformed design-column names, so identify these columns from the
  # reference coxph fit and compare the corresponding positions.
  indicator_idx <- match(c("x_bin", "x_signed"), names(fit_coxph$means))
  expect_false(anyNA(indicator_idx))
  expect_equal(unname(fit_mfp2$means[indicator_idx]), c(0, 0))

  nd <- dat[1:25, c("x_cont", "x_bin", "x_signed"), drop = FALSE]
  got <- predict(
    fit_mfp2,
    newdata = nd,
    type = "lp",
    cox_reference = "zero"
  )
  expected <- predict(
    fit_coxph,
    newdata = nd,
    type = "lp",
    reference = "zero"
  )

  expect_equal(unname(got), unname(expected), tolerance = 1e-8)
})

test_that("explicit nocenter NULL retains the all-columns-centering option", {
  dat <- make_nocenter_cox_data(seed = 15102L)

  fit_mfp2 <- fit_nocenter_mfp2(dat, nocenter_missing = FALSE)
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x_cont + x_bin + x_signed,
    data = dat,
    ties = "breslow",
    nocenter = NULL,
    x = TRUE,
    y = TRUE
  )

  expect_equal(
    unname(stats::coef(fit_mfp2)),
    unname(stats::coef(fit_coxph)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_coxph)),
    tolerance = 1e-8
  )
  expect_equal(unname(fit_mfp2$means), unname(fit_coxph$means), tolerance = 1e-12)

  # NULL has a distinct public meaning: no value-set exemption is applied,
  # so the 0/1 and -1/0/1 columns are allowed to be recentered. As above,
  # use the reference coxph design positions rather than mfp2's transformed
  # column names.
  indicator_idx <- match(c("x_bin", "x_signed"), names(fit_coxph$means))
  expect_false(anyNA(indicator_idx))
  expect_equal(
    unname(fit_mfp2$means[indicator_idx]),
    c(mean(dat$x_bin), mean(dat$x_signed)),
    tolerance = 1e-12
  )
})
