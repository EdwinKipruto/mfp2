# =============================================================================
# Cox strata equivalence tests: mfp2() versus survival::coxph()
# =============================================================================
#
# These regression tests isolate the public `strata =` interface.  MFP
# selection is deliberately disabled as a source of model differences by using
# df = 1 and select = 1, so the fitted mfp2 model is the same ordinary linear
# Cox model as the survival::coxph() reference model.
#
# The suite checks every documented strata representation:
#   * character vector
#   * factor vector
#   * integer vector
#   * numeric vector
#   * logical vector
#   * one-column matrix
#   * one-column data frame
#   * multi-column matrix
#   * multi-column data frame
#
# For every representation, equivalence is checked for fitted coefficients,
# covariance matrices, partial log-likelihoods, and predict() results.  The
# prediction checks include training and new-data linear predictors and risks,
# both zero and stratum-specific Cox references, plus absolute expected-event
# and survival predictions.
# =============================================================================

library(testthat)
library(survival)
library(mfp2)


# Generate a stable Cox problem with enough observations and events in every
# stratum.  The strata are deterministic/balanced so a test failure is not
# caused by an accidentally empty or event-free random stratum.
make_cox_strata_equivalence_data <- function(n = 360L, seed = 31001L) {
  stopifnot(n %% 12L == 0L)
  set.seed(seed)

  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1.5, 1.5),
    s_char = rep(c("A", "B", "C"), length.out = n),
    s_logical = rep(c(FALSE, TRUE), length.out = n),
    s_second = rep(c("X", "X", "Y", "Y"), length.out = n),
    stringsAsFactors = FALSE
  )

  eta <- 0.45 * dat$x1 - 0.30 * dat$x2
  base_mult <-
    c(A = 0.75, B = 1.00, C = 1.35)[dat$s_char] *
    ifelse(dat$s_logical, 1.20, 0.85) *
    c(X = 0.90, Y = 1.15)[dat$s_second]

  event_time <- stats::rexp(n, rate = 0.035 * base_mult * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.012)

  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)
  dat
}


# Compare one forced-linear mfp2 Cox fit with a direct coxph fit.  Names are
# intentionally ignored for coefficient/covariance comparisons because mfp2 may
# retain package-specific transformed-column labels even when df = 1.
expect_cox_strata_fit_equal <- function(fit_mfp2,
                                        fit_coxph,
                                        info,
                                        tolerance = 1e-8) {
  beta_mfp2 <- stats::coef(fit_mfp2)
  beta_coxph <- stats::coef(fit_coxph)
  vcov_mfp2 <- stats::vcov(fit_mfp2)
  vcov_coxph <- stats::vcov(fit_coxph)

  expect_equal(length(beta_mfp2), length(beta_coxph), info = info)
  expect_equal(
    unname(beta_mfp2),
    unname(beta_coxph),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(vcov_mfp2),
    unname(vcov_coxph),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    as.numeric(stats::logLik(fit_mfp2)),
    as.numeric(stats::logLik(fit_coxph)),
    tolerance = tolerance,
    info = info
  )
}


# Compare a predict() result that may be either a numeric vector or the
# list(fit, se.fit) form returned when se.fit = TRUE.
expect_prediction_equal <- function(got,
                                    expected,
                                    info,
                                    tolerance = 1e-8) {
  if (is.list(expected) && all(c("fit", "se.fit") %in% names(expected))) {
    expect_true(
      is.list(got) && all(c("fit", "se.fit") %in% names(got)),
      info = info
    )
    expect_equal(
      unname(got$fit),
      unname(expected$fit),
      tolerance = tolerance,
      info = info
    )
    expect_equal(
      unname(got$se.fit),
      unname(expected$se.fit),
      tolerance = tolerance,
      info = info
    )
  } else {
    expect_equal(
      unname(got),
      unname(expected),
      tolerance = tolerance,
      info = info
    )
  }
}


# Exercise the public Cox prediction paths that depend on correct stratum
# reconstruction.  Relative predictions are checked on both the zero-reference
# and stratum-reference scales.  Absolute predictions use the observed follow-up
# Surv response in newdata, matching predict.coxph(type = "expected"/
# "survival").
expect_cox_strata_predictions_equal <- function(fit_mfp2,
                                                fit_coxph,
                                                newx,
                                                new_strata,
                                                cox_newdata,
                                                new_surv,
                                                info,
                                                tolerance = 1e-8) {
  # Training-data predictions test the strata retained at fit time.
  for (reference in c("zero", "strata")) {
    got_lp <- predict(
      fit_mfp2,
      type = "lp",
      cox_reference = reference,
      se.fit = TRUE
    )
    expected_lp <- predict(
      fit_coxph,
      type = "lp",
      reference = reference,
      se.fit = TRUE
    )
    expect_prediction_equal(
      got_lp,
      expected_lp,
      info = paste(info, "training lp", reference),
      tolerance = tolerance
    )

    got_risk <- predict(
      fit_mfp2,
      type = "risk",
      cox_reference = reference,
      se.fit = TRUE
    )
    expected_risk <- predict(
      fit_coxph,
      type = "risk",
      reference = reference,
      se.fit = TRUE
    )
    expect_prediction_equal(
      got_risk,
      expected_risk,
      info = paste(info, "training risk", reference),
      tolerance = tolerance
    )
  }

  # New-data relative predictions test the explicit `strata =` prediction
  # argument in the same representation that was used at fit time.
  for (reference in c("zero", "strata")) {
    got_lp <- predict(
      fit_mfp2,
      newdata = newx,
      strata = new_strata,
      type = "lp",
      cox_reference = reference,
      se.fit = TRUE
    )
    expected_lp <- predict(
      fit_coxph,
      newdata = cox_newdata,
      type = "lp",
      reference = reference,
      se.fit = TRUE
    )
    expect_prediction_equal(
      got_lp,
      expected_lp,
      info = paste(info, "newdata lp", reference),
      tolerance = tolerance
    )

    got_risk <- predict(
      fit_mfp2,
      newdata = newx,
      strata = new_strata,
      type = "risk",
      cox_reference = reference,
      se.fit = TRUE
    )
    expected_risk <- predict(
      fit_coxph,
      newdata = cox_newdata,
      type = "risk",
      reference = reference,
      se.fit = TRUE
    )
    expect_prediction_equal(
      got_risk,
      expected_risk,
      info = paste(info, "newdata risk", reference),
      tolerance = tolerance
    )
  }

  # Matrix-interface absolute predictions carry their follow-up response as a
  # Surv column in newdata.  Strata remain a separate public predict() argument.
  absolute_newdata <- as.data.frame(newx, check.names = FALSE)
  absolute_newdata$prediction_response <- I(new_surv)

  for (prediction_type in c("expected", "survival")) {
    got <- predict(
      fit_mfp2,
      newdata = absolute_newdata,
      strata = new_strata,
      type = prediction_type,
      se.fit = TRUE
    )
    expected <- predict(
      fit_coxph,
      newdata = cox_newdata,
      type = prediction_type,
      se.fit = TRUE
    )
    expect_prediction_equal(
      got,
      expected,
      info = paste(info, "newdata", prediction_type),
      tolerance = tolerance
    )
  }

  # Absolute training predictions use the fitted response and therefore also
  # verify that the normalized strata stored by mfp2 agree with coxph.
  for (prediction_type in c("expected", "survival")) {
    expect_prediction_equal(
      predict(fit_mfp2, type = prediction_type, se.fit = TRUE),
      predict(fit_coxph, type = prediction_type, se.fit = TRUE),
      info = paste(info, "training", prediction_type),
      tolerance = tolerance
    )
  }
}


# Fit and compare one external-strata representation.  `reference_formula`
# describes the equivalent coxph stratification, while `fit_strata` and
# `new_strata` are passed through mfp2's public `strata =` arguments unchanged.
run_external_strata_equivalence_case <- function(dat,
                                                 rows,
                                                 fit_strata,
                                                 new_strata,
                                                 reference_formula,
                                                 reference_data,
                                                 info) {
  x <- as.matrix(dat[, c("x1", "x2"), drop = FALSE])
  y <- survival::Surv(dat$time, dat$status)

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "cox",
    strata = fit_strata,
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  fit_coxph <- survival::coxph(
    formula = reference_formula,
    data = reference_data,
    ties = "breslow",
    nocenter = NULL,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  # select = 1 and df = 1 are intended to leave both ordinary predictors in a
  # purely linear model.  Assert that premise before comparing against coxph.
  expect_true(
    setequal(get_selected_variable_names(fit_mfp2), c("x1", "x2")),
    info = paste(info, "selected variables")
  )

  expect_cox_strata_fit_equal(
    fit_mfp2,
    fit_coxph,
    info = paste(info, "fit")
  )

  newx <- x[rows, , drop = FALSE]
  cox_newdata <- reference_data[rows, , drop = FALSE]
  new_surv <- y[rows]

  expect_cox_strata_predictions_equal(
    fit_mfp2 = fit_mfp2,
    fit_coxph = fit_coxph,
    newx = newx,
    new_strata = new_strata,
    cox_newdata = cox_newdata,
    new_surv = new_surv,
    info = info
  )
}


# -----------------------------------------------------------------------------
# One-dimensional external strata
# -----------------------------------------------------------------------------

test_that("mfp2 external vector strata match coxph fits and predictions", {
  dat <- make_cox_strata_equivalence_data()
  rows <- c(2L, 7L, 14L, 25L, 38L, 61L, 92L, 137L, 214L, 301L)

  vector_cases <- list(
    character = list(
      fit = dat$s_char,
      new = dat$s_char[rows],
      value = dat$s_char
    ),
    factor = list(
      fit = factor(dat$s_char, levels = c("C", "A", "B")),
      new = factor(dat$s_char[rows], levels = c("C", "A", "B")),
      value = factor(dat$s_char, levels = c("C", "A", "B"))
    ),
    integer = list(
      fit = as.integer(match(dat$s_char, c("B", "C", "A"))),
      new = as.integer(match(dat$s_char[rows], c("B", "C", "A"))),
      value = as.integer(match(dat$s_char, c("B", "C", "A")))
    ),
    numeric = list(
      fit = as.numeric(c(A = 10.5, B = 20.5, C = 40.5)[dat$s_char]),
      new = as.numeric(c(A = 10.5, B = 20.5, C = 40.5)[dat$s_char[rows]]),
      value = as.numeric(c(A = 10.5, B = 20.5, C = 40.5)[dat$s_char])
    ),
    logical = list(
      fit = dat$s_logical,
      new = dat$s_logical[rows],
      value = dat$s_logical
    )
  )

  for (case_name in names(vector_cases)) {
    case <- vector_cases[[case_name]]
    reference_data <- dat
    reference_data$stratum <- case$value

    run_external_strata_equivalence_case(
      dat = dat,
      rows = rows,
      fit_strata = case$fit,
      new_strata = case$new,
      reference_formula = survival::Surv(time, status) ~
        x1 + x2 + strata(stratum),
      reference_data = reference_data,
      info = paste("vector strata:", case_name)
    )
  }
})


# -----------------------------------------------------------------------------
# Matrix and data-frame external strata
# -----------------------------------------------------------------------------

test_that("mfp2 tabular strata match coxph fits and predictions", {
  dat <- make_cox_strata_equivalence_data(seed = 31002L)
  rows <- c(1L, 9L, 18L, 33L, 52L, 79L, 111L, 166L, 243L, 337L)

  # A one-column matrix/data frame should have the same semantics as one strata
  # variable.  Multi-column inputs should have the same semantics as
  # coxph(... + strata(s1, s2)).
  one_matrix <- matrix(dat$s_char, ncol = 1L, dimnames = list(NULL, "s1"))
  one_frame <- data.frame(s1 = dat$s_char, stringsAsFactors = FALSE)

  multi_matrix <- cbind(
    s1 = dat$s_char,
    s2 = ifelse(dat$s_logical, "yes", "no")
  )
  multi_frame <- data.frame(
    s1 = factor(dat$s_char, levels = c("B", "A", "C")),
    s2 = dat$s_logical
  )

  cases <- list(
    one_column_matrix = list(
      fit = one_matrix,
      new = one_matrix[rows, , drop = FALSE],
      reference_data = transform(dat, s1 = dat$s_char),
      formula = survival::Surv(time, status) ~
        x1 + x2 + strata(s1)
    ),
    one_column_data_frame = list(
      fit = one_frame,
      new = one_frame[rows, , drop = FALSE],
      reference_data = transform(dat, s1 = dat$s_char),
      formula = survival::Surv(time, status) ~
        x1 + x2 + strata(s1)
    ),
    multi_column_matrix = list(
      fit = multi_matrix,
      new = multi_matrix[rows, , drop = FALSE],
      reference_data = transform(
        dat,
        s1 = dat$s_char,
        s2 = ifelse(dat$s_logical, "yes", "no")
      ),
      formula = survival::Surv(time, status) ~
        x1 + x2 + strata(s1, s2)
    ),
    multi_column_data_frame = list(
      fit = multi_frame,
      new = multi_frame[rows, , drop = FALSE],
      reference_data = transform(
        dat,
        s1 = factor(dat$s_char, levels = c("B", "A", "C")),
        s2 = dat$s_logical
      ),
      formula = survival::Surv(time, status) ~
        x1 + x2 + strata(s1, s2)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]

    run_external_strata_equivalence_case(
      dat = dat,
      rows = rows,
      fit_strata = case$fit,
      new_strata = case$new,
      reference_formula = case$formula,
      reference_data = case$reference_data,
      info = paste("tabular strata:", case_name)
    )
  }
})


# -----------------------------------------------------------------------------
# Prediction-level reconstruction regression checks
# -----------------------------------------------------------------------------

test_that("mfp2 prediction strata use fitted levels and reject unseen strata", {
  dat <- make_cox_strata_equivalence_data(seed = 31004L)
  x <- as.matrix(dat[, c("x1", "x2"), drop = FALSE])
  y <- survival::Surv(dat$time, dat$status)

  fit <- mfp2(
    x = x,
    y = y,
    family = "cox",
    strata = as.integer(match(dat$s_char, c("A", "B", "C"))),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  rows <- which(dat$s_char == "A")[1:8]
  new_strata <- rep(1L, length(rows))
  prepared <- prepare_newdata_for_predict(
    fit,
    newdata = x[rows, , drop = FALSE],
    strata = new_strata,
    check_binary = FALSE
  )

  fitted_levels <- fit$xlevels[["strata(strata_)"]]
  expect_true(is.factor(prepared$strata_))
  expect_identical(levels(prepared$strata_), fitted_levels)
  expect_identical(as.character(prepared$strata_), rep("1", length(rows)))

  # The source-level validation should fail before survival::model.frame() can
  # emit its less informative "factor ... has new levels" error.
  expect_error(
    predict(
      fit,
      newdata = x[rows, , drop = FALSE],
      strata = rep(99L, length(rows)),
      type = "lp"
    ),
    "not present in the fitted Cox model"
  )
})


# -----------------------------------------------------------------------------
# Formula strata() regression checks
# -----------------------------------------------------------------------------
# The external-strata tests above target the normalization added to mfp2's
# `strata =` argument.  These two checks ensure that the existing formula path
# remains identical to coxph for both one and multiple strata() terms, including
# all native Cox prediction scales.

test_that("mfp2 formula strata terms match coxph fits and predictions", {
  dat <- make_cox_strata_equivalence_data(seed = 31003L)
  rows <- c(3L, 11L, 26L, 47L, 73L, 104L, 158L, 219L, 276L, 352L)

  formula_cases <- list(
    single = survival::Surv(time, status) ~
      x1 + x2 + strata(s_char),
    multiple = survival::Surv(time, status) ~
      x1 + x2 + strata(s_char) + strata(s_second)
  )

  for (case_name in names(formula_cases)) {
    model_formula <- formula_cases[[case_name]]

    fit_mfp2 <- mfp2(
      model_formula,
      data = dat,
      family = "cox",
      df = 1,
      select = 1,
      alpha = 1,
      cycles = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      nocenter = NULL,
      xorder = "original",
      ties = "breslow",
      verbose = FALSE
    )

    fit_coxph <- survival::coxph(
      model_formula,
      data = dat,
      ties = "breslow",
      nocenter = NULL,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )

    info <- paste("formula strata:", case_name)
    expect_cox_strata_fit_equal(fit_mfp2, fit_coxph, info = info)

    nd <- dat[rows, , drop = FALSE]

    for (reference in c("zero", "strata")) {
      expect_prediction_equal(
        predict(
          fit_mfp2,
          newdata = nd,
          type = "lp",
          cox_reference = reference,
          se.fit = TRUE
        ),
        predict(
          fit_coxph,
          newdata = nd,
          type = "lp",
          reference = reference,
          se.fit = TRUE
        ),
        info = paste(info, "newdata lp", reference)
      )
      expect_prediction_equal(
        predict(
          fit_mfp2,
          newdata = nd,
          type = "risk",
          cox_reference = reference,
          se.fit = TRUE
        ),
        predict(
          fit_coxph,
          newdata = nd,
          type = "risk",
          reference = reference,
          se.fit = TRUE
        ),
        info = paste(info, "newdata risk", reference)
      )
    }

    for (prediction_type in c("expected", "survival")) {
      expect_prediction_equal(
        predict(
          fit_mfp2,
          newdata = nd,
          type = prediction_type,
          se.fit = TRUE
        ),
        predict(
          fit_coxph,
          newdata = nd,
          type = prediction_type,
          se.fit = TRUE
        ),
        info = paste(info, "newdata", prediction_type)
      )
    }
  }
})


# Migrated coverage from the former test_mfp2.R

# Exact partial likelihood is deliberately rejected because the low-level Cox
# candidate fitter used during MFP/MFPI selection does not implement it. These
# tests cover both matrix and formula entry points so the unsupported method
# cannot bypass public argument validation.
test_that("exact Cox ties are rejected before MFP-based selection", {
  dat <- data.frame(
    time = c(1, 2, 3, 4, 5, 6),
    status = c(1, 1, 0, 1, 0, 1),
    x = c(0.2, 0.5, 0.1, 0.8, 0.3, 0.7),
    group = c(0, 1, 0, 1, 0, 1)
  )
  y <- survival::Surv(dat$time, dat$status)
  x <- as.matrix(dat[c("x", "group")])
  msg <- "not supported for MFP selection"

  expect_error(
    mfp2(x, y, family = "cox", ties = "exact", verbose = FALSE),
    msg
  )
  expect_error(
    mfp2(
      survival::Surv(time, status) ~ x + group,
      data = dat,
      family = "cox",
      ties = "exact",
      verbose = FALSE
    ),
    msg
  )
  expect_error(
    mfpi(
      x, y,
      group_var = "group",
      cont_vars = "x",
      family = "cox",
      ties = "exact",
      verbose = FALSE
    ),
    msg
  )
  expect_error(
    mfpi(
      survival::Surv(time, status) ~ x + group,
      data = dat,
      group_var = "group",
      cont_vars = "x",
      family = "cox",
      ties = "exact",
      verbose = FALSE
    ),
    msg
  )
})


# Test purpose: 8.2.5 A formula offset without strata should be included in
# training and newdata Cox linear predictors exactly as in coxph().
test_that("8.2.5 Cox formula offset matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8205)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 + offset(off),
    data = dat,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "off"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})


# Test purpose: 8.2.6 A retained categorical predictor without strata or
# offset should match coxph() completely.
test_that("8.2.6 Cox categorical predictor matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8206)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group,
    data = dat,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group,
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})


# Test purpose: 8.2.7 A categorical predictor combined with formula strata
# should preserve both grouped contrasts and stratum-specific risk sets.
test_that("8.2.7 Cox categorical predictor with strata matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8207)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group + strata(stratum1),
    data = dat,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group + strata(stratum1),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group", "stratum1"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})


# Test purpose: 8.2.8 A categorical predictor combined with a formula offset
# should match coxph() in all fitted and predicted quantities.
test_that("8.2.8 Cox categorical predictor with offset matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8208)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group + offset(off),
    data = dat,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group", "off"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})


# Test purpose: 8.2.9 Categorical contrasts, one formula strata term, and a
# formula offset should all be reconstructed together without changing the
# coxph() fit or predictions.
test_that("8.2.9 Cox categorical predictor with strata and offset matches coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8209)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + group +
      strata(stratum1) + offset(off),
    data = dat,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + group +
      strata(stratum1) + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "group", "stratum1", "off"), drop = FALSE]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})


# Test purpose: 8.2.10 Two formula strata terms and an offset should match
# coxph() simultaneously, including the offset contribution to LP and risk.
test_that("8.2.10 Cox multiple strata with offset match coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8210)

  fit_mfp2 <- mfp2(
    survival::Surv(time, status) ~ x1 + x2 +
      strata(stratum1) + strata(stratum2) + offset(off),
    data = dat,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 +
      strata(stratum1) + strata(stratum2) + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[
    1:30,
    c("x1", "x2", "stratum1", "stratum2", "off"),
    drop = FALSE
  ]

  expect_mfp2_cox_predictions_equal(fit_mfp2, fit_coxph, nd)
})


# Test purpose: 8.2.11 Matrix-interface strata and offset vectors should match
# the equivalent coxph() formula model for fitting and newdata prediction.
test_that("8.2.11 Cox matrix-interface strata and offset match coxph completely", {
  dat <- make_cox_equivalence_data(seed = 8211)
  x <- as.matrix(dat[, c("x1", "x2")])
  y <- survival::Surv(dat$time, dat$status)

  fit_mfp2 <- mfp2(
    x = x,
    y = y,
    family = "cox",
    strata = dat$stratum1,
    offset = dat$off,
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
  fit_coxph <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + strata(stratum1) + offset(off),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  nd <- dat[1:30, c("x1", "x2", "stratum1", "off"), drop = FALSE]
  newx <- as.matrix(nd[, c("x1", "x2")])

  expect_mfp2_cox_predictions_equal(
    fit_mfp2,
    fit_coxph,
    newdata = nd,
    mfp2_newdata = newx,
    mfp2_predict_args = list(
      strata = nd$stratum1,
      newoffset = nd$off
    )
  )
})


# Test purpose: External Cox strata supplied as character, numeric, factor, or
# multiple columns must be normalized to one factor before coxph.fit(). This
# prevents character labels from becoming NA through low-level numeric coercion.
test_that("Cox strata inputs are normalized before candidate fitting", {
  set.seed(30031)
  n <- 150L
  x1 <- runif(n, 1, 5)
  x2 <- rnorm(n)
  s1 <- rep(c("A", "B", "C"), length.out = n)
  s2 <- rep(c("X", "Y"), each = ceiling(n / 2L))[seq_len(n)]
  lp <- 0.25 * x1 - 0.15 * x2
  event_time <- rexp(n, rate = 0.03 * exp(lp))
  censor_time <- rexp(n, rate = 0.02)
  y <- survival::Surv(
    pmin(event_time, censor_time),
    as.integer(event_time <= censor_time)
  )
  x <- cbind(x1 = x1, x2 = x2)

  common <- list(
    x = x, y = y, family = "cox", keep = c("x1", "x2"),
    df = 1, select = 1, alpha = 1, cycles = 5, shift = 0, scale = 1,
    center = FALSE, xorder = "original", ties = "breslow", verbose = FALSE
  )

  fit_character <- do.call(mfp2, c(common, list(strata = s1)))
  fit_factor <- do.call(mfp2, c(common, list(strata = factor(s1))))
  fit_numeric <- do.call(mfp2, c(common, list(strata = match(s1, c("A", "B", "C")))))

  expect_equal(unname(coef(fit_character)), unname(coef(fit_factor)), tolerance = 1e-8)
  expect_equal(unname(coef(fit_numeric)), unname(coef(fit_factor)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit_character)), as.numeric(logLik(fit_factor)), tolerance = 1e-8)

  multi <- data.frame(site = s1, period = s2)
  fit_multi <- do.call(mfp2, c(common, list(strata = multi)))
  combined <- do.call(
    survival::strata,
    c(as.list(multi), list(shortlabel = TRUE))
  )
  fit_combined <- do.call(mfp2, c(common, list(strata = combined)))

  expect_equal(unname(coef(fit_multi)), unname(coef(fit_combined)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit_multi)), as.numeric(logLik(fit_combined)), tolerance = 1e-8)
})


# Test purpose: MFPI uses the same Cox-strata normalization as mfp2, including
# raw character labels supplied through the external strata argument.
test_that("MFPI accepts character Cox strata through the external argument", {
  set.seed(30032)
  n <- 120L
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  z <- rnorm(n)
  strata_char <- rep(c("S1", "S2", "S3"), length.out = n)
  lp <- 0.2 * x - 0.1 * z + 0.2 * (trt == "treated")
  event_time <- rexp(n, rate = 0.03 * exp(lp))
  censor_time <- rexp(n, rate = 0.02)

  fit <- mfpi(
    x = data.frame(trt = trt, x = x, z = z),
    y = survival::Surv(
      pmin(event_time, censor_time),
      as.integer(event_time <= censor_time)
    ),
    family = "cox", group_var = "trt", cont_vars = "x",
    cont_var_forms = c(x = "linear"), keep = "z", strata = strata_char,
    df = 1, select = 1, alpha = 1, cycles = 1, shift = 0, scale = 1,
    center = FALSE, p_interact = 1, ties = "breslow", verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
})


# Test purpose: Unsupported or malformed Cox strata fail before reaching
# survival::coxph.fit(), where type/shape checking is intentionally minimal.
test_that("Cox strata normalization rejects malformed inputs", {
  expect_error(
    normalize_cox_strata(as.list(rep(c("A", "B"), 5)), 10L),
    "vector, factor, matrix, or data frame"
  )

  expect_error(
    normalize_cox_strata(array(1:24, dim = c(4, 3, 2)), 4L),
    "vector, factor, matrix, or data frame"
  )

  empty_strata <- data.frame(row.names = seq_len(4L))
  expect_error(
    normalize_cox_strata(empty_strata, 4L),
    "at least one stratification variable"
  )
})
