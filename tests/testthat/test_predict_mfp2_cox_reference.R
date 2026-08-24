# Helpers local to the version-1 Cox reference regression suite.

make_cox_reference_v1_data <- function(n = 480L, seed = 11001L) {
  set.seed(seed)

  stratum <- factor(rep(c("A", "B"), length.out = n))
  x_mean <- ifelse(stratum == "A", 4, 9)
  x <- stats::rnorm(n, mean = x_mean, sd = 1.15)
  off <- stats::rnorm(n, mean = 0.35, sd = 0.12)

  eta <- 0.48 * (x - 6) + 0.30 * (stratum == "B") + off
  event_time <- stats::rexp(n, rate = 0.025 * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.012)

  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    x = x,
    stratum = stratum,
    off = off
  )
}


strip_mfp2_class_v1 <- function(object) {
  class(object) <- setdiff(class(object), "mfp2")
  object
}


prepare_formula_cox_newdata_v1 <- function(object,
                                           newdata,
                                           include_response = FALSE) {
  selected_terms <- get_selected_variable_names(object)

  predictor_data <- reconstruct_formula_newdata(
    object,
    newdata,
    terms = selected_terms
  )

  prediction_strata <- reconstruct_formula_strata_newdata(object, newdata)
  prediction_offset <- reconstruct_formula_offset_newdata(object, newdata)

  prepared <- prepare_newdata_for_predict(
    object,
    predictor_data,
    terms = selected_terms,
    strata = prediction_strata,
    offset = prediction_offset,
    check_binary = FALSE
  )

  if (include_response) {
    response_name <- cox_internal_response_name(object)
    prepared[[response_name]] <- I(
      survival::Surv(newdata$time, newdata$status)
    )
  }

  prepared
}


# The direct coxph fit uses exactly the design scale supplied to the final mfp2
# Cox fit. Matching nocenter and ties is important: otherwise the comparison can
# accidentally test a different Cox centering convention.

# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# The direct coxph fit uses exactly the design scale supplied to the final mfp2
# Cox fit. Matching nocenter and ties is important: otherwise the comparison can
# accidentally test a different Cox centering convention.
test_that("version 1 Cox lp and risk references match an equivalent coxph fit", {
  dat <- make_cox_reference_v1_data(seed = 11002L)
  prediction_x <- data.frame(x = c(2.75, 4.5, 7.25, 10.5))

  for (center_value in c(TRUE, FALSE)) {
    fit_mfp2 <- mfp2(
      x = as.matrix(dat["x"]),
      y = survival::Surv(dat$time, dat$status),
      family = "cox",
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = center_value,
      nocenter = NULL,
      xorder = "original",
      ties = "breslow",
      verbose = FALSE
    )

    fitted_center <- if (center_value) {
      unname(fit_mfp2$centers[[1L]])
    } else {
      0
    }

    dat_cox <- transform(dat, z = x - fitted_center)
    prediction_cox <- data.frame(z = prediction_x$x - fitted_center)

    fit_coxph <- survival::coxph(
      survival::Surv(time, status) ~ z,
      data = dat_cox,
      ties = "breslow",
      nocenter = NULL,
      x = TRUE,
      y = TRUE
    )

    for (reference_value in c("zero", "sample", "strata")) {
      got_lp <- predict(
        fit_mfp2,
        newdata = prediction_x,
        type = "lp",
        se.fit = TRUE,
        cox_reference = reference_value
      )
      expected_lp <- predict(
        fit_coxph,
        newdata = prediction_cox,
        type = "lp",
        se.fit = TRUE,
        reference = reference_value
      )

      expect_equal(
        unname(got_lp$fit),
        unname(expected_lp$fit),
        tolerance = 1e-8
      )
      expect_equal(
        unname(got_lp$se.fit),
        unname(expected_lp$se.fit),
        tolerance = 1e-8
      )

      got_risk <- predict(
        fit_mfp2,
        newdata = prediction_x,
        type = "risk",
        se.fit = TRUE,
        cox_reference = reference_value
      )
      expected_risk <- predict(
        fit_coxph,
        newdata = prediction_cox,
        type = "risk",
        se.fit = TRUE,
        reference = reference_value
      )

      expect_equal(
        unname(got_risk$fit),
        unname(expected_risk$fit),
        tolerance = 1e-8
      )
      expect_equal(
        unname(got_risk$se.fit),
        unname(expected_risk$se.fit),
        tolerance = 1e-8
      )
    }

    lp_zero <- predict(
      fit_mfp2,
      prediction_x,
      type = "lp",
      cox_reference = "zero"
    )
    lp_sample <- predict(
      fit_mfp2,
      prediction_x,
      type = "lp",
      cox_reference = "sample"
    )

    expected_constant <- sum(
      stats::coef(fit_coxph) * fit_coxph$means,
      na.rm = TRUE
    )

    expect_equal(
      unname(lp_zero - lp_sample),
      rep(unname(expected_constant), nrow(prediction_x)),
      tolerance = 1e-8
    )

    if (center_value) {
      expect_equal(
        unname(lp_zero),
        unname(lp_sample),
        tolerance = 1e-8
      )
    } else {
      expect_gt(abs(expected_constant), 1e-4)
    }
  }
})


# Offset centering is independent of the covariate reference. The zero/sample
# change must therefore remain a common covariate constant even when prediction
# offsets differ from row to row.
test_that("version 1 Cox hazard ratios are reference invariant with offsets", {
  dat <- make_cox_reference_v1_data(seed = 11003L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- data.frame(
    x = c(3.5, 8.5),
    off = c(0.10, 0.65)
  )

  lp_zero <- predict(fit, nd, type = "lp", cox_reference = "zero")
  lp_sample <- predict(fit, nd, type = "lp", cox_reference = "sample")
  risk_zero <- predict(fit, nd, type = "risk", cox_reference = "zero")
  risk_sample <- predict(fit, nd, type = "risk", cox_reference = "sample")

  # Use the stripped coxph object as the oracle. This is stronger and safer
  # than recomputing the shift with a positional coef * means expression:
  # predict.coxph() owns the exact model-matrix reconstruction and offset
  # centering rules for the fitted object.
  prepared <- prepare_formula_cox_newdata_v1(fit, nd)
  base <- strip_mfp2_class_v1(fit)

  base_lp_zero <- predict(
    base,
    newdata = prepared,
    type = "lp",
    reference = "zero"
  )
  base_lp_sample <- predict(
    base,
    newdata = prepared,
    type = "lp",
    reference = "sample"
  )
  base_risk_zero <- predict(
    base,
    newdata = prepared,
    type = "risk",
    reference = "zero"
  )
  base_risk_sample <- predict(
    base,
    newdata = prepared,
    type = "risk",
    reference = "sample"
  )

  expect_equal(unname(lp_zero), unname(base_lp_zero), tolerance = 1e-8)
  expect_equal(unname(lp_sample), unname(base_lp_sample), tolerance = 1e-8)
  expect_equal(unname(risk_zero), unname(base_risk_zero), tolerance = 1e-8)
  expect_equal(unname(risk_sample), unname(base_risk_sample), tolerance = 1e-8)

  reference_shift <- unname(base_lp_zero - base_lp_sample)
  expect_equal(
    unname(lp_zero - lp_sample),
    reference_shift,
    tolerance = 1e-8
  )
  expect_equal(
    reference_shift,
    rep(reference_shift[1L], nrow(nd)),
    tolerance = 1e-10
  )
  expect_equal(unname(diff(lp_zero)), unname(diff(lp_sample)), tolerance = 1e-10)
  expect_equal(
    unname(risk_zero[2L] / risk_zero[1L]),
    unname(risk_sample[2L] / risk_sample[1L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(risk_zero / risk_sample),
    exp(reference_shift),
    tolerance = 1e-8
  )
})


# cox_reference = "strata" uses a different weighted training mean in each stratum.
# The resulting shift is constant within a stratum, not across all prediction
# rows. Relative comparisons remain invariant only within the same stratum.
test_that("version 1 Cox strata reference is stratum specific", {
  dat <- make_cox_reference_v1_data(seed = 11004L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x + strata(stratum),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- data.frame(
    x = c(3.25, 7.75, 5.25, 10.75),
    stratum = factor(
      c("A", "A", "B", "B"),
      levels = levels(dat$stratum)
    )
  )

  lp_zero <- predict(fit, nd, type = "lp", cox_reference = "zero")
  lp_strata <- predict(fit, nd, type = "lp", cox_reference = "strata")
  risk_zero <- predict(fit, nd, type = "risk", cox_reference = "zero")
  risk_strata <- predict(fit, nd, type = "risk", cox_reference = "strata")

  prepared <- prepare_formula_cox_newdata_v1(fit, nd)
  base <- strip_mfp2_class_v1(fit)

  expect_equal(
    unname(lp_zero),
    unname(predict(base, prepared, type = "lp", reference = "zero")),
    tolerance = 1e-8
  )
  expect_equal(
    unname(lp_strata),
    unname(predict(base, prepared, type = "lp", reference = "strata")),
    tolerance = 1e-8
  )

  reference_shift <- unname(lp_zero - lp_strata)

  expect_equal(reference_shift[1L], reference_shift[2L], tolerance = 1e-10)
  expect_equal(reference_shift[3L], reference_shift[4L], tolerance = 1e-10)
  expect_gt(abs(reference_shift[1L] - reference_shift[3L]), 1e-4)

  expect_equal(
    unname(lp_zero[2L] - lp_zero[1L]),
    unname(lp_strata[2L] - lp_strata[1L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(lp_zero[4L] - lp_zero[3L]),
    unname(lp_strata[4L] - lp_strata[3L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(risk_zero[2L] / risk_zero[1L]),
    unname(risk_strata[2L] / risk_strata[1L]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(risk_zero[4L] / risk_zero[3L]),
    unname(risk_strata[4L] / risk_strata[3L]),
    tolerance = 1e-10
  )
})


# This is the previously uncovered end-to-end baseline-hazard path. The oracle
# is the same fitted coxph object after only the mfp2 class is removed; newdata
# are independently reconstructed on the stored transformed design scale.
test_that("version 1 Cox expected and survival newdata predictions match coxph", {
  dat <- make_cox_reference_v1_data(seed = 11005L)

  fit <- mfp2(
    survival::Surv(time, status) ~
      x + strata(stratum) + offset(off),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[c(7, 41, 88, 133, 207, 319),
            c("time", "status", "x", "stratum", "off"),
            drop = FALSE]

  prepared <- prepare_formula_cox_newdata_v1(
    fit,
    nd,
    include_response = TRUE
  )
  base <- strip_mfp2_class_v1(fit)

  for (prediction_type in c("expected", "survival")) {
    got <- predict(
      fit,
      newdata = nd,
      type = prediction_type,
      se.fit = TRUE
    )
    expected <- predict(
      base,
      newdata = prepared,
      type = prediction_type,
      se.fit = TRUE
    )

    expect_equal(
      unname(got$fit),
      unname(expected$fit),
      tolerance = 1e-8
    )
    expect_equal(
      unname(got$se.fit),
      unname(expected$se.fit),
      tolerance = 1e-8
    )
  }

  got_expected <- predict(fit, nd, type = "expected")
  got_survival <- predict(fit, nd, type = "survival")
  expect_equal(
    unname(got_survival),
    unname(exp(-got_expected)),
    tolerance = 1e-10
  )
})


test_that("version 1 Cox expected and survival training predictions match coxph", {
  dat <- make_cox_reference_v1_data(seed = 11006L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = TRUE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  base <- strip_mfp2_class_v1(fit)

  for (prediction_type in c("expected", "survival")) {
    expect_equal(
      unname(predict(fit, type = prediction_type)),
      unname(predict(base, type = prediction_type)),
      tolerance = 1e-8
    )
  }
})


# Absolute predictions derive their response from newdata for every interface.
# A matrix-interface caller supplies one Surv column alongside the predictors;
# the column is removed before FP transformation and reattached to the internal
# Cox prediction frame under the response name expected by predict.coxph().
test_that("version 1 matrix-interface absolute Cox predictions derive response from newdata", {
  dat <- make_cox_reference_v1_data(seed = 11007L)
  x <- as.matrix(dat["x"])
  y <- survival::Surv(dat$time, dat$status)

  fit <- mfp2(
    x = x,
    y = y,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  rows <- c(4, 29, 74, 151)
  newdata <- data.frame(x = x[rows, "x"])
  newdata$prediction_response <- I(y[rows])

  prepared <- prepare_newdata_for_predict(
    fit,
    newdata,
    terms = get_selected_variable_names(fit),
    check_binary = FALSE
  )
  response_name <- cox_internal_response_name(fit)
  prepared[[response_name]] <- I(y[rows])

  base <- strip_mfp2_class_v1(fit)

  expect_equal(
    unname(predict(fit, newdata = newdata, type = "expected")),
    unname(predict(base, newdata = prepared, type = "expected")),
    tolerance = 1e-8
  )
  expect_equal(
    unname(predict(fit, newdata = newdata, type = "survival")),
    unname(predict(base, newdata = prepared, type = "survival")),
    tolerance = 1e-8
  )

  expect_error(
    predict(fit, newdata = data.frame(x = newdata$x), type = "expected"),
    "require.*follow-up response information"
  )

  ambiguous <- newdata
  ambiguous$second_response <- I(y[rows])
  expect_error(
    predict(fit, newdata = ambiguous, type = "survival"),
    "more than one Surv column"
  )
})


# The public default remains the reference-zero scale used by mfp2 term
# decomposition. An explicit sample reference is allowed, but it introduces the
# documented common constant when transformed columns have nonzero means.
test_that("version 1 default Cox lp remains aligned with mfp2 terms", {
  dat <- make_cox_reference_v1_data(seed = 11008L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd <- data.frame(x = c(3, 5, 7, 9))
  term_result <- predict(
    fit,
    newdata = nd,
    type = "terms",
    terms_seq = "data"
  )

  term_sum <- Reduce(`+`, lapply(term_result, `[[`, "value"))
  lp_default <- predict(fit, nd, type = "lp")
  lp_null <- predict(fit, nd, type = "lp", cox_reference = NULL)
  lp_zero <- predict(fit, nd, type = "lp", cox_reference = "zero")
  lp_sample <- predict(fit, nd, type = "lp", cox_reference = "sample")

  expect_equal(unname(lp_default), unname(lp_null), tolerance = 1e-10)
  expect_equal(unname(lp_default), unname(lp_zero), tolerance = 1e-10)
  expect_equal(unname(lp_zero), unname(term_sum), tolerance = 1e-8)
  expect_equal(
    unname(lp_zero - lp_sample),
    rep(unname(lp_zero[1L] - lp_sample[1L]), nrow(nd)),
    tolerance = 1e-10
  )
})


test_that("version 1 reports clear Cox reference and response errors", {
  dat <- make_cox_reference_v1_data(seed = 11009L)

  fit <- mfp2(
    survival::Surv(time, status) ~ x,
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    nocenter = NULL,
    xorder = "original",
    ties = "breslow",
    verbose = FALSE
  )

  nd_relative <- data.frame(x = dat$x[1:4])
  nd_absolute <- dat[1:4, c("time", "status", "x"), drop = FALSE]

  # This was formerly the unhelpful duplicate-formal-argument crash.
  expect_no_error(
    predict(fit, nd_relative, type = "lp", cox_reference = "sample")
  )
  expect_no_error(
    predict(fit, type = "lp", cox_reference = "sample")
  )
  base <- strip_mfp2_class_v1(fit)
  expect_equal(
    unname(predict(fit, type = "lp", cox_reference = "sample")),
    unname(predict(base, type = "lp", reference = "sample")),
    tolerance = 1e-8
  )

  expect_error(
    predict(fit, nd_relative, type = "lp", cox_reference = "population"),
    "must be one of 'zero', 'sample', or 'strata'"
  )
  expect_error(
    predict(fit, nd_relative, type = "terms", cox_reference = "zero"),
    "not used for mfp2 term or contrast predictions"
  )
  expect_error(
    predict(fit, nd_absolute, type = "survival", cox_reference = "zero"),
    "does not apply to Cox predictions"
  )
  expect_error(
    predict(fit, nd_relative, type = "survival"),
    "require.*follow-up response"
  )
  expect_error(
    predict(
      fit,
      nd_relative,
      type = "lp",
      newy = survival::Surv(dat$time[1:4], dat$status[1:4])
    ),
    "'newy' has been removed"
  )

  expect_error(
    predict(fit, nd_relative, type = "lp", reference = "zero"),
    "renamed.*cox_reference"
  )

  fit_glm <- mfp2(
    x = as.matrix(dat["x"]),
    # This fit only supplies a valid non-Cox object for the error check below.
    # Use deterministic non-linear variation instead of consuming RNG state.
    y = dat$x + 0.1 * sin(seq_len(nrow(dat))),
    verbose = FALSE
  )
  expect_error(
    predict(fit_glm, nd_relative, cox_reference = "zero"),
    "only available for Cox models"
  )
})
