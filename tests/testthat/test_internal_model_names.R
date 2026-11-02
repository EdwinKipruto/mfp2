# Collision-safe internal model-frame names -----------------------------------
#
# These tests deliberately use predictor names that overlap with helper columns
# created by the formula-based final refit. The MFP search itself remains
# matrix-based; only the final glm()/coxph() model frame needs collision-safe
# response, offset, and strata names.


test_that("internal name allocator avoids preferred and ..mfp2 collisions", {
  used <- c("y", "offset_", "strata_", "..mfp2_response", "..mfp2_response_1")

  expect_identical(
    mfp2_internal_name("response", used, preferred = "y"),
    "..mfp2_response_2"
  )
  expect_identical(
    mfp2_internal_name("offset", used, preferred = "offset_"),
    "..mfp2_offset"
  )
  expect_identical(
    mfp2_internal_name("strata", used, preferred = "strata_"),
    "..mfp2_strata"
  )
})

test_that("grouped-binomial response helpers avoid predictor collisions", {
  set.seed(32004)
  n <- 240L

  x <- cbind(
    `..mfp2_successes` = stats::rnorm(n),
    `..mfp2_failures` = stats::rnorm(n)
  )
  trials <- sample(8:20, n, replace = TRUE)
  eta <- -0.3 + 0.45 * x[, "..mfp2_successes"] -
    0.35 * x[, "..mfp2_failures"]
  successes <- stats::rbinom(n, size = trials, prob = stats::plogis(eta))
  grouped_y <- cbind(successes, trials - successes)

  fit <- mfp2(
    x = x,
    y = grouped_y,
    family = "binomial",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  internal <- fit$mfp2_internal_names
  expect_length(internal$response, 2L)
  expect_false(any(internal$response %in% colnames(x)))

  dat <- as.data.frame(x, check.names = FALSE)
  dat$successes <- grouped_y[, 1L]
  dat$failures <- grouped_y[, 2L]
  ref <- stats::glm(
    cbind(successes, failures) ~ ..mfp2_successes + ..mfp2_failures,
    data = dat,
    family = stats::binomial()
  )

  expect_equal(unname(stats::coef(fit)), unname(stats::coef(ref)), tolerance = 1e-8)
})


test_that("mfp2 GLM offset helper cannot overwrite colliding predictors", {
  set.seed(32001)
  n <- 220L

  x <- cbind(
    y = stats::rnorm(n),
    offset_ = stats::rnorm(n),
    `..mfp2_y` = stats::rnorm(n),
    `..mfp2_offset` = stats::rnorm(n)
  )
  off <- stats::rnorm(n, sd = 0.20)
  response <- 0.4 + 0.35 * x[, "y"] - 0.25 * x[, "offset_"] +
    0.15 * x[, "..mfp2_y"] + 0.20 * x[, "..mfp2_offset"] + off +
    stats::rnorm(n, sd = 0.45)

  fit <- mfp2(
    x = x,
    y = response,
    family = "gaussian",
    offset = off,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  internal <- fit$mfp2_internal_names
  expect_true(is.list(internal))
  expect_false(internal$response %in% colnames(x))
  expect_false(internal$offset %in% colnames(x))
  expect_false(identical(internal$offset, "offset_"))

  dat <- as.data.frame(x, check.names = FALSE)
  dat$response <- response
  dat$off_external <- off
  ref <- stats::glm(
    response ~ y + offset_ + ..mfp2_y + ..mfp2_offset + offset(off_external),
    data = dat,
    family = stats::gaussian()
  )

  expect_equal(unname(stats::coef(fit)), unname(stats::coef(ref)), tolerance = 1e-8)

  rows <- seq_len(25L)
  nd <- as.data.frame(x[rows, , drop = FALSE], check.names = FALSE)
  got <- stats::predict(
    fit,
    newdata = nd,
    newoffset = off[rows],
    type = "link"
  )
  expected <- stats::predict(
    ref,
    newdata = transform(nd, off_external = off[rows]),
    type = "link"
  )
  expect_equal(as.numeric(got), as.numeric(expected), tolerance = 1e-8)

  prepared <- prepare_newdata_for_predict(
    fit,
    nd,
    offset = off[rows],
    check_binary = FALSE
  )
  # Prediction data contain transformed predictor columns plus the one stored
  # internal offset helper. The raw predictor names themselves need not remain
  # after FP/linear transformation; the important invariant is that no column
  # was overwritten and the fitted helper name is present exactly once.
  expect_true("offset_.1" %in% names(prepared))
  expect_true("..mfp2_offset.1" %in% names(prepared))
  expect_true(internal$offset %in% names(prepared))
  expect_false(anyDuplicated(names(prepared)) > 0L)
})


test_that("mfp2 Cox response offset and strata helpers avoid predictor collisions", {
  set.seed(32002)
  n <- 320L

  x <- cbind(
    y = stats::rnorm(n),
    offset_ = stats::rnorm(n),
    strata_ = stats::rnorm(n),
    `..mfp2_response` = stats::rnorm(n),
    `..mfp2_offset` = stats::rnorm(n),
    `..mfp2_strata` = stats::rnorm(n)
  )
  stratum <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  off <- stats::rnorm(n, sd = 0.15)
  eta <- 0.25 * x[, "y"] - 0.20 * x[, "offset_"] +
    0.15 * x[, "strata_"] + 0.10 * x[, "..mfp2_response"] -
    0.12 * x[, "..mfp2_offset"] + 0.08 * x[, "..mfp2_strata"] + off
  event_time <- stats::rexp(n, rate = exp(eta) / 10)
  censor_time <- stats::rexp(n, rate = 0.07)
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  surv_y <- survival::Surv(time, status)

  fit <- mfp2(
    x = x,
    y = surv_y,
    family = "cox",
    offset = off,
    strata = stratum,
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

  internal <- fit$mfp2_internal_names
  expect_true(is.list(internal))
  expect_false(internal$response %in% colnames(x))
  expect_false(internal$offset %in% colnames(x))
  expect_false(internal$strata %in% colnames(x))
  expect_false(identical(internal$response, "y"))
  expect_false(identical(internal$offset, "offset_"))
  expect_false(identical(internal$strata, "strata_"))
  expect_true(inherits(fit$y, "Surv"))

  dat <- as.data.frame(x, check.names = FALSE)
  dat$time <- time
  dat$status <- status
  dat$off_external <- off
  dat$stratum_external <- stratum
  ref <- survival::coxph(
    survival::Surv(time, status) ~ y + offset_ + strata_ +
      ..mfp2_response + ..mfp2_offset + ..mfp2_strata +
      offset(off_external) + strata(stratum_external),
    data = dat,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )

  expect_equal(unname(stats::coef(fit)), unname(stats::coef(ref)), tolerance = 1e-8)
  expect_equal(unname(fit$loglik), unname(ref$loglik), tolerance = 1e-8)

  rows <- seq_len(30L)
  nd <- as.data.frame(x[rows, , drop = FALSE], check.names = FALSE)
  got <- stats::predict(
    fit,
    newdata = nd,
    newoffset = off[rows],
    strata = stratum[rows],
    type = "lp",
    cox_reference = "zero"
  )
  nd_ref <- transform(
    nd,
    off_external = off[rows],
    stratum_external = stratum[rows]
  )
  expected <- stats::predict(
    ref,
    newdata = nd_ref,
    type = "lp",
    reference = "zero"
  )
  expect_equal(as.numeric(got), as.numeric(expected), tolerance = 1e-8)

  prepared <- prepare_newdata_for_predict(
    fit,
    nd,
    strata = stratum[rows],
    offset = off[rows],
    check_binary = FALSE
  )
  expect_true("offset_.1" %in% names(prepared))
  expect_true("strata_.1" %in% names(prepared))
  expect_true(internal$offset %in% names(prepared))
  expect_true(internal$strata %in% names(prepared))
  expect_false(anyDuplicated(names(prepared)) > 0L)
})


test_that("formula-based prediction requires stored internal-name metadata", {
  set.seed(32003)
  n <- 150L
  x <- cbind(x1 = stats::rnorm(n))
  off <- stats::rnorm(n, sd = 0.10)
  stratum <- factor(sample(c("A", "B"), n, replace = TRUE))
  event_time <- stats::rexp(n, rate = exp(0.3 * x[, 1] + off) / 8)
  censor_time <- stats::rexp(n, rate = 0.08)
  y <- survival::Surv(
    pmin(event_time, censor_time),
    as.integer(event_time <= censor_time)
  )

  fit <- mfp2(
    x, y,
    family = "cox",
    offset = off,
    strata = stratum,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    ties = "breslow",
    verbose = FALSE
  )
  fit$mfp2_internal_names <- NULL

  expect_error(
    predict(
      fit,
      newdata = as.data.frame(x[1:8, , drop = FALSE]),
      newoffset = off[1:8],
      strata = stratum[1:8],
      type = "lp"
    ),
    "lacks internal-name metadata"
  )
})


# MFPI uses the same stored-name accessor in its ordinary prediction design.
# Its end-to-end stratified Cox prediction path is exercised in test_mfp2.R;
# keeping this file focused on helper-name collisions avoids coupling the
# namespace regression test to unrelated MFPI adjustment-model selection.
