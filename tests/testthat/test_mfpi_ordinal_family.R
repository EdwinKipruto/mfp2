# End-to-end coverage for ordinal (proportional-odds) MFPI support.

# ---------------------------------------------------------------------------
# Shared data generator
# ---------------------------------------------------------------------------

make_mfpi_ordinal_data <- function(n = 300L, seed = 4201L) {
  set.seed(seed)
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2L)),
    x = stats::runif(n, 0.5, 5),
    z = stats::rnorm(n)
  )
  # Proportional-odds DGP: three ordered categories
  eta <- with(
    dat,
    0.3 * x + 0.4 * (group == "treated") +
      0.35 * x * (group == "treated") - 0.2 * z
  )
  cum_prob2 <- stats::plogis(1.0 + eta)
  cum_prob3 <- stats::plogis(-0.5 + eta)
  u <- stats::runif(n)
  dat$y <- ordered(
    ifelse(u < 1 - cum_prob2, 1L,
      ifelse(u < 1 - cum_prob3, 2L, 3L)
    ),
    levels = 1:3
  )
  dat
}

fit_mfpi_ordinal <- function(dat, link = "logistic", verbose = FALSE) {
  mfpi(
    y ~ group + fp(x, df = 1) + z,
    data = dat,
    family = ordinal_family(link = link),
    group_var = "group",
    interaction_vars = "x",
    interaction_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    p_interact = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = verbose
  )
}

# Reassemble the ordinal class-probability matrix from the flattened
# `fit.<level>` columns returned by predict(type = "response").
ordinal_response_matrix <- function(pred) {
  df <- pred$predictions
  fit_cols <- grep("^fit(\\.|$)", names(df), value = TRUE)
  as.matrix(df[, fit_cols, drop = FALSE])
}


# ---------------------------------------------------------------------------
# 1. Basic fitting
# ---------------------------------------------------------------------------

test_that("ordinal mfpi() fits without error and returns correct class", {
  skip_if_not_installed("rms")

  output <- capture.output(
    fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data(), verbose = TRUE)
  )

  expect_s3_class(fit, "mfpi")
  expect_identical(fit$family_string, "ordinal")
  expect_true(is.character(fit$ordinal_levels))
  expect_true(length(fit$ordinal_levels) >= 3L)
  expect_true(is.character(fit$ordinal_link))
  expect_true(is.integer(fit$n_intercepts))
  expect_equal(fit$n_intercepts, length(fit$ordinal_levels) - 1L)
  expect_match(
    paste(output, collapse = "\n"),
    "Model: Logistic ordinal regression (proportional odds)",
    fixed = TRUE
  )
})


test_that("ordinal mfpi() stores ordinal metadata correctly", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data(), link = "probit")

  expect_identical(fit$ordinal_link, "probit")
  expect_identical(fit$ordinal_levels, c("1", "2", "3"))
  expect_identical(fit$n_intercepts, 2L)
})


# ---------------------------------------------------------------------------
# 2. Coefficient and vcov extraction (accessor API uses `term`, singular)
# ---------------------------------------------------------------------------

test_that("coef() and vcov() work for ordinal mfpi", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())

  cf <- coef(fit, term = "x")
  expect_true(is.numeric(cf))
  expect_true(length(cf) > 0L)

  vc <- vcov(fit, term = "x")
  expect_true(is.matrix(vc))
  expect_equal(nrow(vc), length(cf))
  expect_equal(ncol(vc), length(cf))
  expect_identical(rownames(vc), names(cf))
})


# ---------------------------------------------------------------------------
# 3. Function prediction (fitted group functions / differences)
# ---------------------------------------------------------------------------

test_that("function prediction works for ordinal mfpi", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:20, c("group", "x", "z"), drop = FALSE]

  # type = "function" — group-specific fitted functions (one row per group)
  pred_fun <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "function", se.fit = TRUE
  )
  expect_true(is.data.frame(pred_fun$functions))
  expect_true(all(c("fit", "se.fit") %in% names(pred_fun$functions)))
  # Two groups, so 2 * nrow(nd) fitted-function rows.
  expect_equal(nrow(pred_fun$functions), 2L * nrow(nd))
  expect_true(all(is.finite(pred_fun$functions$fit)))
  expect_true(all(pred_fun$functions$se.fit >= 0))

  # type = "difference" — comparison-group minus reference-group
  pred_diff <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "difference", se.fit = TRUE
  )
  expect_true(is.data.frame(pred_diff$differences))
  expect_true("fit" %in% names(pred_diff$differences))
  expect_equal(nrow(pred_diff$differences), nrow(nd))
  expect_true(all(is.finite(pred_diff$differences$fit)))
})


# ---------------------------------------------------------------------------
# 4. Ordinary prediction — link
# ---------------------------------------------------------------------------

test_that("ordinal link prediction returns the slope linear predictor", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:20, c("group", "x", "z"), drop = FALSE]

  got <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )

  expect_true(is.numeric(got$predictions$fit))
  expect_equal(length(got$predictions$fit), nrow(nd))
  expect_false(got$metadata$used_model_predict)

  # Standard errors present and non-negative.
  expect_true(!is.null(got$predictions$se.fit))
  expect_true(all(got$predictions$se.fit >= 0))
})


# ---------------------------------------------------------------------------
# 5. Ordinary prediction — response (class probabilities)
# ---------------------------------------------------------------------------

test_that("ordinal response prediction returns valid class probabilities", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:20, c("group", "x", "z"), drop = FALSE]

  got <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "response", se.fit = FALSE
  )

  probs <- ordinal_response_matrix(got)
  expect_equal(nrow(probs), nrow(nd))
  expect_equal(ncol(probs), length(fit$ordinal_levels))

  # Probabilities are valid and rows sum to 1.
  expect_true(all(probs >= -1e-10))
  expect_true(all(probs <= 1 + 1e-10))
  expect_equal(rowSums(probs), rep(1, nrow(nd)), tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# 6. Ordinary prediction — mean (expected response)
# ---------------------------------------------------------------------------

test_that("ordinal mean prediction returns numeric expectations", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:20, c("group", "x", "z"), drop = FALSE]

  got_mean <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "mean", se.fit = FALSE
  )

  means <- got_mean$predictions$fit
  expect_true(is.numeric(means))
  expect_equal(length(means), nrow(nd))

  numeric_levels <- as.numeric(fit$ordinal_levels)
  expect_true(all(means >= min(numeric_levels) - 1e-10))
  expect_true(all(means <= max(numeric_levels) + 1e-10))
})


test_that("ordinal mean is consistent with response probabilities", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:20, c("group", "x", "z"), drop = FALSE]

  got_resp <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "response", se.fit = FALSE
  )
  got_mean <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "mean", se.fit = FALSE
  )

  probs <- ordinal_response_matrix(got_resp)
  numeric_levels <- as.numeric(fit$ordinal_levels)
  manual_mean <- as.numeric(probs %*% numeric_levels)

  expect_equal(got_mean$predictions$fit, manual_mean, tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# 7. Probit link
# ---------------------------------------------------------------------------

test_that("ordinal mfpi with probit link returns valid predictions", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data(), link = "probit")
  expect_identical(fit$ordinal_link, "probit")

  nd <- make_mfpi_ordinal_data()[1:15, c("group", "x", "z"), drop = FALSE]

  got_resp <- predict(
    fit, newdata = nd, terms = "x", model = "all",
    type = "response", se.fit = FALSE
  )
  probs <- ordinal_response_matrix(got_resp)
  expect_true(all(probs >= -1e-10))
  expect_true(all(probs <= 1 + 1e-10))
  expect_equal(rowSums(probs), rep(1, nrow(nd)), tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# 8. Matrix interface
# ---------------------------------------------------------------------------

test_that("ordinal mfpi works through the matrix interface", {
  skip_if_not_installed("rms")

  dat <- make_mfpi_ordinal_data()
  fit <- mfpi(
    dat[, c("group", "x", "z")],
    dat$y,
    family = ordinal_family(),
    group_var = "group",
    interaction_vars = "x",
    interaction_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    p_interact = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_identical(fit$family_string, "ordinal")
  expect_true(length(coef(fit, term = "x")) > 0L)
})


# ---------------------------------------------------------------------------
# 9. se.fit warning for non-link types
# ---------------------------------------------------------------------------

test_that("ordinal se.fit warns for response/mean types", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:10, c("group", "x", "z"), drop = FALSE]

  expect_warning(
    predict(
      fit, newdata = nd, terms = "x", model = "all",
      type = "response", se.fit = TRUE
    ),
    "se.fit.*not available.*ordinal"
  )
})


# ---------------------------------------------------------------------------
# 10. type validation
# ---------------------------------------------------------------------------

test_that("ordinal mfpi rejects invalid prediction types", {
  skip_if_not_installed("rms")

  fit <- fit_mfpi_ordinal(make_mfpi_ordinal_data())
  nd <- make_mfpi_ordinal_data()[1:5, c("group", "x", "z"), drop = FALSE]

  expect_error(
    predict(fit, newdata = nd, terms = "x", model = "all", type = "risk"),
    "type.*must be one of"
  )
})
