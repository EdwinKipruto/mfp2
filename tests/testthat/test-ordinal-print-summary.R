# Tests for ordinal threshold intercept printing in print.mfp2() and
# summary.mfp2() / print.summary.mfp2().
#
# Verifies:
#   - print.mfp2()         shows threshold rows above predictor rows by default
#   - print.mfp2()         suppresses thresholds when intercepts = FALSE
#   - print.mfp2()         shows thresholds when intercepts = TRUE even for
#                           models that would default to FALSE (>= 10 thresholds)
#   - summary.mfp2()       includes an ordinal_intercepts data frame
#   - print.summary.mfp2() emits an "Ordinal Intercepts" section
#   - SE and CI are present when vcov is available
#   - Non-ordinal families are unaffected

skip_if_no_rms <- function() skip_if_not_installed("rms")

make_ordinal_data <- function(n = 300, seed = 7, k = 4) {
  set.seed(seed)
  age <- runif(n, 20, 80)
  sex <- rbinom(n, 1, 0.5)
  eta <- 0.04 * age + 0.5 * sex
  # Cut the same vector used to derive the breaks so no observation falls
  # outside them and produces NA.
  z <- eta + rlogis(n)
  y <- as.integer(cut(z, breaks = k, include.lowest = TRUE))
  data.frame(y = ordered(y), age = age, sex = sex)
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

local_ordinal_fit <- function(n = 300, k = 4, link = "logistic") {
  dm <- make_ordinal_data(n = n, k = k)
  mfp2(
    y ~ fp(age, df = 1, select = 1) + fp(sex, df = 1, select = 1),
    data = dm, family = ordinal_family(link = link),
    center = FALSE, verbose = FALSE,
    alpha = 1, select = 1
  )
}

# ---------------------------------------------------------------------------
# print.mfp2() — threshold rows
# ---------------------------------------------------------------------------

test_that("print.mfp2() shows threshold rows by default for ordinal model", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  out <- capture.output(print(fit))

  # At least one "y>=" row must appear in the coefficient section.
  coef_section_start <- grep("Final Model Coefficients", out)
  expect_true(length(coef_section_start) > 0L)
  coef_lines <- out[seq(coef_section_start[1L], length(out))]
  expect_true(any(grepl("y>=", coef_lines, fixed = TRUE)))
  expect_true(any(grepl(
    "Model: Logistic Ordinal Regression (proportional odds)",
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("Frequencies of Responses", out, fixed = TRUE)))
  expect_equal(sum(mfp2:::mfp2_summary_response_frequencies(fit)), fit$nobs)
})

test_that("print.mfp2() intercepts = FALSE suppresses threshold rows", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  out <- capture.output(print(fit, intercepts = FALSE))

  coef_section_start <- grep("Final Model Coefficients", out)
  coef_lines <- out[seq(coef_section_start[1L], length(out))]
  expect_false(any(grepl("y>=", coef_lines, fixed = TRUE)))
})

test_that("print.mfp2() intercepts = TRUE forces threshold rows even if >= 10", {
  skip_on_cran(); skip_if_no_rms()
  # Build a model with many outcome levels so the default would suppress them.
  set.seed(42)
  n <- 480L
  dm <- data.frame(
    y = ordered(sample(rep(seq_len(12L), length.out = n))),
    age = runif(n, 20, 80)
  )
  fit <- mfp2(
    y ~ fp(age, df = 1, select = 1),
    data = dm, family = ordinal_family(),
    center = FALSE, verbose = FALSE, alpha = 1, select = 1
  )
  expect_gte(length(fit$mfp2_ordinal_intercepts), 10L)
  # Default should suppress (>= 10 thresholds).
  out_default <- capture.output(print(fit))
  coef_start  <- grep("Final Model Coefficients", out_default)
  coef_lines  <- out_default[seq(coef_start[1L], length(out_default))]
  expect_false(any(grepl("y>=", coef_lines, fixed = TRUE)))

  # Explicit TRUE must show them.
  out_forced <- capture.output(print(fit, intercepts = TRUE))
  coef_start2 <- grep("Final Model Coefficients", out_forced)
  coef_lines2  <- out_forced[seq(coef_start2[1L], length(out_forced))]
  expect_true(any(grepl("y>=", coef_lines2, fixed = TRUE)))
})

test_that("print.mfp2() threshold rows appear before predictor rows", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  out <- capture.output(print(fit))

  # Find the first threshold and first predictor row in output.
  coef_start   <- grep("Final Model Coefficients", out)[1L]
  coef_lines   <- out[seq(coef_start, length(out))]
  first_thresh <- min(grep("y>=", coef_lines, fixed = TRUE))
  # Predictor rows: lines that are NOT thresholds, NOT the header, and NOT
 # section separators, but DO contain a predictor name at the start.
  pred_names   <- unique(sub("\\.1$", "", names(stats::coef(fit))))
  is_pred_line <- grepl(
    paste0("^\\s*(", paste(pred_names, collapse = "|"), ")\\b"),
    coef_lines
  ) & !grepl("y>=", coef_lines, fixed = TRUE)
  first_pred   <- min(which(is_pred_line))
  expect_lt(first_thresh, first_pred)
})

test_that("print.mfp2() validates the intercepts argument", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()

  expect_error(print(fit, intercepts = NA), "intercepts")
  expect_error(print(fit, intercepts = c(TRUE, FALSE)), "intercepts")
  expect_error(print(fit, intercepts = 1), "intercepts")
})

test_that("print.mfp2() threshold SE is non-NA when vcov is available", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  expect_true(is.matrix(vcov(fit)))

  out <- capture.output(print(fit, intercepts = TRUE))
  # The table should not contain bare "NA" in the SE column next to y>= rows.
  thresh_lines <- grep("y>=", out, fixed = TRUE, value = TRUE)
  expect_true(length(thresh_lines) > 0L)
  # Remove the numeric category embedded in the threshold label before counting
  # estimate and SE tokens.
  has_two_numbers <- vapply(thresh_lines, function(ln) {
    values <- sub("^\\s*y>=\\d+\\s+", "", ln)
    length(regmatches(
      values,
      gregexpr(
        "-?[0-9]+(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?",
        values,
        perl = TRUE
      )
    )[[1L]]) >= 2L
  }, logical(1L))
  expect_true(all(has_two_numbers))
})

# ---------------------------------------------------------------------------
# summary.mfp2() — ordinal_intercepts element
# ---------------------------------------------------------------------------

test_that("summary.mfp2() includes ordinal_intercepts for ordinal family", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  emitted <- capture.output(sm <- summary(fit))

  expect_length(emitted, 0L)
  expect_false(is.null(sm$ordinal_intercepts))
  expect_s3_class(sm$ordinal_intercepts, "data.frame")
  expected_cols <- c("threshold", "estimate", "se", "z", "p", "ci_lower", "ci_upper")
  expect_true(all(expected_cols %in% names(sm$ordinal_intercepts)))

  # Number of rows must equal k - 1.
  n_int <- length(fit$mfp2_ordinal_intercepts)
  expect_equal(nrow(sm$ordinal_intercepts), n_int)
})

test_that("summary.mfp2() ordinal_intercepts estimates match stored intercepts", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  sm  <- summary(fit)

  expect_equal(
    sm$ordinal_intercepts$estimate,
    unname(fit$mfp2_ordinal_intercepts),
    tolerance = 1e-8
  )
})

test_that("summary.mfp2() ordinal_intercepts has finite SE when vcov available", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  sm  <- summary(fit)

  expect_true(all(is.finite(sm$ordinal_intercepts$se)))
  expect_true(all(sm$ordinal_intercepts$se > 0))
})

test_that("summary.mfp2() ordinal_intercepts CI width is positive", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  sm  <- summary(fit)
  oi  <- sm$ordinal_intercepts
  expect_true(all(oi$ci_upper > oi$ci_lower))
})

test_that("summary.mfp2() reports finite ordinal slope inference", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  sm <- summary(fit)

  expect_gt(nrow(sm$linear_terms), 0L)
  expect_true(all(is.finite(sm$linear_terms$se)))
  expect_true(all(sm$linear_terms$se > 0))
  expect_true(all(is.finite(sm$linear_terms$statistic)))
  expect_true(all(is.finite(sm$linear_terms$p)))
})

test_that("ordinal slope inference uses the slope covariance block", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  slope_covariance <- mfp2:::mfp2_ordinal_vcov_slopes(fit)
  expected <- stats::vcov(fit, intercepts = "none")

  expect_equal(slope_covariance, expected)
  expect_true(all(is.finite(diag(slope_covariance))))
  expect_true(all(diag(slope_covariance) > 0))
})

test_that("five-predictor ordinal summaries have finite slope inference", {
  skip_on_cran(); skip_if_no_rms()
  set.seed(99)
  n <- 600L
  predictors <- replicate(5L, stats::rnorm(n))
  colnames(predictors) <- paste0("x", seq_len(5L))
  eta <- drop(predictors %*% c(0.45, -0.35, 0.25, 0.15, -0.2))
  latent <- eta + stats::rlogis(n)
  dat <- data.frame(
    y = ordered(cut(
      latent,
      breaks = stats::quantile(latent, probs = seq(0, 1, length.out = 5L)),
      include.lowest = TRUE
    )),
    predictors,
    check.names = FALSE
  )
  fit <- mfp2(
    y ~ x1 + x2 + x3 + x4 + x5,
    data = dat,
    df = 1,
    center = FALSE,
    family = ordinal_family(link = "logistic"),
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  sm <- summary(fit)

  expect_equal(nrow(sm$linear_terms), 5L)
  expect_true(all(is.finite(sm$linear_terms$se)))
  expect_true(all(is.finite(sm$linear_terms$statistic)))
  expect_true(all(is.finite(sm$linear_terms$p)))
})

test_that("summary.mfp2() excludes every ordinal threshold from model df", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  sm <- summary(fit)
  model_fit <- sm$fit$model_fit

  expected_linear_df <- length(fit$coefficients) -
    length(fit$mfp2_ordinal_intercepts)
  expect_equal(model_fit$df[model_fit$label == "Full linear model"],
               expected_linear_df)
})

test_that("summary.mfp2() ordinal_intercepts is NULL for non-ordinal families", {
  skip_on_cran()
  dm  <- data.frame(y = rbinom(100, 1, 0.5), x = rnorm(100))
  fit <- mfp2(y ~ fp(x, df = 1), data = dm, family = binomial(),
              verbose = FALSE, alpha = 1, select = 1)
  sm  <- summary(fit)
  expect_null(sm$ordinal_intercepts)
})

# ---------------------------------------------------------------------------
# print.summary.mfp2() — Ordinal Intercepts section
# ---------------------------------------------------------------------------

test_that("print.summary.mfp2() emits 'Ordinal Intercepts' section for ordinal", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  out <- capture.output(print(summary(fit)))

  expect_true(any(grepl("Ordinal Intercepts", out, fixed = TRUE)))
  # The section must contain threshold row(s).
  thresh_lines <- grep("y>=", out, fixed = TRUE, value = TRUE)
  expect_true(length(thresh_lines) > 0L)
})

test_that("print.summary.mfp2() 'Ordinal Intercepts' section comes before 'Linear Terms'", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  out <- capture.output(print(summary(fit)))

  ord_line    <- grep("Ordinal Intercepts", out, fixed = TRUE)[1L]
  linear_line <- grep("Linear Terms",       out, fixed = TRUE)[1L]
  expect_lt(ord_line, linear_line)
})

test_that("print.summary.mfp2() 'Ordinal Intercepts' includes 95% CI column", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit()
  out <- capture.output(print(summary(fit)))

  ordinal_start <- grep("Ordinal Intercepts", out, fixed = TRUE)[1L]
  linear_start <- grep("Linear Terms", out, fixed = TRUE)[1L]
  ordinal_section <- out[seq.int(ordinal_start, linear_start - 1L)]
  expect_true(any(grepl("[95% CI]", ordinal_section, fixed = TRUE)))
})

test_that("print.summary.mfp2() labels non-logistic ordinal links correctly", {
  skip_on_cran(); skip_if_no_rms()
  fit <- local_ordinal_fit(link = "probit")
  out <- capture.output(print(summary(fit)))

  expect_true(any(grepl(
    "Thresholds for the probit cumulative-link model.",
    out,
    fixed = TRUE
  )))
  expect_false(any(grepl("log-odds", out, fixed = TRUE)))
  expect_true(any(grepl(
    "Model: Probit Ordinal Regression",
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("Frequencies of Responses", out, fixed = TRUE)))
})

test_that("print.summary.mfp2() does not emit 'Ordinal Intercepts' for GLM", {
  skip_on_cran()
  dm  <- data.frame(y = rpois(100, 2), x = rnorm(100))
  fit <- mfp2(y ~ fp(x, df = 1), data = dm, family = poisson(),
              verbose = FALSE, alpha = 1, select = 1)
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("Ordinal Intercepts", out, fixed = TRUE)))
})
