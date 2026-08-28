# End-to-end test for the full linear reference when a spike term is retained.

# A retained spike implies catzero and zero. The full linear reference must
# therefore contain the positive-part continuous component plus the binary
# structural-zero indicator.
test_that("retained spike uses positive-part plus binary reference design", {
  set.seed(9101)
  n <- 200L
  exposure <- c(rep(0, 60L), stats::rgamma(140L, shape = 2, rate = 1))
  z <- stats::rnorm(n)
  zero_indicator <- as.integer(exposure == 0)
  y <- 1.5 * zero_indicator + 0.6 * exposure + 0.3 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- mfp2(
    x,
    y,
    spike_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  )

  reference <- stats::glm(
    y ~ exposure_positive + exposure_bin + z,
    data = data.frame(
      y = y,
      exposure_positive = pmax(exposure, 0),
      exposure_bin = zero_indicator,
      z = z
    ),
    family = stats::gaussian()
  )

  expect_true(fit$spike[["exposure"]])
  expect_true(fit$catzero[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end test for a spike request that fails eligibility with no explicit
# zero or catzero fallback requested by the user.

# Once the spike is reset, the full linear reference must revert to the ordinary
# one-column linear term. No structural-zero binary component may remain.
test_that("rejected spike reverts full reference to ordinary linear term", {
  set.seed(9102)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- 0 # 1% structural-zero component: below default 10%
  z <- stats::rnorm(n)
  y <- 0.5 * exposure + 0.2 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  ))

  reference <- stats::glm(
    y ~ exposure + z,
    data = data.frame(y = y, exposure = exposure, z = z),
    family = stats::gaussian()
  )

  expect_false(fit$spike[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
  expect_false(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end test for a rejected spike with an explicit zero fallback.

# If the user requested zero handling independently, spike reset must preserve
# that request. The full reference then uses x+ but no binary indicator.
test_that("rejected spike preserves explicit zero reference design", {
  set.seed(9103)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- 0 # too few exact-zero observations for SAZ
  z <- stats::rnorm(n)
  y <- 0.5 * pmax(exposure, 0) + 0.2 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    zero_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  ))

  reference <- stats::glm(
    y ~ exposure_positive + z,
    data = data.frame(
      y = y,
      exposure_positive = pmax(exposure, 0),
      z = z
    ),
    family = stats::gaussian()
  )

  expect_false(fit$spike[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end test for a rejected spike with an explicit catzero fallback.

# Explicit catzero survives spike reset and still implies zero. The full
# reference must therefore retain both x+ and the binary structural-zero term.
test_that("rejected spike preserves explicit catzero binary reference design", {
  set.seed(9104)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- 0 # too few exact-zero observations for SAZ
  z <- stats::rnorm(n)
  exposure_bin <- as.integer(exposure == 0)
  y <- 1.5 * exposure_bin + 0.5 * pmax(exposure, 0) +
    0.2 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    catzero_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  ))

  reference <- stats::glm(
    y ~ exposure_positive + exposure_bin + z,
    data = data.frame(
      y = y,
      exposure_positive = pmax(exposure, 0),
      exposure_bin = exposure_bin,
      z = z
    ),
    family = stats::gaussian()
  )

  expect_false(fit$spike[["exposure"]])
  expect_true(fit$catzero[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end tests for spike-at-zero interaction with generated <term>_bin names.

# A retained spike implies catzero, so it must generate the same <term>_bin
# reference column and therefore reject an existing column with that name.
test_that("retained spike rejects an existing source-term _bin column", {
  set.seed(9105)
  n <- 200L
  exposure <- c(rep(0, 60L), stats::rgamma(140L, shape = 2, rate = 1))
  existing_bin <- rep(c(0, 1), length.out = n)
  x <- cbind(
    exposure = exposure,
    exposure_bin = existing_bin
  )
  y <- 1.2 * as.integer(exposure == 0) + 0.4 * exposure +
    0.2 * existing_bin + stats::rnorm(n)

  expect_error(
    mfp2(
      x,
      y,
      spike_vars = "exposure",
      xorder = "original",
      verbose = FALSE
    ),
    "Generated catzero indicator name.*exposure_bin"
  )
})

# If spike eligibility fails and the user did not request catzero, no binary
# reference column should be generated. An existing exposure_bin predictor must
# therefore remain valid and must not cause the catzero collision error.
test_that("rejected spike does not reserve a source-term _bin name", {
  set.seed(9106)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- 0 # below the default SAZ component threshold
  existing_bin <- rep(c(0, 1), length.out = n)
  x <- cbind(
    exposure = exposure,
    exposure_bin = existing_bin
  )
  y <- 0.4 * exposure + 0.2 * existing_bin + stats::rnorm(n)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    xorder = "original",
    verbose = FALSE
  ))

  expect_s3_class(fit, "mfp2")
  expect_false(fit$spike[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
})


# Migrated coverage from the former test_mfp2.R

# AIC/BIC SAZ selection must search the positive-only and
# positive-plus-binary functional forms independently and compare them in one
# joint table. This mocked FP3 example makes the positive-only FP3 powers win,
# while the full branch has different best powers, so inheriting Stage-1
# powers would fail the test.
test_that("joint SAZ IC selection independently optimizes both FP branches", {
  metric_names <- c(
    "logl", "df", "deviance_rs", "deviance_gaussian",
    "aic", "bic", "df_resid"
  )

  make_metrics <- function(model_names, aic, df) {
    out <- cbind(
      logl = -aic / 2,
      df = df,
      deviance_rs = aic,
      deviance_gaussian = aic,
      aic = aic,
      bic = aic + df,
      df_resid = 100 - df
    )
    colnames(out) <- metric_names
    rownames(out) <- model_names
    out
  }

  testthat::local_mocked_bindings(
    select_ic = function(..., xi, spike) {
      with_binary <- isTRUE(spike[[xi]])
      model_names <- if (with_binary) {
        c("null", "linear + Binary", "FP1 + Binary",
          "FP2 + Binary", "FP3 + Binary")
      } else {
        c("null", "linear", "FP1", "FP2", "FP3")
      }
      aic <- if (with_binary) c(120, 105, 98, 90, 94) else
        c(120, 108, 99, 92, 80)
      powers <- rbind(
        c(NA, NA, NA), c(1, NA, NA), c(2, NA, NA),
        c(-1, 2, NA), c(-1, 0.5, 2)
      )
      rownames(powers) <- model_names
      model_best <- which.min(aic)

      list(
        keep = FALSE, acd = FALSE, powers = powers,
        power_best = powers[model_best, , drop = FALSE],
        metrics = make_metrics(model_names, aic, c(1, 3, 4, 6, 8)),
        model_best = model_best, statistic = NA, pvalue = NA,
        spike = with_binary,
        current_adj_params = list(
          x = list(data_adj = NULL, data_xi = matrix(1, nrow = 8))
        ),
        transform_cache = list(branch = if (with_binary) "BP" else "P")
      )
    },
    fit_model = function(...) structure(list(), class = "mock_fit"),
    calculate_model_metrics = function(...) {
      c(
        logl = -55, df = 2, deviance_rs = 110,
        deviance_gaussian = 110, aic = 110, bic = 112,
        df_resid = 98
      )
    },
    .package = "mfp2"
  )

  out <- select_saz_ic(
    x = matrix(seq_len(8), ncol = 1, dimnames = list(NULL, "x")),
    xi = "x", keep = character(0), degree = 3, acdx = c(x = FALSE),
    y = seq_len(8), powers_current = list(x = c(1, 1, 1)),
    powers = list(x = c(-2, -1, 0, 0.5, 1, 2, 3)),
    criterion = "aic", ftest = FALSE, select = 0.05, alpha = 0.05,
    family = stats::gaussian(), family_string = "gaussian",
    zero = c(x = TRUE),
    catzero = list(x = matrix(rep(c(1, 0), 4), ncol = 1)),
    spike = list(x = TRUE), spike_decision = c(x = 1),
    acd_parameter = list(x = NULL), prev_adj_params = list(x = NULL),
    transform_cache = NULL, force_max_fp = c(x = FALSE),
    has_offset = FALSE, n_obs = 8, term_to_columns = list(x = "x")
  )

  expect_identical(nrow(out$metrics), 2L * 3L + 4L)
  expect_identical(
    rownames(out$metrics),
    c("null", "Binary", "linear", "linear + Binary", "FP1",
      "FP1 + Binary", "FP2", "FP2 + Binary", "FP3", "FP3 + Binary")
  )
  expect_identical(rownames(out$metrics)[out$model_best], "FP3")
  expect_equal(unname(out$power_best), c(-1, 0.5, 2))
  expect_identical(unname(out$spike_decision[["x"]]), 2L)
  expect_identical(out$selection_mode, "joint_ic")
})

# Exact AIC/BIC ties must be resolved by an explicit simplicity hierarchy,
# not by whichever candidate happens to occupy the first row. These are the
# two ordinary-FP pairs that can have equal adjusted SAZ df: binary-only versus
# positive-only linear (1 df), and both-components linear versus positive-only
# FP1 (2 df).
test_that("joint SAZ IC same-df ties prefer the simpler positive form", {
  tie_pairs <- list(
    one_df = list(
      candidates = c("Binary", "linear"),
      df = c(1, 1),
      expected = "Binary"
    ),
    two_df = list(
      candidates = c("linear + Binary", "FP1"),
      df = c(2, 2),
      expected = "linear + Binary"
    )
  )

  for (criterion in c("aic", "bic")) {
    for (tie_pair in tie_pairs) {
      # Check both possible row orders. A row-order tie-break would fail one
      # of these two arrangements for each pair.
      for (candidate_order in list(
        seq_along(tie_pair$candidates),
        rev(seq_along(tie_pair$candidates))
      )) {
        candidate_names <- tie_pair$candidates[candidate_order]
        candidate_df <- tie_pair$df[candidate_order]
        metrics <- cbind(
          df = candidate_df,
          aic = rep(10, length(candidate_names)),
          bic = rep(10, length(candidate_names))
        )
        rownames(metrics) <- candidate_names

        selected <- select_saz_ic_winner(
          metrics = metrics,
          criterion = criterion,
          eligible = seq_len(nrow(metrics)),
          acd = FALSE
        )

        expect_identical(rownames(metrics)[selected], tie_pair$expected)
      }
    }
  }
})

# Adjusted model df remains the first simplicity tie-break. This test uses a
# pair with the same positive-form complexity so that component count cannot
# obscure whether the lower-df rule was applied.
test_that("joint SAZ IC ties first prefer smaller adjusted df", {
  metrics <- cbind(
    df = c(2, 1),
    aic = c(10, 10),
    bic = c(10, 10)
  )
  rownames(metrics) <- c("linear + Binary", "linear")

  for (criterion in c("aic", "bic")) {
    selected <- select_saz_ic_winner(
      metrics = metrics,
      criterion = criterion,
      eligible = seq_len(nrow(metrics)),
      acd = FALSE
    )

    expect_identical(rownames(metrics)[selected], "linear")
  }
})

# The joint SAZ selector also accepts the package's ACD functional-form family.
# Its explicit simplicity metadata must recognize every established ACD label
# with and without the binary component.
test_that("joint SAZ IC simplicity hierarchy recognizes ACD forms", {
  candidate_names <- c(
    "null", "Binary",
    "linear", "linear + Binary",
    "linear(., A(x))", "linear(., A(x)) + Binary",
    "FP1(x, .)", "FP1(x, .) + Binary",
    "FP1(., A(x))", "FP1(., A(x)) + Binary",
    "FP1(x, A(x))", "FP1(x, A(x)) + Binary"
  )

  simplicity <- saz_ic_candidate_simplicity(candidate_names, acd = TRUE)

  expect_identical(
    simplicity$positive_complexity,
    c(0L, 0L, rep(1:5, each = 2L))
  )
  expect_identical(
    simplicity$component_count,
    c(0L, 1L, rep(c(1L, 2L), 5L))
  )
})

# Verbose labels are part of the behavioral distinction: IC selection must not
# imply that the selected powers came from a preceding SAZ Stage 1.
test_that("verbose SAZ output distinguishes joint IC from p-value stages", {
  set.seed(9001)
  x <- c(rep(0, 70), seq(0.2, 7, length.out = 210))
  z <- as.numeric(x == 0)
  y <- 3 * z + 1.5 * x + stats::rnorm(length(x), sd = 0.05)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))

  out_aic <- capture.output(
    suppressMessages(mfp2(
      xmat, y, spike_vars = "x", df = 1, criterion = "aic",
      shift = 0, scale = 1, center = FALSE, verbose = TRUE
    ))
  )
  expect_true(any(grepl("Joint Spike at Zero AIC Selection", out_aic,
                        fixed = TRUE)))
  expect_false(any(grepl("Stage 2 of Spike at Zero", out_aic, fixed = TRUE)))

  out_pvalue <- capture.output(
    suppressMessages(mfp2(
      xmat, y, spike_vars = "x", df = 1, criterion = "pvalue",
      select = 1, alpha = 1, shift = 0, scale = 1, center = FALSE,
      verbose = TRUE
    ))
  )
  expect_true(any(grepl("Stage 1 of Spike at Zero", out_pvalue, fixed = TRUE)))
  expect_true(any(grepl("Stage 2 of Spike at Zero", out_pvalue, fixed = TRUE)))
})

# =============================================================================
# 6. SAZ (spike-at-zero) — eligibility, cascade, and reset
# =============================================================================

# Test purpose: Fits a clear spike-at-zero example and checks that the spike
# flag is retained.
test_that("spike-at-zero basic fitting works", {
  set.seed(123)
  n <- 300
  prop_zero <- 0.25
  x_val <- numeric(n)
  zero_idx <- sample(seq_len(n), size = round(n * prop_zero))
  x_val[zero_idx] <- 0
  x_val[-zero_idx] <- rgamma(n - length(zero_idx), shape = 2, rate = 0.5)

  Z <- ifelse(x_val == 0, 1, 0)
  y_val <- 1.5 * Z + 2 * log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(x_mat, y_val, spike_vars = "exposure", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  # spike should be TRUE because the proportion of zeros is well within threshold
  expect_true(fit$fp_terms["exposure", "spike"])
})


# Test purpose: Checks that SAZ is reset when the zero component fails the
# minimum proportion threshold.
test_that("spike-at-zero is reset when zero proportion is too low", {
  # Almost no zeros
  set.seed(42)
  n <- 200
  x_val <- rgamma(n, shape = 2, rate = 1)
  # Only 2 zeros out of 200 = 1% which is below default 10%
  x_val[1:2] <- 0
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(x_mat, y_val, spike_vars = "exposure", verbose = FALSE),
    "spike"
  )

  # Spike should be reset
  expect_false(fit$fp_terms["exposure", "spike"])
})


# Test purpose: Checks that explicit zero handling survives after an ineligible
# spike request is reset.
test_that("spike cascade restores user-specified zero/catzero on reset", {
  # Build a scenario where spike is reset but user also set zero=TRUE
  set.seed(42)
  n <- 200
  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros => spike reset
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat, y_val,
      spike_vars = "exposure",
      zero_vars = "exposure",
      verbose = FALSE
    )
  })

  # Spike should be reset, but zero should be preserved
  expect_false(fit$fp_terms["exposure", "spike"])
  expect_true(fit$zero["exposure"])
})


# Test purpose: Checks that spike-at-zero handling can be requested inside fp()
# in the formula interface.
test_that("spike formula interface fp(spike = TRUE) works", {
  set.seed(123)
  n <- 300
  prop_zero <- 0.25
  x_val <- numeric(n)
  zero_idx <- sample(seq_len(n), size = round(n * prop_zero))
  x_val[zero_idx] <- 0
  x_val[-zero_idx] <- rgamma(n - length(zero_idx), shape = 2, rate = 0.5)

  Z <- ifelse(x_val == 0, 1, 0)
  y_val <- 1.5 * Z + 2 * log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  dat <- data.frame(y = y_val, exposure = x_val)

  fit <- mfp2(
    y ~ fp(exposure, spike = TRUE, center = FALSE),
    data = dat, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
})


# Test purpose: Checks that changing the SAZ component threshold changes spike
#  eligibility as expected.
test_that("min_saz_prop controls eligibility threshold", {
  set.seed(123)
  n <- 200
  x_val <- rgamma(n, shape = 2, rate = 1)
  # 15% zeros
  x_val[sample(n, 30)] <- 0
  y_val <- 0.5 * x_val + rnorm(n)
  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  # With default threshold 0.10 it should be eligible
  fit_low <- mfp2(x_mat, y_val, spike_vars = "exposure",
                  min_saz_prop = 0.10, verbose = FALSE)
  expect_true(fit_low$fp_terms["exposure", "spike"])

  # With high threshold 0.40 it should be ineligible
  expect_warning(
    fit_high <- mfp2(x_mat, y_val, spike_vars = "exposure",
                     min_saz_prop = 0.40, verbose = FALSE),
    "spike"
  )
  expect_false(fit_high$fp_terms["exposure", "spike"])
})


# Test purpose: Ensures SAZ eligibility requires enough positive-component
# observations, not only enough zero-component observations.
test_that("spike-at-zero is reset when positive component proportion is too low", {
  set.seed(104)
  n <- 200

  x_val <- numeric(n)
  x_val[1:5] <- rgamma(5, shape = 2, rate = 1)  # 2.5% positive component
  y_val <- 1.5 * (x_val == 0) + 0.2 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      verbose = FALSE
    ),
    "spike"
  )

  expect_false(fit$fp_terms["exposure", "spike"])
})


# Test purpose: Checks that negative values are rejected rather than counted in
# the zero component for SAZ eligibility.
test_that("spike-at-zero rejects negative fitting values", {
  set.seed(105)
  n <- 200

  x_val <- numeric(n)
  zero_idx <- sample(seq_len(n), size = 50)
  x_val[zero_idx] <- -1
  x_val[-zero_idx] <- rgamma(n - length(zero_idx), shape = 2, rate = 1)

  y_val <- 1.5 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)
  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_error(
    mfp2(x_mat, y_val, spike_vars = "exposure", verbose = FALSE),
    paste0(
      "require nonnegative covariates(.|[[:space:]])*exposure",
      "(.|[[:space:]])*Recode negative values explicitly"
    )
  )
})


# Test purpose: Ensures that when spike is reset, user-specified catzero handling
# is preserved rather than removed with the spike-implied cascade.
test_that("spike cascade preserves user-specified catzero on reset", {
  set.seed(106)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0  # too few zeros for SAZ eligibility
  y_val <- 1.5 * (x_val == 0) + 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      catzero_vars = "exposure",
      verbose = FALSE
    )
  })

  expect_false(fit$fp_terms["exposure", "spike"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})


# Test purpose: Directly checks reset_spike() for the all-zero case, where
# the positive component is absent.
test_that("reset_spike() resets all-zero variables", {
  x <- matrix(0, nrow = 100, ncol = 1, dimnames = list(NULL, "exposure"))

  spike <- c(exposure = TRUE)
  user_catzero <- c(exposure = FALSE)
  user_zero <- c(exposure = FALSE)

  expect_warning(
    out <- reset_spike(
      x = x,
      spike = spike,
      user_catzero = user_catzero,
      user_zero = user_zero,
      min_saz_prop = 0.10
    ),
    "positive observation proportion"
  )

  expect_false(out$spike["exposure"])
  expect_false(out$catzero["exposure"])
  expect_false(out$zero["exposure"])
})


# Test purpose: Checks that resolve_saz_eligibility() uses exact zeros for the
# structural-zero component without requiring a recoded matrix.
test_that("resolve_saz_eligibility() counts exact zeros", {
  # Only the zero-component proportion matters here. Use deterministic
  # positive values so this validation test is independent of RNG state.
  x <- matrix(
    c(rep(0, 20), seq(0.1, 18, length.out = 180)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  spike <- c(exposure = TRUE)
  catzero <- c(exposure = FALSE)
  zero <- c(exposure = FALSE)

  out <- resolve_saz_eligibility(
    x = x,
    spike = spike,
    catzero = catzero,
    zero = zero,
    min_saz_prop = 0.10
  )

  expect_true(out$spike["exposure"])
  expect_true(out$catzero["exposure"])
  expect_true(out$zero["exposure"])
})


# Test purpose: Exact-zero predicates drive SAZ binary detection and
# eligibility without any negative-to-zero recoding.
test_that("reset_spike exact-zero predicates drive eligibility", {
  x_raw <- cbind(
    exposure = c(rep(0, 4), rep(2, 16)),
    eligible = c(rep(0, 4), rep(1:4, each = 4))
  )

  spike <- c(exposure = TRUE, eligible = TRUE)
  user_catzero <- c(exposure = FALSE, eligible = FALSE)
  user_zero <- c(exposure = FALSE, eligible = FALSE)

  out <- suppressWarnings(reset_spike(
    x = x_raw,
    spike = spike,
    user_catzero = user_catzero,
    user_zero = user_zero,
    min_saz_prop = 0.10
  ))
  # exposure has one effective zero level plus one positive level, so it is
  # binary and must be reset.
  expect_false(out$spike[["exposure"]])
  expect_true(out$spike[["eligible"]])
})


# Test purpose: Checks that structural-zero proportions use finite observations
# and count exact-zero values as the SAZ zero component.
test_that("calculate_saz_prop_zero() reports retained SAZ proportions", {
  x <- cbind(
    exposure = c(0, 0, 1, 2, 3, 4, NA_real_, Inf),
    ordinary = seq_len(8)
  )

  out <- calculate_saz_prop_zero(
    x = x,
    spike = c(exposure = TRUE, ordinary = FALSE)
  )

  expect_equal(out[["exposure"]], 2 / 6)
  expect_true(is.na(out[["ordinary"]]))
})


# Test purpose: Ensures that when spike-only handling is reset, the variable
# returns to ordinary FP handling with no zero or catzero flags.
test_that("spike-only reset restores ordinary FP handling", {
  set.seed(305)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros for SAZ eligibility
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      verbose = FALSE
    ),
    "spike-at-zero option has been reset"
  )

  expect_false(fit$spike["exposure"])
  expect_false(fit$catzero["exposure"])
  expect_false(fit$zero["exposure"])
})


# Test purpose: Ensures that when spike is reset, explicitly requested zero
# handling is preserved rather than removed with the spike-implied cascade.
test_that("spike reset preserves user-specified zero handling", {
  set.seed(306)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros for SAZ eligibility
  y_val <- 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      zero_vars = "exposure",
      verbose = FALSE
    )
  })

  expect_false(fit$spike["exposure"])
  expect_false(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})


# Test purpose: Ensures that when spike is reset, explicitly requested catzero
# handling is preserved and still implies zero handling.
test_that("spike reset preserves user-specified catzero handling", {
  set.seed(307)
  n <- 200

  x_val <- rgamma(n, shape = 2, rate = 1)
  x_val[1:2] <- 0 # too few zeros for SAZ eligibility
  y_val <- 1.5 * (x_val == 0) + 0.5 * x_val + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  suppressWarnings({
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      catzero_vars = "exposure",
      verbose = FALSE
    )
  })

  expect_false(fit$spike["exposure"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})


# Test purpose: Ensures retained SAZ variables satisfy the internal cascade:
# spike implies catzero, and catzero implies zero.
test_that("retained spike variable implies catzero and zero handling", {
  set.seed(308)
  n <- 200

  x_val <- c(rep(0, 60), rgamma(140, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  expect_true(fit$spike["exposure"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
  expect_true(fit$fp_terms["exposure", "spike"])
})


# Test purpose: Ensures formula-interface fp(spike = TRUE) is translated into
# spike, catzero, and zero handling for an eligible SAZ variable.
test_that("formula interface fp(spike = TRUE) activates SAZ handling", {
  set.seed(309)
  n <- 200

  exposure <- c(rep(0, 50), rgamma(150, shape = 2, rate = 1))
  dat <- data.frame(
    y = 2 * (exposure == 0) + log(ifelse(exposure > 0, exposure, 1)) + rnorm(n),
    exposure = exposure,
    z = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ fp(exposure, spike = TRUE) + fp(z),
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})


# Test purpose: Checks that the SAZ algorithm runs under AIC-based selection,
# not only under p-value based closed testing.
test_that("spike-at-zero works with AIC criterion", {
  set.seed(310)
  n <- 220

  x_val <- c(rep(0, 70), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    criterion = "aic",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_true(fit$fp_terms["exposure", "spike"])
  expect_true(fit$spike_dec["exposure"] %in% c(1L, 2L, 3L))
})


# Test purpose: Checks that the SAZ algorithm runs under BIC-based selection.
test_that("spike-at-zero works with BIC criterion", {
  set.seed(311)
  n <- 220

  x_val <- c(rep(0, 70), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    criterion = "bic",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])
  expect_true(fit$fp_terms["exposure", "spike"])
  expect_true(fit$spike_dec["exposure"] %in% c(1L, 2L, 3L))
})


# Test purpose: Ensures retained SAZ variables with at most 3 distinct positive
# values have their maximum FP df forced to 1.
test_that("cap_spike_df() forces df = 1 for at most 3 distinct positive values", {
  x <- matrix(
    c(rep(0, 30), rep(c(1, 2, 3), each = 10)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  df <- c(exposure = 4)
  spike <- c(exposure = TRUE)

  expect_warning(
    out <- cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    "maximum FP df was reduced"
  )

  expect_equal(out[["exposure"]], 1)
})


# Test purpose: Ensures retained SAZ variables with 4 or 5 distinct positive
# values have their maximum FP df capped at FP1, i.e. df = 2.
test_that("cap_spike_df() caps df at 2 for 4 or 5 distinct positive values", {
  x <- matrix(
    c(rep(0, 30), rep(c(1, 2, 3, 4, 5), each = 10)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  df <- c(exposure = 4)
  spike <- c(exposure = TRUE)

  expect_warning(
    out <- cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    "maximum FP df was reduced"
  )

  expect_equal(out[["exposure"]], 2)
})


# Test purpose: Ensures retained SAZ variables with at least 6 distinct positive
# values keep the requested maximum FP df.
test_that("cap_spike_df() keeps df unchanged for at least 6 distinct positive values", {
  x <- matrix(
    c(rep(0, 30), rep(1:6, each = 10)),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  df <- c(exposure = 4)
  spike <- c(exposure = TRUE)

  expect_warning(
    out <- cap_spike_df(
      x = x,
      df = df,
      spike = spike
    ),
    NA
  )

  expect_equal(out[["exposure"]], 4)
})


# Test purpose: GLM SAZ stage 2 must assemble the final intercept-inclusive
# Model 2/Model 3 design once and tell fit_model() not to prepend another
# intercept. This protects the one-allocation GLM path.
test_that("SAZ reduced GLMs pass final intercept-inclusive designs", {
  seen <- new.env(parent = emptyenv())
  seen$x <- list()
  seen$x_has_intercept <- logical()

  testthat::local_mocked_bindings(
    fit_model = function(x, x_has_intercept = FALSE, ...) {
      seen$x[[length(seen$x) + 1L]] <- x
      seen$x_has_intercept <- c(seen$x_has_intercept, x_has_intercept)
      list(logl = -1, df = NCOL(x))
    },
    .package = "mfp2"
  )

  data_xi <- cbind(
    catzero = c(1, 1, 0, 0),
    fp1 = c(0, 0, 1.2, 2.4)
  )
  adjustment <- cbind(z = c(2, 3, 4, 5))
  stage1 <- list(
    current_adj_params = list(
      exposure = list(data_xi = data_xi, data_adj = adjustment)
    )
  )

  out <- fit_saz_reduced_models(
    stage1_selection = stage1,
    xi = "exposure",
    y = rep(0, 4),
    weights = NULL,
    offset = NULL,
    family = stats::gaussian(),
    family_string = "gaussian",
    method = NULL,
    strata = NULL,
    nocenter = NULL,
    control = NULL,
    rownames = NULL,
    has_offset = FALSE
  )

  expected_fit2 <- cbind(
    "(Intercept)" = rep(1, 4),
    fp1 = data_xi[, "fp1"],
    adjustment
  )
  expected_fit3 <- cbind(
    "(Intercept)" = rep(1, 4),
    catzero = data_xi[, "catzero"],
    adjustment
  )

  expect_length(seen$x, 2L)
  expect_equal(seen$x[[1L]], expected_fit2)
  expect_equal(seen$x[[2L]], expected_fit3)
  expect_identical(seen$x_has_intercept, c(TRUE, TRUE))
  expect_false("x" %in% names(out))
  expect_identical(out$data_xi, data_xi)
  expect_identical(out$adjustment_matrix, adjustment)
})


# Test purpose: Cox SAZ stage 2 must preserve the historical no-intercept
# design contract. The GLM allocation optimization must never add an ordinary
# intercept to a matrix sent to the Cox fitter.
test_that("SAZ reduced Cox models remain intercept-free", {
  seen <- new.env(parent = emptyenv())
  seen$x <- list()
  seen$x_has_intercept <- logical()

  testthat::local_mocked_bindings(
    fit_model = function(x, x_has_intercept = FALSE, ...) {
      seen$x[[length(seen$x) + 1L]] <- x
      seen$x_has_intercept <- c(seen$x_has_intercept, x_has_intercept)
      list(logl = -1, df = NCOL(x))
    },
    .package = "mfp2"
  )

  data_xi <- cbind(
    catzero = c(1, 1, 0, 0),
    fp1 = c(0, 0, 1.2, 2.4)
  )
  adjustment <- cbind(z = c(2, 3, 4, 5))
  stage1 <- list(
    current_adj_params = list(
      exposure = list(data_xi = data_xi, data_adj = adjustment)
    )
  )

  fit_saz_reduced_models(
    stage1_selection = stage1,
    xi = "exposure",
    y = rep(0, 4),
    weights = NULL,
    offset = NULL,
    family = NULL,
    family_string = "cox",
    method = "efron",
    strata = NULL,
    nocenter = NULL,
    control = NULL,
    rownames = NULL,
    has_offset = FALSE
  )

  expect_length(seen$x, 2L)
  expect_equal(seen$x[[1L]], cbind(fp1 = data_xi[, "fp1"], adjustment))
  expect_equal(seen$x[[2L]], cbind(catzero = data_xi[, "catzero"], adjustment))
  expect_identical(seen$x_has_intercept, c(FALSE, FALSE))
  expect_false("(Intercept)" %in% colnames(seen$x[[1L]]))
  expect_false("(Intercept)" %in% colnames(seen$x[[2L]]))
})


# Test purpose: The optimized GLM SAZ assembly must be numerically identical
# to the historical path where fit_glm() prepended the intercept itself.
test_that("SAZ reduced GLM fits match historical assembly numerically", {
  data_xi <- cbind(
    catzero = c(1, 1, 0, 0, 0, 0, 0, 0),
    fp1 = c(0, 0, 0.3, 0.8, 1.2, 1.7, 2.1, 2.8)
  )
  adjustment <- cbind(z = c(-1.2, -0.4, 0.1, 0.7, 1.1, 1.8, 2.3, 3.0))
  y <- c(0.4, 0.8, 1.2, 1.7, 2.0, 2.6, 3.0, 3.5)
  family <- stats::gaussian()
  stage1 <- list(
    current_adj_params = list(
      exposure = list(data_xi = data_xi, data_adj = adjustment)
    )
  )

  out <- fit_saz_reduced_models(
    stage1_selection = stage1,
    xi = "exposure",
    y = y,
    weights = NULL,
    offset = NULL,
    family = family,
    family_string = "gaussian",
    method = NULL,
    strata = NULL,
    nocenter = NULL,
    control = NULL,
    rownames = NULL,
    has_offset = FALSE
  )

  old_x2 <- cbind(fp1 = data_xi[, "fp1"], adjustment)
  old_x3 <- cbind(catzero = data_xi[, "catzero"], adjustment)

  expected2 <- fit_model(
    x = old_x2,
    y = y,
    family = family,
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    method = NULL,
    strata = NULL,
    control = NULL,
    rownames = NULL,
    nocenter = NULL,
    has_offset = FALSE
  )
  expected3 <- fit_model(
    x = old_x3,
    y = y,
    family = family,
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    method = NULL,
    strata = NULL,
    control = NULL,
    rownames = NULL,
    nocenter = NULL,
    has_offset = FALSE
  )

  expect_equal(out$fit2$coefficients, expected2$coefficients, tolerance = 1e-12)
  expect_equal(out$fit3$coefficients, expected3$coefficients, tolerance = 1e-12)
  expect_equal(out$fit2$logl, expected2$logl, tolerance = 1e-12)
  expect_equal(out$fit3$logl, expected3$logl, tolerance = 1e-12)
  expect_identical(out$fit2$df, expected2$df)
  expect_identical(out$fit3$df, expected3$df)
})


# Test purpose: Ensures public mfp2() applies the SAZ positive-part df cap
# before final FP model selection metadata are stored.
test_that("mfp2() applies positive-part df cap for retained spike variables", {
  set.seed(312)
  n <- 180

  positive_values <- rep(c(1, 2, 3), each = 40)
  x_val <- c(rep(0, 60), positive_values)
  y_val <- 1.5 * (x_val == 0) + 0.4 * x_val + rnorm(n, sd = 0.2)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = "exposure",
      df = 4,
      verbose = FALSE
    ),
    "maximum FP df was reduced"
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$spike["exposure"])

  # The positive-part cap reduces the continuous FP complexity to one df, but
  # an eligible SAZ term enters selection with its one-df zero indicator as
  # well. Keep these two quantities distinct: df_setting describes the capped
  # continuous component, while df_initial is the total initial SAZ term df.
  expect_equal(fit$fp_terms["exposure", "df_setting"], 1)
  expect_equal(fit$fp_terms["exposure", "df_initial"], 2)
})


# Test purpose: Ensures prediction works for retained SAZ models, including
# new zero and positive values in newdata.
test_that("predict.mfp2() works for retained spike-at-zero models", {
  set.seed(313)
  n <- 200

  x_val <- c(rep(0, 50), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n, sd = 0.3)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  newx <- matrix(
    c(0, 0, 0.5, 1, 2, 4),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  pred <- predict(fit, newdata = newx)

  expect_length(pred, nrow(newx))
  expect_true(all(is.finite(pred)))
})


# Test purpose: Ensures prediction for retained SAZ models rejects negative
# newdata instead of treating it as the zero component.
test_that("predict.mfp2() rejects negative newdata for SAZ models", {
  set.seed(314)
  n <- 200

  x_val <- c(rep(0, 50), rgamma(150, shape = 2, rate = 1))
  y_val <- 2 * (x_val == 0) + log(ifelse(x_val > 0, x_val, 1)) + rnorm(n, sd = 0.3)

  x_mat <- matrix(x_val, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x_mat,
    y_val,
    spike_vars = "exposure",
    verbose = FALSE
  )

  newx <- matrix(
    c(-2, -1, 0, 0.5, 2),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  expect_error(
    predict(fit, newdata = newx),
    "require nonnegative covariates(.|[[:space:]])*exposure"
  )
})


# Test purpose: Ensures SAZ eligibility is resolved independently for multiple
# spike variables, so one reset variable does not reset all spike variables.
test_that("multiple spike variables are reset independently", {
  set.seed(315)
  n <- 240

  eligible <- c(rep(0, 60), rgamma(180, shape = 2, rate = 1))
  ineligible <- rgamma(n, shape = 2, rate = 1)
  ineligible[1:2] <- 0 # too few zeros

  y_val <- 1.5 * (eligible == 0) + log(ifelse(eligible > 0, eligible, 1)) +
    0.2 * ineligible + rnorm(n)

  x_mat <- cbind(
    eligible = eligible,
    ineligible = ineligible
  )

  expect_warning(
    fit <- mfp2(
      x_mat,
      y_val,
      spike_vars = c("eligible", "ineligible"),
      verbose = FALSE
    ),
    "ineligible"
  )

  expect_true(fit$spike["eligible"])
  expect_false(fit$spike["ineligible"])

  expect_true(fit$catzero["eligible"])
  expect_true(fit$zero["eligible"])

  expect_false(fit$catzero["ineligible"])
  expect_false(fit$zero["ineligible"])
})


# Test purpose: Guards the MFP endpoint convention for SAZ Stage 2. With
# alpha = 1, an exact p-value of 1 must not simplify the two-component SAZ
# representation; simplification occurs only when p is strictly greater than
# alpha.
test_that("SAZ Stage 2 retains both components at alpha = 1 and p = 1", {
  metric <- c(
    logl = 0,
    df = 1,
    aic = 0,
    bic = 0,
    deviance_gaussian = 1,
    df_resid = 10
  )

  testthat::local_mocked_bindings(
    calculate_lr_test = function(...) {
      list(statistic = 0, pvalue = 1)
    },
    .package = "mfp2"
  )

  out <- compute_saz_stage2_decision(
    metrics = list(
      metrics1 = metric,
      metrics2 = metric,
      metrics3 = metric
    ),
    criterion = "pvalue",
    alpha = 1,
    n_obs = 10,
    ftest = FALSE
  )

  expect_equal(out$decision, saz_decision_codes[["cont_binary"]])
  expect_equal(unname(out$pvalue), c(1, 1))
})


# Test purpose: When both component-removal tests are non-significant, the two
# reduced models are non-nested and may have different df. Stage 2 must compare
# their already-computed BIC values rather than favoring the model with the
# larger raw likelihood. An exact BIC tie is resolved in favor of binary-only.
test_that("SAZ Stage 2 uses BIC for two non-significant reductions", {
  metric <- c(
    logl = 0,
    df = 1,
    aic = 0,
    bic = 0,
    deviance_gaussian = 1,
    df_resid = 10
  )

  testthat::local_mocked_bindings(
    calculate_lr_test = function(...) {
      list(statistic = 0, pvalue = 0.50)
    },
    .package = "mfp2"
  )

  cases <- list(
    continuous_smaller_bic = list(
      bic2 = 10,
      bic3 = 20,
      # Deliberately make Model 2's likelihood worse so this case would fail
      # under the former larger-log-likelihood tie-break.
      logl2 = -100,
      logl3 = -1,
      expected = saz_decision_codes[["continuous_only"]]
    ),
    binary_smaller_bic = list(
      bic2 = 20,
      bic3 = 10,
      # Deliberately make Model 2's likelihood better for the same reason.
      logl2 = -1,
      logl3 = -100,
      expected = saz_decision_codes[["binary_only"]]
    ),
    exact_bic_tie = list(
      bic2 = 10,
      bic3 = 10,
      logl2 = -1,
      logl3 = -100,
      expected = saz_decision_codes[["binary_only"]]
    )
  )

  for (case in cases) {
    metrics2 <- metric
    metrics3 <- metric
    metrics2[["bic"]] <- case$bic2
    metrics3[["bic"]] <- case$bic3
    metrics2[["logl"]] <- case$logl2
    metrics3[["logl"]] <- case$logl3

    out <- compute_saz_stage2_decision(
      metrics = list(
        metrics1 = metric,
        metrics2 = metrics2,
        metrics3 = metrics3
      ),
      criterion = "pvalue",
      alpha = 0.05,
      n_obs = 10,
      ftest = FALSE
    )

    expect_identical(out$decision, case$expected)
    expect_equal(unname(out$pvalue), c(0.50, 0.50))
  }
})
