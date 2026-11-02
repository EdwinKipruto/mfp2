# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 4.2 mfp2 shift and scale behavior
# =============================================================================
# Scope
# -----
# This section tests the ordinary mfp2 matrix interface first, followed by the
# formula interface. Each test title is numbered so a reviewer can trace the
# expected behavior directly to the corresponding requirement above.

# -----------------------------------------------------------------------------
# 4.2.1 NULL keeps automatic preprocessing
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that leaving both arguments as NULL still delegates shift and scale
#   selection to the package's existing automatic preprocessing logic.
# - Protect against the new normalization helper accidentally treating NULL as
#   an invalid or incomplete user-supplied vector.
# - Confirm that automatic preprocessing produces usable settings for every
#   matrix column before fractional-polynomial fitting begins.
#
# Expected output:
# - fit$transformations contains rows named "age" and "weight".
# - The selected shift and scale cells contain no NA values.
# - Every selected shift and scale value is finite.
# - Both scale values are greater than zero.

test_that("4.2.1 mfp2.default() keeps automatic shift and scale for NULL", {
  fit <- fit_mfp2_settings()
  transformations <- fit$transformations[
    colnames(setting_x), c("shift", "scale"), drop = FALSE
  ]

  expect_false(anyNA(transformations))
  expect_true(all(is.finite(as.matrix(transformations))))
  expect_true(all(transformations$scale > 0))
})


# -----------------------------------------------------------------------------
# 4.2.2 Scalars are recycled to every matrix column
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that a single finite numeric shift remains a global matrix setting.
# - Verify that a single positive numeric scale remains a global matrix setting.
# - Protect the documented scalar interface while tightening validation for
#   vectors containing more than one value.
#
# Expected output:
# - In colnames(setting_x) order, the stored shifts are c(2, 2).
# - In colnames(setting_x) order, the stored scales are c(10, 10).
# - No name matching is required from the user for these scalar inputs.

test_that("4.2.2 mfp2.default() recycles scalar shift and scale", {
  fit <- fit_mfp2_settings(shift = 2, scale = 10)
  transformations <- fit$transformations[
    colnames(setting_x), c("shift", "scale"), drop = FALSE
  ]

  expect_equal(unname(transformations$shift), rep(2, ncol(setting_x)))
  expect_equal(unname(transformations$scale), rep(10, ncol(setting_x)))
})


# -----------------------------------------------------------------------------
# 4.2.2a Named partial settings fix supplied columns and estimate the rest
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that c(age = 20) and c(age = 100) apply only to age rather than
#   being recycled as global shift and scale settings.
# - Verify that an unspecified nonlinear predictor retains internal NA
#   sentinels until the ordinary shift and scale estimation steps.
# - Exercise the public matrix interface through stored transformation metadata.

test_that("4.2.2a mfp2.default() estimates unspecified named partial settings", {
  index <- seq_len(60L)
  partial_x <- cbind(
    age = seq(20, 79),
    weight = -30 + ((17 * index) %% 60)
  )
  expected_weight_shift <- find_shift_factor(partial_x[, "weight"])
  expected_weight_scale <- find_scale_factor(
    partial_x[, "weight"] + expected_weight_shift
  )
  partial_y <-
    1.2 * sqrt(partial_x[, "age"] + 20) -
    0.7 * log(partial_x[, "weight"] + expected_weight_shift) +
    0.05 * sin(index / 4)

  fit <- mfp2(
    x = partial_x,
    y = partial_y,
    shift = c(age = 20),
    scale = c(age = 100),
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = colnames(partial_x),
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  transformations <- fit$transformations[
    colnames(partial_x), c("shift", "scale"), drop = FALSE
  ]
  expect_equal(unname(transformations["age", "shift"]), 20)
  expect_equal(
    unname(transformations["weight", "shift"]),
    expected_weight_shift
  )
  expect_equal(unname(transformations["age", "scale"]), 100)
  expect_equal(
    unname(transformations["weight", "scale"]),
    expected_weight_scale
  )
})


# -----------------------------------------------------------------------------
# 4.2.2b Insufficient explicit shifts fail before FP fitting
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that a supplied named shift is not silently increased.
# - Verify that the existing positivity check reports the affected variable when
#   the shifted nonlinear predictor still contains zero or negative values.

test_that("4.2.2b mfp2.default() rejects insufficient named partial shifts", {
  index <- seq_len(60L)
  insufficient_x <- cbind(
    age = -30 + ((17 * index) %% 60),
    weight = 40 + ((23 * index) %% 60)
  )
  insufficient_y <-
    0.3 * insufficient_x[, "age"] +
    0.2 * sqrt(insufficient_x[, "weight"]) +
    sin(index / 5)

  expect_error(
    mfp2(
      x = insufficient_x,
      y = insufficient_y,
      shift = c(age = 30),
      scale = 1,
      df = 2,
      select = 1,
      alpha = 1,
      force_max_fp_vars = colnames(insufficient_x),
      center = FALSE,
      cycles = 5,
      xorder = "original",
      verbose = FALSE
    ),
    "Problematic variables: age"
  )
})


# -----------------------------------------------------------------------------
# 4.2.2c Partial scale keeps mapped-term automatic defaults
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that an unspecified grouped design block still receives scale = 1.
# - Verify that an explicit partial scale remains authoritative for its column.

test_that("4.2.2c partial scale preserves mapped-term defaults", {
  normalized_scale <- c(group_b = NA_real_, group_c = NA_real_, age = 100)
  term_to_columns <- list(
    group = c("group_b", "group_c"),
    age = "age"
  )

  result <- expand_scale_for_mapped_terms(
    scale = normalized_scale,
    vnames = names(normalized_scale),
    term_to_columns = term_to_columns,
    automatic = is.na(normalized_scale)
  )

  expect_equal(
    result,
    c(group_b = 1, group_c = 1, age = 100)
  )
})


# -----------------------------------------------------------------------------
# 4.2.3 Ordered named vectors are stored against the matching columns
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify the accepted multi-value form: one numeric value for every matrix
#   column, with complete and unique names.
# - Confirm that normalization returns the settings in colnames(setting_x)
#   order and that downstream transformation metadata preserves those values.
# - Establish the reference fit used conceptually by the reordered-input test.
#
# Expected output:
# - The "age" transformation row stores shift = 1 and scale = 10.
# - The "weight" transformation row stores shift = 2 and scale = 100.
# - Reading the two rows in matrix-column order returns shifts c(1, 2) and
#   scales c(10, 100).

test_that("4.2.3 mfp2.default() accepts fully named shift and scale", {
  fit <- fit_mfp2_settings(
    shift = c(age = 1, weight = 2),
    scale = c(age = 10, weight = 100)
  )
  transformations <- fit$transformations[
    colnames(setting_x), c("shift", "scale"), drop = FALSE
  ]

  expect_equal(unname(transformations$shift), c(1, 2))
  expect_equal(unname(transformations$scale), c(10, 100))
})


# -----------------------------------------------------------------------------
# 4.2.4 Named-vector order does not change the fitted model
# -----------------------------------------------------------------------------
# Test purpose:
# - Reproduce the original failure mode by supplying the same named settings in
#   matrix order and then in reverse order.
# - Verify that values are matched by their names rather than silently relabeled
#   according to their supplied positions.
# - Exercise the complete path from normalization through FP transformation,
#   model fitting, stored metadata, and prediction.
#
# Expected output:
# - The ordered and reversed fits have identical transformation tables.
# - Their coefficient vectors and fitted values are numerically identical.
# - Predictions for the first 12 rows of setting_x are numerically identical.
# - In particular, age always receives shift = 1 and scale = 10, while weight
#   always receives shift = 2 and scale = 100, regardless of input order.

test_that("4.2.4 mfp2.default() matches shift and scale by name", {
  ordered <- fit_mfp2_settings(
    shift = c(age = 1, weight = 2),
    scale = c(age = 10, weight = 100)
  )
  reversed <- fit_mfp2_settings(
    shift = c(weight = 2, age = 1),
    scale = c(weight = 100, age = 10)
  )

  expect_equal(ordered$transformations, reversed$transformations)
  expect_equal(unname(stats::coef(ordered)), unname(stats::coef(reversed)))
  expect_equal(unname(stats::fitted(ordered)), unname(stats::fitted(reversed)))

  new_x <- setting_x[1:12, , drop = FALSE]
  expect_equal(
    unname(stats::predict(ordered, newdata = new_x)),
    unname(stats::predict(reversed, newdata = new_x))
  )
})


# -----------------------------------------------------------------------------
# 4.2.5 Invalid multi-value names are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that positional matching is no longer available for multi-value
#   matrix settings.
# - Verify that every accepted multi-value vector identifies only known columns,
#   with no ambiguity or surplus entries.
# - Apply the same naming contract independently to shift and scale.
#
# Expected output for each argument:
# - Unnamed multi-value vectors and vectors containing empty element names
#   error with the single-unnamed-or-named-vector message.
# - Duplicate names error with "`<argument>` names must be unique".
# - Named vectors containing unknown columns error with
#   "`<argument>` contains unknown column name".
# - No model is fitted for any invalid input.

test_that("4.2.5 mfp2.default() rejects invalid shift and scale names", {
  expect_setting_name_errors(mfp2_validation_call, colnames(setting_x))
})


# -----------------------------------------------------------------------------
# 4.2.6 Missing, non-finite, and logical values are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that direct matrix settings are genuine numeric configuration values,
#   not logical vectors that R could silently coerce to 0 and 1.
# - Verify that explicit values are complete and finite before transformation or
#   model fitting is attempted.
# - Apply the same numeric validation to shift and scale.
#
# Expected output for each argument:
# - A vector containing Inf errors with
#   "`<argument>` must contain only finite values".
# - A vector containing NA errors with
#   "`<argument>` must not contain missing values".
# - A logical vector errors with
#   "`<argument>` must contain numeric values".
# - No invalid value reaches the transformation table.

test_that("4.2.6 mfp2.default() rejects invalid numeric settings", {
  expect_setting_numeric_errors(mfp2_validation_call, colnames(setting_x))
})


# -----------------------------------------------------------------------------
# 4.2.7 Scale must be strictly positive
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify the scale-specific constraint that every explicit divisor is greater
#   than zero, in addition to being numeric and finite.
# - Protect transformation calculations from division by zero and from a sign
#   reversal introduced by a negative scale.
#
# Expected output:
# - A named scale vector containing 0 errors with
#   "`scale` must contain only strictly positive values".
# - A named scale vector containing -1 produces the same error.
# - Validation fails before an mfp2 model is fitted.

test_that("4.2.7 mfp2.default() rejects nonpositive scale values", {
  expect_nonpositive_scale_errors(mfp2_validation_call, colnames(setting_x))
})


# -----------------------------------------------------------------------------
# 4.2.8 Formula-level scalar settings remain global
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that the matrix-interface safety change does not impose a new naming
#   requirement on the established top-level formula scalar interface.
# - Verify that formula preprocessing expands one global scalar to each modeled
#   numeric variable before calling the default method.
#
# Expected output:
# - The fitted transformation rows for age and weight both store shift = 2.
# - The fitted transformation rows for age and weight both store scale = 10.
# - The formula fit completes without requiring named top-level vectors.
# - `center = FALSE` is repeated inside each fp() term deliberately: fp() has
#   its own default `center = TRUE`, so the top-level setting alone would not
#   disable centering for these FP terms. Centering is not under test here.

test_that("4.2.8 mfp2.formula() keeps scalar shift and scale compatible", {
  data("prostate", package = "mfp2")

  fit <- mfp2(
    lpsa ~ fp(age, df = 2, center = FALSE, force_max_fp = TRUE) +
      fp(weight, df = 2, center = FALSE, force_max_fp = TRUE),
    data = prostate,
    shift = 2,
    scale = 10,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    unname(fit$transformations[c("age", "weight"), "shift"]),
    c(2, 2)
  )
  expect_equal(
    unname(fit$transformations[c("age", "weight"), "scale"]),
    c(10, 10)
  )
})


# -----------------------------------------------------------------------------
# 4.2.9 Per-variable fp() settings remain supported
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that variable-specific shift and scale values supplied inside fp()
#   remain part of the supported formula interface.
# - Verify that formula preprocessing creates internally named vectors and that
#   the default method matches those values to the correct model-matrix columns.
# - Protect against a regression where formula-generated vectors become unnamed
#   and are rejected by the stricter matrix-interface validation.
#
# Expected output:
# - The age transformation row stores shift = 1 and scale = 10.
# - The weight transformation row stores shift = 2 and scale = 100.
# - Reading rows c("age", "weight") returns shifts c(1, 2) and scales
#   c(10, 100), with no positional reassignment.
# - `center = FALSE` is specified inside each fp() term so this test isolates
#   shift/scale behavior instead of inheriting fp()'s default center = TRUE.

test_that("4.2.9 mfp2.formula() keeps per-variable fp() settings compatible", {
  data("prostate", package = "mfp2")

  fit <- mfp2(
    lpsa ~ fp(
      age,
      df = 2,
      shift = 1,
      scale = 10,
      center = FALSE,
      force_max_fp = TRUE
    ) + fp(
      weight,
      df = 2,
      shift = 2,
      scale = 100,
      center = FALSE,
      force_max_fp = TRUE
    ),
    data = prostate,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    unname(fit$transformations[c("age", "weight"), "shift"]),
    c(1, 2)
  )
  expect_equal(
    unname(fit$transformations[c("age", "weight"), "scale"]),
    c(10, 100)
  )
})


# -----------------------------------------------------------------------------
# 4.2.10 Direct single-shift path matches the historical two-sweep path
# -----------------------------------------------------------------------------
# Test purpose:
# - Exercise the actual Step-16 shift optimization with df > 1. Variables with
#   df = 1 have their shift forced to zero earlier in mfp2.default(), so they do
#   not test this code path meaningfully.
# - Compare the optimized direct-matrix branch with the historical two-sweep
#   branch on the SAME matrix interface. Attaching the private preprocessing
#   matrix makes mfp2.default() take the historical two-matrix branch while
#   keeping the fitting matrix and model specification otherwise identical.
# - This is a direct numerical regression test of the source refactor, not a
#   comparison between matrix and formula coefficient parameterizations.

test_that("4.2.10 direct single-shift path matches historical two-sweep path", {
  n <- 96L
  index <- seq_len(n)
  x_direct <- cbind(
    x1 = seq(-8, 12, length.out = n),
    x2 = -4 + ((17 * index) %% 31)
  )
  x_before <- x_direct
  y <- 1.2 + 0.28 * x_direct[, "x1"] - 0.17 * x_direct[, "x2"] +
    0.015 * x_direct[, "x1"]^2 + 0.03 * sin(index / 5)

  fit_optimized <- mfp2(
    x = x_direct,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  # Force the historical Step-16 branch with a separate preprocessing source
  # containing exactly the same numerical matrix. extract_preprocess_matrix()
  # removes this private attribute before fitting.
  x_historical <- x_direct
  attr(x_historical, "mfp2_preprocess_x") <- x_direct

  fit_historical <- mfp2(
    x = x_historical,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    fit_optimized$transformations,
    fit_historical$transformations,
    tolerance = 0
  )
  expect_equal(
    unname(coef(fit_optimized)),
    unname(coef(fit_historical)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fitted(fit_optimized)),
    unname(fitted(fit_historical)),
    tolerance = 1e-12
  )
  expect_identical(x_direct, x_before)
})


# -----------------------------------------------------------------------------
# 4.2.11 Explicit shift/scale path matches the historical two-sweep path
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that sharing the shifted direct matrix does not alter explicit
#   preprocessing settings or fitted results.
# - Use non-collinear predictors and df > 1 so the supplied shift values are
#   actually active (df = 1 would force shift = 0 before Step 16).

test_that("4.2.11 direct single-shift path preserves explicit settings", {
  n <- 100L
  index <- seq_len(n)
  x_direct <- cbind(
    x1 = seq(-4, 9, length.out = n),
    x2 = 2 + ((19 * index) %% 37)
  )
  y <- 2 + 0.45 * x_direct[, "x1"] + 0.11 * x_direct[, "x2"] +
    0.02 * cos(index / 7)

  shift <- c(x1 = 5, x2 = 1)
  scale <- c(x1 = 10, x2 = 20)

  fit_optimized <- mfp2(
    x = x_direct,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = shift,
    scale = scale,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  x_historical <- x_direct
  attr(x_historical, "mfp2_preprocess_x") <- x_direct

  fit_historical <- mfp2(
    x = x_historical,
    y = y,
    df = 2,
    keep = colnames(x_direct),
    shift = shift,
    scale = scale,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    fit_optimized$transformations,
    fit_historical$transformations,
    tolerance = 0
  )
  expect_equal(
    unname(coef(fit_optimized)),
    unname(coef(fit_historical)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fitted(fit_optimized)),
    unname(fitted(fit_historical)),
    tolerance = 1e-12
  )
})


# -----------------------------------------------------------------------------
# 4.2.12 Direct subset still uses full-data preprocessing
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that the optimized direct branch still estimates shift and scale
#   before applying subset, exactly as the original source does.
# - Use a full-rank retained design; the previous fixture made the retained x1
#   and x2 columns affine functions of each other and therefore correctly
#   triggered mfp2's rank-deficiency validation before this behavior was tested.

test_that("4.2.12 direct subset preserves full-data shift and scale", {
  n <- 122L
  interior_index <- seq_len(n - 2L)
  x1 <- c(-50, seq(-3, 4, length.out = n - 2L), 80)
  x2 <- c(-20, 2 + ((17 * interior_index) %% 53), 70)
  x_direct <- cbind(x1 = x1, x2 = x2)
  y <- 0.7 + 0.22 * x1 - 0.08 * x2 + 0.02 * sin(seq_len(n) / 6)
  rows <- 2:(n - 1L)

  fit_optimized <- mfp2(
    x = x_direct,
    y = y,
    subset = rows,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  # Historical two-sweep reference on the same direct matrix and subset.
  x_historical <- x_direct
  attr(x_historical, "mfp2_preprocess_x") <- x_direct
  fit_historical <- mfp2(
    x = x_historical,
    y = y,
    subset = rows,
    df = 2,
    keep = colnames(x_direct),
    shift = NULL,
    scale = NULL,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expected_shift <- apply(x_direct, 2, find_shift_factor)
  expected_scale <- vapply(
    seq_len(ncol(x_direct)),
    function(j) find_scale_factor(x_direct[, j] + expected_shift[j]),
    numeric(1L)
  )
  names(expected_scale) <- colnames(x_direct)

  stored <- fit_optimized$transformations[
    colnames(x_direct), c("shift", "scale"), drop = FALSE
  ]

  expect_equal(unname(stored$shift), unname(expected_shift), tolerance = 0)
  expect_equal(unname(stored$scale), unname(expected_scale), tolerance = 0)
  expect_equal(
    fit_optimized$transformations,
    fit_historical$transformations,
    tolerance = 0
  )
  expect_equal(
    unname(fitted(fit_optimized)),
    unname(fitted(fit_historical)),
    tolerance = 1e-12
  )

  retained_shift <- apply(x_direct[rows, , drop = FALSE], 2, find_shift_factor)
  expect_false(isTRUE(all.equal(unname(expected_shift), unname(retained_shift))))
})


# =============================================================================
# End of tests
# =============================================================================

# =============================================================================
# Named center settings in default matrix interfaces
# =============================================================================

# Test purpose: A named center vector must be aligned by predictor name rather
# than by the order in which values are supplied.
test_that("mfp2.default() matches named center settings by colnames(x)", {
  fit <- fit_mfp2_settings(
    shift = 0,
    scale = 1,
    center = c(weight = FALSE, age = TRUE)
  )

  expect_true(isTRUE(fit$transformations["age", "center"]))
  expect_false(isTRUE(fit$transformations["weight", "center"]))
})


# Test purpose: A named scalar is a one-variable override, not a global scalar;
# omitted variables retain the public default center = TRUE.
test_that("mfp2.default() fills partial named center settings from the default", {
  fit <- fit_mfp2_settings(
    shift = 0,
    scale = 1,
    center = c(age = FALSE)
  )

  expect_false(isTRUE(fit$transformations["age", "center"]))
  expect_true(isTRUE(fit$transformations["weight", "center"]))
})


# Test purpose: Positional multi-value center vectors are ambiguous and must be
# rejected with guidance to use either a scalar or named values.
test_that("mfp2.default() rejects unnamed multi-value center settings", {
  expect_error(
    mfp2_validation_call("center", c(TRUE, FALSE)),
    "supply a scalar or name each value"
  )
})


# Test purpose: Named center settings use the same duplicate/unknown-name
# validation contract as df, shift, and scale.
test_that("mfp2.default() validates names in center settings", {
  expect_error(
    mfp2_validation_call("center", c(age = TRUE, unknown = FALSE)),
    "unknown column name.*unknown"
  )
  expect_error(
    mfp2_validation_call(
      "center",
      stats::setNames(c(TRUE, FALSE), c("age", "age"))
    ),
    "names must be unique.*age"
  )
})
