# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 4.3 mfpi shift and scale behavior
# =============================================================================
# MFPI-specific review notes
# --------------------------
# The MFPI matrix contains the grouping variable "svi" in addition to ordinary
# covariates. Named shift and scale vectors may specify any subset of its columns;
# unspecified entries remain automatic. Downstream MFPI code continues to treat
# svi as categorical grouping metadata, so its final shift and scale are 0 and 1.
#
# The tests below are intentionally separate from mfp2 tests because MFPI also
# builds an adjustment model, removes or specially handles the grouping column,
# computes interaction metrics, and has its own prediction method.

# -----------------------------------------------------------------------------
# 4.3.1 NULL keeps automatic preprocessing for MFPI covariates
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that NULL retains MFPI's existing automatic preprocessing for the
#   ordinary continuous covariates after the grouping column is identified.
# - Protect the separate MFPI normalization and adjustment-model path from
#   receiving unresolved or invalid settings.
# - Confirm that the new validation does not alter automatic shift/scale choice.
#
# Expected output:
# - fit$shift[c("age", "weight")] contains two non-missing finite values.
# - fit$scale[c("age", "weight")] contains two non-missing finite values.
# - Both continuous-variable scales are strictly greater than zero.
# - The fit completes with the existing categorical handling of svi unchanged.

test_that("4.3.1 mfpi.default() keeps automatic shift and scale for NULL", {
  fit <- fit_mfpi_settings()

  expect_false(anyNA(fit$shift[c("age", "weight")]))
  expect_false(anyNA(fit$scale[c("age", "weight")]))
  expect_true(all(is.finite(fit$shift[c("age", "weight")])))
  expect_true(all(is.finite(fit$scale[c("age", "weight")])))
  expect_true(all(fit$scale[c("age", "weight")] > 0))
})


# -----------------------------------------------------------------------------
# 4.3.2 Scalars are recycled to ordinary MFPI covariates
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that a scalar remains a valid global MFPI setting even though the
#   matrix also contains the grouping variable svi.
# - Confirm that scalar expansion survives grouping-column handling and reaches
#   both ordinary continuous covariates used by the adjustment model.
#
# Expected output:
# - fit$shift[c("age", "weight")] equals c(2, 2).
# - fit$scale[c("age", "weight")] equals c(10, 10).
# - The test intentionally leaves svi to MFPI's existing categorical metadata
#   rules rather than treating it as an ordinary transformed covariate.

test_that("4.3.2 mfpi.default() recycles scalar shift and scale", {
  fit <- fit_mfpi_settings(shift = 2, scale = 10)

  expect_equal(unname(fit$shift[c("age", "weight")]), c(2, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 10))
})


# -----------------------------------------------------------------------------
# 4.3.2a MFPI uses the same named partial-setting contract
# -----------------------------------------------------------------------------

test_that("4.3.2a mfpi.default() estimates unspecified named partial settings", {
  fit <- fit_mfpi_settings(shift = c(age = 2), scale = c(age = 100))

  expect_equal(unname(fit$shift["age"]), 2)
  expect_equal(
    unname(fit$shift["weight"]),
    find_shift_factor(mfpi_setting_x[, "weight"])
  )
  expect_equal(unname(fit$shift["svi"]), 0)
  expect_equal(unname(fit$scale["age"]), 100)
  expect_equal(
    unname(fit$scale["weight"]),
    find_scale_factor(
      mfpi_setting_x[, "weight"] + fit$shift[["weight"]]
    )
  )
  expect_equal(unname(fit$scale["svi"]), 1)
})


# Test purpose: MFPI must neutralize explicitly supplied preprocessing values
# for an ordinary binary adjustment covariate, not only for the group variable.
test_that("4.3.2b mfpi.default() resets supplied scale for binary covariates", {
  binary_adjustment <- rep(c(0, 1), length.out = setting_n)
  x <- cbind(mfpi_setting_x, binary_adjustment = binary_adjustment)
  y <- mfpi_setting_y + 0.3 * binary_adjustment

  fit <- mfpi(
    x = x,
    y = y,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    shift = c(binary_adjustment = 10),
    scale = c(binary_adjustment = 1000),
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = c("age", "weight"),
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(unname(fit$shift["binary_adjustment"]), 0)
  expect_equal(unname(fit$scale["binary_adjustment"]), 1)
})


# -----------------------------------------------------------------------------
# 4.3.3 Fully named vectors include the MFPI grouping column
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that complete named vectors remain accepted at the public matrix
#   boundary, including an explicit value for the grouping column svi.
# - Confirm that, after MFPI applies its special grouping-variable handling, the
#   values for age and weight remain aligned with their names.
#
# Expected output:
# - The complete inputs use names c("svi", "age", "weight").
# - fit$shift[c("age", "weight")] equals c(1, 2).
# - fit$scale[c("age", "weight")] equals c(10, 100).
# - The grouping column does not cause either continuous-variable value to move
#   to the wrong variable.

test_that("4.3.3 mfpi.default() accepts fully named shift and scale", {
  fit <- fit_mfpi_settings(
    shift = c(svi = 0, age = 1, weight = 2),
    scale = c(svi = 1, age = 10, weight = 100)
  )

  expect_equal(unname(fit$shift[c("age", "weight")]), c(1, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 100))
})


# -----------------------------------------------------------------------------
# 4.3.4 Named-vector order does not change MFPI results
# -----------------------------------------------------------------------------
# Test purpose:
# - Reproduce the positional-matching risk in the more complex MFPI path by
#   reversing named settings that include the grouping variable.
# - Verify alignment before and after MFPI handles or removes svi from ordinary
#   transformation processing.
# - Exercise normalized settings, adjustment-model transformations, fitted
#   coefficients, interaction metrics, and MFPI prediction in one regression
#   test.
#
# Expected output:
# - ordered$shift and reversed$shift are identical named vectors.
# - ordered$scale and reversed$scale are identical named vectors.
# - The two adjustment-model transformation tables and coefficient vectors are
#   numerically identical.
# - ordered$all_model_metrics and reversed$all_model_metrics are identical.
# - Link-scale predictions for the age term and all interaction models are
#   identical for the first 12 rows of mfpi_setting_x.
# - Age retains shift = 1 and scale = 10; weight retains shift = 2 and
#   scale = 100, regardless of the supplied order around svi.

test_that("4.3.4 mfpi.default() matches shift and scale by name", {
  ordered <- fit_mfpi_settings(
    shift = c(svi = 0, age = 1, weight = 2),
    scale = c(svi = 1, age = 10, weight = 100)
  )
  reversed <- fit_mfpi_settings(
    shift = c(weight = 2, age = 1, svi = 0),
    scale = c(weight = 100, age = 10, svi = 1)
  )

  expect_equal(ordered$shift, reversed$shift)
  expect_equal(ordered$scale, reversed$scale)
  expect_equal(
    ordered$adjustment_model$transformations,
    reversed$adjustment_model$transformations
  )

  # Stage 1 is fitted after MFPI has already shifted/scaled the working matrix.
  # Its returned mfp2 object must nevertheless restore the original MFPI shifts
  # so that the preprocessing table and direct prediction from raw newdata are
  # both expressed relative to the user's covariate scale.
  expect_equal(
    unname(ordered$adjustment_model$transformations[c("age", "weight"), "shift"]),
    c(1, 2)
  )
  expect_equal(
    unname(ordered$adjustment_model$transformations[c("age", "weight"), "scale"]),
    c(10, 100)
  )
  expect_equal(
    unname(stats::coef(ordered$adjustment_model)),
    unname(stats::coef(reversed$adjustment_model))
  )
  expect_equal(ordered$all_model_metrics, reversed$all_model_metrics)

  new_x <- as.data.frame(mfpi_setting_x[1:12, , drop = FALSE])
  pred_ordered <- predict(
    ordered,
    newdata = new_x,
    terms = "age",
    model = "all",
    type = "link"
  )
  pred_reversed <- predict(
    reversed,
    newdata = new_x,
    terms = "age",
    model = "all",
    type = "link"
  )

  expect_equal(
    pred_ordered$predictions$fit,
    pred_reversed$predictions$fit
  )

  adjustment_newdata <- as.data.frame(
    mfpi_setting_x[1:12, c("age", "weight"), drop = FALSE]
  )
  adjustment_prediction <- predict(
    ordered$adjustment_model,
    newdata = adjustment_newdata
  )
  expect_equal(
    unname(adjustment_prediction),
    unname(ordered$adjustment_model$linear.predictors[1:12]),
    tolerance = 1e-8
  )
})


# -----------------------------------------------------------------------------
# 4.3.5 Invalid MFPI multi-value names are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that MFPI rejects unnamed multi-value vectors, empty element names,
#   duplicate names, and names outside its public matrix columns.
# - Confirm that any supplied names are validated before special treatment of
#   the grouping variable begins.
# - Apply the same partial, unique naming contract to shift and scale.
#
# Expected output for each argument:
# - Unnamed multi-value vectors and empty element names produce the
#   single-unnamed-or-named-vector error.
# - Duplicate names produce the unique-names error.
# - Supplying an unknown name produces the unknown-column error; omitting
#   otherwise valid columns leaves those settings automatic.
# - No adjustment or interaction model is fitted for invalid inputs.

test_that("4.3.5 mfpi.default() rejects invalid shift and scale names", {
  expect_setting_name_errors(
    mfpi_validation_call,
    colnames(mfpi_setting_x)
  )
})


# -----------------------------------------------------------------------------
# 4.3.6 Missing, non-finite, and logical MFPI values are rejected
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that MFPI uses the same finite, nonmissing numeric-value contract as
#   the ordinary mfp2 matrix interface.
# - Ensure invalid values fail before settings are separated into grouping and
#   adjustment-model components.
#
# Expected output for both shift and scale:
# - Inf produces the finite-values error.
# - NA produces the no-missing-values error.
# - TRUE/FALSE vectors produce the numeric-values error.
# - No invalid value is stored in fit$shift, fit$scale, or the adjustment model.

test_that("4.3.6 mfpi.default() rejects invalid numeric settings", {
  expect_setting_numeric_errors(
    mfpi_validation_call,
    colnames(mfpi_setting_x)
  )
})


# -----------------------------------------------------------------------------
# 4.3.7 MFPI scale values must be strictly positive
# -----------------------------------------------------------------------------
# Test purpose:
# - Verify that MFPI enforces the scale-specific positivity rule on each explicit
#   named setting before grouping-variable handling or adjustment fitting.
# - Prevent zero divisors and negative rescaling from entering any MFPI model.
#
# Expected output:
# - A named scale vector containing 0 errors with
#   "`scale` must contain only strictly positive values".
# - A named scale vector containing -1 produces the same error.
# - Neither the adjustment model nor the interaction models are fitted.

test_that("4.3.7 mfpi.default() rejects nonpositive scale values", {
  expect_nonpositive_scale_errors(
    mfpi_validation_call,
    colnames(mfpi_setting_x)
  )
})


# -----------------------------------------------------------------------------
# 4.3.8 Formula-level scalar settings remain global in MFPI
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that the stricter MFPI matrix-vector contract does not alter the
#   documented top-level formula scalar interface.
# - Verify that formula preprocessing applies global scalar values to ordinary
#   continuous variables while continuing to handle svi as the group variable.
#
# Expected output:
# - fit$shift[c("age", "weight")] equals c(2, 2).
# - fit$scale[c("age", "weight")] equals c(10, 10).
# - The formula fit succeeds without a named top-level vector or an explicit
#   scalar setting for svi.

test_that("4.3.8 mfpi.formula() keeps scalar shift and scale compatible", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ svi + fp(age, df = 2, force_max_fp = TRUE) +
      fp(weight, df = 2, force_max_fp = TRUE),
    data = prostate,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    shift = 2,
    scale = 10,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(unname(fit$shift[c("age", "weight")]), c(2, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 10))
})


# -----------------------------------------------------------------------------
# 4.3.9 Per-variable fp() settings remain supported in MFPI formulas
# -----------------------------------------------------------------------------
# Test purpose:
# - Confirm that per-variable settings inside fp() remain supported when the
#   formula also contains the MFPI grouping variable.
# - Verify that formula preprocessing generates correctly named internal
#   settings and that MFPI's grouping-column handling does not shift them.
# - Protect compatibility for existing formula calls that intentionally assign
#   different preprocessing values to different continuous variables.
#
# Expected output:
# - fit$shift[c("age", "weight")] equals c(1, 2).
# - fit$scale[c("age", "weight")] equals c(10, 100).
# - Age receives only the values declared in fp(age, ...), and weight receives
#   only the values declared in fp(weight, ...).

test_that("4.3.9 mfpi.formula() keeps per-variable fp() settings compatible", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ svi +
      fp(
        age,
        df = 2,
        shift = 1,
        scale = 10,
        force_max_fp = TRUE
      ) + fp(
        weight,
        df = 2,
        shift = 2,
        scale = 100,
        force_max_fp = TRUE
      ),
    data = prostate,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 5,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(unname(fit$shift[c("age", "weight")]), c(1, 2))
  expect_equal(unname(fit$scale[c("age", "weight")]), c(10, 100))
})


# Test purpose: MFPI applies the same colnames(x)-based matching while its input
# matrix also contains the grouping column.
test_that("mfpi.default() matches named center settings by colnames(x)", {
  fit <- fit_mfpi_settings(
    shift = 0,
    scale = 1,
    center = c(weight = FALSE, age = TRUE, svi = FALSE)
  )

  expect_true(isTRUE(fit$adjustment_model$transformations["age", "center"]))
  expect_false(isTRUE(fit$adjustment_model$transformations["weight", "center"]))
})


# Test purpose: Partial named MFPI specifications retain center = TRUE for
# omitted continuous predictors.
test_that("mfpi.default() fills partial named center settings from the default", {
  fit <- fit_mfpi_settings(
    shift = 0,
    scale = 1,
    center = c(weight = FALSE)
  )

  expect_true(isTRUE(fit$adjustment_model$transformations["age", "center"]))
  expect_false(isTRUE(fit$adjustment_model$transformations["weight", "center"]))
})


# Test purpose: MFPI rejects positional center vectors and malformed names at
# its public default-method boundary.
test_that("mfpi.default() validates vector center settings", {
  expect_error(
    mfpi_validation_call("center", c(TRUE, FALSE, TRUE)),
    "supply a scalar or name each value"
  )
  expect_error(
    mfpi_validation_call("center", c(age = TRUE, unknown = FALSE)),
    "unknown column name.*unknown"
  )
  expect_error(
    mfpi_validation_call(
      "center",
      stats::setNames(c(TRUE, FALSE), c("age", "age"))
    ),
    "names must be unique.*age"
  )
})
