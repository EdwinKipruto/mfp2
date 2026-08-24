# Shared shift, scale, and center fixtures and expectations.

# =============================================================================
# 4.1 Shift and scale safety contract - shared setup
# =============================================================================
# Review objective
# ----------------
# The tests in Sections 4.2 and 4.3 verify the complete user-input contract for
# shift and scale without changing the existing preprocessing calculations:
#
#   1. NULL keeps automatic shift/scale selection.
#   2. A single finite numeric value is recycled to every matrix column.
#   3. A multi-value vector must be fully and uniquely named.
#   4. Named values are matched to colnames(x), not to supplied position.
#   5. Missing, unknown, duplicate, or empty names are rejected.
#   6. shift values must be finite numeric values with no missing values.
#   7. scale values must additionally be strictly positive.
#   8. Formula-level scalar and per-variable fp() settings remain compatible.
#   9. Correctly specified models retain identical fits and predictions when
#      the same named settings are supplied in a different order.
#
# Sections 4.2 and 4.3 are intentionally separate. MFPI has an additional
# grouping-variable column that must participate in matrix-name validation but
# is subsequently handled as categorical metadata by the fitting procedure.

# -----------------------------------------------------------------------------
# 4.1.1 Deterministic shared data
# -----------------------------------------------------------------------------
# Expected use:
# - setting_x and setting_y exercise the ordinary mfp2 matrix interface.
# - mfpi_setting_x adds the MFPI grouping variable "svi" as a matrix column.
# - The response contains nonlinear signal so stored transformations and
#   predictions are meaningful integration checks rather than validation-only
#   calls.

set.seed(20260717)
setting_n <- 96L
setting_x <- cbind(
  age = runif(setting_n, 20, 80),
  weight = runif(setting_n, 45, 105)
)
setting_y <-
  1.5 * sqrt(setting_x[, "age"] + 2) -
  2.0 * log(setting_x[, "weight"] + 3) +
  rnorm(setting_n, sd = 0.15)

mfpi_setting_x <- cbind(
  svi = rep(c(0, 1), each = setting_n / 2),
  age = setting_x[, "age"],
  weight = setting_x[, "weight"]
)
mfpi_setting_y <-
  setting_y +
  0.4 * mfpi_setting_x[, "svi"] * sqrt(mfpi_setting_x[, "age"] + 2)

# -----------------------------------------------------------------------------
# 4.1.2 Shared fitting helpers
# -----------------------------------------------------------------------------
# Expected use:
# - Keep model specifications identical across ordered and reverse-ordered
#   setting vectors.
# - Limit differences between compared fits to the order in which named shift
#   or scale values were supplied.

fit_mfp2_settings <- function(shift = NULL, scale = NULL, center = FALSE) {
  mfp2(
    x = setting_x,
    y = setting_y,
    shift = shift,
    scale = scale,
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = colnames(setting_x),
    center = center,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )
}

fit_mfpi_settings <- function(shift = NULL, scale = NULL, center = FALSE) {
  mfpi(
    x = mfpi_setting_x,
    y = mfpi_setting_y,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    shift = shift,
    scale = scale,
    df = 2,
    select = 1,
    alpha = 1,
    force_max_fp_vars = c("age", "weight"),
    center = center,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )
}

# -----------------------------------------------------------------------------
# 4.1.3 Shared validation-call helpers
# -----------------------------------------------------------------------------
# Expected use:
# - Reach the public matrix interfaces with the smallest practical fit.
# - Keep validation tests focused on shift/scale errors rather than FP model
#   selection behavior.

mfp2_validation_call <- function(argument, value) {
  args <- list(
    x = setting_x,
    y = setting_y,
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 1,
    verbose = FALSE
  )
  args[[argument]] <- value
  do.call(mfp2, args)
}

mfpi_validation_call <- function(argument, value) {
  args <- list(
    x = mfpi_setting_x,
    y = mfpi_setting_y,
    group_var = "svi",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    cycles = 1,
    verbose = FALSE
  )
  args[[argument]] <- value
  do.call(mfpi, args)
}

# -----------------------------------------------------------------------------
# 4.1.4 Shared malformed-input fixtures
# -----------------------------------------------------------------------------
# Each fixture represents one distinct part of the public input contract. The
# same cases are applied separately to mfp2.default() and mfpi.default().

invalid_named_settings <- function(column_names) {
  n <- length(column_names)
  values <- seq_len(n)

  partially_named <- stats::setNames(values, column_names)
  names(partially_named)[2L] <- ""

  duplicated_names <- stats::setNames(values, column_names)
  names(duplicated_names)[2L] <- column_names[[1L]]

  unknown_replacement <- stats::setNames(values, column_names)
  names(unknown_replacement)[n] <- "replacement_column"

  list(
    # More than one value with no names must never be matched by position.
    unnamed = unname(values),

    # Every element of a multi-value vector must have a non-empty name.
    partially_named = partially_named,

    # Ambiguous duplicate names cannot identify one value per matrix column.
    duplicated = duplicated_names,

    # Replacing a valid name introduces an unknown supplied column. The
    # omitted known column itself is valid under the partial-vector contract.
    unknown_replacement = unknown_replacement,

    # A complete required set plus an extra name must be rejected.
    unknown_extra = stats::setNames(
      c(values, n + 1L),
      c(column_names, "unknown_column")
    ),

    # Extra named values fail because their names are not matrix columns.
    extra_unknown_names = stats::setNames(
      c(values, n + 1L, n + 2L),
      c(column_names, "unknown_column_1", "unknown_column_2")
    ),

    # Explicit numerical settings must be finite and nonmissing.
    non_finite = stats::setNames(
      replace(as.numeric(values), 1L, Inf),
      column_names
    ),
    missing_value = stats::setNames(
      replace(as.numeric(values), 1L, NA_real_),
      column_names
    ),

    # Logical values must not be silently coerced into numeric settings.
    logical = stats::setNames(rep(TRUE, n), column_names)
  )
}

# -----------------------------------------------------------------------------
# 4.1.5 Shared expectation helpers for malformed names
# -----------------------------------------------------------------------------
# These expectations are deliberately explicit. A review can therefore map
# every required naming failure to one assertion while still applying the same
# contract consistently to the two public matrix interfaces.

expect_setting_name_errors <- function(call_setting, column_names) {
  invalid <- invalid_named_settings(column_names)

  # Shift and scale accept named subsets, but never unnamed multi-value
  # vectors, incomplete names, duplicate names, or names outside colnames(x).
  expect_error(
    call_setting("shift", invalid$unnamed),
    "`shift` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("shift", invalid$partially_named),
    "`shift` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("shift", invalid$duplicated),
    "`shift` names must be unique"
  )
  for (invalid_shift in invalid[c(
    "unknown_replacement", "unknown_extra", "extra_unknown_names"
  )]) {
    expect_error(
      call_setting("shift", invalid_shift),
      "`shift` contains unknown column name"
    )
  }

  expect_error(
    call_setting("scale", invalid$unnamed),
    "`scale` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("scale", invalid$partially_named),
    "`scale` must be a single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    call_setting("scale", invalid$duplicated),
    "`scale` names must be unique"
  )
  for (invalid_scale in invalid[c(
    "unknown_replacement", "unknown_extra", "extra_unknown_names"
  )]) {
    expect_error(
      call_setting("scale", invalid_scale),
      "`scale` contains unknown column name"
    )
  }
}

# -----------------------------------------------------------------------------
# 4.1.6 Shared expectation helpers for invalid values
# -----------------------------------------------------------------------------
# Numerical validation is common to shift and scale. Positivity is tested
# separately because it is an additional requirement that applies only to scale.

expect_setting_numeric_errors <- function(call_setting, column_names) {
  invalid <- invalid_named_settings(column_names)

  for (argument in c("shift", "scale")) {
    # Requirement: Inf and -Inf are not valid explicit settings.
    expect_error(
      call_setting(argument, invalid$non_finite),
      paste0("`", argument, "` must contain only finite values")
    )

    # Requirement: direct matrix settings cannot contain NA values.
    expect_error(
      call_setting(argument, invalid$missing_value),
      paste0("`", argument, "` must not contain missing values")
    )

    # Requirement: TRUE/FALSE must not be accepted as 1/0 settings.
    expect_error(
      call_setting(argument, invalid$logical),
      paste0("`", argument, "` must contain numeric values")
    )
  }
}

expect_nonpositive_scale_errors <- function(call_setting, column_names) {
  zero_scale <- stats::setNames(rep(1, length(column_names)), column_names)
  zero_scale[[1L]] <- 0

  negative_scale <- stats::setNames(rep(1, length(column_names)), column_names)
  negative_scale[[1L]] <- -1

  # Requirement: zero is not a valid divisor for scaling, including when
  # supplied through the named-partial interface.
  expect_error(
    call_setting("scale", zero_scale),
    "`scale` must contain only strictly positive values"
  )
  expect_error(
    call_setting("scale", stats::setNames(0, column_names[[1L]])),
    "`scale` must contain only strictly positive values"
  )

  # Requirement: negative scales are not valid transformation settings.
  expect_error(
    call_setting("scale", negative_scale),
    "`scale` must contain only strictly positive values"
  )
  expect_error(
    call_setting("scale", stats::setNames(-1, column_names[[1L]])),
    "`scale` must contain only strictly positive values"
  )
}


# -----------------------------------------------------------------------------
