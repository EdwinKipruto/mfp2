# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# -----------------------------------------------------------------------------
# 4.1.7-4.1.8 Direct normalization regression tests
# -----------------------------------------------------------------------------
# These tests isolate argument normalization from model fitting. They protect
# the precise distinction that caused the original bug: a named length-one
# vector is partial and variable-specific, while an unnamed length-one vector
# is global.

test_that("4.1.7 named scalar settings are partial, not global", {
  columns <- c("cavol", "age", "weight")
  normalized_shift <- normalize_named_numeric_setting(
    value = c(age = 20),
    column_names = columns,
    argument_name = "shift",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = FALSE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )
  normalized_scale <- normalize_named_numeric_setting(
    value = c(age = 10),
    column_names = columns,
    argument_name = "scale",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = TRUE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )

  expect_identical(
    normalized_shift,
    c(cavol = NA_real_, age = 20, weight = NA_real_)
  )
  expect_identical(
    normalized_scale,
    c(cavol = NA_real_, age = 10, weight = NA_real_)
  )
})


test_that("4.1.8 unnamed scalar settings remain global", {
  columns <- c("cavol", "age", "weight")
  normalized_shift <- normalize_named_numeric_setting(
    value = 20,
    column_names = columns,
    argument_name = "shift",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = FALSE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )
  normalized_scale <- normalize_named_numeric_setting(
    value = 10,
    column_names = columns,
    argument_name = "scale",
    allow_null = TRUE,
    scalar_recycle = TRUE,
    strictly_positive = TRUE,
    allow_na = FALSE,
    allow_partial_named = TRUE
  )

  expect_identical(normalized_shift, c(cavol = 20, age = 20, weight = 20))
  expect_identical(normalized_scale, c(cavol = 10, age = 10, weight = 10))
})
