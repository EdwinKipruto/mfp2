library(testthat)
library(mfp2)


test_that("create_z_variables rejects a constant continuous covariate", {
  cont_var <- matrix(
    rep(5, 6),
    ncol = 1L,
    dimnames = list(NULL, "age")
  )
  group_var <- matrix(
    rep(c(0, 1), each = 3L),
    ncol = 1L,
    dimnames = list(NULL, "treatment")
  )

  # Centering would turn the constant FP basis into exact zero columns.
  expect_error(
    mfp2:::create_z_variables(
      cont_var = cont_var,
      group_var = group_var,
      power = 1,
      center = TRUE
    ),
    "must contain at least two distinct values",
    fixed = TRUE
  )

  # Without centering, each constant group-specific column is a scaled group
  # indicator and remains unidentifiable. The guard is therefore unconditional.
  expect_error(
    mfp2:::create_z_variables(
      cont_var = cont_var,
      group_var = group_var,
      power = 1,
      center = FALSE
    ),
    "must contain at least two distinct values",
    fixed = TRUE
  )
})


test_that("create_z_variables retains nonconstant continuous input", {
  cont_var <- matrix(
    c(1, 2, 3, 2, 3, 4),
    ncol = 1L,
    dimnames = list(NULL, "age")
  )
  group_var <- matrix(
    rep(c(0, 1), each = 3L),
    ncol = 1L,
    dimnames = list(NULL, "treatment")
  )

  result <- expect_silent(
    mfp2:::create_z_variables(
      cont_var = cont_var,
      group_var = group_var,
      power = 1,
      center = TRUE
    )
  )

  expect_equal(nrow(result$z), nrow(cont_var))
  expect_false(anyNA(result$z))
  expect_true(all(is.finite(result$z)))
})
