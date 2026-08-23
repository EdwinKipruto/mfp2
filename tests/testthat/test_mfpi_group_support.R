library(testthat)
library(mfp2)


make_mfpi_group_support_data <- function() {
  n <- 20L
  list(
    x = cbind(
      treatment = rep(c(0, 1), each = n / 2L),
      age = seq_len(n)
    ),
    y = seq_len(n) / 10
  )
}


test_that("MFPI rejects a singleton group after subsetting", {
  dat <- make_mfpi_group_support_data()

  expect_error(
    mfpi(
      x = dat$x,
      y = dat$y,
      group_var = "treatment",
      cont_vars = "age",
      subset = c(1:6, 11),
      flex = "flex2",
      family = "gaussian",
      verbose = FALSE
    ),
    "Singleton groups: 1",
    fixed = TRUE
  )
})


test_that("MFPI rejects two-row groups for FP2 interactions", {
  dat <- make_mfpi_group_support_data()

  expect_error(
    mfpi(
      x = dat$x,
      y = dat$y,
      group_var = "treatment",
      cont_vars = "age",
      cont_var_forms = c(age = "fp2"),
      subset = c(1:6, 11:12),
      flex = "flex4",
      family = "gaussian",
      verbose = FALSE
    ),
    "FP2 interactions require at least three observations.*Problematic groups: 1"
  )
})


test_that("group-centered zero handling requires positive group support", {
  x <- cbind(
    treatment = rep(c(0, 1), each = 6L),
    age = c(-4, -3, -2, -1, 0, 1, 1:6)
  )
  y <- seq_len(nrow(x)) / 10

  expect_error(
    mfpi(
      x = x,
      y = y,
      group_var = "treatment",
      cont_vars = "age",
      zero_vars = "age",
      center_type = "group",
      flex = "flex2",
      family = "gaussian",
      verbose = FALSE
    ),
    "age \\(group 0: 1 positive; need 2\\)"
  )
})


test_that("group-centered zero-handled FP2 requires three positive rows", {
  x <- cbind(
    treatment = rep(c(0, 1), each = 6L),
    age = c(-3, -2, -1, 0, 1, 2, 1:6)
  )
  y <- seq_len(nrow(x)) / 10

  expect_error(
    mfpi(
      x = x,
      y = y,
      group_var = "treatment",
      cont_vars = "age",
      cont_var_forms = c(age = "fp2"),
      zero_vars = "age",
      center_type = "group",
      flex = "flex4",
      family = "gaussian",
      verbose = FALSE
    ),
    "age \\(group 0: 2 positive; need 3\\)"
  )
})
