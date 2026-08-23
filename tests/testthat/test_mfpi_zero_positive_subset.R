library(testthat)
library(mfp2)


test_that("zero-handled cont_vars require positive rows after subsetting", {
  n <- 24L
  x <- cbind(
    treatment = rep(c(0, 1), length.out = n),
    age = c(-12:-1, 1:12)
  )
  y <- seq_len(n) / 10
  non_positive_rows <- seq_len(12L)

  # The complete data contain positive age values, but the fitting subset does
  # not. Validation must occur after resolving subset and before dispatching to
  # either flexibility-specific fitting engine.
  for (flex_method in c("flex2", "flex4")) {
    expect_error(
      mfpi(
        x = x,
        y = y,
        group_var = "treatment",
        cont_vars = "age",
        zero_vars = "age",
        subset = non_positive_rows,
        flex = flex_method,
        family = "gaussian",
        verbose = FALSE
      ),
      paste0(
        "Zero-handled variables in `cont_vars` must contain at least one ",
        "positive value.*Problematic variables: age"
      )
    )
  }
})


test_that("zero-handled cont_vars are also validated without subset", {
  n <- 12L
  x <- cbind(
    treatment = rep(c(0, 1), length.out = n),
    age = -seq_len(n)
  )
  y <- seq_len(n) / 10

  expect_error(
    mfpi(
      x = x,
      y = y,
      group_var = "treatment",
      cont_vars = "age",
      zero_vars = "age",
      flex = "flex2",
      family = "gaussian",
      verbose = FALSE
    ),
    "Problematic variables: age",
    fixed = TRUE
  )
})
