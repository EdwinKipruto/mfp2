# MFPI-specific regression tests for the conventional FP1 candidate set.
# These tests intentionally leave ordinary mfp2() closed-test behaviour intact,
# where the fixed linear model is handled separately from nonlinear FP1 search.

library(testthat)

# The old MFPI dispatcher removed p = 1 before every FP1 search. Restricting the
# MFPI candidate set to p = 1 therefore provides a sharp regression test: all
# four flexibility variants must now fit and retain p = 1 as a valid FP1 result.
test_that("mfpi FP1 retains power 1 for all flexibility variants", {
  data("prostate", package = "mfp2")

  for (fl in c("flex1", "flex2", "flex3", "flex4")) {
    fit <- mfpi(
      lpsa ~ fp(age) + svi + fp(cavol, df = 1),
      data = prostate,
      group_var = "svi",
      cont_vars = "cavol",
      cont_var_forms = c(cavol = "fp1"),
      powers = list(cavol = 1),
      flex = fl,
      verbose = FALSE
    )

    expect_s3_class(fit, "mfpi")

    expect_equal(
      fit$cont_var_forms[["cavol"]],
      "fp1",
      info = paste("flex =", fl)
    )

    main_power <- fit$all_model_metrics$fp_powers_main[[1L]]

    expect_equal(
      as.numeric(main_power),
      1,
      info = paste("flex =", fl)
    )
  }
})

# Guard the package-wide contract: ordinary mfp2() keeps its established
# closed-test convention. A df > 1 search space containing only p = 1 is still
# invalid there because the fixed linear model is fitted separately.
test_that("ordinary mfp2 FP1 closed-test handling is unchanged", {
  x <- matrix(
    seq(0.5, 8, length.out = 100),
    ncol = 1,
    dimnames = list(NULL, "x")
  )
  y <- 1 + 2 * x[, 1]

  expect_error(
    mfp2(
      x,
      y,
      powers = list(x = 1),
      df = 2,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    "non-linear|candidate|power",
    ignore.case = TRUE
  )
})

# Retaining the linear FP1 member means preserving p = 1 when it is supplied;
# it must not silently enlarge a user-restricted candidate set that omits p = 1.
test_that("mfpi FP1 does not add power 1 to a restricted candidate set", {
  data("prostate", package = "mfp2")

  allowed <- c(-2, -1, -0.5, 0, 0.5, 2, 3)
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "fp1"),
    powers = list(cavol = allowed),
    flex = "flex1",
    verbose = FALSE
  )

  selected <- as.numeric(fit$all_model_metrics$fp_powers_main[[1L]])
  expect_true(selected %in% allowed)
  expect_false(1 %in% selected)
})
