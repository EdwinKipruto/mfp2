# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Verifies that mfpi.default() uses the same partial named
# override semantics as mfp2.default() for its Stage-1 adjustment model.
test_that("mfpi.default() accepts partial named df select and alpha overrides", {
  set.seed(16001)
  n <- 180L
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), length.out = n)),
    x = stats::runif(n, 1, 8),
    z = stats::runif(n, 1, 5),
    w = stats::runif(n, 2, 7)
  )
  dat$y <- 1 + 0.4 * dat$x + 0.8 * dat$x * (dat$group == "treated") +
    0.5 * dat$z + stats::rnorm(n, sd = 0.3)

  fit <- mfpi(
    dat[, c("group", "x", "z", "w")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 5,
    df = c(z = 1),
    select = c(z = 1),
    alpha = c(w = 1),
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    p_interact = 1,
    verbose = FALSE
  )

  adjustment_terms <- fit$adjustment_model$fp_terms
  expect_equal(as.numeric(adjustment_terms["z", "df_initial"]), 1)
  expect_equal(as.numeric(adjustment_terms["w", "df_initial"]), 4)
  expect_equal(as.numeric(adjustment_terms["z", "select"]), 1)
  expect_equal(as.numeric(adjustment_terms["w", "select"]), 0.05)
  expect_equal(as.numeric(adjustment_terms["z", "alpha"]), 0.05)
  expect_equal(as.numeric(adjustment_terms["w", "alpha"]), 1)
})


# Test purpose: Ensures MFPI no longer accepts positional multi-value settings.
test_that("mfpi.default() rejects unnamed multi-value df select and alpha", {
  set.seed(16002)
  n <- 60L
  x <- data.frame(
    group = factor(rep(c("control", "treated"), length.out = n)),
    x = stats::runif(n, 1, 5),
    z = stats::runif(n, 1, 4)
  )
  y <- stats::rnorm(n)

  common <- list(
    x = x,
    y = y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    verbose = FALSE
  )

  expect_error(
    do.call(mfpi, c(common, list(df = c(1, 4, 4)))),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    do.call(mfpi, c(common, list(select = c(1, 0.05, 0.05)))),
    "single unnamed numeric value or a named numeric vector"
  )
  expect_error(
    do.call(mfpi, c(common, list(alpha = c(1, 0.05, 0.05)))),
    "single unnamed numeric value or a named numeric vector"
  )
})


# Test purpose: Verifies that force_max_fp_vars applies only to MFPI's
# adjustment-model selection and sets both select and alpha to 1 under the
# p-value criterion.
test_that("MFPI force_max_fp_vars forces p-value adjustment variables", {
  dat <- make_mfpi_factor_data()

  # Replace z after generating y so that z is a positive, unassociated
  # adjustment variable. Its retention therefore depends on force_max_fp_vars,
  # while no shift is needed for the FP transformation.
  set.seed(1503)
  dat$z <- runif(nrow(dat), min = 0.5, max = 3)

  fit <- mfpi(
    y ~ trt + x + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex3",
    criterion = "pvalue",
    force_max_fp_vars = "z",
    df = 4,
    select = 0.05,
    alpha = 0.05,
    center = FALSE,
    cycles = 5,
    p_interact = 1,
    verbose = FALSE
  )

  adjustment_terms <- fit$adjustment_model$fp_terms

  expect_equal(as.numeric(adjustment_terms["z", "select"]), 1)
  expect_equal(as.numeric(adjustment_terms["z", "alpha"]), 1)
  expect_true(adjustment_terms["z", "selected"])
  expect_equal(
    sum(!is.na(fit$adjustment_model$fp_powers[["z"]])),
    2
  )
})
