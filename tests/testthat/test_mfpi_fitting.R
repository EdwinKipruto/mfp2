# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Checks that the default MFPI interface fits an interaction-analysis
# object with expected metadata.
test_that("mfpi.default() runs and returns an mfpi object", {
  data("prostate", package = "mfp2")

  x_p <- data.frame(
    svi = prostate$svi,
    age = prostate$age,
    cavol = prostate$cavol,
    pgg45 = prostate$pgg45,
    weight = prostate$weight,
    bph = prostate$bph,
    cp = prostate$cp
  )

  fit <- mfpi(
    x_p, y_prostate,
    group_var = "svi",
    cont_vars = c("cavol", "age"),
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(!is.null(fit$all_model_metrics))
  expect_true(!is.null(fit$adjustment_model))
  expect_equal(fit$group_var, "svi")
  expect_equal(fit$cont_vars, c("cavol", "age"))
})


# Test purpose: Checks that the formula MFPI interface parses fp() terms and
# returns an mfpi object.
test_that("mfpi.formula() runs and returns an mfpi object", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
      fp(bph) + fp(cp),
    data = prostate,
    cont_vars = c("cavol", "age"),
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
})


# Test purpose: Checks that requested interaction functional forms are stored for
# continuous variables.
test_that("mfpi() cont_var_forms specifies functional form correctly", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    cont_var_forms = c(cavol = "fp2", age = "linear"),
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "fp2"))
  expect_equal(fit$cont_var_forms["age"], c(age = "linear"))
})


# Test purpose: Checks that MFPI runs with information-criterion based interaction
# assessment.
test_that("mfpi() with criterion = 'aic' works", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol"),
    group_var = "svi",
    criterion = "aic",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
})


# Test purpose: Checks that all supported MFPI flexibility settings fit without
# error.
test_that("mfpi() flexibility levels run without error", {
  data("prostate", package = "mfp2")

  for (fl in c("flex1", "flex2", "flex3", "flex4")) {
    fit <- mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      cont_var_forms = c(cavol = "fp1"),
      group_var = "svi",
      flex = fl,
      verbose = FALSE
    )
    expect_true(
      inherits(fit, "mfpi"),
      info = paste("flex =", fl)
    )
  }
})


# Test purpose: default behavior
test_that("mfpi() defaults to flex3", {
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    verbose = FALSE
  )
  expect_equal(fit$flex, "flex3")
})


# Test purpose: Checks that omitted cont_var_forms are filled with "linear"
# for every variable listed in cont_vars.
test_that("mfpi() defaults missing cont_var_forms to linear", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    group_var = "svi",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "fp1"))
  expect_equal(fit$cont_var_forms["age"], c(age = "fp1"))
})


# Test purpose: Ensures cont_var_forms may specify only some cont_vars;
# omitted cont_vars are filled with "linear" and ordering follows cont_vars.
test_that("mfpi() fills missing cont_var_forms entries with linear", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    cont_var_forms = c(cavol = "fp2"),
    group_var = "svi",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(names(fit$cont_var_forms), c("cavol", "age"))
  expect_equal(fit$cont_var_forms["cavol"], c(cavol = "fp2"))
  expect_equal(fit$cont_var_forms["age"], c(age = "fp1"))
})


# Test purpose: Verifies that cont_var_forms only accepts "linear", "fp1",
# and "fp2".
test_that("mfpi() rejects invalid cont_var_forms values", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = c("cavol", "age"),
      cont_var_forms = c(cavol = "spline"),
      group_var = "svi",
      verbose = FALSE
    ),
    "Invalid value"
  )
})


# Test purpose: Ensures cont_var_forms entries must be named so each requested
# form is explicitly tied to a variable in cont_vars.
test_that("mfpi() rejects unnamed cont_var_forms", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = c("cavol", "age"),
      cont_var_forms = c("fp2", "linear"),
      group_var = "svi",
      verbose = FALSE
    ),
    "Every entry of `cont_var_forms` must be named"
  )
})


# Test purpose: Checks that cont_var_forms cannot name variables that are not
# being tested as continuous interaction variables.
test_that("mfpi() rejects cont_var_forms names not in cont_vars", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      cont_var_forms = c(age = "fp1"),
      group_var = "svi",
      verbose = FALSE
    ),
    "not in `cont_vars`"
  )
})


# Test purpose: Ensures the grouping variable cannot also be listed as a
# continuous interaction variable.
test_that("mfpi() rejects group_var included in cont_vars", {
  set.seed(204)
  n <- 120

  x <- data.frame(
    group = rep(1:4, length.out = n),
    x1 = runif(n, 1, 10),
    x2 = runif(n, 1, 10)
  )
  y <- rnorm(n)

  expect_error(
    mfpi(
      x,
      y,
      group_var = "group",
      cont_vars = c("group", "x1"),
      verbose = FALSE
    ),
    "must not also appear in `cont_vars`"
  )
})


# Test purpose: Checks that binary variables are not allowed in cont_vars,
# because MFPI requires continuous variables for FP interaction forms.
test_that("mfpi.formula() rejects a binary variable in cont_vars", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "svi",
      group_var = "cavol",
      verbose = FALSE
    ),
    "binary"
  )
})


# Test purpose: Checks that binary variables are not allowed in cont_vars,
# because MFPI requires continuous variables for FP interaction forms.
test_that("mfpi.formula() rejects a simulated binary variable in cont_vars", {
  set.seed(201)
  n <- 100

  dat <- data.frame(
    y = rnorm(n),
    group = rep(0:1, length.out = n),
    binary_x = rep(0:1, length.out = n),
    x = runif(n, 1, 10)
  )

  expect_error(
    mfpi(
      y ~ group + binary_x + fp(x),
      data = dat,
      cont_vars = "binary_x",
      group_var = "group",
      verbose = FALSE
    ),
    "binary"
  )
})


# Test purpose: Ensures MFPI validates multiplicity-adjustment methods against
# stats::p.adjust.methods.
test_that("mfpi() rejects invalid p_adjust_method", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      group_var = "svi",
      p_adjust_method = "not_a_method",
      verbose = FALSE
    ),
    "Invalid `p_adjust_method`"
  )
})


# Test purpose: Checks that a valid multiplicity-adjustment method is accepted
# and stored on the returned mfpi object.
test_that("mfpi() stores valid p_adjust_method", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = c("cavol", "age"),
    group_var = "svi",
    p_adjust_method = "holm",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$p_adjust_method, "holm")
})


# Test purpose: Verifies criterion-specific min_improvement defaults:
# pvalue uses p_interact, while AIC and BIC default to 2.
test_that("mfpi() sets criterion-specific default min_improvement", {
  data("prostate", package = "mfp2")

  fit_p <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "pvalue",
    p_interact = 0.10,
    verbose = FALSE
  )

  fit_aic <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "aic",
    verbose = FALSE
  )

  fit_bic <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "bic",
    verbose = FALSE
  )

  expect_equal(fit_p$min_improvement, 0.10)
  expect_equal(fit_aic$min_improvement, 2)
  expect_equal(fit_bic$min_improvement, 2)
})


# Test purpose: Checks that an explicit min_improvement threshold is respected
# for information-criterion based MFPI selection.
test_that("mfpi() stores explicit min_improvement for AIC", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    criterion = "aic",
    min_improvement = 3.5,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$min_improvement, 3.5)
})


# Test purpose: Ensures min_improvement must be NULL or a single positive
# finite numeric value.
test_that("mfpi() rejects invalid min_improvement", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      cont_vars = "cavol",
      group_var = "svi",
      criterion = "aic",
      min_improvement = 0,
      verbose = FALSE
    ),
    "`min_improvement`"
  )
})


# Test purpose: Checks that include_group_var = TRUE fits successfully and is
# recorded on the returned mfpi object.
test_that("mfpi() accepts include_group_var = TRUE", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    include_group_var = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(isTRUE(fit$include_group_var))
})


# Test purpose: Ensures group-specific centering mode is accepted and stored.
test_that("mfpi() accepts group-specific centering", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    center_type = "group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$center_type, "group")
})


# Test purpose: Ensures mfpi.default() allows group_var to be categorical but
# rejects other categorical predictors in x.
test_that("mfpi.default() rejects non-group categorical predictors", {
  set.seed(202)
  n <- 100

  x <- data.frame(
    group = factor(rep(c("A", "B"), length.out = n)),
    x = runif(n, 1, 10),
    bad_factor = factor(rep(c("low", "high"), length.out = n))
  )
  y <- rnorm(n)

  expect_error(
    mfpi(
      x,
      y,
      group_var = "group",
      cont_vars = "x",
      verbose = FALSE
    ),
    "Only `group_var` may be categorical"
  )
})


# Test purpose: Checks that categorical group labels are retained as metadata
# after internal recoding of group_var.
test_that("mfpi.default() stores original group levels", {
  set.seed(203)
  n <- 120

  x <- data.frame(
    group = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 10),
    z = runif(n, 1, 10)
  )
  y <- 0.2 * x$x + 0.5 * (x$group == "treated") + rnorm(n)

  fit <- mfpi(
    x,
    y,
    group_var = "group",
    cont_vars = "x",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(!is.null(fit$group_levels_original))
  expect_true(all(c("control", "treated") %in% fit$group_levels_original))
})
