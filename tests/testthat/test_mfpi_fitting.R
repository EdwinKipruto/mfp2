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
    interaction_vars = c("cavol", "age"),
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(!is.null(fit$all_model_metrics))
  expect_true(!is.null(fit$adjustment_model))
  expect_equal(fit$group_var, "svi")
  expect_equal(fit$interaction_vars, c("cavol", "age"))
})


# Test purpose: Checks that the formula MFPI interface parses fp() terms and
# returns an mfpi object.
test_that("mfpi.formula() runs and returns an mfpi object", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
      fp(bph) + fp(cp),
    data = prostate,
    interaction_vars = c("cavol", "age"),
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
})


# Formula strata follow the coxph convention and must be written as a formula
# special. The matrix interface continues to accept its strata argument.
test_that("mfpi.formula() rejects a direct strata argument", {
  dat <- data.frame(
    time = seq_len(8L),
    status = rep(c(0L, 1L), 4L),
    age = seq(40, 68, length.out = 8L),
    group = factor(rep(c("A", "B"), 4L)),
    centre = factor(rep(c("C1", "C2"), each = 4L))
  )

  expect_error(
    mfpi(
      survival::Surv(time, status) ~ age + group,
      data = dat,
      family = "cox",
      group_var = "group",
      interaction_vars = "age",
      strata = centre,
      verbose = FALSE
    ),
    paste0(
      "`strata` is not supported as an argument to `mfpi.formula()`. ",
      "Include `strata(...)` in the formula instead."
    ),
    fixed = TRUE
  )
})


# Formula observation arguments follow the same data-mask lookup rule as the
# underlying standard model functions. The bare names below exist only as
# columns of `dat`, so this also guards against accidentally forcing them in the
# test/function environment before formula evaluation.
test_that("mfpi.formula() resolves weights and offset in data", {
  data("prostate", package = "mfp2")
  dat <- prostate
  dat$case_weight <- seq(0.8, 1.2, length.out = nrow(dat))
  dat$external_offset <- seq(-0.1, 0.1, length.out = nrow(dat))

  fit <- mfpi(
    lpsa ~ svi + age + cavol,
    data = dat,
    group_var = "svi",
    interaction_vars = "cavol",
    interaction_forms = c(cavol = "linear"),
    weights = case_weight,
    offset = external_offset,
    flex = "flex1",
    p_interact = 1,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  expect_equal(
    unname(fit$adjustment_model$prior.weights),
    dat$case_weight,
    tolerance = 0
  )
  expect_equal(
    unname(fit$adjustment_model$offset),
    dat$external_offset,
    tolerance = 0
  )
})


# Test purpose: Checks that requested interaction functional forms are stored for
# continuous variables.
test_that("mfpi() interaction_forms specifies functional form correctly", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    interaction_vars = c("cavol", "age"),
    interaction_forms = c(cavol = "fp2", age = "linear"),
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$interaction_forms["cavol"], c(cavol = "fp2"))
  expect_equal(fit$interaction_forms["age"], c(age = "linear"))
})


# Test purpose: Checks that MFPI runs with information-criterion based interaction
# assessment.
test_that("mfpi() with criterion = 'aic' works", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    interaction_vars = c("cavol"),
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
      interaction_vars = "cavol",
      interaction_forms = c(cavol = "fp1"),
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
    interaction_vars = "cavol",
    group_var = "svi",
    verbose = FALSE
  )
  expect_equal(fit$flex, "flex3")
})


# Test purpose: Checks that omitted interaction_forms are filled with "fp1"
# for every variable listed in interaction_vars.
test_that("mfpi() defaults missing interaction_forms to fp1", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    interaction_vars = c("cavol", "age"),
    group_var = "svi",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$interaction_forms["cavol"], c(cavol = "fp1"))
  expect_equal(fit$interaction_forms["age"], c(age = "fp1"))
})


# Test purpose: Ensures interaction_forms may specify only some interaction_vars;
# omitted interaction_vars are filled with "fp1" and ordering follows interaction_vars.
test_that("mfpi() fills missing interaction_forms entries with fp1", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    interaction_vars = c("cavol", "age"),
    interaction_forms = c(cavol = "fp2"),
    group_var = "svi",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(names(fit$interaction_forms), c("cavol", "age"))
  expect_equal(fit$interaction_forms["cavol"], c(cavol = "fp2"))
  expect_equal(fit$interaction_forms["age"], c(age = "fp1"))
})


# Test purpose: Verifies that interaction_forms only accepts "linear", "fp1",
# and "fp2".
test_that("mfpi() rejects invalid interaction_forms values", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      interaction_vars = c("cavol", "age"),
      interaction_forms = c(cavol = "spline"),
      group_var = "svi",
      verbose = FALSE
    ),
    "Invalid value"
  )
})


# Test purpose: Ensures interaction_forms entries must be named so each requested
# form is explicitly tied to a variable in interaction_vars.
test_that("mfpi() rejects unnamed interaction_forms", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      interaction_vars = c("cavol", "age"),
      interaction_forms = c("fp2", "linear"),
      group_var = "svi",
      verbose = FALSE
    ),
    "Every entry of `interaction_forms` must be named"
  )
})


# Test purpose: Checks that interaction_forms cannot name variables that are not
# being tested as continuous interaction variables.
test_that("mfpi() rejects interaction_forms names not in interaction_vars", {
  data("prostate", package = "mfp2")

  expect_error(
    mfpi(
      lpsa ~ fp(age) + svi + fp(cavol),
      data = prostate,
      interaction_vars = "cavol",
      interaction_forms = c(age = "fp1"),
      group_var = "svi",
      verbose = FALSE
    ),
    "not in `interaction_vars`"
  )
})


# Test purpose: Ensures the grouping variable cannot also be listed as a
# interaction variable.
test_that("mfpi() rejects group_var included in interaction_vars", {
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
      interaction_vars = c("group", "x1"),
      verbose = FALSE
    ),
    "must not also appear in `interaction_vars`"
  )
})


# Test purpose: Binary interaction variables use the uncentred linear path.
test_that("mfpi.formula() fits a binary linear interaction", {
  set.seed(201)
  n <- 160

  dat <- data.frame(
    group = rep(0:1, length.out = n),
    binary_x = rep(c(0, 0, 1, 1), length.out = n),
    x = runif(n, 1, 10)
  )
  dat$y <- 0.5 + 0.7 * dat$group - 0.4 * dat$binary_x +
    1.1 * dat$group * dat$binary_x + 0.2 * dat$x + rnorm(n, sd = 0.3)

  fit <- mfpi(
    y ~ group + binary_x + x,
    data = dat,
    interaction_vars = "binary_x",
    interaction_forms = c(binary_x = "linear"),
    group_var = "group",
    select = 1,
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_identical(fit$interaction_specs$binary_x$kind, "binary")
  expect_true(fit$interaction_specs$binary_x$discrete)
  expect_equal(
    fit$var_winners$binary_x$fit$test_results$df$interaction_only,
    1
  )
})

test_that("mfpi() rejects FP forms for binary interaction variables", {
  dat <- data.frame(
    y = rnorm(80),
    group = rep(0:1, each = 40),
    binary_x = rep(c(0, 1), 40),
    x = runif(80, 1, 5)
  )

  expect_error(
    mfpi(
      y ~ group + binary_x + x,
      data = dat,
      group_var = "group",
      interaction_vars = "binary_x",
      interaction_forms = c(binary_x = "fp1"),
      verbose = FALSE
    ),
    "Discrete interaction variable.*must use"
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
      interaction_vars = "cavol",
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
    interaction_vars = c("cavol", "age"),
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
    interaction_vars = "cavol",
    group_var = "svi",
    criterion = "pvalue",
    p_interact = 0.10,
    verbose = FALSE
  )

  fit_aic <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    interaction_vars = "cavol",
    group_var = "svi",
    criterion = "aic",
    verbose = FALSE
  )

  fit_bic <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    interaction_vars = "cavol",
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
    interaction_vars = "cavol",
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
      interaction_vars = "cavol",
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
    interaction_vars = "cavol",
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
    interaction_vars = "cavol",
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
      interaction_vars = "x",
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
    interaction_vars = "x",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true(!is.null(fit$group_levels_original))
  expect_true(all(c("control", "treated") %in% fit$group_levels_original))
})
