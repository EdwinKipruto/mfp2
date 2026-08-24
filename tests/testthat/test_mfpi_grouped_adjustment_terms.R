# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# -----------------------------------------------------------------------------
# 16.1 Grouped categorical adjustment terms
# -----------------------------------------------------------------------------
# These tests verify that factor contrasts and manually specified dummy blocks
# remain one conceptual adjustment term throughout Stage 1 selection and Stage 2
# interaction-model construction.

# Test purpose: Verifies that unordered-factor contrast columns are stored
# and selected as one conceptual MFPI adjustment term.
test_that("MFPI formula groups unordered-factor adjustment columns", {
  dat <- make_mfpi_factor_data()

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true("stage" %in% names(fit$term_to_columns))
  expect_setequal(fit$term_to_columns[["stage"]], c("stageII", "stageIII"))

  expect_true("stage" %in% rownames(fit$adjustment_model$fp_terms))
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_setting"]),
    1
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_initial"]),
    2
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_final"]),
    2
  )
  expect_false(any(c("stageII", "stageIII") %in%
                     rownames(fit$adjustment_model$fp_terms)))
})


# Test purpose: Verifies that ordered-factor polynomial contrasts are accepted
# and grouped as one fixed linear MFPI adjustment term.
test_that("MFPI formula groups ordered-factor polynomial contrasts", {
  dat <- make_mfpi_factor_data(ordered_stage = TRUE)

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_true("stage" %in% names(fit$term_to_columns))
  expect_setequal(fit$term_to_columns[["stage"]], c("stage.L", "stage.Q"))
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_setting"]),
    1
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_initial"]),
    2
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_final"]),
    2
  )
})


# Test purpose: Checks that MFPI uses the source variable name for a simple
# inline factor wrapper while retaining the original contrast-column names.
test_that("MFPI names inline factor adjustments by their source variable", {
  dat <- make_mfpi_factor_data()
  dat$stage_code <- as.integer(dat$stage)

  fit <- mfpi(
    y ~ trt + x + factor(stage_code) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage_code",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_true("stage_code" %in% names(fit$term_to_columns))
  expect_false("factor(stage_code)" %in% names(fit$term_to_columns))
  expect_identical(
    fit$term_to_columns[["stage_code"]],
    c("factor(stage_code)2", "factor(stage_code)3")
  )
  expect_true("stage_code" %in% rownames(fit$adjustment_model$fp_terms))
  expect_false(
    "factor(stage_code)" %in% rownames(fit$adjustment_model$fp_terms)
  )
})


# Test purpose: Verifies that a binary inline factor is retained as a mapped
# one-column MFPI adjustment term rather than reverting to its dummy-column name.
test_that("MFPI retains binary inline factors as source-name mappings", {
  dat <- make_mfpi_factor_data()
  dat$binary_stage <- rep(c(1, 1, 2, 1, 2, 2), length.out = nrow(dat))

  fit <- mfpi(
    y ~ trt + x + factor(binary_stage) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "binary_stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    fit$term_to_columns[["binary_stage"]],
    "factor(binary_stage)2"
  )
  expect_true(fit$adjustment_model$fp_terms["binary_stage", "selected"])

  nd <- dat[1:18, c("trt", "x", "binary_stage", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(fit, "x", fit_result, nd)
  expect_true("factor(binary_stage)2.1" %in% names(design$model_newdata))
})


# Test purpose: Checks manual matrix-interface grouping of dummy columns in
# the Stage 1 MFPI adjustment model.
test_that("MFPI matrix interface accepts manual grouped adjustment columns", {
  dat <- make_mfpi_factor_data()
  stage_mm <- stats::model.matrix(~ stage, dat)[, -1L, drop = FALSE]

  x <- cbind(
    trt = as.numeric(dat$trt) - 1,
    x = dat$x,
    stage_mm,
    z = dat$z
  )

  fit <- mfpi(
    x,
    dat$y,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = c("stageII", "stageIII")),
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    fit$term_to_columns[["stage"]],
    c("stageII", "stageIII")
  )
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
})


# Test purpose: MFPI applies the scalar df default only to ordinary
# continuous columns; a supplied ordered-factor contrast block remains linear.
test_that("MFPI grouped ordered-factor columns are linear under scalar df default", {
  set.seed(51231)
  n <- 192
  trt <- rep(rep(0:1, each = 4L), length.out = n)
  x_cont <- seq(0.5, 9, length.out = n)
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  z <- stats::rnorm(n)

  x <- cbind(
    trt = trt,
    x = x_cont,
    stage_mm,
    z = z
  )
  y <- 0.5 * x_cont + 0.7 * stage_mm[, 1L] -
    0.4 * stage_mm[, 2L] + 0.2 * z + stats::rnorm(n, sd = 0.5)

  fit <- mfpi(
    x,
    y,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_setting"]),
    1
  )
  expect_equal(
    as.numeric(fit$adjustment_model$fp_terms["stage", "df_initial"]),
    ncol(stage_mm)
  )
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
})


# Test purpose: MFPI must preserve a supplied categorical contrast block under
# scale one instead of estimating a separate scale for each member column.
test_that("MFPI grouped matrix columns default to scale one", {
  set.seed(51232)
  n <- 192
  trt <- rep(rep(0:1, each = 4L), length.out = n)
  x_cont <- seq(0.5, 9, length.out = n)
  stage <- ordered(rep(LETTERS[1:4], length.out = n))
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  stage_mm[, 1L] <- 1000 * stage_mm[, 1L]
  z <- stats::rnorm(n)

  x <- cbind(
    trt = trt,
    x = x_cont,
    stage_mm,
    z = z
  )
  y <- 0.5 * x_cont + 0.001 * stage_mm[, 1L] -
    0.4 * stage_mm[, 2L] + 0.2 * z + stats::rnorm(n, sd = 0.5)

  fit <- mfpi(
    x,
    y,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    term_groups = list(stage = colnames(stage_mm)),
    keep = "stage",
    p_interact = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(
    unname(fit$scale[colnames(stage_mm)]),
    rep(1, ncol(stage_mm))
  )
  expect_equal(
    as.numeric(fit$adjustment_model$transformations["stage", "scale"]),
    1
  )
  expect_true(fit$adjustment_model$fp_terms["stage", "selected"])
})


# Test purpose: Ensures grouped categorical terms and their member columns
# cannot be requested as continuous MFPI interaction variables.
test_that("MFPI grouped terms cannot be used as cont_vars", {
  dat <- make_mfpi_factor_data()
  stage_mm <- stats::model.matrix(~ stage, dat)[, -1L, drop = FALSE]

  x <- cbind(
    trt = as.numeric(dat$trt) - 1,
    x = dat$x,
    stage_mm
  )

  expect_error(
    mfpi(
      x,
      dat$y,
      group_var = "trt",
      cont_vars = "stageII",
      term_groups = list(stage = c("stageII", "stageIII")),
      verbose = FALSE
    ),
    "singleton continuous|grouped categorical"
  )
})


# Test purpose: Checks that a multilevel group_var is represented by one
# conceptual dummy block when include_group_var = TRUE.
test_that("include_group_var stores group dummies as one conceptual term", {
  dat <- make_mfpi_factor_data()
  # Use a three-level pattern that is not collinear with the three-level stage
  # factor; otherwise the null and group-dummy models have the same fitted rank.
  dat$trt <- factor(
    rep(c("A", "A", "B", "C", "B", "C"), length.out = nrow(dat)),
    levels = c("A", "B", "C")
  )

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    include_group_var = TRUE,
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_true("trt" %in% names(fit$adjustment_model$term_to_columns))
  expect_length(fit$adjustment_model$term_to_columns[["trt"]], 2L)

  selected_for_stage2 <- fit$univariable_interactions$selected_vars
  if (!is.null(selected_for_stage2)) {
    expect_false("trt" %in% selected_for_stage2)
  }
})
