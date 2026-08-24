# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 24. subset handling across formula and matrix interfaces
# =============================================================================

# Test purpose: The formula interface must rebuild ordered-factor contrasts on
# the fitted rows while retaining full-data scaling for continuous predictors.
test_that("24.1 mfp2 formula rebuilds categorical coding after subset", {
  set.seed(2401)
  n_each <- 35L
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = n_each),
    levels = c("A", "B", "C", "D")
  )
  x <- seq_len(length(stage))
  x[stage == "D"] <- x[stage == "D"] * 10000
  dat <- data.frame(
    y = 1 + 0.03 * seq_len(length(stage)) + rnorm(length(stage), sd = 0.2),
    x = x,
    stage = stage
  )
  fit_rows <- dat$stage != "D"

  fit <- mfp2(
    y ~ x + stage,
    data = dat,
    subset = fit_rows,
    keep = c("x", "stage"),
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    verbose = FALSE
  )

  dat_fit <- droplevels(dat[fit_rows, , drop = FALSE])
  expected_x <- stats::model.matrix(~ x + stage, data = dat_fit)[, -1L, drop = FALSE]

  expect_identical(fit$formula_xlevels$stage, levels(dat_fit$stage))
  expect_identical(fit$formula_model_matrix_columns, colnames(expected_x))
  expect_identical(fit$term_to_columns$stage, c("stage.L", "stage.Q"))
  expect_equal(
    unname(fit$x_original[, c("stage.L", "stage.Q"), drop = FALSE]),
    unname(expected_x[, c("stage.L", "stage.Q"), drop = FALSE]),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(fit$transformations["x", "scale"]),
    find_scale_factor(dat$x)
  )
  expect_false(isTRUE(all.equal(
    find_scale_factor(dat$x),
    find_scale_factor(dat_fit$x)
  )))
  expect_error(
    predict(fit, newdata = dat[dat$stage == "D", , drop = FALSE][1L, ]),
    "new level|new levels|factor"
  )
})


# Test purpose: A formula factor with only one retained level has no estimable
# effect and must fail with a categorical-specific diagnostic.
test_that("24.2 mfp2 formula rejects a one-level factor after subset", {
  set.seed(2402)
  dat <- data.frame(
    y = rnorm(60),
    x = runif(60, 1, 5),
    group = factor(rep(c("A", "B", "C"), each = 20))
  )

  expect_error(
    mfp2(
      y ~ x + group,
      data = dat,
      subset = dat$group == "A",
      df = 1,
      verbose = FALSE
    ),
    "fewer than two observed levels"
  )
})


# Test purpose: The matrix interface cannot regenerate reduced polynomial
# contrasts and therefore rejects a grouped block that loses rank after subset.
test_that("24.3 mfp2 matrix rejects grouped rank loss after subset", {
  set.seed(2403)
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = 25),
    levels = c("A", "B", "C", "D")
  )
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  x <- cbind(z = runif(length(stage), 1, 4), stage_mm)
  y <- rnorm(length(stage))

  expect_error(
    mfp2(
      x,
      y,
      subset = stage != "D",
      term_groups = list(stage = colnames(stage_mm)),
      keep = "stage",
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    "lose estimable dimension"
  )
})


# Test purpose: MFPI formula fitting uses the retained factor levels and
# regenerates the corresponding polynomial-contrast block.
test_that("24.4 MFPI formula rebuilds categorical coding after subset", {
  set.seed(2404)
  n_each <- 35L
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = n_each),
    levels = c("A", "B", "C", "D")
  )
  trt <- factor(rep(c("control", "treated"), length.out = length(stage)))
  x <- runif(length(stage), 1, 7)
  z <- rnorm(length(stage))
  y <- 1 + 0.4 * x + 0.7 * x * (trt == "treated") + 0.2 * z + rnorm(length(stage), sd = 0.3)
  dat <- data.frame(y = y, trt = trt, x = x, stage = stage, z = z)
  fit_rows <- dat$stage != "D"

  fit <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    subset = fit_rows,
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

  expect_identical(fit$formula_xlevels$stage, c("A", "B", "C"))
  expect_identical(fit$term_to_columns$stage, c("stage.L", "stage.Q"))
  expect_false("stage.C" %in% fit$formula_design_columns)
})


# Test purpose: MFPI formula fitting gives a direct diagnostic when a factor is
# reduced to one observed level by subset.
test_that("24.5 MFPI formula rejects a one-level factor after subset", {
  set.seed(2405)
  dat <- make_mfpi_factor_data(n = 120L)

  expect_error(
    mfpi(
      y ~ trt + x + stage + z,
      data = dat,
      subset = dat$stage == "I",
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    "fewer than two observed levels"
  )
})


# Test purpose: The MFPI matrix interface rejects supplied grouped contrast
# columns when subsetting removes an estimable dimension.
test_that("24.6 MFPI matrix rejects grouped rank loss after subset", {
  set.seed(2406)
  n_each <- 30L
  stage <- ordered(
    rep(c("A", "B", "C", "D"), each = n_each),
    levels = c("A", "B", "C", "D")
  )
  stage_mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  trt <- rep(c(0, 1), length.out = length(stage))
  x_cont <- runif(length(stage), 1, 6)
  x <- cbind(trt = trt, x = x_cont, stage_mm, z = rnorm(length(stage)))
  y <- 1 + 0.4 * x_cont + 0.6 * x_cont * trt + rnorm(length(stage), sd = 0.3)

  expect_error(
    mfpi(
      x,
      y,
      subset = stage != "D",
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      term_groups = list(stage = colnames(stage_mm)),
      keep = "stage",
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    "lose estimable dimension"
  )
})


# Test purpose: A raw matrix column that becomes constant after subsetting is
# rejected even when it is not declared through term_groups.
test_that("24.7 MFPI matrix rejects a predictor with no post-subset variation", {
  set.seed(2407)
  n <- 120L
  trt <- rep(c(0, 1), length.out = n)
  x_cont <- runif(n, 1, 5)
  indicator <- rep(c(0, 0, 1, 1), length.out = n)
  x <- cbind(
    trt = trt,
    x = x_cont,
    indicator = indicator,
    z = rnorm(n)
  )
  y <- 1 + 0.3 * x_cont + 0.5 * x_cont * trt + rnorm(n, sd = 0.3)

  expect_error(
    mfpi(
      x,
      y,
      subset = indicator == 0,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    "no variation after applying `subset`"
  )
})


# Test purpose: Removing the original treatment-contrast reference level causes
# the formula interface to rebuild the factor with a retained reference level.
test_that("24.8 mfp2 formula rebuilds treatment contrasts when reference is removed", {
  set.seed(2408)
  group <- factor(
    rep(c("A", "B", "C"), each = 30),
    levels = c("A", "B", "C")
  )
  dat <- data.frame(
    y = rnorm(length(group)),
    x = runif(length(group), 1, 5),
    group = group
  )
  fit_rows <- dat$group != "A"

  fit <- mfp2(
    y ~ x + group,
    data = dat,
    subset = fit_rows,
    keep = "group",
    df = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  dat_fit <- droplevels(dat[fit_rows, , drop = FALSE])
  expected_x <- stats::model.matrix(~ x + group, dat_fit)[, -1L, drop = FALSE]

  expect_identical(fit$formula_xlevels$group, c("B", "C"))
  expect_identical(fit$term_to_columns$group, "groupC")
  expect_equal(
    unname(fit$x_original[, "groupC", drop = FALSE]),
    unname(expected_x[, "groupC", drop = FALSE]),
    tolerance = 1e-12
  )
})


# Test purpose: Formula methods must bypass all subset-specific processing
# when subset is NULL. This guards the ordinary formula path against both the
# model.frame() NSE regression and unnecessary row/factor reconstruction.
test_that("24.9 formula methods bypass subset helpers when subset is NULL", {
  set.seed(2409)
  n <- 90L
  dat <- data.frame(
    y = rnorm(n),
    trt = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 5),
    stage = factor(rep(c("A", "B", "C"), length.out = n)),
    z = rnorm(n)
  )

  testthat::local_mocked_bindings(
    formula_subset_rows = function(...) {
      stop("row resolver must not be called for subset = NULL", call. = FALSE)
    },
    subset_formula_model_frame = function(...) {
      stop("frame helper must not be called for subset = NULL", call. = FALSE)
    },
    .package = "mfp2"
  )

  fit_mfp2 <- mfp2(
    y ~ x + stage,
    data = dat,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  expect_s3_class(fit_mfp2, "mfp2")

  fit_mfpi <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_s3_class(fit_mfpi, "mfpi")
})
