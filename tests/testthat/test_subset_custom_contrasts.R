# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: The retained-frame helper preserves an explicit factor contrast
# when a genuine subset retains every level, including contrasts created by C().
test_that("24.15 formula subset helper preserves custom contrasts when all levels remain", {
  stage <- factor(rep(c("A", "B", "C"), each = 12L))
  contrasts(stage) <- stats::contr.sum(3L)
  # Only factor levels/contrasts are under test; keep the response deterministic.
  dat <- data.frame(y = seq_along(stage), stage = stage)
  mf <- stats::model.frame(y ~ stage, data = dat)
  rows <- c(1:8, 13:20, 25:32)

  retained <- subset_formula_model_frame(mf, rows)

  expect_identical(levels(retained$stage), levels(mf$stage))
  expect_identical(
    attr(retained$stage, "contrasts", exact = TRUE),
    attr(mf$stage, "contrasts", exact = TRUE)
  )

  mf_c <- stats::model.frame(y ~ C(stage, stats::contr.sum), data = dat)
  retained_c <- subset_formula_model_frame(mf_c, rows)
  factor_name <- setdiff(names(mf_c), "y")
  expect_identical(
    attr(retained_c[[factor_name]], "contrasts", exact = TRUE),
    attr(mf_c[[factor_name]], "contrasts", exact = TRUE)
  )
})


# Test purpose: Both formula methods retain factor-specific custom contrasts in
# their fitted prediction metadata when all factor levels survive subsetting.
test_that("24.16 formula methods retain custom contrasts when all levels remain", {
  set.seed(2416)
  n <- 120L
  stage <- factor(rep(c("A", "B", "C"), length.out = n))
  contrasts(stage) <- stats::contr.sum(3L)
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.3 * x + 0.5 * (stage == "B") + rnorm(n, sd = 0.25),
    trt = trt,
    x = x,
    stage = stage,
    z = rnorm(n),
    keep_row = rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)
  )

  expected_mf <- stats::model.frame(
    y ~ x + stage,
    data = dat[dat$keep_row, , drop = FALSE]
  )
  expected_mm <- stats::model.matrix(y ~ x + stage, data = expected_mf)

  fit_mfp2 <- mfp2(
    y ~ x + stage,
    data = dat,
    subset = keep_row,
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_identical(
    fit_mfp2$formula_contrasts$stage,
    attr(expected_mm, "contrasts")$stage
  )
  expect_equal(
    unname(fit_mfp2$x_original[, c("stage1", "stage2"), drop = FALSE]),
    unname(expected_mm[, c("stage1", "stage2"), drop = FALSE]),
    tolerance = 1e-12
  )

  expected_mfpi_mf <- stats::model.frame(
    y ~ trt + x + stage + z,
    data = dat[dat$keep_row, , drop = FALSE]
  )
  expected_mfpi_mm <- stats::model.matrix(
    y ~ trt + x + stage + z,
    data = expected_mfpi_mf
  )

  fit_mfpi <- mfpi(
    y ~ trt + x + stage + z,
    data = dat,
    subset = keep_row,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = c("stage", "z"),
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
    fit_mfpi$formula_contrasts$stage,
    attr(expected_mfpi_mm, "contrasts")$stage
  )
  expect_equal(
    unname(fit_mfpi$x_train_internal[, c("stage1", "stage2"), drop = FALSE]),
    unname(expected_mfpi_mm[, c("stage1", "stage2"), drop = FALSE]),
    tolerance = 1e-12
  )
})


# Test purpose: A custom contrast matrix cannot be reduced uniquely when a
# subset removes one of its factor levels, so both formula methods stop clearly.
test_that("24.17 formula methods reject ambiguous custom contrasts after level removal", {
  set.seed(2417)
  n <- 120L
  stage <- factor(rep(c("A", "B", "C"), length.out = n))
  contrasts(stage) <- stats::contr.sum(3L)
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.3 * x + rnorm(n, sd = 0.25),
    trt = trt,
    x = x,
    stage = stage,
    z = rnorm(n)
  )
  keep_rows <- dat$stage != "C"
  message <- "Custom contrasts for factor `stage` cannot be reconstructed safely"

  expect_error(
    mfp2(
      y ~ x + stage,
      data = dat,
      subset = keep_rows,
      keep = "stage",
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    message,
    fixed = TRUE
  )

  expect_error(
    mfpi(
      y ~ trt + x + stage + z,
      data = dat,
      subset = keep_rows,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      keep = c("stage", "z"),
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    message,
    fixed = TRUE
  )
})


# Test purpose: The custom-contrast inspection is confined to the non-NULL
# subset helper path. An ordinary fit reuses its full model frame unchanged.
test_that("24.18 subset NULL bypasses custom contrast processing", {
  set.seed(2418)
  n <- 90L
  stage <- factor(rep(c("A", "B", "C"), length.out = n))
  contrasts(stage) <- stats::contr.sum(3L)
  dat <- data.frame(y = rnorm(n), x = runif(n, 1, 5), stage = stage)

  testthat::local_mocked_bindings(
    subset_formula_factor_column = function(...) {
      stop("factor subset helper must not be called for subset = NULL", call. = FALSE)
    },
    .package = "mfp2"
  )

  fit <- mfp2(
    y ~ x + stage,
    data = dat,
    keep = "stage",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_identical(
    fit$formula_contrasts$stage,
    attr(stats::model.matrix(y ~ x + stage, data = dat), "contrasts")$stage
  )
})


# Test purpose: Inline factor wrappers, binary factors, and interactions that
# include a categorical variable must all be rebuilt from the retained rows.
test_that("24.19 formula subset rebuilds inline, binary, and interaction coding", {
  set.seed(2419)
  dat <- expand.grid(
    stage_code = 1:3,
    binary_code = 0:1,
    group = c("A", "B", "C"),
    replicate = 1:10,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dat$group <- factor(dat$group, levels = c("A", "B", "C"))
  dat$x <- runif(nrow(dat), 1, 6)
  dat$y <- 0.4 * dat$x +
    0.5 * (dat$stage_code == 2L) -
    0.3 * (dat$stage_code == 3L) +
    0.6 * dat$binary_code +
    0.2 * dat$x * (dat$group == "B") -
    0.15 * dat$x * (dat$group == "C") +
    rnorm(nrow(dat), sd = 0.25)
  dat$keep_row <- seq_len(nrow(dat)) %% 7L != 0L

  fit <- mfp2(
    y ~ x + factor(stage_code) + as.factor(binary_code) +
      group + x:group,
    data = dat,
    subset = keep_row,
    keep = c("x", "stage_code", "binary_code", "group", "x:group"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 5,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  retained <- droplevels(dat[dat$keep_row, , drop = FALSE])
  expected <- stats::model.matrix(
    ~ x + factor(stage_code) + as.factor(binary_code) +
      group + x:group,
    data = retained
  )[, -1L, drop = FALSE]

  expect_identical(fit$formula_model_matrix_columns, colnames(expected))
  expect_equal(
    unname(fit$x_original[, colnames(expected), drop = FALSE]),
    unname(expected),
    tolerance = 1e-12
  )
  expect_true(all(c("stage_code", "binary_code", "group", "x:group") %in%
                    names(fit$term_to_columns)))
})
