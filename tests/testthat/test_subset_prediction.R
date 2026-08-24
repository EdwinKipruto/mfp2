# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Predictions from a subset-fitted formula model must reproduce
# training predictions, accept retained levels, reject removed/unseen levels,
# and retain the same behavior after serialization.
test_that("24.25 subset-fitted mfp2 prediction uses retained factor levels", {
  set.seed(2425)
  n <- 180L
  dat <- data.frame(
    x = runif(n, 1, 5),
    stage = factor(rep(c("A", "B", "C"), length.out = n),
                   levels = c("A", "B", "C"))
  )
  dat$y <- 1 + 0.4 * dat$x + 0.7 * (dat$stage == "B") +
    rnorm(n, sd = 0.25)
  fit_rows <- dat$stage != "C"

  fit <- mfp2(
    y ~ x + stage,
    data = dat,
    subset = fit_rows,
    keep = c("x", "stage"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  retained <- droplevels(dat[fit_rows, , drop = FALSE])

  training <- predict(fit, type = "link", se.fit = TRUE)
  reconstructed <- predict(
    fit,
    newdata = retained,
    type = "link",
    se.fit = TRUE
  )
  expect_equal(unname(reconstructed$fit), unname(training$fit), tolerance = 1e-10)
  expect_equal(unname(reconstructed$se.fit), unname(training$se.fit), tolerance = 1e-10)

  retained_row <- retained[retained$stage == "B", , drop = FALSE][1L, ]
  expect_length(predict(fit, newdata = retained_row), 1L)

  removed_row <- dat[dat$stage == "C", , drop = FALSE][1L, ]
  expect_error(
    predict(fit, newdata = removed_row),
    "new level|new levels|factor",
    ignore.case = TRUE
  )

  unseen_row <- retained_row
  unseen_row$stage <- factor("D", levels = c("A", "B", "D"))
  expect_error(
    predict(fit, newdata = unseen_row),
    "new level|new levels|factor",
    ignore.case = TRUE
  )

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  restored_prediction <- predict(
    restored,
    newdata = retained,
    type = "link",
    se.fit = TRUE
  )
  expect_equal(restored_prediction$fit, reconstructed$fit, tolerance = 1e-12)
  expect_equal(restored_prediction$se.fit, reconstructed$se.fit, tolerance = 1e-12)
})


# Test purpose: MFPI formula metadata, formula-offset reconstruction, and
# ordinary prediction must all use only factor levels retained by subset.
test_that("24.26 MFPI subset prediction uses retained metadata and formula offset", {
  set.seed(2426)
  n <- 180L
  dat <- data.frame(
    trt = factor(rep(c("control", "treated"), length.out = n)),
    x = runif(n, 1, 5),
    stage = factor(rep(c("A", "B", "C"), length.out = n),
                   levels = c("A", "B", "C")),
    z = rnorm(n),
    off = rnorm(n, sd = 0.15)
  )
  dat$y <- 0.7 + 0.45 * dat$x + 0.6 * (dat$trt == "treated") * dat$x +
    0.4 * (dat$stage == "B") + 0.2 * dat$z + dat$off +
    rnorm(n, sd = 0.3)
  fit_rows <- dat$stage != "C"

  fit <- mfpi(
    y ~ trt + x + stage + z + offset(off),
    data = dat,
    subset = fit_rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = c("stage", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    cycles = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(fit$formula_xlevels$stage, c("A", "B"))
  expect_false(is.null(fit$formula_offset_terms))
  expect_true("stageB" %in% fit$formula_design_columns)
  expect_false("stageC" %in% fit$formula_design_columns)

  retained <- droplevels(
    dat[fit_rows, c("trt", "x", "stage", "z", "off"), drop = FALSE]
  )
  reconstructed <- predict(
    fit,
    newdata = retained,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )
  training <- predict(
    fit,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )
  expect_equal(
    reconstructed$predictions,
    training$predictions,
    tolerance = 1e-10
  )

  retained_row <- retained[retained$stage == "B", , drop = FALSE][1L, ]
  retained_prediction <- predict(
    fit,
    newdata = retained_row,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expect_true(all(is.finite(retained_prediction$predictions$fit)))

  removed_row <- dat[
    dat$stage == "C",
    c("trt", "x", "stage", "z", "off"),
    drop = FALSE
  ][1L, ]
  expect_error(
    predict(
      fit,
      newdata = removed_row,
      terms = "x",
      model = "all",
      type = "link",
      se.fit = FALSE
    ),
    "new level|new levels|factor",
    ignore.case = TRUE
  )

  unseen_row <- retained_row
  unseen_row$stage <- factor("D", levels = c("A", "B", "D"))
  expect_error(
    predict(
      fit,
      newdata = unseen_row,
      terms = "x",
      model = "all",
      type = "link",
      se.fit = FALSE
    ),
    "new level|new levels|factor",
    ignore.case = TRUE
  )
})
