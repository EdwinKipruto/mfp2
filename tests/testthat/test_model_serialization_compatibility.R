# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Saving and restoring an mfp2 object must preserve coefficients,
# metadata, training predictions, newdata predictions, and prediction standard
# errors. This protects stored transformations and formula reconstruction.
test_that("23.10 mfp2 serialization preserves prediction behavior", {
  fit <- mfp2(
    lpsa ~ fp(age) + fp(cavol) + svi,
    data = prostate,
    keep = "svi",
    verbose = FALSE
  )
  nd <- prostate[1:15, c("age", "cavol", "svi"), drop = FALSE]
  before <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  after <- predict(restored, newdata = nd, type = "link", se.fit = TRUE)

  expect_s3_class(restored, "mfp2")
  expect_equal(restored$fp_powers, fit$fp_powers)
  expect_equal(coef(restored), coef(fit))
  expect_equal(as.numeric(after$fit), as.numeric(before$fit), tolerance = 1e-12)
  expect_equal(as.numeric(after$se.fit), as.numeric(before$se.fit), tolerance = 1e-12)
})


# Test purpose: Saving and restoring an mfpi object must preserve both ordinary
# subject-level prediction and the manually evaluated fitted-function path.
test_that("23.11 mfpi serialization preserves ordinary and fitted-function predictions", {
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "fp1"),
    flex = "flex1",
    verbose = FALSE
  )
  nd <- prostate[1:12, , drop = FALSE]
  before_link <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  before_fun <- predict(
    fit, terms = "cavol", model = "all", type = "function",
    grid = TRUE, n_grid = 20, se.fit = TRUE
  )

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  after_link <- predict(
    restored, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  after_fun <- predict(
    restored, terms = "cavol", model = "all", type = "function",
    grid = TRUE, n_grid = 20, se.fit = TRUE
  )

  expect_s3_class(restored, "mfpi")
  expect_equal(after_link$predictions, before_link$predictions, tolerance = 1e-12)
  expect_equal(after_fun$functions, before_fun$functions, tolerance = 1e-12)
})


# Test purpose: Continuous-only formula objects created before the explicit
# term and coefficient-column mappings were stored must remain predictable.
# The fallback must also resolve fp()/fp2() formula labels to source variables.
test_that("23.12 legacy continuous formula objects reconstruct prediction mappings", {
  fit <- mfp2(
    lpsa ~ fp(age, center = FALSE) + fp(cavol, center = FALSE) + svi,
    data = prostate,
    keep = c("age", "cavol", "svi"),
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  nd <- prostate[1:18, c("age", "cavol", "svi"), drop = FALSE]
  expected <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
  expected_term <- predict(
    fit, newdata = nd, type = "terms", terms = "age", se.fit = TRUE
  )

  legacy <- fit
  legacy$term_to_columns <- NULL
  legacy$formula_term_to_columns <- NULL
  legacy$transformed_to_model_columns <- NULL
  legacy$formula_prediction_term_names <- NULL

  got <- predict(legacy, newdata = nd, type = "link", se.fit = TRUE)
  got_term <- predict(
    legacy, newdata = nd, type = "terms", terms = "age", se.fit = TRUE
  )

  expect_equal(as.numeric(got$fit), as.numeric(expected$fit), tolerance = 1e-12)
  expect_equal(
    as.numeric(got$se.fit), as.numeric(expected$se.fit), tolerance = 1e-12
  )
  expect_equal(got_term, expected_term, tolerance = 1e-12)
})


# Test purpose: Legacy mapping reconstruction must preserve exact fitted
# coefficient names for matrix columns that require quoting in the final model.
test_that("23.13 legacy matrix objects preserve non-syntactic coefficient names", {
  x <- cbind(
    "age years" = prostate$age,
    "cavol-value" = prostate$cavol
  )
  fit <- mfp2(
    x,
    prostate$lpsa,
    keep = colnames(x),
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  nd <- as.data.frame(x[1:16, , drop = FALSE], check.names = FALSE)
  expected <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
  expected_term <- predict(
    fit, newdata = nd, type = "terms", terms = "age years", se.fit = TRUE
  )

  legacy <- fit
  legacy$term_to_columns <- NULL
  legacy$transformed_to_model_columns <- NULL
  got <- predict(legacy, newdata = nd, type = "link", se.fit = TRUE)
  got_term <- predict(
    legacy, newdata = nd, type = "terms", terms = "age years", se.fit = TRUE
  )

  expect_equal(as.numeric(got$fit), as.numeric(expected$fit), tolerance = 1e-12)
  expect_equal(
    as.numeric(got$se.fit), as.numeric(expected$se.fit), tolerance = 1e-12
  )
  expect_equal(got_term, expected_term, tolerance = 1e-12)
})


# Test purpose: Missing grouped-term metadata cannot be reconstructed as an
# identity mapping. Prediction must fail explicitly rather than silently use a
# conceptual group name as if it were one raw model-matrix column.
test_that("23.14 grouped matrix objects require their member-column mapping", {
  stage <- factor(
    rep(c("I", "II", "III"), length.out = nrow(prostate)),
    levels = c("I", "II", "III")
  )
  mm <- stats::model.matrix(~ stage)[, -1L, drop = FALSE]
  fit <- mfp2(
    mm,
    prostate$lpsa,
    term_groups = list(stage = colnames(mm)),
    keep = "stage",
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  legacy <- fit
  legacy$term_to_columns <- NULL

  expect_error(
    predict(legacy, newdata = mm[1:10, , drop = FALSE]),
    "lacks grouped-term mapping metadata"
  )
})
