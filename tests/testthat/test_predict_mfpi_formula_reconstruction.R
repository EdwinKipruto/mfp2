# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# -----------------------------------------------------------------------------
# 17.2 Grouped categorical adjustment prediction
# -----------------------------------------------------------------------------
# Formula prediction must rebuild the original factor contrasts and pass the
# complete selected adjustment block to the stored interaction model.

# Test purpose: Ensures predict.mfpi() reconstructs a selected factor block
# from ordinary factor-valued formula newdata.
test_that("MFPI prediction reconstructs grouped factor adjustments", {
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

  nd <- dat[1:20, c("trt", "x", "stage", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(
    fit,
    "x",
    fit_result,
    nd
  )

  stored <- fit_result$test_results$interaction_model$fit
  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "link"
  )

  got <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )

  expect_equal(got$predictions$fit, as.numeric(direct), tolerance = 1e-8)
  expect_true(all(c("stageII.1", "stageIII.1") %in%
                    names(design$model_newdata)))
})


# Test purpose: Ensures MFPI formula prediction reconstructs inline factor
# columns while the adjustment model is addressed by the source variable name.
test_that("MFPI prediction reconstructs inline factor adjustments", {
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

  nd <- dat[1:20, c("trt", "x", "stage_code", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(fit, "x", fit_result, nd)

  expect_true(all(c("factor(stage_code)2.1", "factor(stage_code)3.1") %in%
                    names(design$model_newdata)))

  out <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expect_true(all(is.finite(out$predictions$fit)))
})


# Regression: formula ordinary prediction should depend on the selected
# term-specific interaction model, not the complete formula offered to mfpi().
test_that("MFPI formula prediction omits unrelated original predictors", {
  dat <- make_mfpi_factor_data()
  set.seed(17221)
  dat$w <- runif(nrow(dat), 0.5, 4)
  dat$nuisance <- factor(
    rep(
      c("low", "middle", "high", "middle", "high", "low", "high", "low", "middle"),
      length.out = nrow(dat)
    ),
    levels = c("low", "middle", "high")
  )

  fit <- mfpi(
    y ~ trt + x + w + stage + z + nuisance,
    data = dat,
    group_var = "trt",
    cont_vars = c("x", "w"),
    cont_var_forms = c(x = "linear", w = "linear"),
    keep = "stage",
    df = 1,
    select = 0,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  selected <- get_selected_variables(fit$adjustment_model)
  expect_true("stage" %in% selected)
  expect_false(any(c("z", "nuisance") %in% selected))
  expect_identical(
    unname(fit$formula_prediction_term_names[["nuisance"]]),
    "nuisance"
  )

  nd_full <- dat[
    1:18,
    c("trt", "x", "w", "stage", "z", "nuisance"),
    drop = FALSE
  ]
  nd_x <- nd_full[, c("trt", "x", "stage"), drop = FALSE]
  nd_w <- nd_full[, c("trt", "w", "stage"), drop = FALSE]

  pred_x_full <- predict(
    fit, newdata = nd_full, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )
  pred_x_minimal <- predict(
    fit, newdata = nd_x, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )
  pred_w_full <- predict(
    fit, newdata = nd_full, terms = "w", model = "all",
    type = "link", se.fit = TRUE
  )
  pred_w_minimal <- predict(
    fit, newdata = nd_w, terms = "w", model = "all",
    type = "link", se.fit = TRUE
  )

  expect_equal(
    pred_x_minimal$predictions,
    pred_x_full$predictions,
    tolerance = 1e-10
  )
  expect_equal(
    pred_w_minimal$predictions,
    pred_w_full$predictions,
    tolerance = 1e-10
  )

  nd_unseen <- nd_full
  nd_unseen$z <- NA_real_
  nd_unseen$nuisance <- factor(
    rep("unseen", nrow(nd_unseen)),
    levels = c(levels(dat$nuisance), "unseen")
  )
  pred_x_unseen <- predict(
    fit, newdata = nd_unseen, terms = "x", model = "all",
    type = "link", se.fit = TRUE
  )
  expect_equal(
    pred_x_unseen$predictions,
    pred_x_minimal$predictions,
    tolerance = 1e-10
  )
})


# Regression: the formula-label map must allow a selected inline factor block
# to be reconstructed without evaluating an unrelated eliminated inline factor.
test_that("MFPI minimal prediction supports selected inline factor adjustments", {
  dat <- make_mfpi_factor_data()
  dat$stage_code <- as.integer(dat$stage)
  dat$nuisance_code <- rep(
    c(1L, 2L, 3L, 2L, 3L, 1L, 3L, 1L, 2L),
    length.out = nrow(dat)
  )

  fit <- mfpi(
    y ~ trt + x + factor(stage_code) + factor(nuisance_code) + z,
    data = dat,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "stage_code",
    df = 1,
    select = 0,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(
    unname(fit$formula_prediction_term_names[["factor(stage_code)"]]),
    "stage_code"
  )
  expect_identical(
    unname(fit$formula_prediction_term_names[["factor(nuisance_code)"]]),
    "nuisance_code"
  )
  expect_true("stage_code" %in% get_selected_variables(fit$adjustment_model))
  expect_false("nuisance_code" %in% get_selected_variables(fit$adjustment_model))

  nd_full <- dat[
    1:20,
    c("trt", "x", "stage_code", "nuisance_code", "z"),
    drop = FALSE
  ]
  nd_minimal <- nd_full[, c("trt", "x", "stage_code"), drop = FALSE]
  nd_unseen <- nd_full
  nd_unseen$nuisance_code <- 99L

  pred_full <- predict(
    fit, newdata = nd_full, terms = "x", model = "all",
    type = "link", se.fit = FALSE
  )
  pred_minimal <- predict(
    fit, newdata = nd_minimal, terms = "x", model = "all",
    type = "link", se.fit = FALSE
  )
  pred_unseen <- predict(
    fit, newdata = nd_unseen, terms = "x", model = "all",
    type = "link", se.fit = FALSE
  )

  expect_equal(pred_minimal$predictions, pred_full$predictions,
               tolerance = 1e-10)
  expect_equal(pred_unseen$predictions, pred_minimal$predictions,
               tolerance = 1e-10)
  expect_true(all(c("factor(stage_code)2.1", "factor(stage_code)3.1") %in%
                    pred_minimal$metadata$model_newdata_columns))
  expect_false(any(grepl(
    "nuisance_code",
    pred_minimal$metadata$model_newdata_columns,
    fixed = TRUE
  )))
})
