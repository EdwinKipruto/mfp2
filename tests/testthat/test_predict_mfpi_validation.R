# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Ensures predict.mfpi() validates se.fit as a single non-missing
# logical value.
test_that("predict.mfpi() rejects invalid se.fit", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(fit, terms = "cavol", type = "function", model = "all", se.fit = NA),
    "`se.fit` must be"
  )
})


# Test purpose: Ensures predict.mfpi() validates confidence level as a single
# numeric value in (0, 1).
test_that("predict.mfpi() rejects invalid confidence level", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(fit, terms = "cavol", type = "function", model = "all", level = 1),
    "`level` must be"
  )
})


# Test purpose: Ensures predict.mfpi() fails clearly when a requested term has
# no stored MFPI interaction model in the requested model scope.
test_that("predict.mfpi() rejects unknown prediction terms", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(fit, terms = "age", type = "function", model = "all"),
    "requested terms"
  )
})


# Test purpose: Ensures fitted-function prediction checks that newdata contains
# the requested continuous variable.
test_that("predict.mfpi() requires fitted-function newdata to contain requested term", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_error(
    predict(
      fit,
      terms = "cavol",
      type = "function",
      model = "all",
      newdata = data.frame(age = prostate$age[1:10])
    ),
    "must contain a column named `cavol`"
  )
})


# Test purpose: Checks ordinary subject-level link-scale MFPI prediction from
# a term-specific interaction model using supplied newdata.
test_that("predict.mfpi() type = 'link' returns subject-level predictions", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(
    fit,
    terms = "cavol",
    type = "link",
    model = "all",
    newdata = prostate[1:10, ],
    se.fit = FALSE
  )

  expect_s3_class(p, "mfpi_prediction")
  expect_true(!is.null(p$predictions))
  expect_equal(nrow(p$predictions), 10)
  expect_true(all(is.finite(p$predictions$fit)))
})


# Test purpose: Checks ordinary subject-level response-scale MFPI prediction
# using supplied newdata.
test_that("predict.mfpi() type = 'response' returns subject-level predictions", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(
    fit,
    terms = "cavol",
    type = "response",
    model = "all",
    newdata = prostate[1:10, ],
    se.fit = FALSE
  )

  expect_s3_class(p, "mfpi_prediction")
  expect_true(!is.null(p$predictions))
  expect_equal(nrow(p$predictions), 10)
  expect_true(all(is.finite(p$predictions$fit)))
})


# Test purpose: Ensures fitted-function prediction with grid = TRUE returns
# values on the requested evaluation grid.
test_that("predict.mfpi() fitted-function grid uses requested n_grid", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(
    fit,
    terms = "cavol",
    type = "function",
    model = "all",
    grid = TRUE,
    n_grid = 25,
    se.fit = FALSE
  )

  expect_s3_class(p, "mfpi_prediction")
  expect_true(!is.null(p$functions))
  expect_true(length(unique(p$functions$x)) <= 25)
  expect_true(all(is.finite(p$functions$fit)))
})


# Test purpose: Ensures grid = TRUE is ignored with a warning for ordinary
# subject-level MFPI prediction types.
test_that("predict.mfpi() warns when grid is used with link prediction", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  expect_warning(
    predict(
      fit,
      terms = "cavol",
      type = "link",
      model = "all",
      newdata = prostate[1:5, ],
      grid = TRUE,
      se.fit = FALSE
    ),
    "grid"
  )
})


test_that("predict.mfpi enforces family-specific arguments and types", {
  data("prostate", package = "mfp2")
  fit_glm <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )
  nd_glm <- prostate[1:8, , drop = FALSE]

  expect_error(
    predict(fit_glm, newdata = nd_glm, terms = "cavol", model = "all",
            type = "risk"),
    "GLM MFPI models"
  )
  expect_error(
    predict(fit_glm, newdata = nd_glm, terms = "cavol", model = "all",
            type = "link", cox_reference = "zero"),
    "available only for Cox"
  )
  expect_error(
    predict(fit_glm, newdata = nd_glm, terms = "cavol", model = "all",
            type = "link", strata = rep(1, nrow(nd_glm))),
    "available only for Cox"
  )
  expect_error(
    predict(fit_glm, terms = "cavol", model = "all", type = "function",
            strata = rep(1, nrow(prostate))),
    "not used for MFPI fitted-function"
  )

  nd_bad <- nd_glm
  nd_bad$cavol[1] <- Inf
  expect_error(
    predict(fit_glm, newdata = nd_bad, terms = "cavol", model = "all",
            type = "link"),
    "finite"
  )
})


test_that("predict.mfpi Cox references match the retained coxph model", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  dat$inst <- factor(dat$inst)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    ties = "breslow",
    verbose = FALSE
  )

  nd <- dat[1:18, c("age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  strata_new <- reconstruct_formula_strata_newdata(fit, nd)
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = strata_new
  )

  for (prediction_type in c("lp", "risk")) {
    for (reference_value in c("zero", "sample", "strata")) {
      direct <- stats::predict(
        stored,
        newdata = design$model_newdata,
        type = prediction_type,
        se.fit = TRUE,
        reference = reference_value
      )
      got <- predict(
        fit,
        newdata = nd,
        terms = "age",
        model = "all",
        type = prediction_type,
        cox_reference = reference_value,
        se.fit = TRUE
      )
      expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
      expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    }
  }
})


test_that("predict.mfpi validates required predictors but ignores irrelevant extras", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ age + svi + cavol,
    data = prostate,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    keep = "age",
    p_interact = 1,
    verbose = FALSE
  )
  nd <- prostate[1:8, c("age", "svi", "cavol"), drop = FALSE]
  expected <- predict(
    fit, newdata = nd, terms = "cavol", model = "all", type = "link"
  )
  nd_extra <- nd
  nd_extra$unused_na <- NA_real_
  nd_extra$unused_nan <- NaN
  expect_equal(
    predict(
      fit, newdata = nd_extra, terms = "cavol", model = "all", type = "link"
    )$predictions,
    expected$predictions,
    tolerance = 1e-12
  )

  for (column in c("cavol", "age", "svi")) {
    for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
      nd_bad <- nd
      nd_bad[[column]][2] <- bad_value
      expect_error(
        predict(
          fit, newdata = nd_bad, terms = "cavol", model = "all", type = "link"
        ),
        "finite|missing|NA"
      )
    }
  }
})


test_that("formula-derived offsets reject missing and non-finite inputs", {
  data("prostate", package = "mfp2")
  set.seed(21002)
  dat <- prostate
  dat$exposure <- stats::runif(nrow(dat), 0.7, 2.2)

  fit_mfp2 <- mfp2(
    lpsa ~ age + svi + offset(log(exposure)),
    data = dat,
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE
  )
  fit_mfpi <- mfpi(
    lpsa ~ age + svi + cavol + offset(log(exposure)),
    data = dat,
    group_var = "svi",
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    keep = "age",
    p_interact = 1,
    verbose = FALSE
  )
  nd <- dat[1:7, c("age", "svi", "cavol", "exposure"), drop = FALSE]

  for (bad_value in list(NA_real_, NaN, Inf, -Inf)) {
    nd_bad <- nd
    nd_bad$exposure[2] <- bad_value
    expect_error(
      predict(fit_mfp2, newdata = nd_bad, type = "link"),
      "offset|missing|finite"
    )
    expect_error(
      predict(
        fit_mfpi,
        newdata = nd_bad,
        terms = "cavol",
        model = "all",
        type = "link"
      ),
      "offset|missing|finite"
    )
  }
})
