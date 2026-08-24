# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Verifies Cox offsets and every covariate-reference choice by
# comparing MFPI directly with the retained coxph model. No manual prediction
# branch is exercised or retained.
test_that("17.1.9b Cox MFPI offsets and references match stored coxph", {
  set.seed(17109)
  n <- 320

  dat <- data.frame(
    age = stats::rnorm(n, mean = 55, sd = 9),
    sex = factor(sample(c("female", "male"), n, replace = TRUE)),
    exposure = stats::runif(n, 0.5, 3)
  )
  sex_effect <- ifelse(dat$sex == "male", 0.35, 0)
  eta <- 0.025 * (dat$age - 55) + sex_effect + log(dat$exposure)
  event_time <- stats::rexp(n, rate = 0.02 * exp(eta))
  censor_time <- stats::rexp(n, rate = 0.012)
  dat$time <- pmin(event_time, censor_time)
  dat$status <- as.integer(event_time <= censor_time)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    offset = log(dat$exposure),
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

  nd <- dat[1:24, c("age", "sex", "exposure"), drop = FALSE]
  new_offset <- log(nd$exposure)
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, newoffset = new_offset
  )
  stored <- fit_result$test_results$interaction_model$fit

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
        newoffset = new_offset,
        cox_reference = reference_value,
        se.fit = TRUE
      )

      expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
      expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
      expect_identical(got$metadata$cox_reference, reference_value)
    }
  }
})


# Test purpose: Verifies formula-interface strata reconstruction. The caller
# supplies only ordinary newdata; predict.mfpi() must recover the original
# strata variable from stored formula metadata and create model_newdata$strata_
# without converting it to integer codes.
test_that("17.1.10 formula-stratified Cox MFPI reconstructs strata from newdata", {
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

  nd <- dat[1:15, c("age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  reconstructed <- reconstruct_formula_strata_newdata(fit, nd)
  design <- mfpi_build_ordinary_design(
    fit,
    "age",
    fit_result,
    nd,
    strata = reconstructed
  )
  stored <- fit_result$test_results$interaction_model$fit

  expect_identical(design$model_newdata$strata_, nd$inst)

  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "link",
    se.fit = TRUE
  )

  expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
  expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
})


# Test purpose: Exercises actual prediction with two Cox strata columns. Matrix
# or data-frame strata must be combined exactly once with strata(),
# whereas a single vector/factor must remain raw for the stored formula to
# evaluate strata(strata_) itself.
test_that("17.1.11 multiple Cox strata columns are combined once in MFPI prediction", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(
    dat[, c("time", "status", "age", "sex", "inst", "ph.ecog")]
  ), ]
  dat$sex <- factor(dat$sex)

  strata_fit <- data.frame(inst = dat$inst, ecog = dat$ph.ecog)
  fit <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = strata_fit,
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

  nd <- dat[1:15, c("sex", "age", "inst", "ph.ecog"), drop = FALSE]
  strata_new <- nd[, c("inst", "ph.ecog"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(
    fit,
    "age",
    fit_result,
    nd,
    strata = strata_new
  )

  expected_strata <- do.call(
    survival::strata,
    c(as.list(strata_new), list(shortlabel = TRUE))
  )
  stored <- fit_result$test_results$interaction_model$fit
  fitted_strata_levels <- stored$xlevels[["strata(strata_)"]]

  # Prediction values must represent the supplied rows, while the factor keeps
  # the complete fit-time level set required by predict.coxph/model.frame().
  expect_s3_class(design$model_newdata$strata_, "factor")
  expect_identical(
    as.character(design$model_newdata$strata_),
    as.character(expected_strata)
  )
  expect_identical(
    levels(design$model_newdata$strata_),
    fitted_strata_levels
  )

  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    reference = "zero"
  )
  got <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "link",
    strata = strata_new,
    se.fit = FALSE
  )

  expect_equal(got$predictions$fit, as.numeric(direct), tolerance = 1e-8)
})


# Test purpose: Verifies clear validation for stratified Cox prediction. A
# stratified stored interaction model cannot predict supplied rows without one
# stratum value/row per prediction row and without missing stratum values.
test_that("17.1.12 MFPI Cox strata validation rejects missing, short, and NA strata", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = dat$inst,
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  nd <- dat[1:8, c("age", "sex", "inst"), drop = FALSE]

  expect_error(
    predict(fit, newdata = nd, terms = "age", model = "all", type = "link"),
    "stratified|strata"
  )
  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all", type = "link",
      strata = nd$inst[-1]
    ),
    "one value or row per prediction row"
  )
  bad_strata <- nd$inst
  bad_strata[1] <- NA
  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all", type = "link",
      strata = bad_strata
    ),
    "must not contain missing values"
  )
  bad_strata_inf <- as.numeric(nd$inst)
  bad_strata_inf[1] <- Inf
  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all", type = "link",
      strata = bad_strata_inf
    ),
    "strata.*finite"
  )
})


# Test purpose: Verifies fitted-function values using direct matrix algebra.
# The manual fitted-function path should return X_g beta_g for every group and
# x value, using the exact stored basis and coefficient mapping.
test_that("17.1.13 MFPI fitted functions equal direct basis-times-coefficient calculations", {
  set.seed(1723)
  n <- 220
  dat <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    x = runif(n, 1, 8)
  )
  dat$y <- 1 + 0.4 * dat$x + 1.0 * dat$x * (dat$group == "B") +
    rnorm(n, sd = 0.2)

  fit <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  eval_x <- data.frame(x = c(1.5, 3, 6))
  pred <- predict(
    fit,
    newdata = eval_x,
    terms = "x",
    model = "all",
    type = "function",
    grid = FALSE,
    se.fit = FALSE
  )

  fit_result <- fit$var_winners[["x"]]$fit
  prepared <- mfpi_prepare_prediction_data(
    object = fit,
    term = "x",
    newdata = eval_x,
    grid = FALSE,
    n_grid = 200L
  )
  basis <- mfpi_build_function_basis(
    object = fit,
    term = "x",
    fit_result = fit_result,
    cont_var_scaled = prepared$cont_var_scaled,
    x_display = prepared$x_display
  )
  coefficients <- fit_result$test_results$interaction_model$coefficients

  group_internal <- names(basis$coefficient_groups)
  group_display <- mfpi_prediction_group_display_labels(
    fit,
    group_internal
  )
  intercept <- if ("(Intercept)" %in% names(coefficients)) {
    unname(coefficients["(Intercept)"])
  } else {
    0
  }

  expected <- do.call(rbind, lapply(seq_along(group_internal), function(i) {
    g <- group_internal[i]
    cols <- basis$coefficient_groups[[g]]
    dummy_name <- if (i > 1L) paste0(fit$group_var, g) else NULL
    dummy_effect <- if (!is.null(dummy_name)) unname(coefficients[dummy_name]) else 0

    data.frame(
      x = prepared$x_display,
      group = group_display[i],
      fit = intercept +
        as.numeric(basis$x[, cols, drop = FALSE] %*% coefficients[cols]) +
        dummy_effect,
      stringsAsFactors = FALSE
    )
  }))

  observed_key <- paste(pred$functions$x, pred$functions$group, sep = "::")
  expected_key <- paste(expected$x, expected$group, sep = "::")
  expected <- expected[match(observed_key, expected_key), , drop = FALSE]

  expect_false(anyNA(expected$fit))
  expect_equal(pred$functions$fit, expected$fit, tolerance = 1e-10)
})


# Test purpose: Verifies that MFPI ordinary-design prediction reconstructs
# vector and tabular Cox strata against the level set of the corresponding
# fitted Cox model.  The tabular case is fitted with tabular strata as well;
# changing the number of stratification variables at prediction time is not a
# valid prediction operation.
test_that("17.1.14 MFPI ordinary design preserves fitted Cox strata levels", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(
    dat[, c("time", "status", "age", "sex", "inst", "ph.ecog")]
  ), ]
  dat$sex <- factor(dat$sex)

  nd <- dat[1:8, c("sex", "age", "inst", "ph.ecog"), drop = FALSE]

  # One-dimensional external strata.
  fit_vector <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = dat$inst,
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  vector_result <- fit_vector$var_winners[["age"]]$fit
  vector_design <- mfpi_build_ordinary_design(
    fit_vector, "age", vector_result, nd, strata = factor(nd$inst)
  )
  vector_stored <- vector_result$test_results$interaction_model$fit

  expect_identical(
    as.character(vector_design$model_newdata$strata_),
    as.character(nd$inst)
  )
  expect_identical(
    levels(vector_design$model_newdata$strata_),
    vector_stored$xlevels[["strata(strata_)"]]
  )

  # Two-dimensional external strata. Fit and predict with the same two
  # stratification variables, mirroring coxph(... + strata(inst, ph.ecog)).
  strata_fit <- dat[, c("inst", "ph.ecog"), drop = FALSE]
  fit_tabular <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    strata = strata_fit,
    flex = "flex1",
    df = 1,
    select = 1,
    alpha = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  tabular_result <- fit_tabular$var_winners[["age"]]$fit
  tabular_strata <- nd[, c("inst", "ph.ecog"), drop = FALSE]
  tabular_design <- mfpi_build_ordinary_design(
    fit_tabular, "age", tabular_result, nd, strata = tabular_strata
  )
  expected_tabular <- do.call(
    survival::strata,
    c(as.list(tabular_strata), list(shortlabel = TRUE))
  )
  tabular_stored <- tabular_result$test_results$interaction_model$fit

  expect_identical(
    as.character(tabular_design$model_newdata$strata_),
    as.character(expected_tabular)
  )
  expect_identical(
    levels(tabular_design$model_newdata$strata_),
    tabular_stored$xlevels[["strata(strata_)"]]
  )
})


test_that("predict.mfpi Cox expected and survival match stored coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
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

  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit

  for (prediction_type in c("expected", "survival")) {
    direct_training <- stats::predict(
      stored, type = prediction_type, se.fit = FALSE
    )
    got_training <- predict(
      fit, terms = "age", model = "all",
      type = prediction_type, se.fit = FALSE
    )
    expect_equal(
      got_training$predictions$fit,
      as.numeric(direct_training),
      tolerance = 1e-8
    )
  }

  nd <- dat[1:16, c("time", "status", "age", "sex"), drop = FALSE]
  response <- reconstruct_cox_prediction_response(
    object = fit,
    fit_obj = stored,
    newdata = nd
  )
  design <- mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  direct_data <- attach_cox_prediction_response(
    fit_obj = stored,
    newdata = design$model_newdata,
    response = response
  )

  for (prediction_type in c("expected", "survival")) {
    direct <- stats::predict(
      stored,
      newdata = direct_data,
      type = prediction_type,
      se.fit = TRUE
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = prediction_type,
      se.fit = TRUE
    )
    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    expect_true(got$metadata$absolute_cox_prediction)
  }

  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all",
      type = "survival", cox_reference = "zero"
    ),
    "applies only"
  )

  nd_missing_response <- nd[, c("age", "sex"), drop = FALSE]
  expect_error(
    predict(
      fit, newdata = nd_missing_response, terms = "age", model = "all",
      type = "survival"
    ),
    "response information"
  )
})


test_that("predict.mfpi matrix-interface Cox survival accepts one Surv column", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    x = dat[, c("sex", "age")],
    y = survival::Surv(dat$time, dat$status),
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

  rows <- 1:12
  nd <- data.frame(
    sex = dat$sex[rows],
    age = dat$age[rows],
    prediction_response = I(survival::Surv(dat$time[rows], dat$status[rows]))
  )
  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  response <- reconstruct_cox_prediction_response(fit, stored, nd)
  design <- mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  direct_data <- attach_cox_prediction_response(
    stored, design$model_newdata, response
  )

  direct <- stats::predict(
    stored, newdata = direct_data, type = "survival", se.fit = FALSE
  )
  got <- predict(
    fit, newdata = nd, terms = "age", model = "all",
    type = "survival", se.fit = FALSE
  )
  expect_equal(got$predictions$fit, as.numeric(direct), tolerance = 1e-8)
})


test_that("predict.mfpi rejects irrelevant Cox strata", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex")]), ]
  dat$sex <- factor(dat$sex)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    group_var = "sex",
    cont_vars = "age",
    cont_var_forms = c(age = "linear"),
    flex = "flex1",
    p_interact = 1,
    verbose = FALSE
  )
  nd <- dat[1:8, c("age", "sex"), drop = FALSE]

  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all",
      type = "lp", strata = rep(1, nrow(nd))
    ),
    "not stratified"
  )
  expect_error(
    predict(
      fit, terms = "age", model = "all",
      type = "function", strata = rep(1, nrow(dat))
    ),
    "not used for MFPI fitted-function"
  )
})


test_that("predict.mfpi stratified Cox survival matches stored coxph", {
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

  nd <- dat[1:14, c("time", "status", "age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  response <- reconstruct_cox_prediction_response(fit, stored, nd)
  strata_new <- reconstruct_formula_strata_newdata(fit, nd)
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = strata_new
  )
  direct_data <- attach_cox_prediction_response(
    stored, design$model_newdata, response
  )

  for (prediction_type in c("expected", "survival")) {
    direct <- stats::predict(
      stored,
      newdata = direct_data,
      type = prediction_type,
      se.fit = TRUE
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = prediction_type,
      se.fit = TRUE
    )
    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
  }
})


test_that("predict.mfpi rejects replacement offsets for models without offsets", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  nd <- prostate[1:7, , drop = FALSE]
  expect_error(
    predict(
      fit,
      newdata = nd,
      terms = "cavol",
      model = "all",
      type = "link",
      newoffset = rep(0, nrow(nd))
    ),
    "fitted with an offset"
  )
})


test_that("predict.mfpi training reconstruction preserves fitted Cox strata", {
  set.seed(19003)
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  dat$inst <- factor(dat$inst)
  dat$exposure <- stats::runif(nrow(dat), 0.7, 2.4)

  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex + offset(log(exposure)) + strata(inst),
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

  fit_result <- fit$var_winners[["age"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  replacement_offset <- log(dat$exposure)

  design <- mfpi_build_ordinary_design(
    fit,
    "age",
    fit_result,
    newdata = NULL,
    newoffset = replacement_offset
  )
  expect_true("strata_" %in% names(design$model_newdata))
  expect_true("offset_" %in% names(design$model_newdata))

  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got <- predict(
    fit,
    terms = "age",
    model = "all",
    type = "lp",
    newoffset = replacement_offset,
    cox_reference = "zero",
    se.fit = TRUE
  )

  expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
  expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)

  # Reconstructing the training rows solely to resupply their fitted strata
  # must also recover the original, uncentred offset scale. The resulting
  # prediction should therefore equal direct prediction on the stored fit.
  fitted_strata <- stored$strata
  direct_training <- stats::predict(
    stored,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got_strata_only <- predict(
    fit,
    terms = "age",
    model = "all",
    type = "lp",
    strata = fitted_strata,
    cox_reference = "zero",
    se.fit = TRUE
  )
  expect_equal(
    got_strata_only$predictions$fit,
    as.numeric(direct_training$fit),
    tolerance = 1e-8
  )
  expect_equal(
    got_strata_only$predictions$se.fit,
    as.numeric(direct_training$se.fit),
    tolerance = 1e-8
  )
})
