# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 17. predict.mfpi()
# =============================================================================

# Test purpose: Checks that MFPI fitted-function predictions are returned in
# the expected structure.
test_that("predict.mfpi() returns fitted-function predictions", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(fit, terms = "cavol", type = "function", model = "all")

  if (inherits(p, "mfpi_prediction")) {
    expect_true(!is.null(p$functions))
    expect_true(is.data.frame(p$functions))
  } else if (is.list(p)) {
    # May be a list when wrapping
    expect_true(length(p) >= 1)
  }
})


# Test purpose: Checks that MFPI prediction can return both fitted functions
# and group differences.
test_that("predict.mfpi() type = 'both' returns functions and differences", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  p <- predict(fit, terms = "cavol", type = "both", model = "all")

  if (inherits(p, "mfpi_prediction")) {
    expect_true(!is.null(p$functions))
    expect_true(!is.null(p$differences))
  }
})


# Test purpose: Ordinary Gaussian MFPI newdata prediction delegates to the
# stored interaction glm using the reconstructed formula-compatible data frame.
test_that("predict.mfpi ordinary Gaussian prediction matches stored glm", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )
  nd <- prostate[1:12, , drop = FALSE]
  fit_result <- fit$var_winners[["cavol"]]$fit
  design <- mfpi_build_ordinary_design(
    fit, "cavol", fit_result, nd
  )
  direct_link <- stats::predict(
    fit_result$test_results$interaction_model$fit,
    newdata = design$model_newdata,
    type = "link",
    se.fit = TRUE
  )
  got_link <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  expect_equal(got_link$predictions$fit, as.numeric(direct_link$fit))
  expect_equal(got_link$predictions$se.fit, as.numeric(direct_link$se.fit))
  expect_true(got_link$metadata$used_model_predict)

  direct_response <- stats::predict(
    fit_result$test_results$interaction_model$fit,
    newdata = design$model_newdata,
    type = "response",
    se.fit = FALSE
  )
  got_response <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "response", se.fit = FALSE
  )
  expect_equal(got_response$predictions$fit, as.numeric(direct_response))

  # Independent oracle: reconstruct the interaction model matrix and calculate
  # eta = X beta and sqrt(diag(X V X')) directly. This avoids relying solely on
  # predict.glm(), which is also used internally by ordinary MFPI prediction.
  stored <- fit_result$test_results$interaction_model$fit
  reference_terms <- stats::delete.response(stats::terms(stored))
  reference_frame <- stats::model.frame(
    reference_terms,
    data = design$model_newdata,
    xlev = stored$xlevels,
    na.action = stats::na.pass
  )
  manual_x <- stats::model.matrix(
    reference_terms,
    data = reference_frame,
    contrasts.arg = stored$contrasts
  )
  beta <- stats::coef(stored)
  beta_vcov <- stats::vcov(stored)
  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  manual_link <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  expect_equal(got_link$predictions$fit, manual_link, tolerance = 1e-8)
  expect_equal(got_link$predictions$se.fit, manual_se, tolerance = 1e-8)
  expect_equal(got_response$predictions$fit, manual_link, tolerance = 1e-8)
})


# Test purpose: Stratified Cox MFPI prediction reconstructs prediction strata
# with the fitted Cox level set and matches predict.coxph(reference = "zero").
test_that("predict.mfpi stratified Cox prediction matches stored coxph", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]
  dat$sex <- factor(dat$sex)
  fit <- mfpi(
    survival::Surv(time, status) ~ age + sex,
    data = dat,
    family = "cox",
    cont_vars = "age",
    group_var = "sex",
    cont_var_forms = c(age = "linear"),
    strata = dat$inst,
    p_interact = 0.95,
    verbose = FALSE
  )
  nd <- dat[1:10, c("age", "sex", "inst"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(
    fit, "age", fit_result, nd, strata = nd$inst
  )
  stored <- fit_result$test_results$interaction_model$fit
  fitted_strata_levels <- stored$xlevels[["strata(strata_)"]]

  expect_true(is.factor(design$model_newdata$strata_))
  expect_identical(levels(design$model_newdata$strata_), fitted_strata_levels)
  expect_identical(
    as.character(design$model_newdata$strata_),
    as.character(nd$inst)
  )

  direct <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "lp",
    se.fit = TRUE,
    reference = "zero"
  )
  got <- predict(
    fit, newdata = nd, terms = "age", model = "all",
    type = "link", strata = nd$inst, se.fit = TRUE
  )
  expect_equal(got$predictions$fit, as.numeric(direct$fit))
  expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit))
  expect_true(got$metadata$used_model_predict)
})


test_that("predict.mfpi defaults safely and normalizes GLM lp", {
  data("prostate", package = "mfp2")
  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    group_var = "svi",
    flex = "flex1",
    verbose = FALSE
  )

  default_pred <- predict(fit, terms = "cavol", model = "all")
  expect_s3_class(default_pred, "mfpi_prediction")
  expect_identical(default_pred$type, "both")
  expect_true(!is.null(default_pred$functions))
  expect_true(!is.null(default_pred$differences))

  nd <- prostate[1:10, , drop = FALSE]
  link_pred <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "link", se.fit = TRUE
  )
  lp_pred <- predict(
    fit, newdata = nd, terms = "cavol", model = "all",
    type = "lp", se.fit = TRUE
  )
  expect_identical(lp_pred$type, "link")
  expect_equal(lp_pred$predictions, link_pred$predictions)
})
