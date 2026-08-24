# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: In the simplest two-group linear case, MFPI should fit the same
# interaction model as an ordinary Gaussian model y ~ group * x. The comparison
# uses fitted values, log-likelihood, and newdata predictions, avoiding reliance
# on package-specific coefficient names.
test_that("17.1.5 two-group linear MFPI equals an explicit Gaussian interaction model", {
  set.seed(1715)
  n <- 240
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2)),
    x = runif(n, 1, 8)
  )
  dat$y <- 1.2 +
    0.7 * (dat$group == "treated") +
    0.4 * dat$x +
    1.1 * dat$x * (dat$group == "treated") +
    rnorm(n, sd = 0.15)

  fit_mfpi <- mfpi(
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
    xorder = "original",
    p_interact = 1,
    verbose = FALSE
  )

  fit_reference <- stats::glm(y ~ group * x, data = dat)
  stored <- fit_mfpi$var_winners[["x"]]$fit$test_results$interaction_model$fit

  expect_equal(
    unname(stats::fitted(stored)),
    unname(stats::fitted(fit_reference)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(stored)),
    as.numeric(stats::logLik(fit_reference)),
    tolerance = 1e-8
  )

  nd <- dat[c(1, 30, 121, 180), c("group", "x"), drop = FALSE]
  got <- predict(
    fit_mfpi,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = FALSE
  )
  expected <- stats::predict(fit_reference, newdata = nd, type = "link")

  expect_equal(got$predictions$fit, unname(expected), tolerance = 1e-8)
})


# Test purpose: Extends the explicit interaction oracle to three groups. This
# detects incorrect K - 1 dummy construction, swapped group-specific slopes,
# and hard-coded assumptions that only two treatment groups exist.
test_that("17.1.6 three-group linear MFPI equals an explicit Gaussian interaction model", {
  set.seed(1716)
  n_per_group <- 100
  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = n_per_group)),
    x = runif(3 * n_per_group, 1, 8)
  )
  intercept_shift <- c(A = 0, B = 0.6, C = -0.4)[as.character(dat$group)]
  slope_shift <- c(A = 0, B = 0.8, C = -0.5)[as.character(dat$group)]
  dat$y <- 1 + intercept_shift + (0.5 + slope_shift) * dat$x +
    rnorm(nrow(dat), sd = 0.15)

  fit_mfpi <- mfpi(
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
    xorder = "original",
    p_interact = 1,
    verbose = FALSE
  )

  fit_reference <- stats::glm(y ~ group * x, data = dat)
  stored <- fit_mfpi$var_winners[["x"]]$fit$test_results$interaction_model$fit

  expect_equal(
    unname(stats::fitted(stored)),
    unname(stats::fitted(fit_reference)),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(stats::logLik(stored)),
    as.numeric(stats::logLik(fit_reference)),
    tolerance = 1e-8
  )

  nd <- dat[c(1, 101, 201), c("group", "x"), drop = FALSE]
  got <- predict(
    fit_mfpi,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "response",
    se.fit = FALSE
  )
  expected <- stats::predict(fit_reference, newdata = nd, type = "response")

  expect_equal(got$predictions$fit, unname(expected), tolerance = 1e-8)
})


# Test purpose: Verifies that a strong interaction produces the expected MFPI
# test result. The test does not depend on a borderline random p-value: the data
# use a large slope difference and low noise, so failure indicates a structural
# interaction-test regression.
test_that("17.1.7 MFPI detects a strong prespecified linear interaction", {
  set.seed(1717)
  n <- 260
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2)),
    x = runif(n, 1, 8)
  )
  dat$y <- 1 + 0.3 * dat$x + 2.0 * dat$x * (dat$group == "treated") +
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
    p_interact = 0.05,
    verbose = FALSE
  )

  metric <- fit$all_model_metrics[fit$all_model_metrics$variable == "x", , drop = FALSE]
  expect_equal(nrow(metric), 1)
  expect_true(is.finite(metric$pvalue))
  expect_lt(metric$pvalue, 0.05)
  expect_true("x" %in% names(fit$best_interaction_model))
})


# Test purpose: Verifies Poisson MFPI ordinary prediction when the fitted
# interaction model uses an offset. The reconstructed offset_ column must be
# consumed by predict.glm(), and both link and response predictions must match
# direct prediction from the stored interaction model.
test_that("17.1.8 Poisson MFPI offset predictions match the stored glm", {
  set.seed(1718)
  n <- 260
  dat <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    x = runif(n, 1, 6),
    exposure = runif(n, 0.5, 3)
  )
  eta <- 0.2 + 0.12 * dat$x + 0.35 * (dat$group == "B") +
    0.18 * dat$x * (dat$group == "B") + log(dat$exposure)
  dat$y <- rpois(n, exp(eta))

  fit <- mfpi(
    dat[, c("group", "x")],
    dat$y,
    family = "poisson",
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    offset = log(dat$exposure),
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

  nd <- dat[1:18, c("group", "x", "exposure"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  design <- mfpi_build_ordinary_design(
    fit,
    "x",
    fit_result,
    nd,
    newoffset = log(nd$exposure)
  )
  stored <- fit_result$test_results$interaction_model$fit

  expect_true("offset_" %in% names(design$model_newdata))
  expect_equal(design$model_newdata$offset_, log(nd$exposure))

  # Build the interaction design directly from the stored formula. The offset
  # is not a coefficient column; it is added to X beta after multiplication.
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

  manual_offset <- stats::model.offset(reference_frame)
  expect_equal(manual_offset, log(nd$exposure))

  manual_link <- as.numeric(manual_x %*% beta + manual_offset)
  manual_link_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_link_se <- as.numeric(sqrt(pmax(manual_link_variance, 0)))
  manual_response <- as.numeric(stored$family$linkinv(manual_link))

  # predict.glm(type = "response", se.fit = TRUE) applies the delta method:
  # response-scale SE = link-scale SE * abs(d mu / d eta).
  manual_response_se <- manual_link_se * abs(stored$family$mu.eta(manual_link))

  for (prediction_type in c("link", "response")) {
    direct <- stats::predict(
      stored,
      newdata = design$model_newdata,
      type = prediction_type,
      se.fit = TRUE
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "x",
      model = "all",
      type = prediction_type,
      newoffset = log(nd$exposure),
      se.fit = TRUE
    )

    expected_fit <- if (prediction_type == "link") manual_link else manual_response
    expected_se <- if (prediction_type == "link") manual_link_se else manual_response_se

    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    expect_equal(got$predictions$fit, expected_fit, tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, expected_se, tolerance = 1e-8)
  }
})


# Test purpose: Verifies native Cox ordinary prediction types without strata.
# "link" is accepted as an alias for "lp"; "response" is not an alias for
# "risk". Every result is delegated to the stored coxph model.
test_that("17.1.9 unstratified Cox MFPI lp and risk match stored coxph", {
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

  nd <- dat[1:15, c("age", "sex"), drop = FALSE]
  fit_result <- fit$var_winners[["age"]]$fit
  design <- mfpi_build_ordinary_design(fit, "age", fit_result, nd)
  stored <- fit_result$test_results$interaction_model$fit

  cases <- list(
    lp = "lp",
    link = "lp",
    named_link = c(alias = "link"),
    risk = "risk"
  )
  for (case_name in names(cases)) {
    mfpi_type <- cases[[case_name]]
    native_type <- if (identical(unname(mfpi_type), "link")) {
      "lp"
    } else {
      unname(mfpi_type)
    }
    direct <- stats::predict(
      stored,
      newdata = design$model_newdata,
      type = native_type,
      se.fit = TRUE,
      reference = "zero"
    )
    got <- predict(
      fit,
      newdata = nd,
      terms = "age",
      model = "all",
      type = mfpi_type,
      cox_reference = "zero",
      se.fit = TRUE
    )

    expect_equal(got$predictions$fit, as.numeric(direct$fit), tolerance = 1e-8)
    expect_equal(got$predictions$se.fit, as.numeric(direct$se.fit), tolerance = 1e-8)
    expect_true(got$metadata$used_model_predict)
  }

  expect_error(
    predict(
      fit, newdata = nd, terms = "age", model = "all",
      type = "response"
    ),
    "Cox MFPI models"
  )

  # Independent zero-reference oracle for the linear predictor.
  reference_terms <- stats::delete.response(stats::terms(stored))
  manual_x <- stats::model.matrix(
    reference_terms,
    data = design$model_newdata,
    contrasts.arg = stored$contrasts
  )
  manual_x <- manual_x[, colnames(manual_x) != "(Intercept)", drop = FALSE]
  beta <- stats::coef(stored)
  beta_vcov <- stats::vcov(stored)
  manual_lp <- as.numeric(manual_x %*% beta)
  manual_variance <- rowSums((manual_x %*% beta_vcov) * manual_x)
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  got_lp <- predict(
    fit,
    newdata = nd,
    terms = "age",
    model = "all",
    type = "lp",
    cox_reference = "zero",
    se.fit = TRUE
  )
  expect_equal(got_lp$predictions$fit, manual_lp, tolerance = 1e-8)
  expect_equal(got_lp$predictions$se.fit, manual_se, tolerance = 1e-8)
})
