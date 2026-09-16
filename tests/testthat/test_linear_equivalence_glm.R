# Complete forced-linear GLM equivalence --------------------------------------
#
# One canonical formula exercises all likelihood-based families implemented by
# base stats::glm(). Selection and transformation complexity is disabled, so an
# mfp2 df = 1 fit must be the same model as the direct glm() fit. The checks
# cover fitted quantities, likelihood metrics, every public complete-model
# prediction scale (including standard errors), an independent X beta / X V X'
# oracle, and the package-owned term and contrast paths.


expect_complete_glm_linear_equivalence <- function(fit_mfp2,
                                                   fit_glm,
                                                   newdata,
                                                   mfp2_newdata = newdata,
                                                   mfp2_predict_args = list(),
                                                   term_newdata = mfp2_newdata,
                                                   info,
                                                   tolerance = 1e-7) {
  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "glm")
  expect_identical(fit_mfp2$fitter, "base", info = info)
  expect_identical(
    fit_mfp2$family_string,
    fit_glm$family$family,
    info = info
  )
  expect_identical(fit_mfp2$family$link, fit_glm$family$link, info = info)
  expect_true(isTRUE(fit_mfp2$convergence_mfp), info = info)
  expect_true(isTRUE(fit_mfp2$converged), info = info)
  expect_true(isTRUE(fit_glm$converged), info = info)

  expect_mfp2_glm_parameters_equal(
    fit_mfp2,
    fit_glm,
    tolerance = tolerance
  )
  expect_equal(
    unname(fit_mfp2$linear.predictors),
    unname(fit_glm$linear.predictors),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(stats::fitted(fit_mfp2)),
    unname(stats::fitted(fit_glm)),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$rank, fit_glm$rank, info = info)
  expect_equal(fit_mfp2$df.residual, fit_glm$df.residual, info = info)
  expect_equal(fit_mfp2$iter, fit_glm$iter, info = info)
  expect_equal(
    unname(fit_mfp2$prior.weights),
    unname(fit_glm$prior.weights),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    unname(fit_mfp2$weights),
    unname(fit_glm$weights),
    tolerance = tolerance,
    info = info
  )

  for (residual_type in c("deviance", "pearson", "working", "response")) {
    expect_equal(
      unname(stats::residuals(fit_mfp2, type = residual_type)),
      unname(stats::residuals(fit_glm, type = residual_type)),
      tolerance = tolerance,
      info = paste(info, residual_type, "residuals")
    )
  }

  loglik_mfp2 <- stats::logLik(fit_mfp2)
  loglik_glm <- stats::logLik(fit_glm)
  expect_equal(
    as.numeric(loglik_mfp2),
    as.numeric(loglik_glm),
    tolerance = tolerance,
    info = info
  )
  expect_equal(attr(loglik_mfp2, "df"), attr(loglik_glm, "df"), info = info)
  expect_equal(
    attr(loglik_mfp2, "nobs"),
    attr(loglik_glm, "nobs"),
    info = info
  )
  expect_equal(stats::nobs(fit_mfp2), stats::nobs(fit_glm), info = info)
  expect_equal(
    stats::AIC(fit_mfp2),
    stats::AIC(fit_glm),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    stats::deviance(fit_mfp2),
    stats::deviance(fit_glm),
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$null.deviance,
    fit_glm$null.deviance,
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$aic, fit_glm$aic, tolerance = tolerance, info = info)

  # mfp2 metadata must describe the same forced-linear model. For GLMs the
  # package's model deviance is the native residual deviance, while model df is
  # the likelihood df (including an estimated dispersion where applicable).
  expect_equal(
    fit_mfp2$mfp_logl,
    as.numeric(loglik_glm),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$mfp_df, attr(loglik_glm, "df"), info = info)
  expect_equal(
    fit_mfp2$mfp_deviance,
    fit_glm$deviance,
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$null_deviance,
    fit_glm$null.deviance,
    tolerance = tolerance,
    info = info
  )
  expect_equal(
    fit_mfp2$linear_logl,
    as.numeric(loglik_glm),
    tolerance = tolerance,
    info = info
  )
  expect_equal(fit_mfp2$linear_df, attr(loglik_glm, "df"), info = info)
  expect_equal(
    fit_mfp2$linear_deviance,
    fit_glm$deviance,
    tolerance = tolerance,
    info = info
  )
  expect_true(is.na(fit_mfp2$null_logl), info = info)

  check_predictions <- function(reference_data = NULL,
                                mfp2_data = reference_data,
                                label) {
    if (is.null(reference_data)) {
      reference_frame <- stats::model.frame(fit_glm)
      manual_x <- stats::model.matrix(fit_glm)
      mfp2_link <- stats::predict(fit_mfp2, type = "link", se.fit = TRUE)
      glm_link <- stats::predict(fit_glm, type = "link", se.fit = TRUE)
      mfp2_response <- stats::predict(
        fit_mfp2,
        type = "response",
        se.fit = TRUE
      )
      glm_response <- stats::predict(
        fit_glm,
        type = "response",
        se.fit = TRUE
      )
      default_prediction <- stats::predict(fit_mfp2)
      lp_alias <- stats::predict(fit_mfp2, type = "lp")
    } else {
      reference_terms <- stats::delete.response(stats::terms(fit_glm))
      reference_frame <- stats::model.frame(
        reference_terms,
        data = reference_data,
        xlev = fit_glm$xlevels,
        na.action = stats::na.pass
      )
      manual_x <- stats::model.matrix(
        reference_terms,
        data = reference_frame,
        contrasts.arg = fit_glm$contrasts
      )
      mfp2_link <- do.call(
        stats::predict,
        c(
          list(
            object = fit_mfp2,
            newdata = mfp2_data,
            type = "link",
            se.fit = TRUE
          ),
          mfp2_predict_args
        )
      )
      glm_link <- stats::predict(
        fit_glm,
        newdata = reference_data,
        type = "link",
        se.fit = TRUE
      )
      mfp2_response <- do.call(
        stats::predict,
        c(
          list(
            object = fit_mfp2,
            newdata = mfp2_data,
            type = "response",
            se.fit = TRUE
          ),
          mfp2_predict_args
        )
      )
      glm_response <- stats::predict(
        fit_glm,
        newdata = reference_data,
        type = "response",
        se.fit = TRUE
      )
      default_prediction <- do.call(
        stats::predict,
        c(
          list(object = fit_mfp2, newdata = mfp2_data),
          mfp2_predict_args
        )
      )
      lp_alias <- do.call(
        stats::predict,
        c(
          list(object = fit_mfp2, newdata = mfp2_data, type = "lp"),
          mfp2_predict_args
        )
      )
    }

    beta <- stats::coef(fit_glm)
    covariance <- stats::vcov(fit_glm)
    manual_x <- manual_x[, names(beta), drop = FALSE]
    manual_offset <- stats::model.offset(reference_frame)
    if (is.null(manual_offset)) manual_offset <- rep(0, nrow(manual_x))

    manual_link <- as.numeric(manual_x %*% beta + manual_offset)
    manual_variance <- unname(rowSums((manual_x %*% covariance) * manual_x))
    manual_link_se <- sqrt(pmax(manual_variance, 0))
    manual_response <- as.numeric(fit_glm$family$linkinv(manual_link))
    manual_response_se <- abs(fit_glm$family$mu.eta(manual_link)) *
      manual_link_se

    expect_native_prediction_equal(
      mfp2_link,
      glm_link,
      tolerance = tolerance,
      info = paste(info, label, "link")
    )
    expect_native_prediction_equal(
      mfp2_response,
      glm_response,
      tolerance = tolerance,
      info = paste(info, label, "response")
    )
    expect_equal(
      unname(mfp2_link$fit),
      manual_link,
      tolerance = tolerance,
      info = paste(info, label, "manual link")
    )
    expect_equal(
      unname(mfp2_link$se.fit),
      manual_link_se,
      tolerance = tolerance,
      info = paste(info, label, "manual link SE")
    )
    expect_equal(
      unname(mfp2_response$fit),
      manual_response,
      tolerance = tolerance,
      info = paste(info, label, "manual response")
    )
    expect_equal(
      unname(mfp2_response$se.fit),
      manual_response_se,
      tolerance = tolerance,
      info = paste(info, label, "manual response SE")
    )
    expect_equal(
      unname(default_prediction),
      manual_link,
      tolerance = tolerance,
      info = paste(info, label, "default prediction")
    )
    expect_equal(
      unname(lp_alias),
      manual_link,
      tolerance = tolerance,
      info = paste(info, label, "lp alias")
    )
  }

  check_predictions(label = "training")
  check_predictions(newdata, mfp2_newdata, label = "newdata")
  expect_linear_term_and_contrast_equal(
    object = fit_mfp2,
    newdata = term_newdata,
    term = "x1",
    ref = 0.25,
    has_intercept = TRUE,
    tolerance = tolerance
  )
}


test_that("df = 1 likelihood GLMs match base glm completely", {
  set.seed(4101)
  n <- 480L
  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    off = stats::rnorm(n, mean = 0.1, sd = 0.18),
    case_weight = sample(c(1, 2), n, replace = TRUE)
  )
  group_effect <- c(A = 0, B = 0.35, C = -0.25)[as.character(dat$group)]
  eta <- 0.35 + 0.30 * dat$x1 - 0.20 * dat$x2 + group_effect + dat$off
  mu <- exp(eta)

  responses <- list(
    gaussian = eta + stats::rnorm(n, sd = 0.45),
    binomial = stats::rbinom(n, size = 1L, prob = stats::plogis(eta - 0.5)),
    poisson = stats::rpois(n, lambda = mu),
    Gamma = stats::rgamma(n, shape = 8, scale = mu / 8),
    inverse.gaussian = mu * exp(stats::rnorm(n, sd = 0.12))
  )
  families <- list(
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    Gamma = stats::Gamma(link = "log"),
    inverse.gaussian = stats::inverse.gaussian(link = "log")
  )
  newdata <- dat[seq_len(36L), c("x1", "x2", "group", "off"), drop = FALSE]

  for (family_name in names(families)) {
    family <- families[[family_name]]
    dat$y <- responses[[family_name]]

    fit_mfp2 <- mfp2(
      y ~ x1 + x2 + group,
      data = dat,
      family = family,
      weights = case_weight,
      offset = off,
      cycles = 1,
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      xorder = "original",
      verbose = FALSE
    )
    fit_glm <- stats::glm(
      y ~ x1 + x2 + group + offset(off),
      data = dat,
      family = family,
      weights = case_weight,
      x = TRUE,
      y = TRUE,
      model = TRUE
    )

    expect_equal(
      sort(get_selected_variable_names(fit_mfp2)),
      sort(c("x1", "x2", "group")),
      info = family_name
    )
    expect_complete_glm_linear_equivalence(
      fit_mfp2,
      fit_glm,
      newdata = newdata,
      mfp2_predict_args = list(newoffset = newdata$off),
      info = family_name,
      tolerance = 1e-7
    )
  }
})


test_that("df = 1 likelihood GLM matrix fits match base glm completely", {
  set.seed(4105)
  n <- 360L
  dat <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1),
    case_weight = sample(c(1, 2), n, replace = TRUE)
  )
  eta <- 0.30 + 0.25 * dat$x1 - 0.18 * dat$x2
  mu <- exp(eta)
  responses <- list(
    gaussian = eta + stats::rnorm(n, sd = 0.45),
    binomial = stats::rbinom(n, size = 1L, prob = stats::plogis(eta - 0.4)),
    poisson = stats::rpois(n, lambda = mu),
    Gamma = stats::rgamma(n, shape = 8, scale = mu / 8),
    inverse.gaussian = mu * exp(stats::rnorm(n, sd = 0.12))
  )
  families <- list(
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    Gamma = stats::Gamma(link = "log"),
    inverse.gaussian = stats::inverse.gaussian(link = "log")
  )
  x <- as.matrix(dat[, c("x1", "x2"), drop = FALSE])
  rows <- c(2L, 13L, 27L, 54L, 88L, 133L, 191L, 246L, 301L, 349L)
  newx <- x[rows, , drop = FALSE]
  reference_newdata <- dat[rows, c("x1", "x2"), drop = FALSE]

  for (family_name in names(families)) {
    family <- families[[family_name]]
    dat$y <- responses[[family_name]]

    fit_mfp2 <- mfp2(
      x = x,
      y = dat$y,
      family = family,
      weights = dat$case_weight,
      cycles = 1,
      df = 1,
      select = 1,
      alpha = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      xorder = "original",
      verbose = FALSE
    )
    fit_glm <- stats::glm(
      y ~ x1 + x2,
      data = dat,
      family = family,
      weights = case_weight,
      x = TRUE,
      y = TRUE,
      model = TRUE
    )

    expect_equal(
      sort(get_selected_variable_names(fit_mfp2)),
      sort(c("x1", "x2")),
      info = paste(family_name, "matrix")
    )
    expect_complete_glm_linear_equivalence(
      fit_mfp2,
      fit_glm,
      newdata = reference_newdata,
      mfp2_newdata = newx,
      term_newdata = as.data.frame(newx, check.names = FALSE),
      info = paste(family_name, "matrix"),
      tolerance = 1e-7
    )
  }
})
