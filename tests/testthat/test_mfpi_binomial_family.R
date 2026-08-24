# End-to-end coverage for the documented MFPI binomial family.

test_that("binomial MFPI formula predictions match the stored interaction glm", {
  set.seed(2701)
  n <- 260L
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2L)),
    x = stats::runif(n, 0.5, 5),
    # The fit deliberately uses shift = 0, so keep the continuous adjustment
    # variable on a strictly positive scale.
    z = stats::runif(n, 1, 3)
  )
  eta <- with(
    dat,
    -0.8 + 0.2 * x + 0.25 * (group == "treated") +
      0.15 * x * (group == "treated") - 0.15 * z
  )
  dat$y <- stats::rbinom(n, size = 1L, prob = stats::pnorm(eta))

  fit <- mfpi(
    y ~ group + fp(x, df = 1) + z,
    data = dat,
    family = stats::binomial(link = "probit"),
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    p_interact = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_identical(fit$family_string, "binomial")

  nd <- dat[1:24, c("group", "x", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  design <- mfpi_build_ordinary_design(fit, "x", fit_result, nd)

  direct_link <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "link",
    se.fit = TRUE
  )
  got_link <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )

  direct_response <- stats::predict(
    stored,
    newdata = design$model_newdata,
    type = "response",
    se.fit = FALSE
  )
  got_response <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "response",
    se.fit = FALSE
  )

  expect_identical(stored$family$link, "probit")
  expect_true(got_link$metadata$used_model_predict)
  expect_equal(
    got_link$predictions$fit,
    as.numeric(direct_link$fit),
    tolerance = 1e-8
  )
  expect_equal(
    got_link$predictions$se.fit,
    as.numeric(direct_link$se.fit),
    tolerance = 1e-8
  )
  expect_equal(
    got_response$predictions$fit,
    as.numeric(direct_response),
    tolerance = 1e-8
  )
  expect_equal(
    got_response$predictions$fit,
    stats::pnorm(got_link$predictions$fit),
    tolerance = 1e-8
  )
  expect_true(all(got_response$predictions$fit >= 0))
  expect_true(all(got_response$predictions$fit <= 1))
})
