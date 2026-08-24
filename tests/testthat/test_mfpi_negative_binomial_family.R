# End-to-end coverage for the documented MFPI negative-binomial family.

make_mfpi_negbin_data_v <- function(n = 280L, seed = 2801L) {
  set.seed(seed)
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2L)),
    x = stats::runif(n, 0.5, 4.5),
    z = stats::rnorm(n)
  )
  eta <- with(
    dat,
    0.2 + 0.25 * x + 0.3 * (group == "treated") +
      0.22 * x * (group == "treated") - 0.15 * z
  )
  dat$y <- stats::rnbinom(n, mu = exp(eta), size = 3.5)
  dat
}


test_that("negative-binomial MFPI requires the fastglm fitter", {
  dat <- make_mfpi_negbin_data_v(n = 80L, seed = 2802L)

  expect_error(
    mfpi(
      dat[, c("group", "x", "z")],
      dat$y,
      family = "negbin",
      fitter = "base",
      group_var = "group",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      flex = "flex1",
      verbose = FALSE
    ),
    "available only with.*fitter = \"fastglm\""
  )
})


test_that("negative-binomial MFPI reconstructs link and response predictions", {
  skip_if_not_installed("fastglm")

  dat <- make_mfpi_negbin_data_v()
  fit <- mfpi(
    dat[, c("group", "x", "z")],
    dat$y,
    family = "negbin",
    fitter = "fastglm",
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
  expect_identical(fit$family_string, "negbin")
  expect_identical(fit$fitter, "fastglm")

  nd <- dat[1:22, c("group", "x", "z"), drop = FALSE]
  fit_result <- fit$var_winners[["x"]]$fit
  stored <- fit_result$test_results$interaction_model$fit
  design <- mfpi_build_ordinary_design(fit, "x", fit_result, nd)

  link <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )
  response <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "response",
    se.fit = FALSE
  )

  beta <- stored$coefficients
  manual_x <- design$X
  if ("(Intercept)" %in% names(beta) &&
      !"(Intercept)" %in% colnames(manual_x)) {
    manual_x <- cbind("(Intercept)" = 1, manual_x)
  }
  manual_x <- manual_x[, names(beta), drop = FALSE]
  estimable <- !is.na(beta)
  manual_eta <- as.numeric(
    manual_x[, estimable, drop = FALSE] %*% beta[estimable]
  )

  expect_s3_class(stored, "fastglm")
  expect_true(is.numeric(stored$theta))
  expect_true(is.finite(stored$theta) && stored$theta > 0)
  expect_false(link$metadata$used_model_predict)
  expect_equal(link$predictions$fit, manual_eta, tolerance = 1e-8)
  expect_equal(
    response$predictions$fit,
    exp(link$predictions$fit),
    tolerance = 1e-8
  )
  expect_length(link$predictions$se.fit, nrow(nd))
  expect_true(all(is.finite(link$predictions$se.fit)))
  expect_true(all(link$predictions$se.fit >= 0))
  expect_true(all(is.finite(response$predictions$fit)))
  expect_true(all(response$predictions$fit > 0))
})
