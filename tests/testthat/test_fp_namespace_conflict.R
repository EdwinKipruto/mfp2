# Regression coverage for the fp() name shared by mfp and mfp2.

make_shadowed_fp_formula_v <- function(formula) {
  shadow_env <- new.env(parent = environment(formula))
  shadow_env$fp <- function(...) {
    stop("the shadowing fp() was called", call. = FALSE)
  }
  environment(formula) <- shadow_env
  formula
}


test_that("mfp2 prediction retains mfp2 fp semantics under shadowing", {
  set.seed(2601)
  n <- 100L
  dat <- data.frame(
    x = seq(1, 10, length.out = n),
    z = stats::runif(n, 1, 3)
  )
  dat$y <- 0.4 + 0.7 * dat$x - 0.2 * dat$z + stats::rnorm(n, sd = 0.1)

  formula <- make_shadowed_fp_formula_v(
    y ~ fp(x, df = 1, force_max_fp = TRUE) + z
  )
  fit <- mfp2(
    formula,
    data = dat,
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    verbose = FALSE,
    warn_low_information = FALSE
  )

  nd <- dat[1:12, c("x", "z"), drop = FALSE]
  prediction <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)

  expect_s3_class(fit, "mfp2")
  expect_length(prediction$fit, nrow(nd))
  expect_length(prediction$se.fit, nrow(nd))
  expect_true(all(is.finite(prediction$fit)))
  expect_true(all(is.finite(prediction$se.fit)))
})


test_that("mfpi fitting and prediction ignore a shadowing fp function", {
  set.seed(2602)
  n <- 160L
  dat <- data.frame(
    group = factor(rep(c("control", "treated"), each = n / 2L)),
    x = stats::runif(n, 1, 7),
    z = stats::runif(n, 1, 3)
  )
  dat$y <- with(
    dat,
    0.5 + 0.35 * x + 0.8 * (group == "treated") * x +
      0.2 * z + stats::rnorm(n, sd = 0.2)
  )

  formula <- make_shadowed_fp_formula_v(
    y ~ group + fp(x, df = 1, force_max_fp = TRUE) + z
  )
  fit <- mfpi(
    formula,
    data = dat,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    p_interact = 1,
    verbose = FALSE
  )

  nd <- dat[1:14, c("group", "x", "z"), drop = FALSE]
  prediction <- predict(
    fit,
    newdata = nd,
    terms = "x",
    model = "all",
    type = "link",
    se.fit = TRUE
  )

  expect_s3_class(fit, "mfpi")
  expect_length(prediction$predictions$fit, nrow(nd))
  expect_length(prediction$predictions$se.fit, nrow(nd))
  expect_true(all(is.finite(prediction$predictions$fit)))
  expect_true(all(is.finite(prediction$predictions$se.fit)))
})
