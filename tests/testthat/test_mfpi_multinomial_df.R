# MFPI common-power multinomial degree-of-freedom rules ----------------------

test_that("MFPI df formulas multiply coefficients but not shared powers", {
  expect_equal(
    interaction_model_df(n_groups = 3, degree = 0, n_logits = 2,
                         flex = "flex1"),
    list(dfmain = 6L, dfint = 4L, total_df = 10L, n_groups = 3L)
  )
  expect_equal(
    interaction_model_df(n_groups = 3, degree = 1, n_logits = 2,
                         flex = "flex1"),
    list(dfmain = 7L, dfint = 4L, total_df = 11L, n_groups = 3L)
  )
  expect_equal(
    interaction_model_df(n_groups = 3, degree = 2, n_logits = 2,
                         flex = "flex3"),
    list(dfmain = 10L, dfint = 8L, total_df = 18L, n_groups = 3L)
  )
  expect_equal(
    interaction_model_df(n_groups = 3, degree = 2, n_logits = 2,
                         flex = "flex4"),
    list(dfmain = 10L, dfint = 12L, total_df = 22L, n_groups = 3L)
  )
})


test_that("forced-linear multinomial MFPI uses Q-scaled interaction df", {
  set.seed(1615L)
  n <- 240L
  dat <- data.frame(
    group = factor(rep(c("control", "active"), each = n / 2L)),
    x = stats::runif(n, 0.5, 4),
    z = stats::rnorm(n)
  )
  eta_b <- with(dat, -0.3 + 0.25 * x + 0.20 * (group == "active") - 0.2 * z)
  eta_c <- with(dat, 0.2 - 0.15 * x - 0.10 * (group == "active") + 0.3 * z)
  denom <- 1 + exp(eta_b) + exp(eta_c)
  probs <- cbind(A = 1 / denom, B = exp(eta_b) / denom,
                 C = exp(eta_c) / denom)
  draw <- apply(probs, 1L, function(p) sample.int(3L, 1L, prob = p))
  dat$y <- factor(colnames(probs)[draw], levels = colnames(probs))

  fit <- mfpi(
    y ~ group + x + z,
    data = dat,
    family = "multinomial",
    group_var = "group",
    interaction_vars = "x",
    interaction_forms = c(x = "linear"),
    flex = "flex1",
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    p_interact = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_identical(fit$family_string, "multinomial")
  expect_identical(fit$n_logits, 2L)
  expect_equal(fit$all_model_metrics$df_int, 2)
})
