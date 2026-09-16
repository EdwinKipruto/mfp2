# Regression tests for work that must be performed once outside repeated MFP
# candidate fits.

test_that("family-specific controls normalize once per public mfp2 fit", {
  set.seed(5201)
  n <- 72L
  x <- cbind(x1 = stats::rnorm(n))

  original_survreg_normalizer <- mfp2:::normalize_survreg_control
  survreg_normalizations <- 0L
  testthat::local_mocked_bindings(
    normalize_survreg_control = function(control = NULL) {
      survreg_normalizations <<- survreg_normalizations + 1L
      original_survreg_normalizer(control)
    },
    .package = "mfp2"
  )

  y <- survival::Surv(
    exp(1.5 + 0.2 * x[, "x1"] + stats::rnorm(n, sd = 0.35)),
    rep(c(1, 1, 0), length.out = n)
  )
  fit <- mfp2(
    x, y,
    family = survreg_family(),
    control = list(maxiter = 19L),
    cycles = 1,
    df = 1,
    select = 1,
    alpha = 1,
    xorder = "original",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_identical(survreg_normalizations, 1L)
  expect_equal(fit$mfp2_control$iter.max, 19L)
})


test_that("normalize_fit_control resolves partial lists for every fitter family", {
  glm_control <- normalize_fit_control(
    list(epsilon = 1e-7), "poisson", fitter = "base"
  )
  expect_equal(glm_control$epsilon, 1e-7)
  expect_true(all(c("epsilon", "maxit", "trace") %in% names(glm_control)))

  survreg_control <- normalize_fit_control(
    list(maxiter = 11L), "survreg", fitter = "base"
  )
  expect_equal(survreg_control$iter.max, 11L)

  cox_control <- normalize_fit_control(
    list(iter.max = 13L), "finegray", fitter = "base"
  )
  expect_equal(cox_control$iter.max, 13L)

  expect_error(
    normalize_fit_control("invalid", "finegray", fitter = "base"),
    "must be `NULL` or a list"
  )
})


test_that("Fine--Gray preparation caches the invariant expanded offset", {
  n <- 18L
  y <- survival::Surv(
    seq_len(n),
    factor(
      rep(c("censor", "relapse", "death"), length.out = n),
      levels = c("censor", "relapse", "death")
    )
  )
  offset <- seq(-0.3, 0.3, length.out = n)

  family <- prepare_finegray_family(
    finegray_family(etype = "relapse"),
    y = y,
    weights = rep.int(1, n),
    offset = offset
  )

  expect_equal(
    family$prepared$offset_expanded,
    offset[family$prepared$row_map]
  )
})


test_that("offset-aware negative-binomial null fitting stays out of candidates", {
  skip_if_not_installed("fastglm")

  set.seed(5202)
  n <- 96L
  x <- cbind(
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1)
  )
  offset <- stats::runif(n, -0.3, 0.3)
  weights <- sample(c(1, 2), n, replace = TRUE)
  mu <- exp(0.2 + 0.35 * x[, "x1"] - 0.15 * x[, "x2"] + offset)
  y <- stats::rnbinom(n, mu = mu, size = 2.8)

  original_null_deviance <- mfp2:::mfp2_negbin_offset_null_deviance
  null_fit_calls <- 0L
  testthat::local_mocked_bindings(
    mfp2_negbin_offset_null_deviance = function(...) {
      null_fit_calls <<- null_fit_calls + 1L
      original_null_deviance(...)
    },
    .package = "mfp2"
  )

  fit_model(
    x = x,
    y = y,
    family = "negbin",
    family_string = "negbin",
    fitter = "fastglm",
    weights = weights,
    offset = offset,
    control = stats::glm.control(),
    fast = TRUE,
    calculate_fit_statistics = FALSE
  )
  expect_identical(null_fit_calls, 0L)

  fit_model(
    x = x,
    y = y,
    family = "negbin",
    family_string = "negbin",
    fitter = "fastglm",
    weights = weights,
    offset = offset,
    control = stats::glm.control(),
    fast = TRUE,
    calculate_fit_statistics = TRUE
  )
  expect_identical(null_fit_calls, 1L)
})
