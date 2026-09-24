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


test_that("fit_model requires an already-normalized control object", {
  expect_error(
    fit_model(
      x = matrix(1:6, ncol = 1L),
      y = 1:6,
      family = stats::gaussian(),
      family_string = "gaussian"
    ),
    "argument.*control.*missing"
  )
})


test_that("only ordinal hot-path matrices omit a model intercept", {
  expect_true(mfp2_candidate_matrix_has_intercept("gaussian"))
  expect_true(mfp2_candidate_matrix_has_intercept("multinomial"))
  expect_true(mfp2_candidate_matrix_has_intercept("survreg"))
  expect_false(mfp2_candidate_matrix_has_intercept("ordinal"))
  expect_false(mfp2_candidate_matrix_has_intercept("cox"))
  expect_false(mfp2_candidate_matrix_has_intercept("finegray"))

  # Ordinal models still have k - 1 threshold intercepts; only their explicit
  # candidate-matrix intercept is omitted.
  expect_true(mfp2_family_has_intercept("ordinal"))
})


test_that("GLM preparation caches invariant fitter flags", {
  prepared <- prepare_family_for_fit(
    family = stats::gaussian(),
    family_string = "gaussian",
    y = 1:4,
    weights = rep.int(1, 4L)
  )$family

  expect_identical(
    attr(prepared, "mfp2_fit_flags", exact = TRUE),
    list(
      is_gaussian = TRUE,
      is_negbin = FALSE,
      estimates_dispersion = TRUE
    )
  )
  expect_null(attr(mfp2_strip_prepared_family(prepared), "mfp2_fit_flags",
                   exact = TRUE))
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


test_that("survival preparation caches low-level strata codes", {
  n <- 18L
  strata <- factor(rep(c("site1", "site2"), each = n / 2L))
  cox_prepared <- prepare_family_for_fit(
    family = "cox",
    family_string = "cox",
    y = survival::Surv(seq_len(n), rep(c(1, 0), length.out = n)),
    weights = rep.int(1, n),
    strata = strata
  )
  expect_identical(
    attr(cox_prepared$strata, "mfp2_integer_codes", exact = TRUE),
    as.integer(strata)
  )

  y_fg <- survival::Surv(
    seq_len(n),
    factor(
      rep(c("censor", "relapse", "death"), length.out = n),
      levels = c("censor", "relapse", "death")
    )
  )
  fg <- prepare_finegray_family(
    finegray_family(etype = "relapse"),
    y = y_fg,
    weights = rep.int(1, n),
    strata = strata
  )$prepared
  expect_identical(fg$strata_expanded_codes, as.integer(fg$strata_expanded))
})


test_that("multinomial preparation caches candidate response, weights, and offset", {
  y <- factor(c("A", "B", "C", "A"), levels = c("A", "B", "C"))
  weights <- c(1, 2, 3, 4)
  offset <- cbind(
    A = c(0.2, 0.1, -0.1, 0.3),
    B = c(0.4, 0.0, 0.2, -0.2),
    C = c(-0.1, 0.3, 0.1, 0.0)
  )

  family <- prepare_multinomial_family(
    multinomial_family(), y, weights,
    offset = offset, has_offset = TRUE
  )
  prepared <- family$prepared

  expect_equal(prepared$y_fit, prepared$y_matrix)
  expect_equal(prepared$effective_weights, weights)
  expect_equal(prepared$case_weights, weights)
  expect_true(prepared$has_offset)
  expect_equal(
    prepared$offset_matrix,
    sweep(offset, 1L, offset[, "A"], "-")
  )

  no_offset <- prepare_multinomial_family(
    multinomial_family(), y, weights,
    offset = rep.int(0, length(y)), has_offset = FALSE
  )
  expect_null(no_offset$prepared$offset_matrix)
  expect_false(no_offset$prepared$has_offset)

  counts <- rbind(
    c(A = 2, B = 1, C = 0),
    c(A = 0, B = 3, C = 2),
    c(A = 1, B = 1, C = 2)
  )
  count_weights <- c(1, 2, 0.5)
  grouped <- prepare_multinomial_family(
    multinomial_family(), counts, count_weights, has_offset = FALSE
  )$prepared
  totals <- rowSums(counts)
  expect_equal(grouped$y_fit, counts / totals)
  expect_equal(grouped$effective_weights, count_weights * totals)
})


test_that("multinomial fast fits can omit discarded coefficient output", {
  skip_if_not_installed("nnet")
  set.seed(5203)
  n <- 120L
  y <- factor(rep(c("A", "B", "C"), length.out = n),
              levels = c("A", "B", "C"))
  family <- prepare_multinomial_family(
    multinomial_family(), y, rep.int(1, n), has_offset = FALSE
  )
  x <- cbind(`(Intercept)` = 1, x1 = stats::rnorm(n))

  lightweight <- fit_multinomial(
    x = x,
    family = family,
    control = normalize_multinomial_control(),
    fast = TRUE,
    keep_fit = FALSE,
    keep_coefficients = FALSE,
    x_has_intercept = TRUE
  )
  complete <- fit_multinomial(
    x = x,
    family = family,
    control = normalize_multinomial_control(),
    fast = TRUE,
    keep_fit = FALSE,
    keep_coefficients = TRUE,
    x_has_intercept = TRUE
  )

  expect_equal(lightweight$logl, complete$logl, tolerance = 1e-10)
  expect_identical(lightweight$df, complete$df)
  expect_null(lightweight$coefficients)
  expect_null(lightweight$coefficient_matrix)
  expect_true(is.numeric(complete$coefficients))
  expect_true(is.matrix(complete$coefficient_matrix))
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
