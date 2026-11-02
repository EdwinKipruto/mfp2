# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# ACD training-cache and compact-basis invariant regression tests
# =============================================================================

# Test purpose: The explicit fit-time A(x) cache must reproduce the historical
# apply_acd() path exactly for compact ACD candidate generation on the same
# training observations.
test_that("compact ACD generation reuses explicit training A(x) cache", {
  set.seed(411)
  x <- runif(100, 1, 20)
  powers <- c(0, 1)

  acd_fit <- fit_acd(x, powers = powers)
  acd_par <- acd_fit
  acd_par$acd <- NULL

  applied <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  cached <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par,
    acd_training_values = acd_fit$acd
  )

  expect_equal(cached$powers, applied$powers)
  expect_equal(cached$candidate_map, applied$candidate_map)
  expect_equal(cached$basis, applied$basis, tolerance = 1e-12)
})


# Test purpose: The materialized ACD generator shares the same explicit cache
# contract as the compact generator, so legacy/non-compact callers remain
# numerically identical when the training A(x) values are supplied.
test_that("materialized ACD generation reuses explicit training A(x) cache", {
  set.seed(412)
  x <- runif(80, 1, 15)
  powers <- c(0, 1)

  acd_fit <- fit_acd(x, powers = powers)
  acd_par <- acd_fit
  acd_par$acd <- NULL

  applied <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  cached <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par,
    acd_training_values = acd_fit$acd
  )

  expect_equal(cached$powers, applied$powers)
  expect_equal(length(cached$data), length(applied$data))
  for (i in seq_along(applied$data)) {
    expect_equal(cached$data[[i]], applied$data[[i]], tolerance = 1e-12)
  }
})


# Test purpose: A stored $acd component is training-data state, not a general
# prediction cache. Unless it is passed explicitly as acd_training_values, the
# generator must continue applying the stored ACD parameters to the supplied x.
test_that("ACD generator does not implicitly reuse acd_parameter$acd", {
  set.seed(413)
  x_fit <- runif(80, 1, 10)
  x_new <- runif(80, 10, 20)
  powers <- c(0, 1)

  acd_fit <- fit_acd(x_fit, powers = powers)

  with_stored_acd <- generate_transformations_acd_basis(
    x = x_new,
    degree = 1,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_fit
  )

  acd_par_without_values <- acd_fit
  acd_par_without_values$acd <- NULL

  without_stored_acd <- generate_transformations_acd_basis(
    x = x_new,
    degree = 1,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par_without_values
  )

  expect_equal(
    with_stored_acd$basis,
    without_stored_acd$basis,
    tolerance = 1e-12
  )
})


# Test purpose: The explicit training cache is row-aligned with x. Reject a
# wrong-length vector immediately instead of allowing recycling or a later,
# less informative matrix-shape failure.
test_that("ACD training cache validates observation count", {
  x <- seq(1, 10, length.out = 30)

  expect_error(
    generate_transformations_acd_basis(
      x = x,
      degree = 1,
      powers = c(0, 1),
      zero = FALSE,
      acd_parameter = list(
        beta0 = 0,
        beta1 = 1,
        power = 1,
        shift = 0,
        scale = 1
      ),
      acd_training_values = rep(0.5, length(x) - 1L)
    ),
    "one value per observation"
  )
})


# Test purpose: The normal model-search wrapper must pass fit_acd()$acd to the
# compact ACD generator. Invalid application coefficients make this test fail
# if apply_acd() is accidentally called instead of using the training cache.
test_that("transform_data_step() forwards cached training A(x)", {
  x_vec <- seq(1, 5, length.out = 25)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))

  acd_par <- list(
    acd = seq(0.2, 0.8, length.out = length(x_vec)),
    beta0 = NA_real_,
    beta1 = NA_real_,
    power = 1,
    shift = 0,
    scale = 1
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 4,
    powers = list(x = c(0, 1)),
    acdx = c(x = TRUE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = acd_par),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    term_to_columns = list(x = 1L),
    compact_acd = TRUE
  )

  expect_true(all(is.finite(out$acd_basis$basis)))
})


# Test purpose: Ordinary FP and ACD compact representations are mutually
# exclusive for one focal variable. If a future control-flow change violates
# that invariant, find_best_fpm_step() must fail explicitly instead of silently
# preferring fp_basis and potentially fitting the wrong candidate layout.
test_that("find_best_fpm_step() rejects simultaneous FP and ACD compact bases", {
  testthat::local_mocked_bindings(
    transform_data_step = function(...) {
      list(
        data_adj = NULL,
        data_fp = NULL,
        fp_basis = list(source = "fp"),
        acd_basis = list(source = "acd")
      )
    },
    .package = "mfp2"
  )

  x <- matrix(seq_len(8), ncol = 1L, dimnames = list(NULL, "x"))

  expect_error(
    find_best_fpm_step(
      x = x,
      xi = "x",
      degree = 2,
      y = rep(0, nrow(x)),
      powers_current = list(x = c(1, 1)),
      powers = list(x = c(0, 1)),
      acdx = c(x = FALSE),
      family = stats::gaussian(),
      family_string = "gaussian",
      zero = c(x = FALSE),
      catzero = list(x = NULL),
      spike = c(x = FALSE),
      spike_decision = c(x = 2),
      acd_parameter = list(x = NULL),
      prev_adj_params = list(),
      has_offset = FALSE,
      precomputed_adj = NULL,
      n_obs = nrow(x),
      term_to_columns = list(x = 1L)
    ),
    "both `fp_basis` and `acd_basis` are non-NULL"
  )
})
