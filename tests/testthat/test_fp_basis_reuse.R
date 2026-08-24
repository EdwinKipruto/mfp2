# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# Shared focal max-degree basis regression tests
# =============================================================================

# Test purpose: A maximum-degree ordinary FP basis must support lower-degree
# views for a non-default, variable-specific power set without changing the
# existing candidate order or the special FP1 exclusion of power 1.
test_that("shared FP basis reproduces custom-power degree views", {
  x <- seq(1, 12, length.out = 60)
  allowed_powers <- c(-1, 0, 0.5, 1, 2)

  shared <- build_shared_focal_fp_basis(
    x = x,
    max_degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  # find_best_fpm_step() removes power 1 only for FP1. Reproduce that existing
  # rule here and verify that the shared max-degree basis still maps every
  # remaining FP1 candidate correctly.
  fp1_powers <- setdiff(allowed_powers, 1)
  shared_fp1 <- view_shared_focal_fp_basis(
    shared_basis = shared,
    degree = 1,
    powers = fp1_powers
  )
  direct_fp1 <- generate_transformations_fp_basis(
    x = x,
    degree = 1,
    powers = fp1_powers,
    zero = FALSE
  )

  expect_equal(shared_fp1$powers, direct_fp1$powers)
  expect_false(any(shared_fp1$powers == 1))

  for (i in seq_len(nrow(shared_fp1$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(shared_fp1, i)),
      unname(materialize_fp_basis_candidate(direct_fp1, i)),
      tolerance = 1e-12
    )
  }

  # FP2 keeps the user's complete power set, including power 1. Its candidate
  # order and transformed values must also match a separately generated basis.
  shared_fp2 <- view_shared_focal_fp_basis(
    shared_basis = shared,
    degree = 2,
    powers = allowed_powers
  )
  direct_fp2 <- generate_transformations_fp_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  expect_equal(shared_fp2$powers, direct_fp2$powers)
  expect_true(any(shared_fp2$powers == 1))

  for (i in seq_len(nrow(shared_fp2$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(shared_fp2, i)),
      unname(materialize_fp_basis_candidate(direct_fp2, i)),
      tolerance = 1e-12
    )
  }
})


# Test purpose: One joint ACD degree-2 basis must reproduce all three nonlinear
# IC views: FP1(x,.), FP1(.,A(x)), and FP1(x,A(x)), again with custom powers.
test_that("shared ACD basis reproduces all custom-power IC views", {
  x <- seq(1, 10, length.out = 50)
  allowed_powers <- c(-1, 0, 1, 2)
  fp1_powers <- setdiff(allowed_powers, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  shared <- build_shared_focal_acd_basis(
    x = x,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  # FP1(x,.) uses the x component of the joint ACD basis through the ordinary
  # FP view. No ACD-specific candidate semantics are introduced here.
  shared_x <- view_shared_focal_fp_basis(
    shared_basis = shared,
    degree = 1,
    powers = fp1_powers
  )
  direct_x <- generate_transformations_fp_basis(
    x = x,
    degree = 1,
    powers = fp1_powers,
    zero = FALSE
  )

  expect_equal(shared_x$powers, direct_x$powers)
  for (i in seq_len(nrow(shared_x$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(shared_x, i)),
      unname(materialize_fp_basis_candidate(direct_x, i)),
      tolerance = 1e-12
    )
  }

  # FP1(.,A(x)) uses only the A(x) component.
  shared_a <- view_shared_focal_acd_basis(
    shared_basis = shared,
    degree = 1,
    powers = fp1_powers
  )
  direct_a <- generate_transformations_acd_basis(
    x = x,
    degree = 1,
    powers = fp1_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_equal(shared_a$powers, direct_a$powers)
  for (i in seq_len(nrow(shared_a$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(shared_a, i)),
      unname(materialize_acd_basis_candidate(direct_a, i)),
      tolerance = 1e-12
    )
  }

  # FP1(x,A(x)) uses one x and one A(x) column from the same shared basis.
  shared_joint <- view_shared_focal_acd_basis(
    shared_basis = shared,
    degree = 2,
    powers = allowed_powers
  )
  direct_joint <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_equal(shared_joint$powers, direct_joint$powers)
  for (i in seq_len(nrow(shared_joint$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(shared_joint, i)),
      unname(materialize_acd_basis_candidate(direct_joint, i)),
      tolerance = 1e-12
    )
  }
})


# Test purpose: transform_data_step() must use a supplied ordinary shared basis
# rather than regenerating n-length FP transformations for the requested degree.
test_that("transform_data_step reuses supplied shared FP basis", {
  x_vec <- seq(1, 8, length.out = 40)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 0.5, 1, 2)

  shared <- build_shared_focal_fp_basis(
    x = x_vec,
    max_degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  testthat::local_mocked_bindings(
    generate_transformations_fp_basis = function(...) {
      stop("ordinary FP basis was regenerated")
    },
    .package = "mfp2"
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 2,
    powers = list(x = setdiff(allowed_powers, 1)),
    acdx = c(x = FALSE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    focal_basis_cache = shared,
    term_to_columns = list(x = 1L),
    compact_fp = TRUE
  )

  expect_null(out$data_fp)
  expect_equal(out$fp_basis$powers[, 1L], sort(setdiff(allowed_powers, 1)))
})


# Test purpose: transform_data_step() must likewise reuse a supplied joint ACD
# basis without invoking the degree-specific ACD basis generator again.
test_that("transform_data_step reuses supplied shared ACD basis", {
  x_vec <- seq(1, 8, length.out = 40)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 1, 2)
  fp1_powers <- setdiff(allowed_powers, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  shared <- build_shared_focal_acd_basis(
    x = x_vec,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  testthat::local_mocked_bindings(
    generate_transformations_acd_basis = function(...) {
      stop("ACD basis was regenerated")
    },
    .package = "mfp2"
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 2,
    powers = list(x = fp1_powers),
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
    focal_basis_cache = shared,
    term_to_columns = list(x = 1L),
    compact_acd = TRUE
  )

  expect_null(out$data_fp)
  expect_equal(out$acd_basis$powers[, 2L], fp1_powers)
})


# Test purpose: select_ic() should build one maximum-degree basis from the
# complete custom power set and pass the same cache to every FP degree.
test_that("select_ic shares one focal basis across FP degrees", {
  x <- matrix(seq(1, 10, length.out = 40), ncol = 1L,
              dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 0.5, 1, 2)
  sentinel_cache <- list(id = "shared-fp-cache")
  builder_calls <- 0L
  seen_degrees <- integer(0)
  seen_cache <- list()

  fake_metrics <- function(aic) {
    matrix(
      c(0, 1, 0, NA, aic, aic, 38),
      nrow = 1L,
      dimnames = list(NULL, c(
        "logl", "df", "deviance_rs", "deviance_gaussian",
        "aic", "bic", "df_resid"
      ))
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        powers_adj = list(),
        spike_decision_adj = numeric(0),
        data_adj_list = list(),
        data_adj = NULL,
        transform_cache = list()
      )
    },
    fit_null_step = function(...) {
      list(powers = NA, metrics = fake_metrics(10), current_adj_params = list())
    },
    fit_linear_step = function(...) {
      list(powers = 1, metrics = fake_metrics(8), current_adj_params = list())
    },
    build_shared_focal_fp_basis = function(x, max_degree, powers, ...) {
      builder_calls <<- builder_calls + 1L
      expect_equal(max_degree, 2)
      expect_equal(powers, allowed_powers)
      sentinel_cache
    },
    find_best_fpm_step = function(..., degree, focal_basis_cache = NULL) {
      seen_degrees <<- c(seen_degrees, degree)
      seen_cache[[length(seen_cache) + 1L]] <<- focal_basis_cache
      list(
        powers = if (degree == 1) matrix(-1, nrow = 1L) else
          matrix(c(-1, 0), nrow = 1L),
        metrics = fake_metrics(if (degree == 1) 6 else 4),
        model_best = 1L,
        current_adj_params = list()
      )
    },
    .package = "mfp2"
  )

  select_ic(
    x = x,
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = FALSE),
    y = rep(0, nrow(x)),
    powers_current = list(x = c(1, 1)),
    powers = list(x = allowed_powers),
    criterion = "aic",
    ftest = FALSE,
    select = 1,
    alpha = 0.05,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = NULL),
    prev_adj_params = list(),
    transform_cache = list(),
    force_max_fp = c(x = FALSE),
    has_offset = FALSE,
    n_obs = nrow(x),
    term_to_columns = list(x = 1L)
  )

  expect_equal(builder_calls, 1L)
  expect_equal(seen_degrees, c(1, 2))
  expect_length(seen_cache, 2L)
  expect_true(all(vapply(seen_cache, identical, logical(1), sentinel_cache)))
})


# Test purpose: select_ic_acd() should build one joint x/A(x) basis and reuse
# it for FP1(x,.), FP1(.,A(x)), and FP1(x,A(x)).
test_that("select_ic_acd shares one joint focal basis across all nonlinear views", {
  x <- matrix(seq(1, 10, length.out = 40), ncol = 1L,
              dimnames = list(NULL, "x"))
  allowed_powers <- c(-1, 0, 1, 2)
  sentinel_cache <- list(id = "shared-acd-cache")
  builder_calls <- 0L
  seen_degrees <- integer(0)
  seen_cache <- list()

  fake_metrics <- function(aic) {
    matrix(
      c(0, 1, 0, NA, aic, aic, 38),
      nrow = 1L,
      dimnames = list(NULL, c(
        "logl", "df", "deviance_rs", "deviance_gaussian",
        "aic", "bic", "df_resid"
      ))
    )
  }

  testthat::local_mocked_bindings(
    build_adjustment_step = function(...) {
      list(
        powers_adj = list(),
        spike_decision_adj = numeric(0),
        data_adj_list = list(),
        data_adj = NULL,
        transform_cache = list()
      )
    },
    fit_null_step = function(...) {
      list(powers = c(NA, NA), metrics = fake_metrics(12), current_adj_params = list())
    },
    fit_linear_step = function(...) {
      list(powers = c(NA, 1), metrics = fake_metrics(10), current_adj_params = list())
    },
    build_shared_focal_acd_basis = function(x, powers, ...) {
      builder_calls <<- builder_calls + 1L
      expect_equal(powers, allowed_powers)
      sentinel_cache
    },
    find_best_fpm_step = function(..., degree, acdx, focal_basis_cache = NULL) {
      seen_degrees <<- c(seen_degrees, degree)
      seen_cache[[length(seen_cache) + 1L]] <<- focal_basis_cache
      list(
        powers = if (degree == 1) matrix(c(NA, -1), nrow = 1L) else
          matrix(c(-1, 0), nrow = 1L),
        metrics = fake_metrics(if (degree == 1) 8 else 6),
        model_best = 1L,
        current_adj_params = list()
      )
    },
    .package = "mfp2"
  )

  select_ic_acd(
    x = x,
    xi = "x",
    keep = character(0),
    degree = 2,
    acdx = c(x = TRUE),
    y = rep(0, nrow(x)),
    powers_current = list(x = c(1, 1)),
    powers = list(x = allowed_powers),
    criterion = "aic",
    ftest = FALSE,
    select = 1,
    alpha = 0.05,
    family = stats::gaussian(),
    family_string = "gaussian",
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = list(acd = seq(0.1, 0.9, length.out = nrow(x)))),
    prev_adj_params = list(),
    transform_cache = list(),
    force_max_fp = c(x = FALSE),
    has_offset = FALSE,
    n_obs = nrow(x),
    term_to_columns = list(x = 1L)
  )

  expect_equal(builder_calls, 1L)
  expect_equal(seen_degrees, c(1, 1, 2))
  expect_length(seen_cache, 3L)
  expect_true(all(vapply(seen_cache, identical, logical(1), sentinel_cache)))
})
