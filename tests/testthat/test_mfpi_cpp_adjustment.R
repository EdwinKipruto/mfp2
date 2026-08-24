# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Checks that the C++ adjustment-step bridge returns NULL when
# there are no adjustment variables.
test_that("C++ adjustment-step bridge handles no adjustment variables", {
  x <- matrix(c(1, 2, 3), ncol = 1L)
  colnames(x) <- "x1"

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = character(0),
    powers_adj = list(),
    acdx_adj = list(),
    zero_adj = list(),
    catzero = list(),
    spike_adj = list(),
    spike_decision_int_adj = integer(0),
    acd_parameter_adj = list(),
    eliminated = logical(0),
    spike_binary_only_flags = logical(0),
    current_power_keys_adj = list(),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expect_equal(out$data_adj_list, list())
  expect_null(out$data_adj)
})


# Test purpose: Checks that the C++ adjustment-step bridge builds ordinary FP
# adjustment columns on cache miss.
test_that("C++ adjustment-step bridge builds FP adjustment columns", {
  x <- cbind(
    x1 = c(1, 2, 4),
    x2 = c(2, 3, 5)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = c("x1", "x2"),
    powers_adj = list(x1 = c(1), x2 = c(0, 0)),
    acdx_adj = list(x1 = FALSE, x2 = FALSE),
    zero_adj = list(x1 = FALSE, x2 = FALSE),
    catzero = list(x1 = NULL, x2 = NULL),
    spike_adj = list(x1 = FALSE, x2 = FALSE),
    spike_decision_int_adj = c(x1 = 2L, x2 = 2L),
    acd_parameter_adj = list(x1 = NULL, x2 = NULL),
    eliminated = c(x1 = FALSE, x2 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE, x2 = FALSE),
    current_power_keys_adj = list(x1 = c(1), x2 = c(0, 0)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected_x1 <- matrix(x[, "x1"], ncol = 1L)
  colnames(expected_x1) <- "x1_adj1"

  expected_x2 <- cbind(log(x[, "x2"]), log(x[, "x2"])^2)
  colnames(expected_x2) <- c("x2_adj1", "x2_adj2")

  expected <- cbind(expected_x1, expected_x2)

  expect_equal(names(out$data_adj_list), c("x1", "x2"))
  expect_equal(out$data_adj_list$x1, expected_x1, tolerance = 1e-12)
  expect_equal(out$data_adj_list$x2, expected_x2, tolerance = 1e-12)
  expect_equal(out$data_adj, expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ adjustment-step bridge prepends catzero for
# a non-spike variable.
test_that("C++ adjustment-step bridge prepends catzero for non-spike variable", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- cbind(cz[, 1], x[, "x1"])
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj_list$x1, expected)
  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies spike
# decision 1 as catzero plus continuous FP columns.
test_that("C++ adjustment-step bridge handles spike decision 1", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 1L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- cbind(cz[, 1], x[, "x1"])
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies spike
# decision 2 as continuous FP columns only.
test_that("C++ adjustment-step bridge handles spike decision 2", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- matrix(x[, "x1"], ncol = 1L)
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies spike
# decision 3 as catzero only.
test_that("C++ adjustment-step bridge handles spike decision 3", {
  x <- cbind(x1 = c(1, 2, 4, 8))
  cz <- matrix(c(0, 1, 0, 1), ncol = 1L)

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(NA_real_)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 3L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = TRUE),
    spike_binary_only_flags = c(x1 = TRUE),
    current_power_keys_adj = list(x1 = NA_real_),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- matrix(cz[, 1], ncol = 1L)
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge returns n x 0 data_adj
# when adjustment variables exist but all are eliminated.
test_that("C++ adjustment-step bridge handles all-eliminated adjustment variables", {
  x <- cbind(
    x1 = c(1, 2, 4),
    x2 = c(2, 3, 5)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = c("x1", "x2"),
    powers_adj = list(x1 = c(NA_real_), x2 = c(NA_real_)),
    acdx_adj = list(x1 = FALSE, x2 = FALSE),
    zero_adj = list(x1 = FALSE, x2 = FALSE),
    catzero = list(x1 = NULL, x2 = NULL),
    spike_adj = list(x1 = FALSE, x2 = FALSE),
    spike_decision_int_adj = c(x1 = 2L, x2 = 2L),
    acd_parameter_adj = list(x1 = NULL, x2 = NULL),
    eliminated = c(x1 = TRUE, x2 = TRUE),
    spike_binary_only_flags = c(x1 = FALSE, x2 = FALSE),
    current_power_keys_adj = list(x1 = NA_real_, x2 = NA_real_),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expect_true(is.matrix(out$data_adj))
  expect_equal(nrow(out$data_adj), nrow(x))
  expect_equal(ncol(out$data_adj), 0L)
  expect_equal(NCOL(out$data_adj), 0L)
  expect_equal(ncol(out$data_adj_list$x1), 0L)
  expect_equal(ncol(out$data_adj_list$x2), 0L)
})


# Test purpose: The persistent focal-variable cache must drop whole-step
# matrices after a focal fit finishes, while retaining the per-variable blocks
# and metadata needed for fallback reuse in a later cycle.
test_that("persistent adjustment cache drops assembled matrices", {
  cached_block <- matrix(c(1, 2, 3), ncol = 1L)
  assembled <- cbind(cached_block, cached_block)

  params_xi <- list(
    powers_adj = list(x1 = 1),
    spike_decision_adj = c(x1 = 2L),
    data_adj_list = list(x1 = cached_block),
    data_adj = assembled,
    data_xi = matrix(c(4, 5, 6), ncol = 1L)
  )

  out <- compact_prev_adj_cache_entry(params_xi)

  # Keep exact named NULL entries. If `data_adj` were removed entirely, R's
  # partial `$` matching could resolve `out$data_adj` to `out$data_adj_list`.
  expect_true("data_adj" %in% names(out))
  expect_true("data_xi" %in% names(out))
  expect_null(out[["data_adj", exact = TRUE]])
  expect_null(out[["data_xi", exact = TRUE]])
  expect_null(out$data_adj)
  expect_null(out$data_xi)
  expect_identical(out$powers_adj, params_xi$powers_adj)
  expect_identical(out$spike_decision_adj, params_xi$spike_decision_adj)
  expect_identical(out$data_adj_list$x1, cached_block)
})


# Test purpose: build_adjustment_step() must still be able to reuse a historical
# per-focal block when prev_adj_params no longer contains a complete data_adj
# matrix. This protects the compatibility fallback after cache compaction.
test_that("adjustment fallback works without cached assembled matrix", {
  x <- cbind(
    x1 = c(1, 2, 4),
    x2 = c(3, 5, 7)
  )

  cached <- matrix(c(99, 98, 97), ncol = 1L)
  colnames(cached) <- "cached_x1"

  prev_adj_params <- list(
    x1 = NULL,
    x2 = list(
      powers_adj = list(x1 = c(1)),
      spike_decision_adj = c(x1 = 2L),
      data_adj_list = list(x1 = cached)
    )
  )

  out <- build_adjustment_step(
    x = x,
    xi = "x2",
    powers_current = list(x1 = c(1), x2 = c(1)),
    powers = list(x1 = c(-2, -1, 0, 1, 2), x2 = c(-2, -1, 0, 1, 2)),
    acdx = c(x1 = FALSE, x2 = FALSE),
    zero = c(x1 = FALSE, x2 = FALSE),
    catzero = list(x1 = NULL, x2 = NULL),
    spike = c(x1 = FALSE, x2 = FALSE),
    spike_decision = c(x1 = 2L, x2 = 2L),
    acd_parameter = list(x1 = NULL, x2 = NULL),
    prev_adj_params = prev_adj_params,
    transform_cache = list(x1 = NULL, x2 = NULL),
    term_to_columns = list(x1 = "x1", x2 = "x2")
  )

  expected <- cached
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
  expect_equal(out$data_adj_list$x1, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge reuses the cached
# matrix when normalized powers and spike decision are unchanged.
test_that("C++ adjustment-step bridge reuses cache on cache hit", {
  x <- cbind(x1 = c(1, 2, 4))

  cached <- matrix(c(99, 98, 97), ncol = 1L)
  colnames(cached) <- "old_name"

  prev_xi <- list(
    data_adj_list = list(x1 = cached),
    spike_decision_adj = c(x1 = 2L)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = NULL),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = list(x1 = c(1)),
    prev_xi = prev_xi,
    has_prev = TRUE
  )

  expected <- cached
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected)
  expect_equal(out$data_adj_list$x1, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge recomputes the matrix
# when the normalized power key changes.
test_that("C++ adjustment-step bridge recomputes on power-key cache miss", {
  x <- cbind(x1 = c(1, 2, 4))

  cached <- matrix(c(99, 98, 97), ncol = 1L)

  prev_xi <- list(
    data_adj_list = list(x1 = cached),
    spike_decision_adj = c(x1 = 2L)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(0)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = NULL),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(0)),
    prev_power_keys_adj = list(x1 = c(1)),
    prev_xi = prev_xi,
    has_prev = TRUE
  )

  expected <- matrix(log(x[, "x1"]), ncol = 1L)
  colnames(expected) <- "x1_adj1"

  expect_equal(out$data_adj, expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ adjustment-step bridge recomputes the matrix
# when the spike decision changes.
test_that("C++ adjustment-step bridge recomputes on spike-decision cache miss", {
  x <- cbind(x1 = c(1, 2, 4))
  cz <- matrix(c(0, 1, 0), ncol = 1L)

  cached <- matrix(c(99, 98, 97), ncol = 1L)

  prev_xi <- list(
    data_adj_list = list(x1 = cached),
    spike_decision_adj = c(x1 = 2L)
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1)),
    acdx_adj = list(x1 = FALSE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = cz),
    spike_adj = list(x1 = TRUE),
    spike_decision_int_adj = c(x1 = 1L),
    acd_parameter_adj = list(x1 = NULL),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1)),
    prev_power_keys_adj = list(x1 = c(1)),
    prev_xi = prev_xi,
    has_prev = TRUE
  )

  expected <- cbind(cz[, 1], x[, "x1"])
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj, expected)
})


# Test purpose: Checks that the C++ adjustment-step bridge applies stored ACD
# parameters without fitting ACD inside build_adjustment_step().
test_that("C++ adjustment-step bridge builds stored-ACD adjustment columns", {
  x <- cbind(x1 = c(1, 2, 4, 8))

  acd_parameter <- list(
    beta0 = 0,
    beta1 = 1,
    power = 1,
    shift = 0,
    scale = 1
  )

  out <- mfp2_build_adjustment_step_loop(
    x = x,
    vars_adj = "x1",
    powers_adj = list(x1 = c(1, 1)),
    acdx_adj = list(x1 = TRUE),
    zero_adj = list(x1 = FALSE),
    catzero = list(x1 = NULL),
    spike_adj = list(x1 = FALSE),
    spike_decision_int_adj = c(x1 = 2L),
    acd_parameter_adj = list(x1 = acd_parameter),
    eliminated = c(x1 = FALSE),
    spike_binary_only_flags = c(x1 = FALSE),
    current_power_keys_adj = list(x1 = c(1, 1)),
    prev_power_keys_adj = NULL,
    prev_xi = NULL,
    has_prev = FALSE
  )

  expected <- cbind(
    x[, "x1"],
    stats::pnorm(x[, "x1"])
  )
  colnames(expected) <- c("x1_adj1", "x1_adj2")

  expect_equal(out$data_adj, expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ adjustment-step bridge rejects ACD variables
# without stored ACD parameters.
test_that("C++ adjustment-step bridge rejects missing stored ACD parameters", {
  x <- cbind(x1 = c(1, 2, 4, 8))

  expect_error(
    mfp2_build_adjustment_step_loop(
      x = x,
      vars_adj = "x1",
      powers_adj = list(x1 = c(1, 1)),
      acdx_adj = list(x1 = TRUE),
      zero_adj = list(x1 = FALSE),
      catzero = list(x1 = NULL),
      spike_adj = list(x1 = FALSE),
      spike_decision_int_adj = c(x1 = 2L),
      acd_parameter_adj = list(x1 = NULL),
      eliminated = c(x1 = FALSE),
      spike_binary_only_flags = c(x1 = FALSE),
      current_power_keys_adj = list(x1 = c(1, 1)),
      prev_power_keys_adj = NULL,
      prev_xi = NULL,
      has_prev = FALSE
    ),
    "missing stored|require stored|acd_parameter",
    ignore.case = TRUE
  )
})
