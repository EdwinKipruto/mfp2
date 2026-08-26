# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Checks that the C++ FP core returns the expected ordinary
# fractional-polynomial columns for powers 1 and 0.
# Test purpose: Checks that the C++ FP core returns the expected ordinary
# fractional-polynomial columns for powers 1 and 0.
test_that("C++ FP core computes ordinary FP transformations", {
  x <- c(1, 2, 4)
  p <- c(1, 0)

  out <- transform_fp_core(
    x_raw = x,
    power = p,
    shift_val = 0,
    scale_val = 1,
    zero = FALSE
  )

  expected <- cbind(
    x,
    log(x)
  )

  expect_true(is.matrix(out))
  expect_equal(
    unname(out, force = TRUE),
    unname(expected, force = TRUE),
    tolerance = 1e-12
  )
})


# Test purpose: Checks that the C++ FP core implements repeated-power FP rules.
test_that("C++ FP core handles repeated powers", {
  x <- c(1, 2, 4)
  p <- c(0, 0)

  out <- transform_fp_core(
    x_raw = x,
    power = p,
    shift_val = 0,
    scale_val = 1,
    zero = FALSE
  )

  expected <- cbind(
    log(x),
    log(x)^2
  )

  expect_true(is.matrix(out))
  expect_equal(unname(out), expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ FP core applies shift and scale exactly once.
test_that("C++ FP core applies shift and scale", {
  x <- c(1, 3, 5)
  shift <- 1
  scale <- 2
  x_scaled <- (x + shift) / scale

  out <- transform_fp_core(
    x_raw = x,
    power = c(1, 0),
    shift_val = shift,
    scale_val = scale,
    zero = FALSE
  )

  expected <- cbind(
    x_scaled,
    log(x_scaled)
  )

  expect_equal(
    unname(out, force = TRUE),
    unname(expected, force = TRUE),
    tolerance = 1e-12
  )
})


# Test purpose: Checks that the C++ FP core maps exact-zero values to zero
# rows when zero handling is active.
test_that("C++ FP core handles zero mode", {
  x <- c(0, 1, 4)

  out <- transform_fp_core(
    x_raw = x,
    power = c(1, 0),
    shift_val = 0,
    scale_val = 1,
    zero = TRUE
  )

  expected <- rbind(
    c(0, 0),
    c(1, log(1)),
    c(4, log(4))
  )

  expect_equal(unname(out), expected, tolerance = 1e-12)
})


# Test purpose: Checks that the C++ FP core preserves missing and non-finite
# values. Keep this test only if transform_fp_core_internal() has the explicit
# missing/non-finite guard.
test_that("C++ FP core preserves missing and non-finite values", {
  x <- c(1, NA_real_, NaN, Inf, 4)
  p <- c(1, 0)

  out <- transform_fp_core(
    x_raw = x,
    power = p,
    shift_val = 0,
    scale_val = 1,
    zero = FALSE
  )

  expect_equal(out[1, ], c(1, 0))
  expect_true(is.na(out[2, 1]))
  expect_true(is.na(out[2, 2]))
  expect_true(is.nan(out[3, 1]))
  expect_true(is.nan(out[3, 2]))
  expect_true(is.infinite(out[4, 1]))
  expect_true(is.infinite(out[4, 2]))
  expect_equal(out[5, ], c(4, log(4)))
})


# Test purpose: Checks that the R wrapper around the C++ FP core preserves
# variable naming.
test_that("transform_vector_fp keeps expected column names with C++ core", {
  x <- c(1, 2, 4)

  out <- transform_vector_fp(
    x = x,
    power = c(1, 0),
    shift = 0,
    scale = 1,
    name = "x",
    zero = FALSE,
    check_binary = FALSE
  )

  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2L)
  expect_equal(unname(out[, 2]), x)
  expect_equal(unname(out[, 1]), log(x))
  expect_false(is.null(colnames(out)))
})


# Test purpose: Checks that the C++ batch FP generator returns one matrix per
# candidate power row.
test_that("generate_transformations_fp_cpp returns one matrix per power row", {
  x <- c(1, 2, 4)
  powers <- rbind(
    c(1, 1),
    c(0, 0),
    c(1, 0)
  )

  out <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = FALSE,
    catzero = NULL
  )

  expect_type(out, "list")
  expect_length(out, nrow(powers))

  expect_equal(
    unname(out[[1]]),
    unname(cbind(x, x * log(x))),
    tolerance = 1e-12
  )

  expect_equal(
    unname(out[[2]]),
    cbind(log(x), log(x)^2),
    tolerance = 1e-12
  )

  expect_equal(
    unname(out[[3]]),
    unname(cbind(x, log(x))),
    tolerance = 1e-12
  )
})


# Test purpose: Checks that the C++ batch FP generator applies the binary
# shortcut when zero handling is inactive.
test_that("generate_transformations_fp_cpp uses binary shortcut", {
  x <- c(0, 1, 0, 1)
  powers <- rbind(
    c(1, 1),
    c(0, 0)
  )

  out <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = FALSE,
    catzero = NULL
  )

  expect_length(out, 2L)
  expect_equal(unname(out[[1]]), matrix(x, ncol = 1L))
  expect_equal(unname(out[[2]]), matrix(x, ncol = 1L))
})


# Test purpose: Checks that the C++ batch FP generator prepends catzero when
# catzero is supplied.
test_that("generate_transformations_fp_cpp prepends catzero", {
  x <- c(1, 2, 4)
  powers <- rbind(c(1, 0))
  catzero <- matrix(c(0, 1, 0), ncol = 1L)

  out <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = FALSE,
    catzero = catzero
  )

  expect_length(out, 1L)
  expect_equal(ncol(out[[1]]), 3L)
  expect_equal(unname(out[[1]][, 1]), as.numeric(catzero[, 1]))
  expect_equal(unname(out[[1]][, 2]), x)
  expect_equal(unname(out[[1]][, 3]), log(x), tolerance = 1e-12)
  expect_equal(colnames(out[[1]]), c("catzero", "V1", "V2"))
})


# Test purpose: Checks that the compact ordinary-FP basis reconstructs every
# degree-2 candidate exactly while storing only the unique repeated-power terms.
test_that("compact FP basis reconstructs materialized degree-2 candidates", {
  x <- c(1, 2, 4, 8, 16)
  allowed_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  materialized <- generate_transformations_fp(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  compact <- generate_transformations_fp_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE
  )

  expect_equal(nrow(compact$candidate_map), 36L)
  expect_equal(ncol(compact$candidate_map), 2L)
  expect_equal(ncol(compact$basis), 16L)
  expect_equal(compact$powers, materialized$powers)

  for (i in seq_len(nrow(compact$candidate_map))) {
    expect_equal(
      unname(materialize_fp_basis_candidate(compact, i)),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Checks repeated zero powers and structural-zero indicators in
# the compact representation, including the column names required by SAZ stage 2.
test_that("compact FP basis preserves repeated zero powers and catzero", {
  x <- c(1, 2, 4, 8)
  catzero <- matrix(c(0, 1, 0, 1), ncol = 1L)

  compact <- generate_transformations_fp_basis(
    x = x,
    degree = 2,
    powers = 0,
    zero = FALSE,
    catzero = catzero
  )

  out <- materialize_fp_basis_candidate(compact, 1L)

  expect_equal(ncol(compact$basis), 2L)
  expect_equal(unname(compact$basis[, 1L]), log(x), tolerance = 1e-12)
  expect_equal(unname(compact$basis[, 2L]), log(x)^2, tolerance = 1e-12)
  expect_equal(unname(out[, 1L]), as.numeric(catzero[, 1L]))
  expect_equal(unname(out[, 2L]), log(x), tolerance = 1e-12)
  expect_equal(unname(out[, 3L]), log(x)^2, tolerance = 1e-12)
  expect_equal(colnames(out), c("catzero", "V1", "V2"))
})


# Test purpose: Checks exact-zero semantics in the compact basis. Exact-zero
# values must remain zero while positive observations use ordinary FP terms.
test_that("compact FP basis preserves zero-mode transformations", {
  x <- c(0, 0.5, 1, 2, 4)
  powers <- rbind(c(-1, -1), c(0, 0), c(1, 2))

  old <- generate_transformations_fp_cpp(
    x = x,
    powers = powers,
    zero = TRUE,
    catzero = NULL
  )
  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = TRUE
  )

  for (i in seq_len(nrow(powers))) {
    reconstructed <- compact$basis[
      , compact$candidate_map[i, ], drop = FALSE
    ]
    expect_equal(
      unname(reconstructed),
      unname(old[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Checks that compact candidate columns can be copied into a
# reusable design matrix without allocating/materializing the complete candidate.
test_that("compact FP candidate copies into reusable design matrix", {
  x <- c(1, 2, 4, 8)
  powers <- rbind(c(0, 0), c(1, 2))
  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = FALSE
  )

  target <- matrix(-1, nrow = length(x), ncol = 4L)
  target[, 1L] <- 1
  target[, 4L] <- 99

  target <- copy_fp_basis_candidate_cpp(
    target = target,
    basis = compact$basis,
    source_cols = as.integer(compact$candidate_map[1L, ]),
    target_cols = c(2L, 3L)
  )

  expect_equal(target[, 1L], rep(1, length(x)))
  expect_equal(target[, 2L], log(x), tolerance = 1e-12)
  expect_equal(target[, 3L], log(x)^2, tolerance = 1e-12)
  expect_equal(target[, 4L], rep(99, length(x)))
})
