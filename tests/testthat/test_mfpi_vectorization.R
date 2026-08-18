# MFPI compact-basis/vectorization correctness tests
#
# These tests validate the canonical compact MFPI implementation directly. They
# do not depend on a second candidate-construction implementation. Expected
# interaction matrices and centering constants are assembled explicitly from the
# compact basis using the mathematical definition of the MFPI group blocks.


# Build the expected group-specific candidate matrix in simple R code.
#
# `selected_basis` is n x d. The result is group-major:
#   [group 1 term 1..d | group 2 term 1..d | ...].
# Invalid rows represent structural zeros and remain exactly zero.
expected_mfpi_candidate <- function(selected_basis,
                                    group_idx,
                                    n_groups,
                                    valid_rows,
                                    center = TRUE,
                                    group_center = FALSE) {
  selected_basis <- as.matrix(selected_basis)
  n <- nrow(selected_basis)
  n_terms <- ncol(selected_basis)

  out <- matrix(0, nrow = n, ncol = n_groups * n_terms)
  centers <- numeric(n_groups * n_terms)

  if (center) {
    if (group_center) {
      for (g in seq_len(n_groups)) {
        rows_g <- group_idx == g & valid_rows
        expect_true(any(rows_g))
        cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)
        centers[cols_g] <- colMeans(selected_basis[rows_g, , drop = FALSE])
      }
    } else {
      grand <- colMeans(selected_basis[valid_rows, , drop = FALSE])
      centers <- rep(grand, times = n_groups)
    }
  }

  for (g in seq_len(n_groups)) {
    rows_g <- group_idx == g & valid_rows
    cols_g <- seq.int((g - 1L) * n_terms + 1L, g * n_terms)

    if (any(rows_g)) {
      values <- selected_basis[rows_g, , drop = FALSE]
      if (center) {
        values <- sweep(values, 2L, centers[cols_g], "-", check.margin = FALSE)
      }
      out[rows_g, cols_g] <- values
    }
  }

  list(target = out, centers = centers)
}


test_that("MFPI compact scatter computes grand and group centered FP2 candidates", {
  x <- c(1, 2, 4, 8, 16, 32)
  group_idx <- c(1L, 1L, 1L, 2L, 2L, 2L)
  n_groups <- 2L
  powers <- generate_powers_fp(degree = 2L, powers = c(0, 1))

  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = FALSE
  )

  n_terms <- ncol(compact$candidate_map)
  focal_width <- n_groups * n_terms
  valid_rows <- rep(TRUE, length(x))

  for (group_center in c(FALSE, TRUE)) {
    target <- matrix(0, nrow = length(x), ncol = focal_width)

    for (i in seq_len(nrow(compact$candidate_map))) {
      source_cols <- as.integer(compact$candidate_map[i, ])
      selected_basis <- compact$basis[, source_cols, drop = FALSE]

      expected <- expected_mfpi_candidate(
        selected_basis = selected_basis,
        group_idx = group_idx,
        n_groups = n_groups,
        valid_rows = valid_rows,
        center = TRUE,
        group_center = group_center
      )

      built <- fill_mfpi_fp_candidate_cpp(
        target = target,
        basis = compact$basis,
        source_cols = source_cols,
        group_idx = group_idx,
        n_groups = n_groups,
        target_start_col = 1L,
        center = TRUE,
        group_center = group_center,
        valid_rows = valid_rows
      )
      target <- built$target

      expect_equal(unname(target), expected$target, tolerance = 1e-12)
      expect_equal(unname(built$centers), expected$centers, tolerance = 1e-12)
    }
  }
})


test_that("MFPI compact scatter preserves structural zeros", {
  x <- c(0, 1, 2, 0, 4, 8)
  group_idx <- c(1L, 1L, 1L, 2L, 2L, 2L)
  n_groups <- 2L
  valid_rows <- x > 0
  powers <- generate_powers_fp(degree = 2L, powers = c(0, 1))

  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = TRUE
  )

  n_terms <- ncol(compact$candidate_map)
  focal_width <- n_groups * n_terms
  target <- matrix(0, nrow = length(x), ncol = focal_width)

  for (i in seq_len(nrow(compact$candidate_map))) {
    source_cols <- as.integer(compact$candidate_map[i, ])
    selected_basis <- compact$basis[, source_cols, drop = FALSE]

    expected <- expected_mfpi_candidate(
      selected_basis = selected_basis,
      group_idx = group_idx,
      n_groups = n_groups,
      valid_rows = valid_rows,
      center = TRUE,
      group_center = TRUE
    )

    built <- fill_mfpi_fp_candidate_cpp(
      target = target,
      basis = compact$basis,
      source_cols = source_cols,
      group_idx = group_idx,
      n_groups = n_groups,
      target_start_col = 1L,
      center = TRUE,
      group_center = TRUE,
      valid_rows = valid_rows
    )
    target <- built$target

    expect_equal(unname(target), expected$target, tolerance = 1e-12)
    expect_equal(unname(built$centers), expected$centers, tolerance = 1e-12)
    expect_true(all(target[!valid_rows, , drop = FALSE] == 0))
  }
})


test_that("MFPI compact scatter leaves candidates uncentered when requested", {
  x <- c(1, 2, 4, 8)
  group_idx <- c(1L, 1L, 2L, 2L)
  powers <- matrix(c(0, 1), nrow = 1L)

  compact <- generate_transformations_fp_basis_cpp(
    x = x,
    powers = powers,
    zero = FALSE
  )

  source_cols <- as.integer(compact$candidate_map[1, ])
  selected_basis <- compact$basis[, source_cols, drop = FALSE]
  expected <- expected_mfpi_candidate(
    selected_basis = selected_basis,
    group_idx = group_idx,
    n_groups = 2L,
    valid_rows = rep(TRUE, length(x)),
    center = FALSE
  )

  built <- fill_mfpi_fp_candidate_cpp(
    target = matrix(0, nrow = length(x), ncol = 4L),
    basis = compact$basis,
    source_cols = source_cols,
    group_idx = group_idx,
    n_groups = 2L,
    target_start_col = 1L,
    center = FALSE,
    group_center = FALSE,
    valid_rows = rep(TRUE, length(x))
  )

  expect_equal(unname(built$target), expected$target, tolerance = 1e-12)
  expect_equal(unname(built$centers), numeric(4L), tolerance = 1e-12)
})


test_that("mfpi flex2 fits FP2 interactions through the compact path", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "fp2"),
    group_var = "svi",
    flex = "flex2",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$flex, "flex2")
})


test_that("mfpi flex4 supports direct group centering", {
  data("prostate", package = "mfp2")

  fit <- mfpi(
    lpsa ~ fp(age) + svi + fp(cavol),
    data = prostate,
    cont_vars = "cavol",
    cont_var_forms = c(cavol = "fp1"),
    group_var = "svi",
    flex = "flex4",
    center_type = "group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfpi")
  expect_equal(fit$flex, "flex4")
  expect_equal(fit$center_type, "group")
})
