# Tests for catzero indicator naming and reuse in the linear reference design.

# The reference design must derive <term>_bin names directly from the active
# conceptual terms, extend each conceptual block with its binary indicator, and
# leave the stored backfitting blocks unchanged.
test_that("catzero reference indicators use source-term _bin names", {
  x <- cbind(
    exposure = c(0, 1, 2, 3),
    dose = c(0, 0, 2, 4),
    z = c(1, 2, 3, 4)
  )

  exposure_indicator <- matrix(
    c(1, 0, 0, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )
  dose_indicator <- matrix(
    c(1, 1, 0, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )

  blocks <- list(
    exposure = exposure_indicator,
    dose = dose_indicator,
    z = NULL
  )
  mapping <- stats::setNames(as.list(colnames(x)), colnames(x))

  design <- assemble_linear_reference_design(
    x = x,
    term_to_columns = mapping,
    catzero_blocks = blocks,
    intercept = TRUE
  )

  expect_identical(
    colnames(design$x),
    c("(Intercept)", "exposure", "dose", "z", "exposure_bin", "dose_bin")
  )
  expect_identical(
    design$term_to_columns$exposure,
    c("exposure", "exposure_bin")
  )
  expect_identical(
    design$term_to_columns$dose,
    c("dose", "dose_bin")
  )
  expect_identical(design$term_to_columns$z, "z")

  expect_equal(design$x[, "exposure_bin"], exposure_indicator[, 1L])
  expect_equal(design$x[, "dose_bin"], dose_indicator[, 1L])

  # The reusable backfitting blocks must not be renamed in place.
  expect_identical(colnames(blocks$exposure), "catzero")
  expect_identical(colnames(blocks$dose), "catzero")
})

# The helper must use the already-prepared indicator values. This deliberately
# supplies values that cannot be inferred from x; recomputing I(x <= 0) would
# therefore fail the test.
test_that("catzero reference design reuses prepared indicator values", {
  x <- cbind(
    exposure = c(0, 0, 1, 2),
    z = c(4, 3, 2, 1)
  )
  prepared_indicator <- matrix(
    c(0, 1, 1, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )

  design <- assemble_linear_reference_design(
    x = x,
    term_to_columns = list(exposure = "exposure", z = "z"),
    catzero_blocks = list(exposure = prepared_indicator, z = NULL),
    intercept = FALSE
  )

  expect_equal(design$x[, "exposure_bin"], prepared_indicator[, 1L])
})

# Tests for generated catzero <term>_bin name collisions.

# A generated binary indicator name must never silently overwrite or be made
# unique against an existing model-matrix column.
test_that("catzero reference design rejects an existing _bin column", {
  x <- cbind(
    exposure = c(0, 0, 1, 2),
    exposure_bin = c(0, 1, 0, 1)
  )
  indicator <- matrix(
    c(1, 1, 0, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )

  expect_error(
    assemble_linear_reference_design(
      x = x,
      term_to_columns = stats::setNames(as.list(colnames(x)), colnames(x)),
      catzero_blocks = list(exposure = indicator, exposure_bin = NULL),
      intercept = TRUE
    ),
    "Generated catzero indicator name.*exposure_bin"
  )
})

# Exercise the same collision through the public fit path so a future refactor
# cannot bypass the helper-level protection.
test_that("mfp2 rejects catzero _bin name collisions before reference fitting", {
  n <- 60L
  exposure <- c(rep(0, 20L), seq(0.2, 8, length.out = 40L))
  x <- cbind(
    exposure = exposure,
    exposure_bin = rep(c(0, 1), length.out = n)
  )
  y <- 1 + 0.5 * exposure + stats::rnorm(n, sd = 0.2)

  expect_error(
    mfp2(
      x,
      y,
      catzero_vars = "exposure",
      xorder = "original",
      verbose = FALSE
    ),
    "Generated catzero indicator name.*exposure_bin"
  )
})

# Duplicate conceptual source names would generate duplicate <term>_bin names
# and must be rejected explicitly rather than producing an ambiguous mapping.
test_that("catzero reference design rejects duplicate generated indicator names", {
  x <- cbind(exposure = c(0, 0, 1, 2))
  indicator <- matrix(
    c(1, 1, 0, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )
  blocks <- list(indicator, indicator)
  names(blocks) <- c("exposure", "exposure")

  expect_error(
    assemble_linear_reference_design(
      x = x,
      term_to_columns = list(exposure = "exposure"),
      catzero_blocks = blocks,
      intercept = FALSE
    ),
    "duplicate catzero indicator names"
  )
})

# Tests for malformed catzero reference blocks and mappings.

# Named blocks are required because the conceptual source term determines both
# the generated <term>_bin name and the blockwise LRT mapping.
test_that("catzero reference blocks must be named", {
  x <- cbind(exposure = c(0, 0, 1, 2))
  indicator <- matrix(c(1, 1, 0, 0), ncol = 1L)

  expect_error(
    assemble_linear_reference_design(
      x = x,
      term_to_columns = list(exposure = "exposure"),
      catzero_blocks = list(indicator),
      intercept = FALSE
    ),
    "catzero reference blocks must be named"
  )
})

# Every active catzero block must correspond to an existing conceptual term.
test_that("catzero reference blocks must match conceptual terms", {
  x <- cbind(exposure = c(0, 0, 1, 2))
  indicator <- matrix(
    c(1, 1, 0, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )

  expect_error(
    assemble_linear_reference_design(
      x = x,
      term_to_columns = list(exposure = "exposure"),
      catzero_blocks = list(other = indicator),
      intercept = FALSE
    ),
    "do not match conceptual terms"
  )
})

# A structural-zero indicator is exactly one column. Multi-column input must not
# be accepted because it would change the conceptual term df unexpectedly.
test_that("catzero reference block must contain exactly one column", {
  x <- cbind(exposure = c(0, 0, 1, 2))
  bad_block <- cbind(
    catzero = c(1, 1, 0, 0),
    extra = c(0, 0, 1, 1)
  )

  expect_error(
    assemble_linear_reference_design(
      x = x,
      term_to_columns = list(exposure = "exposure"),
      catzero_blocks = list(exposure = bad_block),
      intercept = FALSE
    ),
    "must have one column"
  )
})

# Indicator rows must remain exactly aligned with x.
test_that("catzero reference block must match x row count", {
  x <- cbind(exposure = c(0, 0, 1, 2))
  bad_block <- matrix(
    c(1, 1, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )

  expect_error(
    assemble_linear_reference_design(
      x = x,
      term_to_columns = list(exposure = "exposure"),
      catzero_blocks = list(exposure = bad_block),
      intercept = FALSE
    ),
    "match x rows"
  )
})

# A non-NULL list containing no active indicator blocks should remain a valid
# defensive no-op and should not modify the conceptual mapping.
test_that("empty catzero reference blocks are a no-op", {
  x <- cbind(exposure = 1:4, z = 4:1)
  mapping <- stats::setNames(as.list(colnames(x)), colnames(x))

  design <- assemble_linear_reference_design(
    x = x,
    term_to_columns = mapping,
    catzero_blocks = list(exposure = NULL, z = NULL),
    intercept = FALSE
  )

  expect_identical(design$x, x)
  expect_identical(design$term_to_columns, mapping)
})
