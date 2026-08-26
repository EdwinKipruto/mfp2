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
# supplies values that cannot be inferred from x; recomputing I(x == 0) would
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


# Migrated coverage from the former test_mfp2.R

# Test purpose: When no catzero/SAZ indicator is active, order_variables() must
# stay on the direct x path and must not enter reference-block assembly.
test_that("ordinary and zero-only references bypass catzero assembly", {
  x <- cbind(
    x1 = c(0, 0, 1, 2),
    x2 = c(4, 3, 2, 1)
  )
  mapping <- stats::setNames(as.list(colnames(x)), colnames(x))

  testthat::local_mocked_bindings(
    assemble_linear_reference_design = function(...) {
      stop("catzero assembly must not run when catzero_blocks is NULL")
    },
    fit_full_linear_reference = function(x, ...) {
      list(
        null_deviance = 50,
        model_deviance = 40,
        logl = -20,
        df = 3,
        null_logl = NA_real_
      )
    },
    order_variables_by_significance = function(
    xorder, x, term_to_columns, full_reference, ...) {
      expect_identical(term_to_columns, mapping)
      expect_equal(x[, "x1"], c(0, 0, 1, 2))
      c("x1", "x2")
    },
    .package = "mfp2"
  )

  result <- order_variables(
    xorder = "ascending",
    x = x,
    term_to_columns = mapping,
    catzero_blocks = NULL,
    y = rep(0, nrow(x)),
    family = stats::gaussian(),
    family_string = "gaussian"
  )

  expect_identical(result$variables_ordered, c("x1", "x2"))
  expect_equal(result$linear_deviance, 40)
})


# Test purpose: The full reference and the ordering tests must use the same
# two-column conceptual block for a catzero/retained SAZ term.
test_that("catzero full reference and ordering share one joint block", {
  x <- cbind(
    exposure = c(0, 0, 1, 2),
    z = c(1, 2, 3, 4)
  )
  indicator <- matrix(
    c(1, 1, 0, 0),
    ncol = 1L,
    dimnames = list(NULL, "catzero")
  )
  blocks <- list(exposure = indicator, z = NULL)
  mapping <- list(exposure = "exposure", z = "z")

  full_seen <- NULL
  ordering_seen <- NULL

  testthat::local_mocked_bindings(
    fit_full_linear_reference = function(x, ...) {
      full_seen <<- x
      list(
        null_deviance = 50,
        model_deviance = 30,
        logl = -15,
        df = 4,
        null_logl = NA_real_
      )
    },
    order_variables_by_significance = function(
    xorder, x, term_to_columns, full_reference, ...) {
      ordering_seen <<- list(
        x = x,
        term_to_columns = term_to_columns,
        full_reference = full_reference
      )
      c("exposure", "z")
    },
    .package = "mfp2"
  )

  result <- order_variables(
    xorder = "ascending",
    x = x,
    term_to_columns = mapping,
    catzero_blocks = blocks,
    y = rep(0, nrow(x)),
    family = stats::gaussian(),
    family_string = "gaussian"
  )

  expect_equal(NCOL(full_seen), 4L) # intercept + x + indicator
  expect_identical(ordering_seen$x, full_seen)
  expect_length(ordering_seen$term_to_columns$exposure, 2L)
  expect_identical(ordering_seen$full_reference$logl, -15)
  expect_equal(result$linear_deviance, 30)
  expect_equal(result$linear_df, 4)
})


# Test purpose: The stored full-linear deviance for a zero term must use the
# positive-part reference column rather than the unre-coded covariate.
test_that("zero term changes the full linear reference representation", {
  x <- cbind(
    exposure = c(0, 0, 0.5, 1, 2, 4, 5, 6),
    z = c(-1, 0, 1, 0, -1, 1, 0.5, -0.5)
  )
  y <- c(1.0, 1.2, 0.8, 1.6, 2.1, 2.8, 3.0, 3.5)

  fit <- mfp2(
    x,
    y,
    zero_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  )

  reference_data <- data.frame(
    y = y,
    exposure = pmax(x[, "exposure"], 0),
    z = x[, "z"]
  )
  reference <- stats::glm(
    y ~ exposure + z,
    data = reference_data,
    family = stats::gaussian()
  )

  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})


# Test purpose: The stored full-linear deviance for catzero must include both
# the positive-part continuous column and the structural-zero indicator.
test_that("catzero term changes the full linear reference representation", {
  x <- cbind(
    exposure = c(0, 0, 0.5, 1, 2, 4, 5, 6),
    z = c(-1, 0, 1, 0, -1, 1, 0.5, -0.5)
  )
  y <- c(3.0, 3.2, 2.8, 1.6, 2.1, 2.8, 3.0, 3.5)

  fit <- mfp2(
    x,
    y,
    catzero_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  )

  reference_data <- data.frame(
    y = y,
    exposure = x[, "exposure"],
    zero_indicator = as.integer(x[, "exposure"] == 0),
    z = x[, "z"]
  )
  reference <- stats::glm(
    y ~ exposure + zero_indicator + z,
    data = reference_data,
    family = stats::gaussian()
  )

  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})


# Test purpose: The leave-one-term-out likelihood-ratio test used for
# significance ordering must use the fitted rank contribution of a grouped
# term, not the number of raw columns in its design block. The mocked fits are
# chosen so those two df rules produce opposite orders.
test_that("grouped significance ordering uses fitted rank difference", {
  x <- cbind(
    groupB = c(0, 1, 0, 1),
    groupC = c(0, 0, 1, 1),
    z = c(-1, 0, 1, 2)
  )
  term_to_columns <- list(
    group = c("groupB", "groupC"),
    z = "z"
  )

  testthat::local_mocked_bindings(
    fit_model = function(x, ...) {
      remaining <- colnames(x)
      if (identical(remaining, "z")) {
        # Dropping the two-column group reduces fitted rank by one and gives
        # likelihood-ratio statistic 20.
        return(list(logl = 90, df = 1))
      }
      if (identical(remaining, c("groupB", "groupC"))) {
        # Dropping z also reduces fitted rank by one and gives statistic 18.
        return(list(logl = 91, df = 1))
      }
      stop("Unexpected reduced design in ordering test.")
    },
    .package = "mfp2"
  )

  common_args <- list(
    x = x,
    y = rep(0, nrow(x)),
    family = stats::gaussian(),
    family_string = "gaussian",
    weights = NULL,
    offset = NULL,
    strata = NULL,
    method = NULL,
    control = NULL,
    nocenter = NULL,
    full_reference = list(logl = 100, df = 2),
    term_to_columns = term_to_columns
  )

  ascending <- do.call(
    order_variables_by_significance,
    c(list(xorder = "ascending"), common_args)
  )
  descending <- do.call(
    order_variables_by_significance,
    c(list(xorder = "descending"), common_args)
  )

  expect_identical(ascending, c("group", "z"))
  expect_identical(descending, c("z", "group"))
})


# =============================================================================
# 14. zero_vars and catzero_vars
# =============================================================================

# Test purpose: Checks that zero_vars activates zero-component handling for
# exact-zero values.
test_that("zero_vars preserves exact-zero values", {
  set.seed(1)
  n <- 200
  x_val <- c(rep(0, 40), rgamma(n - 40, shape = 2, rate = 1))
  x_mat <- cbind(exposure = x_val, x2 = runif(n, 1, 10))
  y_val <- 2 * x_val + rnorm(n)

  fit <- mfp2(x_mat, y_val, zero_vars = "exposure", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_true(fit$zero["exposure"])
})


# Test purpose: Checks that catzero_vars adds a zero-component indicator and
# implies zero handling.
test_that("catzero_vars creates binary indicator", {
  set.seed(1)
  n <- 200
  x_val <- c(rep(0, 40), rgamma(n - 40, shape = 2, rate = 1))
  x_mat <- cbind(exposure = x_val, x2 = runif(n, 1, 10))
  y_val <- 2 * x_val + 1.5 * (x_val == 0) + rnorm(n)

  fit <- mfp2(x_mat, y_val, catzero_vars = "exposure", verbose = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_true(fit$catzero["exposure"])
  # zero should also be TRUE (catzero implies zero)
  expect_true(fit$zero["exposure"])
})


# Test purpose: Verifies that fp(x, zero = TRUE) is converted to zero_vars
# and that exact-zero values are handled through the zero component.
test_that("formula interface fp(zero = TRUE) enables zero handling", {
  set.seed(102)
  n <- 200

  exposure <- c(rep(0, 40), rgamma(n - 40, shape = 2, rate = 1))
  dat <- data.frame(
    y = 2 * exposure + rnorm(n),
    exposure = exposure,
    x2 = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ fp(exposure, zero = TRUE) + fp(x2),
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$zero["exposure"])
})


# Test purpose: Verifies that fp(x, catzero = TRUE) creates a zero-component
# indicator and also implies zero handling.
test_that("formula interface fp(catzero = TRUE) enables catzero and zero handling", {
  set.seed(103)
  n <- 200

  exposure <- c(rep(0, 40), rgamma(n - 40, shape = 2, rate = 1))
  dat <- data.frame(
    y = 2 * exposure + 1.5 * (exposure == 0) + rnorm(n),
    exposure = exposure,
    x2 = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ fp(exposure, catzero = TRUE) + fp(x2),
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$catzero["exposure"])
  expect_true(fit$zero["exposure"])
})
