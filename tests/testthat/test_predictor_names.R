# Predictor-name validation ---------------------------------------------------

test_that("predictor-name helper accepts supported non-syntactic names", {
  expect_invisible(
    validate_predictor_names(
      c("age", "age years", "tumour-size", "y", "offset_", "strata_", "..mfp2_predictor")
    )
  )
})

test_that("predictor-name helper rejects missing empty duplicate and backtick names", {
  expect_error(
    validate_predictor_names(NULL),
    "must have column names"
  )

  expect_error(
    validate_predictor_names(c("x1", NA_character_)),
    "must not be missing"
  )

  expect_error(
    validate_predictor_names(c("x1", "")),
    "must not be empty"
  )

  expect_error(
    validate_predictor_names(c("x1", "x1")),
    "must be unique"
  )

  expect_error(
    validate_predictor_names(c("x1", "age`years")),
    "must not contain backticks"
  )
})

test_that("mfp2.default enforces predictor-name contract before fitting", {
  set.seed(5101)
  n <- 40L
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  y <- 1 + 0.5 * x[, 1L] - 0.25 * x[, 2L] + rnorm(n, sd = 0.5)

  x_no_names <- x
  colnames(x_no_names) <- NULL
  expect_error(
    mfp2(x_no_names, y, df = 1, select = 1, verbose = FALSE),
    "must have column names"
  )

  x_na <- x
  colnames(x_na) <- c("x1", NA_character_)
  expect_error(
    mfp2(x_na, y, df = 1, select = 1, verbose = FALSE),
    "must not be missing"
  )

  x_empty <- x
  colnames(x_empty) <- c("x1", "")
  expect_error(
    mfp2(x_empty, y, df = 1, select = 1, verbose = FALSE),
    "must not be empty"
  )

  x_dup <- x
  colnames(x_dup) <- c("x1", "x1")
  expect_error(
    mfp2(x_dup, y, df = 1, select = 1, verbose = FALSE),
    "must be unique"
  )

  x_backtick <- x
  colnames(x_backtick) <- c("x1", "x`2")
  expect_error(
    mfp2(x_backtick, y, df = 1, select = 1, verbose = FALSE),
    "must not contain backticks"
  )
})

test_that("mfp2.default continues to support spaces and hyphens in predictor names", {
  set.seed(5102)
  n <- 50L
  x <- cbind(rnorm(n), rnorm(n))
  colnames(x) <- c("age years", "tumour-size")
  y <- 2 + 0.8 * x[, 1L] - 0.4 * x[, 2L] + rnorm(n, sd = 0.4)

  fit <- mfp2(
    x, y,
    df = 1,
    select = 1,
    keep = c("age years", "tumour-size"),
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
})

test_that("mfpi validates predictor names before resolving group_var", {
  set.seed(5103)
  n <- 40L
  x <- data.frame(
    group = rep(c(0, 1), each = n / 2),
    age = rnorm(n),
    z = rnorm(n),
    check.names = FALSE
  )
  y <- rnorm(n)

  x_dup <- x
  names(x_dup) <- c("group", "age", "age")
  expect_error(
    mfpi(
      x_dup, y,
      group_var = "group",
      cont_vars = "age",
      df = 1, select = 1, verbose = FALSE
    ),
    "must be unique"
  )

  x_backtick <- x
  names(x_backtick)[3L] <- "z`bad"
  expect_error(
    mfpi(
      x_backtick, y,
      group_var = "group",
      cont_vars = "age",
      df = 1, select = 1, verbose = FALSE
    ),
    "must not contain backticks"
  )

  x_empty <- x
  names(x_empty)[3L] <- ""
  expect_error(
    mfpi(
      x_empty, y,
      group_var = "group",
      cont_vars = "age",
      df = 1, select = 1, verbose = FALSE
    ),
    "must not be empty"
  )
})
