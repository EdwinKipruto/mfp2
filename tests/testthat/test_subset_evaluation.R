# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: The mfp2 formula interface evaluates a subset expression once
# in the formula data mask. A data column must take precedence over a same-named
# caller object.
test_that("24.10 mfp2 formula evaluates subset once with data precedence", {
  set.seed(2410)
  n <- 96L
  dat <- data.frame(
    y = 1 + 0.4 * seq_len(n) / n + rnorm(n, sd = 0.2),
    x = runif(n, 1, 5),
    keep = rep(c(TRUE, TRUE, FALSE), length.out = n)
  )
  keep <- rep(TRUE, n)
  evaluations <- 0L
  evaluate_once <- function(value) {
    evaluations <<- evaluations + 1L
    value
  }

  fit <- mfp2(
    y ~ x,
    data = dat,
    subset = evaluate_once(keep),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_identical(evaluations, 1L)
  expect_equal(nrow(fit$x_original), sum(dat$keep))
})


# Test purpose: The MFPI formula interface uses the same exact-once, data-first
# subset lookup as mfp2.formula().
test_that("24.11 MFPI formula evaluates subset once with data precedence", {
  set.seed(2411)
  n <- 120L
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.4 * x + 0.6 * x * (trt == "treated") + rnorm(n, sd = 0.3),
    trt = trt,
    x = x,
    z = rnorm(n),
    keep = rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)
  )
  keep <- rep(FALSE, n)
  evaluations <- 0L
  evaluate_once <- function(value) {
    evaluations <<- evaluations + 1L
    value
  }

  fit <- mfpi(
    y ~ trt + x + z,
    data = dat,
    subset = evaluate_once(keep),
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "z",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )

  expect_identical(evaluations, 1L)
  expect_equal(fit$nobs, sum(dat$keep))
})


# Test purpose: A formula subset need not be a data column. Numeric row indices
# defined in the formula environment remain valid in both formula interfaces.
test_that("24.12 formula methods accept caller-scoped subset row indices", {
  set.seed(2412)
  n <- 120L
  trt <- factor(rep(c("control", "treated"), length.out = n))
  x <- runif(n, 1, 5)
  dat <- data.frame(
    y = 1 + 0.3 * x + 0.5 * x * (trt == "treated") + rnorm(n, sd = 0.3),
    trt = trt,
    x = x,
    z = rnorm(n)
  )
  rows <- which(rep(c(TRUE, TRUE, FALSE), length.out = n))

  fit_mfp2 <- mfp2(
    y ~ x + z,
    data = dat,
    subset = rows,
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  expect_equal(nrow(fit_mfp2$x_original), length(rows))

  fit_mfpi <- mfpi(
    y ~ trt + x + z,
    data = dat,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    keep = "z",
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_equal(fit_mfpi$nobs, length(rows))
})


# Test purpose: Duplicate numeric row positions replicate observations and are
# therefore rejected consistently by all four public interfaces.
test_that("24.13 all subset interfaces reject duplicated numeric row indices", {
  set.seed(2413)
  n <- 120L
  trt_factor <- factor(rep(c("control", "treated"), length.out = n))
  trt_numeric <- as.numeric(trt_factor == "treated")
  x_cont <- runif(n, 1, 5)
  z <- rnorm(n)
  y <- 1 + 0.4 * x_cont + 0.5 * x_cont * trt_numeric + 0.2 * z +
    rnorm(n, sd = 0.3)
  duplicate_rows <- c(seq_len(80L), 80L)
  duplicate_message <- "must not contain duplicated row indices"

  x_matrix <- cbind(trt = trt_numeric, x = x_cont, z = z)
  dat <- data.frame(y = y, trt = trt_factor, x = x_cont, z = z)

  expect_error(
    mfp2(
      x_matrix[, c("x", "z"), drop = FALSE],
      y,
      subset = duplicate_rows,
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    duplicate_message
  )

  expect_error(
    mfp2(
      y ~ x + z,
      data = dat,
      subset = duplicate_rows,
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    ),
    duplicate_message
  )

  expect_error(
    mfpi(
      x_matrix,
      y,
      subset = duplicate_rows,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    duplicate_message
  )

  expect_error(
    mfpi(
      y ~ trt + x + z,
      data = dat,
      subset = duplicate_rows,
      group_var = "trt",
      cont_vars = "x",
      cont_var_forms = c(x = "linear"),
      df = 1,
      shift = 0,
      scale = 1,
      center = FALSE,
      p_interact = 1,
      verbose = FALSE
    ),
    duplicate_message
  )
})


# Test purpose: Unique numeric positions are not sorted. Their supplied order is
# retained by the row resolver and by all four fitting interfaces.
test_that("24.14 unique numeric subset indices preserve supplied order", {
  set.seed(2414)
  n <- 100L
  x <- cbind(
    x = seq_len(n),
    z = rnorm(n)
  )
  y <- 1 + 0.02 * x[, "x"] + rnorm(n, sd = 0.2)
  dat <- data.frame(y = y, x = x[, "x"], z = x[, "z"])
  rows <- c(51:100, 1:50)

  expect_identical(formula_subset_rows(rows, n), as.integer(rows))
  expect_error(
    formula_subset_rows(c(1L, 2L, 2L), n),
    "must not contain duplicated row indices"
  )
  expect_error(
    subset_formula_model_frame(
      data.frame(x = seq_len(4L)),
      c(1L, 2L, 2L)
    ),
    "unique valid integer model-frame positions"
  )

  fit_matrix <- mfp2(
    x,
    y,
    subset = rows,
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  expect_equal(
    unname(fit_matrix$x_original[, "x"]),
    unname(x[rows, "x"]),
    tolerance = 0
  )

  fit_formula <- mfp2(
    y ~ x + z,
    data = dat,
    subset = rows,
    keep = c("x", "z"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    xorder = "original",
    verbose = FALSE
  )
  expect_equal(
    unname(fit_formula$x_original[, "x"]),
    unname(dat$x[rows]),
    tolerance = 0
  )

  trt_factor <- factor(rep(c("control", "treated"), length.out = n))
  trt_numeric <- as.numeric(trt_factor == "treated")
  x_mfpi <- cbind(trt = trt_numeric, x = x[, "x"], z = x[, "z"])
  dat_mfpi <- data.frame(y = y, trt = trt_factor, x = x[, "x"], z = x[, "z"])

  fit_mfpi_matrix <- mfpi(
    x_mfpi,
    y,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_equal(
    unname(fit_mfpi_matrix$x_train_internal[, "x"]),
    unname(x_mfpi[rows, "x"]),
    tolerance = 0
  )

  fit_mfpi_formula <- mfpi(
    y ~ trt + x + z,
    data = dat_mfpi,
    subset = rows,
    group_var = "trt",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    df = 1,
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 1,
    verbose = FALSE
  )
  expect_equal(
    unname(fit_mfpi_formula$x_train_internal[, "x"]),
    unname(dat_mfpi$x[rows]),
    tolerance = 0
  )
})
