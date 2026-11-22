# Regression coverage for spline basis construction, reuse, and prediction.

test_that("small-sample RCS knots match the documented Hmisc rule", {
  old <- options(fractied = 0.05)
  on.exit(options(old), add = TRUE)
  # Hmisc::rcspline.eval(1:8, nk = 3, knots.only = TRUE): initial
  # quantiles 1.7, 4.5, 7.3; replace ends by fifth order statistics,
  # then sort, giving 4, 4.5, 5.
  expect_equal(mfp2:::rcs_default_knots(1:8, 2L), c(4, 4.5, 5))
})

test_that("automatic RCS knots agree with Hmisc for six to eight rows", {
  skip_if_not_installed("Hmisc")
  old <- options(fractied = 0.05)
  on.exit(options(old), add = TRUE)
  for (n in 6:8) {
    expect_equal(mfp2:::rcs_default_knots(seq_len(n), 2L),
      as.numeric(Hmisc::rcspline.eval(seq_len(n), nk = 3, knots.only = TRUE)))
  }
})

test_that("B-spline evaluation computes each distinct knot set once", {
  x <- c(NA_real_, seq(-0.2, 1.2, length.out = 101))
  knots <- c(0.25, 0.5, 0.75)
  boundaries <- c(0, 1)
  terms <- mfp2:::bspline_basis_terms(knots, boundaries)
  # Match the previous per-column evaluator, including NA and extrapolation.
  reference <- cbind(x, do.call(cbind, lapply(terms, function(term) {
    splines::bs(x, knots = term$knots, Boundary.knots = boundaries,
                 intercept = FALSE, warn.outside = FALSE)[, term$column]
  })))
  calls <- 0L
  original_bs <- splines::bs
  testthat::local_mocked_bindings(bs = function(...) {
    calls <<- calls + 1L
    original_bs(...)
  }, .package = "splines")
  actual <- mfp2:::bspline_basis(x, knots, boundaries, terms = terms)
  expect_identical(calls, 2L)
  expect_equal(as.numeric(actual), as.numeric(reference), tolerance = 1e-12)
  # The grouping must also preserve nonconsecutive and reordered descriptors.
  order <- rev(seq_along(terms))
  reordered <- mfp2:::bspline_basis(x, knots, boundaries, terms = terms[order])
  expect_equal(as.numeric(reordered),
               as.numeric(actual[, c(1L, order + 1L)]), tolerance = 1e-12)
})

test_that("spline selection prepares adjustments and fits the null once per step", {
  original_step <- mfp2:::select_spline_step
  original_null <- mfp2:::fit_null_step
  original_adjustment <- mfp2:::build_adjustment_step
  original_linear <- mfp2:::fit_linear_step
  active <- FALSE
  null_calls <- adjustment_calls <- linear_calls <- 0L
  steps <- 0L
  testthat::local_mocked_bindings(
    select_spline_step = function(...) {
      active <<- TRUE
      on.exit({ active <<- FALSE }, add = TRUE)
      null_calls <<- adjustment_calls <<- linear_calls <<- 0L
      result <- original_step(...)
      expect_identical(null_calls, 1L)
      expect_identical(adjustment_calls, 1L)
      expect_identical(linear_calls, 7L)
      steps <<- steps + 1L
      result
    },
    fit_null_step = function(...) {
      if (active) null_calls <<- null_calls + 1L
      original_null(...)
    },
    build_adjustment_step = function(...) {
      if (active) adjustment_calls <<- adjustment_calls + 1L
      original_adjustment(...)
    },
    fit_linear_step = function(...) {
      if (active) linear_calls <<- linear_calls + 1L
      original_linear(...)
    }, .package = "mfp2"
  )
  x <- c(rep(0, 30), seq(0.1, 3, length.out = 100))
  dat <- data.frame(x = x, y = 2 * (x == 0) + sin(3 * x) +
                     0.1 * cos(seq_along(x)))
  for (criterion in c("pvalue", "aic", "bic")) {
    for (form in list(y ~ rcs(x, df = 3, spike = TRUE),
                      y ~ bspline(x, df = 6, spike = TRUE))) {
      fit <- mfp2(form, data = dat, criterion = criterion,
                  ftest = criterion == "pvalue", verbose = FALSE)
      expect_s3_class(fit, "mfp2")
      expect_true(all(is.finite(predict(fit))))
    }
  }
  expect_gt(steps, 0L)
})

test_that("ordinary formula and matrix fits never construct spline bases", {
  forbidden <- function(...) stop("unexpected spline construction")
  testthat::local_mocked_bindings(
    rcs_basis = forbidden, bspline_basis = forbidden,
    rcs_default_knots = forbidden, bspline_default_knots = forbidden,
    prepare_matrix_splines = forbidden, .package = "mfp2"
  )
  x <- seq(1, 3, length.out = 90)
  dat <- data.frame(x = x, y = log(x) + 0.01 * sin(seq_along(x)))
  f <- mfp2(y ~ fp(x, df = 2), data = dat, verbose = FALSE)
  m <- mfp2(cbind(x = x), dat$y, df = 2, verbose = FALSE)
  expect_s3_class(f, "mfp2")
  expect_s3_class(m, "mfp2")
  expect_equal(as.numeric(predict(f)), as.numeric(predict(m)), tolerance = 1e-7)
})

test_that("RCS with fixed knots matches Hmisc on missing rows and tails", {
  skip_if_not_installed("Hmisc")
  knots <- c(-1, -0.5, 0, 0.5, 1)
  x <- c(-5, -2, -1, 0, 1, 2, 5, NA_real_)
  own <- mfp2:::rcs_basis(x, knots)
  ref <- Hmisc::rcspline.eval(x, knots = knots, inclx = TRUE, norm = 2)
  expect_equal(as.numeric(own), as.numeric(ref), tolerance = 1e-10)
})

test_that("B-spline predictions reuse training metadata beyond boundaries", {
  x <- seq(0, 1, length.out = 150)
  dat <- data.frame(x = x, y = sin(8 * x) + 0.02 * cos(31 * x))
  knots <- c(0.25, 0.5, 0.75)
  fit <- mfp2(y ~ bspline(x, knots = knots, boundary_knots = c(0, 1)),
              data = dat, select = 1, alpha = 1, verbose = FALSE)
  ref <- stats::lm(y ~ splines::bs(x, knots = knots,
                       Boundary.knots = c(0, 1)), data = dat)
  nd <- data.frame(x = c(-1, -0.1, 0, 0.3, 1, 1.1, 2))
  expected <- suppressWarnings(predict(ref, nd))
  actual <- predict(fit, nd)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-6)
  expect_equal(as.numeric(predict(fit, nd[1, , drop = FALSE])),
               as.numeric(actual)[1], tolerance = 1e-7)
})
