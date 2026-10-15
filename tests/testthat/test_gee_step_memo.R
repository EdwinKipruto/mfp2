# When a later MFP cycle revisits a GEE term with an unchanged adjustment
# model, the previous (deterministic) power search is reused. These tests check
# that the reuse happens and that it never changes the fitted result.

core_gee_data <- function(K = 30, m = 4, seed = 101) {
  set.seed(seed)
  n <- K * m
  id <- rep(seq_len(K), each = m)
  wave <- rep(seq_len(m), K)
  x1 <- stats::runif(n, 1, 5)
  x2 <- stats::runif(n, 0.5, 3)
  b <- stats::rnorm(K, sd = 0.4)[id]
  tot <- sample(3:7, n, replace = TRUE)
  data.frame(
    id = id, wave = wave, x1 = x1, x2 = x2,
    expo = stats::runif(n, 0.5, 2), w = stats::runif(n, 0.5, 2),
    yg = 1 + 2 / x1 + 0.4 * x2 + b + stats::rnorm(n, sd = 0.5),
    yp = stats::rpois(n, exp(0.2 + 0.5 / x1 + 0.2 * x2 + b)),
    ya = stats::rgamma(n, shape = 5, rate = 5 / exp(0.2 + 0.3 * x2)),
    succ = stats::rbinom(n, tot, stats::plogis(-0.3 + 1 / x1 + b)),
    tot = tot
  )
}

test_that("GEE MFP selection is unchanged by the step memo", {
  skip_if_not_installed("geepack")
  d <- core_gee_data(K = 40, seed = 23)
  fit <- function(criterion) {
    suppressWarnings(mfp2(
      yg ~ fp(x1) + fp(x2), data = d,
      family = gee_family(stats::gaussian(), corstr = "exchangeable"),
      id = id, criterion = criterion, verbose = FALSE
    ))
  }
  for (criterion in c("pvalue", "bic")) {
    fast <- fit(criterion)
    old <- options(mfp2.gee_step_memo = FALSE)
    slow <- fit(criterion)
    options(old)
    expect_identical(fast$fp_powers, slow$fp_powers)
    expect_identical(coef(fast), coef(slow))
    expect_identical(vcov(fast), vcov(slow))
    expect_identical(fast$fp_terms, slow$fp_terms)
    expect_false(is.null(fast$mfp_selection_score))
    expect_identical(fast$mfp_selection_score, slow$mfp_selection_score)
  }
})

test_that("a revisited GEE term with unchanged adjustments reuses its search", {
  skip_if_not_installed("geepack")
  d <- core_gee_data(K = 40, seed = 23)
  calls <- 0L
  real_select <- select_ra2
  testthat::local_mocked_bindings(
    select_ra2 = function(...) {
      calls <<- calls + 1L
      real_select(...)
    },
    .package = "mfp2"
  )
  fit <- function() {
    calls <<- 0L
    m <- suppressWarnings(mfp2(
      yg ~ fp(x1) + fp(x2), data = d,
      family = gee_family(stats::gaussian(), corstr = "exchangeable"),
      id = id, criterion = "pvalue", verbose = FALSE
    ))
    list(model = m, calls = calls)
  }
  memo <- fit()
  old <- options(mfp2.gee_step_memo = FALSE)
  plain <- fit()
  options(old)

  # Two terms, at least two cycles: the final cycle revisits both terms.
  expect_gte(plain$calls, 4L)
  expect_lt(memo$calls, plain$calls)
  expect_identical(memo$model$fp_powers, plain$model$fp_powers)
  expect_identical(coef(memo$model), coef(plain$model))
})
