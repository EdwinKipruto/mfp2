# GEE fitting correctness vs a direct geepack::geeglm() fit.
#
# The binding contract: with df = 1, select = 1, center = FALSE the retained
# mfp2 GEE model is a linear geeglm fit. Coefficients, fitted values, and
# working-correlation parameters reproduce geepack; the sandwich covariance
# additionally carries Stata's K/(K-1) correction. QICu must equal
# geepack::QIC().

make_gee_data <- function(K = 60, m = 4, seed = 11) {
  set.seed(seed)
  n <- K * m
  id <- rep(seq_len(K), each = m)
  wave <- rep(seq_len(m), times = K)
  x1 <- stats::rnorm(n)
  x2 <- stats::runif(n, -1, 1)
  b <- stats::rnorm(K)[id]
  eta <- 1 + 0.8 * x1 - 0.5 * x2 + b
  data.frame(id, wave, x1, x2, eta_lin = eta)
}

fit_both <- function(dat, gfam, corstr, waves_col = NULL) {
  waves_vec <- if (is.null(waves_col)) NULL else dat[[waves_col]]
  m2 <- mfp2(y ~ x1 + x2, data = dat,
             family = gee_family(gfam, corstr = corstr),
             id = dat$id, waves = waves_vec,
             df = 1, select = 1, center = FALSE, verbose = FALSE)
  if (is.null(waves_col)) {
    gg <- geepack::geeglm(y ~ x1 + x2, data = dat, family = gfam,
                          id = id, corstr = corstr)
  } else {
    dat$.wv <- waves_vec
    gg <- geepack::geeglm(y ~ x1 + x2, data = dat, family = gfam,
                          id = id, waves = .wv, corstr = corstr)
  }
  list(m2 = m2, gg = gg)
}

expect_gee_equivalent <- function(m2, gg) {
  expect_true(inherits(m2, "geeglm"))
  expect_equal(unname(coef(m2)), unname(coef(gg)), tolerance = TOL_MED)
  k <- length(gg$geese$clusz)
  expected_vcov <- as.matrix(vcov(gg)) * k / (k - 1)
  expect_equal(unname(as.matrix(vcov(m2))), unname(expected_vcov),
               tolerance = TOL_MED)
  expect_equal(unname(fitted(m2)), unname(fitted(gg)), tolerance = TOL_MED)
  expect_equal(unname(m2$geese$alpha), unname(gg$geese$alpha), tolerance = TOL_MED)
  expect_equal(unname(m2$geese$gamma), unname(gg$geese$gamma), tolerance = TOL_MED)
  expect_false(isTRUE(m2$geese$model$scale.fix))
}

test_that("Gaussian GEE matches geeglm apart from the documented VCE correction", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data()
  dat$y <- dat$eta_lin + stats::rnorm(nrow(dat), sd = 0.5)

  for (cs in c("independence", "exchangeable")) {
    fb <- fit_both(dat, gaussian(), cs)
    expect_gee_equivalent(fb$m2, fb$gg)
  }
  fb <- fit_both(dat, gaussian(), "ar1", waves_col = "wave")
  expect_gee_equivalent(fb$m2, fb$gg)
})

test_that("retained AR1 GEE preserves gaps in factor-valued waves", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data(K = 30, m = 2)
  dat$wave_gap <- factor(rep(c("1", "3"), times = 30),
                         levels = c("1", "3"))
  dat$y <- dat$eta_lin + stats::rnorm(nrow(dat), sd = 0.5)

  m2 <- mfp2(
    y ~ x1 + x2,
    data = dat,
    family = gee_family(gaussian(), corstr = "ar1"),
    id = dat$id,
    waves = dat$wave_gap,
    df = 1,
    select = 1,
    center = FALSE,
    verbose = FALSE
  )

  X <- stats::model.matrix(~ x1 + x2, dat)
  direct <- geepack::geese.fit(
    x = X,
    y = dat$y,
    id = dat$id,
    waves = rep(c(1L, 3L), times = 30),
    family = gaussian(),
    corstr = "ar1"
  )

  expect_identical(m2$geese$waves, rep(c(1L, 3L), times = 30))
  expect_equal(unname(coef(m2)), unname(direct$beta), tolerance = TOL_MED)
  expect_equal(unname(m2$geese$alpha), unname(direct$alpha),
               tolerance = TOL_MED)
})

test_that("Binomial GEE matches geeglm", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data()
  dat$y <- stats::rbinom(nrow(dat), 1, plogis(dat$eta_lin - 1))
  fb <- fit_both(dat, binomial(), "exchangeable")
  expect_gee_equivalent(fb$m2, fb$gg)
  fb <- fit_both(dat, binomial(), "ar1", waves_col = "wave")
  expect_gee_equivalent(fb$m2, fb$gg)
})

test_that("Poisson GEE matches geeglm", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data()
  dat$y <- stats::rpois(nrow(dat), exp(0.4 + 0.3 * dat$x1 - 0.2 * dat$x2))
  fb <- fit_both(dat, poisson(), "exchangeable")
  expect_gee_equivalent(fb$m2, fb$gg)
})

test_that("Gamma GEE matches geeglm", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data()
  mu <- exp(0.5 + 0.2 * dat$x1)
  dat$y <- stats::rgamma(nrow(dat), shape = 2, rate = 2 / mu)
  fb <- fit_both(dat, Gamma(link = "log"), "exchangeable")
  expect_gee_equivalent(fb$m2, fb$gg)
})

test_that("Offset is handled equivalently", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data()
  off <- log(runif(nrow(dat), 1, 3))
  dat$y <- stats::rpois(nrow(dat), exp(0.2 + 0.3 * dat$x1 + off))
  dat$off <- off
  m2 <- mfp2(y ~ x1 + x2, data = dat, family = gee_family(poisson(), corstr = "exchangeable"),
             id = dat$id, offset = dat$off, df = 1, select = 1, center = FALSE, verbose = FALSE)
  gg <- geepack::geeglm(y ~ x1 + x2 + offset(off), data = dat, family = poisson(),
                        id = id, corstr = "exchangeable")
  expect_equal(unname(coef(m2)), unname(coef(gg)), tolerance = TOL_MED)
  expect_equal(unname(fitted(m2)), unname(fitted(gg)), tolerance = TOL_MED)
})

test_that("Reported QICu equals geepack::QIC() for the retained model", {
  skip_if_not_installed("geepack")
  dat <- make_gee_data()
  dat$y <- stats::rpois(nrow(dat), exp(0.4 + 0.3 * dat$x1 - 0.2 * dat$x2))

  m2 <- mfp2(y ~ x1 + x2, data = dat,
             family = gee_family(poisson(), corstr = "exchangeable"),
             id = dat$id, df = 1, select = 1, center = FALSE, verbose = FALSE)

  # Build the reference geeglm here in the test frame with a literal family and
  # `data = dat`. geepack::QIC() refits an independence model by re-evaluating
  # the fitted object's stored `call` in the calling frame, so every symbol in
  # that call (data, family, id, corstr) must resolve where QIC() is invoked.
  gg <- geepack::geeglm(y ~ x1 + x2, data = dat, family = poisson(),
                        id = id, corstr = "exchangeable")
  q <- geepack::QIC(gg)

  # mfp2's retained fit reproduces the direct geeglm fit ...
  expect_equal(unname(fitted(m2)), unname(fitted(gg)), tolerance = TOL_MED)

  # ... and the QICu mfp2 uses for selection matches geepack's QICu exactly.
  mu <- fitted(gg)
  quasi <- sum(dat$y * log(mu) - mu)           # geepack Poisson quasi-likelihood
  p <- length(coef(gg))
  expect_equal(-2 * quasi + 2 * p, unname(q["QICu"]), tolerance = TOL_MED)
})
