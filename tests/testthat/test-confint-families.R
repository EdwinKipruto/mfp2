# confint.mfp2() returns Wald intervals from vcov() for every scalar-coefficient
# family, applying uniformly rather than relying on the profile-likelihood path
# that fails on mfp2's internally transformed GLM design.
#
# Models are fitted with df = 1, select = 1 so every predictor is retained as a
# linear term: this keeps the fits well-conditioned (no extreme FP-power
# candidates, no dropped-to-empty models) and makes the coefficient set that
# confint() reports deterministic.

expect_wald_confint <- function(m, level = 0.95) {
  ci <- confint(m, level = level)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2L)
  expect_gt(nrow(ci), 0L)

  cf <- coef(m)
  V <- as.matrix(vcov(m))
  nm <- intersect(names(cf), rownames(V))
  expect_gt(length(nm), 0L)
  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- sqrt(diag(V))[nm]
  expect_equal(unname(ci[nm, 1]), unname(cf[nm] - z * se), tolerance = 1e-8)
  expect_equal(unname(ci[nm, 2]), unname(cf[nm] + z * se), tolerance = 1e-8)
}

lin <- function(v) sprintf("fp(%s, df = 1, select = 1)", v)

test_that("confint() works for GLM families (gaussian/binomial/poisson/Gamma)", {
  set.seed(101)
  n <- 400
  x1 <- rnorm(n); x2 <- runif(n, 1, 5)

  dg <- data.frame(y = 1 + 0.7 * x1 + 0.5 * x2 + rnorm(n), x1, x2)
  expect_wald_confint(mfp2(y ~ fp(x1, df = 1, select = 1) + fp(x2, df = 1, select = 1),
                           data = dg, verbose = FALSE))

  # moderate signal + large n avoids perfect separation
  db <- data.frame(y = rbinom(n, 1, plogis(-0.2 + 0.6 * x1)), x1, x2)
  expect_wald_confint(mfp2(y ~ fp(x1, df = 1, select = 1) + fp(x2, df = 1, select = 1),
                           data = db, family = "binomial", verbose = FALSE))

  dp <- data.frame(y = rpois(n, exp(0.4 + 0.2 * x1 + 0.1 * x2)), x1, x2)
  expect_wald_confint(mfp2(y ~ fp(x1, df = 1, select = 1) + fp(x2, df = 1, select = 1),
                           data = dp, family = "poisson", verbose = FALSE))

  dga <- data.frame(y = rgamma(n, shape = 4, rate = 4 / exp(0.5 + 0.2 * x1)), x1, x2)
  expect_wald_confint(mfp2(y ~ fp(x1, df = 1, select = 1) + fp(x2, df = 1, select = 1),
                           data = dga, family = "Gamma", verbose = FALSE))
})

test_that("confint() works for Cox and parametric survival models", {
  skip_if_not_installed("survival")
  set.seed(102)
  n <- 400
  x1 <- rnorm(n); x2 <- runif(n, 1, 5)
  tt <- rexp(n, exp(0.5 * x1 + 0.2 * x2)); ev <- rbinom(n, 1, 0.8)
  dc <- data.frame(tt, ev, x1, x2)

  expect_wald_confint(mfp2(survival::Surv(tt, ev) ~ fp(x1, df = 1, select = 1) +
                             fp(x2, df = 1, select = 1),
                           data = dc, family = "cox", verbose = FALSE))
  expect_wald_confint(mfp2(survival::Surv(tt, ev) ~ fp(x1, df = 1, select = 1) +
                             fp(x2, df = 1, select = 1),
                           data = dc, family = survreg_family("weibull"),
                           verbose = FALSE))
})

test_that("confint() honours level and parm", {
  set.seed(103)
  n <- 400
  x1 <- rnorm(n); x2 <- runif(n, 1, 5)
  d <- data.frame(y = 1 + 0.6 * x1 - 0.4 * x2 + rnorm(n), x1, x2)
  m <- mfp2(y ~ fp(x1, df = 1, select = 1) + fp(x2, df = 1, select = 1),
            data = d, verbose = FALSE)

  ci95 <- confint(m, level = 0.95)
  ci90 <- confint(m, level = 0.90)
  expect_true(all(ci90[, 1] > ci95[, 1]))            # 90% narrower
  one <- confint(m, parm = rownames(ci95)[2])
  expect_equal(nrow(one), 1L)
})

test_that("confint() returns an empty matrix when no covariates are retained", {
  skip_if_not_installed("survival")
  set.seed(104)
  n <- 300
  # A Cox model has no intercept, so strict selection on pure-noise predictors
  # can drop every covariate and leave a genuinely empty coefficient vector.
  # confint() must return an empty matrix rather than erroring on the empty fit.
  tt <- rexp(n); ev <- rbinom(n, 1, 0.7)
  d <- data.frame(tt, ev, x1 = rnorm(n), x2 = rnorm(n))
  m <- mfp2(survival::Surv(tt, ev) ~ fp(x1, df = 1) + fp(x2, df = 1),
            data = d, family = "cox", select = 0.001, verbose = FALSE)
  ci <- confint(m)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), length(coef(m)))          # 0 when all dropped
  expect_equal(ncol(ci), 2L)
})
