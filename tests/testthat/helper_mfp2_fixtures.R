# Shared package fixtures used by multiple focused test files.

library(testthat)
library(survival)
# Keep survival formula specials unqualified. Some older survival releases
# used by package-check services do not recognize the namespace-qualified
# form as a formula special, whereas strata() is supported consistently.
library(mfp2)


# =============================================================================
# Test data setup
# =============================================================================

data("prostate", package = "mfp2")

# Shared prostate fixtures used by many Gaussian/default-interface tests.
# x_prostate contains only predictors; y_prostate is the continuous response.
x_prostate <- as.matrix(prostate[, 2:8])
y_prostate <- as.numeric(prostate$lpsa)

# Helper: suppress verbose progress messages in tests that intentionally call
# functions with user-facing output. Keep test assertions outside quiet().
quiet <- function(expr) suppressMessages(capture.output(expr, type = "message"))


# Compare mfp2 and glm parameters positionally. Every GLM-equivalence fit in
# this suite uses xorder = "original", so the coefficient and covariance order
# must match the reference glm() fit directly. mfp2's transformed-column suffix
# is deliberately ignored because these assertions compare numerical values.
expect_mfp2_glm_parameters_equal <- function(fit_mfp2,
                                             fit_glm,
                                             tolerance = 1e-8) {
  coef_mfp2 <- stats::coef(fit_mfp2)
  coef_glm <- stats::coef(fit_glm)
  vcov_mfp2 <- stats::vcov(fit_mfp2)
  vcov_glm <- stats::vcov(fit_glm)

  expect_length(coef_mfp2, length(coef_glm))
  expect_equal(dim(vcov_mfp2), dim(vcov_glm))
  expect_equal(unname(coef_mfp2), unname(coef_glm), tolerance = tolerance)
  expect_equal(unname(vcov_mfp2), unname(vcov_glm), tolerance = tolerance)
  expect_equal(
    unname(sqrt(diag(vcov_mfp2))),
    unname(sqrt(diag(vcov_glm))),
    tolerance = tolerance
  )
}

# Generate stable negative-binomial data for public-interface tests. Keeping the
# signal moderate avoids boundary estimates of theta and makes comparisons with
# an independent MASS::glm.nb() fit reproducible.
make_negbin_test_data <- function(n = 240L, seed = 104L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::runif(n, 0.5, 2.5),
    x2 = stats::rnorm(n),
    x3 = stats::rbinom(n, 1L, 0.4)
  )
  eta <- 0.25 + 0.35 * dat$x1 - 0.25 * dat$x2 + 0.30 * dat$x3
  dat$y <- stats::rnbinom(n, mu = exp(eta), size = 3)
  dat
}

# Cache one independent linear negative-binomial reference problem. Several
# public-method tests use the same fit so that comprehensive NB coverage does
# not repeatedly run the full MFP algorithm.
.negbin_reference_cache <- new.env(parent = emptyenv())

get_negbin_reference_fits <- function() {
  if (!exists("fits", envir = .negbin_reference_cache, inherits = FALSE)) {
    dat <- make_negbin_test_data(n = 320L, seed = 1103L)

    fit_mfp2 <- mfp2(
      y ~ x1 + x2 + x3,
      data = dat,
      family = "negbin",
      fitter = "fastglm",
      cycles = 1,
      df = 1,
      select = 1,
      alpha = 1,
      xorder = "original",
      shift = 0,
      scale = 1,
      center = FALSE,
      verbose = FALSE
    )

    fit_mass <- MASS::glm.nb(
      y ~ x1 + x2 + x3,
      data = dat,
      link = log,
      control = stats::glm.control(maxit = 100L)
    )

    assign(
      "fits",
      list(
        data = dat,
        mfp2 = fit_mfp2,
        mass = fit_mass
      ),
      envir = .negbin_reference_cache
    )
  }

  get("fits", envir = .negbin_reference_cache, inherits = FALSE)
}

expect_negbin_mfp2_mass_equal <- function(fit_mfp2,
                                          fit_mass,
                                          coefficient_tolerance = 1e-4,
                                          likelihood_tolerance = 1e-3) {
  expect_equal(
    unname(stats::coef(fit_mfp2)),
    unname(stats::coef(fit_mass)),
    tolerance = coefficient_tolerance
  )
  expect_equal(
    unname(stats::vcov(fit_mfp2)),
    unname(stats::vcov(fit_mass)),
    tolerance = coefficient_tolerance
  )
  expect_equal(
    unname(fit_mfp2$fitted.values),
    unname(stats::fitted(fit_mass)),
    tolerance = coefficient_tolerance
  )
  expect_equal(
    unname(fit_mfp2$deviance),
    unname(stats::deviance(fit_mass)),
    tolerance = likelihood_tolerance
  )
  expect_equal(
    fit_mfp2$theta,
    fit_mass$theta,
    tolerance = likelihood_tolerance
  )
  expect_equal(
    fit_mfp2$mfp_logl,
    as.numeric(stats::logLik(fit_mass)),
    tolerance = likelihood_tolerance
  )
  expect_equal(
    fit_mfp2$aic,
    stats::AIC(fit_mass),
    tolerance = likelihood_tolerance
  )
}
