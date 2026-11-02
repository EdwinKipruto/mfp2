# Area 2: FP selection correctness (df = 4, criterion = "aic") for each family.

fp_selection_checks <- function(fit, family_string) {
  expect_equal(fit$family_string, family_string)
  expect_true(isTRUE(fit$convergence_mfp))
  # selected FP powers are recorded
  expect_false(is.null(fit$fp_powers))
  # final AIC is no worse than the forced-linear model's AIC (checked by caller)
}

test_that("survreg: AIC-driven FP selection is recorded and not worse than linear", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 4) + fp(nodes, df = 4) + hormon,
              data = gbsg, family = survreg_family(dist = "weibull"),
              criterion = "aic", verbose = FALSE)
  fp_selection_checks(fit, "survreg")

  lin <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1) + fp(nodes, df = 1) + hormon,
              data = gbsg, family = survreg_family(dist = "weibull"),
              criterion = "aic", verbose = FALSE)
  expect_lte(AIC(fit), AIC(lin) + 1e-6)
})

test_that("finegray: AIC-driven FP selection is recorded and not worse than linear", {
  skip_on_cran()
  dat <- make_finegray_data(n = 700, seed = 31)
  fit <- mfp2(Surv(obs, ev) ~ fp(z1, df = 4),
              data = dat, family = finegray_family(etype = "relapse"),
              criterion = "aic", verbose = FALSE)
  fp_selection_checks(fit, "finegray")

  lin <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1),
              data = dat, family = finegray_family(etype = "relapse"),
              criterion = "aic", verbose = FALSE)
  # AIC via -2logLik + 2*df (coxph AIC is well defined for the weighted fit)
  aic <- function(f) -2 * f$loglik[2] + 2 * length(coef(f))
  expect_lte(aic(fit), aic(lin) + 1e-6)
})

test_that("multinomial: AIC-driven FP selection works with default settings", {
  # Per the prompt: df = 4, criterion = 'aic', >= 1 continuous predictor.
  # This uses mfp2()'s default fitting controls, as an ordinary user would.
  skip_on_cran()
  dat <- make_multinomial_data(n = 500, seed = 33)
  fit <- mfp2(y ~ fp(x1, df = 4) + fp(x2, df = 4),
              data = dat, family = multinomial_family(),
              criterion = "aic", verbose = FALSE)
  fp_selection_checks(fit, "multinomial")

  lin <- mfp2(y ~ fp(x1, df = 1) + fp(x2, df = 1),
              data = dat, family = multinomial_family(),
              criterion = "aic", verbose = FALSE)
  expect_lte(AIC(fit), AIC(lin) + 1e-6)
})

test_that("multinomial FP selection succeeds once maxit is raised (root-cause probe)", {
  # Companion to the test above: the multinomial candidate fast path takes its
  # iteration cap from glm.control() (maxit = 25 by default), which is too low
  # for nnet's optimiser on FP-transformed designs. Raising maxit fixes it.
  # If BOTH this test passes and the default-settings test above fails, the
  # defect is the default maxit, not the FP machinery itself.
  skip_on_cran()
  dat <- make_multinomial_data(n = 500, seed = 33)
  expect_error(
    mfp2(y ~ fp(x1, df = 4) + fp(x2, df = 4),
         data = dat, family = multinomial_family(),
         criterion = "aic", control = list(maxit = 300), verbose = FALSE),
    NA
  )
})
