# Area 1c: survreg_family() fitting correctness vs survival::survreg()

test_that("mfp2 survreg (weibull) matches survival::survreg on gbsg", {
  skip_on_cran()
  data(gbsg, package = "mfp2")

  fit <- mfp2(Surv(rectime, censrec) ~
                fp(age,   df = 1, select = 1, scale = 1, shift = 0) +
                fp(nodes, df = 1, select = 1, scale = 1, shift = 0) + hormon,
              data = gbsg, family = survreg_family(dist = "weibull"),
              center = FALSE, verbose = FALSE)
  ref <- survival::survreg(Surv(rectime, censrec) ~ age + nodes + hormon,
                           data = gbsg, dist = "weibull")

  expect_equal(fit$family_string, "survreg")
  expect_true(inherits(fit, "survreg"))

  # slopes match exactly; the intercept is on a different (centred) convention,
  # so it is compared through predictions below rather than directly.
  cf <- coef(fit)
  names(cf) <- sub("\\.1$", "", names(cf))
  for (nm in c("age", "nodes", "hormon"))
    expect_equal(unname(cf[[nm]]), unname(coef(ref)[[nm]]), tolerance = TOL_MED)

  expect_equal(fit$scale, ref$scale, tolerance = TOL_MED)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = TOL_MED)

  # predictions are invariant to the intercept convention -> machine precision
  nd <- gbsg[1:20, ]
  expect_equal(max_abs_diff(predict(fit, newdata = nd, type = "link"),
                            predict(ref, newdata = nd, type = "lp")), 0,
               tolerance = TOL_TIGHT)
  expect_equal(max_abs_diff(predict(fit, newdata = nd, type = "response"),
                            predict(ref, newdata = nd, type = "response")), 0,
               tolerance = TOL_MED)
})

test_that("survreg quantile predictions match predict.survreg", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) +
                fp(nodes, df = 1, select = 1),
              data = gbsg, family = survreg_family(dist = "weibull"),
              verbose = FALSE)
  ref <- survival::survreg(Surv(rectime, censrec) ~ age + nodes,
                           data = gbsg, dist = "weibull")
  nd <- gbsg[1:15, ]
  p <- c(0.1, 0.5, 0.9)
  q_fit <- predict(fit, newdata = nd, type = "quantile", p = p)
  q_ref <- predict(ref, newdata = nd, type = "quantile", p = p)
  expect_equal(max_abs_diff(q_fit, q_ref), 0, tolerance = TOL_MED)
})

test_that("additional survreg distributions fit and match native survreg", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  for (d in c("exponential", "lognormal", "loglogistic")) {
    fit <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) +
                  fp(nodes, df = 1, select = 1),
                data = gbsg, family = survreg_family(dist = d), verbose = FALSE)
    ref <- survival::survreg(Surv(rectime, censrec) ~ age + nodes,
                             data = gbsg, dist = d)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
                 tolerance = TOL_MED, info = d)
    nd <- gbsg[1:10, ]
    expect_equal(max_abs_diff(predict(fit, newdata = nd, type = "link"),
                              predict(ref, newdata = nd, type = "lp")), 0,
                 tolerance = TOL_TIGHT, info = d)
  }
})

test_that("survreg with strata estimates per-stratum scales matching survreg", {
  skip_on_cran()
  dat <- make_survreg_data(n = 600, seed = 5)

  fit <- mfp2(Surv(time, status) ~ fp(x1, df = 1, select = 1) + strata(sx),
              data = dat, family = survreg_family(dist = "weibull"),
              verbose = FALSE)
  ref <- survival::survreg(Surv(time, status) ~ x1 + strata(sx),
                           data = dat, dist = "weibull")

  # two scales, one per stratum, matching native survreg
  expect_length(fit$scale, 2L)
  expect_equal(sort(unname(fit$scale)), sort(unname(ref$scale)), tolerance = TOL_MED)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = TOL_MED)
  expect_equal(fit$mfp2_strata_levels, levels(dat$sx))
})

test_that("fixed scale combined with multiple strata errors", {
  skip_on_cran()
  dat <- make_survreg_data(n = 400, seed = 6)
  expect_error(
    mfp2(Surv(time, status) ~ fp(x1, df = 1, select = 1) + strata(sx),
         data = dat, family = survreg_family(dist = "weibull", scale = 0.8),
         verbose = FALSE),
    regexp = "fixed .*scale|scale strata", ignore.case = TRUE
  )
})

test_that("left-censored survreg response matches native survreg", {
  skip_on_cran()
  set.seed(123); n <- 500
  x1 <- runif(n, 1, 10)
  logt <- 3 + 0.15 * x1 + 0.8 * rnorm(n)
  tt <- exp(logt)
  # per-observation lower detection limit -> proper (well-conditioned) left
  # censoring; status 2 == left censored in Surv(type = "left")
  L <- exp(quantile(logt, 0.15) + 0.3 * rnorm(n))
  event <- ifelse(tt < L, 2L, 1L)
  time  <- ifelse(tt < L, L, tt)
  dat <- data.frame(time = time, event = event, x1 = x1)
  yy <- survival::Surv(dat$time, dat$event, type = "left")

  ref <- survival::survreg(yy ~ x1, data = dat, dist = "lognormal")
  fit <- mfp2(yy ~ fp(x1, df = 1, select = 1), data = dat,
              family = survreg_family(dist = "lognormal"), verbose = FALSE)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = TOL_MED)
})
