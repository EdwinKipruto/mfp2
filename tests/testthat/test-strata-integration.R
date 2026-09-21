# Area 7: cross-cutting strata integration.

find_fit2 <- function(x, cls) {
  if (inherits(x, cls)) return(x)
  if (is.list(x)) for (el in x) { r <- find_fit2(el, cls); if (!is.null(r)) return(r) }
  NULL
}
n_strata <- function(fit) if (is.null(fit$strata)) 0L else length(unique(fit$strata))

test_that("Cox + strata (baseline): strata appear in the final coxph", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) +
                fp(nodes, df = 1, select = 1) + strata(grade),
              data = gbsg, family = "cox", verbose = FALSE)
  expect_true(inherits(fit, "coxph"))
  expect_false(is.null(fit$strata))
  expect_true(any(grepl("strata\\(", as.character(fit$formula))))
})

test_that("Fine-Gray strata_action modes select the correct model structure", {
  skip_on_cran()
  dfg <- make_finegray_data(n = 800, seed = 71, strong_strata = TRUE)

  fit_both <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1) + strata(sx),
                   data = dfg,
                   family = finegray_family(etype = "relapse", strata_action = "both"),
                   verbose = FALSE)
  fit_cens <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1) + strata(sx),
                   data = dfg,
                   family = finegray_family(etype = "relapse", strata_action = "censoring"),
                   verbose = FALSE)
  fit_base <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1) + strata(sx),
                   data = dfg,
                   family = finegray_family(etype = "relapse", strata_action = "baseline"),
                   verbose = FALSE)

  # "both": baseline stratified -> inner Cox has strata
  expect_equal(n_strata(fit_both), nlevels(dfg$sx))
  # "censoring": baseline NOT stratified -> inner Cox has no strata
  expect_equal(n_strata(fit_cens), 0L)
  # "baseline": baseline stratified -> inner Cox has strata
  expect_equal(n_strata(fit_base), nlevels(dfg$sx))

  # censoring/baseline differ in the IPCW expansion (pooled vs stratified),
  # so the number of expanded pseudo-observations differs
  expect_true(fit_base$n != fit_cens$n)

  # coefficients differ across modes when strata affect the processes
  expect_false(isTRUE(all.equal(unname(coef(fit_both)), unname(coef(fit_cens)))))
  expect_false(isTRUE(all.equal(unname(coef(fit_cens)), unname(coef(fit_base)))))
})

test_that("stratified Fine-Gray cumulative-incidence prediction routes strata", {
  skip_on_cran()
  dfg <- make_finegray_data(n = 600, seed = 711, strong_strata = TRUE)
  fit <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1) + strata(sx),
              data = dfg,
              family = finegray_family(etype = "relapse", strata_action = "baseline"),
              verbose = FALSE)
  nd <- data.frame(z1 = c(0, 0), sx = levels(dfg$sx)[1:2])
  cif <- predict(fit, newdata = nd, type = "response",
                 times = stats::median(dfg$obs), se.fit = FALSE)
  expect_length(cif, 2L)
  expect_true(all(is.finite(cif)))
  expect_true(all(cif >= 0 & cif <= 1))
})

test_that("survreg + strata: per-stratum scales match native survreg", {
  skip_on_cran()
  dsr <- make_survreg_data(n = 600, seed = 72)
  fit <- mfp2(Surv(time, status) ~ fp(x1, df = 1, select = 1) + strata(sx),
              data = dsr, family = survreg_family(dist = "weibull"), verbose = FALSE)
  ref <- survival::survreg(Surv(time, status) ~ x1 + strata(sx),
                           data = dsr, dist = "weibull")
  expect_length(fit$scale, 2L)
  expect_equal(sort(unname(fit$scale)), sort(unname(ref$scale)), tolerance = TOL_MED)
})

test_that("mfpi + strata: strata carried into finegray interaction models", {
  skip_on_cran()
  dfg <- make_mfpi_finegray()
  fit <- mfpi(Surv(obs, ev) ~ grp + fp(z1) + strata(sx), data = dfg,
              group_var = "grp", interaction_vars = "z1",
              interaction_forms = c(z1 = "fp1"),
              family = finegray_family(etype = "relapse"),
              p_interact = 0.1, verbose = FALSE)
  int <- find_fit2(fit$best_interaction_model, "coxph")
  expect_false(is.null(int))
  expect_false(is.null(int$strata))
  expect_true(any(grepl("strata\\(", as.character(int$formula))))
})

test_that("mfpi + strata: per-stratum scales carried into survreg interaction models", {
  skip_on_cran()
  dsr <- make_mfpi_survreg()
  fit <- mfpi(Surv(time, status) ~ grp + fp(x1) + strata(sx), data = dsr,
              group_var = "grp", interaction_vars = "x1",
              interaction_forms = c(x1 = "fp1"),
              family = survreg_family(dist = "weibull"),
              p_interact = 0.1, verbose = FALSE)
  int <- find_fit2(fit$best_interaction_model, "survreg")
  expect_false(is.null(int))
  expect_length(int$scale, 2L)
})
