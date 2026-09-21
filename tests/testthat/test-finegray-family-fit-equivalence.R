# Area 1b: finegray_family() fitting correctness vs a manual
# survival::finegray() + weighted survival::coxph() pipeline.

test_that("mfp2 Fine-Gray matches a manual finegray + weighted coxph fit", {
  skip_on_cran()
  dat <- make_finegray_data(n = 600, seed = 9)
  dat$id <- seq_len(nrow(dat))

  fit <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1, scale = 1, shift = 0),
              data = dat, family = finegray_family(etype = "relapse"),
              center = FALSE, verbose = FALSE)

  fg <- survival::finegray(Surv(obs, ev) ~ z1 + id, data = dat, etype = "relapse")
  man <- survival::coxph(Surv(fgstart, fgstop, fgstatus) ~ z1 + cluster(id),
                         data = fg, weights = fgwt)

  expect_equal(fit$family_string, "finegray")
  expect_true(inherits(fit, "coxph"))

  # coefficient, partial log-likelihood, and robust SE should all agree
  expect_equal(unname(coef(fit)["z1.1"]), unname(coef(man)["z1"]),
               tolerance = TOL_TIGHT)
  expect_equal(unname(fit$loglik[2]), unname(man$loglik[2]), tolerance = TOL_TIGHT)
  expect_equal(sqrt(diag(vcov(fit)))[["z1.1"]], sqrt(man$var[1, 1]),
               tolerance = TOL_MED)

  # linear predictors agree up to the shared centring constant. The manual
  # coxph lives on the expanded pseudo-observations, so reconstruct its lp on
  # the ORIGINAL rows from the (matching) coefficient before comparing.
  lp_fit <- predict(fit, type = "lp")
  lp_man <- unname(coef(man)["z1"]) * dat$z1
  expect_length(lp_fit, nrow(dat))
  expect_equal(sd(lp_fit - lp_man), 0, tolerance = TOL_MED)   # differ only by a constant
  expect_equal(cor(lp_fit, lp_man), 1, tolerance = TOL_TIGHT)
})

test_that("Fine-Gray stores expected metadata and type = 'risk' equals exp(lp)", {
  skip_on_cran()
  dat <- make_finegray_data(n = 500, seed = 11)
  fit <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1),
              data = dat, family = finegray_family(etype = "relapse"),
              verbose = FALSE)

  expect_equal(fit$mfp2_finegray_event, "relapse")
  expect_false(is.null(fit$mfp2_finegray_row_map))
  expect_equal(fit$mfp2_finegray_n_original, nrow(dat))
  # prepared slot stripped after fitting
  expect_null(fit$family$prepared)

  lp <- predict(fit, type = "lp")
  rk <- predict(fit, type = "risk")
  expect_equal(max_abs_diff(rk, exp(lp)), 0, tolerance = TOL_TIGHT)
})

test_that("Fine-Gray response prediction returns cumulative incidence at times", {
  skip_on_cran()
  dat <- make_finegray_data(n = 500, seed = 12)
  fit <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1),
              data = dat, family = finegray_family(etype = "relapse"),
              verbose = FALSE)

  times <- as.numeric(stats::quantile(dat$obs, c(0.25, 0.5, 0.75)))
  cif <- predict(fit, type = "response", times = times, se.fit = FALSE)
  expect_equal(dim(cif), c(nrow(dat), length(times)))
  expect_true(all(is.finite(cif)))
  expect_true(all(cif >= 0 & cif <= 1))
  expect_true(all(apply(cif, 1L, function(x) all(diff(x) >= -1e-12))))

  cif_new <- predict(fit, newdata = dat[1:6, "z1", drop = FALSE],
                     type = "response", times = times, se.fit = FALSE)
  expect_equal(unname(cif_new), unname(cif[1:6, , drop = FALSE]),
               tolerance = 1e-8)
  expect_error(predict(fit, type = "response"), "requires `times`")
})

test_that("etype selects the modelled cause (relapse vs death give different fits)", {
  skip_on_cran()
  dat <- make_finegray_data(n = 700, seed = 13)
  f_rel <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1),
                data = dat, family = finegray_family(etype = "relapse"),
                verbose = FALSE)
  f_death <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1),
                  data = dat, family = finegray_family(etype = "death"),
                  verbose = FALSE)
  expect_equal(f_rel$mfp2_finegray_event, "relapse")
  expect_equal(f_death$mfp2_finegray_event, "death")
  expect_false(isTRUE(all.equal(unname(coef(f_rel)), unname(coef(f_death)))))
})
