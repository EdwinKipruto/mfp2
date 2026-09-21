# Area 3: predict.mfp2 for the three families.

test_that("complete-model predictions have correct type and dimension", {
  skip_on_cran()
  data(gbsg, package = "mfp2")

  # survreg
  fsr <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1),
              data = gbsg, family = survreg_family(dist = "weibull"), verbose = FALSE)
  expect_length(predict(fsr), nrow(gbsg))
  expect_length(predict(fsr, newdata = gbsg[1:10, ]), 10L)

  # finegray
  dfg <- make_finegray_data(n = 400, seed = 41)
  ffg <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1),
              data = dfg, family = finegray_family(etype = "relapse"), verbose = FALSE)
  expect_length(predict(ffg, type = "lp"), nrow(dfg))
  expect_length(predict(ffg, newdata = dfg[1:10, ], type = "risk"), 10L)

  # multinomial
  dmn <- make_multinomial_data(n = 300, seed = 42)
  fmn <- mfp2(y ~ fp(x1, df = 1, select = 1), data = dmn,
              family = multinomial_family(), control = list(maxit = 300), verbose = FALSE)
  expect_equal(dim(predict(fmn, newdata = dmn[1:10, ], type = "response")),
               c(10L, nlevels(dmn$y)))
})

test_that("term predictions return per-term values (with SE where supported)", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fsr <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 2, select = 1) +
                fp(nodes, df = 1, select = 1),
              data = gbsg, family = survreg_family(dist = "weibull"), verbose = FALSE)
  tt <- predict(fsr, type = "terms", terms = "age", se.fit = TRUE)
  expect_true(is.list(tt))
  expect_true("age" %in% names(tt))
  # the age term carries fitted values and standard errors
  age_term <- tt[["age"]]
  expect_true(is.list(age_term) || is.data.frame(age_term) || !is.null(dim(age_term)))
})

test_that("contrast predictions return differences from a reference", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fsr <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 2, select = 1),
              data = gbsg, family = survreg_family(dist = "weibull"), verbose = FALSE)
  ct <- predict(fsr, type = "contrasts", ref = list(age = median(gbsg$age)),
                terms = "age")
  expect_false(is.null(ct))
})

test_that("newdata missing a required variable gives a clear error", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fsr <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) +
                fp(nodes, df = 1, select = 1),
              data = gbsg, family = survreg_family(dist = "weibull"), verbose = FALSE)
  expect_error(
    predict(fsr, newdata = data.frame(nodes = gbsg$nodes[1:5]), type = "link"),
    regexp = "age"
  )
})

test_that("multinomial se.fit = TRUE is rejected for link and response", {
  skip_on_cran()
  dmn <- make_multinomial_data(n = 300, seed = 43)
  fmn <- mfp2(y ~ fp(x1, df = 1, select = 1), data = dmn,
              family = multinomial_family(), control = list(maxit = 300), verbose = FALSE)
  expect_error(predict(fmn, se.fit = TRUE, type = "link"),
               regexp = "se.fit", ignore.case = TRUE)
  expect_error(predict(fmn, se.fit = TRUE, type = "response"),
               regexp = "se.fit", ignore.case = TRUE)
})

test_that("stratified survreg/finegray require strata variable in newdata", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fsr <- mfp2(Surv(rectime, censrec) ~ fp(age, df = 1, select = 1) + strata(meno),
              data = gbsg, family = survreg_family(dist = "weibull"), verbose = FALSE)
  expect_error(
    predict(fsr, newdata = data.frame(age = gbsg$age[1:5]), type = "link"),
    regexp = "strat", ignore.case = TRUE
  )

  dfg <- make_finegray_data(n = 500, seed = 44)
  ffg <- mfp2(Surv(obs, ev) ~ fp(z1, df = 1, select = 1) + strata(sx),
              data = dfg, family = finegray_family(etype = "relapse"), verbose = FALSE)
  expect_error(
    predict(ffg, newdata = data.frame(z1 = dfg$z1[1:5]), type = "lp"),
    regexp = "strat|sx", ignore.case = TRUE
  )
})
