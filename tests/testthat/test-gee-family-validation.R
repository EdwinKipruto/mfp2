# gee_family() constructor and GEE input validation

test_that("gee_family() accepts supported families, links, corstr, std.err", {
  f <- gee_family()
  expect_s3_class(f, "mfp2_gee_family")
  expect_s3_class(f, "mfp2_family")
  expect_identical(f$family, "gee")
  expect_identical(f$response_family$family, "gaussian")
  expect_identical(f$corstr, "exchangeable")
  expect_identical(f$std.err, "san.se")
  expect_null(f$scale.fix)
  expect_null(f$scale.value)
  expect_null(f$pvalue_method)

  expect_identical(gee_family(binomial())$response_family$family, "binomial")
  expect_identical(gee_family("poisson")$response_family$family, "poisson")
  expect_identical(gee_family(Gamma(link = "log"))$response_family$family, "Gamma")
  expect_identical(gee_family(corstr = "ar1")$corstr, "ar1")
  expect_identical(gee_family(corstr = "exchangeable")$corstr, "exchangeable")
  expect_identical(gee_family(std.err = "jack")$std.err, "jack")
  expect_identical(gee_family(std.err = "j1s")$std.err, "j1s")
  expect_identical(gee_family(std.err = "fij")$std.err, "fij")
})

test_that("gee_family() rejects invalid settings", {
  expect_error(gee_family("inverse.gaussian"), "must be one of")
  expect_error(gee_family("negbin"), "must be one of")
  expect_error(gee_family(corstr = "unstructured"))
  expect_error(gee_family(std.err = "bogus"))
  expect_error(gee_family(pvalue_method = "stata"), "unused argument")
  expect_error(gee_family(scale.fix = FALSE), "unused argument")
  expect_error(gee_family(scale.value = 1), "unused argument")
})

test_that("mfp2() requires id for GEE and rejects it elsewhere", {
  set.seed(1)
  n <- 40
  d <- data.frame(y = rnorm(n), x1 = rnorm(n), id = rep(1:10, each = 4))

  expect_error(
    mfp2(y ~ x1, data = d, family = gee_family(), df = 1, verbose = FALSE),
    "`id` is required"
  )
  # id supplied to a non-GEE, non-finegray family is rejected
  expect_error(
    mfp2(y ~ x1, data = d, family = "gaussian", id = id, df = 1, verbose = FALSE),
    "only used with"
  )
  # waves without GEE is rejected
  expect_error(
    mfp2(y ~ x1, data = d, family = "gaussian", waves = x1, df = 1, verbose = FALSE),
    "only used with"
  )
})

test_that("GEE rejects non-contiguous clusters and bad waves", {
  set.seed(1)
  n <- 40
  x1 <- rnorm(n)
  y <- rnorm(n)
  good_id <- rep(1:10, each = 4)
  scrambled <- rep(1:10, times = 4)          # not contiguous
  d <- data.frame(y, x1, id = good_id, sid = scrambled)

  expect_error(
    mfp2(y ~ x1, data = d, family = gee_family(), id = sid,
         df = 1, select = 1, center = FALSE, verbose = FALSE),
    "contiguous"
  )

  # duplicate waves within a cluster
  bad_wave <- rep(c(1, 1, 2, 3), times = 10)
  d$bw <- bad_wave
  expect_error(
    mfp2(y ~ x1, data = d, family = gee_family(corstr = "ar1"), id = id,
         waves = bw, df = 1, select = 1, center = FALSE, verbose = FALSE),
    "unique within"
  )
})

test_that("mfpi() rejects GEE families (GEE is mfp2-only)", {
  set.seed(1)
  n <- 120
  d <- data.frame(y = rnorm(n), x = runif(n, 1, 5),
                  g = rbinom(n, 1, 0.5), id = rep(1:30, each = 4))
  expect_error(
    mfpi(y ~ fp(x, df = 4) + g, data = d, family = gee_family(gaussian()),
         group_var = "g", verbose = FALSE),
    "does not support GEE"
  )
})
