test_that("resolve_mfp_criterion() uses family-dependent defaults", {
  expect_identical(mfp2:::resolve_mfp_criterion(NULL, "gee"), "aic")
  for (fam in c("gaussian", "binomial", "poisson", "cox", "negbin")) {
    expect_identical(mfp2:::resolve_mfp_criterion(NULL, fam), "pvalue")
  }
  expect_identical(mfp2:::resolve_mfp_criterion("bic", "gee"), "bic")
  expect_identical(mfp2:::resolve_mfp_criterion("pvalue", "gee"), "pvalue")
  expect_identical(mfp2:::resolve_mfp_criterion("aic", "gaussian"), "aic")
  expect_identical(mfp2:::resolve_mfp_criterion(c("pvalue", "aic", "bic"), "gee"), "pvalue")
  expect_error(mfp2:::resolve_mfp_criterion("qic", "gee"))
  expect_error(mfp2:::resolve_mfp_criterion(1, "gaussian"))
})

test_that("mfp2() defaults to AIC for GEE and p-value otherwise", {
  set.seed(101)
  n_id <- 40; m <- 4
  d <- data.frame(id = rep(seq_len(n_id), each = m), x1 = runif(n_id * m, 1, 5))
  d$y <- 0.5 * log(d$x1) + rep(rnorm(n_id, sd = 0.5), each = m) + rnorm(n_id * m)

  gee_default <- mfp2(y ~ fp(x1, df = 2), data = d, family = gee_family(),
                      id = id, verbose = FALSE)
  expect_identical(gee_default$criterion_mfp, "aic")

  gee_pvalue <- mfp2(y ~ fp(x1, df = 2), data = d, family = gee_family(),
                     id = id, criterion = "pvalue", verbose = FALSE)
  expect_identical(gee_pvalue$criterion_mfp, "pvalue")

  x <- as.matrix(d[, "x1", drop = FALSE])
  gee_matrix <- mfp2(x, d$y, family = gee_family(), id = d$id,
                     df = 2, verbose = FALSE)
  expect_identical(gee_matrix$criterion_mfp, "aic")

  glm_default <- mfp2(y ~ fp(x1, df = 2), data = d, verbose = FALSE)
  expect_identical(glm_default$criterion_mfp, "pvalue")
})
