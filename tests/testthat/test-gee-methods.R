# S3 methods on a fitted GEE mfp2 object: GEE prediction SEs, confidence
# intervals, residuals, the reported GEE parameters, and plotting.
#
# geepack ships no predict()/confint() methods and its residuals()/fitted()
# return one-column matrices; mfp2 must therefore provide correct robust
# behaviour itself. Sandwich quantities are checked after applying Stata's
# K/(K-1) correction to vcov(geeglm).

make_gee_methods_data <- function(K = 60, m = 4, seed = 202) {
  set.seed(seed)
  n <- K * m
  id <- rep(seq_len(K), each = m)
  x1 <- stats::rnorm(n)
  x2 <- stats::runif(n, 1, 5)
  b <- stats::rnorm(K)[id]
  y <- 0.3 + 0.5 * x1 + 2 / x2 + b + stats::rnorm(n, sd = 0.5)
  data.frame(y, x1, x2, id)
}

test_that("confint() gives robust Wald intervals for GEE", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data()
  m2 <- mfp2(y ~ fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
               fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1),
             data = d, family = gee_family(gaussian(), corstr = "exchangeable"),
             id = id, xorder = "original", verbose = FALSE)
  gg <- geepack::geeglm(y ~ x1 + x2, data = d, family = gaussian(),
                        id = id, corstr = "exchangeable")

  z <- stats::qnorm(0.975)
  k <- length(gg$geese$clusz)
  corrected_vcov <- vcov(gg) * k / (k - 1)
  ref <- cbind(coef(gg) - z * sqrt(diag(corrected_vcov)),
               coef(gg) + z * sqrt(diag(corrected_vcov)))
  ci <- confint(m2)
  expect_true(is.matrix(ci))
  expect_equal(unname(ci), unname(ref), tolerance = TOL_MED)
  expect_identical(colnames(ci), c("2.5 %", "97.5 %"))

  # level and parm arguments
  ci90 <- confint(m2, level = 0.90)
  expect_true(all(ci90[, 1] > ci[, 1]))                 # narrower
  expect_equal(nrow(confint(m2, parm = "x1.1")), 1L)
})

test_that("residuals() returns a plain vector for GEE (all types)", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data()
  m2 <- mfp2(y ~ fp(x1) + fp(x2), data = d,
             family = gee_family(gaussian(), corstr = "exchangeable"),
             id = id, verbose = FALSE)
  for (ty in c("deviance", "pearson", "working", "response")) {
    r <- residuals(m2, type = ty)
    expect_true(is.numeric(r))
    expect_null(dim(r))                                 # not a 1-column matrix
    expect_length(r, nrow(d))
  }
})

test_that("GEE deviance residuals match the marginal-mean definition (Stata fracplot)", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data()
  m2 <- mfp2(y ~ fp(x1) + fp(x2), data = d,
             family = gee_family(gaussian(), corstr = "exchangeable"),
             id = id, verbose = FALSE)
  # Deviance residuals are not provided by geepack; the component-plus-residual
  # convention (Stata fracplot) is sign(y - mu) * sqrt(unit deviance), using the
  # GEE marginal fitted mean and the family deviance function.
  fam <- gaussian()
  mu <- as.numeric(fitted(m2))
  y <- as.numeric(m2$y)
  ref <- sign(y - mu) * sqrt(pmax(fam$dev.resids(y, mu, rep(1, length(y))), 0))
  expect_equal(unname(residuals(m2, type = "deviance")), unname(ref),
               tolerance = 1e-10)
  # default type is deviance, as for GLMs
  expect_equal(unname(residuals(m2)), unname(ref), tolerance = 1e-10)
})

test_that("summary() reports the key GEE parameters, matching geepack", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data()
  m2 <- mfp2(y ~ fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
               fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1),
             data = d, family = gee_family(gaussian(), corstr = "exchangeable"),
             id = id, xorder = "original", verbose = FALSE)
  gg <- geepack::geeglm(y ~ x1 + x2, data = d, family = gaussian(),
                        id = id, corstr = "exchangeable")

  g <- summary(m2)$gee
  expect_false(is.null(g))
  expect_identical(g$corstr, "exchangeable")
  expect_identical(g$std_err, "san.se")
  expect_equal(unname(g$correlation), unname(gg$geese$alpha), tolerance = 1e-8)
  expect_equal(unname(g$scale), unname(gg$geese$gamma[1L]), tolerance = 1e-8)
  expect_equal(g$n_clusters, length(gg$geese$clusz))
  expect_equal(g$max_cluster_size, max(gg$geese$clusz))

  # the working-correlation and scale parameters carry their sandwich standard
  # errors, matching geepack (geese$valpha / geese$vgamma)
  expect_equal(unname(g$correlation_se),
               unname(sqrt(diag(as.matrix(gg$geese$valpha)))), tolerance = 1e-8)
  expect_equal(unname(g$scale_se),
               unname(sqrt(diag(as.matrix(gg$geese$vgamma)))[1L]),
               tolerance = 1e-8)

  # printed output labels the model GEE and renders the working-correlation
  # block with the one-line alpha (with SE), the scale (with SE), and the
  # "Number of clusters" label
  out <- capture.output(print(summary(m2)))
  expect_true(any(grepl("Model: Gaussian GEE (identity link)", out, fixed = TRUE)))
  expect_true(any(grepl("GEE Working Correlation", out, fixed = TRUE)))
  expect_true(any(grepl("exchangeable", out, fixed = TRUE)))
  expect_true(any(grepl("SE method:", out, fixed = TRUE) &
                    grepl("san.se", out, fixed = TRUE)))
  expect_true(any(grepl("Estimated correlation:", out, fixed = TRUE) &
                    grepl("alpha =", out, fixed = TRUE) &
                    grepl("(SE =", out, fixed = TRUE)))
  expect_true(any(grepl("Scale (dispersion):", out, fixed = TRUE) &
                    grepl("(SE =", out, fixed = TRUE)))
  expect_true(any(grepl("Number of clusters:", out, fixed = TRUE)))
  expect_false(any(grepl("^Clusters:", out)))

  # the plain print method also surfaces the GEE parameters identically
  out2 <- capture.output(print(m2))
  expect_true(any(grepl("Model: Gaussian GEE (identity link)", out2, fixed = TRUE)))
  expect_true(any(grepl("GEE Working Correlation", out2, fixed = TRUE)))
  expect_true(any(grepl("Number of clusters:", out2, fixed = TRUE)))
  expect_false(any(grepl("\\bNA\\b|NaN", c(out, out2))))
})

test_that("summary term SEs and CIs use the robust covariance", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data()
  m2 <- mfp2(y ~ fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
               fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1),
             data = d, family = gee_family(gaussian(), corstr = "exchangeable"),
             id = id, xorder = "original", verbose = FALSE)
  gg <- geepack::geeglm(y ~ x1 + x2, data = d, family = gaussian(),
                        id = id, corstr = "exchangeable")

  lt <- summary(m2)$linear_terms
  # The intercept is reported as the first linear-terms row, consistently with
  # the other families that estimate one.
  expect_identical(lt$term[[1L]], "(Intercept)")
  expect_identical(lt$variable[[1L]], "(Intercept)")

  # align the mfp2 term names (x1.1) to the geeglm coefficient names (x1),
  # including the intercept, and check the robust SEs match geepack
  k <- length(gg$geese$clusz)
  robust_se <- sqrt(diag(vcov(gg) * k / (k - 1)))
  se_map <- setNames(robust_se, names(robust_se))
  key <- sub("\\.1$", "", lt$term)
  expect_equal(unname(lt$se), unname(se_map[key]), tolerance = TOL_MED)

  # GEE reports the Wald chi-square (z^2) with a Pr(>|W|) column, following
  # geepack's convention rather than a normal-theory z.
  expect_identical(attr(lt, "statistic_label", exact = TRUE), "Wald")
  z <- coef(gg) / robust_se
  expect_equal(unname(lt$statistic), unname((z^2)[key]), tolerance = TOL_MED)

  out <- capture.output(print(summary(m2)))
  expect_true(any(grepl("Wald", out, fixed = TRUE)))
  expect_true(any(grepl("Pr(>|W|)", out, fixed = TRUE)))
})

test_that("selection, vcov, summary, and print use every supported SE method", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data(K = 35)

  covariance_field <- function(geese, method) {
    switch(
      method,
      "san.se" = geese$vbeta,
      "jack" = geese$vbeta.ajs,
      "j1s" = geese$vbeta.j1s,
      "fij" = geese$vbeta.fij
    )
  }
  parameter_covariance_field <- function(geese, method, parameter) {
    if (identical(parameter, "alpha")) {
      switch(
        method,
        "san.se" = geese$valpha,
        "jack" = geese$valpha.ajs,
        "j1s" = geese$valpha.j1s,
        "fij" = geese$valpha.fij
      )
    } else {
      switch(
        method,
        "san.se" = geese$vgamma,
        "jack" = geese$vgamma.ajs,
        "j1s" = geese$vgamma.j1s,
        "fij" = geese$vgamma.fij
      )
    }
  }

  for (method in c("san.se", "jack", "j1s", "fij")) {
    m2 <- mfp2(
      y ~ fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
        fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1),
      data = d,
      family = gee_family(
        gaussian(), corstr = "exchangeable", std.err = method
      ),
      id = id, xorder = "original", verbose = FALSE
    )
    gg <- geepack::geeglm(
      y ~ x1 + x2, data = d, family = gaussian(), id = id,
      corstr = "exchangeable", std.err = method
    )

    expected <- covariance_field(gg$geese, method)
    if (identical(method, "san.se")) {
      k <- length(gg$geese$clusz)
      expected <- expected * k / (k - 1)
    }
    expect_equal(unname(vcov(m2)), unname(expected), tolerance = TOL_MED,
                 info = method)
    gee_summary <- summary(m2)$gee
    expect_identical(gee_summary$std_err, method)
    expect_equal(
      unname(gee_summary$correlation_se),
      unname(sqrt(diag(as.matrix(
        parameter_covariance_field(gg$geese, method, "alpha")
      )))),
      tolerance = TOL_MED,
      info = method
    )
    expect_equal(
      unname(gee_summary$scale_se),
      unname(sqrt(diag(as.matrix(
        parameter_covariance_field(gg$geese, method, "gamma")
      )))[1L]),
      tolerance = TOL_MED,
      info = method
    )
    expect_true(is.finite(m2$mfp_selection_score), info = method)

    printed <- c(capture.output(print(m2)), capture.output(print(summary(m2))))
    expect_false(any(grepl("\\bNA\\b|NaN", printed)), info = method)
  }
})

test_that("nonlinear GEE summaries use finite Stata-style joint tests", {
  skip_if_not_installed("geepack")
  d <- make_gee_methods_data(K = 40)
  m2 <- mfp2(
    y ~ fp(
      x2, df = 4, select = 1, alpha = 1,
      center = FALSE, shift = 0, scale = 1
    ),
    data = d,
    family = gee_family(gaussian(), corstr = "exchangeable"),
    id = id, xorder = "original", verbose = FALSE
  )

  expect_true(is.matrix(m2[["x", exact = TRUE]]))
  expect_identical(colnames(m2[["x", exact = TRUE]]), names(coef(m2)))

  sm <- summary(m2)
  expect_identical(sm$family, "gee")
  expect_gt(nrow(sm$nonlinear_terms), 0L)
  expect_true(all(is.finite(sm$nonlinear_terms$lr_chisq)))
  expect_true(all(is.finite(sm$nonlinear_terms$p)))
  expect_identical(
    unname(sm$nonlinear_terms$df),
    unname(m2$fp_terms[sm$nonlinear_terms$variable, "df_final"])
  )
  expect_identical(
    attr(sm$nonlinear_terms, "test_label", exact = TRUE),
    "Wald-score chi-sq"
  )
  out <- capture.output(print(sm))
  expect_true(any(grepl("Wald-score chi-sq", out, fixed = TRUE)))
  expect_false(any(grepl("\\bNA\\b|NaN", out)))
})

test_that("plot()/fracplot() work for a GEE model with a nonlinear term", {
  skip_if_not_installed("geepack")
  skip_if_not_installed("ggplot2")
  d <- make_gee_methods_data()
  m2 <- mfp2(y ~ fp(x1) + fp(x2), data = d,
             family = gee_family(gaussian(), corstr = "exchangeable"),
             id = id, verbose = FALSE)
  p <- plot(m2)
  expect_true(is.list(p))
  expect_true(length(p) >= 1L)
  # each element is a ggplot
  expect_true(all(vapply(p, function(g) inherits(g, "ggplot"), logical(1))))
})
