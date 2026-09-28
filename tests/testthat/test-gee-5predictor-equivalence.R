# Five-predictor GEE equivalence: mfp2 vs a direct geepack::geeglm() fit.
#
# With fp(x, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) on each
# of five predictors, mfp2 fits a linear GEE model on the raw predictor scale
# and retains every variable. It must therefore reproduce geeglm(), with
# Stata's K/(K-1) correction applied to the sandwich covariance:
#   * coefficients (raw scale; mfp2 names the power-1 column "x.1" and, with
#     xorder = "original", keeps input order, so names align after stripping
#     the ".1" suffix),
#   * point predictions on the link and response scales, and
#   * prediction standard errors from the robust sandwich covariance
#     (geepack ships no predict method; the base predict.glm() SE path is
#     invalid for GEE, so mfp2 computes robust SEs itself).
#
# The equivalence is checked for every geepack response family
# (gaussian/binomial/poisson/Gamma) and every working correlation
# (independence/exchangeable/ar1).

make_gee5_data <- function(family_name, K = 60, m = 4, seed = 101) {
  set.seed(seed)
  n <- K * m
  id <- rep(seq_len(K), each = m)
  wave <- rep(seq_len(m), times = K)
  x1 <- stats::rnorm(n)
  x2 <- stats::runif(n, -1, 1)
  x3 <- stats::rnorm(n)
  x4 <- stats::runif(n, -2, 2)
  x5 <- stats::rnorm(n)
  b <- stats::rnorm(K)[id]                      # cluster random effect
  eta <- 0.3 + 0.5 * x1 - 0.3 * x2 + 0.2 * x3 + 0.25 * x4 - 0.1 * x5 + b
  y <- switch(
    family_name,
    gaussian = eta + stats::rnorm(n, sd = 0.5),
    binomial = stats::rbinom(n, 1, stats::plogis(eta)),
    poisson  = stats::rpois(n, exp(0.3 + 0.2 * x1 - 0.1 * x2 + 0.15 * x3 +
                                     0.1 * x4 + 0.3 * b)),
    Gamma    = stats::rgamma(n, shape = 2,
                             rate = 2 / exp(0.5 + 0.2 * x1 + 0.1 * x3 + 0.2 * b))
  )
  data.frame(y, x1, x2, x3, x4, x5, id, wave)
}

gee5_families <- list(
  gaussian = gaussian(),
  binomial = binomial(),
  poisson  = poisson(),
  Gamma    = Gamma(link = "log")
)

gee5_fp_formula <- y ~
  fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
  fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
  fp(x3, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
  fp(x4, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
  fp(x5, df = 1, select = 1, center = FALSE, shift = 0, scale = 1)

strip_fp1 <- function(nm) sub("\\.1$", "", nm)

for (fam_name in names(gee5_families)) {
  for (corstr in c("independence", "exchangeable", "ar1")) {

    test_that(sprintf("GEE mfp2 == geeglm for %s / %s (5 predictors)",
                      fam_name, corstr), {
      skip_if_not_installed("geepack")
      gfam <- gee5_families[[fam_name]]
      d <- make_gee5_data(fam_name)

      if (identical(corstr, "ar1")) {
        m2 <- mfp2(gee5_fp_formula, data = d,
                   family = gee_family(gfam, corstr = corstr),
                   id = id, waves = wave, xorder = "original", verbose = FALSE)
        gg <- geepack::geeglm(y ~ x1 + x2 + x3 + x4 + x5, data = d,
                              family = gfam, id = id, waves = wave,
                              corstr = corstr)
      } else {
        m2 <- mfp2(gee5_fp_formula, data = d,
                   family = gee_family(gfam, corstr = corstr),
                   id = id, xorder = "original", verbose = FALSE)
        gg <- geepack::geeglm(y ~ x1 + x2 + x3 + x4 + x5, data = d,
                              family = gfam, id = id, corstr = corstr)
      }

      # Retained object is a native geeglm and keeps all five predictors.
      expect_true(inherits(m2, "geeglm"))
      cm <- coef(m2)
      names(cm) <- strip_fp1(names(cm))
      expect_identical(names(cm), names(coef(gg)))          # order preserved
      expect_equal(unname(cm), unname(coef(gg)), tolerance = TOL_MED)

      # Point predictions on both scales, on held-out rows.
      nd <- d[1:20, ]
      expect_equal(unname(predict(m2, nd, type = "link")),
                   unname(predict(gg, nd, type = "link")), tolerance = TOL_MED)
      expect_equal(unname(predict(m2, nd, type = "response")),
                   unname(predict(gg, nd, type = "response")), tolerance = TOL_MED)

      # Robust prediction standard errors: mfp2 must reproduce
      # sqrt(diag(X0 V_robust X0')) on the link scale and its delta-method
      # transform on the response scale. geeglm has no predict(se.fit) of its
      # own, so the reference is computed directly from vcov(gg).
      X0 <- model.matrix(~ x1 + x2 + x3 + x4 + x5, nd)
      k <- length(gg$geese$clusz)
      V <- vcov(gg) * k / (k - 1)
      se_link_ref <- sqrt(diag(X0 %*% V %*% t(X0)))
      eta_ref <- as.numeric(X0 %*% coef(gg))
      se_resp_ref <- se_link_ref * abs(gfam$mu.eta(eta_ref))

      pl <- predict(m2, nd, type = "link", se.fit = TRUE)
      pr <- predict(m2, nd, type = "response", se.fit = TRUE)
      expect_equal(unname(pl$se.fit), unname(se_link_ref), tolerance = TOL_MED)
      expect_equal(unname(pr$se.fit), unname(se_resp_ref), tolerance = TOL_MED)

      # The fit returned alongside se.fit equals the point prediction.
      expect_equal(unname(pl$fit), unname(predict(m2, nd, type = "link")),
                   tolerance = TOL_MED)
      expect_equal(unname(pr$fit), unname(predict(m2, nd, type = "response")),
                   tolerance = TOL_MED)
    })
  }
}
