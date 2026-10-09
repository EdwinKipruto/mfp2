make_reference_gee_data <- function(K = 40L, m = 4L) {
  set.seed(72)
  n <- K * m
  id <- rep(seq_len(K), each = m)
  x1 <- stats::runif(n, 1, 5)
  x2 <- stats::runif(n, 1, 4)
  y <- 2 + 2 / x1 + 0.5 * x2 + stats::rnorm(K)[id] +
    stats::rnorm(n, sd = 0.3)
  data.frame(y, x1, x2, id)
}


test_that("reference dispersion and QBIC count reach cached FP score changes", {
  obj <- list(logl = -12, df = 3, rank = 3,
              coefficients = c(a = 1, b = 1, c = 1),
              reference_dispersion = 2, qbic_n = 40)
  original <- calculate_gee_metrics(obj, n_obs = 10)
  shifted <- mfp2_shift_selection_df(original, 1, n_obs = 10, qbic_n = 40)
  direct <- calculate_gee_metrics(obj, n_obs = 10, df_additional = 1)
  expect_equal(shifted, direct)
  expect_equal(unname(original["aic"]), 12 + 6)
  expect_equal(unname(original["bic"]), 12 + log(40) * 3)
  expect_equal(unname(original["df_resid"]), 7)
})

test_that("a supplied dispersion is used for every GEE information score", {
  skip_if_not_installed("geepack")
  dat <- make_reference_gee_data()
  common <- list(
    stats::as.formula("y ~ x1 + x2"), data = dat,
    family = gee_family(gaussian(), corstr = "exchangeable",
                        reference_dispersion = 2),
    id = dat$id, df = 1, select = 1, center = FALSE,
    xorder = "original", verbose = FALSE
  )
  clusters <- do.call(mfp2, c(common, list(criterion = "bic")))
  observations_args <- common
  observations_args$family <- gee_family(
    gaussian(), corstr = "exchangeable", reference_dispersion = 2,
    qbic_penalty = "observations"
  )
  observations <- do.call(mfp2, c(observations_args, list(criterion = "bic")))
  expect_identical(clusters$reference_dispersion_source, "supplied")
  expect_null(clusters$reference_powers)
  expect_equal(clusters$reference_dispersion, 2)
  expect_equal(clusters$mfp_selection_score,
               -clusters$mfp_logl + log(length(unique(dat$id))) *
                 clusters$mfp_selection_df)
  expect_equal(observations$mfp_selection_score,
               -observations$mfp_logl + log(nrow(dat)) *
                 observations$mfp_selection_df)
  expect_equal(observations$linear_selection_score,
               -observations$linear_logl + log(nrow(dat)) *
                 observations$linear_df)
  expect_equal(unname(coef(clusters)), unname(coef(observations)))
  aic <- do.call(mfp2, c(common, list(criterion = "aic")))
  expect_equal(aic$mfp_selection_score,
               -aic$mfp_logl + 2 * aic$mfp_selection_df)
})

test_that("the GEE reference scale is estimated once from maximal FP terms", {
  skip_if_not_installed("geepack")
  dat <- make_reference_gee_data()
  fit <- mfp2(y ~ x1 + x2, data = dat,
              family = gee_family(gaussian(), corstr = "exchangeable"),
              id = dat$id, df = 2,
              powers = list(x1 = c(1, 2), x2 = c(1, 2)),
              criterion = "aic", xorder = "original", verbose = FALSE)
  expect_identical(fit$reference_dispersion_source, "estimated")
  expect_true(is.finite(fit$reference_dispersion))
  expect_gt(fit$reference_dispersion, 0)
  expect_equal(unname(fit$reference_powers$x1), 2)
  expect_equal(unname(fit$reference_powers$x2), 2)
  expect_true(all(vapply(fit$reference_powers, function(p) p %in% c(1, 2),
                         logical(1L))))
  expect_equal(fit$mfp_selection_score,
               -2 * fit$mfp_logl / fit$reference_dispersion +
                 2 * fit$mfp_selection_df)
  expect_output(print(summary(fit)), "Selection dispersion:")
})

test_that("invalid or inapplicable reference dispersion is rejected", {
  skip_if_not_installed("geepack")
  dat <- make_reference_gee_data()
  args <- list(stats::as.formula("y ~ x1"), data = dat,
               family = gee_family(gaussian()), id = dat$id,
               df = 1, verbose = FALSE)
  expect_error(gee_family(gaussian(), reference_dispersion = 0),
               "positive finite")
  args$family <- gee_family(gaussian(), reference_dispersion = 2)
  expect_error(do.call(mfp2, c(args, list(criterion = "pvalue"))),
               "requires GEE AIC or BIC")
})
