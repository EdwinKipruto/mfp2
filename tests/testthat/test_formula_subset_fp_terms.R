# Regression tests: a formula `subset` must not lose fp()/fp2() term settings.
# fp() stores its settings as attributes of the model-frame column, and row
# subsetting with `[` drops such attributes from plain numeric vectors.

make_subset_fp_data <- function(n = 300, seed = 21) {
  set.seed(seed)
  d <- data.frame(
    g = sample(c("ctl", "trt"), n, TRUE),
    x = runif(n, 0.5, 8),
    w = rexp(n) * 2,
    z = runif(n, 1, 5),
    f = factor(sample(c("p", "q", "r"), n, TRUE))
  )
  d$w[sample(n, 60)] <- 0
  d$y <- log(d$x) * (1 + (d$g == "trt")) + 0.3 * d$z + rnorm(n)
  d$keep <- d$z > 1.5
  d
}

test_that("subset helper keeps fp() term attributes", {
  d <- make_subset_fp_data()
  mf <- stats::model.frame(y ~ fp(x, df = 2, powers = c(0, 1)) + z, data = d)
  sub <- mfp2:::subset_formula_model_frame(mf, which(d$keep))
  fp_col <- sub[["fp(x, df = 2, powers = c(0, 1))"]]
  expect_identical(attr(fp_col, "name"), "x")
  expect_identical(attr(fp_col, "df"), 2)
  expect_identical(attr(fp_col, "powers"), c(0, 1))
  expect_length(fp_col, sum(d$keep))
  expect_identical(as.numeric(fp_col), d$x[d$keep])
})

test_that("mfp2() formula subset works with fp() and matches fitting on subset data", {
  d <- make_subset_fp_data()
  ds <- d[d$keep, ]
  f <- y ~ fp(x, df = 2, shift = 0, scale = 1) + fp(z, shift = 0, scale = 1) + f
  fit_sub <- mfp2(f, data = d, subset = keep, verbose = FALSE)
  fit_ref <- mfp2(f, data = ds, verbose = FALSE)
  expect_equal(fit_sub$fp_powers, fit_ref$fp_powers)
  expect_equal(coef(fit_sub), coef(fit_ref), tolerance = 1e-8)
  expect_lte(length(fit_sub$fp_powers$x), 1L)  # df = 2 is honoured
})

test_that("mfpi() formula subset works with fp() and matches fitting on subset data", {
  d <- make_subset_fp_data()
  ds <- d[d$keep, ]
  f <- y ~ g + fp(x, shift = 0, scale = 1) +
    fp(w, zero = TRUE, shift = 0, scale = 1) + fp(z, shift = 0, scale = 1) + f
  args <- list(group_var = "g", interaction_vars = c("x", "w", "f"),
               p_interact = 1, verbose = FALSE)
  fit_sub <- suppressWarnings(do.call(mfpi, c(list(f, data = d, subset = d$keep), args)))
  fit_ref <- suppressWarnings(do.call(mfpi, c(list(f, data = ds), args)))
  expect_equal(fit_sub$all_model_metrics$pvalue, fit_ref$all_model_metrics$pvalue,
               tolerance = 1e-8)
  expect_equal(fit_sub$nobs, nrow(ds))
  expect_true(isTRUE(fit_sub$zero_vars[["w"]]))  # fp(w, zero = TRUE) is honoured
})

test_that("fp() settings come from the full-data frame", {
  d <- make_subset_fp_data()
  mf_full <- stats::model.frame(y ~ fp(x, df = 2) + z, data = d)
  mf <- mf_full[which(d$keep), , drop = FALSE]  # attributes dropped by `[`
  expect_null(attr(mf[["fp(x, df = 2)"]], "name"))
  settings <- mfp2:::formula_fp_settings_frame(mf, mf_full, 2L)
  expect_identical(attr(settings[[1L]], "name"), "x")
  expect_identical(attr(settings[[1L]], "df"), 2)
})
