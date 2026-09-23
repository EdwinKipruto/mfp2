# coef.mfpi(), vcov.mfpi(), summary.mfpi() and plot.mfpi() now handle the
# multinomial per-logit coefficient matrix. These tests verify that the
# accessors expose one coefficient per (non-reference logit x model column),
# that the values and covariance match the underlying nnet::multinom() fit, and
# that plotting produces one set of curves per logit.

mn_acc_fit <- function() {
  set.seed(7); n <- 900; a <- runif(n, 1, 10)
  g <- factor(sample(c("0", "1"), n, TRUE))
  lp2 <- -0.5 + 0.10 * a + 0.9 * (g == "1") * a
  lp3 <-  0.3 - 0.15 * a + 0.5 * (g == "1")
  y <- factor(vapply(seq_len(n), function(i) {
    p1 <- 1 / (1 + exp(lp2[i]) + exp(lp3[i]))
    sample(c("A", "B", "C"), 1L,
           prob = c(p1, exp(lp2[i]) * p1, exp(lp3[i]) * p1))
  }, character(1)))
  dm <- data.frame(y = y, a = a, g = g)
  mfpi(y ~ g + fp(a), data = dm, group_var = "g", interaction_vars = "a",
       interaction_forms = c(a = "fp1"), family = multinomial_family(),
       control = list(maxit = 300), verbose = FALSE)
}

test_that("coef.mfpi exposes multinomial slopes but not nuisance intercepts", {
  skip_on_cran()
  fit <- mn_acc_fit()
  im <- fit$var_winners[["a"]]$fit$test_results$interaction_model
  cf <- coef(fit, term = "a")

  expect_s3_class(cf, "mfpi_coef")
  slope_columns <- setdiff(colnames(im$coefficient_matrix), "(Intercept)")
  # Q (= 2) non-reference logits x slope columns; class intercepts are omitted.
  expect_length(cf, nrow(im$coefficient_matrix) * length(slope_columns))
  # Values equal the flattened per-logit coefficient matrix (class-major order).
  expect_equal(as.numeric(cf),
               as.vector(t(im$coefficient_matrix[, slope_columns, drop = FALSE])),
               tolerance = 1e-10)
  expect_false(any(grepl("Intercept", names(cf), fixed = TRUE)))
  expect_false(any(grepl(":\\(Intercept\\)$", attr(cf, "raw_names"))))
  # Every coefficient name is prefixed by its non-reference class.
  expect_true(all(grepl("^(B|C) \\| ", names(cf))))
})

test_that("vcov.mfpi returns the per-logit covariance matching nnet::multinom", {
  skip_on_cran()
  fit <- mn_acc_fit()
  im <- fit$var_winners[["a"]]$fit$test_results$interaction_model
  V <- vcov(fit, term = "a")
  cf <- coef(fit, term = "a")

  expect_true(is.matrix(V))
  expect_equal(dim(V), rep(length(cf), 2L))
  expect_identical(rownames(V), names(cf))
  expect_identical(colnames(V), names(cf))
  # Numerically identical to the underlying multinom covariance (reordered to
  # the accessor's class-major layout via the stored raw names).
  raw <- attr(V, "raw_names")
  Vraw <- vcov(im$fit)
  expect_equal(as.numeric(V), as.numeric(Vraw[raw, raw]), tolerance = 1e-12)
})

test_that("summary.mfpi builds readable outcome-specific multinomial displays", {
  skip_on_cran()
  fit <- mn_acc_fit()
  s <- summary(fit)
  out <- capture.output(print(s))
  expect_identical(nrow(s$regression_displays$a$info$adjustment_table), 0L)
  # The regression display names both non-reference outcomes against reference A
  # and separates the coefficient components within each outcome block.
  expect_true(any(grepl("Outcome: B vs A", out, fixed = TRUE)))
  expect_true(any(grepl("Outcome: C vs A", out, fixed = TRUE)))
  expect_true(any(grepl("Group-specific FP terms:", out, fixed = TRUE)))
  expect_true(any(grepl("Group-variable coefficients:", out, fixed = TRUE)))
  expect_false(any(grepl("Adjustment coefficients:", out, fixed = TRUE)))
  expect_true(any(grepl("S.E.", out, fixed = TRUE)))
})

test_that("coef print shows a per-logit block layout", {
  skip_on_cran()
  fit <- mn_acc_fit()
  out <- capture.output(print(coef(fit, term = "a")))
  expect_true(any(grepl("reference class: A", out)))
  expect_true(any(grepl("Logit B vs A", out)))
  expect_true(any(grepl("Logit C vs A", out)))
})

test_that("plot.mfpi returns one set of curves per non-reference logit", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- mn_acc_fit()
  pl <- plot(fit, auto_print = FALSE, plot_type = "fitted")
  expect_true("a" %in% names(pl))
  # One entry per non-reference logit, each a list of per-group-contrast plots.
  expect_setequal(names(pl[["a"]]), c("logit_B", "logit_C"))
  expect_s3_class(pl[["a"]][["logit_B"]][[1]], "ggplot")
  expect_s3_class(pl[["a"]][["logit_C"]][[1]], "ggplot")
})

test_that("plot.mfpi difference plots also resolve per logit", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- mn_acc_fit()
  pl <- plot(fit, auto_print = FALSE, plot_type = "difference", show_ci_diff = TRUE)
  expect_setequal(names(pl[["a"]]), c("logit_B", "logit_C"))
  expect_s3_class(pl[["a"]][["logit_C"]][[1]], "ggplot")
})
