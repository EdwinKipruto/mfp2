# Regression tests: continuous interaction variables that are linear (df = 1)
# in the adjustment model but have an FP interaction form must keep a
# positivity shift. Previously their shift was forced to 0, so the FP
# interaction search took logs/powers of non-positive values and failed.

make_linear_adj_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(g = rbinom(n, 1, 0.5), z = runif(n, 1, 5))
  d$x3 <- sample(c(-1, 0, 1), n, TRUE)  # 3 distinct values -> automatic df = 1
  d$xn <- rnorm(n)                       # continuous with negatives
  d$w <- rexp(n)
  d$w[sample(n, 60)] <- 0
  d$y <- d$xn * (1 + d$g) + d$x3 + log(d$z) + rnorm(n)
  d
}

fit_linear_adj <- function(formula, d, ...) {
  suppressWarnings(mfpi(formula, data = d, group_var = "g",
                        p_interact = 1, verbose = FALSE, ...))
}

test_that("df = 1 interaction variable behaves like the default-df case", {
  # The interaction stage excludes the focal variable from the adjustment
  # block, so with the same shift/scale the interaction test must be identical
  # whether the adjustment term is linear (df = 1) or FP (default df). Any
  # outcome, including data-specific estimability errors, must match.
  d <- make_linear_adj_data()
  outcome <- function(formula, flex, form) {
    tryCatch({
      fit <- fit_linear_adj(formula, d, interaction_vars = "xn",
                            interaction_forms = c(xn = form), flex = flex)
      list(shift = fit$shift[["xn"]],
           pvalue = fit$all_model_metrics$pvalue[1L],
           powers = fit$all_model_metrics$fp_powers_int[[1L]])
    }, error = function(e) conditionMessage(e))
  }
  for (flex in c("flex1", "flex2", "flex3", "flex4")) {
    for (form in c("fp1", "fp2")) {
      linear_adj <- outcome(y ~ g + fp(xn, df = 1) + fp(z), flex, form)
      fp_adj <- outcome(y ~ g + fp(xn) + fp(z), flex, form)
      expect_equal(linear_adj, fp_adj, info = paste(flex, form))
      if (is.list(linear_adj)) {
        expect_true(min(d$xn + linear_adj$shift) > 0, info = paste(flex, form))
      }
    }
  }
})

test_that("automatic df = 1 (3 distinct values) interaction variable works", {
  d <- make_linear_adj_data()
  for (flex in c("flex1", "flex2", "flex3", "flex4")) {
    fit <- fit_linear_adj(y ~ g + x3 + fp(z), d, interaction_vars = "x3",
                          interaction_forms = c(x3 = "fp1"), flex = flex)
    expect_true(min(d$x3 + fit$shift[["x3"]]) > 0, info = flex)
  }
})

test_that("explicit shift is honoured for a df = 1 FP interaction variable", {
  d <- make_linear_adj_data()
  fit <- fit_linear_adj(y ~ g + fp(xn, df = 1, shift = 5) + fp(z), d,
                        interaction_vars = "xn")
  expect_identical(fit$shift[["xn"]], 5)
})

test_that("insufficient explicit shift gives the informative error", {
  d <- make_linear_adj_data()
  expect_error(
    fit_linear_adj(y ~ g + fp(xn, df = 1, shift = 0.1) + fp(z), d,
                   interaction_vars = "xn"),
    "Shifting factors are insufficient"
  )
})

test_that("predictions from newdata match training predictions", {
  d <- make_linear_adj_data()
  fit <- fit_linear_adj(y ~ g + fp(xn, df = 1) + fp(z), d,
                        interaction_vars = "xn", flex = "flex3")
  p_new <- predict(fit, newdata = d, type = "link", model = "all",
                   se.fit = FALSE)$predictions$fit
  p_train <- predict(fit, type = "link", model = "all",
                     se.fit = FALSE)$predictions$fit
  expect_equal(p_new, p_train, tolerance = 1e-10)
})

test_that("linear interaction forms and zero-handled variables are unchanged", {
  d <- make_linear_adj_data()
  fit_lin <- fit_linear_adj(y ~ g + fp(xn, df = 1) + fp(z), d,
                            interaction_vars = "xn",
                            interaction_forms = c(xn = "linear"))
  expect_identical(fit_lin$shift[["xn"]], 0)

  fit_zero <- fit_linear_adj(y ~ g + fp(w, zero = TRUE, df = 1) + fp(z), d,
                             interaction_vars = "w")
  expect_identical(fit_zero$shift[["w"]], 0)
  expect_true(isTRUE(fit_zero$zero_vars[["w"]]))
})

test_that("matrix interface is covered", {
  d <- make_linear_adj_data()
  x <- as.matrix(d[, c("g", "xn", "z")])
  fit <- suppressWarnings(mfpi(x, d$y, group_var = "g", interaction_vars = "xn",
                               df = c(xn = 1), p_interact = 1, verbose = FALSE))
  expect_true(min(d$xn + fit$shift[["xn"]]) > 0)
})
