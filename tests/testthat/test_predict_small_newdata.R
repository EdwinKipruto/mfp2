# Prediction data with one row, or with only two distinct values of a
# predictor, previously triggered the "binary predictor" shortcut in the FP
# transformation: the fitted power was not applied. This corrupted term and
# contrast predictions for every FP term and also full-model predictions for
# ACD terms with an FP component.

small_newdata_fixture <- function() {
  set.seed(12)
  n <- 300
  d <- data.frame(x = stats::rgamma(n, 2, 0.3), x2 = stats::runif(n))
  d$y <- 1 + 0.6 * log(d$x) + 0.3 * d$x2 + stats::rnorm(n, sd = 0.3)
  d
}

expect_small_newdata_consistent <- function(fit, nd, v) {
  full <- predict(fit, newdata = nd)
  terms <- predict(fit, type = "terms", terms = v, newdata = nd,
                   terms_seq = "data", add_intercept = FALSE)[[v]]$value
  ref <- stats::setNames(list(nd[[v]][2]), v)
  contrasts <- predict(fit, type = "contrasts", terms = v, newdata = nd,
                       terms_seq = "data", ref = ref)[[v]]$value
  for (rows in list(1L, c(1L, 2L, 1L))) {
    sub <- nd[rows, , drop = FALSE]
    expect_equal(unname(predict(fit, newdata = sub)), unname(full[rows]))
    expect_equal(
      predict(fit, type = "terms", terms = v, newdata = sub,
              terms_seq = "data", add_intercept = FALSE)[[v]]$value,
      terms[rows]
    )
    expect_equal(
      predict(fit, type = "contrasts", terms = v, newdata = sub,
              terms_seq = "data", ref = ref)[[v]]$value,
      contrasts[rows]
    )
  }
  expect_equal(contrasts[2], 0)
}

test_that("FP term and contrast predictions work for one or two distinct rows", {
  d <- small_newdata_fixture()
  fit <- mfp2(y ~ fp(x, df = 4, select = 1, force_max_fp = TRUE) + x2,
              data = d, select = 1, verbose = FALSE)
  expect_false(identical(unname(fit$fp_powers$x), 1))
  nd <- d[c(3, 10, 50, 80), ]
  expect_small_newdata_consistent(fit, nd, "x")

  # Independent check of a single-row term value from the fitted basis.
  p <- fit$fp_powers$x
  basis <- function(x) {
    x <- x + fit$transformations["x", "shift"]
    b <- if (p[1] == 0) log(x) else x^p[1]
    cbind(b, if (p[2] == p[1]) b * log(x) else if (p[2] == 0) log(x) else x^p[2])
  }
  centred <- basis(nd$x[1]) - colMeans(basis(d$x))
  expected <- unname(drop(centred %*% stats::coef(fit)[c("x.1", "x.2")]))
  expect_equal(
    predict(fit, type = "terms", terms = "x", newdata = nd[1, ],
            terms_seq = "data", add_intercept = FALSE)$x$value,
    expected
  )
})

test_that("ACD predictions are correct for a single new observation", {
  d <- small_newdata_fixture()
  fit <- mfp2(y ~ fp(x, acd = TRUE, select = 1, force_max_fp = TRUE) + x2,
              data = d, select = 1, verbose = FALSE)
  # The bug needs an FP component on x itself, not only on ACD(x).
  expect_false(is.na(fit$fp_powers$x[[1]]))
  expect_small_newdata_consistent(fit, d[c(3, 10, 50, 80), ], "x")
})

test_that("prediction values do not depend on data frame row names", {
  d <- small_newdata_fixture()
  fit <- mfp2(y ~ fp(x, df = 4, select = 1, force_max_fp = TRUE) + x2,
              data = d, select = 1, verbose = FALSE)
  nd <- d[c(3, 10, 50, 80), ]
  nd_default <- nd
  row.names(nd_default) <- NULL
  # Changing row names must not change prediction values, while each result
  # keeps the observation names of its own input.
  pred_df <- predict(fit, newdata = nd)
  pred_default <- predict(fit, newdata = nd_default)
  expect_equal(as.numeric(pred_default), as.numeric(pred_df))
  expect_identical(names(pred_df), row.names(nd))
  expect_identical(names(pred_default), row.names(nd_default))
  se_df <- predict(fit, newdata = nd, se.fit = TRUE)
  se_default <- predict(fit, newdata = nd_default, se.fit = TRUE)
  for (component in c("fit", "se.fit")) {
    expect_equal(as.numeric(se_default[[component]]), as.numeric(se_df[[component]]))
    expect_identical(names(se_df[[component]]), row.names(nd))
    expect_identical(names(se_default[[component]]), row.names(nd_default))
  }
  expect_equal(
    predict(fit, type = "terms", terms = "x", newdata = nd_default, terms_seq = "data"),
    predict(fit, type = "terms", terms = "x", newdata = nd, terms_seq = "data")
  )
  expect_equal(
    predict(fit, type = "contrasts", terms = "x", newdata = nd_default),
    predict(fit, type = "contrasts", terms = "x", newdata = nd)
  )
})

test_that("multinomial term predictions are correct for a single row", {
  d <- small_newdata_fixture()
  set.seed(5)
  d$cls <- factor(sample(c("A", "B", "C"), nrow(d), replace = TRUE))
  fit <- mfp2(cls ~ fp(x, df = 4, select = 1, force_max_fp = TRUE) + x2,
              data = d, family = multinomial_family(), select = 1,
              verbose = FALSE)
  nd <- d[c(3, 10, 50), ]
  many <- predict(fit, type = "terms", terms = "x", newdata = nd,
                  terms_seq = "data", add_intercept = FALSE)$x
  one <- predict(fit, type = "terms", terms = "x", newdata = nd[1, ],
                 terms_seq = "data", add_intercept = FALSE)$x
  expect_equal(one$value, many$value[many$variable == nd$x[1]])
})
