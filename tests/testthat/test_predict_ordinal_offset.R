# Ordinal full-model prediction applies the offset itself. Before the fix,
# predict() with newdata attached the offset as a formula helper column and
# stopped with "lacks internal-name metadata" for any ordinal model with an
# offset.

ordinal_offset_data <- function() {
  set.seed(404)
  n <- 240
  d <- data.frame(x1 = stats::runif(n, 0.5, 20), x2 = stats::runif(n, 0, 3),
                  o = stats::runif(n, -0.4, 0.4))
  latent <- 0.8 * log(d$x1) + 0.3 * d$x2 + d$o + stats::rlogis(n)
  d$y <- cut(latent, stats::quantile(latent, c(0, 0.3, 0.6, 1)),
             include.lowest = TRUE, labels = c("low", "mid", "high"),
             ordered_result = TRUE)
  d
}

test_that("ordinal predictions with a formula offset use newdata's offset", {
  skip_if_not_installed("rms")
  d <- ordinal_offset_data()
  fit <- mfp2(y ~ fp(x1, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
                fp(x2, df = 1, select = 1, center = FALSE, shift = 0, scale = 1) +
                offset(o),
              data = d, family = ordinal_family(), verbose = FALSE)
  native <- rms::orm(y ~ x1 + x2 + offset(o), data = d)
  slopes <- stats::coef(native)[c("x1", "x2")]
  rows <- c(3, 40, 111, 200)
  nd <- d[rows, c("x1", "x2", "o")]
  eta <- drop(as.matrix(nd[, c("x1", "x2")]) %*% slopes) + nd$o

  expect_equal(unname(predict(fit, newdata = nd, type = "link",
                              add_intercept = FALSE)), unname(eta),
               tolerance = 1e-6)
  cum <- stats::plogis(outer(eta, stats::coef(native)[1:2], "+"))
  probs <- cbind(1 - cum[, 1], cum[, 1] - cum[, 2], cum[, 2])
  expect_equal(unname(predict(fit, newdata = nd, type = "response")),
               unname(probs), tolerance = 1e-6)
  # Training-data predictions keep the fitted offset.
  expect_equal(unname(predict(fit, type = "link", add_intercept = FALSE)),
               unname(drop(as.matrix(d[, c("x1", "x2")]) %*% slopes) + d$o),
               tolerance = 1e-6)
  for (type in c("link", "response", "mean")) {
    training <- predict(fit, type = type, se.fit = TRUE)
    supplied <- predict(fit, newdata = nd, type = type, se.fit = TRUE)
    for (component in c("fit", "se.fit")) {
      expected <- if (is.matrix(training[[component]])) {
        training[[component]][rows, , drop = FALSE]
      } else {
        training[[component]][rows]
      }
      expect_equal(supplied[[component]], expected, tolerance = 1e-8)
    }
  }
})

test_that("ordinal matrix-interface predictions use newoffset", {
  skip_if_not_installed("rms")
  d <- ordinal_offset_data()
  x <- as.matrix(d[, c("x1", "x2")])
  fit <- mfp2(x, d$y, offset = d$o, family = ordinal_family(),
              df = c(x1 = 1, x2 = 1), select = 1, center = FALSE,
              shift = 0, scale = 1, verbose = FALSE)
  native <- rms::orm(d$y ~ x + offset(d$o))
  slopes <- stats::coef(native)[-(1:2)]
  rows <- c(7, 90, 150)
  eta <- drop(x[rows, ] %*% slopes) + d$o[rows]
  expect_equal(unname(predict(fit, newdata = x[rows, ], newoffset = d$o[rows],
                              add_intercept = FALSE)),
               unname(eta), tolerance = 1e-6)
  expect_error(predict(fit, newdata = x[rows, ]), "newoffset")
})
