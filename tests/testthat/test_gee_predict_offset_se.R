gee_offset_prediction_data <- function() {
  set.seed(31)
  K <- 40
  id <- rep(seq_len(K), each = 4L)
  n <- length(id)
  d <- data.frame(id = id, x1 = stats::runif(n, 1, 5),
                  x2 = stats::runif(n, 0.5, 3), expo = stats::runif(n, 0.5, 2))
  d$cnt <- stats::rpois(n, d$expo * exp(0.2 + 0.6 / d$x1 + 0.3 * d$x2 +
                                          stats::rnorm(K, sd = 0.3)[id]))
  d
}

expect_gee_prediction_rows <- function(fit, pred_link, pred_resp, rows,
                                       offset) {
  native <- fit
  class(native) <- setdiff(class(native), "mfp2")
  x <- fit$geese$X[rows, , drop = FALSE]
  v <- stats::vcov(native)
  eta <- drop(x %*% stats::coef(native)) + offset
  se_link <- sqrt(rowSums((x %*% v) * x))
  expect_equal(unname(pred_link$fit), unname(eta))
  expect_equal(unname(pred_link$se.fit), unname(se_link))
  expect_equal(unname(pred_resp$fit), unname(exp(eta)))
  expect_equal(unname(pred_resp$se.fit), unname(se_link * exp(eta)))
}

test_that("GEE se.fit predictions with newdata use a formula offset", {
  skip_if_not_installed("geepack")
  d <- gee_offset_prediction_data()
  fit <- suppressWarnings(mfp2(
    cnt ~ fp(x1, df = 2) + x2 + offset(log(expo)), data = d,
    family = gee_family(stats::poisson(), corstr = "exchangeable"),
    id = id, verbose = FALSE
  ))
  rows <- c(2, 7, 30)
  nd <- d[rows, ]
  link <- predict(fit, newdata = nd, se.fit = TRUE)
  resp <- predict(fit, newdata = nd, type = "response", se.fit = TRUE)
  expect_gee_prediction_rows(fit, link, resp, rows, log(d$expo[rows]))
  # Point predictions take the same offset.
  expect_equal(unname(predict(fit, newdata = nd)), unname(link$fit))
})

test_that("GEE se.fit predictions with newdata use newoffset", {
  skip_if_not_installed("geepack")
  d <- gee_offset_prediction_data()
  x <- as.matrix(d[, c("x1", "x2")])
  fit <- suppressWarnings(mfp2(
    x, d$cnt, offset = log(d$expo),
    family = gee_family(stats::poisson(), corstr = "exchangeable"),
    id = d$id, df = c(x1 = 2, x2 = 1), verbose = FALSE
  ))
  rows <- c(1, 4, 9)
  link <- predict(fit, newdata = x[rows, ], newoffset = log(d$expo[rows]),
                  se.fit = TRUE)
  resp <- predict(fit, newdata = x[rows, ], newoffset = log(d$expo[rows]),
                  type = "response", se.fit = TRUE)
  expect_gee_prediction_rows(fit, link, resp, rows, log(d$expo[rows]))
})

test_that("GEE se.fit predictions keep newdata row names", {
  skip_if_not_installed("geepack")
  d <- gee_offset_prediction_data()
  fit <- suppressWarnings(mfp2(
    cnt ~ fp(x1, df = 2) + x2 + offset(log(expo)), data = d,
    family = gee_family(stats::poisson(), corstr = "exchangeable"),
    id = id, verbose = FALSE
  ))
  nd <- d[c(3, 8), ]
  plain <- predict(fit, newdata = nd, type = "response")
  with_se <- predict(fit, newdata = nd, type = "response", se.fit = TRUE)
  expect_identical(names(with_se$fit), names(plain))
  expect_identical(names(with_se$se.fit), names(plain))
})
