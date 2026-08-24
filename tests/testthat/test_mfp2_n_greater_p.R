test_that("matrix and formula interfaces accept fits when n is greater than p", {
  dat <- data.frame(
    y = c(1.0, 2.2, 2.7, 4.1),
    x = c(1.0, 2.0, 3.0, 4.0)
  )

  expect_no_error(
    matrix_fit <- mfp2(
      x = as.matrix(dat["x"]),
      y = dat$y,
      family = "gaussian",
      df = 1,
      xorder = "original",
      verbose = FALSE
    )
  )
  expect_s3_class(matrix_fit, "mfp2")

  expect_no_error(
    formula_fit <- mfp2(
      y ~ fp(x, df = 1),
      data = dat,
      family = "gaussian",
      xorder = "original",
      verbose = FALSE
    )
  )
  expect_s3_class(formula_fit, "mfp2")
})


test_that("matrix and formula interfaces reject fits when n is not greater than p", {
  x <- matrix(
    c(1, 0, 0, 1),
    nrow = 2L,
    dimnames = list(NULL, c("x1", "x2"))
  )
  y <- c(1, 2)

  expect_error(
    mfp2(x, y, df = 1, xorder = "original", verbose = FALSE),
    "n = 2, p = 2",
    fixed = TRUE
  )

  dat <- data.frame(y = y, x1 = x[, 1L], x2 = x[, 2L])
  expect_error(
    mfp2(
      y ~ fp(x1, df = 1) + fp(x2, df = 1),
      data = dat,
      xorder = "original",
      verbose = FALSE
    ),
    "n = 2, p = 2",
    fixed = TRUE
  )
})


test_that("the n greater than p check uses the fitted matrix subset", {
  x <- cbind(
    x1 = seq_len(5),
    x2 = c(1, 0, 1, 0, 1)
  )
  y <- seq_len(5)

  expect_error(
    mfp2(
      x,
      y,
      df = 1,
      subset = 1:2,
      xorder = "original",
      verbose = FALSE
    ),
    "n = 2, p = 2",
    fixed = TRUE
  )
})
