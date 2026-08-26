library(testthat)
library(mfp2)


exact_zero_error <- paste0(
  "require nonnegative covariates(.|[[:space:]])*",
  "Recode negative values explicitly"
)


test_that("matrix mfp2 rejects negatives for zero, catzero, and spike", {
  x <- cbind(exposure = c(-1, 0, seq(0.5, 5, length.out = 38)),
             z = seq_len(40))
  y <- seq_len(nrow(x)) / 10

  for (argument in c("zero_vars", "catzero_vars", "spike_vars")) {
    args <- list(x = x, y = y, verbose = FALSE)
    args[[argument]] <- "exposure"
    expect_error(do.call(mfp2, args), exact_zero_error)
  }
})


test_that("formula mfp2 rejects negatives for zero, catzero, and spike", {
  dat <- data.frame(
    y = seq_len(40) / 10,
    exposure = c(-1, 0, seq(0.5, 5, length.out = 38)),
    z = seq_len(40)
  )

  expect_error(
    mfp2(y ~ fp(exposure, zero = TRUE) + z, data = dat, verbose = FALSE),
    exact_zero_error
  )
  expect_error(
    mfp2(y ~ fp(exposure, catzero = TRUE) + z, data = dat, verbose = FALSE),
    exact_zero_error
  )
  expect_error(
    mfp2(y ~ fp(exposure, spike = TRUE) + z, data = dat, verbose = FALSE),
    exact_zero_error
  )
})


test_that("mfpi matrix and formula interfaces reject negative requested variables", {
  n <- 60L
  dat <- data.frame(
    y = seq_len(n) / 10,
    treatment = rep(c(0, 1), length.out = n),
    age = seq(20, 79),
    exposure = c(-1, rep(0, 15), seq(0.5, 8, length.out = n - 16L))
  )
  x <- as.matrix(dat[c("treatment", "age", "exposure")])

  for (argument in c("zero_vars", "catzero_vars", "spike_vars")) {
    args <- list(
      x = x,
      y = dat$y,
      group_var = "treatment",
      cont_vars = "age",
      flex = "flex1",
      verbose = FALSE
    )
    args[[argument]] <- "exposure"
    expect_error(do.call(mfpi, args), exact_zero_error)
  }

  expect_error(
    mfpi(
      y ~ treatment + fp(age) + fp(exposure, spike = TRUE),
      data = dat,
      group_var = "treatment",
      cont_vars = "age",
      flex = "flex1",
      verbose = FALSE
    ),
    exact_zero_error
  )
})


test_that("one fit error identifies every invalid requested variable", {
  x <- cbind(
    exposure = c(-1, seq(0, 4, length.out = 39)),
    dose = c(-2, seq(0, 8, length.out = 39)),
    z = seq_len(40)
  )
  y <- seq_len(40) / 10

  expect_error(
    mfp2(
      x,
      y,
      zero_vars = "exposure",
      catzero_vars = "dose",
      verbose = FALSE
    ),
    "Affected variable\\(s\\): exposure, dose"
  )
})


test_that("exported FP transformation helper uses exact-zero semantics", {
  expect_equal(
    unname(transform_vector_fp(c(0, 1, 2), power = 1, zero = TRUE)[, 1L]),
    c(0, 1, 2)
  )
})


test_that("exact-zero indicator construction uses equality", {
  x <- matrix(c(0, 0.25, 1, 2), ncol = 1L,
              dimnames = list(NULL, "exposure"))
  transformed <- transform_matrix(
    x = x,
    power_list = list(exposure = 1),
    center = c(exposure = FALSE),
    acdx = c(exposure = FALSE),
    zero = c(exposure = TRUE),
    catzero = c(exposure = TRUE),
    spike = c(exposure = FALSE)
  )$x_transformed

  expect_identical(as.integer(transformed[, "exposure_bin"]), c(1L, 0L, 0L, 0L))
  expect_equal(unname(transformed[, "exposure.1"]), x[, 1L])
})


test_that("mfp2 prediction rejects negatives for every retained zero option", {
  set.seed(8201)
  exposure <- c(rep(0, 50), seq(0.25, 8, length.out = 150))
  z <- stats::rnorm(200)
  y <- 2 * (exposure == 0) + 0.7 * exposure + 0.2 * z + stats::rnorm(200, sd = 0.1)
  x <- cbind(exposure = exposure, z = z)
  bad <- matrix(c(-1, 0), nrow = 1L,
                dimnames = list(NULL, c("exposure", "z")))

  fits <- list(
    zero = mfp2(x, y, zero_vars = "exposure", df = 1, keep = "exposure", verbose = FALSE),
    catzero = mfp2(x, y, catzero_vars = "exposure", df = 1, keep = "exposure", verbose = FALSE),
    spike = mfp2(x, y, spike_vars = "exposure", df = 1, keep = "exposure", verbose = FALSE)
  )

  for (fit in fits) {
    expect_error(predict(fit, newdata = bad), exact_zero_error)
  }
})


test_that("mfp2 contrast references reject negative zero-component values", {
  set.seed(8202)
  exposure <- c(rep(0, 40), seq(0.25, 6, length.out = 120))
  y <- 1.5 * (exposure == 0) + exposure + stats::rnorm(length(exposure), sd = 0.1)
  fit <- mfp2(
    matrix(exposure, ncol = 1L, dimnames = list(NULL, "exposure")),
    y,
    catzero_vars = "exposure",
    df = 1,
    keep = "exposure",
    verbose = FALSE
  )

  expect_error(
    predict(fit, type = "contrasts", terms = "exposure", ref = list(exposure = -1)),
    exact_zero_error
  )
})


test_that("mfpi prediction rejects negatives for binary-only and combined SAZ metadata", {
  newdata <- data.frame(exposure = -1)

  for (decision in c(1L, 3L)) {
    object <- structure(
      list(
        zero_vars = c(exposure = TRUE),
        adjustment_model = list(
          spike_dec = c(exposure = decision)
        )
      ),
      class = "mfpi"
    )

    expect_error(predict(object, newdata = newdata), exact_zero_error)
  }
})


test_that("print and summary label the zero component with x = 0", {
  set.seed(8203)
  exposure <- c(rep(0, 50), seq(0.25, 8, length.out = 150))
  y <- 2 * (exposure == 0) + exposure + stats::rnorm(length(exposure), sd = 0.1)
  fit <- mfp2(
    matrix(exposure, ncol = 1L, dimnames = list(NULL, "exposure")),
    y,
    catzero_vars = "exposure",
    df = 1,
    keep = "exposure",
    center = FALSE,
    verbose = FALSE
  )

  printed <- paste(capture.output(print(fit)), collapse = "\n")
  summarized <- paste(capture.output(print(summary(fit))), collapse = "\n")

  expect_match(printed, "I(exposure = 0)", fixed = TRUE)
  expect_match(summarized, "linear (x > 0) + binary", fixed = TRUE)
  old_zero_label <- paste0("<", "= 0")
  expect_false(grepl(old_zero_label, printed, fixed = TRUE))
  expect_false(grepl(old_zero_label, summarized, fixed = TRUE))
})


test_that("valid nonnegative catzero fitting and prediction match an explicit design", {
  set.seed(8204)
  exposure <- c(rep(0, 40), seq(0.25, 6, length.out = 120))
  indicator <- as.integer(exposure == 0)
  y <- 1 + 1.5 * indicator + 0.8 * exposure + stats::rnorm(length(exposure), sd = 0.1)
  x <- matrix(exposure, ncol = 1L, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    x,
    y,
    catzero_vars = "exposure",
    df = 1,
    keep = "exposure",
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ indicator + exposure)

  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)

  new_exposure <- c(0, 0.5, 2, 5)
  got <- predict(
    fit,
    newdata = matrix(new_exposure, ncol = 1L,
                     dimnames = list(NULL, "exposure")),
    type = "link"
  )
  expected <- stats::predict(
    reference,
    newdata = data.frame(
      indicator = as.integer(new_exposure == 0),
      exposure = new_exposure
    ),
    type = "link"
  )
  expect_equal(unname(got), unname(expected), tolerance = 1e-8)
})


test_that("all-positive zero requests retain ordinary-model behaviour", {
  set.seed(8205)
  x <- cbind(exposure = seq(0.5, 8, length.out = 120), z = stats::rnorm(120))
  y <- 0.7 * x[, "exposure"] + 0.2 * x[, "z"] + stats::rnorm(120)

  ordinary <- mfp2(x, y, df = 1, xorder = "original", verbose = FALSE)
  requested <- suppressWarnings(mfp2(
    x, y, zero_vars = "exposure", df = 1,
    xorder = "original", verbose = FALSE
  ))

  expect_false(requested$zero[["exposure"]])
  expect_equal(stats::coef(requested), stats::coef(ordinary), tolerance = 1e-10)
  expect_equal(stats::fitted(requested), stats::fitted(ordinary), tolerance = 1e-10)
})
