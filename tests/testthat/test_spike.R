# End-to-end test for the full linear reference when a spike term is retained.

# A retained spike implies catzero and zero. The full linear reference must
# therefore contain the positive-part continuous component plus the binary
# structural-zero indicator.
test_that("retained spike uses positive-part plus binary reference design", {
  set.seed(9101)
  n <- 200L
  exposure <- c(rep(0, 60L), stats::rgamma(140L, shape = 2, rate = 1))
  z <- stats::rnorm(n)
  zero_indicator <- as.integer(exposure <= 0)
  y <- 1.5 * zero_indicator + 0.6 * exposure + 0.3 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- mfp2(
    x,
    y,
    spike_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  )

  reference <- stats::glm(
    y ~ exposure_positive + exposure_bin + z,
    data = data.frame(
      y = y,
      exposure_positive = pmax(exposure, 0),
      exposure_bin = zero_indicator,
      z = z
    ),
    family = stats::gaussian()
  )

  expect_true(fit$spike[["exposure"]])
  expect_true(fit$catzero[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end test for a spike request that fails eligibility with no explicit
# zero or catzero fallback requested by the user.

# Once the spike is reset, the full linear reference must revert to the ordinary
# one-column linear term. No structural-zero binary component may remain.
test_that("rejected spike reverts full reference to ordinary linear term", {
  set.seed(9102)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- c(-2, 0) # 1% structural-zero component: below default 10%
  z <- stats::rnorm(n)
  y <- 0.5 * exposure + 0.2 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  ))

  reference <- stats::glm(
    y ~ exposure + z,
    data = data.frame(y = y, exposure = exposure, z = z),
    family = stats::gaussian()
  )

  expect_false(fit$spike[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
  expect_false(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end test for a rejected spike with an explicit zero fallback.

# If the user requested zero handling independently, spike reset must preserve
# that request. The full reference then uses x+ but no binary indicator.
test_that("rejected spike preserves explicit zero reference design", {
  set.seed(9103)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- c(-2, 0) # too few nonpositive observations for SAZ
  z <- stats::rnorm(n)
  y <- 0.5 * pmax(exposure, 0) + 0.2 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    zero_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  ))

  reference <- stats::glm(
    y ~ exposure_positive + z,
    data = data.frame(
      y = y,
      exposure_positive = pmax(exposure, 0),
      z = z
    ),
    family = stats::gaussian()
  )

  expect_false(fit$spike[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end test for a rejected spike with an explicit catzero fallback.

# Explicit catzero survives spike reset and still implies zero. The full
# reference must therefore retain both x+ and the binary structural-zero term.
test_that("rejected spike preserves explicit catzero binary reference design", {
  set.seed(9104)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- c(-2, 0) # too few nonpositive observations for SAZ
  z <- stats::rnorm(n)
  exposure_bin <- as.integer(exposure <= 0)
  y <- 1.5 * exposure_bin + 0.5 * pmax(exposure, 0) +
    0.2 * z + stats::rnorm(n)
  x <- cbind(exposure = exposure, z = z)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    catzero_vars = "exposure",
    df = 1,
    xorder = "original",
    verbose = FALSE
  ))

  reference <- stats::glm(
    y ~ exposure_positive + exposure_bin + z,
    data = data.frame(
      y = y,
      exposure_positive = pmax(exposure, 0),
      exposure_bin = exposure_bin,
      z = z
    ),
    family = stats::gaussian()
  )

  expect_false(fit$spike[["exposure"]])
  expect_true(fit$catzero[["exposure"]])
  expect_true(fit$zero[["exposure"]])
  expect_equal(fit$linear_deviance, stats::deviance(reference), tolerance = 1e-8)
})

# End-to-end tests for spike-at-zero interaction with generated <term>_bin names.

# A retained spike implies catzero, so it must generate the same <term>_bin
# reference column and therefore reject an existing column with that name.
test_that("retained spike rejects an existing source-term _bin column", {
  set.seed(9105)
  n <- 200L
  exposure <- c(rep(0, 60L), stats::rgamma(140L, shape = 2, rate = 1))
  existing_bin <- rep(c(0, 1), length.out = n)
  x <- cbind(
    exposure = exposure,
    exposure_bin = existing_bin
  )
  y <- 1.2 * as.integer(exposure <= 0) + 0.4 * exposure +
    0.2 * existing_bin + stats::rnorm(n)

  expect_error(
    mfp2(
      x,
      y,
      spike_vars = "exposure",
      xorder = "original",
      verbose = FALSE
    ),
    "Generated catzero indicator name.*exposure_bin"
  )
})

# If spike eligibility fails and the user did not request catzero, no binary
# reference column should be generated. An existing exposure_bin predictor must
# therefore remain valid and must not cause the catzero collision error.
test_that("rejected spike does not reserve a source-term _bin name", {
  set.seed(9106)
  n <- 200L
  exposure <- stats::rgamma(n, shape = 2, rate = 1)
  exposure[1:2] <- c(-2, 0) # below the default SAZ component threshold
  existing_bin <- rep(c(0, 1), length.out = n)
  x <- cbind(
    exposure = exposure,
    exposure_bin = existing_bin
  )
  y <- 0.4 * exposure + 0.2 * existing_bin + stats::rnorm(n)

  fit <- suppressWarnings(mfp2(
    x,
    y,
    spike_vars = "exposure",
    xorder = "original",
    verbose = FALSE
  ))

  expect_s3_class(fit, "mfp2")
  expect_false(fit$spike[["exposure"]])
  expect_false(fit$catzero[["exposure"]])
})
