# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: SAZ decision 1 retains both the structural-zero indicator and
# the continuous positive-part effect. Strong independent effects are used so
# AIC has an unambiguous preference for the full two-component representation.
test_that("23.4 SAZ decision 1 retains binary and continuous components", {
  set.seed(2304)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  z <- as.numeric(x == 0)
  y <- 1 + 4.5 * z + 2.2 * x + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ z + x)

  expect_equal(as.integer(fit$spike_dec[["exposure"]]), 1L)
  expect_true(fit$catzero[["exposure"]])
  expect_false(all(is.na(fit$fp_powers[["exposure"]])))
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})


# Test purpose: SAZ decision 2 removes the structural-zero indicator when zero
# observations follow the same continuous relationship as positive values.
test_that("23.5 SAZ decision 2 retains only the continuous component", {
  set.seed(2305)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  y <- 1 + 2.4 * x + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ x)

  expect_equal(as.integer(fit$spike_dec[["exposure"]]), 2L)
  expect_false(fit$catzero[["exposure"]])
  expect_equal(as.numeric(fit$fp_powers[["exposure"]]), 1)
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})


# Test purpose: SAZ decision 3 removes the continuous FP component when only
# membership in the structural-zero group affects the outcome.
test_that("23.6 SAZ decision 3 retains only the binary zero indicator", {
  set.seed(2306)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  z <- as.numeric(x == 0)
  y <- 1 + 4.2 * z + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )
  reference <- stats::glm(y ~ z)

  expect_equal(as.integer(fit$spike_dec[["exposure"]]), 3L)
  expect_true(fit$catzero[["exposure"]])
  expect_true(all(is.na(fit$fp_powers[["exposure"]])))
  expect_equal(unname(stats::fitted(fit)), unname(stats::fitted(reference)),
               tolerance = 1e-8)
})


# Test purpose: Stage-2 SAZ component selection must use alpha rather than
# select. keep forces Stage-1 inclusion (equivalent to select = 1), but the
# binary component should still be removable at alpha = 0.05 when it adds no
# information beyond the positive continuous component. Setting alpha = 1 on
# the same data should retain both components, proving that alpha is the Stage-2
# threshold.
test_that("23.6.1 SAZ Stage 2 uses alpha when select is forced by keep", {
  set.seed(23061)
  positive <- rep(seq(0.5, 8, length.out = 180), each = 2)
  x <- c(rep(0, 140), positive)
  y <- 1 + 2.4 * x + rnorm(length(x), sd = 0.03)
  xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "exposure"))

  fit_alpha_005 <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    keep = "exposure",
    alpha = 0.05,
    criterion = "pvalue",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  fit_alpha_1 <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    keep = "exposure",
    alpha = 1,
    criterion = "pvalue",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  expect_equal(as.integer(fit_alpha_005$spike_dec[["exposure"]]), 2L)
  expect_equal(as.integer(fit_alpha_1$spike_dec[["exposure"]]), 1L)
  expect_true(fit_alpha_005$fp_terms["exposure", "selected"])
  expect_true(fit_alpha_1$fp_terms["exposure", "selected"])
})


# Test purpose: For a retained SAZ model, prediction on zero and positive
# positive new values must match both:
#   1. direct prediction from the final stored GLM; and
#   2. an independent manual calculation using X %*% beta and
#      diag(X %*% vcov(beta) %*% t(X)).
#
# This verifies the complete SAZ prediction contract:
#   - exact-zero values enter the binary structural-zero component;
#   - the continuous component uses the positive part of exposure;
#   - coefficient ordering matches the reconstructed model matrix;
#   - link-scale standard errors use the stored coefficient covariance matrix.
test_that("23.7 SAZ newdata prediction matches stored model and manual matrix calculation", {
  set.seed(2307)

  x <- c(rep(0, 120), seq(0.5, 8, length.out = 300))
  z <- as.numeric(x == 0)

  y <- 1 +
    3.5 * z +
    1.7 * x +
    rnorm(length(x), sd = 0.04)

  xmat <- matrix(
    x,
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  fit <- mfp2(
    xmat,
    y,
    spike_vars = "exposure",
    powers = list(exposure = 1),
    df = 1,
    select = 1,
    criterion = "aic",
    shift = 0,
    scale = 1,
    center = FALSE,
    verbose = FALSE
  )

  # The data-generating model contains both a structural-zero effect and a
  # continuous positive-part effect. Therefore SAZ decision 1 must be retained:
  # binary indicator plus continuous component.
  expect_equal(unname(fit$spike_dec["exposure"]), 1L)

  newx <- matrix(
    c(0, 0.25, 0.5, 2, 6),
    ncol = 1,
    dimnames = list(NULL, "exposure")
  )

  # Public mfp2 prediction.
  got <- predict(
    fit,
    newdata = newx,
    type = "link",
    se.fit = TRUE
  )

  # Because power = 1, shift = 0, scale = 1, and center = FALSE, the final
  # non-intercept design columns are exactly:
  #   exposure_bin = I(exposure == 0)
  #   exposure.1   = exposure
  expected_design <- data.frame(
    exposure_bin = as.numeric(newx[, "exposure"] == 0),
    exposure.1 = newx[, "exposure"],
    check.names = FALSE
  )

  # Remove the mfp2 class so prediction dispatches directly to predict.glm()
  # using the already fitted final model.
  fit_glm <- fit
  class(fit_glm) <- setdiff(class(fit_glm), "mfp2")

  expected <- stats::predict(
    fit_glm,
    newdata = expected_design,
    type = "link",
    se.fit = TRUE
  )

  # Construct the complete model matrix manually, including the intercept.
  manual_x <- cbind(
    `(Intercept)` = 1,
    exposure.1 = expected_design$exposure.1,
    exposure_bin = expected_design$exposure_bin
  )

  beta <- stats::coef(fit_glm)
  beta_vcov <- stats::vcov(fit_glm)

  # The manually constructed design matrix is already in the fitted
  # coefficient order, so the numerical calculation is directly positional.
  expect_true(all(names(beta) %in% colnames(manual_x)))

  expect_equal(ncol(manual_x), length(beta))
  expect_equal(dim(beta_vcov), c(length(beta), length(beta)))

  # Guard against any remaining coefficient-name or ordering mismatch.
  expect_identical(colnames(manual_x), names(beta))
  expect_identical(rownames(beta_vcov), names(beta))
  expect_identical(colnames(beta_vcov), names(beta))

  # Manual link prediction: eta = X beta.
  manual_fit <- as.numeric(manual_x %*% beta)

  # Manual link-scale standard error:
  # se_i = sqrt(x_i' Var(beta) x_i).
  #
  # rowSums((X %*% V) * X) is the diagonal of X V X' without constructing
  # the full prediction covariance matrix.
  manual_variance <- rowSums(
    (manual_x %*% beta_vcov) * manual_x
  )
  manual_se <- as.numeric(sqrt(pmax(manual_variance, 0)))

  # Public mfp2 prediction must equal direct predict.glm().
  expect_equal(
    as.numeric(got$fit),
    as.numeric(expected$fit),
    tolerance = 1e-8
  )

  expect_equal(
    as.numeric(got$se.fit),
    as.numeric(expected$se.fit),
    tolerance = 1e-8
  )

  # Public mfp2 prediction must also equal the independent matrix calculation.
  expect_equal(
    as.numeric(got$fit),
    manual_fit,
    tolerance = 1e-10
  )

  expect_equal(
    as.numeric(got$se.fit),
    manual_se,
    tolerance = 1e-10
  )

  # Direct predict.glm() must agree with the same manual calculation.
  expect_equal(
    as.numeric(expected$fit),
    manual_fit,
    tolerance = 1e-10
  )

  expect_equal(
    as.numeric(expected$se.fit),
    manual_se,
    tolerance = 1e-10
  )
})
