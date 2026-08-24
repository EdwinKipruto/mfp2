# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 7. ACD transformation
# =============================================================================

# Test purpose: Checks that ACD can be requested through the default matrix
# interface.
test_that("ACD transformation via default interface works", {
  fit <- mfp2(x_prostate, y_prostate, acdx = "cavol", verbose = FALSE, warn_low_information = FALSE)

  expect_s3_class(fit, "mfp2")
  expect_true(fit$fp_terms["cavol", "acd"])
})


# Test purpose: Checks that ACD can be requested inside fp() in the formula
# interface.
test_that("ACD transformation via formula interface works", {
  fit <- mfp2(
    lpsa ~ fp(cavol, acdx = TRUE) + fp(age) + fp(svi, df = 1),
    data = prostate, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$fp_terms["cavol", "acd"])
})


# Test purpose: Checks candidate-power matrix dimensions for ACD degrees 0, 1,
# and 2.
test_that("ACD power generation produces correct matrix dimensions", {
  powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  acd0 <- generate_powers_acd(degree = 0, powers = powx)
  expect_equal(ncol(acd0), 2)
  expect_equal(nrow(acd0), 1)

  acd1 <- generate_powers_acd(degree = 1, powers = powx)
  expect_equal(ncol(acd1), 2)
  expect_equal(nrow(acd1), 8)
  expect_true(all(is.na(acd1[, 1]))) # first column all NA

  acd2 <- generate_powers_acd(degree = 2, powers = powx)
  expect_equal(ncol(acd2), 2)
  expect_equal(nrow(acd2), 64)
})


# Test purpose: Ensures reset_acd() turns off ACD for variables with fewer than
# 5 distinct values while preserving ACD for eligible variables.
test_that("reset_acd() resets variables with fewer than five unique values", {
  x <- cbind(
    low_unique = rep(1:4, each = 10),
    enough_unique = seq_len(40)
  )

  acdx <- c(low_unique = TRUE, enough_unique = TRUE)

  expect_warning(
    out <- reset_acd(x, acdx),
    "fewer than 5 unique values"
  )

  expect_false(out["low_unique"])
  expect_true(out["enough_unique"])
  expect_equal(names(out), names(acdx))
})


# Test purpose: Ensures reset_acd() uses variable names, not vector position,
# when acdx is ordered differently from the columns of x.
test_that("reset_acd() aligns acdx by variable name", {
  x <- cbind(
    low_unique = rep(1:4, each = 10),
    enough_unique = seq_len(40)
  )

  acdx <- c(enough_unique = TRUE, low_unique = TRUE)

  expect_warning(
    out <- reset_acd(x, acdx),
    "low_unique"
  )

  expect_true(out["enough_unique"])
  expect_false(out["low_unique"])
  expect_equal(names(out), names(acdx))
})


# Test purpose: Ensures reset_acd() requires acdx to be a named logical vector,
# because ACD variables are aligned by name.
test_that("reset_acd() rejects unnamed acdx vectors", {
  x <- cbind(x1 = seq_len(20), x2 = seq_len(20))

  expect_error(
    reset_acd(x, c(TRUE, FALSE)),
    "`acdx` must be a named logical vector"
  )
})


# Test purpose: Ensures reset_acd() rejects missing ACD flags before model
# fitting starts.
test_that("reset_acd() rejects missing acdx values", {
  x <- cbind(x1 = seq_len(20), x2 = seq_len(20))
  acdx <- c(x1 = TRUE, x2 = NA)

  expect_error(
    reset_acd(x, acdx),
    "`acdx` must not contain missing values"
  )
})


# Test purpose: Ensures public mfp2() applies reset_acd() and records acd = FALSE
# for requested ACD variables with fewer than 5 unique values.
test_that("mfp2() resets ACD for low-cardinality variables", {
  set.seed(401)
  n <- 120

  x <- cbind(
    low_unique = rep(1:4, length.out = n),
    z = runif(n, 1, 10)
  )
  y <- 0.5 * x[, "low_unique"] + rnorm(n)

  expect_warning(
    fit <- mfp2(
      x,
      y,
      acdx = "low_unique",
      verbose = FALSE
    ),
    "fewer than 5 unique values"
  )

  expect_s3_class(fit, "mfp2")
  expect_false(fit$acd["low_unique"])
  expect_false(fit$fp_terms["low_unique", "acd"])
})


# Test purpose: Ensures retained ACD variables are forced to effective df = 4,
# even when the user supplies a smaller df.
test_that("mfp2() forces retained ACD variables to df = 4", {
  set.seed(402)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    df = c(x1 = 2, z = 2),
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$acd["x1"])
  expect_true(fit$fp_terms["x1", "acd"])
  expect_equal(as.numeric(fit$fp_terms["x1", "df_initial"]), 4)
})


# Test purpose: Ensures formula-interface fp(acdx = TRUE) is translated to ACD
# handling and receives the same effective df = 4 treatment as the default interface.
test_that("formula interface fp(acdx = TRUE) forces effective df = 4", {
  set.seed(403)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)

  dat <- data.frame(
    y = 2 * pnorm(scale(x1)) + 0.2 * z + rnorm(n, sd = 0.2),
    x1 = x1,
    z = z
  )

  fit <- mfp2(
    y ~ fp(x1, acdx = TRUE, df = 2) + fp(z, df = 2),
    data = dat,
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$acd["x1"])
  expect_true(fit$fp_terms["x1", "acd"])
  expect_equal(as.numeric(fit$fp_terms["x1", "df_initial"]), 4)
})


# Test purpose: Ensures the formula interface rejects the old acd_vars argument
# and directs users to fp(..., acdx = TRUE).
test_that("formula interface rejects acd_vars argument", {
  data("prostate", package = "mfp2")

  expect_error(
    mfp2(
      lpsa ~ fp(cavol) + fp(age),
      data = prostate,
      acd_vars = "cavol",
      verbose = FALSE
    ),
    "acd_vars.*not supported"
  )
})


# Test purpose: Ensures retained ACD variables store the fitted ACD parameters
# needed for prediction on newdata.
test_that("mfp2() stores ACD parameters for retained ACD variables", {
  set.seed(404)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_true(fit$fp_terms["x1", "acd"])
  expect_true("x1" %in% names(fit$acd_parameter))
  expect_true(is.list(fit$acd_parameter[["x1"]]))

  expect_true(all(c("beta0", "beta1", "power", "shift", "scale") %in%
                    names(fit$acd_parameter[["x1"]])))
  expect_null(fit$acd_parameter[["x1"]]$acd)
})


# Test purpose: Active ACD variables must bypass ordinary MFP scaling in the
# matrix interface, while non-ACD variables retain their requested scales.
test_that("matrix-interface ACD variables use scale one", {
  set.seed(4041)
  n <- 180
  x1 <- runif(n, 10, 2000)
  z <- runif(n, 1, 20)
  x <- cbind(x1 = x1, z = z)
  y <- 1.5 * pnorm(as.numeric(scale(x1))) + 0.15 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    scale = c(x1 = 1000, z = 10),
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  acd_reference <- fit_acd(x1)

  expect_true(fit$acd[["x1"]])
  expect_equal(unname(fit$transformations["x1", "scale"]), 1)
  expect_equal(unname(fit$transformations["z", "scale"]), 10)
  expect_equal(fit$acd_parameter[["x1"]]$scale, 1)
  expect_equal(fit$acd_parameter[["x1"]]$power, acd_reference$power)
  expect_equal(fit$acd_parameter[["x1"]]$beta0, acd_reference$beta0, tolerance = 1e-10)
  expect_equal(fit$acd_parameter[["x1"]]$beta1, acd_reference$beta1, tolerance = 1e-10)
})


# Test purpose: Formula-level scale settings are also overridden only for terms
# that remain eligible for ACD modelling.
test_that("formula-interface ACD variables use scale one", {
  set.seed(4042)
  n <- 180
  dat <- data.frame(
    x1 = runif(n, 10, 2000),
    z = runif(n, 1, 20)
  )
  dat$y <- 1.5 * pnorm(scale(dat$x1)) + 0.15 * dat$z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    y ~ fp(x1, acdx = TRUE, scale = 1000) + fp(z, scale = 10),
    data = dat,
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  expect_true(fit$acd[["x1"]])
  expect_equal(unname(fit$transformations["x1", "scale"]), 1)
  expect_equal(unname(fit$transformations["z", "scale"]), 10)
  expect_equal(fit$acd_parameter[["x1"]]$scale, 1)
})


# Test purpose: A request reset for insufficient distinct values reverts to
# ordinary FP preprocessing and therefore does not have its scale forced to one.
test_that("reset ACD requests retain ordinary FP scaling", {
  set.seed(4043)
  n <- 160
  x1 <- rep(1:4, length.out = n)
  z <- runif(n, 1, 20)
  x <- cbind(x1 = x1, z = z)
  y <- 0.4 * x1 + 0.1 * z + rnorm(n, sd = 0.2)

  expect_warning(
    fit <- mfp2(
      x,
      y,
      acdx = "x1",
      scale = c(x1 = 100, z = 10),
      df = c(x1 = 2, z = 1),
      select = 1,
      alpha = 1,
      verbose = FALSE
    ),
    "fewer than 5 unique values"
  )

  expect_false(fit$acd[["x1"]])
  expect_equal(unname(fit$transformations["x1", "scale"]), 100)
})


# Test purpose: Ensures predict.mfp2() can reuse stored ACD parameters to
# transform newdata for an ACD-fitted model.
test_that("predict.mfp2() works for ACD models with newdata", {
  set.seed(405)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = x[1:10, , drop = FALSE])

  expect_length(pred, 10)
  expect_true(all(is.finite(pred)))
})


# Test purpose: Ensures predict.mfp2() fails clearly if an active ACD variable
# has no stored ACD parameters.
test_that("predict.mfp2() errors when active ACD parameters are missing", {
  set.seed(406)
  n <- 160

  x1 <- runif(n, 1, 20)
  z <- runif(n, 1, 10)
  x <- cbind(x1 = x1, z = z)

  x1_scaled <- as.numeric(scale(x1))
  y <- 2 * pnorm(x1_scaled) + 0.2 * z + rnorm(n, sd = 0.2)

  fit <- mfp2(
    x,
    y,
    acdx = "x1",
    keep = "x1",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  fit$acd_parameter[["x1"]] <- NULL

  expect_error(
    predict(fit, newdata = x[1:10, , drop = FALSE]),
    "Missing stored ACD parameters"
  )
})


# Test purpose: Checks the exported fit_acd() helper returns an ACD-transformed
# vector on the cumulative-probability scale.
test_that("fit_acd() returns ACD values in the unit interval", {
  set.seed(407)
  x <- runif(100, 1, 20)

  acd <- fit_acd(x)

  expect_true(is.list(acd))
  expect_length(acd$acd, length(x))
  expect_true(all(is.finite(acd$acd)))
  expect_true(all(acd$acd >= 0 & acd$acd <= 1))
  expect_true(all(c("beta0", "beta1", "power", "shift", "scale") %in% names(acd)))
})


# Test purpose: Internal GLM candidate searches can request fitted values without
# retaining the complete backend fit object. The lightweight vector must be
# exactly the one produced by the same fitting backend.
test_that("fit_model() can retain fitted values without retaining the GLM fit", {
  set.seed(701)
  x <- cbind(
    "(Intercept)" = 1,
    "x" = stats::runif(40, 0.5, 3)
  )
  y <- 1.5 + 0.8 * x[, "x"] + stats::rnorm(40, sd = 0.2)
  fam <- stats::gaussian()

  retained <- fit_model(
    x = x,
    y = y,
    family = fam,
    family_string = fam$family,
    fitter = "base",
    x_has_intercept = TRUE,
    keep_fit = TRUE,
    keep_fitted_values = FALSE
  )

  lightweight <- fit_model(
    x = x,
    y = y,
    family = fam,
    family_string = fam$family,
    fitter = "base",
    x_has_intercept = TRUE,
    keep_fit = FALSE,
    keep_fitted_values = TRUE
  )

  expect_false("fit" %in% names(lightweight))
  expect_true("fitted_values" %in% names(lightweight))
  expect_identical(lightweight$fitted_values, retained$fit$fitted.values)
  expect_identical(lightweight$coefficients, retained$coefficients)
  expect_identical(lightweight$logl, retained$logl)
})


# Test purpose: The ACD FP1 search must use the lightweight fit_model() contract
# for every candidate and must not accidentally thread outcome-model weights or
# offsets into the distribution-based ACD auxiliary regression.
test_that("find_best_fp1_for_acd() requests only lightweight fitted values", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$keep_fit <- logical()
  calls$keep_fitted_values <- logical()
  calls$weights <- list()
  calls$offset <- list()

  testthat::local_mocked_bindings(
    # This helper should use the compact FP basis path introduced for mfp2
    # candidate fitting, not recreate the legacy list of n x 1 FP matrices.
    generate_transformations_fp = function(...) {
      stop("legacy FP materialization should not be used", call. = FALSE)
    },
    generate_transformations_fp_basis = function(x, degree, powers, zero, ...) {
      expect_equal(degree, 1L)
      expect_false(zero)
      list(
        basis = cbind(x^-1, log(x), x),
        candidate_map = matrix(1:3, ncol = 1L),
        powers = matrix(powers, ncol = 1L),
        catzero = NULL
      )
    },
    fit_model = function(x,
                         y,
                         weights = NULL,
                         offset = NULL,
                         keep_fit,
                         keep_fitted_values,
                         ...) {
      calls$n <- calls$n + 1L
      i <- calls$n
      calls$keep_fit[[i]] <- keep_fit
      calls$keep_fitted_values[[i]] <- keep_fitted_values
      calls$weights[i] <- list(weights)
      calls$offset[i] <- list(offset)

      logl <- c(-5, -1, -3)[[i]]
      list(
        logl = logl,
        coefficients = c("(Intercept)" = i, "fp1" = i + 1),
        fitted_values = rep.int(as.numeric(i), NROW(x))
      )
    },
    .package = "mfp2"
  )

  out <- find_best_fp1_for_acd(
    x = seq(1, 6),
    y = seq(-1, 1, length.out = 6),
    powers = c(-1, 0, 1),
    zero = FALSE,
    fitter = "base"
  )

  expect_equal(calls$n, 3L)
  expect_true(all(!calls$keep_fit))
  expect_true(all(calls$keep_fitted_values))
  expect_true(all(vapply(calls$weights, is.null, logical(1L))))
  expect_true(all(vapply(calls$offset, is.null, logical(1L))))
  expect_equal(out$power, 0)
  expect_equal(out$coefficients, c("(Intercept)" = 2, "fp1" = 3))
  expect_equal(out$fitted_values, rep(2, 6))
})


# Test purpose: Removing the retained backend object must not change the ACD
# candidate selected, its coefficients, or its fitted values relative to the
# previous keep_fit = TRUE implementation.
test_that("find_best_fp1_for_acd() preserves previous fitted results", {
  x <- seq(0.5, 8, length.out = 60)
  y <- stats::qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
  powers <- c(-1, 0, 0.5, 1, 2)
  fam <- stats::gaussian()
  trafo <- generate_transformations_fp(
    x = x,
    degree = 1L,
    powers = powers,
    zero = FALSE
  )$data

  reference <- lapply(seq_along(powers), function(i) {
    design <- cbind(
      "(Intercept)" = 1,
      "fp1" = trafo[[i]][, 1L]
    )
    fit <- fit_model(
      x = design,
      y = y,
      family = fam,
      family_string = fam$family,
      fitter = "base",
      x_has_intercept = TRUE,
      keep_fit = TRUE
    )
    list(
      deviance = -2 * fit$logl,
      coefficients = fit$coefficients,
      fitted_values = fit$fit$fitted.values
    )
  })

  best <- which.min(vapply(reference, `[[`, numeric(1L), "deviance"))
  out <- find_best_fp1_for_acd(
    x = x,
    y = y,
    powers = powers,
    zero = FALSE,
    fitter = "base"
  )

  expect_identical(out$power, powers[[best]])
  expect_identical(out$coefficients, reference[[best]]$coefficients)
  expect_identical(out$fitted_values, reference[[best]]$fitted_values)
})


# Test purpose: Ensures apply_acd() reproduces the training ACD transformation
# when supplied with parameters from fit_acd().
test_that("apply_acd() reproduces fit_acd() values using stored parameters", {
  set.seed(408)
  x <- runif(100, 1, 20)

  acd <- fit_acd(x)

  applied <- apply_acd(
    x = x,
    beta0 = acd$beta0,
    beta1 = acd$beta1,
    power = acd$power,
    shift = acd$shift,
    scale = acd$scale,
    zero = FALSE
  )

  expect_equal(as.numeric(applied), as.numeric(acd$acd), tolerance = 1e-10)
})


# Test purpose: Ensures generate_powers_acd() only accepts supported ACD degrees
# 0, 1, and 2.
test_that("generate_powers_acd() rejects unsupported degrees", {
  expect_error(
    generate_powers_acd(degree = 3),
    "degree.*ACD.*0, 1, or 2"
  )

  expect_error(
    generate_powers_acd(degree = NA),
    "degree.*ACD.*0, 1, or 2"
  )
})


# Test purpose: Ensures ACD power generation uses ordered Cartesian products,
# because the first power applies to x and the second to A(x).
test_that("generate_powers_acd() keeps ordered power pairs", {
  powx <- c(0, 1)

  acd2 <- generate_powers_acd(degree = 2, powers = powx)

  expect_equal(ncol(acd2), 2)
  expect_equal(nrow(acd2), 4)

  expect_true(any(acd2[, 1] == 0 & acd2[, 2] == 1))
  expect_true(any(acd2[, 1] == 1 & acd2[, 2] == 0))
})


# Test purpose: Ensures generate_transformations_acd() can reuse stored ACD
# parameters instead of refitting them.
test_that("generate_transformations_acd() reuses stored ACD parameters", {
  set.seed(409)
  x <- runif(80, 1, 20)
  powers <- c(0, 1)

  acd_par <- fit_acd(x, powers = powers)
  acd_par$acd <- NULL

  out <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_true(is.list(out))
  expect_equal(nrow(out$powers), 4)
  expect_equal(length(out$data), 4)
  expect_true(all(vapply(out$data, nrow, integer(1)) == length(x)))
})


# Test purpose: Ensures ACD transformation generation includes the catzero
# indicator column when catzero is supplied.
test_that("generate_transformations_acd() includes catzero column when supplied", {
  set.seed(410)
  x <- runif(80, 1, 20)
  catzero <- matrix(as.integer(seq_along(x) <= 10), ncol = 1)

  out <- generate_transformations_acd(
    x = x,
    degree = 1,
    powers = c(0, 1),
    zero = FALSE,
    catzero = catzero
  )

  expect_true(is.list(out))
  expect_equal(length(out$data), 2)

  first <- out$data[[1]]
  expect_equal(nrow(first), length(x))
  expect_equal(colnames(first)[1], "catzero")
})


# Test purpose: Checks that compact ACD degree-2 generation reconstructs every
# materialized candidate while storing only the unique x/A(x) transformations.
test_that("compact ACD basis reconstructs materialized degree-2 candidates", {
  x <- seq(1, 8, length.out = 40)
  allowed_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  # Supply fixed ACD parameters so both representations use exactly the same
  # A(x) values and the test isolates candidate storage/reconstruction.
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  materialized <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = allowed_powers,
    zero = FALSE,
    acd_parameter = acd_par
  )

  expect_equal(compact$powers, materialized$powers)
  expect_equal(dim(compact$candidate_map), c(64L, 2L))
  expect_equal(ncol(compact$basis), 16L)

  for (i in seq_len(nrow(compact$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(compact, i)),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Checks the one-term ACD degree-1 representation and verifies
# that the invariant catzero indicator is stored once rather than per candidate.
test_that("compact ACD basis preserves degree-1 catzero candidates", {
  x <- seq(1, 6, length.out = 30)
  catzero <- matrix(as.integer(seq_along(x) <= 5L), ncol = 1L)
  powers <- c(0, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  materialized <- generate_transformations_acd(
    x = x,
    degree = 1,
    powers = powers,
    zero = FALSE,
    catzero = catzero,
    acd_parameter = acd_par
  )

  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 1,
    powers = powers,
    zero = FALSE,
    catzero = catzero,
    acd_parameter = acd_par
  )

  expect_equal(dim(compact$candidate_map), c(2L, 1L))
  expect_equal(ncol(compact$basis), 2L)
  expect_equal(compact$catzero, catzero)

  for (i in seq_len(nrow(compact$candidate_map))) {
    reconstructed <- materialize_acd_basis_candidate(compact, i)
    expect_equal(
      unname(reconstructed),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
    expect_equal(colnames(reconstructed), c("catzero", "V1"))
  }
})


# Test purpose: Checks that compact ACD generation preserves zero-mode
# semantics for nonpositive x values as well as positive transformed values.
test_that("compact ACD basis preserves zero-mode transformations", {
  x <- c(-2, 0, 1, 2, 4, 8)
  powers <- c(0, 1)
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  materialized <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = TRUE,
    acd_parameter = acd_par
  )

  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = TRUE,
    acd_parameter = acd_par
  )

  for (i in seq_len(nrow(compact$candidate_map))) {
    expect_equal(
      unname(materialize_acd_basis_candidate(compact, i)),
      unname(materialized$data[[i]]),
      tolerance = 1e-12
    )
  }
})


# Test purpose: Confirms the model-search transformation wrapper actually uses
# the compact ACD representation instead of returning a materialized data list.
test_that("transform_data_step() returns compact ACD basis when requested", {
  x_vec <- seq(1, 5, length.out = 25)
  x <- matrix(x_vec, ncol = 1L, dimnames = list(NULL, "x"))
  acd_par <- list(
    beta0 = -1,
    beta1 = 0.25,
    power = 1,
    shift = 0,
    scale = 1
  )

  out <- transform_data_step(
    x = x,
    xi = "x",
    powers_current = list(x = c(1, 1)),
    df = 4,
    powers = list(x = c(0, 1)),
    acdx = c(x = TRUE),
    zero = c(x = FALSE),
    catzero = list(x = NULL),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    acd_parameter = list(x = acd_par),
    prev_adj_params = list(),
    precomputed_adj = list(
      powers_adj = list(),
      spike_decision_adj = numeric(0),
      data_adj_list = list(),
      data_adj = NULL
    ),
    term_to_columns = list(x = 1L),
    compact_acd = TRUE
  )

  expect_null(out$data_fp)
  expect_null(out$fp_basis)
  expect_true(is.list(out$acd_basis))
  expect_equal(dim(out$acd_basis$candidate_map), c(4L, 2L))
  expect_equal(ncol(out$acd_basis$basis), 4L)
})


# Test purpose: Checks important fit_acd() input validation branches.
test_that("fit_acd() validates input arguments", {
  expect_error(
    fit_acd(factor(c("a", "b", "c"))),
    "`x` must be a numeric vector"
  )

  expect_error(
    fit_acd(c(1, NA, 3)),
    "missing values"
  )

  expect_error(
    fit_acd(1),
    "at least two values"
  )

  expect_error(
    fit_acd(1:10, scale = 0),
    "`scale` must be a single positive numeric value"
  )

  expect_error(
    fit_acd(1:10, zero = NA),
    "`zero` must be a single logical value"
  )
})
