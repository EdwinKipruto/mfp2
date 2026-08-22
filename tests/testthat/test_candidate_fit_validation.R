test_that("candidate-fit validation accepts finite converged fits", {
  expect_true(invisible(mfp2:::validate_mfp_fit_result(
    logl = -12.5,
    df = 3,
    converged = TRUE,
    family_string = "binomial",
    fast = TRUE
  )))

  expect_true(invisible(mfp2:::validate_mfp_fit_result(
    logl = -8.25,
    df = 4,
    converged = TRUE,
    gaussian_deviance = 15.2,
    require_gaussian_deviance = TRUE,
    family_string = "gaussian",
    fast = TRUE
  )))

  expect_true(invisible(mfp2:::validate_mfp_fit_result(
    logl = -5,
    df = 2,
    converged = 1L,
    family_string = "poisson",
    fast = TRUE
  )))
})


test_that("candidate-fit validation rejects explicit non-convergence", {
  expect_error(
    mfp2:::validate_mfp_fit_result(
      logl = -12.5,
      df = 3,
      converged = FALSE,
      family_string = "binomial",
      fast = TRUE
    ),
    "did not converge"
  )
})


test_that("candidate-fit validation rejects non-finite selection criteria", {
  for (bad_logl in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(
      mfp2:::validate_mfp_fit_result(
        logl = bad_logl,
        df = 2,
        family_string = "poisson",
        fast = TRUE
      ),
      "non-finite log-likelihood"
    )
  }

  for (bad_df in c(NA_real_, NaN, Inf, -Inf, -1)) {
    expect_error(
      mfp2:::validate_mfp_fit_result(
        logl = -10,
        df = bad_df,
        family_string = "poisson",
        fast = TRUE
      ),
      "invalid degrees of freedom"
    )
  }

  expect_error(
    mfp2:::validate_mfp_fit_result(
      logl = -10,
      df = 2,
      gaussian_deviance = NA_real_,
      require_gaussian_deviance = TRUE,
      family_string = "gaussian",
      fast = TRUE
    ),
    "non-finite deviance"
  )
})


test_that("Cox convergence guard converts only iteration exhaustion to an error", {
  expect_error(
    mfp2:::mfp2_with_cox_convergence_guard(
      warning("Ran out of iterations and did not converge"),
      fast = TRUE
    ),
    "did not converge"
  )

  expect_warning(
    value <- mfp2:::mfp2_with_cox_convergence_guard(
      {
        warning("coefficient may be infinite")
        7L
      },
      fast = TRUE
    ),
    "coefficient may be infinite"
  )
  expect_identical(value, 7L)
})


test_that("base GLM candidate non-convergence fails before selection", {
  x <- matrix(seq(-6, 6, length.out = 80), ncol = 1L)
  colnames(x) <- "x"
  y <- as.numeric(x[, 1L] > 0)

  expect_error(
    suppressWarnings(
      mfp2:::fit_model(
        x = x,
        y = y,
        family = stats::binomial(),
        family_string = "binomial",
        control = stats::glm.control(maxit = 1L),
        fast = TRUE,
        fitter = "base"
      )
    ),
    "did not converge"
  )
})
