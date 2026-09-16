# Regression tests for MFPI estimability and rank-aware failure boundaries.


make_estimable_mfpi_design <- function() {
  group <- rep(0:1, each = 4)
  x <- c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5)
  group_dummy <- as.numeric(group == 1)
  centered <- x - mean(x)
  z0 <- ifelse(group == 0, x - mean(x[group == 0]), 0)
  z1 <- ifelse(group == 1, x - mean(x[group == 1]), 0)

  list(
    cont_var = matrix(x, ncol = 1L, dimnames = list(NULL, "x")),
    group_var = matrix(group, ncol = 1L, dimnames = list(NULL, "group")),
    xmain = cbind(group1 = group_dummy, x = centered),
    xinteraction = cbind(group1 = group_dummy, x01 = z0, x11 = z1)
  )
}


test_that("MFPI rejects a constant transformed focal block within a group", {
  group <- rep(0:1, each = 3)
  x <- c(2, 2, 2, 1, 2, 3)
  group_dummy <- as.numeric(group == 1)
  pooled <- x - mean(x)
  z0 <- ifelse(group == 0, x - mean(x[group == 0]), 0)
  z1 <- ifelse(group == 1, x - mean(x[group == 1]), 0)

  expect_error(
    validate_mfpi_design_estimability(
      cont_var = matrix(x, ncol = 1L, dimnames = list(NULL, "exposure")),
      group_var = matrix(group, ncol = 1L, dimnames = list(NULL, "arm")),
      xmain = cbind(group1 = group_dummy, exposure = pooled),
      xinteraction = cbind(group1 = group_dummy, exposure01 = z0, exposure11 = z1),
      degree = 0L,
      flex = "flex0"
    ),
    "exposure.*group '0'.*within-group"
  )
})


test_that("MFPI rejects collinear selected FP columns within a group", {
  group <- rep(0:1, each = 4)
  x <- rep(1:4, 2)
  group_dummy <- as.numeric(group == 1)
  z01 <- ifelse(group == 0, x, 0)
  z02 <- z01
  z11 <- ifelse(group == 1, x, 0)
  z12 <- ifelse(group == 1, x^2, 0)

  expect_error(
    validate_mfpi_design_estimability(
      cont_var = matrix(x, ncol = 1L, dimnames = list(NULL, "biomarker")),
      group_var = matrix(group, ncol = 1L, dimnames = list(NULL, "arm")),
      xmain = cbind(group1 = group_dummy, x1 = x, x2 = x^2),
      xinteraction = cbind(
        group1 = group_dummy,
        biomarker01 = z01,
        biomarker02 = z02,
        biomarker11 = z11,
        biomarker12 = z12
      ),
      degree = 2L,
      flex = "flex1"
    ),
    "biomarker.*group '0'.*rank 2 of 3"
  )
})


test_that("MFPI accepts full-rank selected designs", {
  design <- make_estimable_mfpi_design()

  expect_true(validate_mfpi_design_estimability(
    cont_var = design$cont_var,
    group_var = design$group_var,
    xmain = design$xmain,
    xinteraction = design$xinteraction,
    degree = 0L,
    flex = "flex0"
  ))
})


test_that("MFPI design failure occurs before outcome-model fitting", {
  group <- rep(0:1, each = 3)
  x <- c(2, 2, 2, 1, 2, 3)
  group_dummy <- as.numeric(group == 1)
  calls <- 0L

  testthat::local_mocked_bindings(
    fit_model = function(...) {
      calls <<- calls + 1L
      stop("fit_model should not be called", call. = FALSE)
    },
    .package = "mfp2"
  )

  expect_error(
    test_interaction(
      y = seq_len(6),
      cont_var = matrix(x, ncol = 1L, dimnames = list(NULL, "exposure")),
      group_var = matrix(group, ncol = 1L, dimnames = list(NULL, "arm")),
      xmain = cbind(group1 = group_dummy, exposure = x - mean(x)),
      xinteraction = cbind(
        group1 = group_dummy,
        exposure01 = ifelse(group == 0, 0, 0),
        exposure11 = ifelse(group == 1, x - mean(x[group == 1]), 0)
      ),
      degree = 0L,
      bestfp_main = 1,
      bestfp_interaction = list(exposure01 = 1, exposure11 = 1),
      flex = "flex0",
      use_ftest = FALSE,
      family = stats::gaussian(),
      family_string = "gaussian",
      weights = rep(1, 6),
      offset = rep(0, 6),
      ties = "efron",
      strata = NULL,
      control = NULL,
      nocenter = c(-1, 0, 1),
      has_offset = FALSE
    ),
    "exposure.*group '0'"
  )
  expect_equal(calls, 0L)
})


test_that("public mfpi rejects within-group rank loss", {
  group <- rep(0:1, each = 6)
  exposure <- c(rep(2, 6), 1:6)
  x <- cbind(group = group, exposure = exposure)
  y <- c(1.0, 2.1, 1.6, 2.7, 3.1, 2.4, 1.2, 2.0, 2.8, 3.5, 4.1, 5.0)

  expect_error(
    mfpi(
      x = x,
      y = y,
      group_var = "group",
      cont_vars = "exposure",
      cont_var_forms = c(exposure = "linear"),
      family = "gaussian",
      verbose = FALSE
    ),
    "exposure.*group '0'.*within-group"
  )
})


test_that("MFPI rejects family-specific fitted rank loss", {
  design <- make_estimable_mfpi_design()

  expect_error(
    validate_mfpi_fitted_rank(
      fit = list(rank = ncol(design$xmain)),
      design = design$xmain,
      model_name = "main-effects",
      cont_name = "x",
      family_string = "gaussian"
    ),
    "fitted main-effects model is rank deficient"
  )

  expect_true(validate_mfpi_fitted_rank(
    fit = list(rank = ncol(design$xmain) + 1L),
    design = design$xmain,
    model_name = "main-effects",
    cont_name = "x",
    family_string = "gaussian"
  ))
})
