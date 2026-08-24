# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 2. mfp2.formula() — equivalence and fp() terms
# =============================================================================

# Test purpose: Checks that the formula interface parses fp() terms and returns
# a converged mfp2 object.
test_that("mfp2.formula() returns an mfp2 object", {
  fit <- mfp2(
    lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) + fp(weight) +
      fp(bph) + fp(cp),
    data = prostate, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$convergence_mfp)
  expect_true(is.data.frame(fit$fp_terms))
})


# Test purpose: With all predictors retained and restricted to linear effects,
# verifies the complete NB fit against the independent MASS implementation,
# including covariance, likelihood, AIC, theta, model df, and null likelihood.
test_that("linear negative-binomial mfp2 agrees with MASS::glm.nb()", {
  skip_if_not_installed("fastglm")
  skip_if_not_installed("MASS")

  fits <- get_negbin_reference_fits()
  fit_mfp2 <- fits$mfp2
  fit_mass <- fits$mass

  expect_s3_class(fit_mfp2, "mfp2")
  expect_s3_class(fit_mfp2, "fastglm_nb")
  expect_identical(fit_mfp2$family_string, "negbin")
  expect_identical(fit_mfp2$fitter, "fastglm")
  expect_identical(fit_mfp2$family$link, "log")
  expect_true(all(fit_mfp2$fp_terms[c("x1", "x2", "x3"), "selected"]))
  expect_true(all(vapply(
    fit_mfp2$fp_powers[c("x1", "x2", "x3")],
    function(p) identical(unname(p), 1),
    logical(1L)
  )))

  expect_negbin_mfp2_mass_equal(fit_mfp2, fit_mass)

  # Theta is an estimated nuisance parameter. It contributes to AIC/model df,
  # but not to regression residual degrees of freedom.
  expect_equal(fit_mfp2$mfp_df, fit_mfp2$rank + 1L)
  expect_equal(
    fit_mfp2$aic,
    -2 * fit_mfp2$mfp_logl + 2 * fit_mfp2$mfp_df,
    tolerance = 1e-10
  )
  # With df = 1 and select = 1, the full-linear reference and final MFP model
  # are the same statistical model.
  expect_equal(fit_mfp2$linear_logl, fit_mfp2$mfp_logl, tolerance = 1e-8)
  expect_equal(fit_mfp2$linear_df, fit_mfp2$mfp_df)

  # GLM null_logl is retained as NA for compatibility because Model Fit now
  # reports the already-computed family-specific deviance instead.
  expect_true(is.na(fit_mfp2$null_logl))
  expect_equal(
    fit_mfp2$null_deviance,
    fit_mass$null.deviance,
    tolerance = 1e-3
  )
})


# Test purpose: Compares selected variables from equivalent matrix and formula
# interface fits. Centering is made explicit on both interfaces because fp()
# carries its own center setting (default TRUE) rather than inheriting the
# top-level formula setting.
test_that("default and formula interfaces give consistent selected variables", {
  fit_default <- mfp2(
    x_prostate, y_prostate, center = TRUE, verbose = FALSE
  )
  fit_formula <- mfp2(
    lpsa ~ fp(age, center = TRUE) + fp(svi, df = 1, center = TRUE) +
      fp(pgg45, center = TRUE) + fp(cavol, center = TRUE) +
      fp(weight, center = TRUE) + fp(bph, center = TRUE) +
      fp(cp, center = TRUE),
    data = prostate,
    center = TRUE,
    verbose = FALSE
  )

  sel_default <- sort(get_selected_variable_names(fit_default))
  sel_formula <- sort(get_selected_variable_names(fit_formula))
  expect_equal(sel_default, sel_formula)
})


# Test purpose: Checks that df specified inside fp() is propagated to the fitted
# term metadata.
test_that("fp() applies per-variable df correctly", {
  fit <- mfp2(
    lpsa ~ fp(age, df = 2) + fp(svi, df = 1) + fp(cavol, df = 4),
    data = prostate, verbose = FALSE
  )

  # svi is binary so df should be 1
  expect_equal(as.numeric(fit$fp_terms["svi", "df_initial"]), 1)
})


# Test purpose: Protects formula centering precedence. An fp() term carries its
# own center setting and defaults to TRUE, so it overrides a top-level
# center = FALSE. Ordinary numeric terms continue to use the top-level value.
test_that("fp() center default overrides the global formula center setting", {
  n <- 90L
  i <- seq_len(n)
  dat <- data.frame(
    x1 = seq(1, 9, length.out = n),
    x2 = sin(i / 7) + (i %% 5) / 10
  )
  dat$y <- 0.8 + 0.25 * sqrt(dat$x1) - 0.4 * dat$x2 +
    0.01 * cos(i / 3)

  fit <- mfp2(
    y ~ fp(x1, df = 2, force_max_fp = TRUE) + x2,
    data = dat,
    df = 1,
    keep = c("x1", "x2"),
    select = 1,
    alpha = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    cycles = 2,
    xorder = "original",
    verbose = FALSE
  )

  expect_true(isTRUE(fit$transformations["x1", "center"]))
  expect_false(isTRUE(fit$transformations["x2", "center"]))
})


# Test purpose: Confirms that fp2() can be used as a formula-interface alias for
# fp().
test_that("fp2() is an alias for fp()", {
  fit <- mfp2(
    lpsa ~ fp2(age) + fp2(cavol) + fp2(svi, df = 1),
    data = prostate, verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
})


# Test purpose: Checks that formula-based Cox fitting works with survival responses
# and covariates.
test_that("formula interface works with Cox family and strata", {
  data("gbsg", package = "mfp2")
  fit <- mfp2(
    Surv(rectime, censrec) ~ fp(age) + fp(size) + fp(nodes) +
      fp(er) + meno,
    data = gbsg,
    family = "cox",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_equal(fit$family_string, "cox")
})


# Test purpose: Checks that an unordered factor is represented by one conceptual
# term, while its complete treatment-contrast block is retained through keep.
test_that("formula interface keeps an unordered factor as one grouped term", {
  set.seed(101)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    group = factor(rep(c("A", "B", "C"), length.out = n))
  )

  fit <- mfp2(
    y ~ fp(x) + group,
    data = dat,
    keep = "group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true("group" %in% rownames(fit$fp_terms))
  expect_true(fit$fp_terms["group", "selected"])
  expect_equal(as.numeric(fit$fp_terms["group", "df_setting"]), 1)
  expect_equal(as.numeric(fit$fp_terms["group", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["group", "df_final"]), 2)

  expect_true("group" %in% names(fit$term_to_columns))
  expect_setequal(
    fit$term_to_columns[["group"]],
    c("groupB", "groupC")
  )
  expect_true(
    all(fit$term_to_columns[["group"]] %in% colnames(fit$x_original))
  )
})


# Test purpose: 2.1 Ordered-factor contrast columns are grouped into one fixed
# linear term rather than being treated as separate candidate predictors.
test_that("2.1 Formula interface supports ordered factors as grouped terms", {
  set.seed(2011)
  n <- 120

  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10),
    ordered_group = ordered(
      rep(c("low", "medium", "high"), length.out = n),
      levels = c("low", "medium", "high")
    )
  )

  fit <- mfp2(
    y ~ fp(x) + ordered_group,
    data = dat,
    keep = "ordered_group",
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true("ordered_group" %in% names(fit$term_to_columns))
  expect_length(fit$term_to_columns[["ordered_group"]], 2L)
  expect_true(fit$fp_terms["ordered_group", "selected"])
  expect_equal(as.numeric(fit$fp_terms["ordered_group", "df_setting"]), 1)
  expect_equal(as.numeric(fit$fp_terms["ordered_group", "df_initial"]), 2)
  expect_equal(as.numeric(fit$fp_terms["ordered_group", "df_final"]), 2)

  # Continuous-only extensions are disabled for the complete factor block.
  expect_false(fit$fp_terms["ordered_group", "acd"])
  expect_false(fit$fp_terms["ordered_group", "zero"])
  expect_false(fit$fp_terms["ordered_group", "catzero"])
  expect_false(fit$fp_terms["ordered_group", "spike"])
  expect_true(fit$formula_factor_info[["ordered_group"]]$ordered)
})


# Test purpose: Regression test for final transformation after fitting the ART
# data, which contains both ordinal and nominal three-level predictors.
test_that("formula fit with grouped ART terms completes final transformation", {
  data("art", package = "mfp2")

  # ART stores x4/x9 as numeric level codes. Convert them explicitly so this
  # regression test exercises formula-factor grouping rather than numeric
  # low-cardinality handling.
  art$x4 <- ordered(art$x4)
  art$x9 <- factor(art$x9)

  fit <- mfp2(
    y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = art,
    verbose = FALSE
  )

  expect_s3_class(fit, "mfp2")
  expect_true(fit$convergence_mfp)
  expect_true(is.logical(fit$catzero))
  expect_false(anyNA(fit$catzero))

  expect_true(all(c("x4", "x9") %in% names(fit$term_to_columns)))
  expect_gt(length(fit$term_to_columns[["x4"]]), 1L)
  expect_gt(length(fit$term_to_columns[["x9"]]), 1L)
})


# Test purpose: Ensures formula-interface keep names must match either formula terms
# or expanded model-matrix columns; misspelled names should error.
test_that("formula interface rejects unknown keep variables", {
  expect_error(
    mfp2(
      lpsa ~ fp(age) + fp(cavol),
      data = prostate,
      keep = "does_not_exist",
      verbose = FALSE
    ),
    "Unknown variable"
  )
})


test_that("predict.mfp2 works after fitting with strata()", {
  dat <- survival::lung
  dat$status <- as.integer(dat$status == 2L)
  dat <- dat[stats::complete.cases(dat[, c("time", "status", "age", "sex", "inst")]), ]

  fit <- mfp2(
    survival::Surv(time, status) ~ age + sex + strata(inst),
    data = dat,
    family = "cox",
    df = 1,
    select = 1,
    center = FALSE,
    verbose = FALSE
  )

  nd <- dat[1:5, c("age", "sex", "inst"), drop = FALSE]

  p <- predict(
    fit,
    newdata = nd,
    type = "lp"
  )

  expect_length(p, 5)
  expect_true(all(is.finite(p)))
})
