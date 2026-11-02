# Area 4: interaction analysis with mfpi() for additional families.
#
# Data generators (make_mfpi_*) live in helper-newfamilies.R and produce a
# genuine group x covariate interaction so that the interaction is selected and
# predict()/curve comparisons are meaningful.

# ---- 4a: multinomial + mfpi ------------------------------------------------

test_that("mfpi multinomial returns an mfpi object with interaction metrics", {
  skip_on_cran()
  dm <- make_mfpi_multinomial()
  fit <- mfpi(y ~ g + fp(a), data = dm, group_var = "g",
              interaction_vars = "a", interaction_forms = c(a = "fp1"),
              family = multinomial_family(), control = list(maxit = 300),
              p_interact = 0.1, verbose = FALSE)
  expect_s3_class(fit, "mfpi")
  expect_equal(fit$family_string, "multinomial")
  expect_false(is.null(fit$all_model_metrics))
  expect_true(all(c("pvalue", "df_int") %in% names(fit$best_model_metrics)))
})

test_that("multinomial interaction df equals Q * m (Q non-ref logits, FPm)", {
  skip_on_cran()
  dm <- make_mfpi_multinomial()                 # 3 classes -> Q = 2 logits
  fit <- mfpi(y ~ g + fp(a), data = dm, group_var = "g",
              interaction_vars = "a", interaction_forms = c(a = "fp1"),
              family = multinomial_family(), control = list(maxit = 300),
              p_interact = 0.1, verbose = FALSE)
  # FP1 interaction (m = 1) with Q = 2 -> df_int = 2
  expect_equal(fit$best_model_metrics$df_int, 2L)
})

test_that("multinomial mfpi predict returns per-logit interaction predictions", {
  skip_on_cran()
  dm <- make_mfpi_multinomial()
  fit <- mfpi(y ~ g + fp(a), data = dm, group_var = "g",
              interaction_vars = "a", interaction_forms = c(a = "fp1"),
              family = multinomial_family(), control = list(maxit = 300),
              p_interact = 0.1, verbose = FALSE)
  # predict.mfpi() now has a multinomial branch that consumes the per-logit
  # coefficient matrix and returns one fitted curve / contrast per non-reference
  # logit (Q = 2 here, classes B and C relative to reference A).
  pr <- predict(fit)
  expect_s3_class(pr, "mfpi_prediction")
  expect_true(isTRUE(pr$metadata$multinomial))

  # fitted functions carry an explicit `class` column with both non-ref logits
  expect_true("class" %in% names(pr$functions))
  expect_setequal(unique(pr$functions$class), c("B", "C"))
  # difference table likewise
  expect_true("class" %in% names(pr$differences))
  expect_setequal(unique(pr$differences$class), c("B", "C"))

  # SEs and CI bounds are present, finite, and non-degenerate
  expect_true(all(c("se.fit", "lower", "upper") %in% names(pr$functions)))
  expect_true(all(is.finite(pr$functions$fit)))
  expect_true(all(is.finite(pr$functions$se.fit) & pr$functions$se.fit > 0))
})

# ---- 4b: survreg + mfpi ----------------------------------------------------

test_that("mfpi survreg completes and predicts on the link scale", {
  skip_on_cran()
  dsr <- make_mfpi_survreg()
  fit <- mfpi(Surv(time, status) ~ grp + fp(x1), data = dsr, group_var = "grp",
              interaction_vars = "x1", interaction_forms = c(x1 = "fp1"),
              family = survreg_family(dist = "weibull"),
              p_interact = 0.1, verbose = FALSE)
  expect_s3_class(fit, "mfpi")
  expect_equal(fit$family_string, "survreg")
  expect_error(predict(fit), NA)
})

test_that("mfpi survreg estimates per-stratum scales in both models", {
  skip_on_cran()
  dsr <- make_mfpi_survreg()
  fit <- mfpi(Surv(time, status) ~ grp + fp(x1) + strata(sx), data = dsr,
              group_var = "grp", interaction_vars = "x1",
              interaction_forms = c(x1 = "fp1"),
              family = survreg_family(dist = "weibull"),
              p_interact = 0.1, verbose = FALSE)
  adj <- find_fit(fit$adjustment_model, "survreg")
  int <- find_fit(fit$best_interaction_model, "survreg")
  expect_length(adj$scale, 2L)                  # one scale per stratum
  expect_length(int$scale, 2L)
})
