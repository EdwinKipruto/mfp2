# Area 4: interaction analysis with mfpi() for the three families.
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

# ---- 4b: finegray + mfpi ---------------------------------------------------

test_that("mfpi finegray completes and returns an mfpi object", {
  skip_on_cran()
  dfg <- make_mfpi_finegray()
  fit <- mfpi(Surv(obs, ev) ~ grp + fp(z1), data = dfg, group_var = "grp",
              interaction_vars = "z1", interaction_forms = c(z1 = "fp1"),
              family = finegray_family(etype = "relapse"),
              p_interact = 0.1, verbose = FALSE)
  expect_s3_class(fit, "mfpi")
  expect_equal(fit$family_string, "finegray")
  expect_error(predict(fit), NA)               # returns predictions
})

test_that("mfpi finegray ordinary prediction returns cumulative incidence", {
  skip_on_cran()
  dfg <- make_mfpi_finegray(n = 650, seed = 78)
  fit <- mfpi(Surv(obs, ev) ~ grp + fp(z1), data = dfg, group_var = "grp",
              interaction_vars = "z1", interaction_forms = c(z1 = "fp1"),
              family = finegray_family(etype = "relapse"),
              p_interact = 1, verbose = FALSE)
  prediction_time <- stats::median(dfg$obs)
  training <- predict(fit, terms = "z1", model = "all", type = "response",
                      times = prediction_time, se.fit = FALSE)
  expect_equal(nrow(training$predictions), nrow(dfg))
  expect_true(all(training$predictions$fit >= 0 &
                  training$predictions$fit <= 1))

  nd <- dfg[1:8, c("grp", "z1"), drop = FALSE]
  new_prediction <- predict(
    fit, newdata = nd, terms = "z1", model = "all", type = "response",
    times = prediction_time, se.fit = FALSE
  )
  expect_equal(new_prediction$predictions$fit,
               training$predictions$fit[1:8], tolerance = 1e-8)
})

test_that("mfpi finegray propagates strata into the weighted Cox fits", {
  skip_on_cran()
  dfg <- make_mfpi_finegray()
  fit <- mfpi(Surv(obs, ev) ~ grp + fp(z1) + strata(sx), data = dfg,
              group_var = "grp", interaction_vars = "z1",
              interaction_forms = c(z1 = "fp1"),
              family = finegray_family(etype = "relapse"),
              p_interact = 0.1, verbose = FALSE)

  # adjustment model
  adj <- find_fit(fit$adjustment_model, "coxph")
  expect_false(is.null(adj$strata))
  # interaction model's inner weighted Cox fit carries the strata term
  int <- find_fit(fit$best_interaction_model, "coxph")
  expect_false(is.null(int))
  expect_false(is.null(int$strata))
  expect_true(any(grepl("strata\\(", as.character(int$formula))))
})

# ---- 4c: survreg + mfpi ----------------------------------------------------

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
