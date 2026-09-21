library(testthat)
library(mfp2)


test_that("survreg censoring metadata respects the Surv response type", {
  right <- list(
    family_string = "survreg",
    y = survival::Surv(c(1, 2, 3, 4), c(1, 0, 1, 0), type = "right")
  )
  right_counts <- mfp2:::mfp2_summary_survreg_censoring(right)
  expect_identical(right_counts$exact, 2L)
  expect_identical(right_counts$right_censored, 2L)
  expect_identical(right_counts$left_censored, 0L)

  left <- list(
    family_string = "survreg",
    y = survival::Surv(c(1, 2, 3, 4), c(1, 0, 1, 0), type = "left")
  )
  left_counts <- mfp2:::mfp2_summary_survreg_censoring(left)
  expect_identical(left_counts$exact, 2L)
  expect_identical(left_counts$left_censored, 2L)
  expect_identical(left_counts$right_censored, 0L)

  interval_y <- structure(
    cbind(
      time1 = c(1, 2, 3, 4),
      time2 = c(1, 2, 4, 5),
      status = c(1, 0, 2, 3)
    ),
    class = "Surv",
    type = "interval"
  )
  interval_counts <- mfp2:::mfp2_summary_survreg_censoring(list(y = interval_y))
  expect_identical(interval_counts$exact, 1L)
  expect_identical(interval_counts$right_censored, 1L)
  expect_identical(interval_counts$left_censored, 1L)
  expect_identical(interval_counts$interval_censored, 1L)
})


test_that("survreg metadata identifies distribution and fixed or stratified scale", {
  y <- survival::Surv(c(1, 2, 3), c(1, 0, 1))
  fixed <- list(
    family_string = "survreg",
    family = survreg_family(dist = "weibull", scale = 0.8),
    dist = "weibull",
    scale = 0.8,
    y = y,
    mfp2_survreg_distribution = "weibull",
    mfp2_survreg_scale_fixed = TRUE
  )
  fixed_metadata <- mfp2:::mfp2_summary_model_metadata(fixed)
  expect_identical(fixed_metadata$distribution, "Weibull")
  expect_equal(fixed_metadata$scale, 0.8)
  expect_true(fixed_metadata$scale_fixed)

  stratified <- fixed
  stratified$family <- survreg_family(dist = "weibull")
  stratified$scale <- c(0.7, 1.1)
  stratified$mfp2_survreg_scale_fixed <- FALSE
  stratified$mfp2_strata_levels <- c("A", "B")
  stratified_metadata <- mfp2:::mfp2_summary_model_metadata(stratified)
  expect_false(stratified_metadata$scale_fixed)
  expect_identical(stratified_metadata$scale_strata$stratum, c("A", "B"))
  expect_equal(stratified_metadata$scale_strata$scale, c(0.7, 1.1))

  stratified_print <- paste(
    capture.output(
      mfp2:::mfp2_print_model_specific_parameters(
        "survreg", stratified_metadata, digits = 3L
      )
    ),
    collapse = "\n"
  )
  expect_match(stratified_print, "Scale parameters (estimated):", fixed = TRUE)
  expect_match(stratified_print, "A", fixed = TRUE)
  expect_match(stratified_print, "B", fixed = TRUE)

  distribution_fixed <- list(
    family_string = "survreg",
    family = survreg_family(dist = "exponential"),
    dist = "exponential",
    scale = 1,
    y = y
  )
  expect_true(
    mfp2:::mfp2_summary_model_metadata(distribution_fixed)$scale_fixed
  )

  printed <- capture.output(
    mfp2:::mfp2_print_model_header(
      family_string = "survreg",
      criterion = "p-value",
      converged = TRUE,
      n = 3L,
      metadata = fixed_metadata,
      digits = 3L
    )
  )
  expect_match(printed[[1L]], "Family: survreg | Distribution: Weibull", fixed = TRUE)
  expect_true(any(grepl("Observations: 3 | Events: 2 | Censored: 1", printed, fixed = TRUE)))
  expect_true(any(grepl("Scale: 0.8 (fixed)", printed, fixed = TRUE)))
})


test_that("negative-binomial metadata prints log link and theta", {
  object <- list(
    family_string = "negbin",
    family = list(link = "log"),
    theta = 2.75
  )
  metadata <- mfp2:::mfp2_summary_model_metadata(object)
  expect_identical(metadata$link, "log")
  expect_equal(metadata$theta, 2.75)

  printed <- capture.output(
    mfp2:::mfp2_print_model_header(
      family_string = "negbin",
      criterion = "AIC",
      converged = TRUE,
      n = 120L,
      metadata = metadata,
      digits = 3L
    )
  )
  expect_match(
    printed[[1L]],
    "Family: negative binomial | Link: log | Criterion: AIC | Converged: yes",
    fixed = TRUE
  )
  expect_true(any(grepl("Theta: 2.75 (estimated)", printed, fixed = TRUE)))
})


test_that("mfp2 summaries expose and print fitted survreg metadata", {
  skip_on_cran()
  data(gbsg, package = "mfp2")
  fit <- mfp2(
    survival::Surv(rectime, censrec) ~ fp(age, df = 1, select = 1),
    data = gbsg,
    family = survreg_family(dist = "weibull"),
    verbose = FALSE
  )
  s <- summary(fit)

  expect_identical(s$distribution, "Weibull")
  expect_equal(s$scale, unname(fit$scale))
  expect_false(s$scale_fixed)
  expect_identical(s$censoring$type, "right")
  expect_identical(s$censoring$exact + s$censoring$right_censored, s$n)

  direct <- paste(capture.output(print(fit)), collapse = "\n")
  summarized <- paste(capture.output(print(s)), collapse = "\n")
  for (text in c(direct, summarized)) {
    expect_match(text, "Distribution: Weibull", fixed = TRUE)
    expect_match(text, "Scale:", fixed = TRUE)
    expect_match(text, "(estimated)", fixed = TRUE)
    expect_match(text, "Censored:", fixed = TRUE)
  }
})


test_that("mfp2 summaries expose and print negative-binomial theta", {
  skip_if_not_installed("fastglm")
  fits <- get_negbin_reference_fits()
  s <- summary(fits$mfp2)

  expect_identical(s$link, "log")
  expect_equal(s$theta, fits$mfp2$theta)

  direct <- paste(capture.output(print(fits$mfp2)), collapse = "\n")
  summarized <- paste(capture.output(print(s)), collapse = "\n")
  for (text in c(direct, summarized)) {
    expect_match(text, "Family: negative binomial", fixed = TRUE)
    expect_match(text, "Link: log", fixed = TRUE)
    expect_match(text, "Theta:", fixed = TRUE)
    expect_false(grepl("Dispersion:", text, fixed = TRUE))
  }
})
