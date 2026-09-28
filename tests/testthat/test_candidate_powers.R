# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# =============================================================================
# 5. Candidate-power validation and custom powers
# =============================================================================

# Test purpose: Checks that candidate power vectors are deduplicated and sorted.
test_that("normalize_fp_power_vector() removes duplicates and sorts", {
  result <- normalize_fp_power_vector(c(3, 1, 2, 1, -1))
  expect_equal(result, c(-1, 1, 2, 3))
})


# Test purpose: Checks that invalid power-vector inputs fail validation.
test_that("normalize_fp_power_vector() rejects empty or non-numeric input", {
  expect_error(normalize_fp_power_vector(numeric(0)))
  expect_error(normalize_fp_power_vector("abc"))
  expect_error(normalize_fp_power_vector(c(1, NA, 2)))
})


# Test purpose: Checks that a single non-1 candidate power is allowed and can
# generate repeated-power FP candidates.
test_that("single non-1 candidate power is valid for df > 1", {
  # A single non-unity power should work — repeated-power FP2 candidate
  fit <- mfp2(
    x_prostate[, c("cavol", "age")], y_prostate,
    powers = list(cavol = 2), df = 4,
    verbose = FALSE
  )
  expect_s3_class(fit, "mfp2")
})


# Test purpose: Power 1 remains a valid searched FP1 candidate even when the
# maximum requested model has higher FP degree.
test_that("candidate power set containing only 1 is valid for df > 1", {
  expect_true(
    validate_mfp_candidate_powers(
      powers = list(x = 1),
      df = c(x = 4)
    )
  )
})


# Test purpose: Checks that powers = 1 is allowed when the term is restricted
# to a linear effect.
test_that("candidate power set containing only 1 is valid for df = 1", {
  expect_true(
    validate_mfp_candidate_powers(
      powers = list(x = 1),
      df = c(x = 1)
    )
  )
})


# Test purpose: The numerical power vector c(1) must retain whether it came
# from the fixed-linear row or from the searched FP1 class, because their MFP
# selection degrees of freedom differ.
test_that("searched FP1 power 1 retains its selection identity", {
  metrics <- matrix(0, nrow = 1L, ncol = 1L,
                    dimnames = list("FP1", "deviance_rs"))
  searched <- normalize_selected_step_powers(
    list(power_best = 1, metrics = metrics, model_best = 1L, acd = FALSE),
    xi = "x",
    acdx = c(x = FALSE)
  )

  linear_metrics <- metrics
  rownames(linear_metrics) <- "linear"
  fixed <- normalize_selected_step_powers(
    list(power_best = 1, metrics = linear_metrics,
         model_best = 1L, acd = FALSE),
    xi = "x",
    acdx = c(x = FALSE)
  )

  expect_true(selected_step_is_searched_fp(
    list(metrics = metrics, model_best = 1L)
  ))
  expect_false(selected_step_is_searched_fp(
    list(metrics = linear_metrics, model_best = 1L)
  ))
  expect_null(attr(searched, "mfp_searched_fp", exact = TRUE))
  expect_null(attr(fixed, "mfp_searched_fp", exact = TRUE))
  expect_equal(
    calculate_df(searched, spike_decision = 2, searched_fp = TRUE),
    2
  )
  expect_equal(calculate_df(fixed, spike_decision = 2), 1)
  expect_equal(
    calculate_number_fp_powers(
      list(x = searched),
      searched_fp = c(x = TRUE)
    ),
    1L
  )
  expect_equal(calculate_number_fp_powers(list(x = fixed)), 0L)

  fp_terms <- create_fp_terms(
    fp_powers = list(x = searched),
    acdx = c(x = FALSE),
    df = c(x = 2),
    select = c(x = 0.05),
    alpha = c(x = 0.05),
    criterion = "pvalue",
    zero = c(x = FALSE),
    catzero = c(x = FALSE),
    spike = c(x = FALSE),
    spike_decision = c(x = 2),
    searched_fp = c(x = TRUE)
  )
  expect_true(fp_terms$searched_fp[[1L]])
  expect_equal(fp_terms$df_final[[1L]], 2)
  expect_identical(
    mfp2_summary_form_label(
      powers = 1, acd = FALSE, zero = FALSE, catzero = FALSE,
      spike = FALSE, spike_dec = 2, selected = TRUE,
      searched_fp = TRUE
    ),
    "FP(1)"
  )
})


# Test purpose: Checks that candidate powers supplied inside fp() are accepted
# by the formula interface.
test_that("custom powers via formula fp() work", {
  fit <- mfp2(
    lpsa ~ fp(cavol, powers = c(-1, 0, 1, 2)) + fp(age),
    data = prostate, verbose = FALSE
  )
  expect_s3_class(fit, "mfp2")
})


# Test purpose: Checks FP candidate-generation dimensions, including repeated
# power combinations.
test_that("generate_powers_fp() produces correct number of combinations", {
  powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

  fp1 <- generate_powers_fp(degree = 1, powers = powx)
  expect_equal(nrow(fp1), 8)
  expect_equal(ncol(fp1), 1)

  fp2 <- generate_powers_fp(degree = 2, powers = powx)
  expect_equal(nrow(fp2), 36) # C(8+2-1, 2) = C(9,2) = 36
  expect_equal(ncol(fp2), 2)

  # Single power with degree 2 gives repeated-power pair
  fp_single <- generate_powers_fp(degree = 2, powers = 2)
  expect_equal(nrow(fp_single), 1)
  expect_equal(fp_single[1, ], c(2, 2))
})


# Test purpose: Checks the null-degree FP power matrix used for omitted/null terms.
test_that("generate_powers_fp() degree 0 returns matrix(1)", {
  fp0 <- generate_powers_fp(degree = 0)
  expect_equal(fp0, matrix(1, nrow = 1, ncol = 1))
})
