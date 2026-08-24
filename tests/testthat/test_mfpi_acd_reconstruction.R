# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# 25.2 MFPI adjustment-model ACD power reconstruction

# Test purpose: The canonical fp_powers list must preserve both ACD power
# positions even when the display-oriented fp_terms table shows one power.
test_that("MFPI preserves structural ACD power positions", {
  adjustment_model <- list(
    fp_powers = list(
      age = c(NA_real_, 1),
      cavol = 1
    ),
    fp_terms = data.frame(
      power1 = c(1, 1),
      power2 = c(NA_real_, NA_real_),
      row.names = c("age", "cavol")
    ),
    acd = c(age = TRUE, cavol = FALSE)
  )

  powers <- mfpi_extract_adjustment_powers(
    adjustment_model,
    c("age", "cavol")
  )

  expect_named(powers, c("age", "cavol"))
  expect_length(powers$age, 2L)
  expect_true(is.na(powers$age[[1L]]))
  expect_equal(powers$age[[2L]], 1)
  expect_equal(powers$cavol, 1)
})


# Test purpose: A malformed or legacy ACD object that has lost one structural
# power slot should fail before transform_vector_acd() with a targeted error.
test_that("MFPI rejects malformed one-slot ACD adjustment powers", {
  adjustment_model <- list(
    fp_powers = list(age = 1),
    fp_terms = data.frame(
      power1 = 1,
      power2 = NA_real_,
      row.names = "age"
    ),
    acd = c(age = TRUE)
  )

  expect_error(
    mfpi_extract_adjustment_powers(adjustment_model, "age"),
    "must retain two power positions"
  )
})


# Test purpose: Reproduce the formula-interface failure in which a selected ACD
# adjustment variable is rebuilt while another continuous variable is tested
# for interaction. This previously passed a one-element power vector to
# transform_vector_acd().
test_that("MFPI rebuilds selected ACD adjustment terms during interaction fitting", {
  data("prostate", package = "mfp2")

  fit <- NULL
  expect_error(
    fit <- mfpi(
      lpsa ~ fp(age, acdx = TRUE, select = 1) + svi +
        fp(cavol, select = 1),
      data = prostate,
      cont_vars = "cavol",
      group_var = "svi",
      center = FALSE,
      flex = "flex1",
      winsorize = FALSE,
      include_group_var = TRUE,
      p_adjust_method = "holm",
      criterion = "p",
      show_models = FALSE,
      verbose = FALSE
    ),
    NA
  )

  expect_s3_class(fit, "mfpi")
  expect_length(fit$adjustment_model$fp_powers$age, 2L)
})
