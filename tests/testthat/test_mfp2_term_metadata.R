# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: Checks df down-capping rules for binary, ternary, few-level, and
# continuous variables.
test_that("assign_df() correctly limits df for low-cardinality variables", {
  x <- cbind(
    binary = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    ternary = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1),
    few = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
    continuous = 1:10
  )

  df <- assign_df(x, df_default = 4)
  expect_equal(df[["binary"]], 1)    # <= 3 unique -> 1
  expect_equal(df[["ternary"]], 1)   # <= 3 unique -> 1
  expect_equal(df[["few"]], 2)       # 4-5 unique -> min(2, 4) = 2
  expect_equal(df[["continuous"]], 4) # >= 6 unique -> 4
})


# Test purpose: Checks that the selected-variable accessor returns valid predictor
# names.
test_that("get_selected_variable_names() returns correct names", {
  fit <- mfp2(x_prostate, y_prostate, select = 1, verbose = FALSE, warn_low_information = FALSE)

  sel <- get_selected_variable_names(fit)
  expect_true(is.character(sel))
  expect_true(length(sel) > 0)
  expect_true(all(sel %in% colnames(x_prostate)))
})
