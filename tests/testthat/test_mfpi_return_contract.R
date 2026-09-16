# Regression tests for the public MFPI fitted-model return contract (F-05).


make_f05_contract_data <- function(interaction = FALSE) {
  x_base <- seq(1, 8, length.out = 120)
  residual_pattern <- rep(c(-0.08, 0.08), length.out = length(x_base))
  comparison_slope <- if (interaction) 1.2 else 0.5

  rbind(
    data.frame(
      group = factor("A", levels = c("A", "B")),
      x = x_base,
      y = 1 + 0.5 * x_base + residual_pattern
    ),
    data.frame(
      group = factor("B", levels = c("A", "B")),
      x = x_base,
      y = 1.7 + comparison_slope * x_base + residual_pattern
    )
  )
}


fit_f05_contract_data <- function(interaction = FALSE) {
  dat <- make_f05_contract_data(interaction = interaction)

  mfpi(
    dat[, c("group", "x")],
    dat$y,
    group_var = "group",
    cont_vars = "x",
    cont_var_forms = c(x = "linear"),
    flex = "flex1",
    df = 1,
    select = 1,
    shift = 0,
    scale = 1,
    center = FALSE,
    p_interact = 0.05,
    verbose = FALSE
  )
}


test_that("all_interaction_models stores a direct non-selected model", {
  fit <- fit_f05_contract_data(interaction = FALSE)

  expect_false("x" %in% names(fit$best_interaction_model))
  expect_identical(names(fit$all_interaction_models), "x")

  stored <- fit$all_interaction_models[["x"]]
  canonical <- fit$var_winners[["x"]]$fit$test_results$interaction_model

  expect_identical(stored, canonical)
  expect_true("fit" %in% names(stored))
  expect_false("linear" %in% names(stored))
})


test_that("selected models have the same direct shape in best and all lists", {
  fit <- fit_f05_contract_data(interaction = TRUE)

  expect_true("x" %in% names(fit$best_interaction_model))
  expect_identical(
    fit$all_interaction_models[["x"]],
    fit$best_interaction_model[["x"]]
  )
  expect_identical(
    fit$all_interaction_models[["x"]],
    fit$var_winners[["x"]]$fit$test_results$interaction_model
  )
})


test_that("the flattened all_interaction_models contract survives serialization", {
  fit <- fit_f05_contract_data(interaction = TRUE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  saveRDS(fit, path)
  restored <- readRDS(path)

  restored_model <- restored$all_interaction_models[["x"]]
  original_model <- fit$all_interaction_models[["x"]]

  expect_identical(names(restored$all_interaction_models), "x")
  expect_identical(class(restored_model), class(original_model))
  expect_true("fit" %in% names(restored_model))
  expect_false("linear" %in% names(restored_model))
  expect_equal(
    stats::coef(restored_model$fit),
    stats::coef(original_model$fit),
    tolerance = 0
  )
})
