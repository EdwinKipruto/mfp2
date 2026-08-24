# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

# Test purpose: A balanced dataset with exactly the same x-response slope in
# every group should not be reported as an interaction. Identical residual
# patterns in both groups make the null interaction deterministic.
test_that("23.8 MFPI does not retain a deterministic no-interaction effect", {
  x_base <- seq(1, 8, length.out = 120)
  residual_pattern <- rep(c(-0.08, 0.08), length.out = length(x_base))
  dat <- rbind(
    data.frame(group = factor("A", levels = c("A", "B")), x = x_base,
               y = 1 + 0.5 * x_base + residual_pattern),
    data.frame(group = factor("B", levels = c("A", "B")), x = x_base,
               y = 1.7 + 0.5 * x_base + residual_pattern)
  )

  fit <- mfpi(
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

  metric <- fit$all_model_metrics[fit$all_model_metrics$variable == "x", , drop = FALSE]
  expect_equal(nrow(metric), 1L)
  expect_gt(metric$pvalue, 0.05)
  expect_false("x" %in% names(fit$best_interaction_model))
})
