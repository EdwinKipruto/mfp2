# Structural-zero semantics for ACD and shared FP centering -------------------

test_that("fit_acd() fits the positive conditional distribution in zero mode", {
  x_positive <- c(0.25, 0.5, 1, 2, 4, 8, 16)
  x <- c(0, x_positive[1:3], 0, x_positive[4:7])
  powers <- c(0, 1, 2)

  fit_zero <- fit_acd(x, powers = powers, zero = TRUE)
  fit_positive <- fit_acd(x_positive, powers = powers, zero = FALSE)

  expect_identical(unname(fit_zero$acd[x == 0]), c(0, 0))
  expect_equal(
    unname(fit_zero$acd[x > 0]),
    unname(fit_positive$acd),
    tolerance = 1e-12
  )
  expect_equal(unname(fit_zero$beta0), unname(fit_positive$beta0),
               tolerance = 1e-12)
  expect_equal(unname(fit_zero$beta1), unname(fit_positive$beta1),
               tolerance = 1e-12)
  expect_identical(fit_zero$power, fit_positive$power)
  expect_identical(fit_zero$shift, 0)
})


test_that("apply_acd() preserves exact zeros and transforms only positives", {
  training_x <- c(0, 0, 0.25, 0.5, 1, 2, 4, 8, 16)
  fit <- fit_acd(training_x, powers = c(0, 1, 2), zero = TRUE)
  fit$acd <- NULL
  new_x <- c(0, 0.75, 0, 3)

  applied <- do.call(
    apply_acd,
    utils::modifyList(fit, list(x = new_x, zero = TRUE))
  )
  positive_reference <- do.call(
    apply_acd,
    utils::modifyList(fit, list(x = new_x[new_x > 0], zero = FALSE))
  )

  expect_identical(unname(applied[new_x == 0]), c(0, 0))
  expect_equal(
    unname(applied[new_x > 0]),
    unname(positive_reference),
    tolerance = 1e-12
  )
  expect_identical(
    do.call(
      apply_acd,
      utils::modifyList(fit, list(x = c(0, 0), zero = TRUE))
    ),
    c(0, 0)
  )
})


test_that("fit_acd() validates the positive part required by zero mode", {
  expect_error(
    fit_acd(c(0, 1, -1, 2), zero = TRUE),
    "must be nonnegative"
  )
  expect_error(
    fit_acd(c(0, 0, 1), zero = TRUE),
    "At least two positive values"
  )
  expect_error(
    fit_acd(c(0, 1, 1), zero = TRUE),
    "At least two distinct positive values"
  )
})


test_that("FP of ACD uses the original zero-row mask", {
  x <- c(0, 0.25, 0.5, 1, 2, 4, 8, 0, 16)
  fit <- fit_acd(x, powers = c(0, 1, 2), zero = TRUE)
  fit$acd <- NULL

  transformed <- transform_vector_acd(
    x = x,
    power = c(0, 0),
    acd_parameter = fit,
    name = "x",
    zero = TRUE
  )$acd

  expect_true(all(transformed[x == 0, , drop = FALSE] == 0))
  # x = 1 is positive even though its direct log-FP value equals zero.
  expect_identical(unname(transformed[x == 1, "x.1"]), 0)
  expect_true(is.finite(transformed[x == 1, "A_x.1"]))
})


test_that("all ACD candidate generators preserve structural-zero rows", {
  x <- c(0, 0, 0.25, 0.5, 1, 2, 4, 8, 16)
  powers <- c(0, 1)
  fit <- fit_acd(x, powers = powers, zero = TRUE)
  fit$acd <- NULL

  materialized <- generate_transformations_acd(
    x = x,
    degree = 2,
    powers = powers,
    zero = TRUE,
    acd_parameter = fit
  )
  compact <- generate_transformations_acd_basis(
    x = x,
    degree = 2,
    powers = powers,
    zero = TRUE,
    acd_parameter = fit
  )

  expect_true(all(compact$basis[x == 0, , drop = FALSE] == 0))

  for (i in seq_len(nrow(compact$candidate_map))) {
    compact_i <- materialize_acd_basis_candidate(compact, i)
    materialized_i <- materialized$data[[i]]

    expect_true(all(compact_i[x == 0, , drop = FALSE] == 0))
    expect_true(all(materialized_i[x == 0, , drop = FALSE] == 0))
    expect_equal(unname(compact_i), unname(materialized_i), tolerance = 1e-12)
  }
})


test_that("zero-handled FP columns use ordinary full-column centering", {
  x <- matrix(
    c(0, 1, exp(1), exp(2)),
    ncol = 1,
    dimnames = list(NULL, "x")
  )

  transformed <- transform_matrix(
    x = x,
    power_list = list(x = 0),
    center = c(x = TRUE),
    acdx = c(x = FALSE),
    keep_x_order = TRUE,
    zero = c(x = TRUE),
    catzero = c(x = FALSE),
    spike = c(x = FALSE)
  )

  # The zero-padded log basis is c(0, 0, 1, 2). Its complete-sample mean is
  # 0.75, which is subtracted from every row without a post-centering reset.
  expect_equal(
    unname(transformed$x_transformed[, "x.1"]),
    c(-0.75, -0.75, 0.25, 1.25),
    tolerance = 1e-12
  )
  expect_equal(unname(transformed$centers[["x.1"]]), 0.75, tolerance = 1e-12)
  expect_identical(
    unname(transformed$structural_zero_rows[, "x.1"]),
    c(TRUE, FALSE, FALSE, FALSE)
  )
})


test_that("ordinary FP plus spike retains ordinary continuous centering", {
  x <- matrix(
    c(0, 1, exp(1), exp(2)),
    ncol = 1,
    dimnames = list(NULL, "x")
  )

  transformed <- transform_matrix(
    x = x,
    power_list = list(x = 0),
    center = c(x = TRUE),
    acdx = c(x = FALSE),
    keep_x_order = TRUE,
    zero = c(x = TRUE),
    catzero = c(x = TRUE),
    spike = c(x = TRUE),
    spike_decision = c(x = 1L)
  )

  expect_equal(
    unname(transformed$x_transformed[, "x.1"]),
    c(-0.75, -0.75, 0.25, 1.25),
    tolerance = 1e-12
  )
  expect_identical(
    unname(transformed$x_transformed[, "x_bin"]),
    c(1, 0, 0, 0)
  )
  expect_identical(
    unname(transformed$structural_zero_rows[, "x.1"]),
    c(TRUE, FALSE, FALSE, FALSE)
  )
  expect_false(any(transformed$structural_zero_rows[, "x_bin"]))
})


test_that("ACD final centering uses the complete zero-padded basis", {
  x_vector <- c(0, 0, 0.25, 0.5, 1, 2, 4, 8, 16)
  x <- matrix(x_vector, ncol = 1, dimnames = list(NULL, "x"))

  transformed <- transform_matrix(
    x = x,
    power_list = list(x = c(NA_real_, 1)),
    center = c(x = TRUE),
    acdx = c(x = TRUE),
    keep_x_order = TRUE,
    zero = c(x = TRUE),
    catzero = c(x = FALSE),
    spike = c(x = FALSE)
  )

  uncentered <- transformed$x_trafo[["x"]][, "A_x.1"]
  center <- mean(uncentered)

  expect_equal(
    unname(transformed$x_transformed[, "A_x.1"]),
    unname(uncentered - center),
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed$x_transformed[x_vector == 0, "A_x.1"]),
    rep(-center, sum(x_vector == 0)),
    tolerance = 1e-12
  )
  expect_true(all(transformed$structural_zero_rows[x_vector == 0, "A_x.1"]))
  expect_true(transformed$zero_expanded[["A_x.1"]])
})


test_that("ACD plus spike centers the continuous basis and keeps its indicator", {
  x_vector <- c(0, 0, 0.25, 0.5, 1, 2, 4, 8, 16)
  x <- matrix(x_vector, ncol = 1, dimnames = list(NULL, "x"))

  transformed <- transform_matrix(
    x = x,
    power_list = list(x = c(NA_real_, 1)),
    center = c(x = TRUE),
    acdx = c(x = TRUE),
    keep_x_order = TRUE,
    zero = c(x = TRUE),
    catzero = c(x = TRUE),
    spike = c(x = TRUE),
    spike_decision = c(x = 1L)
  )

  uncentered <- transformed$x_trafo[["x"]][, "A_x.1"]
  center <- mean(uncentered)
  expect_equal(
    unname(transformed$x_transformed[x_vector == 0, "A_x.1"]),
    rep(-center, sum(x_vector == 0)),
    tolerance = 1e-12
  )
  expect_identical(
    unname(transformed$x_transformed[, "x_bin"]),
    as.numeric(x_vector == 0)
  )
  expect_true(transformed$zero_expanded[["A_x.1"]])
  expect_false(transformed$zero_expanded[["x_bin"]])
})


test_that("compiled MFPI adjustment construction matches corrected R ACD", {
  x_vector <- c(0, 0, 0.25, 0.5, 1, 2, 4, 8, 16)
  x <- matrix(x_vector, ncol = 1, dimnames = list(NULL, "x"))
  powers <- c(1, 1)
  fit <- fit_acd(x_vector, powers = c(0, 1, 2), zero = TRUE)
  fit$acd <- NULL

  compiled <- build_adjustment_step_loop_cpp(
    x = x,
    x_col_index = 0L,
    vars_adj = "x",
    powers_adj = list(powers),
    acdx_adj = TRUE,
    zero_adj = TRUE,
    catzero_adj = list(NULL),
    spike_adj = FALSE,
    spike_decision_int_adj = NA_integer_,
    acd_parameter_adj = list(fit),
    eliminated = FALSE,
    spike_binary_only_flags = FALSE,
    current_power_keys_adj = list(powers),
    prev_power_keys_adj = list(NULL),
    prev_data_adj_list = list(NULL),
    prev_spike_decision_int_adj = NA_integer_,
    has_prev = FALSE
  )$data_adj_list[["x"]]

  reference <- transform_vector_acd(
    x = x_vector,
    power = powers,
    acd_parameter = fit,
    zero = TRUE
  )$acd

  expect_true(all(compiled[x_vector == 0, , drop = FALSE] == 0))
  expect_equal(unname(compiled), unname(reference), tolerance = 1e-12)
})


test_that("mfp2 prediction reuses a positive-part ACD transformation with zeros", {
  set.seed(20260914)
  x <- c(rep(0, 24), exp(seq(log(0.25), log(16), length.out = 96)))
  z <- as.numeric(scale(log(pmax(x, 0.25))))
  y <- 1.5 + 2 * (x > 0) * stats::pnorm(z) +
    stats::rnorm(length(x), sd = 0.05)

  fit <- mfp2(
    x = matrix(x, ncol = 1, dimnames = list(NULL, "x")),
    y = y,
    acd_vars = "x",
    zero_vars = "x",
    keep = "x",
    force_max_fp_vars = "x",
    select = 1,
    alpha = 1,
    verbose = FALSE
  )

  fit_uncentered <- mfp2(
    x = matrix(x, ncol = 1, dimnames = list(NULL, "x")),
    y = y,
    acd_vars = "x",
    zero_vars = "x",
    keep = "x",
    force_max_fp_vars = "x",
    select = 1,
    alpha = 1,
    center = FALSE,
    verbose = FALSE
  )

  # Centering is only an intercept reparameterisation. With the same selected
  # powers, it must leave the fitted values and likelihood unchanged.
  expect_equal(
    unname(stats::fitted(fit)),
    unname(stats::fitted(fit_uncentered)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(stats::logLik(fit)),
    as.numeric(stats::logLik(fit_uncentered)),
    tolerance = 1e-10
  )

  mixed_newdata <- matrix(
    c(0, 0.3, 1.5, 8),
    ncol = 1,
    dimnames = list(NULL, "x")
  )
  prepared <- as.matrix(prepare_newdata_for_predict(
    fit,
    newdata = mixed_newdata
  ))
  predictions <- predict(fit, newdata = mixed_newdata, type = "link")

  # The uncentered FP/ACD basis is zero at the exact-zero row. After ordinary
  # centering, that row equals minus the stored fitted-sample center.
  expect_equal(
    unname(prepared[1L, , drop = TRUE]),
    unname(-fit$centers[colnames(prepared)]),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(prepared[-1L, , drop = FALSE])))

  # The public prediction must equal the stored transformed design multiplied
  # by the fitted coefficients; no ACD parameters are re-estimated on newdata.
  model_columns <- prediction_model_column_names(fit, colnames(prepared))
  manual <- as.numeric(
    stats::coef(fit)[["(Intercept)"]] +
      prepared %*% stats::coef(fit)[model_columns]
  )
  expect_equal(unname(predictions), manual, tolerance = 1e-12)

  # Applying prediction to the mixed batch or one row at a time must agree.
  rowwise <- vapply(seq_len(nrow(mixed_newdata)), function(i) {
    unname(predict(
      fit,
      newdata = mixed_newdata[i, , drop = FALSE],
      type = "link"
    ))
  }, numeric(1L))

  expect_true(all(is.finite(predictions)))
  expect_equal(unname(predictions), rowwise, tolerance = 1e-12)
})
