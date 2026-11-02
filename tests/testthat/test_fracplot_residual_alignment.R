library(testthat)
library(mfp2)


make_plot_prediction_frames <- function(n = 3L) {
  list(
    age = data.frame(variable = seq_len(n), value = seq_len(n) / 10),
    weight = data.frame(variable = seq_len(n) + 10, value = seq_len(n) / 5)
  )
}


test_that("plot residuals align to fitted observation identifiers", {
  model <- list(
    x_original = matrix(
      seq_len(3L),
      ncol = 1L,
      dimnames = list(c("row_a", "row_b", "row_c"), "x")
    )
  )
  pred_data <- make_plot_prediction_frames()
  pred_data$age <- pred_data$age[c(3, 2, 1), , drop = FALSE]
  rownames(pred_data$age) <- c("row_c", "row_b", "row_a")
  resid <- c(row_c = 30, row_a = 10, row_b = 20)

  result <- mfp2:::mfp2_plot_attach_residuals(pred_data, resid, model)

  # The age frame carries informative row names and has been reordered, while
  # weight has default row names and follows x_original order. Each receives
  # residuals in its own verified observation order.
  expect_equal(result$age$resid, c(30, 20, 10))
  expect_equal(result$weight$resid, c(10, 20, 30))
  expect_null(names(result$age$resid))
})


test_that("plot residual attachment rejects residual-count mismatch", {
  model <- list(x_original = matrix(seq_len(3L), ncol = 1L))

  expect_error(
    mfp2:::mfp2_plot_attach_residuals(
      make_plot_prediction_frames(),
      resid = c(0.1, 0.2),
      model = model
    ),
    "3 fitted observations but 2 residuals",
    fixed = TRUE
  )
})


test_that("plot residual attachment validates every term frame", {
  model <- list(x_original = matrix(seq_len(3L), ncol = 1L))
  pred_data <- make_plot_prediction_frames()
  pred_data$weight <- pred_data$weight[1:2, , drop = FALSE]

  expect_error(
    mfp2:::mfp2_plot_attach_residuals(
      pred_data,
      resid = c(0.1, 0.2, 0.3),
      model = model
    ),
    "term `weight`.*contains 2 but 3 fitted observations"
  )
})


test_that("plot residual attachment rejects incompatible identifiers", {
  model <- list(
    x_original = matrix(
      seq_len(3L),
      ncol = 1L,
      dimnames = list(c("row_a", "row_b", "row_c"), "x")
    )
  )

  expect_error(
    mfp2:::mfp2_plot_attach_residuals(
      make_plot_prediction_frames(),
      resid = c(row_a = 0.1, row_b = 0.2, other = 0.3),
      model = model
    ),
    "Cannot align model residuals with the fitted observations",
    fixed = TRUE
  )
})
