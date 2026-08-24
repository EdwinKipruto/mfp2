# Direct coverage for exported MFPI prediction and metrics print methods.

make_mfpi_prediction_fixture_v <- function(term = "x") {
  functions <- data.frame(
    term = term,
    x = rep(c(1, 2), 2L),
    group = rep(c("control", "treated"), each = 2L),
    fit = c(0.1, 0.2, 0.3, 0.5),
    se.fit = rep(0.05, 4L),
    lower = c(0.0, 0.1, 0.2, 0.4),
    upper = c(0.2, 0.3, 0.4, 0.6),
    stringsAsFactors = FALSE
  )
  differences <- data.frame(
    term = term,
    x = c(1, 2),
    group = "treated",
    reference = "control",
    fit = c(0.2, 0.3),
    se.fit = c(0.07, 0.08),
    lower = c(0.06, 0.14),
    upper = c(0.34, 0.46),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      term = term,
      type = "both",
      functions = functions,
      differences = differences,
      metadata = list(reference = "control")
    ),
    class = c("mfpi_prediction", "list")
  )
}


test_that("print.mfpi_prediction is readable and returns its object invisibly", {
  prediction <- make_mfpi_prediction_fixture_v("x")

  output <- capture.output(result <- withVisible(print(prediction, digits = 2L)))
  printed <- paste(output, collapse = "\n")

  expect_false(result$visible)
  expect_identical(result$value, prediction)
  expect_match(printed, "MFPI prediction  |  term: x", fixed = TRUE)
  expect_match(printed, "Partial fitted functions:", fixed = TRUE)
  expect_match(printed, "Differences:", fixed = TRUE)
  expect_match(printed, "estimate", fixed = TRUE)
  expect_match(printed, "SE", fixed = TRUE)
  expect_false(grepl("metadata", printed, fixed = TRUE))
  expect_identical(names(prediction$functions),
                   c("term", "x", "group", "fit", "se.fit", "lower", "upper"))
})


test_that("print.mfpi_prediction_list prints every term and returns invisibly", {
  predictions <- structure(
    list(
      x = make_mfpi_prediction_fixture_v("x"),
      z = make_mfpi_prediction_fixture_v("z")
    ),
    class = c("mfpi_prediction_list", "list")
  )

  output <- capture.output(result <- withVisible(print(predictions, digits = 2L)))
  printed <- paste(output, collapse = "\n")

  expect_false(result$visible)
  expect_identical(result$value, predictions)
  expect_match(printed, "term: x", fixed = TRUE)
  expect_match(printed, "term: z", fixed = TRUE)
  expect_equal(length(gregexpr("MFPI prediction", printed, fixed = TRUE)[[1L]]), 2L)
})


test_that("MFPI prediction print methods warn about unused arguments", {
  prediction <- make_mfpi_prediction_fixture_v("x")
  predictions <- structure(
    list(x = prediction),
    class = c("mfpi_prediction_list", "list")
  )

  expect_warning(
    capture.output(print(prediction, unused_argument = TRUE)),
    "Unused arguments in `print.mfpi_prediction(...)`: unused_argument.",
    fixed = TRUE
  )
  expect_warning(
    capture.output(print(predictions, unused_argument = TRUE)),
    "Unused arguments in `print.mfpi_prediction_list(...)`: unused_argument.",
    fixed = TRUE
  )
})


test_that("print.best_model_metrics formats stored FP powers", {
  metrics <- data.frame(
    variable = "x",
    type = "fp2",
    pvalue = 0.0123,
    stringsAsFactors = FALSE
  )
  metrics$fp_powers_main <- I(list(c(1, 2)))
  metrics$fp_powers_int <- I(list(list(control = c(1), treated = c(2))))
  class(metrics) <- c("best_model_metrics", "data.frame")

  output <- capture.output(result <- withVisible(print(metrics)))
  printed <- paste(output, collapse = "\n")

  expect_false(result$visible)
  expect_identical(result$value$fp_powers_main, "(1, 2)")
  expect_identical(result$value$fp_powers_int, "(1), (2)")
  expect_identical(metrics$fp_powers_main[[1L]], c(1, 2))
  expect_match(printed, "(1, 2)", fixed = TRUE)
  expect_match(printed, "(1), (2)", fixed = TRUE)
})
