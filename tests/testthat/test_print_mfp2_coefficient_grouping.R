# Regression coverage for nonadjacent continuous and zero-indicator columns.
grouping_design_info <- function() {
  data.frame(
    variable = c("x1", "x2", "x2", "x4", "x1", "x3"),
    model_column = paste0("b", 1:6),
    basis = c("positive1", "positive2a", "positive2b", "log4", "zero1", "zero3"),
    center = seq_len(6) / 10,
    stringsAsFactors = FALSE
  )
}

grouping_object <- function() {
  list(
    term_to_columns = stats::setNames(as.list(paste0("x", 1:4)), paste0("x", 1:4)),
    fp_terms = data.frame(row.names = c("x4", "x2", "x1", "x3")),
    centers = 1
  )
}

test_that("print metadata groups terms in input order and retains complete rows", {
  info <- grouping_design_info()
  local_mocked_bindings(mfp2_design_column_info = function(object) info,
                        .package = "mfp2")
  object <- grouping_object()
  before <- object
  expected <- c(1L, 5L, 2L, 3L, 6L, 4L)
  expect_identical(mfp2_print_design_info(object), info[expected, , drop = FALSE])
  expect_identical(object, before)

  object$term_to_columns <- NULL
  expect_identical(mfp2_print_design_info(object)$variable,
                   c("x4", "x2", "x2", "x1", "x1", "x3"))
  object$fp_terms <- NULL
  expect_identical(mfp2_print_design_info(object)$variable,
                   c("x1", "x1", "x2", "x2", "x4", "x3"))
  info <- info[FALSE, , drop = FALSE]
  expect_identical(mfp2_print_design_info(object), info)
})

test_that("ordinary coefficient printing keeps estimates and SEs with grouped bases", {
  info <- grouping_design_info()
  local_mocked_bindings(mfp2_design_column_info = function(object) info,
                        .package = "mfp2")
  set.seed(419)
  dat <- as.data.frame(matrix(stats::rnorm(140), ncol = 7))
  names(dat) <- c("y", paste0("b", 1:6))
  fit <- stats::lm(y ~ ., data = dat)
  fit[names(grouping_object())] <- grouping_object()
  fit$family_string <- "gaussian"
  before <- fit
  output <- capture.output(mfp2_print_coefficients(
    fit, digits = 3L, intercepts = FALSE,
    ctx = list(print_section_heading = function(title) cat(title, "\n"))
  ))
  positions <- vapply(info$basis, function(b) which(grepl(b, output, fixed = TRUE)),
                      integer(1L))
  expect_identical(order(positions), c(1L, 5L, 2L, 3L, 6L, 4L))
  expect_match(trimws(output[positions[1L]]), "^x1\\s")
  expect_match(trimws(output[positions[5L]]), "^zero1\\s")
  expect_match(trimws(output[positions[3L]]), "^positive2b\\s")
  estimates <- stats::coef(fit)[info$model_column]
  ses <- sqrt(diag(stats::vcov(fit)))[info$model_column]
  for (i in seq_len(nrow(info))) {
    values <- strsplit(trimws(output[positions[i]]), "[[:space:]]+")[[1L]]
    expect_equal(as.numeric(tail(values, 3L)),
                 unname(round(c(info$center[i], estimates[i], ses[i]), 3L)))
  }
  expect_identical(fit, before)
})

test_that("multinomial printing groups each logit and resets variable labels", {
  info <- grouping_design_info()
  local_mocked_bindings(mfp2_design_column_info = function(object) info,
                        .package = "mfp2")
  object <- grouping_object()
  object$reference_class <- "A"
  object$mfp2_coefficient_matrix <- matrix(
    seq_len(14), nrow = 2, byrow = TRUE,
    dimnames = list(c("B", "C"), c("(Intercept)", info$model_column))
  )
  before <- object
  capture.output(display <- mfp2_print_multinomial_coefficients(object))
  expect_identical(display$Variable,
                   rep(c("(Intercept)", "x1", "", "x2", "", "x3", "x4"), 2L))
  expect_identical(display$Basis,
                   rep(c("", info$basis[c(1, 5, 2, 3, 6, 4)]), 2L))
  expect_equal(as.numeric(display$Estimate),
               c(1, 2, 6, 3, 4, 7, 5, 8, 9, 13, 10, 11, 14, 12))
  expect_equal(as.numeric(display$Center[display$Center != ""]),
               rep(c(0.1, 0.5, 0.2, 0.3, 0.6, 0.4), 2L))
  expect_identical(object, before)
})
