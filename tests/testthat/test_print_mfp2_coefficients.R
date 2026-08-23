data("prostate", package = "mfp2")
# These tests target coefficient extraction and printing. Low-information
# behavior is tested separately, so affected fixture fits opt out explicitly;
# other warning types remain visible to testthat.
x_prostate <- as.matrix(prostate[, 2:8])
y_prostate <- as.numeric(prostate$lpsa)

extract_print_section <- function(output, heading, next_heading) {
  normalized <- trimws(output)
  start <- match(heading, normalized)
  if (is.na(start)) {
    stop(sprintf("Print output is missing the '%s' heading.", heading), call. = FALSE)
  }

  following <- if (start < length(output)) {
    seq.int(start + 1L, length(output))
  } else {
    integer(0L)
  }
  next_relative <- match(next_heading, normalized[following])
  if (is.na(next_relative)) {
    stop(sprintf("Print output is missing the '%s' heading.", next_heading), call. = FALSE)
  }

  end <- following[[next_relative]] - 1L
  output[seq.int(start, end)]
}

make_design_info_object <- function(variables,
                                    components,
                                    powers,
                                    centers = NULL,
                                    zero_handled = NULL,
                                    centered = NULL,
                                    transformations = NULL,
                                    term_to_columns = NULL,
                                    formula_factor_info = NULL) {
  transformed_columns <- paste0("internal_basis_", seq_along(variables))
  model_columns <- paste0("internal_coef_", seq_along(variables))

  if (is.null(zero_handled)) {
    zero_handled <- rep(FALSE, length(variables))
  }
  if (is.null(centered)) {
    centered <- rep(!is.null(centers), length(variables))
  }
  if (is.null(term_to_columns)) {
    term_to_columns <- stats::setNames(
      lapply(unique(variables), identity),
      unique(variables)
    )
  }

  fp_rows <- unique(unname(unlist(term_to_columns, use.names = FALSE)))
  fp_rows <- intersect(fp_rows, names(powers))
  fp_terms <- do.call(rbind, lapply(fp_rows, function(v) {
    p <- powers[[v]]
    length(p) <- max(2L, length(p))
    data.frame(
      selected = TRUE,
      df_final = sum(!is.na(p)),
      power1 = p[[1L]],
      power2 = p[[2L]],
      acd = FALSE,
      zero = FALSE,
      catzero = FALSE,
      spike = FALSE,
      row.names = v,
      check.names = FALSE
    )
  }))

  if (is.null(transformations)) {
    transformations <- data.frame(
      shift = rep(0, nrow(fp_terms)),
      scale = rep(1, nrow(fp_terms)),
      center = rep(!is.null(centers), nrow(fp_terms)),
      row.names = rownames(fp_terms)
    )
  }

  list(
    fp_terms = fp_terms,
    transformations = transformations,
    term_to_columns = term_to_columns,
    transformed_to_model_columns = stats::setNames(
      model_columns,
      transformed_columns
    ),
    transformed_column_to_source = stats::setNames(
      variables,
      transformed_columns
    ),
    transformed_column_component = stats::setNames(
      components,
      transformed_columns
    ),
    transformed_column_zero_handled = stats::setNames(
      zero_handled,
      transformed_columns
    ),
    transformed_column_centered = stats::setNames(
      centered,
      transformed_columns
    ),
    centers = if (is.null(centers)) NULL else stats::setNames(
      centers,
      transformed_columns
    ),
    formula_factor_info = formula_factor_info
  )
}


test_that("ordinary FP basis labels use stored mappings and shifts", {
  object <- make_design_info_object(
    variables = c("age", "age"),
    components = c("fp_basis", "fp_basis"),
    powers = list(age = c(3, 3)),
    centers = c(0.0247, -0.0813),
    transformations = data.frame(
      shift = 3,
      scale = 1,
      center = TRUE,
      row.names = "age"
    )
  )

  info <- mfp2_design_column_info(object)

  expect_identical(info$variable, c("age", "age"))
  expect_identical(
    info$transformed_column,
    c("internal_basis_1", "internal_basis_2")
  )
  expect_identical(
    info$model_column,
    c("internal_coef_1", "internal_coef_2")
  )
  expect_identical(
    info$basis,
    c("(age + 3)^3", "(age + 3)^3 * log(age + 3)")
  )
  expect_equal(info$center, c(0.0247, -0.0813))
})


test_that("ACD columns retain separate FP and ACD basis labels", {
  object <- make_design_info_object(
    variables = c("age", "age"),
    components = c("fp_basis", "acd_basis"),
    powers = list(age = c(1, 0.5)),
    centers = c(1.2, 0.4)
  )
  object$fp_terms["age", "acd"] <- TRUE

  info <- mfp2_design_column_info(object)

  expect_identical(info$basis, c("age", "A(age)^0.5"))
})


test_that("numeric binary variables are displayed as the original variable", {
  object <- make_design_info_object(
    variables = "sex",
    components = "identity_binary",
    powers = list(sex = 1),
    centers = 1
  )

  info <- mfp2_design_column_info(object)

  expect_identical(info$variable, "sex")
  expect_identical(info$basis, "sex")
  expect_equal(info$center, 1)
})


test_that("transform metadata marks unchanged numeric binary columns explicitly", {
  x <- matrix(c(1, 2, 1, 2), ncol = 1L, dimnames = list(NULL, "sex"))
  transformed <- transform_matrix(
    x = x,
    power_list = list(sex = 1),
    center = c(sex = TRUE),
    acdx = c(sex = FALSE),
    zero = c(sex = FALSE),
    catzero = c(sex = FALSE),
    spike = c(sex = FALSE)
  )

  expect_identical(
    unname(transformed$transformed_column_component),
    "identity_binary"
  )
  expect_true(unname(transformed$transformed_column_centered))
})


test_that("zero-handled FP components display their positive-part basis", {
  object <- make_design_info_object(
    variables = "age",
    components = "fp_basis",
    powers = list(age = 0.5),
    centers = 1.83,
    zero_handled = TRUE,
    transformations = data.frame(
      shift = 0,
      scale = 1,
      center = TRUE,
      row.names = "age"
    )
  )

  info <- mfp2_design_column_info(object)

  expect_identical(info$basis, "I(age > 0) * age^0.5")
  expect_true(info$zero_handled)
  expect_equal(info$center, 1.83)
})


test_that("catzero and spike component sets describe the fitted columns", {
  x <- matrix(c(-2, 0, 1, 4), ncol = 1L, dimnames = list(NULL, "age"))
  common <- list(
    x = x,
    power_list = list(age = 0.5),
    center = c(age = TRUE),
    acdx = c(age = FALSE),
    zero = c(age = TRUE),
    reset_zero = FALSE
  )

  zero_fit <- do.call(
    transform_matrix,
    c(common, list(catzero = c(age = FALSE), spike = c(age = FALSE)))
  )
  expect_identical(unname(zero_fit$transformed_column_component), "fp_basis")
  expect_true(unname(zero_fit$transformed_column_zero_handled))

  catzero_fit <- do.call(
    transform_matrix,
    c(common, list(catzero = c(age = TRUE), spike = c(age = FALSE)))
  )
  expect_identical(
    unname(catzero_fit$transformed_column_component),
    c("fp_basis", "zero_indicator")
  )

  spike_1 <- do.call(
    transform_matrix,
    c(common, list(
      catzero = c(age = TRUE),
      spike = c(age = TRUE),
      spike_decision = c(age = saz_decision_codes[["cont_binary"]])
    ))
  )
  expect_identical(
    unname(spike_1$transformed_column_component),
    c("fp_basis", "zero_indicator")
  )

  spike_2 <- do.call(
    transform_matrix,
    c(common, list(
      catzero = c(age = TRUE),
      spike = c(age = TRUE),
      spike_decision = c(age = saz_decision_codes[["continuous_only"]])
    ))
  )
  expect_identical(
    unname(spike_2$transformed_column_component),
    "fp_basis"
  )

  spike_3 <- do.call(
    transform_matrix,
    c(common, list(
      catzero = c(age = TRUE),
      spike = c(age = TRUE),
      spike_decision = c(age = saz_decision_codes[["binary_only"]])
    ))
  )
  expect_identical(
    unname(spike_3$transformed_column_component),
    "zero_indicator"
  )
})


test_that("zero indicator columns display the structural-zero condition", {
  object <- make_design_info_object(
    variables = c("age", "age"),
    components = c("fp_basis", "zero_indicator"),
    powers = list(age = 0.5),
    centers = c(1.83, 0),
    zero_handled = c(TRUE, FALSE),
    centered = c(TRUE, TRUE)
  )

  info <- mfp2_design_column_info(object)

  expect_identical(
    info$basis,
    c("I(age > 0) * age^0.5", "I(age <= 0)")
  )
})


test_that("factor treatment coding is described from design_by_level", {
  term_to_columns <- list(group = c("raw_design_b", "raw_design_c"))
  factor_info <- list(
    group = list(
      variable = "group",
      levels = c("A", "B", "C"),
      ordered = FALSE,
      columns = c("raw_design_b", "raw_design_c"),
      design_by_level = matrix(
        c(
          0, 0,
          1, 0,
          0, 1
        ),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(
          c("A", "B", "C"),
          c("raw_design_b", "raw_design_c")
        )
      )
    )
  )

  object <- make_design_info_object(
    variables = c("raw_design_b", "raw_design_c"),
    components = c("identity_binary", "identity_binary"),
    powers = list(raw_design_b = 1, raw_design_c = 1),
    centers = c(0, 0),
    term_to_columns = term_to_columns,
    formula_factor_info = factor_info
  )
  object$fp_terms <- data.frame(
    selected = TRUE,
    df_final = 2,
    power1 = 1,
    power2 = NA_real_,
    acd = FALSE,
    zero = FALSE,
    catzero = FALSE,
    spike = FALSE,
    row.names = "group",
    check.names = FALSE
  )
  object$transformations <- data.frame(
    shift = 0,
    scale = 1,
    center = TRUE,
    row.names = "group"
  )

  info <- mfp2_design_column_info(object)

  expect_identical(info$variable, c("group", "group"))
  expect_identical(
    info$basis,
    c('I(group = "B")', 'I(group = "C")')
  )
})


test_that("custom factor contrasts report the stored design values", {
  factor_info <- list(
    group = list(
      variable = "group",
      levels = c("A", "B", "C"),
      ordered = FALSE,
      columns = c("contrast_1", "contrast_2"),
      design_by_level = matrix(
        c(
          -1, 0,
          0, -1,
          1, 1
        ),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(
          c("A", "B", "C"),
          c("contrast_1", "contrast_2")
        )
      )
    )
  )

  object <- make_design_info_object(
    variables = c("contrast_1", "contrast_2"),
    components = c("fp_basis", "fp_basis"),
    powers = list(contrast_1 = 1, contrast_2 = 1),
    centers = c(0, 0),
    term_to_columns = list(group = c("contrast_1", "contrast_2")),
    formula_factor_info = factor_info
  )
  object$fp_terms <- data.frame(
    selected = TRUE,
    df_final = 2,
    power1 = 1,
    power2 = NA_real_,
    acd = FALSE,
    zero = FALSE,
    catzero = FALSE,
    spike = FALSE,
    row.names = "group",
    check.names = FALSE
  )
  object$transformations <- data.frame(
    shift = 0,
    scale = 1,
    center = TRUE,
    row.names = "group"
  )

  info <- mfp2_design_column_info(object)

  expect_identical(
    info$basis,
    c(
      'contrast("A"=-1, "B"=0, "C"=1)',
      'contrast("A"=0, "B"=-1, "C"=1)'
    )
  )
})


test_that("print.mfp2 shows basis and center when centering is used", {
  fit <- mfp2(x_prostate, y_prostate, center = TRUE, verbose = FALSE, warn_low_information = FALSE)
  output <- capture.output(print(fit, detailed_settings = FALSE))

  coefficient_lines <- extract_print_section(
    output,
    heading = "Final Model Coefficients",
    next_heading = "Model Fit"
  )
  coefficient_output <- paste(coefficient_lines, collapse = "\n")

  expect_match(coefficient_output, "Estimates are for: Basis - Center", fixed = TRUE)
  expect_match(coefficient_output, "Variable", fixed = TRUE)
  expect_match(coefficient_output, "Basis", fixed = TRUE)
  expect_match(coefficient_output, "Center", fixed = TRUE)
  expect_match(coefficient_output, "Estimate", fixed = TRUE)
  expect_match(coefficient_output, "Std. Error", fixed = TRUE)
  expect_match(coefficient_output, "(Intercept)", fixed = TRUE)
  expect_false(grepl("\\b[A-Za-z][A-Za-z0-9_]*\\.[0-9]+\\b", coefficient_output))
})


test_that("print.mfp2 explains zero-specific centering", {
  set.seed(73103)
  dat <- data.frame(
    y = stats::rnorm(180),
    x = sample(c(0, 0, 0, 1:8), 180, replace = TRUE)
  )

  fit <- mfp2(
    y ~ fp(x, df = 2, select = 1, zero = TRUE),
    data = dat,
    center = TRUE,
    verbose = FALSE
  )
  output <- paste(
    capture.output(print(fit, detailed_settings = FALSE)),
    collapse = "\n"
  )

  expect_match(
    output,
    paste0(
      "Estimates use the displayed Center; for zero-handled terms it is ",
      "subtracted only within the positive component."
    ),
    fixed = TRUE
  )
  expect_match(output, "I(x > 0)", fixed = TRUE)
})


test_that("print.mfp2 omits centering output when centering is not used", {
  fit <- mfp2(x_prostate, y_prostate, center = FALSE, verbose = FALSE, warn_low_information = FALSE)
  output <- capture.output(print(fit, detailed_settings = FALSE))

  coefficient_lines <- extract_print_section(
    output,
    heading = "Final Model Coefficients",
    next_heading = "Model Fit"
  )
  coefficient_output <- paste(coefficient_lines, collapse = "\n")

  expect_match(coefficient_output, "Variable", fixed = TRUE)
  expect_match(coefficient_output, "Basis", fixed = TRUE)
  expect_match(coefficient_output, "Estimate", fixed = TRUE)
  expect_match(coefficient_output, "Std. Error", fixed = TRUE)
  expect_false(grepl("Estimates are for: Basis - Center", coefficient_output, fixed = TRUE))
  expect_false(grepl("Center", coefficient_output, fixed = TRUE))
})



test_that("print.mfp2 standard errors come from the final covariance matrix", {
  fit <- mfp2(x_prostate, y_prostate, center = TRUE, verbose = FALSE, warn_low_information = FALSE)
  covariance <- stats::vcov(fit)
  expected_se <- sqrt(diag(covariance))

  output <- capture.output(
    print(fit, detailed_settings = FALSE, digits = 6)
  )
  coefficient_lines <- extract_print_section(
    output,
    heading = "Final Model Coefficients",
    next_heading = "Model Fit"
  )

  intercept_line <- coefficient_lines[
    grepl("(Intercept)", coefficient_lines, fixed = TRUE)
  ]
  expect_length(intercept_line, 1L)
  expect_match(
    intercept_line,
    format(expected_se[["(Intercept)"]], digits = 6, trim = TRUE),
    fixed = TRUE
  )
})

test_that("coef.mfp2 keeps fitted coefficient names", {
  fit <- mfp2(x_prostate, y_prostate, center = TRUE, verbose = FALSE, warn_low_information = FALSE)
  coefficient_names <- names(stats::coef(fit))

  expect_identical(coefficient_names, names(fit$coefficients))
  expect_true(any(grepl("\\.[0-9]+$", coefficient_names)))
})


test_that("print.mfp2 formats intercept and exact-zero centers distinctly", {
  set.seed(73104)
  n <- 240
  exposure <- stats::rnorm(n, mean = 4, sd = 3)
  dat <- data.frame(
    y = 1.2 * pmax(exposure, 0)^0.5 +
      0.8 * (exposure <= 0) + stats::rnorm(n, sd = 0.5),
    exposure = exposure
  )

  fit <- mfp2(
    y ~ fp(exposure, df = 2, select = 1, catzero = TRUE),
    data = dat,
    center = TRUE,
    verbose = FALSE
  )

  output <- capture.output(print(fit, detailed_settings = FALSE))
  coefficient_lines <- extract_print_section(
    output,
    heading = "Final Model Coefficients",
    next_heading = "Model Fit"
  )
  coefficient_output <- paste(coefficient_lines, collapse = "\n")

  intercept_line <- coefficient_lines[grepl("(Intercept)", coefficient_lines, fixed = TRUE)]
  expect_length(intercept_line, 1L)
  expect_false(grepl("\\bNA\\b", intercept_line))

  indicator_line <- coefficient_lines[
    grepl("I(exposure <= 0)", coefficient_lines, fixed = TRUE)
  ]
  expect_length(indicator_line, 1L)
  expect_match(indicator_line, "I\\(exposure <= 0\\)\\s+0\\s+[-+0-9]", perl = TRUE)
  expect_false(grepl("I\\(exposure <= 0\\)\\s+0\\.0+", indicator_line, perl = TRUE))

  expect_false(grepl("\\bNA\\b", coefficient_output))
})


# An intercept-only selected model has no transformed predictor metadata.
# Printing must still render the intercept rather than treating the empty
# final design as incomplete metadata.
test_that("print.mfp2 handles an intercept-only final model", {
  set.seed(2406)
  n <- 160L
  dat <- data.frame(
    y = rnorm(n),
    x = runif(n, 1, 10)
  )

  fit <- mfp2(
    y ~ x,
    data = dat,
    select = 0,
    verbose = FALSE
  )

  expect_length(get_selected_variable_names(fit), 0L)

  design_info <- mfp2_design_column_info(fit)
  expect_s3_class(design_info, "data.frame")
  expect_equal(nrow(design_info), 0L)
  expect_identical(
    names(design_info),
    c(
      "variable", "transformed_column", "model_column", "basis",
      "center", "centered", "zero_handled", "component"
    )
  )

  output <- capture.output(print(fit))
  coefficient_section <- extract_print_section(
    output,
    "Final Model Coefficients",
    "Model Fit"
  )

  expect_true(any(grepl("\\(Intercept\\)", coefficient_section)))
  expect_false(any(grepl("Final design-column metadata is incomplete", output)))
})

# An empty design map is valid only for a model with no fitted predictor
# coefficients. Do not silently omit predictor coefficients when metadata are
# genuinely incomplete.
test_that("mfp2_design_column_info rejects an empty map with predictor coefficients", {
  object <- list(
    coefficients = c("(Intercept)" = 1, x = 0.5),
    transformed_to_model_columns = stats::setNames(character(0L), character(0L)),
    transformed_column_to_source = NULL,
    transformed_column_component = NULL,
    transformed_column_zero_handled = NULL,
    transformed_column_centered = NULL,
    term_to_columns = list(x = "x")
  )

  expect_error(
    mfp2_design_column_info(object),
    "Final design-column metadata is incomplete"
  )
})
