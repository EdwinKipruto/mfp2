# MFPI coef() and vcov() accessors

# Build a small hand-crafted MFPI object around a standard lm fit. This keeps
# accessor tests focused on MFPI storage, labelling, and S3 behavior rather than
# repeating the full interaction-selection algorithm in every test.
make_mfpi_accessor_fixture <- function(selected = TRUE) {
  set.seed(901)
  n <- 80L
  dat <- data.frame(
    age01 = rnorm(n),
    age02 = rnorm(n),
    age11 = rnorm(n),
    age12 = rnorm(n),
    age21 = rnorm(n),
    age22 = rnorm(n),
    age31 = rnorm(n),
    age32 = rnorm(n),
    trt1 = rbinom(n, 1L, 0.3),
    trt2 = rbinom(n, 1L, 0.3),
    trt3 = rbinom(n, 1L, 0.3),
    hx = rnorm(n)
  )
  dat$y <- with(
    dat,
    0.4 * age01 - 0.2 * age02 + 0.3 * age11 + 0.1 * age12 +
      0.2 * age21 - 0.1 * age22 + 0.15 * age31 + 0.05 * age32 +
      0.25 * trt1 - 0.15 * trt2 + 0.1 * trt3 + 0.35 * hx + rnorm(n)
  )

  fit_lm <- stats::lm(
    y ~ age01 + age02 + age11 + age12 + age21 + age22 + age31 + age32 +
      trt1 + trt2 + trt3 + hx,
    data = dat
  )

  interaction_model <- list(
    fit = fit_lm,
    coefficients = stats::coef(fit_lm),
    transformed_to_model_columns = stats::setNames(
      names(stats::coef(fit_lm)),
      names(stats::coef(fit_lm))
    )
  )

  fit_result <- list(
    bestfp_interaction = list(
      `0` = c(-2, -2),
      `1` = c(0, 1),
      `2` = c(-1, 0.5),
      `3` = c(0.5, 1)
    ),
    coefficient_groups = list(
      `0` = c("age01", "age02"),
      `1` = c("age11", "age12"),
      `2` = c("age21", "age22"),
      `3` = c("age31", "age32")
    ),
    test_results = list(interaction_model = interaction_model)
  )

  out <- list(
    group_var = "trt",
    flex = "flex4",
    shift = c(age = 0),
    cont_var_forms = c(age = "fp2"),
    group_levels_new = 0:3,
    group_levels_original = c("Placebo", "Drug A", "Drug B", "Drug C"),
    group_level_map = data.frame(
      original = c("Placebo", "Drug A", "Drug B", "Drug C"),
      internal = 0:3,
      stringsAsFactors = FALSE
    ),
    var_winners = list(
      age = list(fit = fit_result, type = "fp2")
    ),
    best_interaction_model = if (selected) {
      list(age = interaction_model)
    } else {
      list()
    }
  )
  class(out) <- "mfpi"
  out
}


test_that("coef.mfpi returns readable group-specific coefficients", {
  fit <- make_mfpi_accessor_fixture()
  b <- stats::coef(fit)

  expect_s3_class(b, "mfpi_coef")
  expect_true(is.numeric(b))
  expect_identical(
    names(b),
    c(
      "(Intercept)",
      "Placebo: age^-2",
      "Placebo: age^-2 * log(age)",
      "Drug A: log(age)",
      "Drug A: age",
      "Drug B: age^-1",
      "Drug B: age^0.5",
      "Drug C: age^0.5",
      "Drug C: age",
      "trt [Drug A]",
      "trt [Drug B]",
      "trt [Drug C]",
      "hx"
    )
  )

  raw_names <- attr(b, "raw_names", exact = TRUE)
  expect_identical(
    raw_names,
    names(stats::coef(fit$var_winners$age$fit$test_results$interaction_model$fit))
  )
  expect_equal(
    as.numeric(b),
    as.numeric(stats::coef(fit$var_winners$age$fit$test_results$interaction_model$fit))
  )
})


test_that("coef.mfpi prints powers and separates group-variable coefficients", {
  fit <- make_mfpi_accessor_fixture()
  txt <- capture.output(print(stats::coef(fit), digits = 4))
  txt <- paste(txt, collapse = "\n")

  expect_match(txt, "Interaction: age", fixed = TRUE)
  expect_match(txt, "Group variable: trt", fixed = TRUE)
  expect_match(txt, "FLEX: FLEX4", fixed = TRUE)
  expect_match(txt, "FP form: FP2", fixed = TRUE)
  expect_match(txt, "Interaction powers:", fixed = TRUE)
  expect_match(txt, "Placebo: (-2, -2)", fixed = TRUE)
  expect_match(txt, "Drug A:  (0, 1)", fixed = TRUE)
  expect_match(txt, "Drug B:  (-1, 0.5)", fixed = TRUE)
  expect_match(txt, "Drug C:  (0.5, 1)", fixed = TRUE)
  expect_match(txt, "Group-specific FP terms:", fixed = TRUE)
  expect_match(txt, "Group level", fixed = TRUE)
  expect_match(txt, "age^-2 * log(age)", fixed = TRUE)
  expect_match(txt, "Group-variable coefficients:", fixed = TRUE)
  expect_match(txt, "Reference level: Placebo", fixed = TRUE)
  expect_match(txt, "Drug A", fixed = TRUE)
  expect_false(grepl("vs Placebo", txt, fixed = TRUE))
  expect_match(txt, "Model intercept:", fixed = TRUE)
  expect_match(txt, "Adjustment coefficients:", fixed = TRUE)
  expect_match(txt, "hx", fixed = TRUE)
})


test_that("coef.mfpi prints powers by group even when FLEX uses common powers", {
  fit <- make_mfpi_accessor_fixture()
  fit$flex <- "flex1"
  fit$var_winners$age$fit$bestfp_interaction <- list(
    `0` = c(3, 3),
    `1` = c(3, 3),
    `2` = c(3, 3),
    `3` = c(3, 3)
  )

  txt <- paste(capture.output(print(stats::coef(fit), digits = 4)), collapse = "\n")

  expect_match(txt, "FLEX: FLEX1", fixed = TRUE)
  expect_match(txt, "Interaction powers:", fixed = TRUE)
  expect_match(txt, "Placebo: (3, 3)", fixed = TRUE)
  expect_match(txt, "Drug A:  (3, 3)", fixed = TRUE)
  expect_match(txt, "Drug B:  (3, 3)", fixed = TRUE)
  expect_match(txt, "Drug C:  (3, 3)", fixed = TRUE)
})


test_that("vcov.mfpi matches the underlying model and coef labels", {
  fit <- make_mfpi_accessor_fixture()
  b <- stats::coef(fit, term = "age")
  V <- stats::vcov(fit, term = "age")
  V_raw <- stats::vcov(fit$var_winners$age$fit$test_results$interaction_model$fit)

  expect_true(is.matrix(V))
  expect_identical(rownames(V), names(b))
  expect_identical(colnames(V), names(b))
  expect_true("trt [Drug A]" %in% rownames(V))
  expect_false(any(grepl("vs Placebo", rownames(V), fixed = TRUE)))
  expect_equal(as.numeric(V), as.numeric(V_raw))
  expect_identical(attr(V, "raw_names", exact = TRUE), rownames(V_raw))
  expect_identical(rownames(V_raw), colnames(V_raw))
})


test_that("model = all exposes evaluated but unselected MFPI interactions", {
  fit <- make_mfpi_accessor_fixture(selected = FALSE)

  expect_error(
    stats::coef(fit),
    "No retained MFPI interaction models"
  )
  expect_error(
    stats::vcov(fit),
    "No retained MFPI interaction models"
  )

  b <- stats::coef(fit, term = "age", model = "all")
  V <- stats::vcov(fit, term = "age", model = "all")
  expect_s3_class(b, "mfpi_coef")
  expect_true(is.matrix(V))
})


test_that("MFPI accessors reject unknown terms and unused arguments", {
  fit <- make_mfpi_accessor_fixture()

  expect_error(
    stats::coef(fit, term = "unknown"),
    "do not have best MFPI interaction models"
  )
  expect_error(stats::coef(fit, foo = 1), "Unused argument")
  expect_error(stats::vcov(fit, foo = 1), "Unused argument")
})


test_that("MFPI accessors return named collections for multiple interactions", {
  fit <- make_mfpi_accessor_fixture()
  second <- fit$var_winners$age
  fit$var_winners$hg <- second
  fit$best_interaction_model$hg <- second$fit$test_results$interaction_model
  fit$cont_var_forms <- c(age = "fp2", hg = "fp2")
  fit$shift <- c(age = 0, hg = 0)

  b <- stats::coef(fit)
  V <- stats::vcov(fit)

  expect_s3_class(b, "mfpi_coef_list")
  expect_identical(names(b), c("age", "hg"))
  expect_s3_class(V, "mfpi_vcov_list")
  expect_identical(names(V), c("age", "hg"))

  txt <- paste(capture.output(print(V)), collapse = "\n")
  expect_match(txt, "MFPI variance-covariance matrices", fixed = TRUE)
  expect_match(txt, "age", fixed = TRUE)
  expect_match(txt, "hg", fixed = TRUE)
  expect_match(txt, "vcov(object, term", fixed = TRUE)
})


test_that("summary.mfpi builds readable regression displays", {
  fit <- make_mfpi_accessor_fixture()
  s <- summary(fit)

  expect_s3_class(s, "summary.mfpi")
  expect_true("age" %in% names(s$regression_displays))

  display <- s$regression_displays$age
  expect_identical(display$info$reference_level, "Placebo")
  expect_identical(display$info$group_variable_table$level,
                   c("Drug A", "Drug B", "Drug C"))
  expect_identical(
    vapply(display$info$power_table$powers, mfpi_format_power_vector, character(1L)),
    c("(-2, -2)", "(0, 1)", "(-1, 0.5)", "(0.5, 1)")
  )

  raw_summary <- summary(
    fit$var_winners$age$fit$test_results$interaction_model$fit
  )$coefficients
  stats <- display$statistics
  expect_equal(
    stats$std_error,
    unname(raw_summary[stats$raw_name, "Std. Error"])
  )
  expect_identical(attr(stats, "statistic_label", exact = TRUE), "t")
})


test_that("MFPI regression summary print is readable and omits raw model output", {
  fit <- make_mfpi_accessor_fixture()
  display <- summary(fit)$regression_displays$age

  txt <- paste(
    capture.output(mfpi_print_summary_regression_block(display, digits = 3L)),
    collapse = "\n"
  )

  expect_match(txt, "Interaction: age", fixed = TRUE)
  expect_match(txt, "Interaction: age\n----------------", fixed = TRUE)
  expect_false(grepl("Group variable: trt", txt, fixed = TRUE))
  expect_false(grepl("FLEX: FLEX4", txt, fixed = TRUE))
  expect_match(txt, "FP form: FP2", fixed = TRUE)
  expect_match(txt, "Interaction powers:", fixed = TRUE)
  expect_match(txt, "Placebo: (-2, -2)", fixed = TRUE)
  expect_match(txt, "Drug A:  (0, 1)", fixed = TRUE)
  expect_match(txt, "Drug B:  (-1, 0.5)", fixed = TRUE)
  expect_match(txt, "Drug C:  (0.5, 1)", fixed = TRUE)
  expect_match(txt, "Group-specific FP terms:", fixed = TRUE)
  expect_match(txt, "Std. Error", fixed = TRUE)
  expect_match(txt, "p-value", fixed = TRUE)
  expect_match(txt, "Group-variable coefficients:", fixed = TRUE)
  expect_match(txt, "Reference level: Placebo", fixed = TRUE)
  expect_match(txt, "Adjustment coefficients:", fixed = TRUE)
  expect_match(txt, "Model intercept:", fixed = TRUE)

  expect_false(grepl("age01", txt, fixed = TRUE))
  expect_false(grepl("vs Placebo", txt, fixed = TRUE))
  expect_false(grepl("Call:", txt, fixed = TRUE))
  expect_false(grepl("exp(coef)", txt, fixed = TRUE))
  expect_false(grepl("Signif. codes", txt, fixed = TRUE))
  expect_false(grepl("Concordance", txt, fixed = TRUE))
  expect_false(grepl("Likelihood ratio", txt, fixed = TRUE))
  expect_false(grepl("Wald test", txt, fixed = TRUE))
  expect_false(grepl("Score", txt, fixed = TRUE))
})


test_that("MFPI summary interaction heading underline follows the term name", {
  fit <- make_mfpi_accessor_fixture()
  display <- summary(fit)$regression_displays$age
  display$info$term <- "haemoglobin"

  txt <- paste(
    capture.output(mfpi_print_summary_regression_block(display, digits = 3L)),
    collapse = "\n"
  )

  heading <- "Interaction: haemoglobin"
  expect_match(
    txt,
    paste0(heading, "\n", strrep("-", nchar(heading))),
    fixed = TRUE
  )
})


test_that("MFPI summary shift note is local to shifted interaction blocks", {
  fit <- make_mfpi_accessor_fixture()
  display_zero <- summary(fit)$regression_displays$age

  expect_false(mfpi_summary_display_has_nonzero_shift(display_zero))
  expect_length(
    capture.output(mfpi_print_summary_shift_note(display_zero)),
    0L
  )

  fit$shift["age"] <- 1
  display_shifted <- summary(fit)$regression_displays$age
  expect_true(mfpi_summary_display_has_nonzero_shift(display_shifted))
  expect_equal(display_shifted$info$shift, 1)
  expect_match(
    display_shifted$info$group_table$transformation[[1L]],
    "(age + 1)^-2",
    fixed = TRUE
  )

  txt <- paste(
    capture.output(mfpi_print_summary_regression_block(display_shifted, digits = 3)),
    collapse = "\n"
  )
  expect_match(
    txt,
    "The constant added inside the FP transformation is the shifting factor",
    fixed = TRUE
  )
  expect_match(
    txt,
    "applied to all observations before constructing the group-specific FP",
    fixed = TRUE
  )
  expect_match(
    txt,
    "structural zeros therefore remain untransformed",
    fixed = TRUE
  )

  group_pos <- regexpr("Group-specific FP terms:", txt, fixed = TRUE)[1L]
  shift_pos <- regexpr("Note: The constant added inside", txt, fixed = TRUE)[1L]
  group_var_pos <- regexpr("Group-variable coefficients:", txt, fixed = TRUE)[1L]
  expect_true(group_pos < shift_pos)
  expect_true(shift_pos < group_var_pos)
})


test_that("MFPI summary prints the shift explanation only once", {
  fit <- make_mfpi_accessor_fixture()

  display_zero <- summary(fit)$regression_displays$age

  fit$shift["age"] <- 1
  display_shifted_1 <- summary(fit)$regression_displays$age
  display_shifted_1$info$term <- "sz"
  display_shifted_1$info$group_table$transformation <- sub(
    "age + 1",
    "sz + 1",
    display_shifted_1$info$group_table$transformation,
    fixed = TRUE
  )

  display_shifted_2 <- display_shifted_1
  display_shifted_2$info$term <- "weight"
  display_shifted_2$info$shift <- 2
  display_shifted_2$info$group_table$transformation <- sub(
    "sz + 1",
    "weight + 2",
    display_shifted_2$info$group_table$transformation,
    fixed = TRUE
  )

  displays <- list(
    age_zero = display_zero,
    sz_shifted = display_shifted_1,
    weight_shifted = display_shifted_2
  )

  txt <- paste(
    capture.output(mfpi_print_summary_regression_displays(displays, digits = 3L)),
    collapse = "\n"
  )

  note_text <- "Note: The constant added inside the FP transformation is the shifting factor"
  note_hits <- gregexpr(note_text, txt, fixed = TRUE)[[1L]]
  expect_equal(sum(note_hits > 0L), 1L)

  first_shifted_pos <- regexpr("Interaction: sz\n---------------", txt)[1L]
  note_pos <- regexpr(note_text, txt, fixed = TRUE)[1L]
  second_shifted_pos <- regexpr("Interaction: weight", txt, fixed = TRUE)[1L]

  expect_true(first_shifted_pos < note_pos)
  expect_true(note_pos < second_shifted_pos)
  expect_match(txt, "(weight + 2)^-2", fixed = TRUE)
})


test_that("MFPI adjustment FP labels use compact transformation notation", {
  fit <- make_mfpi_accessor_fixture()

  # Mimic a selected FP2 adjustment term carried into an interaction model.
  # These are the source-column names produced by the MFP adjustment fit; the
  # MFPI display should reconstruct the selected FP basis rather than expose
  # the parse-like labels used by the generic mfp2 summary formatter.
  adj <- list(
    fp_terms = data.frame(
      selected = TRUE,
      df_final = 4,
      power1 = 3,
      power2 = 3,
      row.names = "age",
      check.names = FALSE
    ),
    x = structure(
      matrix(numeric(0), nrow = 0L, ncol = 2L),
      dimnames = list(NULL, c("age.1", "age.2"))
    ),
    coefficients = c(age.1 = -0.1, age.2 = 0.02),
    fp_powers = list(age = c(3, 3)),
    transformations = data.frame(
      shift = 0,
      scale = 1,
      center = TRUE,
      row.names = "age"
    ),
    term_to_columns = list(age = "age"),
    transformed_to_model_columns = c(
      age.1 = "age.1",
      age.2 = "age.2"
    ),
    transformed_column_to_source = c(
      age.1 = "age",
      age.2 = "age"
    ),
    transformed_column_component = c(
      age.1 = "fp_basis",
      age.2 = "fp_basis"
    ),
    transformed_column_zero_handled = c(
      age.1 = FALSE,
      age.2 = FALSE
    ),
    transformed_column_centered = c(
      age.1 = TRUE,
      age.2 = TRUE
    )
  )
  class(adj) <- "mfp2"
  fit$adjustment_model <- adj

  interaction_model <- list(
    transformed_to_model_columns = c(
      age.1 = "age.1",
      age.2 = "age.2"
    )
  )

  labels <- mfpi_adjustment_coefficient_labels(
    object = fit,
    interaction_model = interaction_model,
    coefficient_names = c("age.1", "age.2")
  )

  expect_identical(labels, c("age^3", "age^3 * log(age)"))
  expect_false(any(grepl("((age))", labels, fixed = TRUE)))

  # Genuine MFPI shifts remain visible, but with the same compact notation.
  fit$adjustment_model$transformations["age", "shift"] <- 2
  shifted <- mfpi_adjustment_coefficient_labels(
    object = fit,
    interaction_model = interaction_model,
    coefficient_names = c("age.1", "age.2")
  )
  expect_identical(
    shifted,
    c("(age + 2)^3", "(age + 2)^3 * log(age + 2)")
  )
})
