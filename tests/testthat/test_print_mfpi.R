# Focused tests split from the former test_mfp2.R monolith.


# Migrated coverage from the former test_mfp2.R

test_that("small interaction-power sets print inline and large sets summarize", {
  small <- format_mfpi_powers(
    make_mfpi_print_metrics(2L),
    max_inline_int_powers = 2L
  )
  expect_identical(small$fp_powers_int, "A: (1); B: (2)")

  large <- format_mfpi_powers(
    make_mfpi_print_metrics(3L),
    max_inline_int_powers = 2L
  )
  expect_identical(large$fp_powers_int, "3 group-specific FPs")

  singleton <- format_mfpi_powers(
    make_mfpi_print_metrics(1L),
    max_inline_int_powers = 0L
  )
  expect_identical(singleton$fp_powers_int, "1 group-specific FP")
})


test_that("inline interaction-power threshold is validated", {
  metrics <- make_mfpi_print_metrics(2L)

  expect_error(
    format_mfpi_powers(metrics, max_inline_int_powers = -1L),
    "single non-negative integer",
    fixed = TRUE
  )
  expect_error(
    format_mfpi_powers(metrics, max_inline_int_powers = 1.5),
    "single non-negative integer",
    fixed = TRUE
  )
})


test_that("p-value Step 1 prints the cutoff, adjustment method, select, and alpha", {
  x <- make_mfpi_adjustment_object("pvalue")
  x$p_adjust_method <- "none"
  x$p_interact <- 0.075
  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")
  info <- mfpi_prepare_adjustment_display(x, digits = 3L)
  expect_true(all(c("df_initial", "df_final") %in% names(info$display)))
  expect_false("df_setting" %in% names(info$display))
  expect_true(all(c("select", "alpha") %in% names(info$display)))
  expect_false(grepl("df_setting", printed, fixed = TRUE))
  expect_match(printed, "df_initial", fixed = TRUE)
  expect_match(printed, "df_final", fixed = TRUE)

  expect_match(printed, " criterion             : p-value", fixed = TRUE)
  expect_match(
    printed,
    " interaction selection : p-value < 0.075",
    fixed = TRUE
  )
  expect_match(
    printed,
    " p-value adjustment    : none",
    fixed = TRUE
  )
  expect_match(printed, "select", fixed = TRUE)
  expect_match(printed, "alpha", fixed = TRUE)
})


test_that("p-value Step 1 prints the active adjustment method dynamically", {
  x <- make_mfpi_adjustment_object("pvalue")
  x$p_adjust_method <- "hochberg"
  x$p_interact <- 0.01
  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(
    printed,
    " interaction selection : adjusted p-value < 0.01",
    fixed = TRUE
  )
  expect_match(
    printed,
    " p-value adjustment    : hochberg",
    fixed = TRUE
  )
})


test_that("AIC and BIC Step 1 use dynamic human-readable thresholds", {
  for (criterion in c("aic", "bic")) {
    x <- make_mfpi_adjustment_object(criterion)
    x$min_improvement <- 3.5
    output <- capture.output(
      print_adjustment_step(x, ruler = "-----", digits = 3L)
    )
    printed <- paste(output, collapse = "\n")
    info <- mfpi_prepare_adjustment_display(x, digits = 3L)
    criterion_label <- toupper(criterion)

    expect_match(
      printed,
      sprintf(" criterion             : %s", criterion_label),
      fixed = TRUE
    )
    expect_match(
      printed,
      sprintf(" interaction selection : %s reduction > 3.5", criterion_label),
      fixed = TRUE
    )
    expect_false(grepl("min_improvement", printed, fixed = TRUE))
    expect_false(grepl("df_setting", printed, fixed = TRUE))
    expect_match(printed, "df_initial", fixed = TRUE)
    expect_match(printed, "df_final", fixed = TRUE)
    expect_match(printed, "power1", fixed = TRUE)
    expect_match(printed, "power2", fixed = TRUE)
    expect_false("select" %in% names(info$display))
    expect_false("alpha" %in% names(info$display))
    expect_match(printed, "Selected adjustment variables (2): hx, age", fixed = TRUE)
  }
})


test_that("MFPI Step 1 prints grouped initial and final df only", {
  x <- make_minimal_mfpi_print_object()
  x$adjust_terms <- data.frame(
    df_setting = 1,
    df_initial = 2,
    select = 0.05,
    alpha = 0.05,
    selected = TRUE,
    df_final = 2,
    power1 = 1,
    power2 = NA_real_,
    row.names = "stage",
    check.names = FALSE
  )

  info <- mfpi_prepare_adjustment_display(x, digits = 3L)
  expect_false("df_setting" %in% names(info$display))
  expect_equal(info$display[["df_initial"]], 2)
  expect_equal(info$display[["df_final"]], 2)

  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")
  expect_false(grepl("df_setting", printed, fixed = TRUE))
  expect_match(printed, "df_initial", fixed = TRUE)
  expect_match(printed, "df_final", fixed = TRUE)
})


test_that("candidate output omits interaction powers shown in the detail table", {
  x <- make_minimal_mfpi_print_object()
  x$all_model_metrics <- make_mfpi_print_metrics(2L)

  output <- capture.output(
    print_candidates_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "fp_powers_main", fixed = TRUE)
  expect_false(grepl("fp_powers_int", printed, fixed = TRUE))
  expect_false(grepl("A: (1); B: (2)", printed, fixed = TRUE))
})


test_that("interaction-power detail output has no underline", {
  output <- capture.output(
    print_interaction_power_details_step(
      metrics = make_mfpi_print_metrics(2L),
      group_var = "group",
      title = "  FP powers by group level"
    )
  )

  expect_match(paste(output, collapse = "\n"), "FP powers by group level:", fixed = TRUE)
  expect_false(any(grepl("^-+$", trimws(output))))
})


test_that("main header ruler matches the displayed header width", {
  x <- make_minimal_mfpi_print_object()
  x$nevents <- 7L
  output <- capture.output(print(x))

  expect_identical(output[[1L]], output[[3L]])
  expect_identical(
    nchar(output[[1L]], type = "width"),
    nchar(output[[2L]], type = "width")
  )
})


test_that("criterion labels are shown in settings, not the main header", {
  for (criterion in c("aic", "bic", "pvalue")) {
    x <- make_minimal_mfpi_print_object()
    x$criterion <- criterion
    output <- capture.output(print(x))
    expected <- if (criterion == "pvalue") "p-value" else toupper(criterion)

    # The main header identifies the fit only; criterion details are printed
    # once in the settings block below the ruler.
    expect_false(grepl("criterion:", output[[2L]], fixed = TRUE))
    expect_match(
      paste(output, collapse = "\n"),
      sprintf("criterion             : %s", expected),
      fixed = TRUE
    )

    # The criterion should not be duplicated elsewhere in the printed output.
    criterion_lines <- grep(
      "^\\s*criterion\\s*:",
      output,
      value = TRUE
    )
    expect_length(criterion_lines, 1L)
    expect_false(grepl("p-adjust:", output[[2L]], fixed = TRUE))
  }
})


test_that("MFPI methods warn consistently about unused dots", {
  x <- make_minimal_mfpi_print_object()

  expect_warning(
    capture.output(print(x, unused_argument = TRUE)),
    "Unused arguments in `print.mfpi(...)`: unused_argument.",
    fixed = TRUE
  )
  expect_warning(
    warn_unused_mfpi_dots(list(TRUE), method = "summary.mfpi"),
    "Unused arguments in `summary.mfpi(...)`: <unnamed>.",
    fixed = TRUE
  )
})


test_that("Step 3 reports Yes and No and embeds the p-value rule in its heading", {
  x <- make_mfpi_interaction_summary_object("pvalue")
  output <- capture.output(
    print_interaction_summary_step(
      x,
      ruler = "-----",
      digits = 3L,
      step_no = 3L
    )
  )
  printed <- paste(output, collapse = "\n")

  expect_match(
    printed,
    "Step 3 - Interaction Summary (selected when p_adjusted < 0.05):",
    fixed = TRUE
  )
  expect_match(printed, "selected", fixed = TRUE)
  expect_match(printed, "Yes", fixed = TRUE)
  expect_match(printed, "No", fixed = TRUE)
  expect_false(grepl("* = selected", printed, fixed = TRUE))
  expect_false(grepl("p_adjust_method =", printed, fixed = TRUE))
})


test_that("Step 3 embeds the AIC and BIC rules and uses Yes and No", {
  for (criterion in c("aic", "bic")) {
    x <- make_mfpi_interaction_summary_object(criterion)
    output <- capture.output(
      print_interaction_summary_step(
        x,
        ruler = "-----",
        digits = 3L,
        step_no = 3L
      )
    )
    printed <- paste(output, collapse = "\n")
    label <- if (criterion == "aic") "dAIC" else "dBIC"

    expect_match(
      printed,
      sprintf(
        "Step 3 - Interaction Summary (selected when %s > 2):",
        label
      ),
      fixed = TRUE
    )
    expect_match(printed, "Yes", fixed = TRUE)
    expect_match(printed, "No", fixed = TRUE)
    expect_false(grepl("* = selected", printed, fixed = TRUE))
  }
})


test_that("shared adjustment display contains selected variables only", {
  x <- make_mfpi_adjustment_object("pvalue")
  dropped <- data.frame(
    df_setting = 4,
    df_initial = 4,
    select = 0.05,
    alpha = 0.05,
    df_final = 0,
    power1 = NA_real_,
    power2 = NA_real_,
    selected = FALSE,
    row.names = "hg",
    check.names = FALSE
  )
  x$adjust_terms <- rbind(x$adjust_terms, dropped)

  info <- mfpi_prepare_adjustment_display(x, digits = 3L)

  expect_identical(rownames(info$display), c("hx", "age"))
  expect_identical(info$selected_names, c("hx", "age"))
  expect_false("selected" %in% names(info$display))
  expect_false("df_setting" %in% names(info$display))
  expect_false("hg" %in% rownames(info$display))
})


test_that("MFPI Step 1 uses mfp2 SAZ display names and labels", {
  x <- make_minimal_mfpi_print_object()
  x$p_adjust_method <- "none"
  x$adjust_terms <- data.frame(
    df_setting = c(4, 4, 4),
    df_initial = c(4, 4, 4),
    select = c(0.05, 1.00, 1.00),
    alpha = c(0.05, 0.05, 0.05),
    acd = c(FALSE, FALSE, TRUE),
    zero = c(FALSE, TRUE, FALSE),
    catzero = c(FALSE, TRUE, FALSE),
    spike = c(FALSE, TRUE, FALSE),
    prop_zero = c(NA_real_, 0.25, NA_real_),
    spike_dec = c(2L, 1L, 2L),
    selected = c(TRUE, TRUE, TRUE),
    df_final = c(2, 2, 1),
    power1 = c(0, 1, 1),
    power2 = c(NA_real_, NA_real_, NA_real_),
    row.names = c("cavol", "pgg45", "age"),
    check.names = FALSE
  )

  info <- mfpi_prepare_adjustment_display(x, digits = 3L)

  expect_true("catzero_final" %in% names(info$display))
  expect_true("prop_zero" %in% names(info$display))
  expect_true("saz_decision" %in% names(info$display))
  expect_false("catzero" %in% names(info$display))
  expect_false("spike_dec" %in% names(info$display))
  expect_identical(info$display$prop_zero, c(".", "0.250", "."))
  expect_identical(
    info$display$saz_decision,
    c("not SAZ", "continuous + binary", "not SAZ")
  )

  output <- capture.output(
    print_adjustment_step(x, ruler = "-----", digits = 3L)
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "catzero_final", fixed = TRUE)
  expect_match(printed, "prop_zero", fixed = TRUE)
  expect_match(printed, "0.250", fixed = TRUE)
  expect_match(printed, "saz_decision", fixed = TRUE)
  expect_match(printed, "not SAZ", fixed = TRUE)
  expect_match(printed, "continuous + binary", fixed = TRUE)
  expect_false(grepl("spike_dec", printed, fixed = TRUE))
})


test_that("verbose Step 3 reports final p-value decisions without stars", {
  winners <- list(
    age = list(
      fit = list(ok = TRUE),
      metric = data.frame(pvalue = 0.0273),
      type = "fp2",
      score = 0.0273
    ),
    wt = list(
      fit = list(ok = TRUE),
      metric = data.frame(pvalue = 0.663),
      type = "fp1",
      score = 0.663
    )
  )

  output <- capture.output(
    print_interaction_step3_summary(
      var_winners = winners,
      cont_vars = c("age", "wt"),
      mode = "pvalue",
      p_interact = 0.05,
      p_adjust_method = "none",
      digits = 4L
    )
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "p_raw", fixed = TRUE)
  expect_match(printed, "p_adjusted", fixed = TRUE)
  expect_match(printed, "selected", fixed = TRUE)
  expect_match(printed, "Yes", fixed = TRUE)
  expect_match(printed, "No", fixed = TRUE)
  expect_false(grepl("* = selected", printed, fixed = TRUE))
  expect_false(grepl("p_interact", printed, fixed = TRUE))
})


test_that("verbose Step 3 reports final AIC decisions without stars", {
  winners <- list(
    age = list(
      fit = list(ok = TRUE),
      metric = data.frame(AIC_main_minus_int = 3.681),
      type = "fp2",
      score = 3.681
    ),
    wt = list(
      fit = list(ok = TRUE),
      metric = data.frame(AIC_main_minus_int = -1.824),
      type = "fp1",
      score = -1.824
    )
  )

  output <- capture.output(
    print_interaction_step3_summary(
      var_winners = winners,
      cont_vars = c("age", "wt"),
      mode = "ic",
      ic_col = "AIC_main_minus_int",
      ic_label = "dAIC",
      min_improvement = 2,
      digits = 3L
    )
  )
  printed <- paste(output, collapse = "\n")

  expect_match(printed, "dAIC", fixed = TRUE)
  expect_match(printed, "selected", fixed = TRUE)
  expect_match(printed, "Yes", fixed = TRUE)
  expect_match(printed, "No", fixed = TRUE)
  expect_false(grepl("* = selected", printed, fixed = TRUE))
})


test_that("verbose candidate evaluation contains no provisional decision text", {
  body_text <- paste(
    deparse(body(evaluate_interaction_for_variable)),
    collapse = "\n"
  )

  expect_match(body_text, "p-value =", fixed = TRUE)
  expect_false(grepl("provisional", body_text, fixed = TRUE))
  expect_false(grepl("No significant interaction retained", body_text, fixed = TRUE))
  expect_false(grepl("Not selected", body_text, fixed = TRUE))
  expect_false(grepl("Selected  (", body_text, fixed = TRUE))
})
