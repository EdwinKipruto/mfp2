
#' Print the Coefficient Block for a Multinomial `mfp2` Model
#'
#' Prints the fitted multinomial coefficient matrix carried by an `"mfp2"`
#' object. Rows correspond to non-reference outcomes and columns to design
#' terms. When the coefficient matrix is absent or empty, `(none)` is
#' printed so the print method never silently reports nothing.
#'
#' @param object An `"mfp2"` model object with a multinomial family.
#' @param digits Integer scalar controlling the number of digits used when
#'   formatting coefficient values. Default `3`.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of printing
#'   to the console.
#'
#' @keywords internal
#' @noRd
mfp2_print_multinomial_coefficients <- function(object, digits = 3L) {
  coefficients <- object$mfp2_coefficient_matrix
  if (!is.matrix(coefficients) || length(coefficients) == 0L) {
    cat("(none)\n")
    return(invisible(NULL))
  }
  base_object <- object
  class(base_object) <- setdiff(class(base_object), "mfp2")
  standard_errors <- tryCatch(
    summary(base_object)$standard.errors,
    error = function(e) matrix(
      NA_real_, nrow(coefficients), ncol(coefficients),
      dimnames = dimnames(coefficients)
    )
  )
  if (is.null(dim(standard_errors))) {
    standard_errors <- matrix(
      standard_errors,
      nrow = nrow(coefficients),
      dimnames = dimnames(coefficients)
    )
  }

  reference <- object$reference_class
  design_info <- mfp2_design_column_info(object)
  model_rows <- match(design_info$model_column, colnames(coefficients))
  if (anyNA(model_rows)) {
    stop("Final multinomial design metadata is not aligned with coefficients.",
         call. = FALSE)
  }
  show_center <- !is.null(object$centers)
  rows <- vector("list", nrow(coefficients))
  for (i in seq_len(nrow(coefficients))) {
    outcome <- rownames(coefficients)[[i]]
    intercept_position <- match("(Intercept)", colnames(coefficients))
    variable <- c(
      if (!is.na(intercept_position)) "(Intercept)",
      design_info$variable
    )
    basis <- c(if (!is.na(intercept_position)) "", design_info$basis)
    estimate <- c(
      if (!is.na(intercept_position)) coefficients[i, intercept_position],
      coefficients[i, model_rows]
    )
    se <- c(
      if (!is.na(intercept_position)) standard_errors[i, intercept_position],
      standard_errors[i, model_rows]
    )
    block <- data.frame(
      Logit = c(paste0(outcome, " vs ", reference),
                rep("", max(0L, length(variable) - 1L))),
      Variable = variable,
      Basis = basis,
      Estimate = unname(estimate),
      SE = unname(se),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (show_center) {
      block$Center <- c(
        if (!is.na(intercept_position)) NA_real_,
        design_info$center
      )
      block <- block[, c("Logit", "Variable", "Basis", "Center",
                         "Estimate", "SE")]
    }
    rows[[i]] <- block
  }
  display <- do.call(rbind, rows)
  if (show_center) {
    display$Center <- format_print_decimal(
      display$Center,
      digits,
      na_string = ""
    )
  }
  display <- format_model_print_table(display, digits)
  print.data.frame(
    display,
    row.names = FALSE,
    right = FALSE,
    quote = FALSE,
    na.print = ""
  )
  invisible(display)
}

#' Prepare Display Data for `print.mfp2()`
#'
#' Derives variable selection, functional-form labels, SAZ status, and display
#' helpers once so each output section uses the same interpretation of the fit.
#' The returned context is for internal printing only; it does not modify `x`.
#'
#' @param x A fitted `mfp2` object.
#' @param digits Number of decimal places for displayed values.
#' @return A named list of prepared display values and formatting functions.
#' @keywords internal
#' @noRd
mfp2_print_prepare <- function(x, digits) {
  output_width <- 78L

  make_boundary <- function(character = "-") {
    paste(rep(character, output_width), collapse = "")
  }

  # Full-width heading, used for top-level sections.
  print_section_heading <- function(title) {
    cat(make_boundary("-"), "\n", sep = "")
    cat(title, "\n")
    cat(make_boundary("-"), "\n\n", sep = "")
  }

  # Light heading, used for the Standard MFP / ACD / SAZ sub-tables nested
  # inside "Summary of Function Selection": an underline matching the title's
  # own width, so it reads as "one topic, several angles" rather than three
  # unrelated top-level sections.
  print_subsection_heading <- function(title) {
    cat(title, "\n", sep = "")
    cat(paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")
  }

  print_variable_list <- function(label, variables) {
    variable_text <- if (length(variables) > 0L) {
      paste(variables, collapse = ", ")
    } else {
      "none"
    }

    complete_text <- paste0(label, ": ", variable_text)

    wrapped_text <- strwrap(
      complete_text,
      width = output_width,
      exdent = nchar(label) + 2L
    )

    cat(paste(wrapped_text, collapse = "\n"), "\n", sep = "")
  }

  format_criterion <- function(value) {
    if (length(value) != 1L || is.na(value)) {
      return(NULL)
    }

    value <- as.character(value)
    criterion_key <- tolower(gsub("[^[:alnum:]]", "", value))

    switch(
      criterion_key,
      "pvalue" = "p-value",
      "aic" = "AIC",
      "bic" = "BIC",
      value
    )
  }

  # Build a single human-readable functional-form label from a variable's
  # selected powers together with its zero/catzero status. `acd_prefix`
  # should be FALSE whenever the calling table already has dedicated
  # power-on-x / power-on-A(x) columns (the ACD table), and TRUE only when
  # the table has no such columns but the row is nonetheless ACD-transformed
  # (an ACD variable that also went through the SAZ table).
  format_function_label <- function(powers, zero, catzero,
                                    searched_fp = FALSE,
                                    acd_prefix = FALSE,
                                    selected = TRUE,
                                    show_positive_domain = TRUE,
                                    binary_only_label = "binary indicator only") {
    # An unselected variable is always "out", regardless of any lingering
    # catzero/zero flag values -- those describe what was *requested*, not
    # what survived selection, and the two can disagree for eliminated
    # spike variables depending on how spike_decision was left set.
    if (!isTRUE(selected)) {
      return("out")
    }

    has_continuous <- length(powers) > 0L

    if (has_continuous) {
      base_label <- if (length(powers) == 1L && isTRUE(powers == 1) &&
                        !isTRUE(searched_fp)) {
        "linear"
      } else {
        sprintf("FP(%s)", paste(powers, collapse = ", "))
      }

      if (isTRUE(acd_prefix)) {
        base_label <- paste0("ACD ", base_label)
      }

      if (isTRUE(zero) && isTRUE(show_positive_domain)) {
        base_label <- paste0(base_label, " (x > 0)")
      }

      if (isTRUE(catzero)) {
        base_label <- paste0(base_label, " + binary")
      }

      return(base_label)
    }

    if (isTRUE(catzero)) {
      return(binary_only_label)
    }

    "out"
  }

  # Extract a variable's non-NA selected powers from a single fp_terms row,
  # given the names of the power1, power2, ... columns.
  extract_powers <- function(row, power_cols) {
    p <- as.numeric(row[power_cols])
    p[!is.na(p)]
  }

  # Sort a data.frame so that selected rows (selected_status = TRUE) come
  # first, preserving original relative order within each group.
  sort_selected_first <- function(df, selected_status) {
    df[order(!selected_status), , drop = FALSE]
  }


  # ---------------------------------------------------------------------------
  # Step 2: Prepare the detailed MFP table and selection indicators
  # ---------------------------------------------------------------------------

  # Work on a print-only copy. Do not modify `x$fp_terms`, because downstream
  # package code may depend on its original column names and numeric SAZ codes.
  fp_terms <- x$fp_terms

  variable_names <- rownames(fp_terms)

  if (is.null(variable_names) &&
      !is.null(x$fp_powers) &&
      length(x$fp_powers) == nrow(fp_terms)) {
    variable_names <- names(x$fp_powers)
  }

  if (is.null(variable_names) ||
      length(variable_names) != nrow(fp_terms)) {
    variable_names <- paste0("V", seq_len(nrow(fp_terms)))
  }

  rownames(fp_terms) <- variable_names

  if ("selected" %in% names(fp_terms)) {
    selected_status <- fp_terms[["selected"]]

    if (!is.logical(selected_status)) {
      selected_status <- tolower(as.character(selected_status)) %in%
        c("true", "t", "yes", "y", "1")
    }
  } else if ("df_final" %in% names(fp_terms)) {
    selected_status <- !is.na(fp_terms[["df_final"]]) &
      fp_terms[["df_final"]] > 0
  } else {
    selected_status <- rep(TRUE, nrow(fp_terms))
  }

  if (anyNA(selected_status)) {
    selected_fallback <- if ("df_final" %in% names(fp_terms)) {
      !is.na(fp_terms[["df_final"]]) &
        fp_terms[["df_final"]] > 0
    } else {
      rep(FALSE, nrow(fp_terms))
    }

    selected_status[is.na(selected_status)] <-
      selected_fallback[is.na(selected_status)]
  }

  acd_flag <- if ("acd" %in% names(fp_terms)) {
    as.logical(fp_terms[["acd"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }

  zero_flag <- if ("zero" %in% names(fp_terms)) {
    as.logical(fp_terms[["zero"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }

  catzero_flag <- if ("catzero" %in% names(fp_terms)) {
    as.logical(fp_terms[["catzero"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }

  spike_flag <- if ("spike" %in% names(fp_terms)) {
    sp <- fp_terms[["spike"]]
    if (!is.logical(sp)) {
      tolower(as.character(sp)) %in% c("true", "t", "yes", "y", "1")
    } else {
      sp
    }
  } else {
    rep(FALSE, nrow(fp_terms))
  }

  searched_fp_flag <- if ("searched_fp" %in% names(fp_terms)) {
    as.logical(fp_terms[["searched_fp"]])
  } else {
    rep(FALSE, nrow(fp_terms))
  }
  searched_fp_flag[is.na(searched_fp_flag)] <- FALSE

  decision_code <- if ("spike_dec" %in% names(fp_terms)) {
    suppressWarnings(as.integer(fp_terms[["spike_dec"]]))
  } else {
    rep(NA_integer_, nrow(fp_terms))
  }

  prop_zero_all <- if ("prop_zero" %in% names(fp_terms)) {
    suppressWarnings(as.numeric(fp_terms[["prop_zero"]]))
  } else {
    rep(NA_real_, nrow(fp_terms))
  }

  format_prop_zero <- function(value) {
    ifelse(
      is.na(value),
      ".",
      formatC(value, format = "f", digits = digits)
    )
  }

  power_cols <- grep("^power[0-9]+$", names(fp_terms), value = TRUE)
  powers_by_row <- lapply(seq_len(nrow(fp_terms)), function(i) {
    extract_powers(fp_terms[i, , drop = FALSE], power_cols)
  })

  df_initial_all <- if ("df_initial" %in% names(fp_terms)) fp_terms[["df_initial"]] else rep(NA, nrow(fp_terms))
  df_final_all   <- if ("df_final" %in% names(fp_terms)) fp_terms[["df_final"]] else rep(NA, nrow(fp_terms))
  df_change_all  <- sprintf("%s -> %s", df_initial_all, df_final_all)

  selected_variables <- variable_names[selected_status]
  excluded_variables <- variable_names[!selected_status]

  # SAZ decision text, one label per variable (only meaningful for
  # spike-eligible variables; used both in the SAZ table and in
  # Detailed Settings' saz_decision column).
  saz_decision_text <- rep("not SAZ", nrow(fp_terms))
  saz_decision_text[spike_flag & !selected_status] <- "not selected"
  saz_decision_text[spike_flag & selected_status] <- saz_decision_label(
    decision_code[spike_flag & selected_status],
    style = "print",
    unknown = "unknown"
  )
  # saz_decision_label(style = "print") returns "cont + binary" for the
  # combined decision; standardize to the fuller wording used elsewhere in
  # this print method.
  saz_decision_text[saz_decision_text == "cont + binary"] <- "continuous + binary"


  # ---------------------------------------------------------------------------
  # Step 3: Determine the model-selection criterion
  # ---------------------------------------------------------------------------

  criterion_label <- format_criterion(x$criterion_mfp)

  if (is.null(criterion_label)) {
    criterion_values <- character(0L)

    if ("select" %in% names(fp_terms)) {
      criterion_values <- c(criterion_values, as.character(fp_terms[["select"]]))
    }

    if ("alpha" %in% names(fp_terms)) {
      criterion_values <- c(criterion_values, as.character(fp_terms[["alpha"]]))
    }

    criterion_values <- unique(toupper(trimws(criterion_values)))
    criterion_values <- criterion_values[!is.na(criterion_values) & nzchar(criterion_values)]

    if ("AIC" %in% criterion_values) {
      criterion_label <- "AIC"
    } else if ("BIC" %in% criterion_values) {
      criterion_label <- "BIC"
    } else {
      criterion_label <- "p-value"
    }
  }

  is_pvalue_criterion <- identical(criterion_label, "p-value")



  in_saz <- spike_flag
  in_acd_only <- acd_flag & !spike_flag
  in_standard <- !acd_flag & !spike_flag

  list(
    output_width = output_width,
    make_boundary = make_boundary,
    print_section_heading = print_section_heading,
    print_subsection_heading = print_subsection_heading,
    print_variable_list = print_variable_list,
    format_function_label = format_function_label,
    sort_selected_first = sort_selected_first,
    fp_terms = fp_terms,
    variable_names = variable_names,
    selected_status = selected_status,
    acd_flag = acd_flag,
    zero_flag = zero_flag,
    catzero_flag = catzero_flag,
    spike_flag = spike_flag,
    searched_fp_flag = searched_fp_flag,
    prop_zero_all = prop_zero_all,
    format_prop_zero = format_prop_zero,
    powers_by_row = powers_by_row,
    df_change_all = df_change_all,
    excluded_variables = excluded_variables,
    saz_decision_text = saz_decision_text,
    criterion_label = criterion_label,
    is_pvalue_criterion = is_pvalue_criterion,
    in_saz = in_saz,
    in_acd_only = in_acd_only,
    in_standard = in_standard
  )
}

#' Format Covariate Preprocessing for Printing
#'
#' Displays integral shifts and scales without decimal places while retaining
#' significant fractional values. The stored transformation data is unchanged.
#'
#' @param transformations Data frame containing preprocessing settings.
#' @return A display copy with formatted `shift` and `scale` columns.
#' @keywords internal
#' @noRd
mfp2_format_preprocessing <- function(transformations) {
  display <- transformations
  format_value <- function(value) {
    if (is.na(value)) return("NA")
    if (!is.finite(value)) return(as.character(value))
    if (value == trunc(value)) {
      return(format(value, trim = TRUE, scientific = FALSE))
    }
    format(value, digits = 15L, trim = TRUE, scientific = FALSE)
  }
  for (column in intersect(c("shift", "scale"), names(display))) {
    if (is.numeric(display[[column]])) {
      display[[column]] <- vapply(
        display[[column]], format_value, character(1L), USE.NAMES = FALSE
      )
    }
  }
  display
}

#' Print Model Identity and Selected Functions
#'
#' Prints the call and model header, the selected and excluded variables,
#' preprocessing information, and the Standard MFP, ACD, and SAZ tables.
#'
#' @param x A fitted `mfp2` object.
#' @param digits Number of decimal places for displayed values.
#' @param ctx The display context returned by `mfp2_print_prepare()`.
#' @return Called for its printing side effects; its return value is unused.
#' @keywords internal
#' @noRd
mfp2_print_selection <- function(x, digits, ctx) {
  with(ctx, {
    # ---------------------------------------------------------------------------
    # Step 4: Print the main output banner
    # ---------------------------------------------------------------------------

    cat(make_boundary("="), "\n", sep = "")
    cat("MFP Model Fit\n")
    cat(make_boundary("="), "\n\n", sep = "")


    # ---------------------------------------------------------------------------
    # Step 5: Print the original mfp2 call and one-line meta block
    # ---------------------------------------------------------------------------
    #
    # A plain "Call:" label (no dashed rule) opens the output, followed by the
    # deparsed mfp2() call, a one-line model / criterion / convergence block,
    # and an observation/event line. This matches the opening produced by
    # print.summary.mfp2(); the two methods therefore present the call and
    # top-level metadata identically. The dashed-rule section headings resume
    # only from Selection Summary onward, where they bracket the more
    # substantial tabular sections.

    cat("Call:\n")
    if (!is.null(x$call_mfp)) {
      print(x$call_mfp)
    } else {
      cat("(original mfp2 call unavailable)\n")
    }
    cat("\n")

    # Model-specific metadata is extracted and rendered by the same helpers used
    # by print.summary.mfp2(), keeping direct and summary output synchronized.
    n_obs <- if (!is.null(x$nobs)) {
      x$nobs
    } else {
      tryCatch(NROW(x$y_original), error = function(e) NA_integer_)
    }
    n_events <- if (mfp2_family_uses_event_count(x$family_string) &&
                    !is.null(x$nevents)) x$nevents else NA_integer_
    model_metadata <- mfp2_summary_model_metadata(x)
    response_frequencies <- mfp2_summary_response_frequencies(x)
    mfp2_print_model_header(
      family_string = x$family_string,
      criterion = criterion_label,
      converged = x$convergence_mfp,
      n = n_obs,
      nevents = n_events,
      metadata = model_metadata,
      digits = digits,
      response_frequencies = response_frequencies
    )
    if (identical(x$family_string, "multinomial")) {
      cat(sprintf(
        "\nReference class: %s | Non-reference logits: %d\n",
        x$reference_class,
        x$n_logits
      ))
      cat("FP powers: common across logits\n")
    }

    cat("\n")


    # ---------------------------------------------------------------------------
    # Step 6: Print the model-selection summary
    # ---------------------------------------------------------------------------
    #
    # The Selection Summary section previously opened with "Converged:" and
    # "Criterion:" lines. Both fields are now shown in the top meta block above,
    # so they are omitted here to avoid duplication. The section retains its
    # dashed-rule heading and its three "Selected variables" lines.

    print_section_heading("Selection Summary")

    # A variable counts as "linear" if its sole continuous term is the fixed
    # linear row with power exactly 1 (and it is not ACD-transformed), or if it is
    # binary-only spike variable (no continuous component at all, just a 0/1
    # indicator). Every other selected variable is "nonlinear": any FP2
    # (including FP(1, 1), which contains a log(x) term), any FP1 with power
    # != 1, a searched FP1(power 1), and any ACD variable regardless of its
    # specific powers. Whether a zero/catzero indicator is additionally present
    # does not, by itself,
    # change this classification.
    is_plain_linear <- !acd_flag & !searched_fp_flag &
      vapply(powers_by_row, function(p) length(p) == 1L && isTRUE(p == 1), logical(1L))
    is_binary_only <- vapply(powers_by_row, length, integer(1L)) == 0L & catzero_flag

    linear_status <- selected_status & (is_plain_linear | is_binary_only)
    nonlinear_status <- selected_status & !linear_status

    print_variable_list(
      label = "Selected variables (linear)",
      variables = variable_names[linear_status]
    )

    print_variable_list(
      label = "Selected variables (nonlinear)",
      variables = variable_names[nonlinear_status]
    )

    print_variable_list(
      label = "Excluded variables",
      variables = excluded_variables
    )

    cat("\n")


    # ---------------------------------------------------------------------------
    # Step 7: Print covariate preprocessing information
    # ---------------------------------------------------------------------------

    print_section_heading("Covariate Preprocessing")

    if (!is.null(x$transformations)) {
      transformation_display <- mfp2_format_preprocessing(x$transformations)
      print.data.frame(transformation_display, right = FALSE)
    } else {
      cat("(preprocessing information unavailable)\n")
    }

    cat("\n")


    # ---------------------------------------------------------------------------
    # Step 8: Summary of Function Selection (Standard MFP / ACD / SAZ)
    # ---------------------------------------------------------------------------

    print_section_heading("Summary of Function Selection")

    function_label_plain <- mapply(
      FUN = format_function_label,
      powers = powers_by_row,
      zero = zero_flag,
      catzero = catzero_flag,
      searched_fp = searched_fp_flag,
      selected = selected_status,
      MoreArgs = list(acd_prefix = FALSE),
      SIMPLIFY = TRUE
    )

    # Partition: every variable belongs to exactly one of the three tables.
    # Spike-eligible variables go to SAZ regardless of ACD status (SAZ is the
    # more consequential decision for them); among the rest, ACD variables go
    # to the ACD table; everyone else goes to Standard MFP.
    # --- Standard MFP ----------------------------------------------------------

    if (any(in_standard)) {
      standard_table <- data.frame(
        Variable = variable_names[in_standard],
        Selected = ifelse(selected_status[in_standard], "yes", "no"),
        `df (init->final)` = df_change_all[in_standard],
        Function = unname(function_label_plain[in_standard]),
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      standard_table <- sort_selected_first(
        standard_table,
        selected_status[in_standard]
      )

      print_subsection_heading("Standard MFP")
      print.data.frame(standard_table, row.names = FALSE, right = FALSE)
      cat("\n")
    }

    # --- ACD (non-spike) --------------------------------------------------------

    if (any(in_acd_only)) {
      power1_all <- if ("power1" %in% names(fp_terms)) fp_terms[["power1"]] else rep(NA, nrow(fp_terms))
      power2_all <- if ("power2" %in% names(fp_terms)) fp_terms[["power2"]] else rep(NA, nrow(fp_terms))

      format_na_dot <- function(v) ifelse(is.na(v), ".", format(v, trim = TRUE))

      acd_table <- data.frame(
        Variable = variable_names[in_acd_only],
        Selected = ifelse(selected_status[in_acd_only], "yes", "no"),
        `df (init->final)` = df_change_all[in_acd_only],
        `Power on x` = format_na_dot(power1_all[in_acd_only]),
        `Power on A(x)` = format_na_dot(power2_all[in_acd_only]),
        Function = unname(function_label_plain[in_acd_only]),
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      acd_table <- sort_selected_first(acd_table, selected_status[in_acd_only])

      print_subsection_heading(
        "Approximate Cumulative Distribution (ACD) -- non-spike variables"
      )
      print.data.frame(acd_table, row.names = FALSE, right = FALSE)
      cat("\n")
    }

    # --- Spike-at-Zero (SAZ) -----------------------------------------------------

    if (any(in_saz)) {
      function_label_saz <- mapply(
        FUN = format_function_label,
        powers = powers_by_row[in_saz],
        zero = zero_flag[in_saz],
        catzero = catzero_flag[in_saz],
        searched_fp = searched_fp_flag[in_saz],
        acd_prefix = acd_flag[in_saz],
        selected = selected_status[in_saz],
        MoreArgs = list(
          show_positive_domain = FALSE,
          binary_only_label = "binary"
        ),
        SIMPLIFY = TRUE
      )

      # Build the SAZ table from columns that are always relevant. The ACD
      # indicator is intentionally omitted by default: an all-"no" column can
      # imply that ACD played a role even when no SAZ variable actually used it.
      saz_table <- data.frame(
        Variable = variable_names[in_saz],
        prop_zero = format_prop_zero(prop_zero_all[in_saz]),
        `df (init->final)` = df_change_all[in_saz],
        saz_decision = saz_decision_text[in_saz],
        Function = unname(function_label_saz),
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      # Add ACD only for an actual ACD/SAZ overlap. Restricting the check to
      # in_saz prevents ACD use elsewhere from adding an irrelevant SAZ column.
      if (any(acd_flag[in_saz])) {
        saz_table$ACD <- ifelse(acd_flag[in_saz], "yes", "no")
        saz_table <- saz_table[, c(
          "Variable", "prop_zero", "df (init->final)", "ACD",
          "saz_decision", "Function"
        ), drop = FALSE]
      }

      saz_table <- sort_selected_first(saz_table, selected_status[in_saz])

      print_subsection_heading(
        "Spike-at-Zero (SAZ; x is transformed only for x > 0, binary is I(x = 0))"
      )
      print.data.frame(saz_table, row.names = FALSE, right = FALSE)
      cat("\n")
    }


  })
}

#' Print Detailed Selection Settings
#'
#' Prints the variable settings table when requested. Its setting notes and
#' column explanations appear only when both settings and notes are enabled.
#'
#' @param detailed_settings Whether to print the settings table.
#' @param notes Whether to print explanatory notes.
#' @param digits Number of decimal places for displayed values.
#' @param ctx The display context returned by `mfp2_print_prepare()`.
#' @return Called for its printing side effects; its return value is unused.
#' @keywords internal
#' @noRd
mfp2_print_detailed_settings <- function(detailed_settings, notes, digits, ctx) {
  with(ctx, {
    # ---------------------------------------------------------------------------
    # Step 9: Detailed Settings (complete raw fp_terms table)
    # ---------------------------------------------------------------------------

    if (isTRUE(detailed_settings)) {

      print_section_heading("Detailed Settings")

      detailed_table <- data.frame(
        Selected = ifelse(selected_status, "yes", "no"),
        `df (init->final)` = df_change_all,
        row.names = variable_names,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      if (is_pvalue_criterion) {
        if ("select" %in% names(fp_terms)) detailed_table[["select"]] <- fp_terms[["select"]]
        if ("alpha" %in% names(fp_terms)) detailed_table[["alpha"]] <- fp_terms[["alpha"]]
      }
      detailed_table[["acd"]] <- ifelse(acd_flag, "yes", "no")
      detailed_table[["zero"]] <- ifelse(zero_flag, "yes", "no")
      detailed_table[["catzero_final"]] <- ifelse(catzero_flag, "yes", "no")
      detailed_table[["spike"]] <- ifelse(spike_flag, "yes", "no")
      detailed_table[["saz_decision"]] <- saz_decision_text
      if ("power1" %in% names(fp_terms)) detailed_table[["power1"]] <- fp_terms[["power1"]]
      if ("power2" %in% names(fp_terms)) detailed_table[["power2"]] <- fp_terms[["power2"]]

      detailed_table <- sort_selected_first(detailed_table, selected_status)
      detailed_table <- format_model_print_table(detailed_table, digits)

      print.data.frame(detailed_table, right = FALSE)
      cat("\n")

      # --- Notes and column definitions ---------------------------------------
      #
      # Everything here is auxiliary explanatory text, not part of the data
      # itself, so it is all controlled by the `notes` argument together.

      if (isTRUE(notes)) {

        note_lines <- character(0L)

        if (any(catzero_flag)) {
          note_lines <- c(
            note_lines,
            "catzero implies zero.",
            paste0(
              "catzero_final = yes means a binary indicator for the ",
              "x = 0 component is included in the final model."
            )
          )
        }

        if (is_pvalue_criterion && all(c("select", "alpha") %in% names(fp_terms))) {
          select_num <- suppressWarnings(as.numeric(fp_terms[["select"]]))
          alpha_num  <- suppressWarnings(as.numeric(fp_terms[["alpha"]]))

          mode_value <- function(v) {
            v <- v[!is.na(v)]
            if (length(v) == 0L) return(NA_real_)
            tab <- table(v)
            as.numeric(names(tab)[which.max(tab)])
          }

          # Forced inclusion and alpha = 1 are both deliberate variable-specific
          # modelling settings. Report them together so the note is concise while
          # still explaining that alpha = 1 permits the maximum function complexity.
          setting_note_parts <- character(0L)

          forced <- !is.na(select_num) & select_num == 1
          if (any(forced)) {
            forced_vars <- variable_names[forced]
            setting_note_parts <- c(
              setting_note_parts,
              sprintf(
                "%s %s forced into the model (select = 1)",
                paste(forced_vars, collapse = ", "),
                if (length(forced_vars) == 1L) "is" else "are"
              )
            )
          }

          alpha_one <- !is.na(alpha_num) & alpha_num == 1
          alpha_one_is_exception <- alpha_one & any(!is.na(alpha_num) & !alpha_one)
          if (any(alpha_one_is_exception)) {
            alpha_one_vars <- variable_names[alpha_one_is_exception]
            setting_note_parts <- c(
              setting_note_parts,
              sprintf(
                "%s use%s alpha = 1, allowing the maximum permitted function complexity",
                paste(alpha_one_vars, collapse = ", "),
                if (length(alpha_one_vars) == 1L) "s" else ""
              )
            )
          }

          if (length(setting_note_parts) > 0L) {
            note_lines <- c(
              note_lines,
              paste0(paste(setting_note_parts, collapse = "; "), ".")
            )
          }

          # Preserve reporting of other non-default selection settings. Forced
          # inclusion is handled above and is therefore excluded here.
          selectable <- !forced & !is.na(select_num)
          select_common <- mode_value(select_num[selectable])
          select_diff <- selectable & !is.na(select_common) &
            select_num != select_common
          if (any(select_diff)) {
            vars <- variable_names[select_diff]
            note_lines <- c(note_lines, sprintf(
              "%s use%s select = %s (all others use select = %s).",
              paste(vars, collapse = ", "),
              if (length(vars) == 1L) "s" else "",
              paste(unique(format(select_num[select_diff], trim = TRUE)), collapse = "/"),
              format(select_common, trim = TRUE)
            ))
          }

          # alpha = 1 has the informative wording above. Keep the established
          # generic note only for other variable-specific alpha values.
          alpha_common <- mode_value(alpha_num[!alpha_one])
          alpha_diff <- !alpha_one & !is.na(alpha_num) & !is.na(alpha_common) &
            alpha_num != alpha_common
          if (any(alpha_diff)) {
            vars <- variable_names[alpha_diff]
            note_lines <- c(note_lines, sprintf(
              "%s use%s alpha = %s.",
              paste(vars, collapse = ", "),
              if (length(vars) == 1L) "s" else "",
              paste(unique(format(alpha_num[alpha_diff], trim = TRUE)), collapse = "/")
            ))
          }
        }

        for (note in note_lines) {
          cat("Note: ", note, "\n", sep = "")
        }

        if (length(note_lines) > 0L) cat("\n")

        # --- Column definitions ---------------------------------------------
        # catzero_final's meaning is covered in the Notes above (together with
        # the catzero-implies-zero relationship), not repeated here.

        if (any(in_saz) && any(acd_flag[in_saz])) {
          cat(
            "Note: In the SAZ table, ACD = yes indicates that the variable also\n",
            "underwent an ACD transformation.\n\n",
            sep = ""
          )
        }

        if (any(zero_flag)) {
          cat(
            "zero:\n",
            "  yes means FP transformations are applied only to the positive ",
            "component\n  of the variable.\n\n",
            sep = ""
          )
        }

        if (any(spike_flag)) {
          cat(
            "spike:\n",
            "  yes means SAZ modelling was requested and remained eligible after ",
            "the\n  eligibility checks.\n\n",
            sep = ""
          )

          cat(
            "saz_decision:\n",
            "  Final SAZ status for each variable.\n\n",
            sep = ""
          )
        }

      }

    }


  })
}

#' Print Final-Model Coefficients
#'
#' Matches fitted coefficients to design-column labels and standard errors.
#' Handles multinomial logits, optional ordinal thresholds, and centering.
#'
#' @param x A fitted `mfp2` object.
#' @param digits Number of decimal places for displayed values.
#' @param intercepts Whether to display ordinal threshold intercepts.
#' @param ctx The display context returned by `mfp2_print_prepare()`.
#' @return Called for its printing side effects; its return value is unused.
#' @keywords internal
#' @noRd
mfp2_print_coefficients <- function(x, digits, intercepts, ctx) {
  with(ctx, {
    # ---------------------------------------------------------------------------
    # Step 10: Print final-model coefficients
    # ---------------------------------------------------------------------------

    print_section_heading("Final Model Coefficients")

    if (identical(x$family_string, "multinomial")) {
      mfp2_print_multinomial_coefficients(x, digits = digits)
    } else {
      coefficients <- stats::coef(x)

      if (length(coefficients) > 0L) {
        design_info <- mfp2_design_column_info(x)
        coefficient_names <- names(coefficients)

        if (is.null(coefficient_names)) {
          stop("Final model coefficients must be named.", call. = FALSE)
        }

        model_rows <- match(design_info$model_column, coefficient_names)
        if (anyNA(model_rows)) {
          stop(
            "Final design-column metadata is not aligned with model coefficients.",
            call. = FALSE
          )
        }

        # Standard errors come from the final fitted model covariance matrix and
        # are matched by exact coefficient names. Missing covariance entries remain
        # NA rather than being aligned by position or reconstructed from names.
        standard_errors <- stats::setNames(
          rep(NA_real_, length(coefficients)),
          coefficient_names
        )
        covariance <- if (mfp2_family_is_ordinal(x$family_string)) {
          mfp2_ordinal_vcov_slopes(x)
        } else {
          tryCatch(stats::vcov(x), error = function(e) NULL)
        }
        if (is.matrix(covariance) && !is.null(rownames(covariance)) &&
            !is.null(colnames(covariance))) {
          covariance_names <- intersect(
            coefficient_names,
            intersect(rownames(covariance), colnames(covariance))
          )
          if (length(covariance_names) > 0L) {
            covariance_rows <- match(covariance_names, rownames(covariance))
            covariance_cols <- match(covariance_names, colnames(covariance))
            variances <- covariance[cbind(covariance_rows, covariance_cols)]
            valid_variance <- !is.na(variances) & variances >= 0
            standard_errors[covariance_names[valid_variance]] <- sqrt(
              variances[valid_variance]
            )
          }
        }

        display <- data.frame(
          Variable = design_info$variable,
          Basis = design_info$basis,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        if (!is.null(x$centers)) {
          display[["Center"]] <- design_info$center
        }
        display[["Estimate"]] <- unname(coefficients[model_rows])
        display[["SE"]] <- unname(
          standard_errors[coefficient_names[model_rows]]
        )

        # Repeated rows for one variable represent additional basis columns. Blank
        # the repeated label while retaining the exact basis, center, and estimate.
        if (nrow(display) > 1L) {
          repeated_variable <- c(FALSE, display$Variable[-1L] == display$Variable[-nrow(display)])
          display$Variable[repeated_variable] <- ""
        }

        intercept_position <- match("(Intercept)", coefficient_names)
        if (!is.na(intercept_position)) {
          intercept <- data.frame(
            Variable = "(Intercept)",
            Basis = "",
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          if (!is.null(x$centers)) intercept[["Center"]] <- NA_real_
          intercept[["Estimate"]] <- unname(coefficients[[intercept_position]])
          intercept[["SE"]] <- unname(
            standard_errors[["(Intercept)"]]
          )
          display <- rbind(intercept, display)
        }

        # Ordinal models: prepend threshold intercepts when requested.
        # Mirrors rms::orm's print method: intercepts appear as the first rows of
        # the same coefficient table, above the predictor slopes.
        if (isTRUE(intercepts) && mfp2_family_is_ordinal(x$family_string)) {
          ord_int <- x$mfp2_ordinal_intercepts
          if (!is.null(ord_int) && length(ord_int) > 0L) {
            # SEs from vcov: rms::orm's vcov() defaults to intercepts = "mid"
            # (middle intercept only). Request all intercepts explicitly.
            # The mfp2 object IS the orm object (class c("mfp2", "orm")), so
            # vcov() dispatches directly on x, not on a $fit slot.
            covariance_full <- mfp2_ordinal_vcov_all(x)
            variances <- diag(
              covariance_full[names(ord_int), names(ord_int), drop = FALSE]
            )
            ord_se <- sqrt(variances)
            threshold_rows <- data.frame(
              Variable   = names(ord_int),
              Basis      = "",
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
            if (!is.null(x$centers)) threshold_rows[["Center"]] <- NA_real_
            threshold_rows[["Estimate"]]   <- unname(ord_int)
            threshold_rows[["SE"]] <- unname(ord_se)
            display <- rbind(threshold_rows, display)
          }
        }

        if (!is.null(x$centers)) {
          cat("Estimates are for: Basis - Center\n\n")
        }

        if (!is.null(x$centers)) {
          # Center is numeric while the table is assembled. Convert it to display
          # text only at the print boundary: a missing center is not applicable and
          # is left blank; actual centering constants use the common fixed-decimal
          # display convention.
          format_center <- function(value) {
            if (length(value) != 1L || is.na(value)) {
              return("")
            }
            format_print_decimal(value, digits)
          }

          display[["Center"]] <- vapply(
            display[["Center"]],
            FUN = format_center,
            FUN.VALUE = character(1L),
            USE.NAMES = FALSE
          )
        }

        display <- format_model_print_table(display, digits)

        print.data.frame(
          display,
          row.names = FALSE,
          right = FALSE,
          quote = FALSE,
          na.print = ""
        )
      } else {
        cat("(none)\n")
      }
    }

    cat("\n")


  })
}

#' Print Model-Fit Information
#'
#' Prints GEE working-correlation parameters when applicable, then renders
#' the family-specific fit statistics using the shared summary helper.
#'
#' @param x A fitted `mfp2` object.
#' @param digits Number of decimal places for displayed values.
#' @param notes Whether to print the model-df explanation.
#' @param ctx The display context returned by `mfp2_print_prepare()`.
#' @return Called for its printing side effects; its return value is unused.
#' @keywords internal
#' @noRd
mfp2_print_fit <- function(x, digits, notes, ctx) {
  with(ctx, {
    # ---------------------------------------------------------------------------
    # Step 11: GEE working-correlation parameters, then the Model Fit block
    # ---------------------------------------------------------------------------
    #
    # For marginal GEE models, surface the key GEE parameters (working
    # correlation and its estimate, robust SE type, scale, cluster structure)
    # using the same shared renderer as print.summary.mfp2().
    if (identical(x$family_string, "gee")) {
      mfp2_print_gee_parameters(
        mfp2_summary_gee_parameters(x),
        digits,
        print_section_heading
      )
    }

    # The Model Fit block (formerly "Model Deviances") is produced by
    # mfp2_format_model_fit_block(), the same helper called by
    # print.summary.mfp2(). This guarantees the two methods display the same
    # two family-specific fit statistics (GLM deviance or Cox -2 log L; df on
    # the MFP-adjusted scale) and the same optional explanatory note, so there is
    # exactly one place in the package where the block's format lives.
    mfp2_format_model_fit_block(
      values          = mfp2_summary_model_fit_values(x),
      digits          = digits,
      heading_printer = print_section_heading,
      notes           = notes
    )

    cat("\n")
    cat(make_boundary("="), "\n", sep = "")

  })
}

#' Print a Fitted MFP Model
#'
#' Displays the choices made by [mfp2()] and the resulting fitted model.
#' Use `print(fit)` for the full report, or set `detailed_settings = FALSE`
#' for a shorter view.
#'
#' @details
#' The report shows the original call and model details, selected and excluded
#' variables, preprocessing, selected functions, final coefficients with
#' standard errors, and model-fit statistics.
#' Preprocessing shifts and scales display without unnecessary trailing
#' zeros, independently of the `digits` print setting.
#'
#' Function choices are grouped by type: standard MFP, approximate cumulative
#' distribution (ACD), and spike-at-zero (SAZ). Selected variables appear
#' before excluded variables in each table. In the SAZ table, `prop_zero`
#' is the proportion of fitting observations with a structural zero;
#' `binary` denotes an indicator for `x = 0`. For SAZ variables, the
#' continuous component is transformed only when `x > 0`.
#'
#' The optional Detailed Settings table shows each variable's initial and
#' final degrees of freedom, selection settings, and transformation flags.
#' The coefficients table labels each fitted basis term and its standard
#' error. If the fit was centered, it also displays the centering constant
#' and explains how to read the estimate.
#'
#' Model-fit measures depend on the family. GLMs report deviances, Cox models
#' report minus twice the partial log likelihood, and other families display
#' their applicable statistics. Multinomial fits show one coefficient block
#' per non-reference outcome. Ordinal fits may include threshold intercepts.
#'
#' @param x A fitted `mfp2` object.
#' @param detailed_settings Show the Detailed Settings table. Default `TRUE`;
#'   set to `FALSE` for a shorter report.
#' @param notes Show explanatory notes beneath the settings and model-fit
#'   sections. Default `TRUE`.
#' @param intercepts Show ordinal threshold intercepts. By default (`NULL`),
#'   they are shown when there are fewer than 10 thresholds; use `TRUE` or
#'   `FALSE` to override this. Has no effect for other model families.
#' @param ... Additional print arguments. Use `digits` to choose the number
#'   of decimal places for displayed values (default 3); FP powers and
#'   degrees of freedom retain their own formatting.
#'
#' @return Invisibly returns `x`, allowing the fitted object to be reused.
#'
#' @examples
#' data("prostate")
#' fit <- mfp2(lpsa ~ fp(age) + svi, data = prostate,
#'             select = 1, alpha = 1, verbose = FALSE)
#' print(fit, detailed_settings = FALSE)
#'
#' @seealso [mfp2()], [summary.mfp2()]
#' @export
print.mfp2 <- function(x, detailed_settings = TRUE, notes = TRUE,
                       intercepts = NULL, ...) {
  dots <- list(...)
  digits <- dots$digits
  if (is.null(intercepts)) {
    n_int <- length(x$mfp2_ordinal_intercepts)
    intercepts <- mfp2_family_is_ordinal(x$family_string) && n_int > 0L && n_int < 10L
  } else {
    validate_logical_vector(intercepts, "intercepts", allowed_lengths = 1L)
  }
  validate_logical_vector(notes, "notes", allowed_lengths = 1L)
  if (is.null(digits)) digits <- 3L

  ctx <- mfp2_print_prepare(x, digits)
  mfp2_print_selection(x, digits, ctx)
  mfp2_print_detailed_settings(detailed_settings, notes, digits, ctx)
  mfp2_print_coefficients(x, digits, intercepts, ctx)
  mfp2_print_fit(x, digits, notes, ctx)
  invisible(x)
}
