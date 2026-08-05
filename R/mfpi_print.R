# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
#' Format one fractional polynomial power vector
#'
#' Converts a single fractional polynomial power vector into a compact character
#' label for printing. Missing, empty, or entirely `NA` power vectors are shown
#' as `"."`.
#'
#' @param p Numeric vector, list-like power specification, `NULL`, or `NA`.
#'
#' @return A character scalar such as `"(1, 2)"`, `"(0.5)"`, or `"."`.
#'
#' @keywords internal
#' @noRd
format_power_vector <- function(p) {
  # Treat absent, empty, or entirely missing power specifications as unavailable.
  if (is.null(p) || length(p) == 0L || all(is.na(p))) {
    return(".")
  }

  # If a one-element list is supplied, unwrap it. This supports future storage
  # such as list(age = c(1, 2)) while still printing only the power vector.
  if (is.list(p) && length(p) == 1L) {
    p <- p[[1L]]
  }

  # If the object is still list-like after unwrapping, collapse the unlisted
  # values. This is defensive and mainly protects older or irregular objects.
  if (is.list(p)) {
    p <- unlist(p, use.names = FALSE)
  }

  # Re-check after unlisting because the list may have contained only NULL/NA.
  if (length(p) == 0L || all(is.na(p))) {
    return(".")
  }

  paste0("(", paste(p, collapse = ", "), ")")
}


#' Format main-effect fractional polynomial powers
#'
#' Formats the stored main-effect FP powers for one metric row. The function
#' supports both the current storage style, where powers may be stored directly
#' as a numeric vector, and the proposed named-list style, for example
#' `list(age = c(1, 2))`.
#'
#' @param p Main-effect FP power specification for one metric row.
#'
#' @return A character scalar suitable for console printing.
#'
#' @keywords internal
#' @noRd
format_main_power_label <- function(p) {
  # For the main effect, there should normally be exactly one power vector.
  # Delegate the actual vector formatting to format_power_vector().
  format_power_vector(p)
}

#' Format interaction fractional polynomial powers compactly
#'
#' Formats the stored interaction FP powers for one metric row when a print
#' table explicitly requests that column. Small sets are displayed inline, while
#' larger sets use a compact count. The standard Step 2 candidate table omits
#' this column because the detailed group-level mapping follows immediately.
#'
#' @param p_list Interaction FP power specification for one metric row.
#' @param max_inline_int_powers Integer scalar. Maximum number of group-specific
#'   interaction power vectors to print inline.
#'
#' @return A character scalar suitable for console printing.
#'
#' @keywords internal
#' @noRd
format_interaction_power_label <- function(p_list,
                                           max_inline_int_powers = 4L) {
  # No interaction powers available.
  if (is.null(p_list) || length(p_list) == 0L) {
    return(".")
  }

  # A plain numeric vector represents one interaction function. Normalise it to
  # a list so singular/plural handling and inline formatting use one code path.
  if (!is.list(p_list)) {
    p_list <- list(p_list)
  }

  n_powers <- length(p_list)

  # Small sets are readable inline. Preserve complete group-level names when
  # available so the compact label remains interpretable.
  if (n_powers <= max_inline_int_powers) {
    power_labels <- vapply(p_list, format_power_vector, character(1L))
    group_labels <- names(p_list)

    if (!is.null(group_labels) &&
        length(group_labels) == n_powers &&
        all(!is.na(group_labels) & nzchar(group_labels))) {
      power_labels <- paste0(group_labels, ": ", power_labels)
    }

    return(paste(power_labels, collapse = "; "))
  }

  # Larger sets are summarised wherever this compact label is requested. The
  # complete mapping remains available through the dedicated detail table.
  sprintf(
    "%d group-specific FP%s",
    n_powers,
    if (n_powers == 1L) "" else "s"
  )
}

#' Format fractional polynomial power columns for printing
#'
#' Converts fractional polynomial power columns in an MFPI metrics data frame
#' from their stored numeric/list-column representation into compact character
#' labels suitable for console printing.
#'
#' The stored metric objects keep `fp_powers_main` and `fp_powers_int` as
#' list-columns so that the underlying powers remain machine-readable. This
#' helper is used only for display. Main-effect powers are printed as compact
#' labels such as `"(1, 2)"`. When interaction powers are requested by a
#' caller, small sets are printed inline and larger sets use a compact summary
#' such as `"9 group-specific FPs"`.
#'
#' @param d A data frame-like object containing MFPI model metrics. May include
#'   columns named `fp_powers_main` and/or `fp_powers_int`.
#' @param max_inline_int_powers Integer scalar. Maximum number of group-specific
#'   interaction power vectors to print inline when that column is retained.
#'
#' @return A base data frame with the same columns as `d`, except that
#'   `fp_powers_main` and `fp_powers_int`, when present, are converted to
#'   character vectors for printing.
#'
#' @keywords internal
#' @noRd
format_mfpi_powers <- function(d, max_inline_int_powers = 4L) {
  # Validate the inline threshold once before applying it row by row.
  if (!is.numeric(max_inline_int_powers) ||
      length(max_inline_int_powers) != 1L ||
      is.na(max_inline_int_powers) ||
      !is.finite(max_inline_int_powers) ||
      max_inline_int_powers < 0L ||
      max_inline_int_powers != floor(max_inline_int_powers)) {
    stop(
      "`max_inline_int_powers` must be a single non-negative integer.",
      call. = FALSE
    )
  }
  max_inline_int_powers <- as.integer(max_inline_int_powers)

  # Preserve NULL and empty inputs so callers can safely pass optional metric
  # tables without checking them first.
  if (is.null(d) || nrow(d) == 0L) {
    return(d)
  }

  # Drop tibble or other data-frame subclasses before printing. This keeps
  # console output based on standard data.frame behavior.
  d <- as.data.frame(d)

  # Format main-effect FP powers. Supports both old storage, e.g. c(1, 2), and
  # future named-list storage, e.g. list(age = c(1, 2)).
  if ("fp_powers_main" %in% names(d)) {
    d$fp_powers_main <- vapply(
      d$fp_powers_main,
      format_main_power_label,
      character(1L)
    )
  }

  # Format interaction FP powers. Large multi-level interactions are deliberately
  # summarised here; the detailed group-level table will show the full mapping.
  if ("fp_powers_int" %in% names(d)) {
    d$fp_powers_int <- vapply(
      d$fp_powers_int,
      format_interaction_power_label,
      character(1L),
      max_inline_int_powers = max_inline_int_powers
    )
  }

  d
}

#' Print Group-Level Interaction FP Powers
#'
#' Prints a supplementary table showing group-specific interaction-side
#' fractional polynomial powers for MFPI interaction models.
#'
#' This helper is used for both all-candidate output and retained-model output.
#' It is intentionally formatted as a detail table rather than as a separate
#' numbered model-building step. For example, candidate powers are printed below
#' the candidate-model table, and retained powers are printed before the
#' regression summaries in \code{print.summary.mfpi()}.
#'
#' The input metric table is expected to contain an \code{fp_powers_int}
#' list-column. Each row represents one interaction candidate or retained
#' interaction model. The helper expands that list-column into one display row
#' per group level.
#'
#' @param metrics Data frame of MFPI model metrics. Usually
#'   \code{x$all_model_metrics} for candidate output or
#'   \code{x$best_model_metrics} for retained-model output.
#' @param group_var Character scalar naming the grouping variable. Used only for
#'   display.
#' @param title Character scalar giving the table title.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
print_interaction_power_details_step <- function(metrics,
                                                 group_var = NULL,
                                                 title     = "  FP powers by group level") {
  tab <- mfpi_interaction_power_table(
    d = metrics,
    group_var = group_var
  )

  if (is.null(tab) || nrow(tab) == 0L) {
    return(invisible(NULL))
  }

  # Keep the display compact and predictable. The main candidate or retained
  # metric table already shows the main-effect FP powers, so this detail table
  # prints only the interaction-side powers by group level.
  cols <- c(
    "variable",
    "type",
    "group_var",
    "group_level",
    "fp_powers_int"
  )
  tab <- tab[, cols[cols %in% names(tab)], drop = FALSE]

  cat(sprintf("\n%s:\n", title))
  print(tab, row.names = FALSE)

  invisible(NULL)
}

# Warn consistently when an MFPI method receives unused arguments through `...`.
warn_unused_mfpi_dots <- function(dots, method) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("<unnamed>", length(dots))
  } else {
    missing_names <- is.na(dot_names) | !nzchar(dot_names)
    dot_names[missing_names] <- "<unnamed>"
  }

  warning(
    "Unused arguments in `",
    method,
    "(...)`: ",
    paste(dot_names, collapse = ", "),
    ".",
    call. = FALSE
  )

  invisible(NULL)
}

# Resolve the number of decimal places used by print.mfpi().
# User-supplied digits takes precedence over x$digits. If neither is available,
# use 3 for backward-compatible display.
#' Resolve the number of digits used for printing
#'
#' Determines the number of digits to use when printing an MFPI object or one
#' of its summary components. If `digits` is supplied, that value is validated
#' and returned. If `digits` is `NULL`, the function first tries to use
#' `x$digits`; if that is also unavailable, it falls back to `3`.
#'
#' This helper centralises validation of print precision so that all print and
#' summary methods use the same rules.
#'
#' @param x An object that may contain a `digits` element, typically an MFPI
#'   model or summary object.
#' @param digits Integer or `NULL`. Optional number of digits to use for
#'   printed numeric output. Must be a single non-negative integer when
#'   supplied.
#'
#' @return A single integer giving the number of digits to use for printing.
#'
#' @keywords internal
#' @noRd
resolve_print_digits <- function(x, digits = NULL) {
  # If the caller did not supply digits explicitly, try to recover the value
  # stored on the object. Fall back to 3 digits when neither source is available.
  if (is.null(digits)) {
    digits <- if (!is.null(x$digits)) x$digits else 3L
  }

  # Validate that digits is exactly one finite, non-negative integer-valued
  # numeric value. This rejects vectors, NA, Inf, negative values, and decimals.
  if (!is.numeric(digits) ||
      length(digits) != 1L ||
      is.na(digits) ||
      !is.finite(digits) ||
      digits < 0L ||
      digits != floor(digits)) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }

  # Return a clean integer value for downstream formatting helpers.
  as.integer(digits)
}

# Round numeric columns for console display only.
# This function must never be used upstream for stored metrics or selection.
round_print_numeric <- function(d, digits) {
  if (is.null(d) || nrow(d) == 0L) return(d)

  numeric_cols <- vapply(d, is.numeric, logical(1L))
  if (any(numeric_cols)) {
    d[numeric_cols] <- lapply(d[numeric_cols], round, digits = digits)
  }

  d
}

# Format scalar thresholds printed in explanatory footnotes.
format_print_number <- function(x, digits) {
  if (length(x) == 0L || is.na(x)) return("NA")
  format(round(x, digits = digits), trim = TRUE, scientific = FALSE)
}

#' Select metric columns for printing
#'
#' Selects the subset of MFPI metric columns that should be displayed in print
#' and summary output for a given model-selection criterion.
#'
#' The stored metric tables contain all available comparison quantities, such as
#' likelihood-ratio test statistics, AIC values, BIC values, and fractional
#' polynomial powers. This helper reduces those tables to the columns relevant
#' for the active selection criterion. For AIC and BIC output, it also creates
#' display aliases `dAIC` and `dBIC` from the stored improvement columns
#' `AIC_main_minus_int` and `BIC_main_minus_int` when needed.
#'
#' @param d A data frame-like object containing MFPI model metrics.
#' @param criterion Character scalar or `NULL`. Model-selection criterion used
#'   to choose the displayed metric columns. Supported values are `"pvalue"`,
#'   `"aic"`, and `"bic"`. If `NULL`, `"pvalue"` is used.
#' @param include_adjusted Logical scalar. If `TRUE` and `criterion` is
#'   `"pvalue"`, include the adjusted p-value column `p_adjusted` when present.
#' @param include_interaction_powers Logical scalar. If `TRUE`, include the
#'   `fp_powers_int` column when present. If `FALSE`, suppress it, which is
#'   useful when interaction powers are printed separately in a detailed
#'   group-level table.
#'
#' @return A data frame containing only the columns selected for printing.
#'   If `d` is `NULL` or has zero rows, `d` is returned unchanged.
#'
#' @keywords internal
#' @noRd
select_print_metric_columns <- function(d,
                                        criterion,
                                        include_adjusted = FALSE,
                                        include_interaction_powers = TRUE) {
  # Preserve NULL and empty inputs so callers can safely pass optional metric
  # tables without checking them first.
  if (is.null(d) || nrow(d) == 0L) {
    return(d)
  }

  # Use display labels for interaction type in printed output.
  if ("type" %in% names(d)) {
    d$type <- format_type_label(d$type)
  }

  # Default to p-value based output when no criterion is stored on the object.
  criterion <- if (!is.null(criterion)) criterion else "pvalue"

  # For AIC-based printing, expose the stored AIC improvement under the compact
  # display name dAIC. Positive values favour the interaction model because the
  # stored value is AIC_main - AIC_interaction.
  if (criterion == "aic" &&
      !("dAIC" %in% names(d)) &&
      "AIC_main_minus_int" %in% names(d)) {
    d$dAIC <- d$AIC_main_minus_int
  }

  # For BIC-based printing, expose the stored BIC improvement under the compact
  # display name dBIC. Positive values favour the interaction model because the
  # stored value is BIC_main - BIC_interaction.
  if (criterion == "bic" &&
      !("dBIC" %in% names(d)) &&
      "BIC_main_minus_int" %in% names(d)) {
    d$dBIC <- d$BIC_main_minus_int
  }

  # Select the criterion-specific statistic columns.
  criterion_cols <- switch(
    criterion,
    "pvalue" = c(
      "deviance_int",
      "deviance_diff",
      "df_int",
      "pvalue",
      if (isTRUE(include_adjusted)) "p_adjusted"
    ),
    "aic" = c("AIC_main", "AIC_interaction", "dAIC"),
    "bic" = c("BIC_main", "BIC_interaction", "dBIC"),
    character(0L)
  )

  # Always try to include identifying columns and main-effect FP powers.
  # Interaction FP powers can be suppressed when a separate group-level detail
  # table is printed.
  cols <- c(
    "variable",
    "type",
    criterion_cols,
    "fp_powers_main",
    if (isTRUE(include_interaction_powers)) "fp_powers_int"
  )
  cols <- unique(cols[cols %in% names(d)])

  d[, cols, drop = FALSE]
}

#' Build reader-facing MFPI interaction-selection settings
#'
#' Returns the criterion label and compact human-readable selection rule used by
#' both the standard print method and the verbose fitting output.
#'
#' @keywords internal
#' @noRd
mfpi_interaction_settings <- function(criterion = "pvalue",
                                      p_interact = 0.05,
                                      min_improvement = 2,
                                      p_adjust_method = "none",
                                      digits = 3L) {
  crit <- tolower(as.character(criterion)[1L])

  if (identical(crit, "pvalue")) {
    method <- as.character(p_adjust_method)[1L]
    if (is.na(method) || !nzchar(method)) method <- "none"
    no_adjustment <- identical(tolower(method), "none")
    selection <- sprintf(
      "%sp-value < %s",
      if (no_adjustment) "" else "adjusted ",
      format_print_number(p_interact, digits)
    )
    lines <- c(
      "criterion             : p-value",
      sprintf("interaction selection : %s", selection),
      sprintf("p-value adjustment    : %s", method)
    )
    return(list(
      criterion = crit,
      criterion_label = "p-value",
      selection = selection,
      p_adjust_method = method,
      lines = lines
    ))
  }

  if (crit %in% c("aic", "bic")) {
    label <- toupper(crit)
    selection <- sprintf(
      "%s reduction > %s",
      label,
      format_print_number(min_improvement, digits)
    )
    return(list(
      criterion = crit,
      criterion_label = label,
      selection = selection,
      p_adjust_method = NULL,
      lines = c(
        sprintf("criterion             : %s", label),
        sprintf("interaction selection : %s", selection)
      )
    ))
  }

  list(
    criterion = crit,
    criterion_label = as.character(criterion)[1L],
    selection = "",
    p_adjust_method = NULL,
    lines = sprintf("criterion             : %s", as.character(criterion)[1L])
  )
}

#' Prepare the selected adjustment-model table for printing
#'
#' Centralizes the selected-row filtering and criterion-specific column rules so
#' that `print.mfpi()` and verbose fitting output cannot drift apart.
#'
#' @keywords internal
#' @noRd
mfpi_prepare_adjustment_display <- function(x, digits) {
  settings <- mfpi_interaction_settings(
    criterion = if (!is.null(x$criterion)) x$criterion else "pvalue",
    p_interact = if (!is.null(x$p_interact)) x$p_interact else 0.05,
    min_improvement = if (!is.null(x$min_improvement)) x$min_improvement else 2,
    p_adjust_method = if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none",
    digits = digits
  )

  fp <- x$adjust_terms
  if (is.null(fp)) {
    return(list(
      settings = settings,
      display = NULL,
      selected_names = character(0L),
      model_fitted = FALSE
    ))
  }

  selected_rows <- fp[!is.na(fp$selected) & fp$selected, , drop = FALSE]
  selected_names <- rownames(selected_rows)

  if (nrow(selected_rows) == 0L) {
    selected_rows[["df_setting"]] <- NULL
    return(list(
      settings = settings,
      display = selected_rows,
      selected_names = selected_names,
      model_fitted = TRUE
    ))
  }

  display <- selected_rows
  display[["selected"]] <- NULL
  # `df_setting` is retained in the fitted object for internal diagnostics,
  # but the public table preserves the established df_initial/df_final display.
  display[["df_setting"]] <- NULL

  if (settings$criterion %in% c("aic", "bic")) {
    for (column in intersect(c("select", "alpha"), names(display))) {
      display[[column]] <- NULL
    }
  }

  logical_cols <- vapply(display, is.logical, logical(1L))
  if (any(logical_cols)) {
    display[logical_cols] <- lapply(
      display[logical_cols],
      function(col) ifelse(col, "yes", "no")
    )
  }

  if ("catzero" %in% names(display)) {
    names(display)[names(display) == "catzero"] <- "catzero_final"
  }

  spike_active <- if ("spike" %in% names(display)) {
    display[["spike"]] == "yes"
  } else {
    rep(FALSE, nrow(display))
  }

  if ("spike_dec" %in% names(display)) {
    dec_labels <- saz_decision_label(
      as.integer(display[["spike_dec"]]),
      style = "print",
      unknown = "unknown"
    )
    dec_labels[dec_labels == "cont + binary"] <- "continuous + binary"
    dec_labels[!spike_active] <- "not SAZ"
    display[["spike_dec"]] <- dec_labels
    names(display)[names(display) == "spike_dec"] <- "saz_decision"
  }

  if ("prop_zero" %in% names(display)) {
    prop_zero <- suppressWarnings(as.numeric(display[["prop_zero"]]))

    if (any(spike_active)) {
      display[["prop_zero"]] <- ifelse(
        spike_active & !is.na(prop_zero),
        formatC(prop_zero, format = "f", digits = digits),
        "."
      )
    } else {
      display[["prop_zero"]] <- NULL
    }
  }

  flag_cols <- intersect(
    c("acd", "zero", "catzero_final", "spike", "saz_decision"),
    names(display)
  )
  for (fc in flag_cols) {
    col_vals <- display[[fc]]
    if (all(col_vals %in% c("no", FALSE, NA, ".", "not SAZ"))) {
      display[[fc]] <- NULL
    }
  }

  list(
    settings = settings,
    display = round_print_numeric(display, digits),
    selected_names = selected_names,
    model_fitted = TRUE
  )
}

#' Print the adjustment-model step
#'
#' Prints only the adjustment variables retained by MFP, using the same
#' criterion-specific settings and columns as verbose fitting.
#'
#' @param x An MFPI fit object or print-ready object containing `adjust_terms`.
#' @param ruler Character scalar printed below the section heading.
#' @param digits Integer scalar controlling displayed numeric precision.
#' @param show_settings Logical scalar. Whether to print the interaction-selection
#'   criterion and threshold before the adjustment table. The top-level print
#'   method sets this to `FALSE` because it prints the same settings globally.
#'
#' @return Invisibly `NULL`.
#' @keywords internal
#' @noRd
print_adjustment_step <- function(x, ruler, digits, show_settings = TRUE) {
  cat("\nStep 1 - Selected Adjustment Model (MFP):\n")
  cat(ruler, "\n")

  info <- mfpi_prepare_adjustment_display(x, digits)

  # State the interaction-selection criterion alongside the adjustment-model
  # results. Keep this delegated to `mfpi_interaction_settings()` (via
  # `mfpi_prepare_adjustment_display()`) so p-value adjustment methods and
  # AIC/BIC thresholds are rendered identically in every output path.
  if (isTRUE(show_settings)) {
    cat("\n", paste0(" ", info$settings$lines, "\n"), "\n", sep = "")
  }

  if (!isTRUE(info$model_fitted)) {
    cat("  No adjustment model fitted.\n")
    return(invisible(NULL))
  }

  if (!is.null(info$display) && nrow(info$display) > 0L) {
    print(info$display)
  } else {
    cat("  No adjustment variables selected.\n")
  }

  if (length(info$selected_names) > 0L) {
    cat(sprintf(
      "\n  Selected adjustment variables (%d): %s\n",
      length(info$selected_names),
      paste(info$selected_names, collapse = ", ")
    ))
  }

  invisible(NULL)
}

#' Print all interaction candidates
#'
#' Prints the second step of an MFPI fit: the complete set of interaction
#' candidates evaluated after the adjustment model has been fitted.
#'
#' The function retrieves `all_model_metrics` from the MFPI object, formats
#' fractional polynomial power columns for display, selects the metric columns
#' relevant to the active selection criterion, rounds numeric columns for
#' printing, and writes the resulting table to the console.
#'
#' This helper is used internally by MFPI print methods. It writes directly to
#' the console using `cat()` and `print()` and returns `NULL` invisibly.
#'
#' @param x An MFPI fit object or print-ready MFPI object containing an
#'   `all_model_metrics` component and, optionally, a `criterion` component.
#' @param ruler Character scalar. Separator line printed below the section
#'   heading.
#' @param digits Integer scalar. Number of digits used when printing numeric
#'   columns.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_candidates_step <- function(x, ruler, digits) {
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")

  all_m <- x$all_model_metrics
  if (is.null(all_m) || nrow(all_m) == 0L) {
    cat("  No interaction candidates computed.\n")
    return(invisible(NULL))
  }

  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"

  # Remove interaction powers before display formatting. Their complete
  # group-level mapping is printed in the dedicated table immediately below.
  all_m <- select_print_metric_columns(
    all_m,
    criterion = crit,
    include_adjusted = FALSE,
    include_interaction_powers = FALSE
  )
  all_m <- format_mfpi_powers(all_m)

  print(round_print_numeric(as.data.frame(all_m), digits), row.names = FALSE)

  invisible(NULL)
}

#' Print the interaction-summary step
#'
#' Prints the interaction-selection summary for an MFPI fit. The summary shows,
#' for each continuous variable, the best interaction candidate found during the
#' MFPI interaction search and whether that candidate was selected for the final
#' model.
#'
#' The displayed columns depend on the model-selection criterion. For
#' p-value-based selection, the function prints raw and adjusted p-values. For
#' AIC- or BIC-based selection, it prints the corresponding information-criterion
#' improvement score. The final column reports the decision explicitly as
#' `"Yes"` or `"No"`. The exact criterion-specific selection rule is included in
#' the section heading, so no explanatory marker note is required after the
#' table.
#'
#' This helper is used internally by MFPI print methods. It writes directly to
#' the console using `cat()` and `print()` and returns `NULL` invisibly.
#'
#' @param x An MFPI fit object or print-ready MFPI object containing
#'   `var_winners`, `best_model_metrics`, and optionally `criterion`,
#'   `p_adjust_method`, `p_interact`, and `min_improvement` components.
#' @param ruler Character scalar. Minimum separator line printed below the
#'   section heading. The line is extended when necessary to match the heading.
#' @param digits Integer scalar. Number of digits used when printing numeric
#'   columns and the selection threshold in the heading.
#' @param step_no Integer scalar. Step number to display in the section heading.
#'   Defaults to `3L`.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_interaction_summary_step <- function(x, ruler, digits, step_no = 3L) {
  # Resolve the criterion used for interaction selection. Fall back to p-value
  # based output for older or partial objects.
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"

  # Put the exact decision rule in the heading. P-value selection uses the
  # generic displayed field `p_adjusted`, regardless of the adjustment method;
  # AIC/BIC selection uses the stored minimum-improvement threshold.
  if (crit == "pvalue") {
    threshold <- format_print_number(x$p_interact, digits)
    rule <- sprintf("selected when p_adjusted < %s", threshold)
  } else {
    ic_label <- if (crit == "aic") "dAIC" else "dBIC"
    min_imp <- if (!is.null(x$min_improvement)) x$min_improvement else 2
    rule <- sprintf(
      "selected when %s > %s",
      ic_label,
      format_print_number(min_imp, digits)
    )
  }

  heading <- sprintf(
    "Step %d - Interaction Summary (%s):",
    step_no,
    rule
  )
  heading_ruler <- strrep("-", max(nchar(heading), nchar(ruler)))

  cat("\n", heading, "\n", sep = "")
  cat(heading_ruler, "\n")

  # Extract the per-variable winner list and the final selected interaction
  # metrics.
  vw <- x$var_winners
  best_m <- x$best_model_metrics

  # Determine which variables were selected into the final model.
  selected_vars <- if (!is.null(best_m) && nrow(best_m) > 0L) {
    best_m$variable
  } else {
    character(0L)
  }

  # Older or partial objects may not contain `var_winners`. In that case, fall
  # back to printing `best_model_metrics` directly when available. Every row in
  # that retained-only table is selected by definition.
  if (is.null(vw) || length(vw) == 0L) {
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      best_m <- format_mfpi_powers(best_m)
      best_m <- select_print_metric_columns(best_m, criterion = crit)
      best_m$selected <- "Yes"

      print(round_print_numeric(as.data.frame(best_m), digits), row.names = FALSE)
    } else {
      cat("  No interactions selected.\n")
    }

    return(invisible(NULL))
  }

  # Names of the continuous variables considered for interaction selection.
  cont_names <- names(vw)

  if (crit == "pvalue") {
    # Build one summary row per continuous variable for p-value based selection.
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]

      # If no valid winner was fitted for this variable, report the unavailable
      # metrics and an explicit non-selection decision.
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(
          variable   = vn,
          type       = "---",
          p_raw      = NA_real_,
          p_adjusted = NA_real_,
          selected   = "No",
          check.names = FALSE
        )
      } else {
        # Prefer adjusted p-values when present. If they are absent, use the raw
        # p-value as the displayed adjusted value.
        p_adj <- if ("p_adjusted" %in% names(w$metric)) {
          w$metric$p_adjusted[1L]
        } else {
          w$metric$pvalue[1L]
        }

        data.frame(
          variable   = vn,
          type       = format_type_label(w$type),
          p_raw      = w$metric$pvalue[1L],
          p_adjusted = p_adj,
          selected   = if (vn %in% selected_vars) "Yes" else "No",
          check.names = FALSE
        )
      }
    }))

    print(round_print_numeric(tab, digits), row.names = FALSE)
  } else {
    # Resolve the displayed information-criterion column. The display aliases
    # are dAIC and dBIC, but stored metric rows may contain only the original
    # improvement columns.
    ic_label <- if (crit == "aic") "dAIC" else "dBIC"
    stored_ic_col <- if (crit == "aic") {
      "AIC_main_minus_int"
    } else {
      "BIC_main_minus_int"
    }

    # Build one summary row per continuous variable for AIC/BIC selection.
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]

      # If no valid winner was fitted for this variable, report the unavailable
      # score and an explicit non-selection decision.
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(
          variable = vn,
          type     = "---",
          score    = NA_real_,
          selected = "No",
          check.names = FALSE
        )
      } else {
        # Prefer the display alias if already present. Otherwise use the stored
        # improvement column.
        val <- if (ic_label %in% names(w$metric)) {
          w$metric[[ic_label]][1L]
        } else if (stored_ic_col %in% names(w$metric)) {
          w$metric[[stored_ic_col]][1L]
        } else {
          NA_real_
        }

        data.frame(
          variable = vn,
          type     = format_type_label(w$type),
          score    = val,
          selected = if (vn %in% selected_vars) "Yes" else "No",
          check.names = FALSE
        )
      }
    }))

    # Rename the generic score column to the criterion-specific display label.
    names(tab)[names(tab) == "score"] <- ic_label

    print(round_print_numeric(tab, digits), row.names = FALSE)
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# print.mfpi()
# -----------------------------------------------------------------------------


#' Print an \code{"mfpi"} Object
#'
#' Displays the adjustment model, criterion-specific interaction model metrics,
#' detailed group-level interaction FP powers, and the final interaction
#' decision summary. The printed candidate table shows the metrics relevant
#' to the selected criterion and the main-effect FP powers. Interaction FP
#' powers are omitted from that table because the complete group-level mapping
#' is printed immediately below it. Full diagnostics remain available in the
#' returned object.
#'
#' Step 1 reports the interaction-selection rule in reader-facing language and
#' displays `df_initial` and `df_final`. For a grouped adjustment term,
#' `df_initial` counts all member design columns in the initial model.
#' P-value fits print the p-value cutoff and the stored adjustment method,
#' including \code{none}, while retaining the variable-specific \code{select} and
#' \code{alpha} columns. AIC/BIC fits print the required criterion reduction
#' and omit only those two inapplicable columns.
#'
#' The final interaction summary reports selection explicitly as \code{Yes} or
#' \code{No}. Its heading contains the exact criterion-specific rule:
#' \code{p_adjusted < p_interact} for \code{criterion = "pvalue"}, or the
#' minimum \code{dAIC}/\code{dBIC} improvement for information-criterion
#' selection.
#'
#' Numeric values are rounded for display only. The underlying metrics stored in
#' the \code{"mfpi"} object are not modified, and model-selection decisions are
#' based on the unrounded values computed during fitting.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param digits Optional non-negative integer controlling the number of digits
#'   used when printing numeric output. If \code{NULL}, \code{x$digits} is used
#'   when available; otherwise 3 is used.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @seealso [mfpi()], [summary.mfpi()]
#'
#' @method print mfpi
#' @export
print.mfpi <- function(x, digits = NULL, ...) {
  warn_unused_mfpi_dots(list(...), method = "print.mfpi")

  # Resolve display metadata from the fitted object. Defaults are used for
  # backward compatibility with older or partially constructed objects.
  crit <- if (!is.null(x$criterion)) tolower(x$criterion) else "pvalue"
  criterion_label <- if (crit %in% c("aic", "bic")) {
    toupper(crit)
  } else if (crit == "pvalue") {
    "p-value"
  } else {
    crit
  }

  # Resolve and validate the number of digits used for display rounding.
  digits <- resolve_print_digits(x, digits)

  # Build the sample-size label. For survival models, include
  # the number of events in brackets so the user sees both n and events.
  n_label <- if (!is.null(x$nevents)) {
    sprintf("n = %d (events = %d)", x$nobs, x$nevents)
  } else {
    sprintf("n = %d", x$nobs)
  }

  # Build the complete header before printing it so the surrounding ruler has
  # exactly the same display width as the text. This avoids an oversized fixed
  # ruler and remains aligned when model metadata changes.
  header_text <- sprintf(
    "MFPI  |  group: '%s'  |  %s  |  %s  |  criterion: %s",
    x$group_var,
    n_label,
    x$flex,
    criterion_label
  )

  header <- strrep("=", nchar(header_text, type = "width"))
  cat(header, "\n", header_text, "\n", header, "\n", sep = "")

  # Print the global analysis settings (criterion, selection rule, adjustment
  # method) between the header and the numbered steps so the reader sees the
  # ground rules before any results appear.
  settings <- mfpi_interaction_settings(
    criterion       = crit,
    p_interact      = if (!is.null(x$p_interact)) x$p_interact else 0.05,
    min_improvement = if (!is.null(x$min_improvement)) x$min_improvement else 2,
    p_adjust_method = if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none",
    digits          = digits
  )
  cat("\n", paste0(" ", settings$lines, "\n"), "\n", sep = "")

  # Use a shorter ruler inside individual sections.
  ruler2 <- strrep("-", 41)
  ruler3 <- strrep("-", 36)
  ruler4 <- strrep("-", 29)

  # Step 1: print adjustment-model MFP selection results.
  print_adjustment_step(x, ruler2, digits, show_settings = FALSE)

  # Step 2: print all interaction candidates and criterion-specific metrics.
  print_candidates_step(x, ruler3, digits)

  # Step 2b: print the detailed group-level mapping of interaction FP powers.
  # These powers are deliberately omitted from the candidate table above to
  # avoid duplicating the same information in adjacent sections.
  print_interaction_power_details_step(
    metrics   = x$all_model_metrics,
    group_var = x$group_var,
    title     = "  FP powers by group level"
  )

  # Step 3: print the final interaction summary. For p-value based selection,
  # the summary includes both raw and adjusted p-values, so a separate
  # adjustment-only section would be redundant.
  print_interaction_summary_step(x, ruler4, digits, step_no = 3L)

  # Follow the standard print-method convention: return the original object
  # invisibly so it can still be assigned or piped if needed.
  invisible(x)
}