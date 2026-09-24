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

  # Keep storage names stable while presenting reader-facing labels.
  display_names <- c(
    group_var = "Group",
    group_level = "Group level"
  )
  matched <- intersect(names(display_names), names(tab))
  names(tab)[match(matched, names(tab))] <- unname(display_names[matched])

  cat(sprintf("\n%s:\n", title))
  print(tab, row.names = FALSE)

  invisible(NULL)
}


#' Rename MFPI Metric-Table Columns for Public Display
#'
#' Renames only the metric columns that are meant for public display in the
#' printed MFPI selection tables. The stored metric-table names on the
#' fitted object remain unchanged so that programmatic access continues to
#' use stable identifiers.
#'
#' @param d Data frame of MFPI metrics, or `NULL`.
#'
#' @return A data frame with display-friendly column names, or `d`
#'   unchanged when the input is `NULL` or has no rows.
#'
#' @keywords internal
#' @noRd
rename_mfpi_metric_display_columns <- function(d) {
  if (is.null(d) || nrow(d) == 0L) return(d)
  if ("pvalue" %in% names(d)) {
    names(d)[names(d) == "pvalue"] <- "p-value"
  }
  d
}

#' Warn About Unused `...` Arguments to an MFPI Method
#'
#' Emits a single consistent warning when an MFPI method receives arguments
#' through `...` that it does not use. Centralising the warning here keeps
#' the messages identical across accessors, so users see one uniform
#' complaint regardless of which method they called.
#'
#' @param dots List of trailing arguments captured by the calling method.
#' @param method Character scalar naming the calling method, used to
#'   qualify the warning message.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
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
#' @param digits Integer or `NULL`. Optional number of decimal places to use for
#'   printed numeric output other than fractional-polynomial powers and degrees
#'   of freedom. Must be a single non-negative integer when supplied.
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

#' Round Numeric Columns for Console Display
#'
#' Rounds every numeric column of a data frame to the requested number of
#' digits for console display only. Must never be used upstream for stored
#' metrics or selection: rounding the stored numbers would change what
#' downstream comparisons and selection decisions see.
#'
#' @param d Data frame of values to display.
#' @param digits Integer scalar; number of significant digits.
#'
#' @return `d` with numeric columns rounded, or `d` unchanged when it is
#'   `NULL` or has no rows.
#'
#' @keywords internal
#' @noRd
round_print_numeric <- function(d, digits) {
  if (is.null(d) || nrow(d) == 0L) return(d)

  numeric_cols <- vapply(d, is.numeric, logical(1L))
  if (any(numeric_cols)) {
    d[numeric_cols] <- lapply(d[numeric_cols], round, digits = digits)
  }

  d
}

#' Format Numeric Values With a Fixed Number of Decimal Places
#'
#' Formats a numeric vector to a fixed number of decimal places for
#' console tables. Missing or non-finite values are replaced by `na_string`
#' so the printed column stays aligned.
#'
#' @param x Numeric vector to format.
#' @param digits Integer scalar; number of decimal places.
#' @param na_string Character scalar used for `NA` and non-finite values.
#'
#' @return Character vector of formatted values.
#'
#' @keywords internal
#' @noRd
format_print_decimal <- function(x, digits, na_string = "NA") {
  vapply(x, function(value) {
    if (length(value) != 1L || is.na(value)) return(na_string)
    formatC(value, format = "f", digits = digits)
  }, character(1L), USE.NAMES = FALSE)
}

#' Format P-Values Without Rounding a Positive Value to Zero
#'
#' Formats a numeric p-value vector to the requested decimal precision
#' without ever displaying a positive value as `0`. Values below the
#' display resolution are printed as `"<10^{-digits}"` so the printed
#' column always tells the truth about a positive test statistic.
#'
#' @param x Numeric vector of p values.
#' @param digits Integer scalar; number of decimal places.
#' @param na_string Character scalar used for `NA` values.
#'
#' @return Character vector of formatted p values.
#'
#' @keywords internal
#' @noRd
format_print_pvalue <- function(x, digits, na_string = "NA") {
  cutoff <- 10^(-digits)
  cutoff_label <- formatC(cutoff, format = "f", digits = digits)

  vapply(x, function(value) {
    if (length(value) != 1L || is.na(value)) return(na_string)
    if (is.finite(value) && value >= 0 && value < cutoff) {
      return(paste0("<", cutoff_label))
    }
    formatC(value, format = "f", digits = digits)
  }, character(1L), USE.NAMES = FALSE)
}

#' Format Model-Output Table Columns For Display
#'
#' Formats numeric columns of a model-output table with a fixed number of
#' decimals while preserving the natural notation for degrees of freedom
#' and fractional-polynomial power vectors. Character and integer columns
#' are passed through unchanged.
#'
#' @param d Data frame of model-output values, or `NULL`.
#' @param digits Integer scalar; number of decimal places.
#'
#' @return `d` with numeric columns formatted for printing.
#'
#' @keywords internal
#' @noRd
format_model_print_table <- function(d, digits) {
  if (is.null(d) || nrow(d) == 0L) return(d)

  numeric_columns <- names(d)[vapply(d, is.numeric, logical(1L))]
  for (column in numeric_columns) {
    key <- tolower(gsub("[^[:alnum:]]+", "_", column))
    is_df <- grepl("^df($|_)", key) || grepl("_df$", key)
    is_power <- grepl("power", key, fixed = TRUE)
    if (is_df || is_power) next

    is_pvalue <- key %in% c("p", "pvalue", "p_value", "p_raw", "p_adjusted") ||
      grepl("^pr_", key)
    d[[column]] <- if (is_pvalue) {
      format_print_pvalue(d[[column]], digits)
    } else {
      format_print_decimal(d[[column]], digits)
    }
  }

  d
}

#' Format Scalar Thresholds for Footnotes
#'
#' Formats a scalar threshold or comparison value printed in the
#' explanatory footnotes of MFPI output. Empty and `NA` inputs return
#' `"NA"` so the footnote remains legible.
#'
#' @param x Numeric scalar to format.
#' @param digits Integer scalar; number of decimal places.
#'
#' @return Character scalar with the formatted value.
#'
#' @keywords internal
#' @noRd
format_print_number <- function(x, digits) {
  if (length(x) == 0L || is.na(x)) return("NA")
  format_print_decimal(x[1L], digits)
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
    display = format_model_print_table(display, digits),
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
#' @param digits Integer scalar controlling displayed decimal places.
#' @param show_settings Logical scalar. Whether to print the interaction-selection
#'   criterion and threshold before the adjustment table. The top-level print
#'   method sets this to `FALSE` because it prints the same settings in Step 2.
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
#' @param digits Integer scalar. Number of decimal places used when printing
#'   numeric columns other than powers and degrees of freedom.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_candidates_step <- function(x, ruler, digits) {
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")

  settings <- mfpi_interaction_settings(
    criterion = if (!is.null(x$criterion)) x$criterion else "pvalue",
    p_interact = if (!is.null(x$p_interact)) x$p_interact else 0.05,
    min_improvement = if (!is.null(x$min_improvement)) x$min_improvement else 2,
    p_adjust_method = if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none",
    digits = digits
  )

  cat("\n")
  cat(sprintf("  Grouping variable     : %s\n", x$group_var))
  cat(sprintf("  Selection criterion   : %s\n", settings$criterion_label))
  cat(sprintf("  Selection rule        : %s\n", settings$selection))
  if (identical(settings$criterion, "pvalue")) {
    cat(sprintf("  P-value adjustment    : %s\n", settings$p_adjust_method))
  }
  cat("\n")

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
  all_m <- rename_mfpi_metric_display_columns(all_m)

  print(format_model_print_table(as.data.frame(all_m), digits), row.names = FALSE)

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
#' p-value-based selection, the function prints one p-value when no multiplicity
#' adjustment is used, and raw and adjusted p-values otherwise. For
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
#' @param digits Integer scalar. Number of decimal places used when printing
#'   numeric columns and the selection threshold in the heading.
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

  # The exact decision rule is printed with the interaction settings in Step 2,
  # so repeating it in this heading would add noise without information.
  heading <- sprintf("Step %d - Interaction Summary:", step_no)
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
      adjustment_method <- if (!is.null(x$p_adjust_method)) {
        tolower(as.character(x$p_adjust_method)[1L])
      } else {
        "none"
      }
      best_m <- format_mfpi_powers(best_m)
      best_m <- select_print_metric_columns(
        best_m,
        criterion = crit,
        include_adjusted = identical(crit, "pvalue") &&
          !identical(adjustment_method, "none")
      )
      if (identical(crit, "pvalue") &&
          identical(adjustment_method, "none")) {
        best_m[["p_adjusted"]] <- NULL
      }
      best_m <- rename_mfpi_metric_display_columns(best_m)
      best_m$selected <- "Yes"

      print(format_model_print_table(as.data.frame(best_m), digits), row.names = FALSE)
    } else {
      cat("  No interactions selected.\n")
    }

    return(invisible(NULL))
  }

  # Names of the continuous variables considered for interaction selection.
  cont_names <- names(vw)

  if (crit == "pvalue") {
    adjustment_method <- if (!is.null(x$p_adjust_method)) {
      tolower(as.character(x$p_adjust_method)[1L])
    } else {
      "none"
    }
    adjusted <- !identical(adjustment_method, "none")

    # Build one summary row per continuous variable for p-value based selection.
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]

      # If no valid winner was fitted for this variable, report the unavailable
      # metrics and an explicit non-selection decision.
      if (is.null(w$fit) || is.null(w$metric)) {
        if (adjusted) {
          data.frame(
            variable   = vn,
            type       = "---",
            p_raw      = NA_real_,
            p_adjusted = NA_real_,
            selected   = "No",
            check.names = FALSE
          )
        } else {
          data.frame(
            variable = vn,
            type = "---",
            `p-value` = NA_real_,
            selected = "No",
            check.names = FALSE
          )
        }
      } else {
        # Prefer adjusted p-values when present. If they are absent, use the raw
        # p-value as the displayed adjusted value.
        p_adj <- if ("p_adjusted" %in% names(w$metric)) {
          w$metric$p_adjusted[1L]
        } else {
          w$metric$pvalue[1L]
        }

        if (adjusted) {
          data.frame(
            variable   = vn,
            type       = format_type_label(w$type),
            p_raw      = w$metric$pvalue[1L],
            p_adjusted = p_adj,
            selected   = if (vn %in% selected_vars) "Yes" else "No",
            check.names = FALSE
          )
        } else {
          data.frame(
            variable = vn,
            type = format_type_label(w$type),
            `p-value` = w$metric$pvalue[1L],
            selected = if (vn %in% selected_vars) "Yes" else "No",
            check.names = FALSE
          )
        }
      }
    }))

    print(format_model_print_table(tab, digits), row.names = FALSE)
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

    print(format_model_print_table(tab, digits), row.names = FALSE)
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# print.mfpi()
# -----------------------------------------------------------------------------


#' Concise Model Description for MFPI Print Output
#'
#' Converts the normalised family identifier into the concise model
#' description shown immediately below the printed MFPI call. This mirrors
#' the mfp2 model label but is deliberately short so the MFPI print output
#' remains compact.
#'
#' @param x An MFPI object (or a compatible list with `family_string`).
#'
#' @return Character scalar with the model description.
#'
#' @keywords internal
#' @noRd
mfpi_print_model_label <- function(x) {
  family_string <- if (!is.null(x$family_string)) {
    as.character(x$family_string)[1L]
  } else {
    NA_character_
  }

  link <- NULL
  if (identical(family_string, "ordinal") && !is.null(x$ordinal_link)) {
    link <- as.character(x$ordinal_link)[1L]
  } else if (is.list(x$family) && !is.null(x$family$link)) {
    link <- as.character(x$family$link)[1L]
  }

  switch(
    family_string,
    cox = "Cox proportional hazards",
    survreg = "Parametric survival regression",
    multinomial = "Multinomial logistic regression",
    ordinal = if (is.null(link) || is.na(link) || !nzchar(link)) {
      "Ordinal regression"
    } else {
      switch(
        link,
        logistic = "Logistic ordinal regression (proportional odds)",
        probit = "Probit ordinal regression",
        loglog = "Log-log ordinal regression",
        cloglog = "Complementary log-log ordinal regression",
        cauchit = "Cauchit ordinal regression",
        "Ordinal regression"
      )
    },
    negbin = paste0(
      "Negative-binomial GLM",
      if (!is.null(link) && !is.na(link) && nzchar(link)) {
        paste0(" (", link, " link)")
      } else ""
    ),
    gaussian = "Gaussian GLM",
    binomial = "Binomial GLM",
    poisson = "Poisson GLM",
    Gamma = "Gamma GLM",
    inverse.gaussian = "Inverse-Gaussian GLM",
    if (!is.na(family_string) && nzchar(family_string)) family_string else "Unknown model"
  )
}


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
#' The header reports the original call and fitted model type. Step 1 displays
#' `df_initial` and `df_final`. For a grouped adjustment term,
#' `df_initial` counts all member design columns in the initial model.
#' Step 2 reports the grouping variable, interaction-selection criterion and
#' rule, and the p-value adjustment method when applicable. The candidate table
#' retains the stored metric names except that \code{pvalue} is printed as
#' \code{p-value}. The group-level power table uses the display labels
#' \code{Group} and \code{Group level}.
#'
#' The final interaction summary reports selection explicitly as \code{Yes} or
#' \code{No} without repeating the rule already shown in Step 2. With no
#' multiplicity adjustment it prints one \code{p-value} column; otherwise it
#' prints \code{p_raw} and \code{p_adjusted}.
#'
#' Numeric values use fixed decimal places for display only, except FP powers
#' and degrees of freedom, which retain their natural notation. The underlying
#' metrics stored in
#' the \code{"mfpi"} object are not modified, and model-selection decisions are
#' based on the unrounded values computed during fitting.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param digits Optional non-negative integer controlling the number of decimal
#'   places used when printing numeric output other than FP powers and degrees
#'   of freedom. If \code{NULL}, \code{x$digits} is used when available;
#'   otherwise 3 is used.
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
  header_text <- sprintf("MFPI  |  %s  |  %s", n_label, x$flex)

  header <- strrep("=", nchar(header_text, type = "width"))
  cat(header, "\n", header_text, "\n", header, "\n", sep = "")

  cat("\nCall:\n")
  if (!is.null(x$call)) {
    print(x$call)
  } else {
    cat("(original mfpi call unavailable)\n")
  }

  cat("\nModel: ", mfpi_print_model_label(x), "\n", sep = "")
  if (identical(x$family_string, "multinomial") &&
      !is.null(x$reference_class)) {
    cat("Reference outcome: ", x$reference_class, "\n", sep = "")
  }
  cat("\n")

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
  # the displayed p-value columns follow the adjustment method reported in
  # Step 2, so a separate adjustment-only section would be redundant.
  print_interaction_summary_step(x, ruler4, digits, step_no = 3L)

  # Follow the standard print-method convention: return the original object
  # invisibly so it can still be assigned or piped if needed.
  invisible(x)
}
# -----------------------------------------------------------------------------
# print.mfpi_prediction()
# -----------------------------------------------------------------------------

#' Format One MFPI Prediction Table for Console Display
#'
#' Formats a single MFPI prediction table for console display without
#' altering the stored prediction object. The public object retains its
#' original column names (`fit`, `se.fit`, and `x`) for backwards
#' compatibility; only the printed labels are made more descriptive.
#' For fitted-function output, `x` is displayed using the actual term name
#' (for example `age`).
#'
#' @param d Prediction data frame produced by `predict.mfpi()`, or `NULL`.
#' @param digits Integer scalar; number of decimal places.
#' @param term Optional character scalar naming the term, used to relabel
#'   the `x` column.
#'
#' @return A data frame with display-friendly column names, or `NULL` when
#'   `d` is `NULL`.
#'
#' @keywords internal
#' @noRd
format_mfpi_prediction_table <- function(d, digits, term = NULL) {
  if (is.null(d)) {
    return(NULL)
  }

  d <- as.data.frame(d)

  # A single mfpi_prediction is term-specific, so repeating the term in every
  # printed row adds noise. Keep it in the object but show it once in the header.
  if ("term" %in% names(d)) {
    d$term <- NULL
  }

  if (nrow(d) == 0L) {
    return(d)
  }

  display_names <- names(d)
  display_names[display_names == "fit"] <- "estimate"
  display_names[display_names == "se.fit"] <- "SE"
  if (!is.null(term) && length(term) == 1L && !is.na(term) && nzchar(term)) {
    display_names[display_names == "x"] <- term
  }
  names(d) <- display_names

  round_print_numeric(d, digits)
}


#' Print a Compact Ruler Around an MFPI Prediction Heading
#'
#' Prints a compact dynamic ruler around an MFPI prediction section
#' heading. The ruler width adapts to the header text so headings for
#' short and long term names line up consistently in the console.
#'
#' @param title Character scalar heading text.
#' @param term Character scalar term name shown in the heading.
#' @param detail Optional character scalar with extra qualifying text.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_mfpi_prediction_header <- function(title, term, detail = NULL) {
  pieces <- c(title, paste0("term: ", term))
  if (!is.null(detail) && length(detail) == 1L && !is.na(detail) && nzchar(detail)) {
    pieces <- c(pieces, detail)
  }

  header_text <- paste(pieces, collapse = "  |  ")
  header <- strrep("=", nchar(header_text, type = "width"))
  cat(header, "\n", header_text, "\n", header, "\n", sep = "")

  invisible(NULL)
}


#' Print an MFPI Prediction
#'
#' Prints the user-facing prediction results from an \code{"mfpi_prediction"}
#' object without exposing internal prediction metadata by default.
#'
#' For fitted-function predictions, the displayed curves are group-specific
#' partial fitted functions on the linear-predictor scale. Partial fitted
#' functions comprise the intercept (when present), group main effect, and
#' group-specific FP function estimated from the adjusted interaction model.
#'
#' The stored prediction object is not modified. In particular, the underlying
#' result tables retain the programmatic column names \code{fit},
#' \code{se.fit}, and \code{x}; these are printed as \code{estimate},
#' \code{SE}, and the actual term name (for example, \code{age}) only for
#' readability. Internal components such as \code{metadata} remain accessible
#' explicitly from the returned object.
#'
#' @param x An object of class \code{"mfpi_prediction"}, as returned by
#'   \code{predict.mfpi()}.
#' @param digits Optional non-negative integer controlling printed numeric
#'   precision. If \code{NULL}, 3 digits are used.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print mfpi_prediction
#' @export
print.mfpi_prediction <- function(x, digits = NULL, ...) {
  warn_unused_mfpi_dots(list(...), method = "print.mfpi_prediction")
  digits <- resolve_print_digits(x, digits)

  term <- if (!is.null(x$term) && length(x$term) > 0L) {
    paste(x$term, collapse = ", ")
  } else {
    "<unknown>"
  }
  type <- if (!is.null(x$type) && length(x$type) == 1L) x$type else "<unknown>"

  if (type %in% c("function", "difference", "both")) {
    print_mfpi_prediction_header(
      title = "MFPI prediction",
      term = term,
      detail = "scale: linear predictor"
    )

    if (type %in% c("function", "both")) {
      cat("\nPartial fitted functions:\n")
      tab <- format_mfpi_prediction_table(x$functions, digits, term = term)
      if (is.null(tab) || nrow(tab) == 0L) {
        cat("  No fitted-function results.\n")
      } else {
        print(tab, row.names = FALSE)
      }
    }

    if (type %in% c("difference", "both")) {
      cat("\nDifferences:\n")
      tab <- format_mfpi_prediction_table(x$differences, digits, term = term)
      if (is.null(tab) || nrow(tab) == 0L) {
        cat("  No fitted-function differences.\n")
      } else {
        print(tab, row.names = FALSE)
      }
    }

    if (type == "function") {
      cat(
        "\nNote: Partial fitted functions comprise the intercept (when present), ",
        "group main effect,\n",
        "and group-specific FP function estimated from the adjusted interaction ",
        "model.\n",
        sep = ""
      )
    } else if (type == "difference") {
      cat(
        "\nNote: Differences are comparison-group minus reference-group partial ",
        "fitted functions.\n",
        "Partial fitted functions comprise the intercept (when present), group ",
        "main effect, and\n",
        "group-specific FP function estimated from the adjusted interaction ",
        "model.\n",
        sep = ""
      )
    } else if (type == "both") {
      cat(
        "\nNote: Partial fitted functions comprise the intercept (when present), ",
        "group main effect,\n",
        "and group-specific FP function estimated from the adjusted interaction ",
        "model.\n",
        "Differences are comparison-group minus reference-group partial fitted ",
        "functions.\n",
        sep = ""
      )
    }

    return(invisible(x))
  }

  # Ordinary subject-level prediction. These predictions use the complete
  # term-specific interaction model, so the prediction type itself identifies
  # the relevant output scale (for example, link or response).
  print_mfpi_prediction_header(
    title = "MFPI prediction",
    term = term,
    detail = paste0("type: ", type)
  )
  cat("\n")

  tab <- format_mfpi_prediction_table(x$predictions, digits)
  if (is.null(tab) || nrow(tab) == 0L) {
    cat("  No prediction results.\n")
  } else {
    print(tab, row.names = FALSE)
  }

  invisible(x)
}


#' Print a List of MFPI Predictions
#'
#' Prints each term-specific \code{"mfpi_prediction"} in an
#' \code{"mfpi_prediction_list"} using the compact prediction print method.
#'
#' @param x An object of class \code{"mfpi_prediction_list"}.
#' @param digits Optional non-negative integer controlling printed numeric
#'   precision. If \code{NULL}, 3 digits are used.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print mfpi_prediction_list
#' @export
print.mfpi_prediction_list <- function(x, digits = NULL, ...) {
  warn_unused_mfpi_dots(list(...), method = "print.mfpi_prediction_list")
  digits <- resolve_print_digits(x, digits)

  if (length(x) == 0L) {
    cat("MFPI prediction list: no results.\n")
    return(invisible(x))
  }

  for (i in seq_along(x)) {
    if (i > 1L) cat("\n")
    print.mfpi_prediction(x[[i]], digits = digits)
  }

  invisible(x)
}
