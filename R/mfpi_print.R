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
#' Formats the stored interaction FP powers for one metric row. The main
#' interaction-candidate table intentionally shows only a compact summary,
#' because the detailed group-level mapping is printed separately.
#'
#' @param p_list Interaction FP power specification for one metric row.
#'
#' @return A character scalar suitable for console printing.
#'
#' @keywords internal
#' @noRd
format_interaction_power_label <- function(p_list) {
  # No interaction powers available.
  if (is.null(p_list) || length(p_list) == 0L) {
    return(".")
  }
  
  # A plain numeric vector represents one interaction function.
  if (!is.list(p_list)) {
    return("1 group-specific FP")
  }
  
  sprintf("%d group-specific FPs", length(p_list))
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
#' labels such as `"(1, 2)"`. Interaction powers are printed inline only when
#' the number of group-specific functions is small; otherwise a compact summary
#' such as `"9 group-specific functions"` is shown.
#'
#' @param d A data frame-like object containing MFPI model metrics. May include
#'   columns named `fp_powers_main` and/or `fp_powers_int`.
#' @param max_inline_int_powers Integer scalar. Maximum number of group-specific
#'   interaction power vectors to print inline in the main table.
#'
#' @return A base data frame with the same columns as `d`, except that
#'   `fp_powers_main` and `fp_powers_int`, when present, are converted to
#'   character vectors for printing.
#'
#' @keywords internal
#' @noRd
format_mfpi_powers <- function(d, max_inline_int_powers = 4L) {
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
      character(1L)
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
#' @param ruler Character scalar used as the horizontal separator printed below
#'   the title.
#' @param title Character scalar giving the table title.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
print_interaction_power_details_step <- function(metrics,
                                                 group_var = NULL,
                                                 ruler,
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
  
  cat(sprintf("\n%s:\n\n", title))
  
  print(tab, row.names = FALSE)
  
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

#' Print the adjustment-model step
#'
#' Prints the first step of an MFPI fit: the adjustment model selected by the
#' MFP procedure. The function displays selected adjustment variables and, when
#' applicable, reports adjustment variables that were dropped during MFP
#' selection.
#'
#' This helper is used internally by MFPI print methods. It writes directly to
#' the console using `cat()` and `print()` and returns `NULL` invisibly.
#'
#' @param x An MFPI fit object or print-ready MFPI object containing an
#'   `adjust_terms` component.
#' @param ruler Character scalar. Separator line printed below the section
#'   heading.
#' @param digits Integer scalar. Number of digits used when printing numeric
#'   columns.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_adjustment_step <- function(x, ruler, digits) {
  # Print the section heading and visual separator.
  cat("\nStep 1 - Adjustment Model (MFP):\n")
  cat(ruler, "\n")
  
  # Extract the adjustment-term table created by the initial MFP selection step.
  fp <- x$adjust_terms
  
  # If no adjustment model was fitted, report this explicitly and stop printing
  # this section.
  if (is.null(fp)) {
    cat("  No adjustment model fitted.\n")
    return(invisible(NULL))
  }
  
  # Split adjustment terms into variables selected by MFP and variables dropped
  # by MFP. NA values in `selected` are ignored in both subsets.
  selected_rows <- fp[!is.na(fp$selected) &  fp$selected, , drop = FALSE]
  dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, , drop = FALSE]
  
  # Print selected adjustment variables when available. Numeric columns are
  # rounded only for display; the stored object is not modified.
  if (nrow(selected_rows) > 0L) {
    print(round_print_numeric(selected_rows, digits))
  } else {
    cat("  No adjustment variables selected.\n")
  }
  
  # Report variables eliminated by MFP. Row names are used here because the MFP
  # adjustment table stores variable names as row names.
  if (nrow(dropped_rows) > 0L) {
    cat(sprintf(
      "\n  Dropped (%d variables eliminated by MFP): %s\n",
      nrow(dropped_rows),
      paste(rownames(dropped_rows), collapse = ", ")
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
  
  all_m <- format_mfpi_powers(all_m)
  all_m <- select_print_metric_columns(
    all_m,
    criterion = crit,
    include_adjusted = FALSE,
    include_interaction_powers = FALSE
  )
  
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
#' improvement score. Selected interactions are marked with `"*"`.
#'
#' This helper is used internally by MFPI print methods. It writes directly to
#' the console using `cat()` and `print()` and returns `NULL` invisibly.
#'
#' @param x An MFPI fit object or print-ready MFPI object containing
#'   `var_winners`, `best_model_metrics`, and optionally `criterion`,
#'   `p_interact`, and `min_improvement` components.
#' @param ruler Character scalar. Separator line printed below the section
#'   heading.
#' @param digits Integer scalar. Number of digits used when printing numeric
#'   columns.
#' @param step_no Integer scalar. Step number to display in the section heading.
#'   Defaults to `3L`.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_interaction_summary_step <- function(x, ruler, digits, step_no = 3L) {
  # Print the section heading and visual separator. The step number may differ
  # depending on whether a p-value adjustment section is printed separately.
  cat(sprintf("\nStep %d - Interaction Summary:\n", step_no))
  cat(ruler, "\n")
  
  # Resolve the criterion used for interaction selection. Fall back to p-value
  # based output for older or partial objects.
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  
  # Extract the per-variable winner list and the final selected interaction
  # metrics.
  vw <- x$var_winners
  best_m <- x$best_model_metrics
  
  # Determine which variables were selected into the final model. These are
  # later marked with "*" in the summary table.
  selected_vars <- if (!is.null(best_m) && nrow(best_m) > 0L) {
    best_m$variable
  } else {
    character(0L)
  }
  
  # Older or partial objects may not contain `var_winners`. In that case, fall
  # back to printing `best_model_metrics` directly when available.
  if (is.null(vw) || length(vw) == 0L) {
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      # Format FP power list-columns for display and keep only criterion-specific
      # print columns.
      best_m <- format_mfpi_powers(best_m)
      best_m <- select_print_metric_columns(best_m, criterion = crit)
      
      # Round numeric columns for display and print without row names.
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
      
      # If no valid winner was fitted for this variable, print an empty row.
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(
          variable   = vn,
          type       = "---",
          p_raw      = NA_real_,
          p_adjusted = NA_real_,
          selected   = "",
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
          selected   = if (vn %in% selected_vars) "*" else "",
          check.names = FALSE
        )
      }
    }))
    
    # Print the p-value summary table.
    print(round_print_numeric(tab, digits), row.names = FALSE)
    
    # Print the p-value adjustment method and explain the selection marker.
    method <- if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none"
    cat(sprintf("\n  p_adjust_method = %s", method))
    cat(sprintf(
      "\n  * = selected at p_interact = %s\n",
      format_print_number(x$p_interact, digits)
    ))
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
      
      # If no valid winner was fitted for this variable, print an empty row.
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(
          variable = vn,
          type     = "---",
          score    = NA_real_,
          selected = "",
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
          selected = if (vn %in% selected_vars) "*" else "",
          check.names = FALSE
        )
      }
    }))
    
    # Rename the generic score column to the criterion-specific display label.
    names(tab)[names(tab) == "score"] <- ic_label
    
    # Print the AIC/BIC summary table.
    print(round_print_numeric(tab, digits), row.names = FALSE)
    
    # Resolve and print the minimum improvement threshold used for selection.
    min_imp <- if (!is.null(x$min_improvement)) x$min_improvement else 2
    cat(sprintf(
      "\n  * = selected (%s > %s)\n",
      ic_label,
      format_print_number(min_imp, digits)
    ))
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
#' decision summary. The printed candidate table shows only the metrics relevant
#' to the selected criterion. Full diagnostics remain available in the returned
#' object.
#'
#' For \code{criterion = "pvalue"}, the final interaction summary includes
#' both raw and adjusted p-values together with the selection marker and
#' adjustment method.
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
#' @method print mfpi
#' @export
print.mfpi <- function(x, digits = NULL, ...) {
  # Prepare fixed-width separator lines used throughout the printed output.
  ruler  <- strrep("-", 85)
  header <- strrep("=", 85)
  
  # Resolve display metadata from the fitted object. Defaults are used for
  # backward compatibility with older or partially constructed objects.
  padj <- if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none"
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  
  # Resolve and validate the number of digits used for display rounding.
  digits <- resolve_print_digits(x, digits)
  
  # Print the main header with key model settings.
  cat(header, "\n")
  cat(sprintf(
    "MFPI  |  group: '%s'  |  n = %d  |  %s  |  criterion: %s",
    x$group_var,
    x$nobs,
    x$flex,
    crit
  ))
  
  # For p-value based selection, also show the p-value adjustment method in the
  # header because it affects the final interaction decisions.
  if (crit == "pvalue") {
    cat(sprintf("  |  p-adjust: %s", padj))
  }
  
  cat("\n")
  cat(header, "\n")
  
  # Use a shorter ruler inside individual sections.
  ruler2 <- strrep("-", 37)
  
  # Step 1: print adjustment-model MFP selection results.
  print_adjustment_step(x, ruler2, digits)
  
  # Step 2: print all interaction candidates and criterion-specific metrics.
  print_candidates_step(x, ruler2, digits)
  
  # Step 2b: print the detailed group-level mapping of interaction FP powers.
  # The main candidate table remains compact, while this section preserves the
  # interpretable mapping from group level to interaction FP powers.
  print_interaction_power_details_step(
    metrics   = x$all_model_metrics,
    group_var = x$group_var,
    ruler     = ruler2,
    title     = "  FP powers by group level"
  )
  
  # Step 3: print the final interaction summary. For p-value based selection,
  # the summary includes both raw and adjusted p-values, so a separate
  # adjustment-only section would be redundant.
  print_interaction_summary_step(x, ruler2, digits, step_no = 3L)
  
  # Follow the standard print-method convention: return the original object
  # invisibly so it can still be assigned or piped if needed.
  invisible(x)
}