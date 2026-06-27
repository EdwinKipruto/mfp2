# S3 print method for mfpi objects
#
# print.mfpi() gives a structured output:
#   Step 1 - Adjustment model (fp_terms for selected variables)
#   Step 2 - All interaction models with their criterion-specific metrics
#   Step 3 - P-value adjustment when criterion = "pvalue"
#   Step 4 - Interaction summary when criterion = "pvalue"
# For AIC/BIC, Step 3 is the interaction summary.
#
# Numeric rounding in this file is display-only. Stored model metrics remain
# unrounded so that selection decisions are not affected by print precision.

# -----------------------------------------------------------------------------
# Helpers ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Convert FP power columns to compact character labels before printing.
format_mfpi_powers <- function(d) {
  if (is.null(d) || nrow(d) == 0L) return(d)
  
  if ("fp_powers_main" %in% names(d)) {
    d$fp_powers_main <- vapply(d$fp_powers_main, function(p) {
      if (is.null(p) || length(p) == 0L || all(is.na(p))) "."
      else paste0("(", paste(p, collapse = ", "), ")")
    }, character(1L))
  }
  
  if ("fp_powers_int" %in% names(d)) {
    d$fp_powers_int <- vapply(d$fp_powers_int, function(p_list) {
      if (is.null(p_list) || length(p_list) == 0L) return(".")
      if (!is.list(p_list)) {
        return(paste0("(", paste(p_list, collapse = ", "), ")"))
      }
      paste(vapply(p_list, function(p) {
        if (is.null(p) || length(p) == 0L || all(is.na(p))) "."
        else paste0("(", paste(p, collapse = ", "), ")")
      }, character(1L)), collapse = ", ")
    }, character(1L))
  }
  
  d
}

# Resolve the number of decimal places used by print.mfpi().
# User-supplied digits takes precedence over x$digits. If neither is available,
# use 3 for backward-compatible display.
resolve_print_digits <- function(x, digits = NULL) {
  if (is.null(digits)) {
    digits <- if (!is.null(x$digits)) x$digits else 3L
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0L || digits != floor(digits)) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }
  
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

# Select only the columns relevant to the criterion used for interaction
# selection. Full diagnostics remain available in the object itself.
select_print_metric_columns <- function(d, criterion, include_adjusted = FALSE) {
  if (is.null(d) || nrow(d) == 0L) return(d)
  
  criterion <- if (!is.null(criterion)) criterion else "pvalue"
  
  if (criterion == "aic" && !("dAIC" %in% names(d)) &&
      "AIC_main_minus_int" %in% names(d)) {
    d$dAIC <- d$AIC_main_minus_int
  }
  if (criterion == "bic" && !("dBIC" %in% names(d)) &&
      "BIC_main_minus_int" %in% names(d)) {
    d$dBIC <- d$BIC_main_minus_int
  }
  
  criterion_cols <- switch(
    criterion,
    "pvalue" = c(
      "deviance_int", "deviance_diff", "df_int", "pvalue",
      if (isTRUE(include_adjusted)) "p_adjusted"
    ),
    "aic" = c("AIC_main", "AIC_interaction", "dAIC"),
    "bic" = c("BIC_main", "BIC_interaction", "dBIC"),
    character(0L)
  )
  
  cols <- c(
    "type", "variable",
    criterion_cols,
    "fp_powers_main", "fp_powers_int"
  )
  cols <- unique(cols[cols %in% names(d)])
  
  d[, cols, drop = FALSE]
}

print_adjustment_step <- function(x, ruler, digits) {
  cat("\nStep 1 - Adjustment Model (MFP):\n")
  cat(ruler, "\n")
  
  fp <- x$adjust_terms
  if (is.null(fp)) {
    cat("  No adjustment model fitted.\n")
    return(invisible(NULL))
  }
  
  selected_rows <- fp[!is.na(fp$selected) &  fp$selected, , drop = FALSE]
  dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, , drop = FALSE]
  
  if (nrow(selected_rows) > 0L) {
    print(round_print_numeric(selected_rows, digits))
  } else {
    cat("  No adjustment variables selected.\n")
  }
  
  if (nrow(dropped_rows) > 0L) {
    cat(sprintf(
      "\n  Dropped (%d variables eliminated by MFP): %s\n",
      nrow(dropped_rows), paste(rownames(dropped_rows), collapse = ", ")
    ))
  }
  
  invisible(NULL)
}

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
    include_adjusted = FALSE
  )
  
  print(round_print_numeric(as.data.frame(all_m), digits), row.names = FALSE)
  invisible(NULL)
}

print_candidate_pvalue_step <- function(x, ruler, digits) {
  cat("\nStep 3 - P-value Adjustment:\n")
  cat(ruler, "\n")
  
  all_m <- x$all_model_metrics
  if (is.null(all_m) || nrow(all_m) == 0L) {
    cat("  No candidate p-values available.\n")
    return(invisible(NULL))
  }
  
  method <- if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none"
  type_label <- format_type_label(all_m$type)
  padj <- if ("p_adjusted" %in% names(all_m)) all_m$p_adjusted else all_m$pvalue
  
  tab <- data.frame(
    Variable   = all_m$variable,
    Type       = type_label,
    p_raw      = all_m$pvalue,
    p_adjusted = padj,
    check.names = FALSE
  )
  
  print(round_print_numeric(tab, digits), row.names = FALSE)
  cat(sprintf("\n  p_adjust_method = %s\n", method))
  
  invisible(NULL)
}

print_interaction_summary_step <- function(x, ruler, digits, step_no = 3L) {
  cat(sprintf("\nStep %d - Interaction Summary:\n", step_no))
  cat(ruler, "\n")
  
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  vw <- x$var_winners
  best_m <- x$best_model_metrics
  selected_vars <- if (!is.null(best_m) && nrow(best_m) > 0L) {
    best_m$variable
  } else {
    character(0L)
  }
  
  if (is.null(vw) || length(vw) == 0L) {
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      best_m <- format_mfpi_powers(best_m)
      best_m <- select_print_metric_columns(best_m, criterion = crit)
      print(round_print_numeric(as.data.frame(best_m), digits), row.names = FALSE)
    } else {
      cat("  No interactions selected.\n")
    }
    return(invisible(NULL))
  }
  
  cont_names <- names(vw)
  
  if (crit == "pvalue") {
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(
          Variable = vn,
          Type = "---",
          p_raw = NA_real_,
          p_adjusted = NA_real_,
          Selected = "",
          check.names = FALSE
        )
      } else {
        p_adj <- if ("p_adjusted" %in% names(w$metric)) {
          w$metric$p_adjusted[1L]
        } else {
          w$metric$pvalue[1L]
        }
        
        data.frame(
          Variable = vn,
          Type = format_type_label(w$type),
          p_raw = w$metric$pvalue[1L],
          p_adjusted = p_adj,
          Selected = if (vn %in% selected_vars) "*" else "",
          check.names = FALSE
        )
      }
    }))
    
    print(round_print_numeric(tab, digits), row.names = FALSE)
    cat(sprintf(
      "\n  * = selected at p_interact = %s\n",
      format_print_number(x$p_interact, digits)
    ))
  } else {
    ic_col   <- if (crit == "aic") "dAIC" else "dBIC"
    ic_label <- if (crit == "aic") "dAIC" else "dBIC"
    
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(
          Variable = vn,
          Type = "---",
          score = NA_real_,
          Selected = "",
          check.names = FALSE
        )
      } else {
        val <- if (ic_col %in% names(w$metric)) w$metric[[ic_col]][1L] else NA_real_
        data.frame(
          Variable = vn,
          Type = format_type_label(w$type),
          score = val,
          Selected = if (vn %in% selected_vars) "*" else "",
          check.names = FALSE
        )
      }
    }))
    
    names(tab)[names(tab) == "score"] <- ic_label
    print(round_print_numeric(tab, digits), row.names = FALSE)
    
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
# print.mfpi() ----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print an \code{"mfpi"} Object
#'
#' Displays the adjustment model, criterion-specific interaction model metrics,
#' and the final interaction decision summary. The printed candidate table shows
#' only the metrics relevant to the selected criterion. Full diagnostics remain
#' available in the returned object.
#'
#' For \code{criterion = "pvalue"}, an additional p-value adjustment step is
#' printed before the final summary so that raw and adjusted p-values can be
#' inspected separately.
#'
#' Numeric values are rounded for display only. The underlying metrics stored in
#' the \code{"mfpi"} object are not modified, and model-selection decisions are
#' based on the unrounded values computed during fitting.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param digits Optional non-negative integer controlling the number of decimal
#'   places used when printing numeric output. If \code{NULL}, \code{x$digits}
#'   is used when available; otherwise 3 is used.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print mfpi
#' @export
print.mfpi <- function(x, digits = NULL, ...) {
  ruler  <- strrep("-", 85)
  header <- strrep("=", 85)
  padj   <- if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none"
  crit   <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  digits <- resolve_print_digits(x, digits)
  
  cat(header, "\n")
  cat(sprintf(
    "MFPI  |  group: '%s'  |  n = %d  |  %s  |  criterion: %s",
    x$group_var, x$nobs, x$flex, crit
  ))
  if (crit == "pvalue") cat(sprintf("  |  p-adjust: %s", padj))
  cat("\n")
  cat(header, "\n")
  
  ruler2  <- strrep("-", 37)
  print_adjustment_step(x, ruler2, digits)
  print_candidates_step(x, ruler2, digits)
  
  if (crit == "pvalue") {
    print_candidate_pvalue_step(x, ruler2, digits)
    print_interaction_summary_step(x, ruler2, digits, step_no = 4L)
  } else {
    print_interaction_summary_step(x, ruler2, digits, step_no = 3L)
  }
  
  invisible(x)
}