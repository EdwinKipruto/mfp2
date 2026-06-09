# S3 print method for mfpi objects
#
# print.mfpi() gives a structured output:
#   Step 1 - Adjustment model (fp_terms for selected variables)
#   Step 2 - All interaction models with their metrics
#   Step 3 - P-value adjustment when criterion = "pvalue"
#   Step 4 - Interaction summary when criterion = "pvalue"
# For AIC/BIC, Step 3 is the interaction summary.

# -----------------------------------------------------------------------------
# Helpers ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

format_mfpi_powers <- function(d) {
  if ("fp_powers_main" %in% names(d)) {
    d$fp_powers_main <- vapply(d$fp_powers_main, function(p) {
      if (is.null(p) || length(p) == 0L || all(is.na(p))) "."
      else paste0("(", paste(p, collapse = ", "), ")")
    }, character(1L))
  }
  if ("fp_powers_int" %in% names(d)) {
    d$fp_powers_int <- vapply(d$fp_powers_int, function(p_list) {
      if (is.null(p_list) || length(p_list) == 0L) return(".")
      if (!is.list(p_list)) return(paste0("(", paste(p_list, collapse = ", "), ")"))
      paste(vapply(p_list, function(p) {
        if (is.null(p) || length(p) == 0L || all(is.na(p))) "."
        else paste0("(", paste(p, collapse = ", "), ")")
      }, character(1L)), collapse = ", ")
    }, character(1L))
  }
  d
}

print_adjustment_step <- function(x, ruler) {
  cat("\nStep 1 - Adjustment Model (MFP):\n")
  cat(ruler, "\n")
  fp <- x$adjust_terms
  if (!is.null(fp)) {
    selected_rows <- fp[!is.na(fp$selected) &  fp$selected, , drop = FALSE]
    dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, , drop = FALSE]
    if (nrow(selected_rows) > 0L) print(selected_rows) else cat("  No adjustment variables selected.\n")
    if (nrow(dropped_rows) > 0L) {
      cat(sprintf("\n  Dropped (%d variables eliminated by MFP): %s\n",
                  nrow(dropped_rows), paste(rownames(dropped_rows), collapse = ", ")))
    }
  } else {
    cat("  No adjustment model fitted.\n")
  }
}

print_candidates_step <- function(x, ruler) {
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")
  all_m <- x$all_model_metrics
  if (is.null(all_m) || nrow(all_m) == 0L) {
    cat("  No interaction candidates computed.\n")
    return(invisible(NULL))
  }
  all_m <- format_mfpi_powers(all_m)
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  if (crit == "pvalue") {
    all_m <- all_m[, setdiff(names(all_m), "p_adjusted"), drop = FALSE]
  }
  print(as.data.frame(all_m), row.names = FALSE)
  invisible(NULL)
}

print_candidate_pvalue_step <- function(x, ruler) {
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
  print(tab, row.names = FALSE)
  cat(sprintf("\n  p_adjust_method = %s\n", method))
  invisible(NULL)
}

print_interaction_summary_step <- function(x, ruler, step_no = 3L) {
  cat(sprintf("\nStep %d - Interaction Summary:\n", step_no))
  cat(ruler, "\n")
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  vw <- x$var_winners
  best_m <- x$best_model_metrics
  selected_vars <- if (!is.null(best_m) && nrow(best_m) > 0L) best_m$variable else character(0L)
  if (is.null(vw) || length(vw) == 0L) {
    if (!is.null(best_m) && nrow(best_m) > 0L) print(best_m, row.names = FALSE) else cat("  No interactions selected.\n")
    return(invisible(NULL))
  }
  cont_names <- names(vw)
  if (crit == "pvalue") {
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(Variable = vn, Type = "---", p_raw = NA_real_, p_adjusted = NA_real_, Selected = "", check.names = FALSE)
      } else {
        data.frame(Variable = vn,
                   Type = format_type_label(w$type),
                   p_raw = w$metric$pvalue[1L],
                   p_adjusted = if ("p_adjusted" %in% names(w$metric)) w$metric$p_adjusted[1L] else w$metric$pvalue[1L],
                   Selected = if (vn %in% selected_vars) "*" else "",
                   check.names = FALSE)
      }
    }))
    print(tab, row.names = FALSE)
    cat(sprintf("\n  * = selected at p_interact = %g\n", x$p_interact))
  } else {
    ic_col   <- if (crit == "aic") "dAIC" else "dBIC"
    ic_label <- if (crit == "aic") "dAIC" else "dBIC"
    tab <- do.call(rbind, lapply(cont_names, function(vn) {
      w <- vw[[vn]]
      if (is.null(w$fit) || is.null(w$metric)) {
        data.frame(Variable = vn, Type = "---", score = NA_real_, Selected = "", check.names = FALSE)
      } else {
        val <- if (ic_col %in% names(w$metric)) w$metric[[ic_col]][1L] else NA_real_
        data.frame(Variable = vn,
                   Type = format_type_label(w$type),
                   score = val,
                   Selected = if (vn %in% selected_vars) "*" else "",
                   check.names = FALSE)
      }
    }))
    names(tab)[names(tab) == "score"] <- ic_label
    print(tab, row.names = FALSE)
    min_imp <- if (!is.null(x$min_improvement)) x$min_improvement else 2
    cat(sprintf("\n  * = selected (%s > %g)\n", ic_label, min_imp))
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# print.mfpi() ----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print an \code{"mfpi"} Object
#'
#' Displays the adjustment model, interaction model metrics, and the final
#' interaction decision summary. For \code{criterion = "pvalue"}, an additional
#' p-value adjustment step is printed before the final summary so that raw and
#' adjusted p-values can be inspected separately.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print mfpi
#' @export
print.mfpi <- function(x, ...) {
  ruler  <- strrep("-", 65)
  header <- strrep("=", 65)
  padj   <- if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none"
  crit   <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  
  cat(header, "\n")
  cat(sprintf("MFPI  |  group: '%s'  |  n = %d  |  %s  |  criterion: %s",
              x$group_var, x$nobs, x$flex, crit))
  if (crit == "pvalue") cat(sprintf("  |  p-adjust: %s", padj))
  cat("\n")
  cat(header, "\n")
  
  print_adjustment_step(x, ruler)
  print_candidates_step(x, ruler)
  if (crit == "pvalue") {
    print_candidate_pvalue_step(x, ruler)
    print_interaction_summary_step(x, ruler, step_no = 4L)
  } else {
    print_interaction_summary_step(x, ruler, step_no = 3L)
  }
  invisible(x)
}