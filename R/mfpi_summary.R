# S3 summary method for mfpi objects
#
# summary.mfpi()       - collects all results into a "summary.mfpi" list
# print.summary.mfpi() - four-step formatted display:
#   Step 1 - Adjustment model
#   Step 2 - All interaction candidates (best-per-variable flagged with +)
#   Step 3 - Interaction summary (all cont_vars, * for selected)
#   Step 4 - Regression output for winning models


# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarise an \code{"mfpi"} Object
#'
#' Returns a structured summary of an \code{"mfpi"} fit. The returned object
#' has a dedicated print method that displays four sections:
#' \enumerate{
#'   \item Adjustment model selected by MFP.
#'   \item All interaction candidates (linear, FP1, FP2) with the best
#'     candidate per variable flagged.
#'   \item Interaction summary for all tested variables, with p-values
#'     (raw and adjusted when \code{p_adjust_method != "none"}) or
#'     information criterion differences.
#'   \item Regression output (coefficients) for each winning model.
#' }
#'
#' @param object An object of class \code{"mfpi"}, as returned by
#'   \code{mfpi()}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.mfpi"} - a named list with:
#' \describe{
#'   \item{\code{group_var}}{Name of the grouping variable.}
#'   \item{\code{nobs}}{Number of observations.}
#'   \item{\code{family}}{Regression family.}
#'   \item{\code{flex}}{Flexibility level used.}
#'   \item{\code{criterion}}{Selection criterion used.}
#'   \item{\code{p_adjust_method}}{Method used for p-value adjustment.}
#'   \item{\code{p_interact}}{Significance threshold for interaction selection.}
#'   \item{\code{min_improvement}}{Threshold for AIC/BIC selection.}
#'   \item{\code{group_levels_original}}{Original levels of \code{group_var}.}
#'   \item{\code{adjust_terms}}{Full \code{fp_terms} table from the MFP
#'     adjustment model.}
#'   \item{\code{all_model_metrics}}{Metrics for all candidates.}
#'   \item{\code{best_model_metrics}}{Metrics for selected (winning) models.}
#'   \item{\code{var_winners}}{Named list of per-variable best candidates
#'     (regardless of significance).}
#'   \item{\code{best_interaction_model}}{Named list of fitted interaction
#'     model objects, one per significant variable.}
#'   \item{\code{all_interaction_models}}{Named list of all candidate
#'     interaction model objects.}
#'   \item{\code{best_fitted_functions}}{Named list of fitted-function
#'     matrices for the winning model per significant variable.}
#'   \item{\code{all_fitted_functions}}{Named list of fitted-function
#'     matrices for all candidates.}
#'   \item{\code{model_summaries}}{Named list of \code{summary(fit$fit)}
#'     objects, one per significant variable.}
#' }
#'
#' @seealso \code{print.summary.mfpi()}, \code{mfpi()},
#'   [stats::summary.glm()], [survival::summary.coxph()]
#'
#' @method summary mfpi
#' @export
summary.mfpi <- function(object, ...) {
  
  model_summaries <- lapply(object$best_interaction_model, function(fit_obj) {
    if (!is.null(fit_obj$fit)) summary(fit_obj$fit) else NULL
  })
  
  structure(
    list(
      group_var              = object$group_var,
      nobs                   = object$nobs,
      family                 = object$family,
      flex                   = object$flex,
      criterion              = object$criterion,
      p_adjust_method        = object$p_adjust_method,
      p_interact             = object$p_interact,
      min_improvement        = object$min_improvement,
      group_levels_original  = object$group_levels_original,
      adjust_terms           = object$adjust_terms,
      all_model_metrics      = object$all_model_metrics,
      best_model_metrics     = object$best_model_metrics,
      var_winners            = object$var_winners,
      best_interaction_model = object$best_interaction_model,
      all_interaction_models = object$all_interaction_models,
      best_fitted_functions  = object$best_fitted_functions,
      all_fitted_functions   = object$all_fitted_functions,
      model_summaries        = model_summaries
    ),
    class = "summary.mfpi"
  )
}


# -----------------------------------------------------------------------------
# print.summary.mfpi() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print a \code{"summary.mfpi"} Object
#'
#' Displays a four-step structured summary. Steps 1-3 match
#' \code{print.mfpi()}; Step 4 adds regression output for each winning
#' interaction model.
#'
#' @param x An object of class \code{"summary.mfpi"}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print summary.mfpi
#' @export
print.summary.mfpi <- function(x, ...) {
  
  ruler  <- strrep("-", 65)
  header <- strrep("=", 65)
  padj   <- if (!is.null(x$p_adjust_method)) x$p_adjust_method else "none"
  crit   <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  adjusting <- padj != "none" && crit == "pvalue"
  
  cat(header, "\n")
  cat(sprintf(
    "MFPI Summary  |  group: '%s'  |  n = %d  |  %s  |  p-adjust: %s\n",
    x$group_var, x$nobs, x$flex, padj
  ))
  cat(header, "\n")
  
  # ---------------------------------------------------------------------------
  # Step 1: Adjustment model
  # ---------------------------------------------------------------------------
  cat("\nStep 1 - Adjustment Model (MFP):\n")
  cat(ruler, "\n")
  
  fp <- x$adjust_terms
  if (!is.null(fp)) {
    selected_rows <- fp[!is.na(fp$selected) &  fp$selected, , drop = FALSE]
    dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, , drop = FALSE]
    
    if (nrow(selected_rows) > 0L) {
      print(selected_rows)
    } else {
      cat("  No adjustment variables selected.\n")
    }
    if (nrow(dropped_rows) > 0L) {
      cat(sprintf(
        "\n  Dropped (%d variables eliminated by MFP): %s\n",
        nrow(dropped_rows),
        paste(rownames(dropped_rows), collapse = ", ")
      ))
    }
  } else {
    cat("  No adjustment model fitted.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 2: All interaction candidates with best-per-variable flagged
  # ---------------------------------------------------------------------------
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")
  
  all_m <- x$all_model_metrics
  vw    <- x$var_winners
  
  if (!is.null(all_m) && nrow(all_m) > 0L) {
    
    if (!is.null(vw)) {
      winner_keys <- vapply(names(vw), function(vn) {
        wtype <- vw[[vn]]$type
        if (is.null(wtype)) NA_character_
        else paste(wtype, vn, sep = "__")
      }, character(1L))
      winner_keys <- winner_keys[!is.na(winner_keys)]
    } else {
      best_m <- x$best_model_metrics
      if (!is.null(best_m) && nrow(best_m) > 0L) {
        winner_keys <- paste(best_m$type, best_m$variable, sep = "__")
      } else {
        winner_keys <- character(0L)
      }
    }
    
    if ("fp_powers_main" %in% names(all_m)) {
      all_m$fp_powers_main <- vapply(all_m$fp_powers_main, function(p)
        paste0("(", paste(p, collapse = ", "), ")"), character(1L))
    }
    if ("fp_powers_int" %in% names(all_m)) {
      all_m$fp_powers_int <- vapply(all_m$fp_powers_int, function(p_list)
        paste(vapply(p_list, function(p)
          paste0("(", paste(p, collapse = ", "), ")"),
          character(1L)), collapse = ", "), character(1L))
    }
    
    row_keys    <- paste(all_m$type, all_m$variable, sep = "__")
    all_m$best  <- ifelse(row_keys %in% winner_keys, "+", "")
    
    # Insert best after the last criterion-specific column, before powers
    insert_after <- switch(crit,
                           pvalue = "BIC_main_minus_int",
                           aic    = "AIC_main_minus_int",
                           bic    = "BIC_main_minus_int"
    )
    cols <- names(all_m)
    best_pos <- which(cols == "best")
    after_pos <- which(cols == insert_after)
    if (length(after_pos) == 1L && length(best_pos) == 1L) {
      cols_no_best <- cols[-best_pos]
      insert_at <- which(cols_no_best == insert_after)
      all_m <- all_m[, append(cols_no_best, "best", after = insert_at),
                     drop = FALSE]
    }
    
    print(as.data.frame(all_m), row.names = FALSE)
    
    if (crit == "pvalue") {
      cat("\n  + = best candidate (by pvalue)\n")
    }
    
  } else {
    cat("  No interaction candidates computed.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: Interaction summary (all cont_vars)
  # ---------------------------------------------------------------------------
  cat("\nStep 3 - Interaction Summary:\n")
  cat(ruler, "\n")
  
  best_m <- x$best_model_metrics
  selected_vars <- if (!is.null(best_m) && nrow(best_m) > 0L)
    best_m$variable else character(0L)
  
  if (!is.null(vw) && length(vw) > 0L) {
    cont_names <- names(vw)
    raw_pvals  <- vapply(vw, function(w) {
      if (is.null(w$fit)) NA_real_ else w$metric$pvalue[1L]
    }, numeric(1L))
    
    if (crit == "pvalue") {
      if (adjusting) {
        adj_pvals <- stats::p.adjust(raw_pvals, method = padj)
        cat(sprintf("  %-12s %-8s %12s %12s  %s\n",
                    "Variable", "Type", "p_raw", "p_adjusted", ""))
        cat(strrep("-", 55), "\n")
        for (i in seq_along(cont_names)) {
          vn <- cont_names[i]
          w  <- vw[[vn]]
          if (is.null(w$fit)) {
            cat(sprintf("  %-12s %-8s %12s %12s\n", vn, "---", "NA", "NA"))
          } else {
            type_label <- toupper(gsub("fp", "FP", w$type))
            sel_flag   <- if (vn %in% selected_vars) " *" else ""
            cat(sprintf("  %-12s %-8s %12s %12s%s\n",
                        vn, type_label,
                        formatC(raw_pvals[i], format = "g", digits = 4),
                        formatC(adj_pvals[i], format = "g", digits = 4),
                        sel_flag))
          }
        }
        cat(strrep("-", 55), "\n")
        cat(sprintf("  * = selected at p_interact = %g\n", x$p_interact))
      } else {
        cat(sprintf("  %-12s %-8s %12s  %s\n",
                    "Variable", "Type", "pvalue", ""))
        cat(strrep("-", 42), "\n")
        for (i in seq_along(cont_names)) {
          vn <- cont_names[i]
          w  <- vw[[vn]]
          if (is.null(w$fit)) {
            cat(sprintf("  %-12s %-8s %12s\n", vn, "---", "NA"))
          } else {
            type_label <- toupper(gsub("fp", "FP", w$type))
            sel_flag   <- if (vn %in% selected_vars) " *" else ""
            cat(sprintf("  %-12s %-8s %12s%s\n",
                        vn, type_label,
                        formatC(raw_pvals[i], format = "g", digits = 4),
                        sel_flag))
          }
        }
        cat(strrep("-", 42), "\n")
        cat(sprintf("  * = selected at p_interact = %g\n", x$p_interact))
      }
    } else {
      ic_col   <- if (crit == "aic") "AIC_main_minus_int" else "BIC_main_minus_int"
      ic_label <- if (crit == "aic") "dAIC" else "dBIC"
      min_imp  <- if (!is.null(x$min_improvement)) x$min_improvement else 2
      scores   <- vapply(vw, function(w) {
        if (is.null(w$fit)) NA_real_ else w$metric[[ic_col]][1L]
      }, numeric(1L))
      
      cat(sprintf("  %-12s %-8s %12s  %s\n",
                  "Variable", "Type", ic_label, ""))
      cat(strrep("-", 42), "\n")
      for (i in seq_along(cont_names)) {
        vn <- cont_names[i]
        w  <- vw[[vn]]
        if (is.null(w$fit)) {
          cat(sprintf("  %-12s %-8s %12s\n", vn, "---", "NA"))
        } else {
          type_label <- toupper(gsub("fp", "FP", w$type))
          sel_flag   <- if (vn %in% selected_vars) " *" else ""
          cat(sprintf("  %-12s %-8s %12s%s\n",
                      vn, type_label,
                      formatC(scores[i], format = "f", digits = 2),
                      sel_flag))
        }
      }
      cat(strrep("-", 42), "\n")
      cat(sprintf("  * = selected (%s > %g)\n", ic_label, min_imp))
    }
    
  } else {
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      disp <- best_m
      if ("fp_powers_main" %in% names(disp)) {
        disp$fp_powers_main <- vapply(disp$fp_powers_main, function(p)
          paste0("(", paste(p, collapse = ", "), ")"), character(1L))
      }
      if ("fp_powers_int" %in% names(disp)) {
        disp$fp_powers_int <- vapply(disp$fp_powers_int, function(p_list)
          paste(vapply(p_list, function(p)
            paste0("(", paste(p, collapse = ", "), ")"),
            character(1L)), collapse = ", "), character(1L))
      }
      print(as.data.frame(disp), row.names = FALSE)
    } else {
      cat("  No significant interactions found.\n")
    }
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Regression output for winning models
  # ---------------------------------------------------------------------------
  if (length(x$model_summaries) > 0L) {
    
    cat("\nStep 4 - Regression Output (Winning Interaction Models):\n")
    cat(ruler, "\n")
    cat("  Note: p-values below are from glm/coxph and do not account for\n")
    cat("        FP power estimation df. Use Step 3 p-values for the correct\n")
    cat("        interaction test results.\n")
    
    for (var_name in names(x$model_summaries)) {
      sm <- x$model_summaries[[var_name]]
      form_lbl <- ""
      if (!is.null(vw) && !is.null(vw[[var_name]]$type)) {
        form_lbl <- toupper(gsub("fp", "FP", vw[[var_name]]$type))
      } else if (!is.null(best_m) && "type" %in% names(best_m)) {
        idx <- which(best_m$variable == var_name)
        if (length(idx) > 0L)
          form_lbl <- toupper(gsub("fp", "FP", best_m$type[idx[1L]]))
      }
      cat(sprintf("\n  '%s' x '%s'  [%s]\n", x$group_var, var_name, form_lbl))
      cat(strrep("-", 45), "\n")
      if (!is.null(sm)) print(sm) else cat("  (model summary unavailable)\n")
    }
  }
  
  invisible(x)
}