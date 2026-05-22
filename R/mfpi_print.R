# -----------------------------------------------------------------------------
# print.mfpi() ----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print an \code{mfpi} Object
#'
#' Displays a structured summary of an \code{"mfpi"} object, including the
#' adjustment variables selected by the MFP algorithm and the interaction test
#' results for each continuous variable in \code{cont_vars}. Optionally prints
#' the full regression output for each final interaction model when
#' \code{show_models = TRUE} was set in \code{mfpi()}.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @section Output sections:
#' \enumerate{
#'   \item **Adjustment model** — the \code{fp_terms} table from the MFP
#'     algorithm, showing which variables were selected and their estimated FP
#'     powers.
#'   \item **Interaction test** — for each variable in \code{cont_vars}, the
#'     best-selected functional form (linear, FP1, or FP2) and its evaluation
#'     metrics: deviance, degrees of freedom, p-value, AIC, and BIC. Columns:
#'     \itemize{
#'       \item \code{type}: selected functional form (\code{linear}, \code{fp1},
#'         or \code{fp2}).
#'       \item \code{fp_powers_main}: FP powers in the main-effects model.
#'       \item \code{fp_powers_int}: FP powers in the interaction model, one set per
#'         group. For \code{flex4} these may differ across groups.
#'       \item \code{deviance_diff}: likelihood-ratio test statistic
#'         \eqn{T = -2\ell_\text{main} - (-2\ell_\text{int})}.
#'       \item \code{df_interaction}: degrees of freedom for the LRT
#'         (\eqn{= (K-1)m} for flex1/2/3; \eqn{= 2(K-1)m} for flex4).
#'       \item \code{pvalue}: \eqn{\Pr[\chi^2(\text{df\_interaction}) > T]}.
#'       \item \code{AIC_main_minus_int}, \code{BIC_main_minus_int}: improvement of interaction
#'         model over main-effects model; positive values favour interaction.
#'     }
#'   \item **Interaction model coefficients** (only when
#'     \code{show_models = TRUE}) — the regression summary for the final
#'     interaction model for each variable.
#' }
#'
#' @method print mfpi
#' @export
print.mfpi <- function(x, ...) {
  
  # ---------------------------------------------------------------------------
  # Section 1: Adjustment model
  # ---------------------------------------------------------------------------
  cat("Adjustment Variables Selected by MFP Algorithm:\n")
  cat(strrep("-", 65), "\n")
  
  fp <- x$adjust_terms
  if (!is.null(fp)) {
    # Show all fp_terms columns for selected rows only.
    # Dropped variables are summarised in a compact one-line list.
    selected_rows <- fp[!is.na(fp$selected) &  fp$selected, , drop = FALSE]
    dropped_rows  <- fp[!is.na(fp$selected) & !fp$selected, , drop = FALSE]
    
    if (nrow(selected_rows) > 0L) {
      print(selected_rows)
    } else {
      cat("No adjustment variables selected.\n")
    }
    if (nrow(dropped_rows) > 0L) {
      cat(sprintf(
        "\nDropped (%d variables eliminated by MFP): %s\n",
        nrow(dropped_rows),
        paste(rownames(dropped_rows), collapse = ", ")
      ))
    }
  } else {
    cat("No adjustment model fitted.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Section 2: Interaction test results
  # ---------------------------------------------------------------------------
  cat("\n", strrep("-", 65), "\n", sep = "")
  cat("Interaction Test\n")
  cat(strrep("-", 65), "\n")
  cat(sprintf(
    "\nInteraction with '%s'  |  %d observations  |  %s strategy\n\n",
    x$group_var, x$nobs, x$flex
  ))
  
  print(x$model_evaluation_metrics)
  
  cat(
    "\nNote: df_interaction = LRT degrees of freedom for interaction;",
    "df_total = total df in interaction model.\n"
  )
  
  # For flex4, each group has its own FP powers so fp_powers_int shows K power sets
  if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
    cat(sprintf(
      "\nFor flex4, 'fp_powers_int' shows separate FP powers for each level of '%s': %s.\n",
      x$group_var,
      paste(x$group_levels_original, collapse = ", ")
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Section 3: Interaction model summaries (optional)
  # ---------------------------------------------------------------------------
  if (isTRUE(x$show_models) && length(x$interaction_models) > 0L) {
    
    cat("\n", strrep("-", 65), "\n", sep = "")
    cat("Regression Output for Final Interaction Models:\n")
    cat(strrep("-", 65), "\n")
    
    # model_evaluation_metrics has one row per variable with a 'type' column
    # (the selected functional form) — use it to label the output
    metrics <- x$model_evaluation_metrics
    type_lookup <- if (!is.null(metrics) && "type" %in% names(metrics)) {
      setNames(
        toupper(gsub("fp", "FP", metrics$type)),  # "linear"->"LINEAR" / "fp1"->"FP1"
        metrics$variable
      )
    } else {
      setNames(rep("", length(x$interaction_models)),
               names(x$interaction_models))
    }
    
    for (var_name in names(x$interaction_models)) {
      fit_obj  <- x$interaction_models[[var_name]]
      form_lbl <- if (!is.null(type_lookup[[var_name]])) type_lookup[[var_name]] else ""
      
      cat(sprintf(
        "\n'%s' x '%s'  [%s]\n",
        x$group_var, var_name, form_lbl
      ))
      cat(strrep("-", 45), "\n")
      
      # fit_obj is the object returned by mfp2:::fit_model (fast = FALSE).
      # It carries a $fit component holding the raw glm / coxph object.
      if (!is.null(fit_obj$fit)) {
        print(summary(fit_obj$fit))
      } else {
        # Fallback: print whatever is available
        print(fit_obj)
      }
    }
  }
  
  invisible(x)
}



# S3 summary method for mfpi objects
#
# summary.mfpi() returns a structured list containing the key results from an
# mfpi fit. Unlike print.mfpi(), which formats output for the console,
# summary.mfpi() returns an object that can be inspected programmatically.
# A print method for the returned object formats it for display.
#
# Note on degrees of freedom: the p-values in the regression summary tables
# produced by summary.glm() and summary.coxph() do NOT account for the
# degrees of freedom consumed by estimating FP powers. They should be treated
# as approximate. The correct interaction test p-values, with proper df, are
# in the $model_evaluation_metrics component.


# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarise an \code{"mfpi"} Object
#'
#' Produces a structured summary of an \code{"mfpi"} model fit, including the
#' adjustment model, interaction test metrics, and the regression output for
#' each final interaction model. The returned object has its own print method
#' for formatted console display.
#'
#' @param object An object of class \code{"mfpi"}, as returned by
#'   \code{mfpi()}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.mfpi"}, which is a named list
#'   with the following components:
#' \describe{
#'   \item{\code{group_var}}{Name of the grouping variable.}
#'   \item{\code{nobs}}{Number of observations.}
#'   \item{\code{family}}{Regression family (GLM object or \code{"cox"}).}
#'   \item{\code{flex}}{Flexibility level used (\code{"flex1"} to
#'     \code{"flex4"}).}
#'   \item{\code{group_levels_original}}{Original levels of \code{group_var}.}
#'   \item{\code{adjust_terms}}{The full \code{fp_terms} data frame from the
#'     MFP adjustment model, as returned by \code{mfp2:::fit_mfp()}.}
#'   \item{\code{model_evaluation_metrics}}{Tibble of interaction test results
#'     (one row per significant variable): functional form selected, FP powers,
#'     deviance, df, p-value, AIC, and BIC. These p-values correctly account
#'     for the degrees of freedom of the interaction test.}
#'   \item{\code{all_model_metrics}}{Tibble of metrics for all three candidate
#'     models (linear, FP1, FP2) for every variable tested, regardless of
#'     significance.}
#'   \item{\code{model_summaries}}{Named list of regression summaries, one per
#'     variable in \code{cont_vars} that had a significant interaction. Each
#'     element is the output of \code{summary(fit$fit)}, where \code{fit} is
#'     the raw \code{glm} or \code{coxph} object. \strong{Note:} p-values in
#'     these summaries are from the standard \code{glm}/\code{coxph} machinery
#'     and do not account for the degrees of freedom used to estimate FP
#'     powers; use \code{model_evaluation_metrics$pvalue} for the correct
#'     interaction test p-values.}
#' }
#'
#' @seealso \code{print.summary.mfpi()}, \code{mfpi()},
#'   [stats::summary.glm()], [survival::summary.coxph()]
#'
#' @method summary mfpi
#' @export
summary.mfpi <- function(object, ...) {
  
  # Summarise each final interaction model (flat named list by var_name)
  model_summaries <- lapply(object$interaction_models, function(fit_obj) {
    if (!is.null(fit_obj$fit)) {
      summary(fit_obj$fit)
    } else {
      NULL
    }
  })
  
  structure(
    list(
      group_var                = object$group_var,
      nobs                     = object$nobs,
      family                   = object$family,
      flex                     = object$flex,
      group_levels_original    = object$group_levels_original,
      adjust_terms             = object$adjust_terms,
      model_evaluation_metrics = object$model_evaluation_metrics,
      all_model_metrics        = object$all_model_metrics,
      model_summaries          = model_summaries
    ),
    class = "summary.mfpi"
  )
}


# -----------------------------------------------------------------------------
# print.summary.mfpi() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print a \code{"summary.mfpi"} Object
#'
#' Formats and prints the structured summary produced by
#' \code{summary.mfpi()}.
#'
#' @param x An object of class \code{"summary.mfpi"}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print summary.mfpi
#' @export
print.summary.mfpi <- function(x, ...) {
  
  cat(strrep("=", 65), "\n")
  cat(sprintf(
    "MFPI Summary  |  group: '%s'  |  n = %d  |  %s\n",
    x$group_var, x$nobs, x$flex
  ))
  cat(strrep("=", 65), "\n")
  
  # --- Adjustment model ------------------------------------------------------
  cat("\nAdjustment Variables (MFP):\n")
  cat(strrep("-", 65), "\n")
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
        "\n  Dropped (%d): %s\n",
        nrow(dropped_rows),
        paste(rownames(dropped_rows), collapse = ", ")
      ))
    }
  } else {
    cat("  No adjustment model fitted.\n")
  }
  
  # --- Interaction test results ----------------------------------------------
  cat("\n", strrep("-", 65), "\n", sep = "")
  cat("Interaction Test Results:\n")
  cat(strrep("-", 65), "\n\n")
  print(x$model_evaluation_metrics)
  cat("\nNote: df_interaction = LRT df for interaction; df_total = total df in interaction model.\n")
  cat("      p-values correctly account for interaction test degrees of freedom.\n")
  
  if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
    cat(sprintf(
      "\nFor flex4, 'fp_powers_int' shows per-group FP powers (levels: %s).\n",
      paste(x$group_levels_original, collapse = ", ")
    ))
  }
  
  # --- Regression summaries --------------------------------------------------
  if (length(x$model_summaries) > 0L) {
    cat("\n", strrep("-", 65), "\n", sep = "")
    cat("Regression Output for Final Interaction Models:\n")
    cat("(Note: p-values below are from glm/coxph and do not account for\n")
    cat(" FP power estimation df. Use 'model_evaluation_metrics$pvalue'\n")
    cat(" for the correct interaction test p-values.)\n")
    cat(strrep("-", 65), "\n")
    
    metrics <- x$model_evaluation_metrics
    type_lookup <- if (!is.null(metrics) && "type" %in% names(metrics)) {
      setNames(toupper(gsub("fp", "FP", metrics$type)), metrics$variable)
    } else {
      setNames(rep("", length(x$model_summaries)), names(x$model_summaries))
    }
    
    for (var_name in names(x$model_summaries)) {
      sm       <- x$model_summaries[[var_name]]
      form_lbl <- if (!is.null(type_lookup[[var_name]])) type_lookup[[var_name]] else ""
      cat(sprintf("\n'%s' x '%s'  [%s]\n", x$group_var, var_name, form_lbl))
      cat(strrep("-", 45), "\n")
      if (!is.null(sm)) print(sm) else cat("  (model summary unavailable)\n")
    }
  }
  
  invisible(x)
}