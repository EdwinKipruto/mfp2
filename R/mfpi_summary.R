# S3 summary method for mfpi objects
#
# summary.mfpi() returns all printed components plus regression summaries for
# retained interaction models. print.summary.mfpi() reuses print.mfpi() for the
# model-building output and appends the regression output as the final step.

# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarize an \code{"mfpi"} Object
#'
#' Produces a structured summary of an \code{"mfpi"} fit. The printed summary
#' contains the same model-building output as \code{print.mfpi()} and adds
#' regression output for retained interaction models.
#'
#' @param object An object of class \code{"mfpi"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.mfpi"}.
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
#' Displays the structured MFPI output and appends regression summaries for
#' retained interaction models. For p-value selection, the regression output is
#' printed after Step 4; for AIC/BIC selection, it is printed after Step 3.
#'
#' @param x An object of class \code{"summary.mfpi"}.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print summary.mfpi
#' @export
print.summary.mfpi <- function(x, ...) {
  obj <- x
  class(obj) <- "mfpi"
  print.mfpi(obj)
  
  crit <- if (!is.null(x$criterion)) x$criterion else "pvalue"
  step_no <- if (crit == "pvalue") 5L else 4L
  ruler <- strrep("-", 65)
  
  cat(sprintf("\nStep %d - Regression Output for Retained Interaction Models:\n", step_no))
  cat(ruler, "\n")
  
  ms <- x$model_summaries
  if (is.null(ms) || length(ms) == 0L || all(vapply(ms, is.null, logical(1L)))) {
    cat("  No retained interaction models.\n")
    return(invisible(x))
  }
  
  for (vn in names(ms)) {
    sm <- ms[[vn]]
    if (is.null(sm)) next
    w <- if (!is.null(x$var_winners)) x$var_winners[[vn]] else NULL
    type_label <- if (!is.null(w$type)) format_type_label(w$type) else "?"
    cat(sprintf("\n  Interaction model (%s, %s):\n", vn, type_label))
    cat(strrep("-", 50), "\n")
    print(sm)
    cat("\n")
  }
  
  cat("  Note: The regression-table Wald p-values do not account for FP power\n")
  cat("        selection. Use the interaction summary p-values or IC metrics\n")
  cat("        for the MFPI interaction decision.\n")
  
  invisible(x)
}