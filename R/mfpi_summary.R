# S3 summary method for mfpi objects
#
# summary.mfpi() creates a structured summary object for an MFPI fit.
# The summary object keeps the same model-building components used by
# print.mfpi(), and adds ordinary regression summaries for retained interaction
# models.
#
# print.summary.mfpi() deliberately reuses print.mfpi() for the model-building
# output. This keeps the ordinary print output and the summary output aligned.
# The summary print method then appends retained-only details and the regression
# summaries for the final retained interaction models.
#
# Important interpretation note:
# The regression summaries printed here are summaries of the final fitted
# interaction models. Their coefficient-level Wald p-values do not account for
# the fractional-polynomial power-selection process. The MFPI interaction
# decision should therefore be based on the interaction summary table and the
# criterion-specific metrics stored in the object.

# -----------------------------------------------------------------------------
# Internal summary helpers -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Build Regression Summaries for Retained MFPI Interaction Models
#'
#' Extracts ordinary model summaries from the retained interaction models stored
#' in an \code{"mfpi"} object.
#'
#' The retained models are stored in \code{object$best_interaction_model}. Each
#' element is expected to be a fitted-model wrapper produced by the MFPI fitting
#' path. The actual fitted model is stored in its \code{fit} component. This
#' helper calls \code{summary()} on that fitted model and returns \code{NULL}
#' for missing or incomplete retained model objects.
#'
#' The returned regression summaries are intended for inspection of the final
#' retained models only. They are not used for MFPI model selection.
#'
#' @param object Object of class \code{"mfpi"}.
#'
#' @return A named list. Each element is either an ordinary model-summary object
#'   returned by \code{summary()} or \code{NULL} when no fitted model is
#'   available for that retained interaction.
#'
#' @keywords internal
#' @noRd
build_mfpi_model_summaries <- function(object) {
  retained_models <- object$best_interaction_model
  
  if (is.null(retained_models) || length(retained_models) == 0L) {
    return(list())
  }
  
  lapply(retained_models, function(fit_obj) {
    # A retained-model entry should be a list-like object with a fitted model in
    # its $fit component. If that structure is absent, there is nothing useful
    # to summarise for this variable.
    if (is.null(fit_obj) || is.null(fit_obj$fit)) {
      return(NULL)
    }
    
    # summary() dispatches to the appropriate model method, for example
    # summary.glm() or summary.coxph(). If a fitted object is malformed, return
    # NULL rather than breaking the whole summary printout.
    tryCatch(
      summary(fit_obj$fit),
      error = function(e) NULL
    )
  })
}


#' Check Whether an MFPI Summary Contains Retained Model Summaries
#'
#' Determines whether a \code{"summary.mfpi"} object contains at least one
#' non-\code{NULL} regression summary for a retained interaction model.
#'
#' @param model_summaries A list, usually \code{x$model_summaries}.
#'
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
has_retained_model_summaries <- function(model_summaries) {
  !is.null(model_summaries) &&
    length(model_summaries) > 0L &&
    any(!vapply(model_summaries, is.null, logical(1L)))
}


#' Resolve the Regression-Output Step Number for an MFPI Summary
#'
#' Computes the numbered step used for the regression-output block appended by
#' \code{print.summary.mfpi()}.
#'
#' The retained interaction-power table printed immediately before the
#' regression output is intentionally unnumbered. It is a detail table for the
#' retained models, not an additional model-building step.
#'
#' For all interaction-selection criteria, \code{print.mfpi()} ends with:
#'
#' \preformatted{
#' Step 3 - Interaction Summary
#' }
#'
#' so the regression output is printed as Step 4.
#'
#' @param criterion Character scalar or \code{NULL}. The selection criterion.
#'
#' @return Integer scalar.
#'
#' @keywords internal
#' @noRd
mfpi_summary_regression_step_number <- function(criterion) {
  4L
}


# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarize an \code{"mfpi"} Object
#'
#' Produces a structured summary of an \code{"mfpi"} fit.
#'
#' The summary object contains the model-building information printed by
#' \code{print.mfpi()}, together with ordinary regression summaries for the
#' retained interaction models. It is designed to support a two-layer display:
#'
#' \enumerate{
#'   \item the MFPI model-building output, including the adjustment model,
#'     candidate interaction models, optional p-value adjustment, and final
#'     interaction-selection summary;
#'   \item regression summaries for the retained interaction models.
#' }
#'
#' The regression summaries are useful for inspecting the final fitted models,
#' but their coefficient-level Wald p-values do not account for fractional
#' polynomial power selection. The MFPI interaction decision should be based on
#' the interaction summary table and the relevant model-selection metric rather
#' than on the regression table alone.
#'
#' @param object An object of class \code{"mfpi"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.mfpi"}. The object is a list
#'   containing the main printed MFPI components, group-level metadata, retained
#'   model objects, and ordinary regression summaries for retained interaction
#'   models.
#'
#' @method summary mfpi
#' @export
summary.mfpi <- function(object, ...) {
  dots <- list(...)
  
  if (length(dots) > 0L) {
    warning(
      "Unused arguments in `summary.mfpi(...)`: ",
      paste(names(dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  if (!inherits(object, "mfpi")) {
    stop("`object` must be an object of class \"mfpi\".", call. = FALSE)
  }
  
  # Build ordinary regression summaries for retained interaction models. These
  # are appended by print.summary.mfpi(); they are not used for interaction
  # selection.
  model_summaries <- build_mfpi_model_summaries(object)
  
  structure(
    list(
      # Basic model metadata ---------------------------------------------------
      call                    = object$call,
      group_var               = object$group_var,
      nobs                    = object$nobs,
      family                  = object$family,
      flex                    = object$flex,
      criterion               = object$criterion,
      p_adjust_method         = object$p_adjust_method,
      p_interact              = object$p_interact,
      min_improvement         = object$min_improvement,
      digits                  = object$digits,
      
      # Group-level metadata --------------------------------------------------
      # group_levels_new are the internal 0, 1, ..., K - 1 codes used by the
      # fitting and prediction internals. group_levels_original are the
      # user-facing labels. group_level_map links the two representations when
      # available.
      group_levels_new        = object$group_levels_new,
      group_levels_original   = object$group_levels_original,
      group_level_map         = object$group_level_map,
      
      # Model-building outputs ------------------------------------------------
      adjust_terms            = object$adjust_terms,
      all_model_metrics       = object$all_model_metrics,
      best_model_metrics      = object$best_model_metrics,
      var_winners             = object$var_winners,
      
      # Fitted interaction model objects --------------------------------------
      best_interaction_model  = object$best_interaction_model,
      all_interaction_models  = object$all_interaction_models,
      
      # Ordinary regression summaries for retained models ---------------------
      model_summaries         = model_summaries
    ),
    class = "summary.mfpi"
  )
}


# -----------------------------------------------------------------------------
# print.summary.mfpi() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print a \code{"summary.mfpi"} Object
#'
#' Displays the structured MFPI summary.
#'
#' The printed summary first reuses \code{print.mfpi()} to display the
#' model-building output. This includes the adjustment model, all interaction
#' candidates, the candidate interaction FP powers by group level, optional
#' p-value adjustment, and the final interaction summary.
#'
#' After the model-building output, the summary print method appends:
#'
#' \enumerate{
#'   \item an unnumbered retained-only detail table showing the group-specific
#'     interaction-side FP powers for retained interaction models;
#'   \item ordinary regression summaries for retained interaction models.
#' }
#'
#' The retained-power table is intentionally not printed as a numbered step. It
#' is a detail table that explains the retained models before the regression
#' output. This avoids confusing labels such as \code{"Step 2b"} when there is
#' no visible \code{"Step 2a"}.
#'
#' The regression output is printed after \code{print.mfpi()}'s Step 3 and is
#' labelled Step 4 for all interaction-selection criteria.
#'
#' @param x An object of class \code{"summary.mfpi"}.
#' @param digits Optional non-negative integer controlling the number of decimal
#'   places used when printing numeric output. If \code{NULL}, \code{x$digits}
#'   is used when available; otherwise the print helper default is used.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @method print summary.mfpi
#' @export
print.summary.mfpi <- function(x, digits = NULL, ...) {
  dots <- list(...)
  
  if (length(dots) > 0L) {
    warning(
      "Unused arguments in `print.summary.mfpi(...)`: ",
      paste(names(dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  if (!inherits(x, "summary.mfpi")) {
    stop("`x` must be an object of class \"summary.mfpi\".", call. = FALSE)
  }
  
  # Reuse print.mfpi() for the model-building output. This keeps print(fit) and
  # print(summary(fit)) consistent. The temporary class switch is safe because
  # summary.mfpi stores the fields consumed by print.mfpi().
  obj <- x
  class(obj) <- "mfpi"
  print.mfpi(obj, digits = digits)
  
  # Use the same display precision as print.mfpi() for any summary-only numeric
  # tables printed below.
  digits <- resolve_print_digits(x, digits)
  
  ruler <- strrep("-", 65)
  
  # Print retained-only interaction FP powers as an unnumbered detail table.
  # This complements the all-candidate power table printed by print.mfpi() and
  # places the final selected group-specific FP powers next to the regression
  # summaries.
  print_interaction_power_details_step(
    metrics   = x$best_model_metrics,
    group_var = x$group_var,
    ruler     = ruler,
    title     = "  Retained FP powers by group level"
  )
  
  # The retained-power table above is unnumbered. Since print.mfpi() ends with
  # Step 3 - Interaction Summary for all criteria, the regression-output block is
  # printed as Step 4.
  regression_step_no <- mfpi_summary_regression_step_number(x$criterion)
  
  cat(sprintf(
    "\nStep %d - Regression Output for Retained Interaction Models:\n",
    regression_step_no
  ))
  cat(ruler, "\n")
  
  ms <- x$model_summaries
  
  # No retained models: report explicitly. This can occur when no continuous
  # variable meets the interaction-selection criterion.
  if (!has_retained_model_summaries(ms)) {
    cat("  No retained interaction models.\n")
    return(invisible(x))
  }
  
  # Print one regression summary per retained interaction model. The variable
  # name comes from the model_summaries list name. The interaction type is taken
  # from var_winners when available so that users can see whether the retained
  # model was linear, FP1, or FP2.
  for (vn in names(ms)) {
    sm <- ms[[vn]]
    
    if (is.null(sm)) {
      next
    }
    
    w <- if (!is.null(x$var_winners)) x$var_winners[[vn]] else NULL
    
    type_label <- if (!is.null(w) && !is.null(w$type)) {
      format_type_label(w$type)
    } else {
      "?"
    }
    
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