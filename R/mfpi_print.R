# S3 print method for mfpi objects
#
# print.mfpi() gives a structured three-section output:
#   Step 1 - Adjustment model (all fp_terms columns for selected variables)
#   Step 2 - All interaction candidates (linear/FP1/FP2) with winner flagged
#   Step 3 - Selected interactions (variables that cleared the threshold)


# -----------------------------------------------------------------------------
# print.mfpi() ----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print an \code{"mfpi"} Object
#'
#' Displays a structured three-step summary of an \code{"mfpi"} fit:
#' \enumerate{
#'   \item The adjustment model selected by MFP, showing all \code{fp_terms}
#'     columns for retained variables and a compact list of dropped variables.
#'   \item All interaction candidates (linear, FP1, FP2) for every variable in
#'     \code{cont_vars}, with the selected model flagged by \code{*}. This
#'     allows the user to judge whether the evidence is consistent across
#'     functional forms or driven by a single degree.
#'   \item The final selected interactions - one row per variable that cleared
#'     the selection threshold - or a message if none were found.
#' }
#'
#' Regression output for the winning models is available via
#' \code{summary.mfpi()}.
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
  
  cat(header, "\n")
  cat(sprintf(
    "MFPI  |  group: '%s'  |  n = %d  |  %s\n",
    x$group_var, x$nobs, x$flex
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
  # Step 2: All interaction candidates with winner flagged
  # ---------------------------------------------------------------------------
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")
  
  all_m <- x$all_model_metrics
  
  if (!is.null(all_m) && nrow(all_m) > 0L) {
    
    # Identify winners: match on type + variable
    best_m <- x$best_model_metrics
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      winner_keys <- paste(best_m$type, best_m$variable, sep = "__")
    } else {
      winner_keys <- character(0L)
    }
    
    # Format power list-columns as strings
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
    
    # Add winner flag as last column; type is already first from assembly
    row_keys     <- paste(all_m$type, all_m$variable, sep = "__")
    all_m$winner <- ifelse(row_keys %in% winner_keys, "*", "")
    
    print(as.data.frame(all_m), row.names = FALSE)
    cat("\n  * = selected as best functional form for that variable\n")
    
  } else {
    cat("  No interaction candidates computed.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: Selected interactions
  # ---------------------------------------------------------------------------
  cat("\nStep 3 - Selected Interactions:\n")
  cat(ruler, "\n")
  
  best_m <- x$best_model_metrics
  
  if (!is.null(best_m) && nrow(best_m) > 0L) {
    
    # Format power list-columns
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
    
    if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
      cat(sprintf(
        "\n  For flex4, 'fp_powers_int' shows per-group FP powers (levels: %s).\n",
        paste(x$group_levels_original, collapse = ", ")
      ))
    }
    cat(
      paste0(
        "\n  df_interaction: degrees of freedom for the likelihood-ratio test\n",
        "    (extra parameters in interaction model vs main-effects model;\n",
        "    used for the p-value or information-criterion comparison).\n",
        "  df_total: total parameters in the interaction model excluding\n",
        "    the intercept and adjustment variables; includes group dummies\n",
        "    and all group-specific FP terms.\n"
      )
    )
    
  } else {
    cat("  No significant interactions found.\n")
  }
  
  invisible(x)
}
# S3 summary method for mfpi objects
#
# summary.mfpi()       - collects all results into a "summary.mfpi" list
# print.summary.mfpi() - four-step formatted display:
#   Step 1 - Adjustment model
#   Step 2 - All interaction candidates (with winner flagged)
#   Step 3 - Selected interactions
#   Step 4 - Regression output for winning models (coefficients + vcov)


# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarise an \code{"mfpi"} Object
#'
#' Returns a structured summary of an \code{"mfpi"} fit. The returned object
#' has a dedicated print method that displays four sections:
#' \enumerate{
#'   \item Adjustment model selected by MFP.
#'   \item All interaction candidates (linear, FP1, FP2) with the winner
#'     flagged.
#'   \item Final selected interactions.
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
#'   \item{\code{group_levels_original}}{Original levels of \code{group_var}.}
#'   \item{\code{adjust_terms}}{Full \code{fp_terms} table from the MFP
#'     adjustment model.}
#'   \item{\code{all_model_metrics}}{Metrics for all three candidates (linear,
#'     FP1, FP2) for every variable in \code{cont_vars}.}
#'   \item{\code{best_model_metrics}}{Metrics for the selected (winning) model
#'     per significant variable.}
#'   \item{\code{best_interaction_model}}{Named list of fitted interaction model
#'     objects, one per significant variable.}
#'   \item{\code{all_interaction_models}}{Named list of all candidate interaction
#'     model objects, one element per variable.}
#'   \item{\code{best_fitted_functions}}{Named list of fitted-function matrices
#'     for the winning model per significant variable. See \code{mfpi()} for
#'     column definitions.}
#'   \item{\code{all_fitted_functions}}{Named list of fitted-function matrices
#'     for all candidates, one element per variable.}
#'   \item{\code{model_summaries}}{Named list of \code{summary(fit$fit)}
#'     objects, one per significant variable. \strong{Note:} p-values here
#'     are from \code{glm}/\code{coxph} and do not account for FP power
#'     estimation df. Use \code{best_model_metrics$pvalue} for the
#'     correct interaction test p-values.}
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
      group_var                = object$group_var,
      nobs                     = object$nobs,
      family                   = object$family,
      flex                     = object$flex,
      group_levels_original    = object$group_levels_original,
      adjust_terms             = object$adjust_terms,
      all_model_metrics        = object$all_model_metrics,
      best_model_metrics = object$best_model_metrics,
      best_interaction_model       = object$best_interaction_model,
      all_interaction_models   = object$all_interaction_models,
      best_fitted_functions         = object$best_fitted_functions,
      all_fitted_functions     = object$all_fitted_functions,
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
#' Displays a four-step structured summary. Steps 1-3 match \code{print.mfpi()};
#' Step 4 adds regression output for each final interaction model.
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
  
  cat(header, "\n")
  cat(sprintf(
    "MFPI Summary  |  group: '%s'  |  n = %d  |  %s\n",
    x$group_var, x$nobs, x$flex
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
  # Step 2: All interaction candidates with winner flagged
  # ---------------------------------------------------------------------------
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")
  
  all_m <- x$all_model_metrics
  
  if (!is.null(all_m) && nrow(all_m) > 0L) {
    
    # Identify winners: match on type + variable
    best_m <- x$best_model_metrics
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      winner_keys <- paste(best_m$type, best_m$variable, sep = "__")
    } else {
      winner_keys <- character(0L)
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
    
    # Add winner flag as last column; type is already first from assembly
    row_keys       <- paste(all_m$type, all_m$variable, sep = "__")
    all_m$winner   <- ifelse(row_keys %in% winner_keys, "*", "")
    
    print(as.data.frame(all_m), row.names = FALSE)
    cat("\n  * = selected as best functional form for that variable\n")
    
  } else {
    cat("  No interaction candidates computed.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: Selected interactions
  # ---------------------------------------------------------------------------
  cat("\nStep 3 - Selected Interactions:\n")
  cat(ruler, "\n")
  
  best_m <- x$best_model_metrics
  
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
    
    if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
      cat(sprintf(
        "\n  For flex4, 'fp_powers_int' shows per-group FP powers (levels: %s).\n",
        paste(x$group_levels_original, collapse = ", ")
      ))
    }
    cat(
      paste0(
        "\n  df_interaction: degrees of freedom for the likelihood-ratio test\n",
        "    (extra parameters in interaction model vs main-effects model;\n",
        "    used for the p-value or information-criterion comparison).\n",
        "  df_total: total parameters in the interaction model excluding\n",
        "    the intercept and adjustment variables; includes group dummies\n",
        "    and all group-specific FP terms.\n"
      )
    )
    cat("        p-values correctly account for interaction test degrees of freedom.\n")
    
  } else {
    cat("  No significant interactions found.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Regression output for winning models
  # ---------------------------------------------------------------------------
  if (length(x$model_summaries) > 0L) {
    
    cat("\nStep 4 - Regression Output (Winning Interaction Models):\n")
    cat(ruler, "\n")
    cat("  Note: p-values below are from glm/coxph and do not account for\n")
    cat("        FP power estimation df. Use Step 3 pvalues for the correct\n")
    cat("        interaction test p-values.\n")
    
    best_m <- x$best_model_metrics
    type_lookup <- if (!is.null(best_m) && "type" %in% names(best_m)) {
      setNames(toupper(gsub("fp", "FP", best_m$type)), best_m$variable)
    } else {
      setNames(rep("", length(x$model_summaries)), names(x$model_summaries))
    }
    
    for (var_name in names(x$model_summaries)) {
      sm       <- x$model_summaries[[var_name]]
      form_lbl <- if (!is.null(type_lookup[[var_name]])) type_lookup[[var_name]] else ""
      cat(sprintf("\n  '%s' x '%s'  [%s]\n", x$group_var, var_name, form_lbl))
      cat(strrep("-", 45), "\n")
      if (!is.null(sm)) print(sm) else cat("  (model summary unavailable)\n")
    }
  }
  
  invisible(x)
}

# -----------------------------------------------------------------------------
# print.summary.mfpi() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print a \code{"summary.mfpi"} Object
#'
#' Displays a four-step structured summary. Steps 1-3 match \code{print.mfpi()};
#' Step 4 adds regression output for each final interaction model.
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
  
  cat(header, "\n")
  cat(sprintf(
    "MFPI Summary  |  group: '%s'  |  n = %d  |  %s\n",
    x$group_var, x$nobs, x$flex
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
  # Step 2: All interaction candidates with winner flagged
  # ---------------------------------------------------------------------------
  cat("\nStep 2 - All Interaction Candidates:\n")
  cat(ruler, "\n")
  
  all_m <- x$all_model_metrics
  
  if (!is.null(all_m) && nrow(all_m) > 0L) {
    
    # Identify winners: match on type + variable
    best_m <- x$best_model_metrics
    if (!is.null(best_m) && nrow(best_m) > 0L) {
      winner_keys <- paste(best_m$type, best_m$variable, sep = "__")
    } else {
      winner_keys <- character(0L)
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
    
    # Add winner flag as last column; type is already first from assembly
    row_keys       <- paste(all_m$type, all_m$variable, sep = "__")
    all_m$winner   <- ifelse(row_keys %in% winner_keys, "*", "")
    
    print(as.data.frame(all_m), row.names = FALSE)
    cat("\n  * = selected as best functional form for that variable\n")
    
  } else {
    cat("  No interaction candidates computed.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: Selected interactions
  # ---------------------------------------------------------------------------
  cat("\nStep 3 - Selected Interactions:\n")
  cat(ruler, "\n")
  
  best_m <- x$best_model_metrics
  
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
    
    if (identical(x$flex, "flex4") && !is.null(x$group_levels_original)) {
      cat(sprintf(
        "\n  For flex4, 'fp_powers_int' shows per-group FP powers (levels: %s).\n",
        paste(x$group_levels_original, collapse = ", ")
      ))
    }
    cat("\n  Note: df_interaction = LRT df; df_total = total df in interaction model.\n")
    cat("        p-values correctly account for interaction test degrees of freedom.\n")
    
  } else {
    cat("  No significant interactions found.\n")
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Regression output for winning models
  # ---------------------------------------------------------------------------
  if (length(x$model_summaries) > 0L) {
    
    cat("\nStep 4 - Regression Output (Winning Interaction Models):\n")
    cat(ruler, "\n")
    cat("  Note: p-values below are from glm/coxph and do not account for\n")
    cat("        FP power estimation df. Use Step 3 pvalues for the correct\n")
    cat("        interaction test p-values.\n")
    
    best_m <- x$best_model_metrics
    type_lookup <- if (!is.null(best_m) && "type" %in% names(best_m)) {
      setNames(toupper(gsub("fp", "FP", best_m$type)), best_m$variable)
    } else {
      setNames(rep("", length(x$model_summaries)), names(x$model_summaries))
    }
    
    for (var_name in names(x$model_summaries)) {
      sm       <- x$model_summaries[[var_name]]
      form_lbl <- if (!is.null(type_lookup[[var_name]])) type_lookup[[var_name]] else ""
      cat(sprintf("\n  '%s' x '%s'  [%s]\n", x$group_var, var_name, form_lbl))
      cat(strrep("-", 45), "\n")
      if (!is.null(sm)) print(sm) else cat("  (model summary unavailable)\n")
    }
  }
  
  invisible(x)
}