# S3 summary method for mfpi objects
#
# summary.mfpi() creates a structured summary object for an MFPI fit.
# print.summary.mfpi() reuses print.mfpi() for Steps 1-3 and appends a
# user-facing regression display for retained interaction models in Step 4.
#
# Step 4 is deliberately MFPI-specific rather than a verbatim printout of
# summary.glm(), summary.coxph(), or another underlying model summary. This
# avoids exposing internal design-column names and separates group-specific FP
# terms, group-variable coefficients, and adjustment coefficients.
#
# Coefficient-level tests in Step 4 are conditional on the selected model and
# FP transformations. They are descriptive summaries of the retained fitted
# regression models and are not used for the MFPI interaction decision.

# -----------------------------------------------------------------------------
# Internal summary helpers -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Build Regression Displays for Retained MFPI Interaction Models
#'
#' Creates the structured coefficient information used by
#' \code{print.summary.mfpi()} for each retained interaction model.
#'
#' @param object Object of class \code{"mfpi"}.
#'
#' @return A named list. Each non-\code{NULL} element contains the readable
#'   coefficient metadata produced by \code{mfpi_coefficient_info()} and the
#'   corresponding coefficient-level standard errors, test statistics, and
#'   p-values from the retained fitted model.
#'
#' @keywords internal
#' @noRd
build_mfpi_regression_displays <- function(object) {
  retained <- object$best_interaction_model
  if (is.null(retained) || length(retained) == 0L) {
    return(list())
  }

  terms <- names(retained)
  if (is.null(terms)) {
    return(list())
  }
  keep <- !vapply(retained, is.null, logical(1L)) & !is.na(terms) & nzchar(terms)
  terms <- terms[keep]
  if (length(terms) == 0L) {
    return(list())
  }

  out <- stats::setNames(vector("list", length(terms)), terms)

  for (term in terms) {
    winner <- if (!is.null(object$var_winners)) object$var_winners[[term]] else NULL
    fit_result <- if (!is.null(winner)) winner$fit else NULL

    if (is.null(fit_result) || is.null(fit_result$test_results) ||
        is.null(fit_result$test_results$interaction_model) ||
        is.null(fit_result$test_results$interaction_model$fit)) {
      out[[term]] <- NULL
      next
    }

    info <- mfpi_coefficient_info(object, term, fit_result)
    fit_obj <- fit_result$test_results$interaction_model$fit
    statistics <- mfpi_summary_coefficient_statistics(
      fit_obj = fit_obj,
      raw_names = info$raw_names
    )
    out[[term]] <- list(
      info = info,
      statistics = statistics,
      family_string = object$family_string,
      model_metadata = mfp2_summary_model_metadata(
        fit_obj,
        family_string = object$family_string
      )
    )
  }

  out
}


#' Extract Coefficient-Level Statistics for an MFPI Regression Display
#'
#' Uses the fitted coefficient vector and variance-covariance matrix for the
#' estimate and standard error, and uses the underlying model summary when
#' available to recover the model-specific test statistic and p-value.
#'
#' @param fit_obj Underlying fitted regression model.
#' @param raw_names Character vector of fitted coefficient names in the order
#'   required by the MFPI display.
#'
#' @return A data frame with one row per coefficient and columns
#'   \code{raw_name}, \code{estimate}, \code{std_error}, \code{statistic}, and
#'   \code{p_value}. The model-specific statistic label (usually \code{"z"} or
#'   \code{"t"}) is stored in the \code{"statistic_label"} attribute.
#'
#' @keywords internal
#' @noRd
mfpi_summary_coefficient_statistics <- function(fit_obj, raw_names) {
  beta <- stats::coef(fit_obj)
  V <- stats::vcov(fit_obj)

  # Multinomial fits expose a per-logit coefficient matrix. Flatten it into the
  # same class-major, single-colon-named vector used by vcov.multinom() so the
  # shared alignment below works unchanged.
  if (is.matrix(beta)) {
    flat_names <- as.vector(t(outer(
      rownames(beta), colnames(beta), function(a, b) paste0(a, ":", b)
    )))
    beta <- stats::setNames(as.vector(t(beta)), flat_names)
  }

  if (!is.numeric(beta) || is.null(names(beta))) {
    stop("The retained model does not expose named coefficients.", call. = FALSE)
  }
  if (!is.matrix(V) || is.null(rownames(V)) || is.null(colnames(V))) {
    stop("The retained model does not expose a named covariance matrix.", call. = FALSE)
  }

  missing_beta <- setdiff(raw_names, names(beta))
  missing_vcov <- setdiff(raw_names, intersect(rownames(V), colnames(V)))
  if (length(missing_beta) > 0L || length(missing_vcov) > 0L) {
    stop("Retained-model coefficients could not be aligned for summary output.", call. = FALSE)
  }

  estimate <- unname(beta[raw_names])
  V <- V[raw_names, raw_names, drop = FALSE]
  std_error <- sqrt(pmax(0, diag(V)))
  statistic <- estimate / std_error
  p_value <- rep(NA_real_, length(raw_names))
  statistic_label <- "Statistic"

  sm <- tryCatch(summary(fit_obj), error = function(e) NULL)
  coef_table <- if (!is.null(sm) && !is.null(sm$coefficients)) {
    sm$coefficients
  } else if (!is.null(sm)) {
    # summary.survreg() calls its coefficient matrix `table` and appends scale
    # rows; the raw-name matching below selects only regression coefficients.
    sm$table
  } else {
    NULL
  }

  if (is.data.frame(coef_table)) {
    coef_table <- as.matrix(coef_table)
  }

  if (is.matrix(coef_table) && !is.null(rownames(coef_table))) {
    rows <- match(raw_names, rownames(coef_table))
    if (all(!is.na(rows))) {
      cn <- colnames(coef_table)
      if (is.null(cn)) cn <- rep("", ncol(coef_table))

      stat_idx <- mfpi_summary_statistic_column(cn)
      p_idx <- mfpi_summary_pvalue_column(cn)

      if (!is.na(stat_idx)) {
        candidate <- suppressWarnings(as.numeric(coef_table[rows, stat_idx]))
        usable <- is.finite(candidate) | is.na(candidate)
        if (all(usable)) statistic <- candidate
        statistic_label <- mfpi_summary_statistic_label(cn[stat_idx])
      }

      if (!is.na(p_idx)) {
        candidate <- suppressWarnings(as.numeric(coef_table[rows, p_idx]))
        if (length(candidate) == length(p_value)) p_value <- candidate
      }
    }
  }

  # If an underlying summary does not expose p-values, provide the standard
  # model-conditional reference calculation for the common fitted-model classes.
  if (all(is.na(p_value))) {
    if (identical(statistic_label, "t") && !is.null(fit_obj$df.residual) &&
        length(fit_obj$df.residual) == 1L && is.finite(fit_obj$df.residual)) {
      p_value <- 2 * stats::pt(
        abs(statistic),
        df = fit_obj$df.residual,
        lower.tail = FALSE
      )
    } else if (inherits(fit_obj, "glm") || inherits(fit_obj, "coxph") ||
               inherits(fit_obj, "survreg") ||
               inherits(fit_obj, "fastglm") ||
               inherits(fit_obj, "multinom")) {
      statistic_label <- if (identical(statistic_label, "Statistic")) "z" else statistic_label
      p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    }
  }

  out <- data.frame(
    raw_name = raw_names,
    estimate = estimate,
    std_error = unname(std_error),
    statistic = unname(statistic),
    p_value = unname(p_value),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  attr(out, "statistic_label") <- statistic_label
  out
}


#' Locate the Test-Statistic Column in a Coefficient Table
#'
#' Returns the position of the test-statistic column in an underlying
#' model summary coefficient table, matching a small set of common
#' spellings after normalising case and whitespace. Used when the MFPI
#' summary has to extract a Wald-style z or t statistic from tables
#' produced by different fitting backends.
#'
#' @param column_names Character vector of coefficient-table column names.
#'
#' @return Integer scalar column position, or `NA_integer_` when no
#'   candidate label matches.
#'
#' @keywords internal
#' @noRd
mfpi_summary_statistic_column <- function(column_names) {
  if (length(column_names) == 0L) return(NA_integer_)
  normalized <- tolower(trimws(column_names))
  candidates <- c("z", "z value", "t", "t value", "statistic")
  idx <- match(candidates, normalized, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx) == 0L) NA_integer_ else idx[1L]
}


#' Locate the P-Value Column in a Coefficient Table
#'
#' Returns the position of the p-value column in an underlying model
#' summary coefficient table, matching a small set of common spellings
#' after normalising case and whitespace.
#'
#' @param column_names Character vector of coefficient-table column names.
#'
#' @return Integer scalar column position, or `NA_integer_` when no
#'   candidate label matches.
#'
#' @keywords internal
#' @noRd
mfpi_summary_pvalue_column <- function(column_names) {
  if (length(column_names) == 0L) return(NA_integer_)
  idx <- grep(
    "^Pr\\(|^p$|p[-._ ]?value|pvalue",
    column_names,
    ignore.case = TRUE
  )
  if (length(idx) == 0L) NA_integer_ else idx[1L]
}


#' Concise Test-Statistic Column Label
#'
#' Converts an underlying coefficient-table test-statistic heading into the
#' concise label used in the printed MFPI summary (`"z"`, `"t"`, or
#' `"Statistic"`).
#'
#' @param column_name Character scalar column name from an underlying
#'   coefficient table.
#'
#' @return Character scalar with the concise label.
#'
#' @keywords internal
#' @noRd
mfpi_summary_statistic_label <- function(column_name) {
  x <- tolower(trimws(column_name))
  if (grepl("^z($| )", x)) return("z")
  if (grepl("^t($| )", x)) return("t")
  "Statistic"
}


#' Are Any Retained Regression Displays Available?
#'
#' Reports whether at least one retained regression-display block was
#' assembled for the MFPI summary. Callers use this to decide whether to
#' print the regression-block section header at all.
#'
#' @param regression_displays List of regression-display blocks, possibly
#'   empty.
#'
#' @return `TRUE` when at least one non-empty display is present, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
has_retained_regression_displays <- function(regression_displays) {
  !is.null(regression_displays) &&
    length(regression_displays) > 0L &&
    any(!vapply(regression_displays, is.null, logical(1L)))
}


#' Resolve the Step Number for Retained Regression Output
#'
#' Returns the numbered MFPI step under which retained regression output
#' is printed. The step number depends on the active selection criterion,
#' because closed-testing p-value selection and information-criterion
#' selection produce different numbers of preceding steps.
#'
#' @param criterion Character scalar selection criterion.
#'
#' @return Integer scalar step number.
#'
#' @keywords internal
#' @noRd
mfpi_summary_regression_step_number <- function(criterion) {
  4L
}


#' Does One Retained Interaction Use a Nonzero MFPI Shift?
#'
#' Reports whether one retained interaction display uses a nonzero
#' covariate-wide MFPI shift. When it does, the shift-explanation note is
#' printed inside that interaction's block, immediately after the table
#' whose transformation labels contain the shift, rather than as a
#' summary-wide footer.
#'
#' @param display One retained interaction display.
#'
#' @return `TRUE` when the display uses a nonzero shift, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfpi_summary_display_has_nonzero_shift <- function(display) {
  if (is.null(display) || is.null(display$info$shift)) return(FALSE)

  value <- display$info$shift
  length(value) == 1L && !is.na(value) && is.finite(value) && value != 0
}


#' Print the Shift-Explanation Note for One Interaction Block
#'
#' Prints the explanatory note describing the covariate-wide MFPI shift
#' used by one retained interaction display. The note is placed next to
#' the group-specific FP table for that interaction, so it never suggests
#' that the shift is a global feature of every retained interaction.
#'
#' @param display One retained interaction display.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
mfpi_print_summary_shift_note <- function(display) {
  if (!mfpi_summary_display_has_nonzero_shift(display)) {
    return(invisible(FALSE))
  }

  cat("\nNote: The constant added inside the FP transformation is the shifting factor\n")
  cat("      applied to all observations before constructing the group-specific FP\n")
  cat("      terms; structural zeros therefore remain untransformed.\n")
  invisible(TRUE)
}


#' Build One User-Facing Coefficient Table for `print.summary.mfpi()`
#'
#' Assembles one user-facing coefficient statistics table for the MFPI
#' summary layout, combining the term metadata (variable, basis,
#' centering) with the fitted regression statistics into the printed
#' column order.
#'
#' @param metadata_table Data frame of per-coefficient metadata.
#' @param statistics Data frame of fitted regression statistics (estimate,
#'   standard error, test statistic, p value).
#' @param leading_columns Character vector of metadata columns to place at
#'   the front of the display table, in order.
#'
#' @return A data frame with columns arranged for direct printing.
#'
#' @keywords internal
#' @noRd
mfpi_summary_display_table <- function(metadata_table, statistics, leading_columns) {
  if (is.null(metadata_table) || nrow(metadata_table) == 0L) {
    return(NULL)
  }

  idx <- match(metadata_table$raw_name, statistics$raw_name)
  if (anyNA(idx)) {
    stop("Summary coefficient metadata do not align with fitted coefficients.", call. = FALSE)
  }

  out <- metadata_table[, leading_columns, drop = FALSE]
  out$Estimate <- statistics$estimate[idx]
  out[["S.E."]] <- statistics$std_error[idx]
  stat_label <- attr(statistics, "statistic_label", exact = TRUE)
  if (is.null(stat_label) || !nzchar(stat_label)) stat_label <- "Statistic"
  out[[stat_label]] <- statistics$statistic[idx]
  out[["p-value"]] <- statistics$p_value[idx]
  out
}


#' Format Stored Centering Constants at the Print Boundary
#'
#' Formats a vector of stored centering constants for the printed MFPI
#' summary. Missing values become blank cells because no centering
#' constant applies to that row; all actual constants use the same
#' fixed-decimal convention as the coefficient statistics so the columns
#' line up.
#'
#' @param x Numeric vector of centering constants (possibly with `NA`).
#' @param digits Integer scalar; number of decimal places.
#'
#' @return Character vector of formatted constants.
#'
#' @keywords internal
#' @noRd
mfpi_summary_format_centers <- function(x, digits) {
  format_print_decimal(x, digits, na_string = "")
}


#' Does an Interaction Display Contain Any Centering Constants?
#'
#' Reports whether a retained interaction display carries any stored
#' centering constants in either its group-specific or adjustment basis
#' metadata. Callers use this to decide whether the `Center` column should
#' appear in the printed table at all.
#'
#' @param info Interaction display information list.
#'
#' @return `TRUE` when at least one stored centering constant is present,
#'   `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfpi_summary_info_has_centers <- function(info) {
  tables <- list(info$group_table, info$adjustment_table)
  any(vapply(tables, function(tab) {
    !is.null(tab) && "center" %in% names(tab) && any(!is.na(tab$center))
  }, logical(1L)))
}


#' Format P-Values for a Regression-Display Table
#'
#' Formats a numeric p-value vector for one MFPI regression-display table
#' without appending significance stars. Significance stars are avoided in
#' MFPI output because model selection has already been performed and the
#' printed p values are conditional on the retained model.
#'
#' @param p Numeric vector of p values.
#' @param digits Integer scalar; number of decimal places.
#'
#' @return Character vector of formatted p values.
#'
#' @keywords internal
#' @noRd
mfpi_summary_format_pvalues <- function(p, digits) {
  format_print_pvalue(p, digits)
}


#' Print One Coefficient Table in the MFPI Summary Layout
#'
#' Prints one coefficient statistics table using the shared MFPI summary
#' layout: fixed-decimal numeric columns, left-aligned metadata columns,
#' and a compact header rule.
#'
#' @param tab Data frame of coefficient statistics with layout produced by
#'   `mfpi_summary_display_table()`.
#' @param digits Integer scalar; number of decimal places.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
mfpi_print_summary_coefficient_table <- function(tab, digits) {
  if (is.null(tab) || nrow(tab) == 0L) return(invisible(NULL))
  if ("Center" %in% names(tab)) {
    tab[["Center"]] <- mfpi_summary_format_centers(tab[["Center"]], digits)
  }
  tab[["p-value"]] <- mfpi_summary_format_pvalues(tab[["p-value"]], digits)
  tab <- format_model_print_table(tab, digits)
  print(tab, row.names = FALSE)
  invisible(NULL)
}


#' Print One Retained MFPI Interaction Block
#'
#' Prints one retained interaction model in the MFPI-specific summary
#' layout: interaction identifier, selected FP powers per group,
#' group-specific and adjustment coefficient tables, and any per-block
#' shift-explanation note.
#'
#' @param display One retained interaction display.
#' @param digits Integer scalar; number of decimal places.
#' @param print_shift_note Logical. If `TRUE`, print the shift-explanation
#'   note for this interaction when a nonzero shift was used.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
mfpi_print_summary_regression_block <- function(display, digits, print_shift_note = TRUE) {
  info <- display$info
  statistics <- display$statistics

  if (isTRUE(info$multinomial)) {
    return(mfpi_print_summary_multinomial_block(display, digits))
  }

  heading <- paste0("Interaction: ", info$term)
  cat(heading, "\n", sep = "")
  cat(strrep("-", nchar(heading)), "\n", sep = "")
  if (identical(info$form, "Linear")) {
    cat("Form: Linear\n")
  } else {
    cat("FP form: ", info$form, "\n", sep = "")
  }

  family_string <- display$family_string
  metadata <- display$model_metadata
  if (identical(family_string, "survreg") &&
      !is.null(metadata$distribution)) {
    cat("Distribution: ", metadata$distribution, "\n", sep = "")
  }
  if (identical(family_string, "negbin") && !is.null(metadata$link)) {
    cat("Link: ", metadata$link, "\n", sep = "")
  }
  mfp2_print_model_specific_parameters(
    family_string = family_string,
    metadata = metadata,
    digits = digits
  )

  mfpi_print_interaction_powers(info$power_table)

  if (mfpi_summary_info_has_centers(info)) {
    cat("\nEstimates are for: Basis - Center\n")
  }

  group_leading <- c("group_level", "transformation")
  if ("center" %in% names(info$group_table) &&
      any(!is.na(info$group_table$center))) {
    group_leading <- c(group_leading, "center")
  }
  group_tab <- mfpi_summary_display_table(
    metadata_table = info$group_table,
    statistics = statistics,
    leading_columns = group_leading
  )
  if (!is.null(group_tab)) {
    names(group_tab)[1:2] <- c("Group level", "Basis")
    if ("center" %in% names(group_tab)) {
      names(group_tab)[names(group_tab) == "center"] <- "Center"
    }
    cat("\nGroup-specific FP terms:\n")
    mfpi_print_summary_coefficient_table(group_tab, digits)

    # A nonzero covariate-wide shift is visible inside the FP labels above.
    # The caller controls whether the explanatory note is emitted so a summary
    # with several shifted retained interactions can explain the convention once
    # at its first occurrence without repeating the same note in later blocks.
    if (isTRUE(print_shift_note)) {
      mfpi_print_summary_shift_note(display)
    }
  }

  group_var_tab <- mfpi_summary_display_table(
    metadata_table = info$group_variable_table,
    statistics = statistics,
    leading_columns = "level"
  )
  if (!is.null(group_var_tab)) {
    names(group_var_tab)[1L] <- "Level"
    cat("\nGroup-variable coefficients:\n")
    cat("Reference level: ", info$reference_level, "\n\n", sep = "")
    mfpi_print_summary_coefficient_table(group_var_tab, digits)
  }

  intercept_tab <- mfpi_summary_display_table(
    metadata_table = info$intercept_table,
    statistics = statistics,
    leading_columns = "term"
  )
  if (!is.null(intercept_tab)) {
    names(intercept_tab)[1L] <- "Term"
    cat("\nModel intercept:\n")
    mfpi_print_summary_coefficient_table(intercept_tab, digits)
  }

  adjustment_leading <- if (all(c("variable", "basis") %in%
                                names(info$adjustment_table))) {
    c("variable", "basis")
  } else {
    "term"
  }
  if ("center" %in% names(info$adjustment_table) &&
      any(!is.na(info$adjustment_table$center))) {
    adjustment_leading <- c(adjustment_leading, "center")
  }
  adjustment_tab <- mfpi_summary_display_table(
    metadata_table = info$adjustment_table,
    statistics = statistics,
    leading_columns = adjustment_leading
  )
  if (!is.null(adjustment_tab)) {
    names(adjustment_tab)[names(adjustment_tab) == "variable"] <- "Variable"
    names(adjustment_tab)[names(adjustment_tab) == "basis"] <- "Basis"
    names(adjustment_tab)[names(adjustment_tab) == "center"] <- "Center"
    names(adjustment_tab)[names(adjustment_tab) == "term"] <- "Term"
    cat("\nAdjustment coefficients:\n")
    mfpi_print_summary_coefficient_table(adjustment_tab, digits)
  }

  invisible(NULL)
}


#' Print One Retained Multinomial MFPI Interaction Block
#'
#' Prints the summary layout for a multinomial MFPI interaction model:
#' one coefficient statistics table per non-reference logit, each against
#' the common reference class named in the header.
#'
#' @param display One retained multinomial interaction display.
#' @param digits Integer scalar; number of decimal places.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
mfpi_print_summary_multinomial_block <- function(display, digits) {
  info <- display$info
  statistics <- display$statistics

  heading <- paste0("Interaction: ", info$term)
  cat(heading, "\n", sep = "")
  cat(strrep("-", nchar(heading)), "\n", sep = "")
  if (identical(info$form, "Linear")) {
    cat("Form: Linear\n")
  } else {
    cat("FP form: ", info$form, "\n", sep = "")
  }
  cat("Reference outcome: ", info$reference_class, "\n", sep = "")

  mfpi_print_interaction_powers(info$power_table)

  if (mfpi_summary_info_has_centers(info)) {
    cat("\nEstimates are for: Basis - Center\n")
  }

  for (cls in info$classes) {
    outcome_heading <- paste0("Outcome: ", cls, " vs ", info$reference_class)
    cat("\n", outcome_heading, "\n", sep = "")
    cat(strrep("-", nchar(outcome_heading)), "\n", sep = "")

    if (!is.null(info$group_table)) {
      group_meta <- info$group_table[
        info$group_table$class == cls, , drop = FALSE
      ]
      group_leading <- c("group_level", "transformation")
      if ("center" %in% names(group_meta) && any(!is.na(group_meta$center))) {
        group_leading <- c(group_leading, "center")
      }
      group_tab <- mfpi_summary_display_table(
        metadata_table = group_meta,
        statistics = statistics,
        leading_columns = group_leading
      )
      if (!is.null(group_tab)) {
        names(group_tab)[1:2] <- c("Group level", "Basis")
        names(group_tab)[names(group_tab) == "center"] <- "Center"
        cat("\nGroup-specific FP terms:\n")
        mfpi_print_summary_coefficient_table(group_tab, digits)
      }

      group_var_meta <- info$group_variable_table[
        info$group_variable_table$class == cls, , drop = FALSE
      ]
      group_var_tab <- mfpi_summary_display_table(
        metadata_table = group_var_meta,
        statistics = statistics,
        leading_columns = "level"
      )
      if (!is.null(group_var_tab)) {
        names(group_var_tab)[1L] <- "Level"
        cat("\nGroup-variable coefficients:\n")
        cat("Reference level: ", info$reference_level, "\n\n", sep = "")
        mfpi_print_summary_coefficient_table(group_var_tab, digits)
      }

      adjustment_meta <- info$adjustment_table[
        info$adjustment_table$class == cls, , drop = FALSE
      ]
      adjustment_leading <- c("variable", "basis")
      if ("center" %in% names(adjustment_meta) &&
          any(!is.na(adjustment_meta$center))) {
        adjustment_leading <- c(adjustment_leading, "center")
      }
      adjustment_tab <- mfpi_summary_display_table(
        metadata_table = adjustment_meta,
        statistics = statistics,
        leading_columns = adjustment_leading
      )
      if (!is.null(adjustment_tab)) {
        names(adjustment_tab)[names(adjustment_tab) == "variable"] <- "Variable"
        names(adjustment_tab)[names(adjustment_tab) == "basis"] <- "Basis"
        names(adjustment_tab)[names(adjustment_tab) == "center"] <- "Center"
        cat("\nAdjustment coefficients:\n")
        mfpi_print_summary_coefficient_table(adjustment_tab, digits)
      }
    } else {
      # Backward-compatible rendering for older summary objects that contain
      # only the original flat per-logit metadata table.
      meta <- info$logit_table[info$logit_table$class == cls, , drop = FALSE]
      tab <- mfpi_summary_display_table(
        metadata_table = meta,
        statistics = statistics,
        leading_columns = "term"
      )
      if (!is.null(tab)) {
        names(tab)[1L] <- "Term"
        mfpi_print_summary_coefficient_table(tab, digits)
      }
    }
  }

  invisible(NULL)
}


#' Print Every Retained MFPI Regression Display
#'
#' Prints every retained regression display in order. The
#' shift-explanation note is shown at most once, immediately after the
#' first retained interaction whose focal variable uses a nonzero MFPI
#' shift; later shifted interactions still show the shift explicitly in
#' their transformation labels, so repeating the note adds no information.
#'
#' @param displays List of retained interaction displays.
#' @param digits Integer scalar; number of decimal places.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
mfpi_print_summary_regression_displays <- function(displays, digits) {
  first <- TRUE
  shift_note_printed <- FALSE

  for (term in names(displays)) {
    display <- displays[[term]]
    if (is.null(display)) next

    if (!first) cat("\n")

    show_shift_note <- !shift_note_printed &&
      mfpi_summary_display_has_nonzero_shift(display)

    mfpi_print_summary_regression_block(
      display,
      digits = digits,
      print_shift_note = show_shift_note
    )

    if (show_shift_note) shift_note_printed <- TRUE
    first <- FALSE
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# summary.mfpi() --------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Summarize an \code{"mfpi"} Object
#'
#' Produces a structured summary of an \code{"mfpi"} fit.
#'
#' @details
#' The printed summary has four steps. Steps 1--3 are the same model-building
#' display produced by [print.mfpi()]: the selected adjustment model, all
#' interaction candidates, and the final interaction-selection summary. Step 4
#' reports regression results for retained interaction models using MFPI-specific
#' labels rather than the raw output of the underlying regression model.
#'
#' The main MFPI header reports the FLEX setting, while Step 2 reports the
#' grouping variable and interaction-selection settings. Step 4 avoids
#' repeating those fit-level properties within each retained model. Instead,
#' each retained interaction is introduced by an
#' underlined interaction heading and reports its FP form and interaction powers
#' for every group level. Regression coefficients are then separated into:
#'
#' \enumerate{
#'   \item group-specific FP terms, displayed using the original group labels
#'     and readable FP transformations;
#'   \item group-variable coefficients, with the fitted reference level stated
#'     explicitly;
#'   \item a model intercept when the underlying model contains one; and
#'   \item adjustment coefficients.
#' }
#'
#' Centered basis terms report their stored centering constants using the same
#' \code{Basis - Center} convention as \code{print.mfp2()}. Multinomial models
#' use one outcome-versus-reference block per non-reference outcome, retaining
#' the same coefficient-table separation within every block. Each coefficient
#' table includes the estimate, standard error, model-specific
#' test statistic (for example, \code{z} or \code{t}), and p-value. Significance
#' stars, exponentiated coefficients, confidence-interval tables, the raw model
#' call, and model-level likelihood-ratio/Wald/score tests are intentionally not
#' printed. In particular, the usual global model tests condition on the
#' selected regression specification and their standard degrees of freedom do
#' not represent the additional model-selection complexity introduced by FP
#' power selection.
#'
#' Group-variable coefficients are coefficients of the fitted regression
#' parameterization. When a group-by-continuous interaction is present, they are
#' not overall comparisons between group levels; group comparisons depend on
#' the interacting variable and should be obtained from [predict.mfpi()].
#'
#' Coefficient-level p-values are also conditional on the selected functional
#' forms and do not account for uncertainty introduced by adjustment-model, FP,
#' or interaction selection. They should therefore not be used for the MFPI
#' interaction decision. The criterion reported in Step 3 is the relevant MFPI
#' selection result.
#'
#' @param object An object of class \code{"mfpi"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.mfpi"}. In addition to the fields
#'   needed to reproduce the MFPI Steps 1--3 display, the object contains
#'   \code{regression_displays}, a named list of structured coefficient displays
#'   for retained interaction models. Each display contains readable coefficient
#'   metadata, the corresponding model-conditional coefficient statistics, and
#'   model-specific nuisance metadata. Retained negative-binomial models report
#'   their own theta, while retained `survreg` models report their own scale(s).
#'   The normalized \code{family_string} is retained. Multinomial summaries also
#'   retain \code{class_levels}, \code{reference_class}, and \code{n_logits};
#'   ordinal summaries retain \code{ordinal_levels}, \code{ordinal_link}, and
#'   \code{n_intercepts}.
#'   Printing the object produces the four-step summary described above.
#'
#' @examples
#' \donttest{
#' data("prostate", package = "mfp2")
#'
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(cavol),
#'   data = prostate,
#'   group_var = "svi",
#'   interaction_vars = "cavol",
#'   interaction_forms = c(cavol = "fp2"),
#'   flex = "flex4",
#'   p_interact = 1,
#'   verbose = FALSE
#' )
#'
#' summary(fit)
#' }
#'
#' @seealso [mfpi()], [print.mfpi()], [coef.mfpi()], [vcov.mfpi()],
#'   [predict.mfpi()], [plot.mfpi()]
#'
#' @method summary mfpi
#' @export
summary.mfpi <- function(object, ...) {
  warn_unused_mfpi_dots(list(...), method = "summary.mfpi")

  if (!inherits(object, "mfpi")) {
    stop("`object` must be an object of class \"mfpi\".", call. = FALSE)
  }

  regression_displays <- build_mfpi_regression_displays(object)

  structure(
    list(
      # Basic model metadata ---------------------------------------------------
      call                    = object$call,
      group_var               = object$group_var,
      nobs                    = object$nobs,
      nevents                 = object$nevents,
      family                  = object$family,
      family_string           = object$family_string,
      flex                    = object$flex,
      criterion               = object$criterion,
      p_adjust_method         = object$p_adjust_method,
      p_interact              = object$p_interact,
      min_improvement         = object$min_improvement,
      digits                  = object$digits,

      # Group-level metadata --------------------------------------------------
      group_levels_new        = object$group_levels_new,
      group_levels_original   = object$group_levels_original,
      group_level_map         = object$group_level_map,

      # Ordinal response metadata --------------------------------------------
      ordinal_levels          = object$ordinal_levels,
      ordinal_link            = object$ordinal_link,
      n_intercepts            = object$n_intercepts,
      class_levels            = object$class_levels,
      reference_class         = object$reference_class,
      n_logits                = object$n_logits,

      # Model-building outputs ------------------------------------------------
      adjust_terms            = object$adjust_terms,
      all_model_metrics       = object$all_model_metrics,
      best_model_metrics      = object$best_model_metrics,
      var_winners             = object$var_winners,

      # Fitted interaction model objects --------------------------------------
      best_interaction_model  = object$best_interaction_model,
      all_interaction_models  = object$all_interaction_models,

      # MFPI-specific retained regression displays ----------------------------
      regression_displays     = regression_displays
    ),
    class = "summary.mfpi"
  )
}


# -----------------------------------------------------------------------------
# print.summary.mfpi() --------------------------------------------------------
# -----------------------------------------------------------------------------

#' Print a \code{"summary.mfpi"} Object
#'
#' Displays the four-step MFPI summary produced by [summary.mfpi()].
#'
#' @details
#' Steps 1--3 are delegated to [print.mfpi()] so that \code{print(fit)} and
#' \code{print(summary(fit))} use the same adjustment-model and interaction-
#' selection output. Step 4 adds a compact regression display for each retained
#' interaction model.
#'
#' Step 4 uses readable group labels and FP transformations rather than internal
#' design-column names. The grouping variable is shown in Step 2 and the FLEX
#' setting in the main MFPI header, so neither is repeated for every retained
#' interaction.
#' Each retained interaction is introduced by an underlined heading. Interaction
#' powers are always printed by group level, including when the selected FLEX
#' strategy gives all groups the same powers. Group-variable coefficients and
#' adjustment coefficients are displayed in separate tables. For models with an
#' intercept, it is displayed separately. Centered terms include a \code{Center}
#' column and use the \code{Basis - Center} convention. Multinomial models use
#' a separate outcome-versus-reference block for each non-reference outcome.
#'
#' When one or more retained interaction variables have a nonzero MFPI shift,
#' Step 4 prints one clarification immediately below the first shifted
#' interaction's group-specific FP table. The constant shown inside an FP
#' transformation is the shifting factor applied to all observations before the
#' group-specific FP terms are constructed; structural zeros remain
#' untransformed. Later shifted interactions retain the explicit constant in
#' their transformation labels, but the note is not repeated. No shift note is
#' printed when all retained interactions have zero shift.
#'
#' The original \code{mfpi()} call is printed with Steps 1--3. Step 4 does not
#' repeat the calls of the underlying retained regression models or print
#' exponentiated-coefficient tables, significance stars, concordance, or global
#' likelihood-ratio/Wald/score tests. Those quantities remain available from the
#' stored fitted regression model when required, but they are not part of the
#' MFPI interaction-selection result.
#'
#' @param x An object of class \code{"summary.mfpi"}.
#' @param digits Optional non-negative integer controlling displayed decimal
#'   places other than fractional-polynomial powers and degrees of freedom. If
#'   \code{NULL}, \code{x$digits} is used when available; otherwise the MFPI
#'   print-helper default is used.
#' @param ... Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @seealso [summary.mfpi()], [print.mfpi()], [coef.mfpi()], [vcov.mfpi()]
#'
#' @method print summary.mfpi
#' @export
print.summary.mfpi <- function(x, digits = NULL, ...) {
  warn_unused_mfpi_dots(list(...), method = "print.summary.mfpi")

  if (!inherits(x, "summary.mfpi")) {
    stop("`x` must be an object of class \"summary.mfpi\".", call. = FALSE)
  }

  # Reuse print.mfpi() for Steps 1-3. The summary object stores the fields
  # consumed by print.mfpi(), so a temporary class switch is sufficient.
  obj <- x
  class(obj) <- "mfpi"
  print.mfpi(obj, digits = digits)

  digits <- resolve_print_digits(x, digits)
  ruler <- strrep("-", 65)
  regression_step_no <- mfpi_summary_regression_step_number(x$criterion)

  cat(sprintf(
    "\nStep %d - Regression Output for Retained Interaction Models:\n",
    regression_step_no
  ))
  cat(ruler, "\n\n")

  displays <- x$regression_displays

  if (!has_retained_regression_displays(displays)) {
    cat("  No retained interaction models.\n")
    return(invisible(x))
  }

  mfpi_print_summary_regression_displays(displays, digits = digits)

  cat("\nNote: Each selected interaction is fitted in its own model. Use the results\n")
  cat("      in Step 3 for the MFPI interaction decision.\n")

  invisible(x)
}
