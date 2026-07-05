# Helpers for extracting information from mfp2 model objects
#
# get_selected_variables() -- returns names of variables retained by MFP
# get_fp_powers()          -- returns the FP power vector(s) for named variables


# -----------------------------------------------------------------------------
# get_selected_variables() ----------------------------------------------------
# -----------------------------------------------------------------------------

#' Extract Selected Variables from an MFP Model
#'
#' Returns the names of variables selected by an object of class \code{"mfp2"}.
#'
#' The \code{fp_terms} component is expected to be the data frame created by
#' the internal \code{create_fp_terms()} helper. Its row names are the original
#' variable names, and its logical \code{selected} column records whether each
#' variable was retained by the MFP algorithm.
#'
#' @param object Object of class \code{"mfp2"}.
#'
#' @return Character vector containing the row names of selected variables. If
#'   the \code{selected} column is absent, returns \code{character(0)}.
#'
#' @keywords internal
#' @noRd
get_selected_variables <- function(object) {
  
  if (!inherits(object, "mfp2")) {
    stop(
      "! `object` must inherit from class \"mfp2\".",
      call. = FALSE
    )
  }
  
  if (!"fp_terms" %in% names(object)) {
    stop(
      "! `object` must contain a component named \"fp_terms\".",
      call. = FALSE
    )
  }
  
  fp_terms <- object$fp_terms
  
  if (!is.data.frame(fp_terms)) {
    stop("! `fp_terms` must be a data frame.", call. = FALSE)
  }
  
  if (is.null(colnames(fp_terms)) || anyNA(colnames(fp_terms))) {
    stop("! `fp_terms` must have non-missing column names.", call. = FALSE)
  }
  
  # If the 'selected' column is absent, no variable survived selection.
  if (!"selected" %in% colnames(fp_terms)) {
    return(character(0L))
  }
  
  selected <- fp_terms[, "selected"]
  
  if (!is.logical(selected)) {
    stop("! The 'selected' column in `fp_terms` must be logical.", call. = FALSE)
  }
  
  rn <- rownames(fp_terms)
  
  if (is.null(rn) || length(rn) != nrow(fp_terms) || anyNA(rn) ||
      any(!nzchar(rn))) {
    stop(
      "! `fp_terms` must have non-empty row names giving variable names.",
      call. = FALSE
    )
  }
  
  rn[selected]
}


# -----------------------------------------------------------------------------
# get_fp_powers() -------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Extract FP Powers for Named Variables from an \code{fp_terms} Table
#'
#' Retrieves the fractional polynomial powers (\code{power1}, \code{power2},
#' …) for one or more variables from the \code{fp_terms} data frame of an
#' \code{"mfp2"} object. \code{NA} values in the power columns indicate that
#' a variable was eliminated during model selection.
#'
#' @param varname Character vector of variable names. Each name must appear as
#'   a row name in \code{df}.
#' @param df A data frame, typically \code{object$fp_terms} from an
#'   \code{"mfp2"} object. Must contain at least one column named
#'   \code{power1}; additional power columns (\code{power2}, etc.) are
#'   included if present.
#' @param warn_if_not_selected Logical. If \code{TRUE}, a warning is issued
#'   for any variable whose power columns are all \code{NA} (i.e. the variable
#'   was not selected). Default is \code{FALSE} because \code{get_fp_powers()}
#'   is typically called after \code{get_selected_variables()} has already
#'   filtered to retained variables, making the warning redundant. Set to
#'   \code{TRUE} when calling on an unfiltered variable list.
#'
#' @return A named list with one element per variable in \code{varname}. Each
#'   element is a numeric vector of FP powers with \code{NA}s removed. For
#'   a linear term the vector is \code{1}; for FP1 it has length 1 (e.g.
#'   \code{0.5}); for FP2 it has length 2 (e.g. \code{c(-1, 0.5)}). If all
#'   powers are \code{NA} (variable not selected), the element is
#'   \code{numeric(0)}.
#'
#' @examples
#' \dontrun{
#' fit    <- mfp2(y ~ fp(age) + fp(bmi), data = mydata)
#' powers <- get_fp_powers(c("age", "bmi"), fit$fp_terms)
#'
#' # Warn when a variable was not selected
#' powers <- get_fp_powers("age", fit$fp_terms, warn_if_not_selected = TRUE)
#' }
#'
#' @keywords internal
#' @noRd
get_fp_powers <- function(varname, df, warn_if_not_selected = FALSE) {
  
  if (!is.character(varname) || length(varname) == 0L) {
    stop("! `varname` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.data.frame(df)) {
    stop("! `df` must be a data frame.", call. = FALSE)
  }
  
  missing_vars <- setdiff(varname, rownames(df))
  if (length(missing_vars) > 0L) {
    stop(
      paste0("! The following variable(s) are not found as row names in `df`: ",
             paste(missing_vars, collapse = ", "), "."),
      call. = FALSE
    )
  }
  
  power_cols <- grep("^power\\d+$", names(df), value = TRUE)
  if (length(power_cols) == 0L) {
    stop(
      "! No power columns (e.g. 'power1', 'power2') found in `df`.",
      call. = FALSE
    )
  }
  
  result <- lapply(varname, function(v) {
    raw <- unlist(df[v, power_cols, drop = FALSE], use.names = FALSE)
    
    if (all(is.na(raw))) {
      if (warn_if_not_selected) {
        warning(
          paste0("i Variable '", v, "' was not selected in the model ",
                 "(all FP powers are NA)."),
          call. = FALSE
        )
      }
      return(numeric(0L))
    }
    
    # Drop trailing NAs: FP1 has power1 only; power2 is NA and not needed
    raw[!is.na(raw)]
  })
  
  setNames(result, varname)
}

#' Bind metric rows into a base data frame
#'
#' Combines a list of metric data frames into a single base data frame while
#' preserving list-columns such as fractional polynomial power columns. This is
#' an internal base-R replacement for `dplyr::bind_rows()` used for MFPI model
#' metric objects.
#'
#' The function is designed for metric rows that may not all have exactly the
#' same columns. Missing columns are added before binding. Missing scalar
#' columns are filled with `NA`, while missing list-columns are filled with
#' empty list elements of the appropriate length.
#'
#' @param rows A list of data frames. `NULL` elements, non-data-frame elements,
#'   and zero-row data frames are ignored.
#' @param metric_class Character scalar or `NULL`. Optional class to prepend to
#'   the returned data frame. The default is `"best_model_metrics"`. If `NULL`,
#'   no custom class is added.
#'
#' @return A base data frame containing all valid rows in `rows`. If
#'   `metric_class` is not `NULL`, the returned object has class
#'   `c(metric_class, "data.frame")`. List-columns are preserved.
#'
#' @keywords internal
#' @noRd
bind_metric_rows <- function(rows, metric_class = "best_model_metrics") {
  # Keep only valid, non-empty data frames.
  rows <- Filter(function(x) {
    !is.null(x) && is.data.frame(x) && nrow(x) > 0L
  }, rows)
  
  # If there are no valid rows, return an empty data frame with the requested
  # custom class. This gives downstream code a predictable object type.
  if (length(rows) == 0L) {
    out <- data.frame()
    
    if (!is.null(metric_class)) {
      class(out) <- c(metric_class, "data.frame")
    }
    
    return(out)
  }
  
  # Strip custom classes before calling base rbind(). This avoids method
  # dispatch or class-specific behavior while combining rows.
  rows <- lapply(rows, function(x) {
    class(x) <- "data.frame"
    x
  })
  
  # Compute the union of all column names, preserving the order in which column
  # names first appear across the input data frames.
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  
  # Identify columns that should be treated as list-columns. A column is treated
  # as a list-column if it is a list in at least one input data frame, excluding
  # nested data-frame columns.
  list_cols <- vapply(all_names, function(nm) {
    any(vapply(rows, function(d) {
      nm %in% names(d) && is.list(d[[nm]]) && !is.data.frame(d[[nm]])
    }, logical(1L)))
  }, logical(1L))
  
  # Add missing columns to each input data frame, reorder columns consistently,
  # and mark list-columns with I() so base rbind() preserves them.
  rows <- lapply(rows, function(d) {
    missing <- setdiff(all_names, names(d))
    
    # Fill missing columns. List-columns need list placeholders; ordinary
    # columns can be filled with NA.
    if (length(missing) > 0L) {
      for (nm in missing) {
        if (isTRUE(list_cols[[nm]])) {
          d[[nm]] <- I(vector("list", nrow(d)))
        } else {
          d[[nm]] <- NA
        }
      }
    }
    
    # Reorder columns so every data frame has the same column layout before
    # binding.
    d <- d[, all_names, drop = FALSE]
    
    # Preserve list-columns during base rbind().
    for (nm in names(d)) {
      if (isTRUE(list_cols[[nm]])) {
        d[[nm]] <- I(d[[nm]])
      }
    }
    
    d
  })
  
  # Bind all rows using base R.
  out <- do.call(rbind, rows)
  
  # Remove row names introduced by rbind().
  rownames(out) <- NULL
  
  # Restore the requested custom class on the final combined data frame.
  if (!is.null(metric_class)) {
    class(out) <- c(metric_class, "data.frame")
  }
  
  out
}

#' Construct a long interaction-power display table
#'
#' Converts the interaction fractional polynomial powers stored in an MFPI
#' metrics table into a long-format data frame with one row per variable and
#' group level. This table is intended for printing only; it does not modify the
#' stored MFPI object.
#'
#' The helper supports both current and future storage formats. Interaction
#' powers may be stored as unnamed lists, named lists, or simple numeric vectors.
#' When names are available on `fp_powers_int`, they are used as group-level
#' labels. Otherwise generic labels such as `"level_1"` are used.
#'
#' @param d A data frame-like object containing MFPI model metrics. Usually
#'   `x$all_model_metrics` or `x$best_model_metrics`.
#' @param group_var Character scalar or `NULL`. Name of the grouping variable.
#'   If supplied, it is printed in the `group_var` column.
#'
#' @return A base data frame with columns `variable`, `type`, `group_var`,
#'   `group_level`, and `fp_powers_int`. If no interaction
#'   powers are available, an empty data frame is returned.
#'
#' @keywords internal
#' @noRd
mfpi_interaction_power_table <- function(d, group_var = NULL) {
  # Return an empty table if no metrics or interaction powers are available.
  if (is.null(d) ||
      nrow(d) == 0L ||
      !("fp_powers_int" %in% names(d))) {
    return(data.frame())
  }
  
  # Work with a base data frame for predictable indexing and printing.
  d <- as.data.frame(d)
  
  # Build one or more display rows for each metric row.
  rows <- lapply(seq_len(nrow(d)), function(i) {
    # Extract the interaction power specification for this metric row.
    p_int <- d$fp_powers_int[[i]]
    
    # Skip rows with no interaction powers.
    if (is.null(p_int) || length(p_int) == 0L) {
      return(NULL)
    }
    
    # Normalise a simple numeric vector to a one-element list. This keeps the
    # rest of the code list-based.
    if (!is.list(p_int)) {
      p_int <- list(p_int)
    }
    
    # Use group-level names when present. If unavailable, create generic labels.
    group_levels <- names(p_int)
    if (is.null(group_levels) ||
        length(group_levels) != length(p_int) ||
        any(!nzchar(group_levels))) {
      group_levels <- paste0("level_", seq_along(p_int))
    }
    
    # Extract the variable name, if present.
    variable <- if ("variable" %in% names(d)) {
      as.character(d$variable[i])
    } else {
      NA_character_
    }
    
    # Extract and format the interaction candidate type, if present.
    type <- if ("type" %in% names(d)) {
      format_type_label(d$type[i])
    } else {
      NA_character_
    }
    
    # Create one row per group-specific interaction power vector.
    data.frame(
      type          = type,
      variable      = variable,
      group_var     = if (!is.null(group_var)) group_var else NA_character_,
      group_level   = group_levels,
      fp_powers_int = vapply(p_int, format_power_vector, character(1L)),
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  })
  
  # Drop skipped rows.
  rows <- Filter(Negate(is.null), rows)
  
  # If all rows were skipped, return an empty data frame.
  if (length(rows) == 0L) {
    return(data.frame())
  }
  
  # Bind the display rows using base R.
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  
  out
}


#' Normalise main-effect FP powers for metric storage
#'
#' Converts a main-effect fractional polynomial power specification into the
#' standard named-list representation used in MFPI metric tables.
#'
#' The target representation is a named list with one element:
#' `list(variable = powers)`. This keeps the stored metric object
#' machine-readable while preserving the association between the continuous
#' variable and its selected fractional polynomial powers.
#'
#' @param powers Main-effect FP powers. Usually a numeric vector, but older or
#'   intermediate objects may already contain a list.
#' @param variable Character scalar. Name of the continuous variable.
#'
#' @return A named list with one element named by `variable`.
#'
#' @keywords internal
#' @noRd
normalise_main_fp_powers <- function(powers, variable) {
  # Store absent powers explicitly under the variable name.
  if (is.null(powers) || length(powers) == 0L) {
    return(stats::setNames(list(NULL), variable))
  }
  
  # If the object is already stored as list(variable = powers), keep it.
  if (is.list(powers) &&
      length(powers) == 1L &&
      !is.null(names(powers)) &&
      identical(names(powers), variable)) {
    return(powers)
  }
  
  # If a one-element list is supplied, unwrap it to the underlying power vector.
  # This handles objects such as list(c(1, 2)).
  if (is.list(powers) && length(powers) == 1L) {
    powers <- powers[[1L]]
  }
  
  # Store the main-effect powers as a one-element named list.
  stats::setNames(list(powers), variable)
}


#' Normalise interaction FP powers for metric storage
#'
#' Converts an interaction fractional polynomial power specification into the
#' standard named-list representation used in MFPI metric tables.
#'
#' The target representation is a named list with one element per group level.
#' Names should ideally be the original group labels. If the number of power
#' vectors suggests a non-reference-only representation, the reference level is
#' omitted from the labels. If original group labels cannot be matched, existing
#' valid names are preserved. If no valid names are available, generic labels are
#' generated.
#'
#' @param powers Interaction FP powers. Usually a list of numeric vectors, one
#'   per group level. Simpler or older objects may contain a plain numeric
#'   vector.
#' @param group_levels Optional vector of original group levels.
#'
#' @return A named list of FP power vectors.
#'
#' @keywords internal
#' @noRd
normalise_interaction_fp_powers <- function(powers, group_levels = NULL) {
  # No interaction powers available.
  if (is.null(powers) || length(powers) == 0L) {
    return(list())
  }
  
  # A plain numeric vector represents one interaction power specification.
  # Wrap it in a list so all downstream code can use list-based logic.
  if (!is.list(powers)) {
    powers <- list(powers)
  }
  
  # Resolve one readable label per interaction power vector. Prefer original
  # group labels, but tolerate non-reference-only or already named inputs.
  labels <- resolve_interaction_power_labels(
    powers       = powers,
    group_levels = group_levels
  )
  
  # Apply labels to the list-column entry.
  names(powers) <- labels
  
  powers
}


#' Normalise FP power columns in an MFPI metric row
#'
#' Ensures that the `fp_powers_main` and `fp_powers_int` list-columns of an MFPI
#' metric row use named-list storage. This helper is intended to be applied
#' after an interaction candidate has been fitted and before its metric row is
#' stored in `all_model_metrics`, `best_model_metrics`, or `var_winners`.
#'
#' Only the reporting-oriented metric row is modified. Fitted model internals,
#' design matrices, coefficient names, and prediction-related objects are left
#' unchanged.
#'
#' @param metrics A one-row MFPI metric data frame.
#' @param variable Character scalar. Name of the continuous variable.
#' @param group_levels Optional vector of original group levels.
#'
#' @return The input metric data frame with normalised FP power list-columns.
#'
#' @keywords internal
#' @noRd
normalise_metric_fp_powers <- function(metrics, variable, group_levels = NULL) {
  # Leave invalid or empty metric objects unchanged.
  if (is.null(metrics) || !is.data.frame(metrics) || nrow(metrics) == 0L) {
    return(metrics)
  }
  
  # Normalise the main-effect FP power column, if present.
  if ("fp_powers_main" %in% names(metrics)) {
    main_powers <- metrics$fp_powers_main[[1L]]
    
    metrics$fp_powers_main <- I(list(
      normalise_main_fp_powers(
        powers   = main_powers,
        variable = variable
      )
    ))
  }
  
  # Normalise the interaction FP power column, if present.
  if ("fp_powers_int" %in% names(metrics)) {
    int_powers <- metrics$fp_powers_int[[1L]]
    
    metrics$fp_powers_int <- I(list(
      normalise_interaction_fp_powers(
        powers       = int_powers,
        group_levels = group_levels
      )
    ))
  }
  
  metrics
}

#' Resolve group labels for interaction FP powers
#'
#' Determines the group-level labels to use for a list of interaction fractional
#' polynomial power vectors. The preferred labels are the original group levels,
#' but the function also supports non-reference-only power lists and already
#' named power lists.
#'
#' If the number of original group levels matches the number of power vectors,
#' all original group levels are used. If there is one fewer power vector than
#' original group levels, the first group level is treated as the reference and
#' omitted from the labels. If neither case applies, existing valid names are
#' preserved; otherwise generic labels are generated.
#'
#' @param powers A list of interaction FP power vectors.
#' @param group_levels Optional vector of original group levels.
#'
#' @return A character vector with one label per element of `powers`.
#'
#' @keywords internal
#' @noRd
resolve_interaction_power_labels <- function(powers, group_levels = NULL) {
  # Number of group-specific power vectors to label.
  n_powers <- length(powers)
  
  # Empty input has no labels.
  if (n_powers == 0L) {
    return(character(0L))
  }
  
  # Prefer original group labels when there is one power vector per group level.
  # This is the expected representation for the current MFPI interaction design.
  if (!is.null(group_levels) && length(group_levels) == n_powers) {
    return(as.character(group_levels))
  }
  
  # Support a possible non-reference-only representation. In that case the first
  # group level is treated as the reference and omitted from the interaction
  # power labels.
  if (!is.null(group_levels) &&
      length(group_levels) > 1L &&
      length(group_levels) - 1L == n_powers) {
    return(as.character(group_levels[-1L]))
  }
  
  # Preserve existing names if they are complete and non-empty.
  existing_names <- names(powers)
  if (!is.null(existing_names) &&
      length(existing_names) == n_powers &&
      all(nzchar(existing_names))) {
    return(existing_names)
  }
  
  # Fall back to generic labels.
  paste0("level_", seq_len(n_powers))
}

#' Check for ordered-factor variables in formula interfaces
#'
#' Ordered factors are unsafe in the formula interfaces because model.matrix()
#' expands them using polynomial contrasts by default, creating columns such as
#' .L, .Q, and .C. The MFP/MFPI engines would then treat those contrast columns
#' as independent variables rather than one conceptual predictor.
#'
#' @param vars Character vector of variable names to check.
#' @param data Data frame containing the variables.
#' @param interface Character label used in the error message.
#'
#' @keywords internal
#' @noRd
check_ordered_factor_variables <- function(vars,
                                           data,
                                           interface = "model") {
  if (is.null(data)) {
    return(invisible(NULL))
  }
  
  vars <- unique(vars)
  vars <- vars[!is.na(vars) & nzchar(vars)]
  vars <- intersect(vars, names(data))
  
  if (length(vars) == 0L) {
    return(invisible(NULL))
  }
  
  is_ordered <- vapply(
    vars,
    function(v) {
      is.ordered(data[[v]])
    },
    logical(1L)
  )
  
  ordered_vars <- vars[is_ordered]
  
  if (length(ordered_vars) == 0L) {
    return(invisible(NULL))
  }
  
  stop(
    sprintf(
      paste0(
        "`%s()` does not support ordered-factor variables: %s.\n",
        "Ordered factors are expanded by `model.matrix()` into polynomial ",
        "contrast columns such as `.L`, `.Q`, and `.C`. The MFP/MFPI engine ",
        "would then treat those contrast columns as separate variables, which ",
        "breaks variable-wise selection, degree testing, and group-specific ",
        "interaction construction.\n",
        "Convert ordered factors to unordered factors, numeric scores, or ",
        "explicit dummy/ordinal indicator variables before fitting."
      ),
      interface,
      paste(sQuote(ordered_vars), collapse = ", ")
    ),
    call. = FALSE
  )
  
  invisible(NULL)
}


#' Extract a Named Adjustment-Model Flag
#'
#' Extracts a named logical or numeric flag vector from the fitted adjustment
#' model and aligns it to the selected adjustment variables.
#'
#' This is used after the adjustment model has been fitted because fit_mfp()
#' may update zero, catzero, and spike flags during preprocessing or
#' spike-at-zero checks. The interaction-evaluation stage must therefore use
#' the final post-fit flags stored in the adjustment model, not the original
#' pre-fit user inputs.
#'
#' @param object Fitted adjustment model object.
#' @param field Character scalar naming the field to extract.
#' @param selected_vars Character vector of selected adjustment variables.
#' @param default Default value used when the field is absent.
#'
#' @return A named vector aligned to `selected_vars`.
#'
#' @keywords internal
#' @noRd
mfpi_extract_adjustment_flag <- function(object,
                                         field,
                                         selected_vars,
                                         default = FALSE) {
  if (!is.character(field) || length(field) != 1L ||
      is.na(field) || !nzchar(field)) {
    stop("`field` must be a single non-empty character string.", call. = FALSE)
  }
  
  selected_vars <- as.character(selected_vars)
  
  flag <- object[[field]]
  
  if (is.null(flag)) {
    return(stats::setNames(
      rep(default, length(selected_vars)),
      selected_vars
    ))
  }
  
  if (is.null(names(flag))) {
    stop(
      paste0(
        "! Internal error: `adjustment_model$", field,
        "` must be named."
      ),
      call. = FALSE
    )
  }
  
  missing_flag <- setdiff(selected_vars, names(flag))
  
  if (length(missing_flag) > 0L) {
    stop(
      paste0(
        "! Internal error: `adjustment_model$", field,
        "` is missing selected variable(s): ",
        paste(missing_flag, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  out <- flag[selected_vars]
  names(out) <- selected_vars
  out
}
