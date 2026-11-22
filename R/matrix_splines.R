#' Specify splines for the matrix interface
#'
#' These constructors describe a raw numeric column in `mfp2(x, y, splines = ...)`.
#' They do not evaluate a basis. Knots omitted here are chosen from the rows
#' used for fitting, then held fixed during selection and prediction.
#'
#' @param df Maximum number of spline columns. If omitted, inferred from
#'   supplied knots, or defaults to 6 for `bs_spec()` and 4 for `rcs_spec()`.
#' @param knots For `bs_spec()`, strictly increasing interior knots; for
#'   `rcs_spec()`, the entire knot vector including its two boundary knots.
#' @param boundary_knots Optional pair of B-spline boundary knots; defaults
#'   to the range of the fitted predictor.
#' @return A specification for one entry of `splines` in `mfp2.default()`.
#' @export
bs_spec <- function(df = NULL, knots = NULL, boundary_knots = NULL) {
  structure(list(type = "bspline", df = df, knots = knots,
                 boundary_knots = boundary_knots), class = "mfp2_spline_spec")
}

#' @rdname bs_spec
#' @export
rcs_spec <- function(df = NULL, knots = NULL) {
  structure(list(type = "rcs", df = df, knots = knots,
                 boundary_knots = NULL), class = "mfp2_spline_spec")
}

#' Resolve matrix spline requests and create the full fixed-knot design
#'
#' Expands only the specified raw columns. All other matrix columns keep their
#' original names and per-column options. The first spline column is the raw
#' predictor, so the existing term prediction and plotting code can display
#' values on their original scale.
#' @param x Raw numeric predictor matrix.
#' @param splines Named character vector or named list of spline specifications.
#' @param subset Optional fitted-row index or logical mask.
#' @param df Existing matrix-interface FP df setting.
#' @param df_missing Whether the caller supplied `df` explicitly.
#' @param catzero_vars,spike_vars Existing zero-indicator requests.
#' @param term_groups Existing grouped linear-term specification.
#' @param min_saz_prop Minimum fitted-row proportion for both spike groups.
#' @return Expanded design, fixed spline specifications, and column mapping.
#' @keywords internal
#' @noRd
prepare_matrix_splines <- function(x, splines, subset, df, df_missing,
                                   catzero_vars, spike_vars, term_groups,
                                   min_saz_prop) {
  if (!is.matrix(x) || !is.numeric(x) || is.null(colnames(x)) ||
      anyNA(x) || any(!is.finite(x))) {
    stop("`splines` requires a finite numeric matrix with column names.", call. = FALSE)
  }
  original <- colnames(x)
  validate_predictor_names(original, object = "`x`")
  validate_variable_names(catzero_vars, "catzero_vars", original)
  validate_variable_names(spike_vars, "spike_vars", original)
  if (is.character(splines)) {
    splines <- lapply(unname(splines), function(type) {
      if (identical(type, "bspline")) bs_spec() else if (identical(type, "rcs"))
        rcs_spec() else stop("`splines` types must be 'bspline' or 'rcs'.", call. = FALSE)
    }) |> stats::setNames(names(splines))
  }
  if (!is.list(splines) || !length(splines) || is.null(names(splines)) ||
      anyNA(names(splines)) || any(!nzchar(names(splines))) ||
      anyDuplicated(names(splines)) ||
      !all(vapply(splines, inherits, logical(1L), "mfp2_spline_spec"))) {
    stop("`splines` must be a named vector of spline types or a named list of bs_spec()/rcs_spec() values.", call. = FALSE)
  }
  if (!is.null(term_groups)) normalize_term_groups(original, term_groups)
  unknown <- setdiff(names(splines), original)
  if (length(unknown)) stop("Unknown spline column(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  if (length(intersect(names(splines), unlist(term_groups, use.names = FALSE)))) {
    stop("A spline column cannot also appear in `term_groups`.", call. = FALSE)
  }
  if (is.null(subset)) {
    rows <- seq_len(nrow(x))
  } else if (is.logical(subset) && length(subset) == nrow(x) && !anyNA(subset)) {
    rows <- which(subset)
  } else if (is.numeric(subset) && !anyNA(subset) && all(is.finite(subset)) &&
             all(subset == as.integer(subset)) && all(subset >= 1L) &&
             all(subset <= nrow(x)) && !anyDuplicated(subset)) {
    rows <- as.integer(subset)
  } else stop("Invalid `subset` for matrix spline fitting.", call. = FALSE)
  if (!length(rows)) stop("`subset` selects no observations.", call. = FALSE)
  df_setting <- normalize_named_override_setting(
    df, original, default = 4L, argument_name = "df"
  )
  specs <- list()
  groups <- term_groups
  expanded <- list()
  sources <- character(0L)
  for (term in original) {
    if (!term %in% names(splines)) {
      expanded[[term]] <- x[, term]
      sources <- c(sources, term)
      next
    }
    request <- splines[[term]]
    type <- request$type
    if (!type %in% c("bspline", "rcs")) stop("Invalid spline type.", call. = FALSE)
    supplied <- request$knots
    if (!is.null(supplied) && (!is.numeric(supplied) || anyNA(supplied) ||
                              any(!is.finite(supplied)))) {
      stop("Spline knots must be finite numeric values.", call. = FALSE)
    }
    inferred_df <- if (is.null(supplied)) NULL else if (type == "bspline")
      length(supplied) + 3L else length(supplied) - 1L
    user_df <- if (!is.null(request$df)) request$df else if (!df_missing &&
      (isTRUE(df_setting$global_scalar) || isTRUE(df_setting$supplied[[term]])))
        df_setting$value[[term]] else if (!is.null(inferred_df)) inferred_df else
          if (type == "bspline") 6L else 4L
    min_df <- if (type == "bspline") 4L else 2L
    if (!is.numeric(user_df) || length(user_df) != 1L || is.na(user_df) ||
        !is.finite(user_df) || user_df < min_df || user_df != as.integer(user_df)) {
      stop("Invalid maximum spline `df` for `", term, "`.", call. = FALSE)
    }
    if (!is.null(supplied) && inferred_df != user_df) {
      stop("Spline knot count and `df` disagree for `", term, "`.", call. = FALSE)
    }
    fitted <- x[rows, term]
    requested_catzero <- term %in% catzero_vars
    requested_spike <- term %in% spike_vars
    if ((requested_catzero || requested_spike) && any(x[, term] < 0)) {
      stop("Zero-indicator spline `", term, "` requires nonnegative values.", call. = FALSE)
    }
    has_zero <- any(fitted == 0)
    binary <- length(unique(fitted)) <= 2L
    eligible <- requested_spike && !binary && has_zero &&
      mean(fitted == 0) >= min_saz_prop && mean(fitted > 0) >= min_saz_prop
    if (requested_spike && !eligible) {
      warning("spike was reset for spline `", term, "`; explicit catzero is preserved.",
              call. = FALSE)
    }
    if (requested_catzero && (!has_zero || binary)) {
      warning("catzero was reset for spline `", term, "`.", call. = FALSE)
    }
    use_indicator <- (requested_catzero && has_zero && !binary) || eligible
    knot_values <- if (use_indicator) fitted[fitted > 0] else fitted
    if (type == "bspline") {
      boundaries <- request$boundary_knots
      if (is.null(boundaries)) boundaries <- range(fitted)
      bspline_validate_boundaries(boundaries)
      knots <- if (is.null(supplied)) bspline_default_knots(knot_values, user_df, boundaries) else supplied
      bspline_validate_knots(knots, length(knots) + 3L, boundaries)
      effective_df <- length(knots) + 3L
      terms <- bspline_basis_terms(knots, boundaries)
      candidate_df <- lengths(bspline_knot_sets(knots)) + 3L
      spec <- list(type = type, df = effective_df, knots = knots,
                   boundary_knots = boundaries, basis_terms = terms,
                   candidate_df = candidate_df)
    } else {
      knots <- if (is.null(supplied)) rcs_default_knots(knot_values, user_df) else supplied
      rcs_validate_knots(knots, length(knots) - 1L)
      effective_df <- length(knots) - 1L
      spec <- list(type = type, df = effective_df, knots = knots,
                   boundary_knots = NULL, basis_terms = NULL,
                   candidate_df = seq.int(2L, effective_df))
    }
    if (effective_df < user_df) {
      warning(type, " spline `", term, "` has fewer distinct knots than requested; ",
              "maximum df reduced from ", user_df, " to ", effective_df, ".",
              call. = FALSE)
    }
    spec$catzero <- requested_catzero && has_zero && !binary
    spec$spike <- eligible
    spec$indicator <- use_indicator
    spec$label <- term
    specs[[term]] <- spec
    basis <- spline_basis_from_spec(x[, term], spec)
    colnames(basis) <- c(term, paste0(term, "__spline", seq.int(2L, effective_df)),
                         if (use_indicator) paste0(term, "__zero"))
    if (any(colnames(basis)[-1L] %in% original)) {
      stop("Generated spline column names collide with `x`.", call. = FALSE)
    }
    expanded <- c(expanded, as.list(as.data.frame(basis, check.names = FALSE)))
    sources <- c(sources, rep(term, ncol(basis)))
    groups[[term]] <- colnames(basis)
  }
  design <- as.matrix(as.data.frame(expanded, check.names = FALSE))
  storage.mode(design) <- "double"
  rownames(design) <- rownames(x)
  attr(design, "mfp2_spline_specs") <- specs
  if (anyDuplicated(colnames(design))) {
    stop("Generated spline column names are not unique.", call. = FALSE)
  }
  list(x = design, specs = specs, groups = groups,
       sources = stats::setNames(sources, colnames(design)),
       original_columns = original)
}

#' Translate matrix column settings to the expanded spline design
#' @param prepared Output from `prepare_matrix_splines()`.
#' @param df,alpha,select,center,shift,scale,powers Existing matrix settings.
#' @return Named per-column settings suitable for `mfp2.default()`.
#' @keywords internal
#' @noRd
matrix_spline_settings <- function(prepared, df, alpha, select, center,
                                   shift, scale, powers) {
  raw <- prepared$original_columns
  sources <- prepared$sources
  spline <- names(prepared$specs)
  expand <- function(values) stats::setNames(unname(values[sources]), names(sources))
  df_setting <- normalize_named_override_setting(df, raw, 4L, "df")
  df_values <- df_setting$value
  # Preserve the original matrix interface's cardinality limits for FP
  # predictors. Only the spline columns are subsequently forced to df=1 in
  # the expanded design; their actual maximum is stored in the spline spec.
  raw_design <- prepared$x[, raw, drop = FALSE]
  if (isTRUE(df_setting$global_scalar)) {
    df_values <- assign_df(raw_design, df_default = df_values)
  } else {
    defaulted <- !df_setting$supplied
    if (any(defaulted)) {
      defaults <- assign_df(raw_design, df_default = df_values)
      df_values[defaulted] <- defaults[defaulted]
    }
    df_values[vapply(seq_along(raw), function(i) {
      length(unique(raw_design[, i])) <= 3L
    }, logical(1L))] <- 1L
  }
  df_out <- expand(df_values)
  df_out[sources %in% spline] <- 1L
  # A pre-existing term_groups block is already a fixed linear term. A new
  # spline beside it must not turn its global FP default into an explicit df=4.
  for (group in setdiff(names(prepared$groups), spline)) {
    df_out[prepared$groups[[group]]] <- 1L
  }
  alpha_out <- expand(normalize_named_override_setting(alpha, raw, 0.05, "alpha")$value)
  select_out <- expand(normalize_named_override_setting(select, raw, 0.05, "select")$value)
  center_out <- expand(normalize_named_logical_setting(center, raw, TRUE, "center"))
  expand_numeric <- function(value, argument, spline_value) {
    normalized <- normalize_named_numeric_setting(
      value, raw, argument, allow_partial_named = TRUE
    )
    result <- expand(normalized)
    result[sources %in% spline] <- spline_value
    result
  }
  if (length(intersect(names(powers), spline))) {
    stop("`powers` cannot be specified for a spline column.", call. = FALSE)
  }
  list(df = df_out, alpha = alpha_out, select = select_out,
       center = center_out,
       shift = expand_numeric(shift, "shift", 0),
       scale = expand_numeric(scale, "scale", 1), powers = powers)
}

#' Recreate matrix spline columns from raw prediction values
#' @param object Fitted matrix-interface spline model.
#' @param newdata New raw matrix or data frame.
#' @param terms Active conceptual terms.
#' @return Design columns needed for prediction.
#' @keywords internal
#' @noRd
reconstruct_matrix_spline_newdata <- function(object, newdata, terms) {
  lookup <- prediction_term_to_columns(object, terms = terms)
  expected <- unlist(lookup, use.names = FALSE)
  newdata <- as.data.frame(newdata, check.names = FALSE)
  if (!length(expected)) return(data.frame(row.names = row.names(newdata)))
  if (all(expected %in% names(newdata))) return(newdata[, expected, drop = FALSE])
  result <- vector("list", length(lookup))
  names(result) <- names(lookup)
  for (term in names(lookup)) {
    columns <- lookup[[term]]
    spec <- object$matrix_spline_specs[[term]]
    if (is.null(spec)) {
      absent <- setdiff(columns, names(newdata))
      if (length(absent)) stop("Missing required predictor column(s) in newdata: ",
                               paste(absent, collapse = ", "), call. = FALSE)
      result[[term]] <- newdata[, columns, drop = FALSE]
      next
    }
    if (all(columns %in% names(newdata))) {
      result[[term]] <- newdata[, columns, drop = FALSE]
      next
    }
    if (!term %in% names(newdata)) {
      stop("Missing raw spline column `", term, "` in newdata.", call. = FALSE)
    }
    values <- newdata[[term]]
    mfp2_validate_prediction_newdata(newdata[, term, drop = FALSE])
    if (isTRUE(spec$indicator) && any(values < 0)) {
      stop("Spline `", term, "` with a zero indicator requires nonnegative newdata.",
           call. = FALSE)
    }
    basis <- spline_basis_from_spec(values, spec)
    colnames(basis) <- object$formula_term_to_columns[[term]]
    result[[term]] <- as.data.frame(basis[, columns, drop = FALSE], check.names = FALSE)
  }
  design <- data.frame(row.names = row.names(newdata))
  for (term in names(result)) {
    for (column in names(result[[term]])) {
      design[[column]] <- result[[term]][[column]]
    }
  }
  design[, expected, drop = FALSE]
}
