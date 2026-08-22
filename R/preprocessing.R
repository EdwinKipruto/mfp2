#' Calculate a scaling factor for a predictor variable
#' 
#' @details
#' This function exposes one preprocessing calculation used by [mfp2()].
#' Most users do not need to call it directly because shifting and scaling are
#' handled automatically during model fitting.
#' 
#' For details on why scaling is useful, see the corresponding section in the
#' documentation of \code{mfp2()}.
#' 
#' The determination of the scaling factor is independent of (i.e. not affected 
#' by) shifts in the input data, as it depends only on the range of the 
#' input data.
#' 
#' Note that the estimation of powers is unaffected by scaling; the same powers 
#' are found for scaled input data. In extreme cases, scaling is necessary to 
#' preserve accuracy; see Royston and Sauerbrei (2008).
#' This function uses the scaling formula from Section 4.11.1 of 
#' Royston and Sauerbrei (2008). Further information can also be found in the 
#' Stata manual for mfp at https://www.stata.com/manuals/rfp.pdf.
#'
#' @param x A numeric vector already shifted to positive values (see 
#' \code{find_shift_factor()}). Requires at least two distinct values.
#' 
#' @examples
#' x = 1:1000
#' find_scale_factor(x)
#' 
#' @return 
#' An integer that can be used to scale `x` to a reasonable range. For binary
#' variables, 1 is returned.
#' 
#' @references 
#' Royston, P., and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @keywords internal
#' @noRd
find_scale_factor <- function(x) {
  
  n_unique <- length(unique(x))
  
  if (n_unique == 1)
    stop("! Input data must not be constant.", 
         "i All values of x are identical, hence log(max(x)-min(x)) = log(0) is not defined.")
  
  if (n_unique == 2)
    return(1)
  
  lrange <- log10(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  
  10^(sign(lrange) * floor(abs(lrange)))
}

#' Calculate a shift factor for a predictor variable
#' 
#' @details
#' This function exposes one preprocessing calculation used by [mfp2()].
#' Most users do not need to call it directly because shifting and scaling are
#' handled automatically during model fitting.
#' For details on why shifting is necessary, see the corresponding section in the
#' documentation of \code{mfp2()}.
#' 
#' This function implements the formula in Section 4.7 of Royston and 
#' Sauerbrei (2008).
#' 
#' @param x A numeric vector.
#' 
#' @examples
#' x = 1:1000
#' find_shift_factor(x)
#' 
#' @return 
#' A numeric value that can be used to shift `x` to positive values. 
#' If all values are positive, or if `x` is binary, then 0 is returned.
#' 
#' @references 
#' Royston, P., and Sauerbrei, W., 2008. \emph{Multivariable Model - Building: 
#' A Pragmatic Approach to Regression Analysis based on Fractional Polynomials 
#' for Modelling Continuous Variables. John Wiley & Sons.}
#' 
#' @keywords internal
#' @noRd
find_shift_factor <- function(x) {
  
  n_unique <- length(unique(x))
  
  if (all(x > 0) || n_unique <= 2) 
    return(0)
  
  difx <- diff(sort(x))
  eps <- min(difx[difx != 0], na.rm = TRUE)
  
  eps - min(x, na.rm = TRUE)
}

#' Shift and scale a predictor vector
#' 
#' Shifts `x` to positive values if it contains negative or zero values. If all
#' values of `x` are already positive, the original values are returned without
#' shifting but scaled if the scaling factor is not equal to 1. If `x` has
#' already been shifted and scaled, the function returns it unchanged.
#' 
#' @param x A numeric vector of predictor values.
#' @param scale Scaling factor for `x`. Must be a positive numeric value.
#' Default is `NULL`, in which case the scaling factor is estimated
#' automatically using `find_scale_factor()`. Set `scale = 1` to disable
#' scaling.
#' @param shift Shift factor for `x`. Default is `NULL`, in which case the
#' shift is estimated automatically using `find_shift_factor()`.
#' @examples
#' x = 1:1000
#' apply_shift_scale(x)
#' 
#' @returns 
#' A numeric vector of the same length as `x`, shifted and scaled.
#'  
#' @keywords internal
#' @noRd
apply_shift_scale <- function(x, scale = NULL, shift = NULL) {
  # restrict x to be a vector not matrix
  if (is.matrix(x))
    stop("x must be a vector not a matrix", call. = FALSE)
  
  N <- length(x)
  # If adjustment factors are NULL we use R&S formula to shift x to positive values
  if (is.null(shift)) {
    # Check whether all x values are positive. If true-No need of shifting
    if (all(x > 0)) {
      x <- x
    } else {
      # estimate adjustment factors
      shift <- find_shift_factor(x)
      # Shift x to positive
      x <- x + shift
    }
    # use the adjustment factors supplied by the user to shift x to positive values
  } else {
    # term + a
    x <- x + shift
    # check whether all x are now positive
    if (!all(x > 0))
      stop(
        "The minimum value of x after shifting x is ",
        min(x, na.rm = TRUE),
        " which is not > 0. Check your adjustment factors",
        call. = FALSE
      )
  }
  
  # if scale is NULL then scale x for computational stability using R&S formula
  if (is.null(scale)) { # No scaling
    x <- x / find_scale_factor(x)
  } else {
    x <- x / scale
  }
  return(x)
}


#' Validate rank and estimability of a default-interface design matrix
#'
#' Checks that a numeric design matrix supplied to the default matrix interface
#' has complete column names, contains no constant or non-informative columns,
#' and is full column rank after adding an intercept when appropriate.
#'
#' This is an internal defensive check for the matrix/default interface. Unlike
#' the formula interface, the default interface receives only a numeric design
#' matrix and cannot know whether columns came from raw variables, dummy coding,
#' ordered-factor polynomial contrasts, spline bases, or user-defined features.
#' Therefore this helper validates the actual estimability of the supplied
#' columns rather than trying to infer their origin from column names.
#'
#' For Cox models, `intercept` should be `FALSE` because Cox partial-likelihood
#' models do not include an intercept. For ordinary Gaussian, binomial, poisson,
#' and other GLM-style fits, `intercept` should usually be `TRUE`.
#'
#' @param x A numeric design matrix, or an object coercible to a matrix. Columns
#'   are candidate variables supplied to the default interface.
#' @param intercept Logical scalar. If `TRUE`, an intercept column is included
#'   when checking matrix rank. If `FALSE`, rank is checked on `x` alone.
#'
#' @return Invisibly returns `TRUE` if the design passes validation.
#'
#' @details
#' A column is treated as non-informative when, after removing non-finite values,
#' it has at most one unique value. Rank deficiency is detected using `qr()`.
#' If the design is rank-deficient, the helper reports columns identified by the
#' QR pivot as aliased or not estimable.
#'
#' This check is intentionally conservative. It is meant to fail early for input
#' matrices that cannot support valid model-comparison tests in the MFP
#' selection cycle. Step-level rank checks are still needed because a candidate
#' can become non-estimable in a particular adjustment model even if the initial
#' design passes this validation.
#'
#' @keywords internal
#' @noRd
validate_default_design_rank <- function(x, intercept = TRUE) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
    stop("`x` must have complete, non-empty column names.", call. = FALSE)
  }
  
  bad_constant <- vapply(
    seq_len(ncol(x)),
    function(j) {
      z <- x[, j]
      z <- z[is.finite(z)]
      length(unique(z)) <= 1L
    },
    logical(1L)
  )
  
  if (any(bad_constant)) {
    stop(
      sprintf(
        "The following columns in `x` are constant or non-informative: %s.",
        paste(colnames(x)[bad_constant], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  X <- if (intercept) cbind("(Intercept)" = 1, x) else x
  
  qr_x <- qr(X)
  full_rank <- qr_x$rank == ncol(X)
  
  if (!full_rank) {
    pivot <- qr_x$pivot
    # `!full_rank` guarantees at least one aliased column, so the bounds are
    # ordered. `seq.int()` makes that invariant explicit and avoids a fragile
    # programmatic colon expression.
    aliased_positions <- pivot[seq.int(qr_x$rank + 1L, ncol(X))]
    aliased_names <- colnames(X)[aliased_positions]
    aliased_names <- setdiff(aliased_names, "(Intercept)")
    
    stop(
      sprintf(
        paste0(
          "The supplied design matrix is rank-deficient. ",
          "The following column(s) are aliased or not estimable in the initial design: %s. ",
          "Remove redundant columns or provide a design in which each candidate variable ",
          "adds an estimable degree of freedom."
        ),
        paste(aliased_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}
#' Resolve Fitting Rows for Formula Interfaces
#'
#' Normalize a formula method's already validated `subset` argument to integer
#' row positions. The helper deliberately performs no expression evaluation.
#' Public methods validate the user-facing argument first; this helper also
#' rejects duplicated resolved positions defensively.
#'
#' @details
#' Formula methods need the selected row positions before they create the fitted
#' model frame. The full model frame remains available for documented full-data
#' preprocessing, while the returned positions identify the observations used
#' to rebuild factor levels, contrasts, responses, weights, offsets, and strata.
#'
#' A logical subset is converted with `which()`. Numeric positions are converted
#' to integer without reordering. `NULL` returns all row positions for common
#' sample-size checks and downstream alignment; callers must nevertheless bypass
#' model-frame and observation-vector subsetting when the public argument is
#' `NULL`. Callers must ensure that logical masks have length `nobs` and numeric
#' positions are otherwise valid. Numeric order is preserved, but duplicated
#' resolved positions are rejected here as an internal safety check.
#'
#' @param subset `NULL`, a logical mask, or numeric row positions that have
#'   already passed the formula method's public validation.
#' @param nobs Positive integer giving the number of rows in the full model
#'   frame.
#'
#' @return Integer vector containing the row positions used for fitting.
#' @keywords internal
#' @noRd
formula_subset_rows <- function(subset, nobs) {
  # Return concrete positions for shared size checks and alignment. Formula
  # methods still branch on the original public argument so they do not subset
  # frames or observation-level vectors when `subset` is NULL.
  if (is.null(subset)) {
    return(seq_len(nobs))
  }
  
  # Public validation has already ruled out NA values and length mismatches.
  rows <- if (is.logical(subset)) {
    which(subset)
  } else {
    # Preserve the user-supplied numeric ordering; do not sort here.
    as.integer(subset)
  }
  
  if (anyDuplicated(rows)) {
    stop(
      "! `subset` must not contain duplicated row indices.",
      call. = FALSE
    )
  }
  
  rows
}


#' Subset One Formula Factor Column
#'
#' Retain selected observations from a factor while preserving or safely
#' rebuilding its fitted contrast basis.
#'
#' @details
#' This helper is used only by the non-`NULL` formula-subset path. If all
#' original factor levels remain represented, the original level set and any
#' factor-specific `contrasts` attribute are preserved exactly. If levels are
#' removed and no factor-specific contrasts were supplied, unused levels are
#' dropped so a later `model.matrix()` call can regenerate the configured
#' default treatment or ordered-factor contrasts.
#'
#' A factor-specific contrast matrix or contrast specification generally cannot
#' be reduced uniquely after levels disappear. For predictor-side factor columns,
#' such a case is rejected rather than silently replacing the requested coding.
#' A factor response does not contribute contrast columns to the predictor design,
#' so its unused levels may be dropped without this predictor-side diagnostic.
#'
#' @param x Full-data factor column.
#' @param rows Unique retained row positions.
#' @param variable Display name used in diagnostics.
#' @param predictor Whether the factor contributes to the predictor-side formula
#'   model frame and therefore requires contrast reconstruction.
#'
#' @return The retained factor column.
#' @keywords internal
#' @noRd
subset_formula_factor_column <- function(x, rows, variable, predictor = TRUE) {
  if (!is.factor(x)) {
    stop("Internal error: `x` must be a factor.", call. = FALSE)
  }
  if (!is.character(variable) || length(variable) != 1L || is.na(variable)) {
    stop("Internal error: `variable` must be one non-missing character value.", call. = FALSE)
  }
  if (!is.logical(predictor) || length(predictor) != 1L || is.na(predictor)) {
    stop("Internal error: `predictor` must be TRUE or FALSE.", call. = FALSE)
  }
  
  retained <- x[as.integer(rows)]
  original_levels <- levels(x)
  observed_values <- unique(as.character(retained[!is.na(retained)]))
  retained_levels <- original_levels[original_levels %in% observed_values]
  levels_removed <- !identical(retained_levels, original_levels)
  factor_contrasts <- attr(x, "contrasts", exact = TRUE)
  
  if (!levels_removed) {
    # `[.factor` normally retains contrasts, but restore them explicitly so the
    # helper contract does not depend on that implementation detail.
    if (!is.null(factor_contrasts)) {
      attr(retained, "contrasts") <- factor_contrasts
    }
    return(retained)
  }
  
  if (predictor && !is.null(factor_contrasts)) {
    removed_levels <- original_levels[!original_levels %in% retained_levels]
    stop(
      paste0(
        "! Custom contrasts for factor `", variable,
        "` cannot be reconstructed safely after `subset` removed level(s): ",
        paste(removed_levels, collapse = ", "), ".\n",
        "i Use default contrasts, recode the retained-level factor before fitting, ",
        "or supply an explicit numeric design through the matrix interface."
      ),
      call. = FALSE
    )
  }
  
  droplevels(retained)
}


#' Build a Retained-Row Formula Model Frame
#'
#' Subset an already evaluated full-data model frame and remove unused factor
#' levels without re-evaluating the formula or the `subset` expression.
#'
#' @details
#' Formula methods first evaluate a complete model frame because automatic
#' preprocessing for continuous predictors is intentionally based on the full
#' data. They need a second frame only when the public `subset` argument is not
#' `NULL`. With no subset, callers reuse the complete frame directly and must not
#' call this helper.
#'
#' Calling `model.frame(..., subset = fit_rows)` from inside a formula method is
#' unsafe when `fit_rows` is a local variable. `model.frame()` captures the
#' expression and evaluates it in the formula/data environment, where that local
#' symbol is not necessarily visible. Directly subsetting the evaluated model
#' frame avoids that non-standard-evaluation boundary and also avoids evaluating
#' formula expressions, offsets, or strata terms a second time.
#'
#' Factor handling follows a preserve-or-reject policy. If all levels remain,
#' factor-specific contrasts are preserved exactly. If levels are removed from a
#' factor with no factor-specific contrasts, unused levels are dropped and a
#' subsequent `model.matrix()` call regenerates default treatment or ordered
#' polynomial contrasts. If a predictor factor loses levels while carrying a
#' custom contrast specification, the helper stops because no unique reduced
#' contrast basis can be inferred safely.
#'
#' The `terms` attribute is restored explicitly because it is part of the model
#' frame contract used by downstream formula parsing and model-matrix generation.
#' No model-matrix columns are copied or subset by this helper.
#'
#' @param model_frame Full-data model frame returned by `stats::model.frame()`.
#' @param rows Integer row positions selected by a non-`NULL` subset. The caller
#'   must have validated their range, uniqueness, and minimum sample size.
#'
#' @return A model frame containing only `rows`, with unused factor levels
#'   removed and the original `terms` attribute retained.
#' @keywords internal
#' @noRd
subset_formula_model_frame <- function(model_frame, rows) {
  if (!is.data.frame(model_frame)) {
    stop("Internal error: `model_frame` must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) ||
      any(rows != as.integer(rows)) || any(rows < 1L) ||
      any(rows > nrow(model_frame)) || anyDuplicated(rows)) {
    stop(
      "Internal error: `rows` must contain unique valid integer model-frame positions.",
      call. = FALSE
    )
  }
  
  # Preserve the formula contract before ordinary data-frame subsetting and
  # droplevels() manipulate row names and factor columns.
  terms_object <- attr(model_frame, "terms", exact = TRUE)
  
  retained <- model_frame[as.integer(rows), , drop = FALSE]
  
  factor_columns <- which(vapply(model_frame, is.factor, logical(1L)))
  response_position <- if (is.null(terms_object)) {
    0L
  } else {
    response_value <- attr(terms_object, "response")
    if (is.null(response_value)) 0L else response_value
  }
  
  for (column_position in factor_columns) {
    column_name <- names(model_frame)[column_position]
    retained[[column_position]] <- subset_formula_factor_column(
      x = model_frame[[column_position]],
      rows = rows,
      variable = column_name,
      predictor = column_position != response_position
    )
  }
  
  if (!is.null(terms_object)) {
    attr(retained, "terms") <- terms_object
  }
  
  retained
}


#' Align the Full-Data Preprocessing Design to the Fitting Design
#'
#' Create a full-row preprocessing matrix whose columns exactly match the
#' subset-specific fitting matrix. This is the bridge between the two data
#' sources used by formula methods:
#'
#' * the full data, which retain the package's documented automatic continuous
#'   preprocessing source; and
#' * the retained observations, which determine factor levels, contrasts,
#'   fitting columns, model selection, and prediction metadata.
#'
#' @details
#' Subsetting can change a factor's generated columns. For example, reducing a
#' four-level ordered factor to three levels changes `contr.poly(4)` to
#' `contr.poly(3)`. Consequently, `full_x` and `x` can have different
#' categorical columns even though their continuous columns retain stable names.
#'
#' The returned matrix has `nrow(full_x)` rows and the columns of `x`, in the
#' same order. Same-named columns are copied from `full_x`; columns that exist
#' only in the fitted design are filled with zero. The zero fill is an internal
#' placeholder only. Such columns belong to mapped categorical blocks whose
#' preprocessing is fixed (`df = 1`, shift 0, scale 1, and no continuous-only
#' extensions), so their full-row placeholder values must never affect fitting
#' or automatic continuous preprocessing.
#'
#' This helper does not transfer factor levels, contrasts, or other prediction
#' metadata. Those must always come from the subset-specific model frame and
#' fitting design.
#'
#' @param x Numeric subset-specific fitting design matrix after public column
#'   renaming and intercept removal.
#' @param full_x Numeric full-data design matrix after the same public column
#'   renaming and intercept removal.
#'
#' @return Numeric matrix with the rows of `full_x` and the columns of `x`.
#' @keywords internal
#' @noRd
align_formula_preprocess_matrix <- function(x, full_x) {
  if (!is.matrix(x) || !is.matrix(full_x)) {
    stop("Internal error: formula preprocessing designs must be matrices.", call. = FALSE)
  }
  if (is.null(colnames(x)) || is.null(colnames(full_x))) {
    stop("Internal error: formula preprocessing designs must have column names.", call. = FALSE)
  }
  
  # Start from a deterministic placeholder matrix. Only same-named columns are
  # eligible to carry full-data values into automatic preprocessing.
  aligned <- matrix(
    0,
    nrow = nrow(full_x),
    ncol = ncol(x),
    dimnames = list(rownames(full_x), colnames(x))
  )
  
  # Continuous model-matrix columns normally appear in both designs with stable
  # names. Factor columns may also share names, but their values are harmless
  # because mapped categorical settings prevent their automatic preprocessing.
  common <- intersect(colnames(x), colnames(full_x))
  if (length(common) > 0L) {
    aligned[, common] <- full_x[, common, drop = FALSE]
  }
  
  aligned
}


#' Attach the Full-Data Preprocessing Design to a Fitting Matrix
#'
#' Store the aligned full-data preprocessing source as a private matrix
#' attribute before a formula method delegates to the default method.
#'
#' @details
#' This is an internal transport mechanism, not fitted-object metadata and not a
#' public API. `extract_preprocess_matrix()` removes the attribute immediately
#' on entry to the default method so it cannot leak into later matrix operations,
#' fitting objects, or user-visible output.
#'
#' Direct matrix-interface calls do not carry this attribute and therefore use
#' their supplied matrix for both preprocessing and fitting, preserving existing
#' matrix behavior.
#'
#' @inheritParams align_formula_preprocess_matrix
#'
#' @return `x` with a private `mfp2_preprocess_x` attribute containing the
#'   aligned full-data preprocessing matrix.
#' @keywords internal
#' @noRd
attach_formula_preprocess_matrix <- function(x, full_x) {
  # Align before attaching so the default method can require exact column-name
  # compatibility and avoid positional matching.
  attr(x, "mfp2_preprocess_x") <- align_formula_preprocess_matrix(x, full_x)
  x
}


#' Extract and Validate the Preprocessing Source
#'
#' Separate the matrix used for fitting from the matrix used to estimate
#' automatic preprocessing quantities.
#'
#' @details
#' For direct matrix calls, no private attribute is present and the supplied
#' matrix is returned for both roles. For formula calls, the fitting matrix is
#' subset-specific, while the attached preprocessing matrix contains full-data
#' values aligned to the fitted columns.
#'
#' The private attribute is removed before returning. The preprocessing matrix
#' must be numeric, finite, and have exactly the same column-name set as the
#' fitting matrix; it is reordered to the fitting-column order. Row counts may
#' differ by design because preprocessing can use all observations while fitting
#' uses only the retained subset.
#'
#' @param x Numeric fitting matrix, optionally carrying the private
#'   `mfp2_preprocess_x` attribute set by a formula method.
#'
#' @return Named list with:
#'   \describe{
#'     \item{`x`}{The fitting matrix with the private attribute removed.}
#'     \item{`preprocess_x`}{The full-data preprocessing source, or `x` for a
#'       direct matrix call.}
#'   }
#' @keywords internal
#' @noRd
extract_preprocess_matrix <- function(x) {
  preprocess_x <- attr(x, "mfp2_preprocess_x", exact = TRUE)
  
  # Remove the transport attribute immediately. Downstream code should work
  # with two explicit matrices rather than relying on hidden matrix attributes.
  attr(x, "mfp2_preprocess_x") <- NULL
  
  if (is.null(preprocess_x)) {
    return(list(x = x, preprocess_x = x))
  }
  if (!is.matrix(preprocess_x) || !is.numeric(preprocess_x)) {
    stop("Internal error: the preprocessing source must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(colnames(preprocess_x)) ||
      !setequal(colnames(preprocess_x), colnames(x))) {
    stop(
      "Internal error: preprocessing and fitting matrices must have identical column names.",
      call. = FALSE
    )
  }
  
  # Name-based reordering is required because formula processing can change the
  # physical column order while preserving conceptual column identities.
  preprocess_x <- preprocess_x[, colnames(x), drop = FALSE]
  if (anyNA(preprocess_x) || any(!is.finite(preprocess_x))) {
    stop(
      "Internal error: the full-data preprocessing matrix contains non-finite values.",
      call. = FALSE
    )
  }
  
  list(x = x, preprocess_x = preprocess_x)
}


#' Validate Grouped Matrix Terms After Subsetting
#'
#' Ensure that every explicitly mapped matrix block retains the estimable
#' dimension it had before rows were removed.
#'
#' @details
#' The matrix interface receives numeric columns but no factor levels, contrast
#' function, or model-frame recipe. It therefore cannot regenerate a categorical
#' basis when `subset` removes a represented level. The safe behavior is to keep
#' the supplied columns unchanged and reject a mapped block whose within-block
#' rank decreases.
#'
#' Rank is measured after adjoining a constant column and subtracting that
#' constant's contribution. This measures the block's estimable dimension in an
#' intercept-based model and detects, for example, an all-zero treatment dummy,
#' a binary indicator that becomes constant, or retained rows of `contr.poly(k)`
#' after one of the original levels disappears.
#'
#' Only explicitly mapped terms are checked. Cross-term collinearity and the
#' rank of the complete fitted design remain the responsibility of the existing
#' global design validation. The helper does not modify or drop columns.
#'
#' @param x_full Numeric design matrix before applying `subset`.
#' @param x_fit Numeric design matrix containing only the retained rows.
#' @param term_to_columns Complete conceptual-term mapping from term names to
#'   raw design-column names.
#'
#' @return Invisibly returns `TRUE`. Raises an error naming all mapped terms
#'   whose estimable dimension decreases.
#' @keywords internal
#' @noRd
validate_grouped_subset_rank <- function(x_full, x_fit, term_to_columns) {
  mapped_terms <- names(term_to_columns)[mapped_term_flags(term_to_columns)]
  if (length(mapped_terms) == 0L) {
    return(invisible(TRUE))
  }
  
  block_rank <- function(z) {
    if (ncol(z) == 0L || nrow(z) == 0L) return(0L)
    
    # Subtract the rank contributed by the constant so a one-column 0/1 dummy
    # contributes one degree of freedom only when both values remain observed.
    qr(cbind("(Constant)" = 1, z))$rank - 1L
  }
  
  lost <- mapped_terms[vapply(mapped_terms, function(term) {
    cols <- term_to_columns[[term]]
    block_rank(x_fit[, cols, drop = FALSE]) <
      block_rank(x_full[, cols, drop = FALSE])
  }, logical(1L))]
  
  if (length(lost) > 0L) {
    stop(
      paste0(
        "The following grouped term(s) lose estimable dimension after applying `subset`: ",
        paste(lost, collapse = ", "), ". The matrix interface cannot rebuild ",
        "factor levels or contrasts from numeric columns. Construct the model ",
        "matrix from the intended fitting rows, or use the formula interface."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Validate Predictor Variation After Subsetting
#'
#' Check that selected numeric predictor columns still contain at least two
#' distinct values after the matrix interface applies `subset`.
#'
#' @details
#' This check is complementary to `validate_grouped_subset_rank()`. It reports
#' ordinary singleton columns that become constant, while grouped blocks are
#' validated by their joint estimable dimension. Input matrices are expected to
#' have passed the package's earlier finite-value checks, so missing and infinite
#' values are not handled here.
#'
#' `exclude` is used by MFPI for `group_var`, whose post-subset group-count and
#' group-size requirements are checked separately with more informative errors.
#'
#' @param x Numeric fitting matrix after subsetting.
#' @param exclude Optional character vector of column names omitted from this
#'   check because they have dedicated validation elsewhere.
#'
#' @return Invisibly returns `TRUE`. Raises an error listing constant columns.
#' @keywords internal
#' @noRd
validate_subset_predictor_variation <- function(x, exclude = NULL) {
  columns <- setdiff(colnames(x), exclude)
  if (length(columns) == 0L) return(invisible(TRUE))
  
  bad <- vapply(columns, function(column) {
    length(unique(x[, column])) <= 1L
  }, logical(1L))
  
  if (any(bad)) {
    stop(
      paste0(
        "The following predictor column(s) have no variation after applying `subset`: ",
        paste(columns[bad], collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Validate Formula Factor Levels in the Fitting Frame
#'
#' Confirm that every factor- or character-valued predictor represented in the
#' formula has at least two observed levels after `subset` and unused-level
#' dropping.
#'
#' @details
#' Formula methods rebuild their model matrix from the retained observations.
#' A categorical predictor with one remaining level has no estimable effect and
#' cannot receive valid treatment or polynomial contrasts. Detecting it before
#' `model.matrix()` yields a term-specific error rather than a lower-level
#' contrast or rank failure.
#'
#' Predictor expressions are obtained from the terms incidence matrix and then
#' matched to evaluated columns in `model_frame`. This includes ordinary factor
#' variables and evaluated inline expressions such as `factor(stage)` or
#' `ordered(stage)` when they appear as model-frame columns. The response is not
#' checked.
#'
#' @param terms_object Predictor terms object after removing non-predictor
#'   specials such as `strata()` and `offset()`.
#' @param model_frame Subset-specific model frame with unused levels already
#'   dropped.
#'
#' @return Invisibly returns `TRUE`. Raises an error listing categorical
#'   predictors with fewer than two fitted levels.
#' @keywords internal
#' @noRd
validate_formula_factor_levels <- function(terms_object, model_frame) {
  incidence <- attr(terms_object, "factors")
  if (is.null(incidence) || nrow(incidence) == 0L) {
    return(invisible(TRUE))
  }
  
  # Restrict the check to expressions that participate in predictor terms.
  variables <- rownames(incidence)[rowSums(incidence != 0) > 0L]
  factor_variables <- variables[vapply(variables, function(variable) {
    variable %in% names(model_frame) &&
      (is.factor(model_frame[[variable]]) || is.character(model_frame[[variable]]))
  }, logical(1L))]
  
  bad <- factor_variables[vapply(factor_variables, function(variable) {
    value <- model_frame[[variable]]
    observed_levels <- if (is.factor(value)) {
      nlevels(value)
    } else {
      length(unique(value))
    }
    observed_levels < 2L
  }, logical(1L))]
  
  if (length(bad) > 0L) {
    stop(
      paste0(
        "The following categorical variable(s) have fewer than two observed levels after applying `subset`: ",
        paste(bad, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}
