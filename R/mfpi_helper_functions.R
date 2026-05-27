# Helpers for extracting information from mfp2 model objects
#
# Both functions are exported because they are useful to callers who build
# custom workflows on top of mfpi() and need to inspect the adjustment model
# selected by the MFP algorithm.
#
# get_selected_variables() -- returns names of variables retained by MFP
# get_fp_powers()          -- returns the FP power vector(s) for named variables


# -----------------------------------------------------------------------------
# get_selected_variables() ----------------------------------------------------
# -----------------------------------------------------------------------------

#' Extract Selected Variable Names from an MFP Model
#'
#' Returns the names of variables that were retained by the MFP algorithm,
#' as recorded in the \code{fp_terms} component of an \code{"mfp2"} object.
#' Names refer to the original variable names before any FP transformation.
#' If no variables were selected the function returns an empty character vector.
#'
#' @param object An object of class \code{"mfp2"}, as returned by
#'   \code{mfp2::mfp2()} or \code{mfp2:::fit_mfp()}.
#'
#' @return A character vector of variable names (row names of
#'   \code{object$fp_terms}) for which the \code{selected} column is
#'   \code{TRUE}. Returns \code{character(0)} if no variables were selected.
#'
#' @examples
#' \dontrun{
#' fit          <- mfp2::mfp2(y ~ fp(x1) + fp(x2), data = mydata)
#' selected     <- get_selected_variables(fit)
#' }
#'
#' @export
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
  
  if (!is.matrix(fp_terms) && !is.data.frame(fp_terms)) {
    stop("! `fp_terms` must be a matrix or data frame.", call. = FALSE)
  }
  
  # If the 'selected' column is absent, no variable survived selection
  if (!"selected" %in% colnames(fp_terms)) {
    return(character(0L))
  }
  
  selected <- fp_terms[, "selected"]
  
  if (!is.logical(selected)) {
    stop("! The 'selected' column in `fp_terms` must be logical.", call. = FALSE)
  }
  
  rownames(fp_terms)[selected]
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
#' @export
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
