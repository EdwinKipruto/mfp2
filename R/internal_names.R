# Internal formula-column names ------------------------------------------------

#' Allocate a Collision-Safe Internal Column Name
#'
#' Formula-based final refits need temporary data-frame columns for package-
#' created quantities such as the response, offset, and Cox strata. User
#' predictor names are left unchanged; only these temporary names are allocated
#' here. The helper is called once per final refit, not during repeated MFP
#' candidate fitting, so it has no meaningful effect on search performance.
#'
#' @param base Character scalar used to form the internal name.
#' @param used Character vector of names that must not be reused.
#' @param preferred Optional preferred helper name. It is used when available;
#'   otherwise a collision-safe `..mfp2_` name is allocated.
#'
#' @return A collision-free syntactic character scalar.
#'
#' @keywords internal
#' @noRd
mfp2_internal_name <- function(base, used = character(), preferred = NULL) {
  if (!is.character(base) || length(base) != 1L || is.na(base) || !nzchar(base)) {
    stop("Internal error: `base` must be a non-empty character scalar.", call. = FALSE)
  }

  used <- as.character(used)

  if (!is.null(preferred)) {
    if (!is.character(preferred) || length(preferred) != 1L ||
        is.na(preferred) || !nzchar(preferred)) {
      stop(
        "Internal error: `preferred` must be NULL or a non-empty character scalar.",
        call. = FALSE
      )
    }
    if (!preferred %in% used) {
      return(preferred)
    }
  }

  stem <- paste0("..mfp2_", base)
  candidate <- stem
  suffix <- 0L

  while (candidate %in% used) {
    suffix <- suffix + 1L
    candidate <- paste0(stem, "_", suffix)
  }

  candidate
}


#' Recover a Stored Internal Formula-Column Name
#'
#' Formula-based final fits store the exact collision-safe names embedded in
#' their model formulas. Prediction must reuse those names rather than guessing
#' from historical conventions or parsing the formula. Fitted objects that lack
#' this metadata are deliberately rejected and must be refitted.
#'
#' @param fit_obj Fitted model object.
#' @param component Stored scalar component to retrieve, such as `"offset"` or
#'   `"strata"`.
#'
#' @return The stored non-empty character scalar.
#'
#' @keywords internal
#' @noRd
mfp2_internal_fit_name <- function(fit_obj, component) {
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    stop(
      "Internal error: `component` must be a non-empty character scalar.",
      call. = FALSE
    )
  }

  internal <- fit_obj$mfp2_internal_names
  if (!is.list(internal) || !component %in% names(internal)) {
    stop(
      "The fitted model lacks internal-name metadata required for prediction. ",
      "Refit the model with the current package version.",
      call. = FALSE
    )
  }

  value <- internal[[component]]
  if (!is.character(value) || length(value) != 1L ||
      is.na(value) || !nzchar(value)) {
    stop(
      "The fitted model contains invalid internal-name metadata for `",
      component, "`. Refit the model with the current package version.",
      call. = FALSE
    )
  }

  value
}
