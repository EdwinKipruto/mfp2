#' Resolve Interaction Forms from Variable Types
#'
#' Continuous variables default to FP1. Binary and categorical variables use
#' linear interactions, with a warning when an explicit FP form is overridden.
#'
#' @param interaction_forms NULL or a named character vector of requested forms.
#' @param interaction_specs Named list describing the tested interaction terms.
#'
#' @return A named character vector in interaction-term order.
#' @keywords internal
#' @noRd
resolve_mfpi_interaction_forms <- function(interaction_forms,
                                           interaction_specs) {
  interaction_vars <- names(interaction_specs)
  discrete <- vapply(interaction_specs, function(spec) {
    isTRUE(spec$discrete)
  }, logical(1L))
  resolved <- stats::setNames(
    ifelse(discrete, "linear", "fp1"), interaction_vars
  )

  if (is.null(interaction_forms)) return(resolved)

  if (!is.character(interaction_forms)) {
    stop(
      "! `interaction_forms` must be a named character vector or NULL.",
      call. = FALSE
    )
  }
  form_names <- names(interaction_forms)
  if (is.null(form_names) || anyNA(form_names) || any(!nzchar(form_names))) {
    stop(
      paste0(
        "! Every entry of `interaction_forms` must be named with a variable ",
        "from `interaction_vars`. Unnamed entries are not allowed.\n",
        "  Example: interaction_forms = c(age = \"fp2\", bmi = \"fp1\")"
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(form_names)) {
    stop("! Names in `interaction_forms` must be unique.", call. = FALSE)
  }

  unknown_names <- setdiff(form_names, interaction_vars)
  if (length(unknown_names) > 0L) {
    stop(
      paste0(
        "! The following name(s) in `interaction_forms` are not in `interaction_vars`: ",
        paste(unknown_names, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  valid_forms <- c("linear", "fp1", "fp2")
  bad_values <- interaction_forms[!interaction_forms %in% valid_forms]
  if (length(bad_values) > 0L) {
    stop(
      paste0(
        "! Invalid value(s) in `interaction_forms`: ",
        paste(unique(bad_values), collapse = ", "), ".\n",
        "  Allowed values are: \"linear\", \"fp1\", \"fp2\"."
      ),
      call. = FALSE
    )
  }

  resolved[form_names] <- interaction_forms
  overridden <- form_names[
    discrete[form_names] & interaction_forms != "linear"
  ]
  if (length(overridden) > 0L) {
    warning(
      paste0(
        "Interaction forms for binary or categorical variables have been reset ",
        "to 'linear': ", paste(overridden, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  resolved[discrete] <- "linear"

  resolved
}
