#' Evaluate a fitted spline specification at original-scale values
#' @param x Numeric predictor values.
#' @param spec Fitted RCS or cubic B-spline specification.
#' @return A complete spline matrix including the optional zero indicator.
#' @keywords internal
#' @noRd
spline_basis_from_spec <- function(x, spec) {
  if (identical(spec$type, "bspline")) {
    return(bspline_basis(x, spec$knots, spec$boundary_knots,
                         indicator = spec$indicator, terms = spec$basis_terms))
  }
  rcs_basis(x, spec$knots, indicator = spec$indicator)
}

#' Describe either selected spline form
#' @param df Selected coefficient count, or `NA` for exclusion.
#' @param indicator Whether the zero indicator is present.
#' @param spec Fitted spline specification.
#' @return A human-readable model label.
#' @keywords internal
#' @noRd
spline_form_label <- function(df, indicator, spec) {
  if (identical(spec$type, "bspline")) {
    bspline_form_label(df, indicator)
  } else {
    rcs_form_label(df, indicator)
  }
}

#' Describe a fitted spline basis column
#' @param term Original predictor name.
#' @param position One-based position in the full spline design.
#' @param spec Fitted spline specification.
#' @return A human-readable coefficient label.
#' @keywords internal
#' @noRd
spline_column_label <- function(term, position, spec) {
  if (identical(spec$type, "bspline")) {
    bspline_column_label(term, position, spec)
  } else {
    rcs_column_label(term, position, spec$knots)
  }
}

#' Retrieve spline metadata, including from older RCS-only fits
#' @param object A fitted `mfp2` object.
#' @param term Original predictor name.
#' @return Fitted term specification or `NULL`.
#' @keywords internal
#' @noRd
prediction_spline_spec <- function(object, term) {
  spline_specs_from_fit(object)[[term]]
}

#' Retrieve selected spline complexity
#' @param object A fitted `mfp2` object.
#' @param term Original predictor name.
#' @return Selected coefficient count, or `NA` when excluded.
#' @keywords internal
#' @noRd
prediction_spline_selected_df <- function(object, term) {
  spline_selected_df_from_fit(object)[[term]]
}

#' Retrieve all spline specifications from current and older fits
#' @param object A fitted `mfp2` object.
#' @return Named list of spline specifications.
#' @keywords internal
#' @noRd
spline_specs_from_fit <- function(object) {
  if (!is.null(object$spline_specs)) object$spline_specs else object$rcs_specs
}

#' Retrieve all selected spline coefficient counts
#' @param object A fitted `mfp2` object.
#' @return Named integer vector; missing values indicate exclusion.
#' @keywords internal
#' @noRd
spline_selected_df_from_fit <- function(object) {
  if (!is.null(object$spline_selected_df)) {
    object$spline_selected_df
  } else object$rcs_selected_df
}

#' Label the selected continuous spline form for a plot title
#'
#' Selection can reduce a spline to a linear term or a zero indicator. In
#' those cases, a title naming the requested maximum spline would misstate
#' what was fitted.
#' @param spec Fitted spline specification, including its type.
#' @param selected_df Number of continuous basis columns retained.
#' @return A concise label for the selected continuous form.
#' @keywords internal
#' @noRd
spline_plot_form <- function(spec, selected_df) {
  if (selected_df == 0L) return("zero indicator")
  if (selected_df == 1L) return("linear")
  form <- if (identical(spec$type, "bspline")) "bs" else "rcs"
  sprintf("%s(df=%d)", form, selected_df)
}
