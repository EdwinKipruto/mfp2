#' Select a fixed-knot cubic B-spline form
#'
#' Candidate coefficient counts correspond to centered, nested subsets of
#' the same interior knots. Their basis columns were assembled into prefixes
#' before selection, so the existing grouped fitter can score each candidate.
#' @param args Arguments prepared by `find_best_fp_step()`.
#' @param spec Fitted B-spline specification with nested candidate counts.
#' @return Selected spline df, model metrics, and adjustment-cache metadata.
#' @keywords internal
#' @noRd
select_bspline_step <- function(args, spec) {
  select_spline_step(args, spec$df, c(1L, spec$candidate_df))
}
