#' Select a restricted cubic spline term
#' @param args Arguments prepared by `find_best_fp_step()`.
#' @param max_df Maximum number of spline columns.
#' @return Selected form and candidate fit metrics.
#' @keywords internal
#' @noRd
select_rcs_step <- function(args, max_df) {
  select_spline_step(args, max_df, seq_len(max_df))
}
