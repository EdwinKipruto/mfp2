# Run-time checks for packages listed in Suggests.
#
# ggplot2 (plotting) and geepack (GEE models) are optional. Every code path
# that needs one of them calls the matching helper first so that the user gets
# an informative error instead of a failure deep inside a computation.

#' Is a suggested package available?
#'
#' Thin wrapper around [requireNamespace()] so that the availability checks can
#' be mocked in tests. Loading the namespace also registers its S3 methods.
#'
#' @param pkg A single package name.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @keywords internal
#' @noRd
mfp2_has_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Stop unless ggplot2 is available
#'
#' @return Invisibly `TRUE`; otherwise an error.
#'
#' @keywords internal
#' @noRd
mfp2_require_ggplot2 <- function() {
  if (!mfp2_has_package("ggplot2")) {
    stop(
      "! Plotting requires the 'ggplot2' package.\n",
      "i Install it with install.packages(\"ggplot2\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless geepack is available
#'
#' @return Invisibly `TRUE`; otherwise an error.
#'
#' @keywords internal
#' @noRd
mfp2_require_geepack <- function() {
  if (!mfp2_has_package("geepack")) {
    stop(
      "! GEE models require the 'geepack' package.\n",
      "i Install it with install.packages(\"geepack\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless geepack is available for a fitted GEE model
#'
#' Methods for a fitted GEE model (for example `predict()`, `summary()` or
#' `residuals()`) dispatch to geepack's `geeglm` methods. A model restored with
#' [readRDS()] in a new session does not load geepack, so check, and load its
#' namespace, before delegating.
#'
#' @param object A fitted `mfp2` object.
#'
#' @return Invisibly `TRUE`; otherwise an error.
#'
#' @keywords internal
#' @noRd
mfp2_require_geepack_for_fit <- function(object) {
  if (identical(object$family_string, "gee")) {
    mfp2_require_geepack()
  }
  invisible(TRUE)
}
