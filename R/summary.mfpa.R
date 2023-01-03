#' Summarizing `mfpa` model fits
#' 
#' This function is a method for the generic [base::summary()] function for
#' objects of class `mfpa`.
#' 
#' @param object an object of class `mfpa`, usually, a result of a call to
#' [mfpa()].
#' @param ... further arguments passed to the summary functions for `glm()` 
#' ([stats::summary.glm()], i.e. families supported by `glm()`) or `coxph()` 
#' ([survival::summary.coxph()], if `object$family = "cox"`).
#' 
#' @return 
#' An object returned from [stats::summary.glm()] or
#' [survival::summary.coxph()], depending on the family parameter of `object`.
#' 
#' @seealso 
#' [mfpa()], [stats::glm()], [stats::summary.glm()], [survival::coxph()],
#' [survival::summary.coxph()]
#' 
#' @export
summary.mfpa <- function(object, ...) {
  NextMethod("summary", object)
}
