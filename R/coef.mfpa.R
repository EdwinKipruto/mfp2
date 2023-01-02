#' Extract coefficients from a `mfpa` object 
#' 
#' This function is a method for the generic [stats::coef()] function for objects 
#' of class `mfpa`. 
#' 
#' @param object an object of class `mfpa`, usually, a result of a call to
#' [mfpa()].
#' 
#' @return 
#' Named numeric vector of coefficients extracted from the model `object`.
#' 
#' @export
coef.mfpa <- function(x) {
  x$coefficients
}
