#' Helper to assign attributes to a variable undergoing FP-transformation
#' 
#' Used in formula interface to `mfp2()`.
#' 
#' @param x a vector representing a continuous variable undergoing 
#' fp-transformation.
#' @param df,alpha,select,shift,scale,center,acd See [mfp2::mfp2()]) for details. 
#' @param powers a vector of powers to be evaluated for `x`. Default is `NULL` 
#' and `powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)` will be used.
#' 
#' @return 
#' The vector `x` with new attributes relevant for fp-transformation. All 
#' arguments passed to this function will be stored as attributes. 
#' 
#' @export
fp <- function(x, 
               df = 4, 
               alpha = 0.05,
               select = 0.05, 
               shift = NULL, 
               scale=NULL,
               center = TRUE, 
               acd = FALSE, 
               powers = NULL) {

  name <- deparse(substitute(x))
  
  # Assert that a factor variable must not be subjected to fp transformation
  if(is.factor(x))
    stop(name," is a factor variable and should not be passed to the fp() function.")
  
  attr(x, "df") <- df
  attr(x, "alpha") <- alpha
  attr(x, "select") <- select
  attr(x, "shift") <- shift
  attr(x, "scale") <- scale
  attr(x, "center") <- center
  attr(x, "acd") <- acd
  attr(x, "powers") <- powers
  attr(x, "name") <- name
  
  x
}
