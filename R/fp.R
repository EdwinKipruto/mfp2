#' A function that assign attribute to a variable that undergo fp transformation
#' 
#' @param x a vector of variable to be assigned attributes
#' @param df,alpha,select,shift,scale,center,acd See [mfp2::mfp2()]) for details. 
#' @return 
#' a variable with its attributes.
#' 
#' @export
fp <- function(x, df = 4, alpha = 0.05,select = 0.05, shift = NULL, scale=NULL,center = TRUE, acd = FALSE)
{
  name <- deparse(substitute(x))
  # attributes of x
  attr(x, "df") <- df
  attr(x, "alpha") <- alpha
  attr(x, "select") <- select
  attr(x, "shift") <- shift
  attr(x, "scale") <- scale
  attr(x, "center") <- center
  attr(x, "acd") <- acd
  attr(x, "name") <- name
  x
}

