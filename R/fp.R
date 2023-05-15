#' A function that assign attribute to a continuous variable that undergo fp transformation
#' 
#' @param x a vector of continuous variable to be assigned attributes
#' @param df,alpha,select,shift,scale,center,acd See [mfp2::mfp2()]) for details. 
#' @param pow a vector of powers to be evaluated for x. Default is NULL and 
#' powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3) will be used
#' @return 
#' a variable with its attributes.
#' 
#' @export
fp <- function(x, df = 4, alpha = 0.05,select = 0.05, shift = NULL, scale=NULL,center = TRUE, acd = FALSE, pow = NULL)
{

  name <- deparse(substitute(x))
  # Assert that a factor variable must not be subjected to fp transformation
  if(is.factor(x))
    stop(name," is a factor variable and should not be passed to the fp() function.")
  # attributes of x
  attr(x, "df") <- df
  attr(x, "alpha") <- alpha
  attr(x, "select") <- select
  attr(x, "shift") <- shift
  attr(x, "scale") <- scale
  attr(x, "center") <- center
  attr(x, "acd") <- acd
  attr(x, "pow") <- pow
  attr(x, "name") <- name
  x
}

