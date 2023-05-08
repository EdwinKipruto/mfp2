#' Helper function to assign attribute to a variable
#' 
#' @return 
#' a variable with its attributes.
#' 
#' @export
fp <- function(x, df = 4, select = 0.05, alpha = 0.05, scale=TRUE, acd = FALSE)
{
  #
  name <- deparse(substitute(x))
  attr(x, "df") <- df
  attr(x, "alpha") <- alpha
  attr(x, "select") <- select
  attr(x, "scale") <- scale
  attr(x, "acd") <- acd
  attr(x, "name") <- name
  x
}

