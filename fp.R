#' Helper function to assign attribute to a variable
#' 
#' @return 
#' The variable with its attributes.
#' 
#' @export
fp <- function(x, df = 4, select = NA, alpha = NA, scale=TRUE, acd = FALSE)
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
ut=data.frame(y=1:10)
kk=fp(ut)
attributes(kk)
