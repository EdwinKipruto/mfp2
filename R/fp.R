#' A function that assign attribute to a continuous variable that undergo fp transformation
#' 
#' @param x a vector of continuous variable to be assigned attributes

#' @param df a numeric value specifying the degrees of freedom (df) for
#'  predictor x. Default is 4 (FP2)
#' @param select a numeric value that sets the nominal significance levels for 
#' variable selection by backward elimination. The default nominal significance
#' level is 0.05. Setting the nominal significance level to be 1 forces the into
#' the model.
#' @param alpha a numeric value that sets the significance levels for 
#' function selection. The default nominal significance level is 0.05. 
#' @param shift a numeric value specifying a shifting value for x. Default is
#'  `NULL`which lets the program estimate the shifting factor. If shifting is 
#'  not required, set `shift = 0` to disable it.
#' @param scale a numeric value specifying  scaling factors for x. Default 
#' is `NULL` which lets the program estimate the scaling factors.
#' If scaling is not required set `scale = 1` to disable it.
#' @param center a logical determining whether transformed variable x 
#' should be centered before fitting the final model. The default `TRUE` 
#' which implies mean centering, except for binary covariates, where the 
#' covariate is centered using the lower of the two distinct values of the
#' covariate. 
#' @param acd logical indicating whether variable x should undergo 
#' the approximate cumulative distribution (ACD) transformation.
#' It also invokes the function-selection procedure to determine the 
#' best-fitting FP1(p1, p2) model. 
#' @param pow a distinct vector of powers to be evaluated for covariate x. Default is NULL
#' and powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3) will be used. 
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

