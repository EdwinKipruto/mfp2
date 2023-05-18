#' A function that assign attribute to a continuous variable that undergo fp transformation
#' 
#' @param x a vector of continuous variable to be assigned attributes
#' @param df,alpha,select,shift,scale,center,acd See [mfp2::mfp2()]) for details. 
#' @param scale a numeric vector of length nvars or single numeric specifying 
#' scaling factors. If a single numeric, then the value will be replicated as
#' necessary. Default is `NULL` which lets the program estimate the scaling 
#' factors (see Details section).
#' If scaling is not required set `scale = 1` to disable it.
#' @param shift a numeric vector of length nvars or a single numeric specifying
#' shift terms. If a single numeric, then the value will be replicated as
#' necessary. Default is `NULL` which lets the program estimate the shifts
#' (see Details section).
#' If shifting is not required, set `shift = 0` to disable it.
#' @param df a numeric vector of length nvars or a single numeric that sets the 
#' (default) degrees of freedom (df) for each predictor. If a single numeric, 
#' then the value will be replicated as necessary. The df (not counting 
#' the intercept) are twice the degree of a fractional polynomial (FP). 
#' For example, an FP2 has 4 df, while FP3 has 6 df. 
#' The program overrides default df based on the number of distinct (unique) 
#' values for a variable as follows: 
#' 2-3 distinct values are assigned `df = 1` (linear), 4-5 distinct values are
#' assigned `df = min(2, default)` and >= 6 distinct values are assigned  
#' `df = default`. NOT in mfp2.formula()...df = 1 makes sense in formula to avoid warnings from binary variables
#' @param center a logical determining whether variables are centered before 
#' model fit. The default `TRUE` implies mean centering, except for binary 
#' covariates, where the covariate is centered using the lower of the two 
#' distinct values of the covariate. See Details section below.
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

