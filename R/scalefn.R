#' A function that calculates an integer used to scale x. Estimation of powers
#' is unaffected by scaling, the same powers are found for x/scale. In extreme
#' cases scaling is necessary to preserve accuracy (Royston and Sauerbrei 2008).
#' We used the formula in section 4.11.1 of Royston and Sauerbrei (2008) book
#'
#'@param x a vector of predictor variable where x is already shifted to positive values.
#'@return An integer
#'@export
scalefn <- function(x){
  # see https://www.stata.com/manuals/rfp.pdf
  mx <- unique(x)
  if(length(mx)==1) stop("all values of x are identical hence log(max(x)-min(x)) = log(0) not allowed")
  lrange <- log10(max(x, na.rm = T)-min(x, na.rm = T))
  #scales <- 10^(sign(lrange)*as.integer(lrange))
  scales <- 10^(sign(lrange)*floor(abs(lrange)))
  return(scales)
}

