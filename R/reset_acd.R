#' Function to reset acd transformation for variables with few values
#' 
#' To be used in [fit_mfp()].
#' This function resets the `acdx` parameter (logical vector) of variables with
#' less than 5 distinct values to `FALSE`.
#' 
#' @param x a design matrix of dimension nobs x nvars where nvars is the number 
#' of predictors excluding an intercept.  
#' @param acdx a named logical vector of length nvars indicating continuous
#' variables to undergo the approximate cumulative distribution (ACD) 
#' transformation. May be ordered differently than the columns of `x`.
#' 
#' @return 
#' Logical vector of same length as `acdx`.
reset_acd <- function(x, acdx) {
  
  names_acd <- names(acdx)[which(acdx == TRUE)]
  
  # number of unique values of each column in acdx
  n_unique <- apply(x[, names_acd, drop = FALSE], 2, 
                    function(col) length(unique(col)))
  ind_reset <- which(n_unique < 5)
  
  if (length(ind_reset) > 0) {
    acdx[names_acd][ind_reset] <- FALSE
    warning("i For any variable with fewer than 5 unique values no acd transformation can be performed.\n", 
            sprintf("i The requested acd transform has been reset to FALSE for the following variables: %s.", 
                    paste0(names_acd[ind_reset], collapse = ", ")))
  }
  
  acdx
}
