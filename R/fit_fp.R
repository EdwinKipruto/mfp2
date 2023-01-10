#' Function to fit univariable FP1 models for acd transformation
#' 
#' To be used in [fit_acd()].
#' 
#' @return 
#' The best FP power with smallest deviance.
find_best_fp1_for_acd <- function(x, 
                                  y, 
                                  powers) {
  
    if (!is.null(dim(x))) 
    stop("! `x` must be a vector.")
  
  # generate all possible FP1 transformations
  trafo <- generate_transformations_fp(x = x, degree = 1, powers = powers)$data
  
  # Fit linear models for each FP1 function
  n_powers <- length(powers)
  # store deviance and model object
  devs <- vector("numeric", n_powers)
  fits <- vector("list", n_powers)
  
  for (i in seq_len(n_powers)) {
    # fit linear model
    fit <- fit_model(x = trafo[[i]], y = y, family = "gaussian")
    
    devs[i] <- -2 * fit$logl
    fits[[i]] <- fit$fit
  }
  
  # find best model
  index <- which.min(devs)
  
  list(
    power = powers[[index]], 
    fit = fits[[index]]
  )
}
