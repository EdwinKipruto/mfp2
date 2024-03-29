#' Calculates the total number of fractional polynomial powers in adjustment variables.
#'
#' This function takes a list \code{x} containing fractional polynomial powers for all variables
#' and calculates the total number of powers across the variables.
#'
#' @param x A list of fractional polynomial powers for all variables.
#' 
#' @return Numeric value denoting total number of fractional polynomial powers in the adjustment
#' variables.
calculate_number_fp_powers <- function(x) {
  # Remove NA-variables eliminated in usual mfp. In mfpa, a variable can have
  # powers c(1,NA) or c(NA,1) and its unaffected
  x <- x[!is.na(x)]
  
  if (!is.null(x)) {
    # Remove NAs in powers (2,NA) or (NA,1) etc. Happens because of acd
    x <- lapply(x, function(x) x[!is.na(x)])
    
    # Find variables having powers equal to 1
    ft <- unlist(lapply(x, function(x) identical(x, 1)))
    
    # If all ft are TRUE then no fp power in adjustment variables
    if (all(ft)) {
      # Total number of FP powers in adjustment are 0
      tFP <- 0
    } else {
      # Powers of variables having FP powers
      fp.ok <- x[!ft]
      
      # Total estimated FP powers in adjustment variables
      tFP <- sum(unlist(lapply(fp.ok, length)), na.rm = TRUE)
    }
    # if all adjustment variables were eliminated then tFP = 0
  } else {
    tFP <- 0
  }
  return(tFP)
}
