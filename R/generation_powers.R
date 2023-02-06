#' Function that generates a matrix of FP powers for any degree
#'
#' @param degree The degree of fractional polynomial; degree = 1 is FP1 and
#' returns 8 powers; degree 2 is FP2 and returns 36 pairs of powers; degree 3 is
#' FP3 and returns 120 triples of powers etc. If the ACD transformation is used,
#' this degree is assumed to be 2. #'
#' @param powers the set of allowed powers for the fractional polynomials. 
#' Default is `NULL` and the set {-2, -1, -0.5, 0, 0.5, 1, 2, 3} is used.
#' 
#' @details 
#' For FP powers, this function returns all combinations of the powers of 
#' length `degree`, that is all pairs in which each entry is taken from the 
#' set `powers`, but no pair is repeated (i.e. the order of the entries does 
#' not matter).
#' Thus, for the default set of powers and degree 2, this function returns
#' 36 combinations.
#' 
#' For ACD powers, this function simply returns all possible tuples of 
#' powers of length n. 
#' Thus, for the default set of powers, this function returns 8 possible
#' powers, and for degree 2 it returns 64 pairs of powers. Higher degrees
#' are not supported by the function. 
#' 
#' @return 
#' A matrix of powers with degree columns and rows depending on the `degree`.
#' 
#' @import arrangements
#' @export
generate_powers_fp <- function(degree = NULL, 
                               powers = NULL) {
  if (is.null(degree)) {
    degree <- 2
  }
  if (is.null(powers)) {
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  
  # using replacement because powers may be repeated e.g. (0,0) or (1,1)
  arrangements::combinations(x = powers, k = degree, replace = TRUE)
}

#' @describeIn generate_powers_fp Function to generate acd powers.
#' @export
generate_powers_acd <- function(degree = NULL, 
                                powers = NULL) {
  if (is.null(degree)) {
    degree <- 2
  }
  if (is.null(powers)) {
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  
  if (degree == 1) {
    return(as.matrix(powers, ncol = 1))
  }
  
  matrix(as.matrix(expand.grid(powers, powers)), ncol = 2)
}
