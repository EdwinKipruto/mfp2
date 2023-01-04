#' Function that generates a matrix of FP powers for any degree
#'
#' @param degree The degree of fractional polynomial; degree = 1 is FP1 and
#' returns 8 powers; degree 2 is FP2 and returns 36 pairs of powers; degree 3 is
#'  FP3 and returns 120 triples of powers etc.
#' @param powers the set of allowed powers for the fractional polynomials. 
#' Default is `NULL` and the set {-2, -1, -0.5, 0, 0.5, 1, 2, 3} is used.
#' 
#' @return 
#' A matrix of powers with degree columns and rows depending on the degree.
#' 
#' @importFrom arrangements  combinations
#' @export
generate_fp_powers <- function(degree = NULL, 
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
