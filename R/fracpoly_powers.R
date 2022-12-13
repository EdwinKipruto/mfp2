#' Function that generates a matrix of FP powers for any degree
#'
#' @param degree The degree of fractional polynomial; degree = 1 is FP1 and
#' returns 8 powers; degree 2 is FP2 and returns 36 pairs of powers; degree 3 is
#'  FP3 and returns 120 triple of powers etc.
#' @param powers the set of fractional polynomial. Default is NULL in which the set
#' s = {-2, -1, -0.5, 0, 0.5, 1, 2, 3} is used
#' @return A matrix of powers depending on the degree
#' @importFrom arrangements  combinations
#' @export
generate_fp_powers <- function(degree = NULL, powers = NULL) {
  # set FP2 as the default maximum degree permitted
  if (is.null(degree)) {
    degree <- 2
  }
  # set default FP set when not provided
  if (is.null(powers)) {
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  # Find all possible combination of set s. Here we use combinations with
  # repetition because we allow powers to be repeated e.g (0,0) or (1,1) etc
  all_powers <- arrangements::combinations(x = powers, k = degree, replace = T)
  return(all_powers)
}
