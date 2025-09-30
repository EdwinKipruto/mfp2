#' Function that generates a matrix of FP powers for any degree
#'
#' @param degree The degree of fractional polynomial. For example,
#' degree = 1 is FP1 and returns 8 powers; degree 2 is FP2 and 
#' returns 36 pairs of powers; degree 3 is FP3 and returns 120 
#' triples of powers, and so on. If the ACD transformation is used,
#' this degree is assumed to be 2.
#' @param powers the set of allowed powers for the fractional polynomials. 
#' Default is `NULL` and the set \eqn{(-2, -1, -0.5, 0, 0.5, 1, 2, 3)} is used.
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
#' are not supported by the function. In case that `degree = 0` or `degree = 1`, 
#' the first column of the matrix representing untransformed data are set to 
#' `NA` to indicate that the normal data do not play a role. Higher degrees
#' than two are not supported. 
#' 
#' @examples
#' powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
#' generate_powers_fp(degree = 2, powers = powx)
#' generate_powers_acd(degree = 2, powers = powx)
#' 
#' @return 
#' A matrix of powers with degree columns and rows depending on the `degree`.
#' For ACD powers always a matrix with two columns. For normal fps each row
#' will be sorted in increasing order (in alignment with
#' how \code{transform_vector_fp()} processes the data).
#' 
#' @export
generate_powers_fp <- function(degree = NULL, 
                               powers = NULL) {
  if (is.null(degree)) {
    degree <- 2
  }
  
  if (is.null(powers)) {
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  
  if (degree == 0) {
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
   # combination below does not work if x is of length 1
     if (length(powers) == 1){
    return(matrix(powers, nrow = 1, ncol = 1))
    }
  
  # using replacement because powers may be repeated e.g. (0,0) or (1,1)
  generate_combinations_with_replacement(powers, degree)
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
  
  if (degree == 0)
    return(matrix(c(NA, 1), ncol = 2))
  
  if (degree == 1) 
    return(matrix(c(rep(NA, length(powers)), powers), ncol = 2))
  
  matrix(as.matrix(expand.grid(powers, powers)), ncol = 2)
}

#' Helper function to generate combinations with replacement
#' 
#' This very simple helper generates combinations with replacement. 
#' 
#' @param x vector of elements to choose from.
#' @param k number of elements to choose. 
#' 
#' @details 
#' This is replicating the functionality from `arrangements::combinations` with
#' `replace = TRUE`. Note that base R function `utils::combn` only returns
#' combinations without replacement, thus pairs like (0, 0) are not in the 
#' output. 
#' 
#' Note that this function is extremely inefficient and only intended to be
#' used with small use cases, i.e. small k. This is typically the case in the
#' context of MFP, but a warning is given if this is not the case since the 
#' algorithm may take a while to compute the combinations, and even longer
#' to do model selection.
#' 
#' @return 
#' A m x k matrix, where m is the number of combinations.
generate_combinations_with_replacement <- function(x, 
                                                   k) {
  
  if (k > 5) {
    warning("i FP degree higher than 5, the MFP algorithm may take a while to do model selection.")
  }
  
  # generate all possible pairs
  pairs <- expand.grid(rep(list(x), k))
  
  # and remove duplicates
  # reverse column order to conform with ordering of arrangements::combinations
  # i.e. we have identical results
  unname(as.matrix(pairs[!duplicated(t(apply(pairs, 1, sort))), rev(1:k)]))
}
