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
  
  
  # Internal contract: ACD supports degree 0, 1, or 2 only.
  if (length(degree) != 1L || is.na(degree) || degree < 0L || degree > 2L) {
    stop("Internal error: `degree` for ACD must be 0, 1, or 2.", call. = FALSE)
  }
  
  
  if (degree == 0)
    return(matrix(c(NA, 1), ncol = 2))
  
  if (degree == 1) 
    return(matrix(c(rep(NA, length(powers)), powers), ncol = 2))
  
  matrix(as.matrix(expand.grid(powers, powers)), ncol = 2)
}

#' Helper function to generate combinations with replacement
#'
#' This helper generates combinations with replacement.
#'
#' @param x Vector of elements to choose from.
#' @param k Number of elements to choose.
#'
#' @details
#' This function replicates the functionality of
#' arrangements::combinations(x, k, replace = TRUE). Unlike
#' utils::combn(x, k), which returns combinations without replacement, this
#' function allows repeated elements, so combinations such as (0, 0) are
#' included in the output.
#'
#' Internally, the function uses a stars-and-bars indexing construction. It
#' first applies utils::combn() to the sequence 1:(length(x) + k - 1) and
#' then shifts the resulting indices by 1:k. This produces the indices needed
#' to select nondecreasing combinations from the sorted input vector x.
#'
#' This avoids generating all ordered tuples with expand.grid() and then
#' removing duplicates. However, the number of combinations can still grow
#' quickly with increasing k. In the MFP context, high FP degrees correspond
#' to a large number of possible FP power combinations, and the subsequent model
#' selection step may therefore be computationally intensive. A warning is
#' issued for k > 5.
#'
#' @return
#' A matrix with one row per combination and k columns.
generate_combinations_with_replacement <- function(x,
                                                   k) {
  
  if (k > 5) {
    warning("FP degree higher than 5; the MFP algorithm may take a while to do model selection.")
  }
  
  # Sort input so that returned combinations are ordered consistently.
  x <- sort(x)
  n <- length(x)
  
  # Generate combinations of positions using the stars-and-bars representation.
  # This gives combinations with replacement without constructing all ordered
  # tuples first.
  idx <- utils::combn(seq_len(n + k - 1), k)
  
  # Shift the indices by 1:k to map stars-and-bars positions back to indices
  # of the original sorted vector x.
  idx <- idx - seq_len(k)
  
  # Extract the selected elements. The matrix is first formed with one
  # combination per column, then transposed so that each row is one combination.
  out <- matrix(x[idx + 1], nrow = k)
  
  t(out)
  
}
