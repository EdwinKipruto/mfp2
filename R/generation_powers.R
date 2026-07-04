#' Generate fractional-polynomial power combinations
#'
#' Generates the candidate power matrix used for ordinary fractional-polynomial
#' transformations.
#'
#' @param degree Integer degree of the fractional polynomial. For example,
#'   `degree = 1` generates FP1 powers, `degree = 2` generates FP2 pairs,
#'   and `degree = 3` generates FP3 triples. If `NULL`, the default is `2`.
#' @param powers Numeric vector of allowed fractional-polynomial powers. If
#'   `NULL`, the default set `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)` is used.
#'
#' @details
#' `generate_powers_fp()` returns all combinations with replacement of length
#' `degree` from the supplied power set. The order of powers within a row does
#' not matter, and each row is sorted in increasing order. Repeated powers are
#' allowed, for example `(0, 0)` or `(1, 1)`, corresponding to repeated-power
#' fractional-polynomial terms.
#'
#' With the default set of eight powers, `degree = 1` gives 8 candidate powers,
#' `degree = 2` gives 36 candidate pairs, and `degree = 3` gives 120 candidate
#' triples.
#'
#' If `degree = 0`, the function returns a one-row, one-column matrix containing
#' `1`, representing the null transformation.
#'
#' @return
#' A numeric matrix. For `degree > 0`, the matrix has `degree` columns and one
#' row for each candidate power combination. For `degree = 0`, the matrix has
#' one column containing `1`.
#'
#' @examples
#' powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
#'
#' generate_powers_fp(degree = 1, powers = powx)
#' generate_powers_fp(degree = 2, powers = powx)
#'
#' # A single supplied power with degree 2 gives the repeated-power FP2 basis.
#' generate_powers_fp(degree = 2, powers = 2)
#'
#' @keywords internal
#' @noRd
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
  
  # using replacement because powers may be repeated e.g. (0,0) or (1,1)
  generate_combinations_with_replacement(powers, degree)
}

#' Generate ACD power combinations
#'
#' Generates the candidate power matrix used for approximate cumulative
#' distribution (ACD) transformations.
#'
#' @param degree Integer ACD degree. Supported values are `0`, `1`, and `2`.
#'   If `NULL`, the default is `2`.
#' @param powers Numeric vector of allowed powers. If `NULL`, the default set
#'   `c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)` is used.
#'
#' @details
#' `generate_powers_acd()` returns candidate power pairs for ACD transformations.
#' Unlike ordinary fractional-polynomial powers, ACD powers are generated as
#' Cartesian products, so order matters.
#'
#' The returned matrix always has two columns. The first column corresponds to
#' the untransformed-data component and the second column corresponds to the ACD
#' component. For `degree = 0` and `degree = 1`, entries in the first column are
#' set to `NA` to indicate that the untransformed-data component is not used.
#'
#' For the default set of eight powers, `degree = 1` gives 8 candidate rows and
#' `degree = 2` gives 64 candidate pairs. Degrees greater than 2 are not
#' supported.
#'
#' @return
#' A two-column numeric matrix of ACD power specifications. Rows containing
#' `NA` in the first column indicate that the untransformed-data component is
#' not included.
#'
#' @examples
#' powx <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
#'
#' generate_powers_acd(degree = 0, powers = powx)
#' generate_powers_acd(degree = 1, powers = powx)
#' generate_powers_acd(degree = 2, powers = powx)
#'
#' @keywords internal
#' @noRd
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
  
  if (degree == 0) {
    return(matrix(c(NA, 1), ncol = 2))
  }
  
  if (degree == 1) {
    return(matrix(c(rep(NA, length(powers)), powers), ncol = 2))
  }
  
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
#' @keywords internal
#' @noRd
generate_combinations_with_replacement <- function(x,
                                                   k) {
  
  # if (k > 5) {
  #   warning("FP degree higher than 5; the MFP algorithm may take a while to do model selection.")
  # }
  
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
