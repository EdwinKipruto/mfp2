#' Function to generate all requested FP transformations for a single variable
#' 
#' @param x a numeric vector of length nobs assumed to have been shifted and 
#' scaled.
#' @param degree numeric indicating the degree of FPs. Assumed to be 2 for acd
#' transformation.
#' @param powers a vector of allowed FP powers.
#' 
#' @details 
#' Any FP transformation is given by a vector of powers, e.g. (p1, p2) for 
#' degree 2. These correspond to powers x^p1 and x^p2. Thus, we only need to 
#' consider combinations of all values in `powers`, since order of the entries
#' does not matter. See [generate_powers_fp()]. 
#' A special case are repeated powers, i.e. p1 = p2. In this case, the repeated 
#' entries are multiplied by log(x) (see [transform_vector_fp()]).
#' 
#' When the ACD transformation is requested, then all pairs of length 2
#' are considered, i.e. 64. See [generate_powers_acd()].
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `data`: a list with length equal to the number of possible FPs for the 
#' variable of interest. Each entry is a matrix with degree many columns, 
#' and nobs observations comprising the FP transformed input variable. 
#' For example, for degree = 2 and nobs = 10, each entry is a 10 x 2 matrix.
#' * `powers`: the associated FP powers for each entry in data. 
generate_transformations_fp <- function(x, 
                                        degree, 
                                        powers) {
  
  # all possible combination of powers given degree
  combs <- generate_powers_fp(degree = degree, powers = powers)
  nfp <- dim(combs)[1]
  
  # save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in 1:nfp) {
    fpdt[[i]] <- transform_vector_fp(
      x = x, power = combs[i, ], scale = 1, shift = 0, center = FALSE
    )
  }

  list(
    data = fpdt,
    powers = combs
  )
}

#' @describeIn generate_transformations_fp Function to generate acd transformations.
generate_transformations_acd <- function(x, 
                                         powers) {

  # all possible pairs of powers
  combs <- generate_powers_acd(powers = powers)
  nfp <- dim(combs)[1L]
  
  # Save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in seq_len(nfp)) {
    fpdt[[i]] <- transform_vector_acd(
      x = x, power = combs[i, ], scale = 1, shift = 0, center = FALSE
    )
  }
  
  list(
    data = fpdt, 
    powers = combs
  )
}
