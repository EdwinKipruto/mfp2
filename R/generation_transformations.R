#' Function to generate all requested FP transformations for a single variable
#' 
#' @param x a numeric vector of length nobs assumed to have been shifted and 
#' scaled.
#' @param degree numeric indicating the degree of FPs. Assumed to be 2 for acd
#' transformation.
#' @param powers a vector of allowed FP powers.
#' 
#' @return 
#' A list with two entries: 
#' 
#' * `data`: a list with length equal to the number of possible FPs for the 
#' variable of interest. Each entry is a matrix with as degree many columns, 
#' and nobs observations comprising the FP transformed input variable. 
#' For example, for degree = 2 and nobs = 10, each entry is a 10 x 2 matrix.
#' * `powers`: the associated FP powers for each entry in data. 
generate_fp_transformations <- function(x, 
                                        degree, 
                                        powers) {
  
  # all possible combination of powers given degree
  combs <- generate_fp_powers(degree = degree, powers = powers)
  nfp <- dim(combs)[1]
  
  # save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in 1:nfp) {
    fpdt[[i]] <- transform_fp_vector(
      x = x, power = combs[i, ], scale = 1,
      shift = 0, center = FALSE
    )
  }
  # Return a list of length of tranformed variables for instance in FP2 we have
  # a list of length 36, each list having a n*2 matrix, when default s is used
  
  list(
    data = fpdt,
    powers = combs
  )
}

#' @describeIn generate_fp_transformations Function to generate acd transformations.
generate_acd_transformations <- function(x, powers) {
  # Possible combination of powers given degree = 2 i.e p1,p2
  # we set degree = 2 because we are interested in two powers for x and ax = acd(x)
  combs <- generate_fp_powers(degree = 2, powers = powers)
  # interchange the columns of combs, select unrepeated powers and rbind
  combs1 <- combs[, 2:1]
  keep <- apply(combs1, 1, function(x) length(unique(x[!is.na(x)])) != 1)
  # Final pairs of powers.if default FP set is used then its length is 64
  combs2 <- rbind(combs, combs1[keep, ])
  nfp <- dim(combs2)[1L]
  # Save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in seq_len(nfp)) {
    fpdt[[i]] <- transform_acd_vector(
      x = x, power = combs2[i, ], scale = 1,
      shift = 0, center = F
    )
  }
  
  list(
    data = fpdt, 
    powers = combs2
  )
}
