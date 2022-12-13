# Generates FP data for x of interest.
# x = a vector of a predictor in which FP data is required. x is assumed to have
# shifted and scaled
# powers = set of FP powers
# degree = degree of FP. 1 = FP1, 2 = FP2 and so on
generate_fp_data <- function(x, degree, powers) {
  # Possible combination of powers given degree
  combs <- generate_fp_powers(degree = degree, powers = powers)
  nfp <- dim(combs)[1]
  # Save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in 1:nfp) {
    fpdt[[i]] <- generate_fp(
      x = x, power = combs[i, ], scale = 1,
      shift = 0, center = F
    )
  }
  # Return a list of length of tranformed variables for instance in FP2 we have
  # a list of length 36, each list having a n*2 matrix, when default s is used
  return(list(data = fpdt, powers = combs))
}
