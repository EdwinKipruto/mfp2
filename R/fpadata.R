# Generates FPa data for x of interest. It utilizes acd(x) function
# x = a vector of a predictor in which FP data is required. x is assumed to have
# shifted and scaled
# powers = set of FP powers
fpadata <- function(x, powers) {
  # Possible combination of powers given degree = 2 i.e p1,p2
  # we set degree = 2 because we are interested in two powers for x and ax = acd(x)
  combs <- fracpoly_powers(degree = 2, powers = powers)
  # interchange the columns of combs, select unrepeated powers and rbind
  combs1 <- combs[, 2:1]
  keep <- apply(combs1, 1, function(x) length(unique(x[!is.na(x)])) != 1)
  # Final pairs of powers.if default FP set is used then its length is 64
  combs2 <- rbind(combs, combs1[keep, ])
  nfp <- dim(combs2)[1L]
  # Save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in seq_len(nfp)) {
    fpdt[[i]] <- afracgen(
      x = x, power = combs2[i, ], scale = 1,
      shift = 0, center = F
    )
  }
  return(list(data = fpdt, powers = combs2))
}
