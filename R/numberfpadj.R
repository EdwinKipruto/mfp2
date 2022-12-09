# Calculates total number of fp powers in adjustment variables
# x is a list of fp powers for all variables
numberfpadj <- function(x) {
  # Remove NA-variables eliminated in usual mfp. In MFPA, a variable can have
  # powers c(1,NA) or c(NA,1) and its unaffected
  x <- x[!is.na(x)]
  if (!is.null(x)) {
    # Remove NAs in powers (2,NA) or (NA,1) etc. Happens because of acd
    x <- lapply(x, function(x) x[!is.na(x)])
    # Find variables having powers equal to 1
    ft <- unlist(lapply(x, function(x) identical(x, 1)))
    # If all ft are TRUE then no fp power in adjustment variables
    if (all(ft)) {
      # Total number of FP powers in adjustment are 0
      tFP <- 0
    } else {
      # Powers of variables having FP powers
      fp.ok <- x[!ft]
      # Total estimated FP powers in adjustment variables
      tFP <- sum(unlist(lapply(fp.ok, length)), na.rm = T)
    }
    # if all adjustment variables were eliminated then tFP = 0
  } else {
    tFP <- 0
  }
  return(tFP)
}
#
#
#
# x = list(a = c(NA,NA),b = 3, c = c(NA,1), d = 1, e = 0, f = c(1,1))
# numberfpadj(x)
#
#
# lapply(x, function(x) x[!is.na(x)])
# ax = list(a)
# pw[!is.na(pw)]
# xx =c(NA,1)
# xx[!is.na(xx)]
