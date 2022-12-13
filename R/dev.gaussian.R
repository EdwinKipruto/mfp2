deviance_gaussian <- function(RSS, weights, n) {
  # calculate lognormalized weights
  if (length(unique(weights)) == 1) {
    meanwts <- 0
  } else {
    logwts <- log(weights)
    mlogwts <- mean(logwts)
    # lognormalize weights
    lognormweights <- logwts / mlogwts
    # average of lognormalized weights
    meanwts <- mean(lognormweights)
  }
  # use the formula found here https://www.stata.com/manuals/rfp.pdf
  k <- log((2 * pi * RSS) / n)
  devx <- n * (1 - meanwts + k)
  return(devx)
}
