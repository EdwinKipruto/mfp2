# standard error of a fitted line
# vcox is the variance-covariance of an mfpa model
# xnames are the names of transformed variables of x of interest
# family error distribution
# X = subsetted Xtransformed for variable of interest
vcalc <- function(model, X) { # Add complete = F, maybe
  # The first column is the intercept in glm
  vcovx <- vcov(object = model)
  # Get rid of variance and covariance of intercept if any
  xnames <- colnames(X)
  ind <- match(xnames, colnames(vcovx))
  vcovx1 <- vcovx[ind, ind, drop = FALSE] # vcovx[xnames, xnames, drop = FALSE]
  # variance of the fitted line
  nx <- length(xnames)
  v1 <- 0
  for (i in 1:nx) {
    for (j in 1:i) {
      v1 <- v1 + ((j < i) + 1) * vcovx1[i, j] * X[, i] * X[, j]
    }
  }
  # deal with intercept
  family <- model$family
  if (is.character(family)) { # cox comes as character
    family <- family
  } else {
    family <- family$family # glm comes as functions
  }
  if (family != "cox") {
    ind <- c(1, ind)
    vcovx2 <- vcovx[ind, ind, drop = F]
    v2 <- 0
    for (i in seq_along(xnames)) {
      # we assume intercept is in column 1 of vcov which is always the case
      v2 <- v2 + 2 * vcovx2[1, i + 1] * X[, i]
    }
    v1 <- vcovx2[1, 1] + v2 + v1
  }
  return(sqrt(v1))
}
