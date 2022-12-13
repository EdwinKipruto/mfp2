# calculates p-values for Chi-square distribution
# dev = a vector of deviance for models i.e Null, linear, FP1,...FPm
calculate_chisquare_test <- function(dev, acd) {
  if (acd) {
    df <- c(4, 3, 2, 2, 1)
    # we have dev = (NULL, M4, M2, M3, M1, M5)
    # calculate p-values for the functions selection procedure:
    # a). M1 vs Null     b). M1 vs M4     c). M1 vs M2    d) M1 vs M3
    pvalues <- dev.diff <- numeric(5)
    for (i in 1:4) {
      dev.diff[i] <- dev[i] - dev[5]
      pvalues[i] <- pchisq(dev.diff[i], df = df[i], lower.tail = F)
    }
    # e). M3 vs M5
    dev.diff[5] <- dev[6] - dev[4]
    pvalues[5] <- pchisq(dev.diff[5], df = df[5], lower.tail = F)
  } else {
    # we have dev = c(NULL, Linear, FP1, FP2,....FPm)
    nn <- length(dev)
    # Maximum permitted degree: nn-2, the first and second position is null and linear
    m <- nn - 2
    # calculate degrees of freedom for testing. Note that:
    # FPm vs Null = 2m; FPm vs linear = 2m-1; FPm vs FP1 = 2m-2; FPm vs FP2 = 2m-4
    # in general: FPm vs FPk = 2(m-k) for m>k
    k <- seq(m - 1, 1)
    df <- c(2 * m, 2 * m - 1, if (m >= 2) {
      2 * (k)
    }) # same as 2*(m-k) where k = seq(1,m-1)
    # calculate p-values for the functions selection procedure:
    # FPm vs Null; FPm vs linear; FPm vs FP1; FPm vs FP2 etc.
    pvalues <- dev.diff <- numeric(nn - 1)
    for (i in 1:(nn - 1)) {
      dev.diff[i] <- dev[i] - dev[nn]
      pvalues[i] <- pchisq(dev.diff[i], df = df[i], lower.tail = F)
    }
  }
  return(list(pvalues = pvalues, dev.diff = dev.diff))
}
