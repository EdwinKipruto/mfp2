# calculates p-values for F distribution based on Royston formula on page 23 of
# stata manual found here: https://www.stata.com/manuals/rfp.pdf
# dev = a vector of deviance for models i.e Null, linear, FP1,...FPm
# resid.df = a vector of residual degrees of freedom for models i.e Null, linear, FP1,...FPm
# n = sample size/number of observations
# acd = TRUE/FALSE
calculate_f_test_royston <- function(dev, resid.df, n, acd) {
  if (acd) {
    df <- c(4, 3, 2, 2, 1)
    # CHECK ON THIS SECTION, THE OTHER SECTION HAS BEEN CORRECTED
    # we have dev = (NULL, M4, M2, M3, M1, M5)
    # calculate p-values for the functions selection procedure:
    # a). M1 vs Null     b). M1 vs M4     c). M1 vs M2    d) M1 vs M3
    # number of FP powers estimated in NULL = 0, M4 = 0, M2 = 1, M3 = 1, M1 = 2, M5 = 0
    pwrs <- c(0, 0, 1, 1, 2, 0)
    pvalues <- dev.diff <- fstatistic <- numeric(5)
    for (i in 1:4) {
      stats <- calculcate_f_statistic_stata(
        dev.reduced = dev[i], dev.full = dev[5], d1 = df[i],
        d2 = resid.df[5], n = n
      )
      fstatistic[i] <- stats$fstatistic
      pvalues[i] <- stats$pval
      dev.diff[i] <- stats$dev.diff
    }
    # e). M3 vs M5
    stats2 <- calculcate_f_statistic_stata(
      dev.reduced = dev[6], dev.full = dev[4], d1 = df[5],
      d2 = resid.df[4], n = n
    )
    fstatistic[5] <- stats2$fstatistic
    pvalues[5] <- stats2$pval
    dev.diff[5] <- stats2$dev.diff
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
    # number of FP powers estimated in null = 0, linear = 0, fp1 =1, fp2 = 2,...
    pwrs <- c(0, 0, seq_len(m))
    # calculate p-values for the functions selection procedure:
    # FPm vs Null; FPm vs linear; FPm vs FP1; FPm vs FP2 etc.
    pvalues <- dev.diff <- fstatistic <- numeric(nn - 1)
    for (i in 1:(nn - 1)) {
      stats <- calculcate_f_statistic_stata(
        dev.reduced = dev[i], dev.full = dev[nn], d1 = df[i],
        d2 = resid.df[nn], n = n
      )
      pvalues[i] <- stats$pval
      fstatistic[i] <- stats$fstatistic
      dev.diff[i] <- stats$dev.diff
    }
  }
  return(list(pvalues = pvalues, dev.diff = dev.diff, fstatistic = fstatistic))
}
