# calculates p-values for f distribution
# sse = a vector of residual sum of squares for null, linear, fp1, fp2,...fpm
# df  = a vector of degrees of freedom for each sse
calculate_f_test <- function(sse, df, acd) {
  if (acd) {
    # we have sse = (NULL, M4, M2, M3, M1, M5)
    # calculate p-values for the functions selection procedure:
    # a). M1 vs Null     b). M1 vs M4     c). M1 vs M2    d) M1 vs M3
    pvalues <- fstatistic <- numeric(5)
    for (i in 1:4) {
      fsttas <- calculcate_f_statistic(sse.reduced = sse[i], sse.full = sse[5], df.reduced = df[i], df.full = df[5])
      fstatistic[i] <- fsttas$fstatistic
      pvalues[i] <- fsttas$pval
    }
    # e). M3 vs M5
    fsttas <- calculcate_f_statistic(sse.reduced = sse[6], sse.full = sse[4], df.reduced = df[6], df.full = df[4])
    fstatistic[5] <- fsttas$fstatistic
    pvalues[5] <- fsttas$pval
  } else {
    # we have sse = c(NULL, Linear, FP1, FP2,....FPm)---usual mfp approach
    nn <- length(sse)
    # Maximum permitted degree: nn-2, the first and second position is null and lin
    m <- nn - 2
    # calculate p-values for the functions selection procedure:
    # FPm vs Null; FPm vs linear; FPm vs FP1; FPm vs FP2...FPm vs FPm-1 etc.
    pvalues <- fstatistic <- numeric(nn - 1)
    for (i in 1:(nn - 1)) {
      fsttas <- calculcate_f_statistic(sse.reduced = sse[i], sse.full = sse[nn], df.reduced = df[i], df.full = df[nn])
      fstatistic[i] <- fsttas$fstatistic
      pvalues[i] <- fsttas$pval
    }
  }

  return(list(pvalues = pvalues, fstatistic = fstatistic))
}

# model.pvalues.ftest(sse = c(230, 240, 240), df = c(20,10,10), acd = F)
# pf(fstatistic, df1 = df.reduced-df.full, df2 = df.full, lower.tail = F)
#
# num <- (sse.reduced-sse.full)/(df.reduced-df.full)
# denom <- sse.full/(df.full)
# fstatistic <- num/denom
# pval <- pf(0, df1 = 0, df2 = 240, lower.tail = F)
# return(list(fstatistic = fstatistic, pval = pval))
