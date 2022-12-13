# A function that calculates f statistic and p-value
# mfp in stata uses different formula see https://www.stata.com/manuals/rfp.pdf
calculcate_f_statistic <- function(sse.reduced, sse.full, df.reduced, df.full) {
  # see page 297 or equation next to 7.4 in http://users.stat.ufl.edu/~winner/sta4211/ALSM_5Ed_Kutner.pdf
  if (df.reduced == df.full) { # it can happen that the best fp1 is linear and we are testing between the same models
    fstatistic <- 0
    pval <- 1
  } else {
    num <- (sse.reduced - sse.full) / (df.reduced - df.full)
    denom <- sse.full / (df.full)
    fstatistic <- num / denom
    pval <- stats::pf(fstatistic, df1 = df.reduced - df.full, df2 = df.full, lower.tail = F)
  }
  return(list(fstatistic = fstatistic, pval = pval))
}
