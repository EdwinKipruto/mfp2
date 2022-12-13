calculcate_f_statistic_stata <- function(dev.reduced, dev.full, d1, d2, n) {
  # Find deviance difference
  devdiff <- dev.reduced - dev.full
  # use formula on page 23 from here: https://www.stata.com/manuals/rfp.pdf
  aa <- exp(devdiff / n)
  fstatx <- (d2 / d1) * (aa - 1)
  # pvalue
  if (devdiff == 0) {
    pval <- 1
  } else {
    pval <- stats::pf(fstatx, df1 = d1, df2 = d2, lower.tail = F)
  }
  return(list(fstatistic = fstatx, pval = pval, dev.diff = devdiff))
}
