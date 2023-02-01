#' Function to calculate p-values for likelihood-ratio test
#' 
#' @param logl a numeric vector of length 2 with log-likelihoods. Typically 
#' ordered in increasing order (i.e. null model first, then full model) and 
#' used to test the ratio `logl[1] / logl[2]`.
#' @param dfs a numeric vector with degrees of freedom.
#'  
#' @details 
#' Uses Wilk's theorem that -2log(LR) (LR = likelihood ratio) asymptotically
#' approaches a Chi-square distribution under the null hypothesis that both 
#' likelihoods are equal. 
#' 
#' Model likelihoods can then be compared by computing 
#' D = -2 log(likelihood reduced model / likelihood full model), and then 
#' use a Chi-square distribution with df_full - df_reduced many degrees
#' of freedom to derive a p-value.
#' 
#' This is basically the same way as [stats::anova()] implements the 
#' likelihood ratio test.
#' 
#' @return
#' The p-value for the likelihood ratio test for the ratio `logl[1] / logl[2]`.
calculate_lr_test <- function(logl, dfs) {
  pchisq(2 * (logl[2] - logl[1]) , 
         df = dfs[2] - dfs[1], 
         lower.tail = FALSE)
}

#' Function to calculate p-values for Chi-square distribution
#' 
#' @param dev a vector of deviance for models i.e Null, linear, FP1,...FPm.
#' @param acd logical indicating the use of acd transformation.
calculate_chisquare_test <- function(dev, 
                                     acd) {
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
  
  list(
    pvalues = pvalues, 
    dev.diff = dev.diff
  )
}

#' Function to compute F-statistic and p-value
#' 
#' @details 
#' `mfp` in Stata uses different formula, see
#'  https://www.stata.com/manuals/rfp.pdf. That formula is implemented in 
#'  [calculate_f_statistic_stata()].
#'  
#'  This functions uses page 297 or equation next to 7.4 in 
#'  http://users.stat.ufl.edu/~winner/sta4211/ALSM_5Ed_Kutner.pdf.
#'  
#'  @seealso 
#'  [calculate_f_statistic_stata()]
calculate_f_statistic <- function(sse_reduced, 
                                   sse_full,
                                   df_reduced, 
                                   df_full) {
  
  # it can happen that the best fp1 is linear
  # so we are testing between the same models
  if (df_reduced == df_full) { 
    fstatistic <- 0
    pval <- 1
  } else {
    num <- (sse_reduced - sse_full) / (df_reduced - df_full)
    denom <- sse_full / (df_full)
    fstatistic <- num / denom
    pval <- stats::pf(fstatistic, df1 = df_reduced - df_full,
                      df2 = df_full, lower.tail = FALSE)
  }
  
  list(
    fstatistic = fstatistic,
    pval = pval
  )
}

#' Function to compute F-statistic as defined in `mfp` in Stata 
#' 
#' Alternative to [calculate_f_statistic()].
#' 
#' @details 
#' Uses formula on page 23 from here: https://www.stata.com/manuals/rfp.pdf.
#' 
#' @seealso 
#' [calculate_f_statistic()]
calculate_f_statistic_stata <- function(dev_reduced, 
                                         dev_full,
                                         d1, 
                                         d2,
                                         n) {
  devdiff <- dev_reduced - dev_full
  aa <- exp(devdiff / n)
  fstatx <- (d2 / d1) * (aa - 1)
  if (devdiff == 0) {
    pval <- 1
  } else {
    pval <- stats::pf(fstatx, df1 = d1, df2 = d2, lower.tail = FALSE)
  }
  
  list(
    fstatistic = fstatx, 
    pval = pval, 
    dev.diff = devdiff
  )
}

#' Function to calculates p-values for F-distribution
#' 
#' @param sse a vector of residual sum of squares for null, linear, fp1,
#' fp2 to fpm.
#' @param df  a vector of degrees of freedom for each sse.
#' @param acd logical indicating the use of acd transformation.
calculate_f_test <- function(sse, 
                             df, 
                             acd) {
  if (acd) {
    # we have sse = (NULL, M4, M2, M3, M1, M5)
    # calculate p-values for the functions selection procedure:
    # a). M1 vs Null     b). M1 vs M4     c). M1 vs M2    d) M1 vs M3
    pvalues <- fstatistic <- numeric(5)
    for (i in 1:4) {
      fsttas <- calculate_f_statistic(sse_reduced = sse[i], 
                                       sse_full = sse[5],
                                       df_reduced = df[i],
                                       df_full = df[5])
      fstatistic[i] <- fsttas$fstatistic
      pvalues[i] <- fsttas$pval
    }
    # e). M3 vs M5
    fsttas <- calculate_f_statistic(sse_reduced = sse[6],
                                     sse_full = sse[4], 
                                     df_reduced = df[6],
                                     df_full = df[4])
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
      fsttas <- calculate_f_statistic(sse_reduced = sse[i], 
                                       sse_full = sse[nn],
                                       df_reduced = df[i], 
                                       df_full = df[nn])
      fstatistic[i] <- fsttas$fstatistic
      pvalues[i] <- fsttas$pval
    }
  }
  
  list(
    pvalues = pvalues,
    fstatistic = fstatistic
  )
}

#' Function to calculate p-values for F-distribution based on Royston formula
#' 
#' @details
#' Uses formula on page 23 of Stata manual at 
#' https://www.stata.com/manuals/rfp.pdf. 
#' 
#' @param dev a vector of deviance for models i.e Null, linear, FP1,...FPm.
#' @param resid.df a vector of residual degrees of freedom for models.
#' @param n sample size/number of observations.
#' @param acd logical indicating use of acd transformation.
calculate_f_test_royston <- function(dev, 
                                     resid.df,
                                     n, 
                                     acd) {
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
      stats <- calculate_f_statistic_stata(
        dev_reduced = dev[i], dev_full = dev[5], d1 = df[i],
        d2 = resid.df[5], n = n
      )
      fstatistic[i] <- stats$fstatistic
      pvalues[i] <- stats$pval
      dev.diff[i] <- stats$dev.diff
    }
    # e). M3 vs M5
    stats2 <- calculate_f_statistic_stata(
      dev_reduced = dev[6], dev_full = dev[4], d1 = df[5],
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
      stats <- calculate_f_statistic_stata(
        dev_reduced = dev[i], dev_full = dev[nn], d1 = df[i],
        d2 = resid.df[nn], n = n
      )
      pvalues[i] <- stats$pval
      fstatistic[i] <- stats$fstatistic
      dev.diff[i] <- stats$dev.diff
    }
  }
  
  list(
    pvalues = pvalues, 
    dev.diff = dev.diff, 
    fstatistic = fstatistic
  )
}
