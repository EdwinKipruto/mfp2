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
calculate_lr_test <- function(logl, 
                              dfs) {
  list(
    statistic = 2 * (logl[2] - logl[1]), 
    pvalue = pchisq(2 * (logl[2] - logl[1]), 
                    df = dfs[2] - dfs[1], 
                    lower.tail = FALSE)
  )
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

#' Function to compute F-statistic and p-value from deviances
#' 
#' Alternative to likelihood ratio tests in normal / Gaussian error models. 
#' 
#' @param deviances a numeric vector of length 2 with deviances. Typically 
#' ordered in increasing order (i.e. null model first, then full model) and 
#' used to test the difference `deviances[1] - deviances[2]`.
#' @param dfs_resid a numeric vector with residual degrees of freedom.
#' @param n_obs a numeric value with the number of observations.
#' @param d1 a numeriv value giving `d1` in the formula below directly as 
#' the number of additional degrees of freedom in model 2 compared to model 1.
#' In this case `dfs_resid` must be a single numeric value giving the residual
#' df for model 2. This interface is sometimes more convenient than to specify
#' both residual dfs. 
#' 
#' @details 
#' Uses formula on page 23 from here: https://www.stata.com/manuals/rfp.pdf:
#' \deqn{F = \frac{d_2}{d_1} (exp(\frac{D_2 - D_1}{n}) - 1),}
#' where \eqn{D} refers to deviances of two models 1 and 2. 
#' \eqn{d1} is the number of additional parameters used in in model 2 as 
#' compared to model 1, i.e. `dfs_resid[1] - dfs_resid[2]`. 
#' \eqn{d2} is the number of residual degrees of freedom minus the number of
#' estimated powers for model 2, i.e. `dfs_resid[2]`.
#' #' The p-value then results from the use of a F-distribution with 
#' (d1, d2) degrees of freedom.
#' 
#' Note that this computation is completely equivalent to the computation
#' of a F-test using sum of squared errors as in e.g. Kutner at al. (2004), 
#' p 263. The formula there is given as 
#' \deqn{F = \frac{SSE(R) - SSE(F)}{df_R - df_F} / \frac{SSE(F)}{df_F},}
#' where the \eqn{df} terms refer to residual degrees of freedom, and \eqn{R}
#' and \eqn{F} to the reduced (model 1) and full model (model 2), respectively.
#' 
#' @return
#' The p-value for the F-test for the comparison of `deviance[1]` to 
#' `deviance[2]`.
#' 
#' @references 
#' Kutner, M.H., et al., 2004. \emph{Applied linear statistical models. 
#' McGraw-Hill Irwin.}
#' 
#' @seealso 
#' [calculate_f_test_royston()]
calculate_f_test <- function(deviances,
                             dfs_resid,
                             n_obs, 
                             d1 = NULL) {
  
  dev_diff <- deviances[1] - deviances[2]
  
  # number of additional parameters in model 2 compared to model 1
  if (is.null(d1)) {
    # since we use residual dfs, we subtract in changed order
    d1 <- dfs_resid[1] - dfs_resid[2]
    d2 <- dfs_resid[2]
  } else {
    d2 <- dfs_resid
  }
  
  statistic <- (d2 / d1) * (exp(dev_diff / n_obs) - 1)
  
  p_value <- 1
  if (dev_diff != 0) {
    p_value <- stats::pf(statistic, df1 = d1, df2 = d2, lower.tail = FALSE)
  }
  
  list(
    statistic = statistic, 
    p_value = p_value, 
    dev_diff = dev_diff
  )
}

#' Function to calculate p-values for F-distribution based on Royston formula
#' 
#' Alternative to [calculate_f_test()].
#' 
#' @details
#' Uses formula on page 23 of Stata manual at 
#' https://www.stata.com/manuals/rfp.pdf. 
#' 
#' @param dev a vector of deviance for models i.e Null, linear, FP1,...FPm.
#' @param resid.df a vector of residual degrees of freedom for models.
#' @param n sample size/number of observations.
#' @param acd logical indicating use of acd transformation.
#' 
#' @seealso 
#' [calculate_f_test()]
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
      stats <- calculate_f_test(
        deviances = c(dev[i], dev[5]), 
        d1 = df[i], dfs_resid = resid.df[5], n_obs = n
      )
      fstatistic[i] <- stats$statistic
      pvalues[i] <- stats$p_value
      dev.diff[i] <- stats$dev_diff
    }
    # e). M3 vs M5
    stats2 <- calculate_f_test(
      deviances = c(dev[6], dev[4]),
      d1 = df[5], dfs_resid = resid.df[4], n_obs = n
    )
    fstatistic[5] <- stats2$statistic
    pvalues[5] <- stats2$p_value
    dev.diff[5] <- stats2$dev_diff
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
      stats <- calculate_f_test(
        deviances = c(dev[i], dev[nn]), 
        d1 = df[i], dfs_resid = resid.df[nn], n_obs = n
      )
      pvalues[i] <- stats$p_value
      fstatistic[i] <- stats$statistic
      dev.diff[i] <- stats$dev_diff
    }
  }
  
  list(
    pvalues = pvalues, 
    dev.diff = dev.diff, 
    fstatistic = fstatistic
  )
}
