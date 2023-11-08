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
#' use a Chi-square distribution with df_full - df_reduced degrees
#' of freedom to derive a p-value.
#' 
#' This is basically the same way as [stats::anova()] implements the 
#' likelihood ratio test.
#' 
#' @return
#' A list with two entries for the likelihood ratio test for the ratio 
#' `logl[1] / logl[2]`.
#' 
#' * `statistic`: test statistic.
#' * `pvalue`: p-value
calculate_lr_test <- function(logl, 
                              dfs) {
  statistic <- 2 * (logl[2] - logl[1])
  
  list(
    statistic = statistic, 
    pvalue = pchisq(statistic, 
                    df = dfs[2] - dfs[1], 
                    lower.tail = FALSE)
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
#' @param d1 a numeric value giving `d1` in the formula below directly as 
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
#' A list with three entries giving the test statistic and p-value for the F-test 
#' for the comparison of `deviance[1]` to `deviance[2]`.
#' 
#' * `statistic`: test statistic.
#' * `pvalue`: p-value. 
#' * `dev_diff`: difference in deviances tested.
#' 
#' @references 
#' Kutner, M.H., et al., 2004. \emph{Applied linear statistical models. 
#' McGraw-Hill Irwin.}
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
  
  pvalue <- 1
  if (dev_diff != 0) {
    pvalue <- stats::pf(statistic, df1 = d1, df2 = d2, lower.tail = FALSE)
  }
  
  list(
    statistic = statistic, 
    pvalue = pvalue, 
    dev_diff = dev_diff
  )
}
