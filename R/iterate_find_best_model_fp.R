# A function that estimates the best FP functions for all predictors (x) in the
# ith cycle.
iterate_find_best_model_fp <- function(x, y, allpowers, df, weights, offset, family, criterion,
                      select, alpha, keep, powers, method, strata, verbose, ftest, control,
                      rownames, nocenter, acdx) {
  # FP powers have the same names as variables in design matrix x, but the former is ordered.
  varnames <- names(allpowers)
  # For printing on the screen: fpgen distinguishes between p-value and AIC/BIC output
  if (verbose) {
      print_mfp_summary(criterion, ftest = ftest)
  }
  # Initialization, starting from the first variable
  i <- 1
  # Evaluate the variables according to how order_variables() ordered them. If they are
  # ascending, they are evaluated from most to least significant.
  while (i <= ncol(x)) {
    # Estimate xi's best FP power. It can also be NA (variable not significant),
    # linear, FP1, FP2, and so on.
    fitx <- find_best_model_fp(
      x = x, y = y,
      xi = varnames[i],
      allpowers = allpowers,
      weights = weights,
      offset = offset,
      df = df[i], # df for the ith variable
      select = select[i], # select for the ith variable
      alpha = alpha[i],
      keep = keep,
      family = family,
      criterion = criterion,
      powers = powers,
      method = method,
      strata = strata,
      ftest = ftest,
      control = control,
      rownames = rownames,
      nocenter = nocenter,
      acdx = acdx,
      verbose = verbose
    )
    # update the FP power of the ith variable
    allpowers[[i]] <- fitx$overall.best.fn
    # Move to the next variable
    i <- i + 1
  }
  # we will return the deviances, aic, bic etc. in the future
  return(allpowers)
}
