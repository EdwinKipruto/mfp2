# If the maximum allowed degree is 5, this function will calculate metrics for
# FP2 through FP5, which will then be combined with metrics for FP1 estimated by
# the find_best_fp1() function. We should have started with FP1 and worked our way up
# to FPm, but the linear function makes the formula more complicated, so there
# are separate functions for FP1 and FPm.
calculate_metrics_fpm <- function(y, x, xi, allpowers, powers, family, weights, offset, strata,
                        control, method, rownames, nocenter, degree, acdx) {
  # we output best fpm parameters from degree 2 to m for degree 1 see bestfp1
  mm <- seq(2, degree)
  out <- vector(mode = "list", length = length(mm))
  for (k in seq_along(mm)) {
    out[[k]] <- find_best_fpm(
      y = y, x = x, xi = xi, allpowers = allpowers,
      powers = powers, family = family, weights = weights,
      offset = offset, strata = strata, control = control,
      method = method, rownames = rownames,
      nocenter = nocenter, degree = mm[k], acdx = acdx
    )
  }
  # save AIC, BIC, DEV and SSE for m = 2,3,...
  dev <- unlist(lapply(out, `[[`, 1), use.names = F) # deviances are in position 1 of nested list
  aic <- unlist(lapply(out, `[[`, 2), use.names = F) # aic are in position 2 of nested list
  bic <- unlist(lapply(out, `[[`, 3), use.names = F) # bic are in position 3 of nested list
  sse <- unlist(lapply(out, `[[`, 4), use.names = F) # sse are in position 4 of nested list
  dev.roy <- unlist(lapply(out, `[[`, 11), use.names = F) # royston deviances are in position 11 of nested list

  # degree of freedom: regression coefficients plus estimated FP powers
  df.best.fpm.sse <- unlist(lapply(out, `[[`, 9), use.names = F) # df are in position 9 of nested list
  names(dev) <- names(dev.roy) <- names(aic) <- names(bic) <- names(sse) <- names(df.best.fpm.sse) <- paste0("FP", mm)
  # save corresponding fp powers
  fn.bestfpm <- lapply(out, `[[`, 10) # functions are in position 10 of nested list
  fun.dev <- lapply(fn.bestfpm, `[[`, 1)
  fun.aic <- lapply(fn.bestfpm, `[[`, 2)
  fun.bic <- lapply(fn.bestfpm, `[[`, 3)
  fun.sse <- lapply(fn.bestfpm, `[[`, 4)
  fun.dev.roy <- lapply(fn.bestfpm, `[[`, 5)

  return(list(
    dev = dev, aic = aic, bic = bic, sse = sse, df.best.fpm.sse = df.best.fpm.sse,
    fun = list(dev = fun.dev, aic = fun.aic, bic = fun.bic, sse = fun.sse, dev.roy = fun.dev.roy),
    dev.roy = dev.roy
  ))
}
