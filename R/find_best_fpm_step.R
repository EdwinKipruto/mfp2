# similar with find_best_fp1_step() but for degree>1
# degree =  degree of fp. m = 2 is fp2, m = 3 is fp3 etc
find_best_fpm_step <- function(y, x, xi, allpowers, powers, family, weights, offset, strata,
                    control, method, rownames, nocenter, degree, acdx) {
  #
  if (degree < 2) stop("Degree must be >= 2. Here we are interest on best FPm where
                    m>=2. For m = 1 see find_best_fp1_step()")
  # Generate FP data for x of interest (xi) and adjustment variables
  m <- degree
  df1 <- extract_adjustment_data(
    x = x, xi = xi, allpowers = allpowers,
    df = 2 * m, powers = powers, acdx = acdx
  )
  # Matrix of adjustment data
  adjdata <- df1$adjdata
  # List of FP data for xi (continuous) of interest.
  fpdata <- df1$fpdata
  # Total FP powers different from 1 estimated in adjustment model
  # tFP <- calculate_number_fp_powers(df1$adjustpowers)
  # The length of generated FP variables for x of interest.
  nv <- length(fpdata)
  # Generate all possible FPm powers
  fpmpowers <- generate_fp_powers(degree = m, powers = powers)
  # log(n) for bic calculation
  n <- nrow(x)
  logn <- log(n)
  devs <- devs.royston <- sse <- aic <- bic <- dfpm <- dfx <- numeric(nv)
  xnames <- paste0("newx", seq_along(1:m))
  for (i in seq_len(nv)) {
    # combine FP variables for x of interest with adjustment variables
    xout <- cbind(fpdata[[i]], adjdata)
    colnames(xout)[1:m] <- xnames
    fit1 <- fit_model(
      x = xout, y = y, family = family, method = method,
      weights = weights, offset = offset, strata = strata,
      control = control, rownames = rownames,
      nocenter = nocenter
    )
    # Deviance of the fitted model
    devs[i] <- -2 * fit1$logl
    # number of regression coefficients (plus variance if gaussian)
    dfx[i] <- fit1$df
    # save degrees of freedom of sse for fpm i.e n-#parameters in fpm model
    # useful for F test. we add m because of m fp powers and tFP because of
    # estimated FP powers in adjustment model
    # dfpm[i] = n-(fit1$df + m + tFP)
    dfpm[i] <- n - ((fit1$df - 1) + m) # subtract scale parameter
    # AIC and BIC of the fitted model. Add m because of the m FPm powers.
    # aic[i] = devs[i] + 2*(fit1$df + m + tFP) #
    # bic[i] = devs[i] + logn*(fit1$df + m + tFP)
    aic[i] <- devs[i] + 2 * (fit1$df + m) #
    bic[i] <- devs[i] + logn * (fit1$df + m)
    # sse and deviance for gaussian family.
    sse[i] <- fit1$SSE
    devs.royston[i] <- deviance_gaussian(RSS = sse[i], weights = weights, n = dim(x)[1L])
  }
  # Best FPm function based on dev (calculated using loglik), aic, bic, sse and dev(calculated using SSE)
  fn.bestfpm <- list(
    dev = fpmpowers[which.min(devs), ],
    aic = fpmpowers[which.min(aic), ],
    bic = fpmpowers[which.min(bic), ],
    sse = fpmpowers[which.min(sse), ],
    dev.r = fpmpowers[which.min(devs.royston), ]
  )

  mt <- list(
    dev.best.fpm = devs[which.min(devs)], # Deviance of the best FPm function
    aic.best.fpm = aic[which.min(aic)], # aic of the best FPm function
    bic.best.fpm = bic[which.min(bic)],
    sse.best.fpm = sse[which.min(sse)],
    dev.all.FPm = devs,
    aic.all.FPm = aic,
    bic.all.FPm = bic,
    sse.all.FPm = sse,
    df.best.fpm.sse = dfpm[which.min(sse)],
    fn.bestfpm = fn.bestfpm,
    dev.roy.best.fpm = devs.royston[which.min(devs.royston)]
  ) # Deviance of the best FPm function based on royston formula

  return(mt)
}
