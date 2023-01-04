# Fits univariable Fractional Polynomial models. No testing is conducted
# x is a vector or matrix (with one column) of predictor
# y is the outcome variable
# returns the best FP power with smallest deviance
fit_fp <- function(y, x, powers = NULL, family, method, weights, shift = NULL,
                     scalex = NULL, degree, offset, strata, control, rownames,
                     nocenter) {
  # Check whether x has more than one variables
  x <- as.matrix(x)
  if (dim(x)[2L] > 1) stop("Only one predictor is allowed")
  # estimate shifting factor from data when not supplied
  if (is.null(shift)) {
    shift <- find_shift_factor(x)
  }
  # shift x to positive values if any(x)<= 0
  x <- x + shift
  # scale x for computational stability
  if (is.null(scalex)) {
    scalex <- find_scale_factor(x)
  }
  x <- x / scalex
  # Set default powers when missing
  if (is.null(powers)) {
    powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  }
  nv <- length(powers)
  # Generate FP data
  dt <- generate_fp_transformations(x = x[, 1, drop = T], degree = degree, powers = powers)$data
  devs <- fits <- vector(mode = "list", length = nv)
  # Fit linear models for each FP1 function
  xnames <- paste0("newx", seq_along(1:degree))
  for (i in seq_len(nv)) {
    xout <- dt[[i]]
    colnames(xout)[1:degree] <- xnames
    # Fit the  model which can be glm or cox depending on the family chosen
    fit1 <- fit_model(
      x = xout, y = y, family = family, method = method,
      weights = weights, offset = offset, strata = strata,
      control = control, rownames = rownames,
      nocenter = nocenter
    )
    # Deviance of the fitted model
    devx <- -2 * fit1$logl
    devs[i] <- devx
    # save also each model
    fits[[i]] <- fit1$fit
  }
  # Index of a model with the smallest deviance
  index <- which.min(devs)
  # The  FP power that yielded the smallest deviance
  s <- powers[[index]]
  # The model that produced the smallest deviance
  bestmodel <- fits[[index]]
  return(list(bestpower = s, bestmodel = bestmodel))
}
