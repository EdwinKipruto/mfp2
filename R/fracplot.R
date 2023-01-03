#' @title plot data and fit from fractional polynomial model
#'
#' @description plots the data and fit, with 95% confidence limits, from the fractional
#' polynomial (FP) model. The data and fit are plotted against the selected covariate of
#' interest.
#'

#' @param model mfpa model
#' @param x  variable name to be plotted.
#' @param plotype selects which graphics engine to use to generate the
#' fractional polynomial plot. Either 'ggplot' or 'rplot'. Default is "ggplot"
#' @param partials produces component-plus-residual plots of the continuous variables when set to TRUE otherwise
#' partial predictor would be produced. Default is TRUE
#' @param col.points The color of the data points when component-plus-residual plot is generated
#' @param col.line The color of the partial predictor line
#' @param title The title of the plot
#' @param size,cex line width of the partial predictor and the size of the points
#' when component-plus-residual plot is generated
#' @param ylab,xlab,pch other graphical parameters
#' @param fill the color for confidence interval
#' @param alpha refers to the opacity of the confidence interval.
#' Values of alpha range from 0 to 1, with lower values corresponding to more
#' transparent colors. See ggplot2

#' ## Details on `plotype` option
#'
#' Selects which graphics engine to use to generate the
#' fractional polynomial plot. Either 'ggplot' or 'rplot'. Default is "ggplot"
#'
#' For Cox models, the response should preferably be a \code{Surv} object,
#' created by the \code{Surv()} function in \pkg{survival} package and the \code{family = cox}. Only
#' right-censored data are currently supported. To fit stratified Cox
#' models, strata option can be used.


## @method plot mfpa
#' @import graphics
#' @import ggplot2
#' @export
## plot.mfpa
fracplot <- function(model, x, plotype = c("ggplot", "rplot"), partials = F, col.points = "gray",
                     title = NULL, ylim = NULL,
                     col.line = "red", pch = 16, size = 2, lty = 1, cex = 1, # lty.true = 5,
                     ylab = NULL, xlab = NULL, fill = "gray", alpha = 0.3, ...) {
  plotype <- match.arg(plotype)
  # if (!inherits(x, "mfpa"))
  #   stop("use only with \"mfpa\" objects")
  # Does a check if ggplot2 is available
  # It should be as it is in the imports section but in CRAN checks some systems don't have it!
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  # cox comes as character while glm is a function
  family <- ifelse(is.character(model$family), model$family, model$family$family)
  acdx <- model$acd
  # Assert that the model must be specified
  if (missing(model)) stop("Enter the model", call. = F)
  # Assert that x must be a character
  if (!is.character(x)) stop("x must be a character", call. = F)
  # Get rid of irrelevant variables denoted by NA in the FP powers
  px <- sapply(model$fp_powers, function(x) all(is.na(x)))
  # In acd we can have one FP power being NA and the other is not
  pwrs <- model$fp_powers[!px]
  #  Check whether x is in the final model.
  if (!(x %in% names(pwrs))) stop("The variable ", x, " not in model", call. = F)
  # strip the names of the variables by getting rid of dot extension
  varnamesx <- unlist(lapply(strsplit(c(names(model$coefficients)), "[.]"), `[[`, 1))
  # Get rid of "A" in acd variables such that FP2 variables or FP1 with acd are put together
  varnamesx <- unlist(strsplit(varnamesx, "[()]"))
  varnamesx <- varnamesx[varnamesx != ""]
  varnamesx <- varnamesx[varnamesx != "A"]
  # Rename coefficients such that FP2 have the same name
  coefx <- setNames(model$coefficients, varnamesx)
  # converts coefx to a list for easy identification of coefficients of x of interest
  sepcoef <- lapply(split(coefx, names(coefx)), unname)
  # select the intercept and coefficients of x of interest
  x.coefs <- unlist(sepcoef[if (family == "cox") {
    x
  } else {
    c("Intercept", x)
  }])
  nx <- length(sepcoef[[x]])
  # select transformed x values for x of interest. Note that names of transformed
  # variables are x.1, x.2 if FP2 and x.1 and (A)x.1 if acd. It can happen that
  # in acd one of the variables has FP power NA so its transformed values are missing
  kk <- if (acdx[x]) {
    paste0(c(x, paste0("(A)", x)), ".", c(1, 1))
  } else {
    paste(rep(x, nx), 1:nx, sep = ".")
  }
  k2 <- kk[which(!is.na(pwrs[[x]]))]
  X <- model$X[, k2, drop = F]
  # calculate partial predictor: beta0 + beta^T*X or beta^T*X in case of cox
  if (family == "cox") {
    partial <- X %*% x.coefs
  } else {
    partial <- cbind(1, X) %*% x.coefs
  }
  # standard error of the partial predictors
  stder <- calculate_standard_error(model = model, X = X)
  # lower and upper interval for partial
  pp <- qnorm(0.975)
  lowerci <- partial - pp * stder
  upperci <- partial + pp * stder
  # GLM:Deviance residuals while cox is martingale residuals
  residx <- if (family == "cox") {
    model$residuals
  } else {
    residuals.glm(model, type = "deviance")
  }
  # Component + residual.
  compresid <- residx + partial
  # Combine original x, partial, component, lowerci and upperci
  xx <- cbind(model$x[, x], partial, compresid, lowerci, upperci)
  colnames(xx) <- c(x, "partial", "component", "lower", "upper")
  # order xx
  xx <- xx[order(xx[, x]), ]
  # for title
  xpower <- unlist(model$fp_powers[x], use.names = F)

  # Add title
  # x label
  if (plotype == "rplot") {
    if (missing(xlab)) xlab <- x
    if (partials) {
      if (is.null(ylab)) ylab <- "Partial predictor"
      if (is.null(ylim)) ylim <- c(min(c(xx[, 2], xx[, 4])), max(c(xx[, 2], xx[, 5])))
      plot(
        x = xx[, x], y = xx[, "partial"], xlab = xlab, ylab = ylab, type = "l",
        col = col.line, lwd = size, ylim = ylim, ...
      )
    } else {
      if (is.null(ylab)) ylab <- "Component + residual"
      if (is.null(ylim)) ylim <- c(min(c(xx[, 3], xx[, 4])), max(c(xx[, 3], xx[, 5])))
      plot(
        x = xx[, x], y = xx[, "component"], xlab = xlab, pch = pch, ylab = ylab,
        type = "n", ylim = ylim, ...
      )
      # lines(x=xx[,x],y = xx[,"lower"], col = col.line, lwd = lwd, lty = lty)
      # lines(x=xx[,x],y = xx[,"upper"], col = col.line, lwd = lwd, lty = lty)
      polygon(c(xx[, x], rev(x = xx[, x])), c(xx[, "lower"], rev(xx[, "upper"])), col = fill, border = NA)
      points(x = xx[, x], y = xx[, "component"], cex = cex, pch = pch, col = col.points, ...)
      lines(x = xx[, x], y = xx[, "partial"], col = col.line, lwd = size, lty = lty)
    }
    if (is.null(title)) {
      # title(paste0("FP(",paste0(xpower, collapse = ","), ")"),line = 0.5)
      title(if (acdx[x]) {
        paste0("FP1(", paste0(xpower, collapse = ","), ")")
      } else {
        paste0("FP(", paste0(xpower, collapse = ","), ")")
      }, line = 0.5)
    }
  } else {
    xx <- data.frame(xx)
    nn <- colnames(xx)
    if (partials) {
      gp <- ggplot(xx, aes_string(x = nn[1], y = nn[3])) +
        geom_line(aes_string(x = nn[1], y = nn[2]), size = size, color = col.line) +
        geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = alpha, fill = fill)
    } else {
      gp <- ggplot(xx, aes_string(x = nn[1], y = nn[3])) +
        geom_point(aes_string(x = nn[1], y = nn[3]), size = cex, color = col.points, shape = pch) +
        geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = alpha, fill = fill) +
        geom_line(aes_string(x = nn[1], y = nn[2]), size = size, color = col.line)
    }
    tt <- if (acdx[x]) {
      paste0("FP1(", paste0(xpower, collapse = ","), ")")
    } else {
      paste0("FP(", paste0(xpower, collapse = ","), ")")
    }
    if (is.null(title)) gp <- gp + ggtitle(tt)
    gp
  }
}

# alias to keep things from breaking
plot_mfp <- function(...) {
    fracplot(...)
}
