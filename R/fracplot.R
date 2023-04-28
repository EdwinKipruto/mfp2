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
fracplot <- function(model, 
                     terms = NULL, 
                     partial_only = FALSE, 
                     terms_seq = "data",
                     col.points = "gray",
                     title = NULL, 
                     ylim = NULL,
                     col.line = "red", 
                     pch = 16, 
                     size = 2, 
                     lty = 1, 
                     cex = 1, 
                     ylab = NULL, 
                     xlab = NULL, 
                     fill = "gray",
                     alpha = 0.3, 
                     ...) {
  
  # Does a check if ggplot2 is available
  # It should be as it is in the imports section but in CRAN checks some systems don't have it!
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  pred <- predict(model, type = "terms", terms = terms, terms_seq = terms_seq)
  
  # for points also need the data point predictions
  pred_data <- pred
  if (!partial_only && terms_seq != "data")
    pred_data <- predict(model, type = "terms", 
                         terms = terms, 
                         terms_seq = "data")
  
  if (length(pred) == 0) {
    warning("i Variables specified in terms not used in final model.")
    return()
  }
  
  if (!partial_only) {
    # compute residuals to plot points
    # for glm, deviance residuals are required
    # while for cox martingale residuals
    resid <- if (model$family_string == "cox") {
      model$residuals
    } else {
      residuals.glm(model, type = "deviance")
    }
  }
  
  plots <- list()
  for (v in names(pred)) {
    plots[[v]] <- ggplot2::ggplot(data = pred[[v]], 
                         aes(x = variable, y = contrast)) + 
      ggplot2::geom_line() + 
      ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper), 
                  alpha = 0.1) + 
      ggplot2::ggtitle(sprintf("FP(%s)%s", 
                      paste0(model$fp_powers[[v]], collapse = ", "),
                      ifelse(model$fp_terms[v, "acd"], ", ACD", "")
                      )) +
      ggplot2::xlab(v) +
      ggplot2::theme_bw()
    
    if (!partial_only) 
      plots[[v]] <- plots[[v]] + 
        ggplot2::geom_point(data = pred_data[[v]], aes(y = contrast + resid))
  }
  
  plots
}

# alias to keep things from breaking
plot_mfp <- function(...) {
    fracplot(...)
}
