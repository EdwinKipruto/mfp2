#' Plot response functions from a fitted `mfpa` object
#'
#' Plots partial linear predictors against the selected covariate(s) of 
#' interest, including a confidence interval.
#'
#' @param model fitted `mfpa` model.
#' @param terms character vector with variable names to be plotted.
#' @param partial_only a logical value indicating whether only the partial 
#' predictor is drawn (`TRUE`), or also observed data points (`FALSE`, the 
#' default). See below for details. 
#' @param terms_seq `terms_seq` argument of [predict.mfpa()].
#' @param alpha `alpha` argument of [predict.mfpa()].
#' @param shape,size_points,color_points `ggplot2` properties of drawn 
#' data points.
#' @param color_line,linetype,linewidth `ggplot2` properties of line for 
#' partial predictor.
#' @param color_fill,alpha_fill `ggplot2` properties of ribbon for confidence
#' interval.
#' 
#' @details 
#' The data points are computed as the partial linear predictors plus residuals
#' extracted from the fitted model (deviance residuals for glms, and martingale
#' residuals for Cox models), following the mfp implementation in Stata.
#' 
#' @return 
#' A list of `ggplot2` plot objects, one for each term requested. Can be 
#' drawn as individual plots of facetted / combined easily using e.g. 
#' `patchwork::wrap_plots` and further customized. 
#' 
#' @seealso 
#' [predict.mfpa()]
#' 
#' @import ggplot2
#' @export
fracplot <- function(model, 
                     terms = NULL, 
                     partial_only = FALSE, 
                     terms_seq = "data",
                     alpha = 0.05,
                     color_points = "#AAAAAA",
                     color_line = "#000000", 
                     color_fill = "#000000",
                     shape = 1,
                     size_points = 1,
                     linetype = "solid", 
                     linewidth = 1,
                     alpha_fill = 0.1) {
  
  # Does a check if ggplot2 is available
  # It should be as it is in the imports section but in CRAN checks some systems don't have it!
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  pred <- predict(model, 
                  type = "terms", 
                  terms = terms, 
                  terms_seq = terms_seq,
                  alpha = alpha)
  
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
                         aes(x = variable, y = value)) + 
      ggplot2::geom_line(linewidth = linewidth,
                         linetype = linetype,
                         color = color_line) + 
      ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper), 
                           alpha = alpha_fill) + 
      ggplot2::ggtitle(sprintf("FP(%s)%s", 
                      paste0(model$fp_powers[[v]], collapse = ", "),
                      ifelse(model$fp_terms[v, "acd"], ", ACD", "")
                      )) +
      ggplot2::xlab(v) +
      ggplot2::theme_bw()
    
    if (!partial_only) 
      plots[[v]] <- plots[[v]] + 
        ggplot2::geom_point(data = pred_data[[v]], aes(y = value + resid),
                            color = color_points,
                            size = size_points, 
                            shape = shape)
  }
  
  plots
}

#' @describeIn fracplot Alias for fracplot.
plot_mfp <- function(...) {
    fracplot(...)
}
