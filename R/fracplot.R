#' Plot response functions from a fitted `mfp2` object
#'
#' Produces partial predictor plots (or contrasts) with confidence intervals
#' against selected covariates from a fitted \code{`mfp2`} model. If requested, 
#' component-plus-residual plots are also supported. 
#'
#' @param model A fitted \code{`mfp2`} model.
#' @param terms Character vector with variable names to be plotted.  If `NULL`,
#' all fractional polynomial terms in the model are plotted.
#' @param partial_only Logical. If `TRUE`, only the partial predictor (model
#' component) is plotted. If `FALSE` (default), component-plus-residual plots
#' are drawn. Only used if `type = "terms"`. See below for details. 
#' @param type Character, one of `"terms"` or `"contrasts"`. Passed to
#' \code{predict.mfp2()}. `"terms"` plots partial predictors (with or without
#' residuals), while `"contrasts"` plots contrasts relative to a reference
#' value. 
#' @param ref Reference list passed to \code{predict.mfp2()} when `type =
#'  "contrasts"`. Ignored otherwise.
#' @param terms_seq Character, one of `"data"` or `"equidistant"`. Passed to
#'  \code{predict.mfp2()}. `"data"` uses the observed values, `"equidistant"` creates
#' a grid over the covariate range.
#' @param alpha Confidence level for intervals. Passed to \code{predict.mfp2()}.
#' @param shape Numeric value. Shape of points used when residuals 
#' are displayed.
#' @param size_points Numeric value. Size of points used when residuals 
#' are displayed.
#' @param size_points_spike Numeric value. Size of the point drawn at zero 
#' when the covariate includes a spike-at-zero component (`spike_decision = 1`).
#' @param color_points Character value. Color of points used when residuals 
#'   are displayed.
#' @param color_line Character value. Color of the line representing 
#'   the partial predictor.
#' @param color_line_spike Character value. Color of the point drawn at zero 
#'   when the covariate includes a spike-at-zero component (`spike_decision = 1`).
#' @param linetype Character value. Line type for the partial predictor. 
#'   See [ggplot2::geom_line()] for options.
#' @param linewidth Numeric value. Width of the line representing 
#' the partial predictor.
#' @param color_fill Character value. Fill color of the confidence interval ribbon.
#' @param alpha_fill Numeric value between 0 and 1. Transparency of the 
#' confidence interval ribbon.
#' @param ... Further arguments passed from the alias `plot_mfp()`. 
#' 
#' @details 
#' Confidence intervals are based on the variance–covariance matrix of the 
#' final fitted model. They reflect uncertainty in the regression coefficients 
#' but not in the selection of fractional polynomial powers. Intervals may 
#' therefore be too narrow. A bootstrap approach (not yet implemented) is 
#' recommended for more realistic intervals (see Royston & Sauerbrei, 2008, 
#' Section 4.9.2).
#' 
#' Component-plus-residual plots are available if `type = "terms"`. Deviance 
#' residuals are used for generalized linear models, while martingale residuals 
#' are used for Cox regression. This matches the behavior of the Stata `mfp` 
#' program.
#' 
#' Spike-at-zero covariates are handled according to the `spike_decision` code:
#' * `1` – Include both the transformed FP function for positive values and the binary 
#'       spike-at-zero indicator.
#' * `2` – Ignore the spike; treat the variable as continuous (usual FP plot).
#' * `3` – Show only the binary spike-at-zero indicator.
#'
#' Plot behavior for each decision:
#' * If `spike_decision == 1`, the plot shows the FP function for positive values 
#'   and includes the binary spike-at-zero indicator. The term 
#'   \eqn{\hat{\beta}_0 + \hat{\beta}} for observations equal to zero is also 
#'   displayed with a vertical error bar. The plot title includes 
#'   `+ z` to indicate the presence of the spike-at-zero component. The FP power 
#'   for the positive part is enclosed in parentheses. For example, `FP(0) + z` 
#'   indicates an FP power of 0 (log) for the positive values.
#' * If `spike_decision == 3`, the plot shows the binary indicator alone (`z only` in 
#'   the title). Mean values at 0 and 1 are connected with a line, and a ribbon showing 
#'   confidence intervals is displayed.
#' * If `spike_decision == 2` (or not specified), the covariate is plotted as a 
#'   continuous FP function in the usual way.
#'  See \code{fracplot} for details on partial predictors
#' @examples
#'
#' # Gaussian response
#' data("prostate")
#' x = as.matrix(prostate[,2:8])
#' y = as.numeric(prostate$lpsa)
#' # default interface
#' fit = mfp2(x, y, verbose = FALSE)
#' fracplot(fit) # generate plots
#'
#' @return 
#' A list of `ggplot2` plot objects, one for each term requested. Can be 
#' drawn as individual plots or facetted / combined easily using e.g. 
#' `patchwork::wrap_plots` and further customized. 
#' 
#' @seealso 
#' \code{predict.mfp2()}
#' 
#' @import ggplot2
#' @importFrom ggplot2 .data
#' @export
fracplot <- function(model, 
                     terms = NULL, 
                     partial_only = FALSE, 
                     type = c("terms","contrasts"),
                     ref = NULL,
                     terms_seq = c("data", "equidistant"),
                     alpha = 0.05,
                     color_points = "#AAAAAA",
                     color_line = "red", 
                     color_line_spike = "red",
                     color_fill = "#000000",
                     shape = 1,
                     size_points = 1,
                     size_points_spike = 2,
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
  # assert that the object must be mfp2
  if (!inherits(model, "mfp2")) 
    stop("The model entered is not an mfp2 object.", call. = FALSE)
  
 # check if ref is a list (only if it is not NULL)
  if (!is.null(ref) && !is.list(ref))
    stop("ref must be a list", call. = FALSE)
  
 # set defaults depending on type
  type <- match.arg(type)
  if (type == "contrasts") {
    partial_only <- TRUE
    if (missing(terms_seq))
      terms_seq <- "equidistant"
  }
  
  terms_seq <- match.arg(terms_seq)
  
  pred <- predict(model, 
                  type = type, 
                  terms = terms,
                  ref = ref,
                  terms_seq = terms_seq,
                  alpha = alpha)
  # for points also need the data point predictions
  pred_data <- pred
  if (!partial_only && terms_seq != "data") {
    pred_data <- predict(model,
                         type = "terms", 
                         terms = terms, 
                         terms_seq = "data")
    }
  
  if (length(pred) == 0) {
    warning("i Variables specified in terms not used in final model.")
    return()
  }
  
  ylab <- "Partial Predictor"
  
  if (!partial_only) {
    # compute residuals to plot points
    # for glm, deviance residuals are required
    # while for cox martingale residuals
    resid <- if (model$family_string == "cox") {
      model$residuals
    } else {
      residuals.glm(model, type = "deviance")
    }
    # add residuals to the data
    #pred_data <- lapply(pred_data, function(v) transform(v, resid = resid))
    pred_data <- lapply(pred_data, function(v) {
      v$resid <- resid
    v})
    
    # y label for the plot
    ylab <- "Partial Predictor + residuals"
  }
  
  # Preallocate the list for plots
  plots <- setNames(vector("list", length(names(pred))), names(pred))  
  
  # in the calls to ggplot2::aes use the .data pronoun to avoid notes 
  # generated by R CMD CHECK about missing bindings of global variables
  for (v in names(pred)) {
    df <- pred[[v]]
    # Check if variable has a spike-at-zero binary indicator
    is_spike <- !is.null(model$catzero_list[[v]])
    
    # Title with FP powers and spike label
    title_spike <- ifelse(model$spike_decision[v] == 1, " + z",
                          ifelse(model$spike_decision[v] == 3, " (z only)", ""))
    
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = .data$variable, y = .data$value)) +
      ggplot2::ggtitle(sprintf("FP%s(%s)%s%s", 
                               ifelse(model$fp_terms[v, "acd"], "1", ""),
                               paste0(model$fp_powers[[v]], collapse = ", "),
                               ifelse(model$fp_terms[v, "acd"], ":ACD", ""),
                               title_spike)
      ) +
      ggplot2::xlab(v) + ggplot2::ylab(ylab) +
      ggplot2::theme_bw()
    
    # Residuals go first
    if (!partial_only) {
      p <- p + ggplot2::geom_point(data = pred_data[[v]],
                                   ggplot2::aes(y = .data$value + .data$resid),
                                   color = color_points,
                                   size = size_points,
                                   shape = shape)
    }
    
    # Then fitted line/ribbon on top
    # Add line for positive values only if spike-at-zero exists
    if (is_spike && model$spike_decision[v] != 2) {
      pos_df <- df[df$variable > 0, , drop = FALSE]
      zero_df <- df[df$variable == 0, , drop = FALSE]
      
      if (model$spike_decision[v] == 3) {
        # Expect df to have two rows: one for 0 and one for 1
        p <- p + ggplot2::geom_line(linewidth = linewidth,
                                    linetype = linetype,
                                    color = color_line) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                               alpha = alpha_fill,
                               fill = color_fill)
      } else {
      p <- p + ggplot2::geom_line(data = pos_df, 
                                  ggplot2::aes(x = .data$variable, y = .data$value),
                                  linewidth = linewidth,
                                  linetype = linetype,
                                  color = color_line) +
        ggplot2::geom_ribbon(data = pos_df,
                             ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             alpha = alpha_fill, fill = color_fill)
      # Add point at zero for the spike
      zero_df <- df[df$variable == 0, , drop = FALSE]
      if (nrow(zero_df) > 0) {
        p <- p + ggplot2::geom_point(data = zero_df,
                                     ggplot2::aes(y = .data$value),
                                     color = color_line_spike,
                                     size = size_points_spike) + 
          ggplot2::geom_errorbar(data = zero_df,
                                 ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                 width = 0.2,   # small width so it's vertical
                                 color = color_line)
      }
      }
      
    } else {
      # regular continuous plot
      p <- p + ggplot2::geom_line(linewidth = linewidth,
                                  linetype = linetype,
                                  color = color_line) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             alpha = alpha_fill, fill = color_fill)
    }
    
    # Add residual points if needed
    # if (!partial_only) {
    #   p <- p + ggplot2::geom_point(data = pred_data[[v]],
    #                                ggplot2::aes(y = .data$value + .data$resid),
    #                                color = color_points,
    #                                size = size_points,
    #                                shape = shape)
    # }
    plots[[v]] <- p
  }

  plots
}

#' @describeIn fracplot Alias for fracplot.
plot_mfp <- function(...) {
    fracplot(...)
}
