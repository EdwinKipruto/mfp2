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
#' @param show_titles Logical scalar. If \code{TRUE}, add automatically
#'   generated plot titles describing the FP/ACD and spike-at-zero
#'   representation. If \code{FALSE}, omit plot titles.
#' @param shape Numeric value. Shape of points used when residuals 
#' are displayed.
#' @param size_points Numeric value. Size of points used when residuals 
#' are displayed.
#' @param size_points_spike Numeric value. Size of the point drawn for the
#'   zero component of spike-at-zero covariates.
#' @param color_points Character value. Color of points used when residuals 
#'   are displayed.
#' @param color_line Character value. Color of the line representing 
#'   the partial predictor.
#' @param color_line_spike Character value. Color of the point drawn for the
#'   zero component of spike-at-zero covariates.
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
#' Spike-at-zero covariates are plotted according to the representation retained
#' in the final model:
#' * both components: the transformed positive-value FP/ACD component plus the
#'   structural-zero binary indicator.
#' * positive-part only: the transformed positive-value FP/ACD component, with
#'   the structural-zero binary indicator suppressed.
#' * zero indicator only: the structural-zero binary indicator without a
#'   continuous FP/ACD component.
#'
#' Plot titles use explicit spike-at-zero labels:
#' * `spike: both components FP(<powers>) + zero indicator` when both the
#'   positive-value FP component and the structural-zero indicator are retained.
#' * `spike: positive-part only FP(<powers>)` when only the positive-value FP
#'   component is retained.
#' * `spike: zero indicator only` when only the structural-zero indicator is
#'   retained.
#'
#' For ACD terms, the continuous component is labelled as
#' `ACD FP(<powers>)`, for example
#' `spike: both components ACD FP(0.5) + zero indicator`.
#'
#' Set \code{show_titles = FALSE} to suppress these automatically generated
#' plot titles. This is useful when plots are combined in multi-panel figures
#' or when custom titles are added later with ggplot2.
#'
#' For spike-at-zero covariates with a continuous component, the curve is drawn
#' only over `x > 0`. If a fitted value at `x = 0` is available, it is displayed
#' as a separate point with a vertical error bar, consistent with the convention
#' that the FP/ACD function is not evaluated at zero.
#'
#' See \code{predict.mfp2()} for details on partial predictors.
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
                     show_titles = TRUE,
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
  
  if (!is.logical(show_titles) || length(show_titles) != 1L || is.na(show_titles)) {
    stop("show_titles must be a single non-missing logical value.", call. = FALSE)
  }
  
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
    # Check if variable was assessed by the spike-at-zero algorithm.
    # Use stored spike metadata rather than catzero_list: for positive-part-only
    # spike terms, the binary indicator is intentionally suppressed, but the
    # variable is still a spike-at-zero covariate and must be plotted as a
    # positive-value function.
    is_spike <- FALSE
    if (!is.null(model$spike) && v %in% names(model$spike)) {
      is_spike <- isTRUE(model$spike[[v]])
    } else if (!is.null(model$fp_terms) &&
               v %in% rownames(model$fp_terms) &&
               "spike" %in% colnames(model$fp_terms)) {
      is_spike <- isTRUE(model$fp_terms[v, "spike"])
    }
    
    spike_dec_v <- saz_decision_codes[["continuous_only"]]
    if (is_spike &&
        !is.null(model$spike_dec) &&
        v %in% names(model$spike_dec) &&
        !is.na(model$spike_dec[[v]])) {
      spike_dec_v <- as.integer(model$spike_dec[[v]])
    }
    
    is_acd <- !is.null(model$fp_terms) &&
      v %in% rownames(model$fp_terms) &&
      "acd" %in% colnames(model$fp_terms) &&
      isTRUE(model$fp_terms[v, "acd"])
    
    power_label <- paste0(model$fp_powers[[v]], collapse = ", ")
    continuous_label <- if (is_acd) {
      sprintf("ACD FP(%s)", power_label)
    } else {
      sprintf("FP(%s)", power_label)
    }
    
    # Title with FP/ACD powers and explicit spike-at-zero representation.
    # The internal spike_dec code is not exposed in the plot title.
    plot_title <- if (is_spike) {
      saz_decision_label(
        spike_dec_v,
        style = "plot_title",
        continuous_label = continuous_label
      )
    } else {
      continuous_label
    }
    
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = .data$variable, y = .data$value)) +
      ggplot2::xlab(v) + ggplot2::ylab(ylab) +
      ggplot2::theme_bw()
    
    if (show_titles) {
      p <- p + ggplot2::ggtitle(plot_title)
    }
    
    # Residuals go first
    if (!partial_only) {
      p <- p + ggplot2::geom_point(data = pred_data[[v]],
                                   ggplot2::aes(y = .data$value + .data$resid),
                                   color = color_points,
                                   size = size_points,
                                   shape = shape)
    }
    
    # Then fitted line/ribbon on top.
    # For active spike-at-zero covariates, the positive-value FP/ACD curve must
    # not extend through x = 0. The retained representation is handled as one
    # of: both components, positive-part only, or zero-indicator only.
    if (is_spike && spike_dec_v == saz_decision_codes[["binary_only"]]) {
      # Binary-only SAZ:
      # predict.mfp2() represents the x-axis as the structural-zero indicator
      # I(x <= 0), so the prediction data may contain many repeated 0/1 rows
      # when terms_seq = "data". Collapse to one row per level before drawing;
      # otherwise geom_line() may connect duplicated values in data order.
      bin_df <- df[!duplicated(df$variable), , drop = FALSE]
      bin_df <- bin_df[order(bin_df$variable), , drop = FALSE]
      
      p <- p + ggplot2::geom_line(
        data = bin_df,
        ggplot2::aes(x = .data$variable, y = .data$value),
        linewidth = linewidth,
        linetype = linetype,
        color = color_line
      ) +
        ggplot2::geom_point(
          data = bin_df,
          ggplot2::aes(x = .data$variable, y = .data$value),
          color = color_line,
          size = size_points_spike
        ) +
        ggplot2::geom_errorbar(
          data = bin_df,
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2,
          color = color_line
        )
      
    } else if (is_spike && spike_dec_v == saz_decision_codes[["cont_binary"]]) {
      # Cont + binary: FP curve for x > 0, point at x = 0
      pos_df <- df[df$variable > 0, , drop = FALSE]
      pos_df <- pos_df[order(pos_df$variable), , drop = FALSE]
      zero_df <- df[df$variable == 0, , drop = FALSE]
      
      p <- p + ggplot2::geom_line(data = pos_df,
                                  ggplot2::aes(x = .data$variable, y = .data$value),
                                  linewidth = linewidth,
                                  linetype = linetype,
                                  color = color_line) +
        ggplot2::geom_ribbon(data = pos_df,
                             ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             alpha = alpha_fill, fill = color_fill)
      
      if (nrow(zero_df) > 0) {
        p <- p + ggplot2::geom_point(data = zero_df,
                                     ggplot2::aes(y = .data$value),
                                     color = color_line_spike,
                                     size = size_points_spike) +
          ggplot2::geom_errorbar(data = zero_df,
                                 ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                 width = 0.2,
                                 color = color_line)
      }
      
    } else if (is_spike && spike_dec_v == saz_decision_codes[["continuous_only"]]) {
      # Cont only: FP curve for x > 0, separate reference point at x = 0
      pos_df <- df[df$variable > 0, , drop = FALSE]
      pos_df <- pos_df[order(pos_df$variable), , drop = FALSE]
      zero_df <- df[df$variable == 0, , drop = FALSE]
      
      p <- p + ggplot2::geom_line(data = pos_df,
                                  ggplot2::aes(x = .data$variable, y = .data$value),
                                  linewidth = linewidth,
                                  linetype = linetype,
                                  color = color_line) +
        ggplot2::geom_ribbon(data = pos_df,
                             ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             alpha = alpha_fill, fill = color_fill)
      
      if (nrow(zero_df) > 0) {
        p <- p + ggplot2::geom_point(data = zero_df,
                                     ggplot2::aes(y = .data$value),
                                     color = color_line_spike,
                                     size = size_points_spike) +
          ggplot2::geom_errorbar(data = zero_df,
                                 ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                 width = 0.2,
                                 color = color_line)
      }
      
    } else {
      # Regular continuous plot (non-spike variable)
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

#' @describeIn fracplot S3 plot method for \code{mfp2} objects.
#'
#' @param x A fitted \code{mfp2} object.
#' @param ... Further arguments passed to \code{fracplot()}.
#'
#' @examples
#' data("prostate")
#' x <- as.matrix(prostate[, 2:8])
#' y <- as.numeric(prostate$lpsa)
#' fit <- mfp2(x, y, verbose = FALSE)
#' plot(fit)
#'
#' @method plot mfp2
#' @export
plot.mfp2 <- function(x, ...) {
  fracplot(model = x, ...)
}