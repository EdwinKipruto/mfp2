#' Plot response functions from a fitted `mfp2` object
#'
#' Produces partial-predictor or contrast plots with confidence intervals for
#' selected covariates from a fitted \code{mfp2} model. For term plots, the
#' fitted component can be shown alone or combined with residuals.
#'
#' @param x A fitted \code{mfp2} object.
#' @param terms Character vector naming variables to plot. If \code{NULL}, all
#'   eligible terms retained in the final model are plotted.
#' @param partial_only Logical scalar. If \code{TRUE}, plot only the fitted
#'   partial predictor. If \code{FALSE}, the default, add residuals to produce
#'   component-plus-residual plots. This argument is used only when
#'   \code{type = "terms"}; contrast plots always use \code{partial_only = TRUE}.
#' @param type Character scalar, either \code{"terms"} or \code{"contrasts"}.
#'   Term plots display partial predictors; contrast plots display differences
#'   from the reference supplied through \code{ref}.
#' @param ref Optional named list of reference values used when
#'   \code{type = "contrasts"}. Ignored for term plots.
#' @param terms_seq Character scalar, either \code{"data"} or
#'   \code{"equidistant"}. The first evaluates the fitted function at observed
#'   covariate values; the second uses an equally spaced grid over the observed
#'   range. For contrast plots, the default is \code{"equidistant"} when this
#'   argument is omitted.
#' @param alpha Numeric scalar giving the significance level used to construct
#'   confidence intervals; for example, \code{alpha = 0.05} gives 95 percent
#'   intervals.
#' @param show_titles Logical scalar. If \code{TRUE}, add titles describing the
#'   selected FP or ACD representation and any retained spike-at-zero component.
#' @param color_points Colour used for residual points.
#' @param color_line Colour used for the fitted continuous curve.
#' @param color_line_spike Colour used for the structural-zero point.
#' @param color_fill Fill colour used for confidence ribbons.
#' @param shape Point shape used for residual observations.
#' @param size_points Point size used for residual observations.
#' @param size_points_spike Point size used for structural-zero observations.
#' @param linetype Line type used for fitted curves.
#' @param linewidth Line width used for fitted curves.
#' @param alpha_fill Numeric scalar between 0 and 1 controlling confidence-band
#'   transparency.
#' @param ... Reserved for compatibility with the base \code{plot()} generic.
#'
#' @details
#' Confidence intervals are based on the variance-covariance matrix of the final
#' fitted model. They account for uncertainty in the regression coefficients but
#' not for uncertainty introduced by model and power selection. Consequently,
#' they may be narrower than intervals obtained from a full bootstrap procedure.
#'
#' Component-plus-residual plots are available for \code{type = "terms"}.
#' Deviance residuals are used for generalized linear models and martingale
#' residuals for Cox regression.
#'
#' Spike-at-zero covariates are displayed according to the representation
#' retained by the final model:
#' \itemize{
#'   \item both components: the positive-value FP or ACD function and the
#'     structural-zero indicator;
#'   \item positive-part only: the positive-value function without the binary
#'     zero indicator;
#'   \item zero-indicator only: the binary structural-zero component without a
#'     continuous function.
#' }
#'
#' For spike-at-zero terms with a continuous component, the curve is drawn only
#' for \code{x > 0}. A fitted value at \code{x = 0}, when available, is shown as
#' a separate point with a vertical confidence interval because an FP or ACD
#' function is not evaluated at zero.
#'
#' @return A named list of \code{ggplot2} objects, one for each plotted term.
#'   Individual plots can be printed directly or combined with packages such as
#'   \pkg{patchwork}.
#'
#' @examples
#' data("prostate")
#' x <- as.matrix(prostate[, 2:8])
#' y <- as.numeric(prostate$lpsa)
#' fit <- mfp2(x, y, verbose = FALSE)
#'
#' plots <- plot(fit)
#' plots[[1]]
#'
#' @seealso \code{\link{predict.mfp2}}, \code{\link{fracplot}}
#' @importFrom ggplot2 .data
#' @method plot mfp2
#' @export
plot.mfp2 <- function(x,
                      terms = NULL,
                      partial_only = FALSE,
                      type = c("terms", "contrasts"),
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
                      alpha_fill = 0.1,
                      ...) {
  terms_seq_missing <- missing(terms_seq)
  
  plot_mfp2_impl(
    model = x,
    terms = terms,
    partial_only = partial_only,
    type = type,
    ref = ref,
    terms_seq = terms_seq,
    alpha = alpha,
    show_titles = show_titles,
    color_points = color_points,
    color_line = color_line,
    color_line_spike = color_line_spike,
    color_fill = color_fill,
    shape = shape,
    size_points = size_points,
    size_points_spike = size_points_spike,
    linetype = linetype,
    linewidth = linewidth,
    alpha_fill = alpha_fill,
    terms_seq_missing = terms_seq_missing
  )
}

# Internal implementation shared by plot.mfp2() and fracplot().
#
# Keeping the plotting logic here prevents the preferred plot() method from
# passing through the deprecated fracplot() wrapper and emitting a warning.
plot_mfp2_impl <- function(model,
                           terms = NULL,
                           partial_only = FALSE,
                           type = c("terms", "contrasts"),
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
                           alpha_fill = 0.1,
                           terms_seq_missing = FALSE) {
  
  
  # Fail clearly if ggplot2 is unavailable at runtime.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  # Validate the fitted object before prediction.
  if (!inherits(model, "mfp2")) 
    stop("The model entered is not an mfp2 object.", call. = FALSE)
  
  # Contrast references must be supplied as a named list.
  if (!is.null(ref) && !is.list(ref))
    stop("`ref` must be a list or NULL.", call. = FALSE)
  
  if (!is.logical(show_titles) || length(show_titles) != 1L || is.na(show_titles)) {
    stop("show_titles must be a single non-missing logical value.", call. = FALSE)
  }
  
  # set defaults depending on type
  type <- match.arg(type)
  if (type == "contrasts") {
    partial_only <- TRUE
    if (terms_seq_missing)
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
    warning("Variables specified in `terms` were not retained in the final model.", call. = FALSE)
    return(list())
  }
  
  ylab <- "Partial Predictor"
  
  if (!partial_only) {
    # compute residuals to plot points
    # for glm, deviance residuals are required
    # while for cox martingale residuals
    resid <- if (model$family_string == "cox") {
      model$residuals
    } else {
      stats::residuals(model, type = "deviance")
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

#' Deprecated fractional-polynomial plotting interface
#'
#' @description
#' \code{fracplot()} is deprecated. Use \code{plot()} on a fitted \code{mfp2}
#' object instead. The deprecated function remains available temporarily for
#' backward compatibility and returns the same plot list as \code{plot()}.
#'
#' @inheritParams plot.mfp2
#' @param model A fitted \code{mfp2} object.
#'
#' @return A named list of \code{ggplot2} objects, identical to the result from
#'   \code{plot(model, ...)}.
#'
#' @examples
#' data("prostate")
#' x <- as.matrix(prostate[, 2:8])
#' y <- as.numeric(prostate$lpsa)
#' fit <- mfp2(x, y, verbose = FALSE)
#'
#' # Preferred interface
#' plots <- plot(fit)
#'
#' \dontrun{
#' # Deprecated compatibility interface
#' old_plots <- fracplot(fit)
#' }
#'
#' @seealso \code{\link{plot.mfp2}}
#' @export
fracplot <- function(model,
                     terms = NULL,
                     partial_only = FALSE,
                     type = c("terms", "contrasts"),
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
  terms_seq_missing <- missing(terms_seq)
  
  .Deprecated(
    new = "plot",
    package = "mfp2",
    old = "fracplot",
    msg = paste0(
      "`fracplot()` is deprecated and will be removed in a future version ",
      "of mfp2. Use `plot()` with a fitted `mfp2` object instead, for example ",
      "`plot(fit)`."
    )
  )
  
  plot_mfp2_impl(
    model = model,
    terms = terms,
    partial_only = partial_only,
    type = type,
    ref = ref,
    terms_seq = terms_seq,
    alpha = alpha,
    show_titles = show_titles,
    color_points = color_points,
    color_line = color_line,
    color_line_spike = color_line_spike,
    color_fill = color_fill,
    shape = shape,
    size_points = size_points,
    size_points_spike = size_points_spike,
    linetype = linetype,
    linewidth = linewidth,
    alpha_fill = alpha_fill,
    terms_seq_missing = terms_seq_missing
  )
}
