# -----------------------------------------------------------------------------
# plot.mfpi() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Plot MFPI Fitted Functions and Fitted-Function Differences
#'
#' Produces diagnostic plots for selected or requested continuous variables in an
#' object of class \code{"mfpi"}. The method can display group-specific fitted
#' functions, fitted-function differences between a comparison group and the
#' reference group, or both plot types side by side.
#'
#' The plotting method is deliberately prediction-based. It does not rely on
#' precomputed fitted-function matrices stored in the fitted \code{"mfpi"}
#' object. Instead, plotting data are generated on demand by
#' \code{predict.mfpi()} using the same prediction path used for ordinary MFPI
#' prediction. This keeps plotting consistent with prediction and avoids storing
#' large plotting grids inside the fitted model object.
#'
#' @section Plot types:
#' Three plot types are available through \code{plot_type}.
#'
#' \describe{
#'   \item{\code{"fitted"}}{
#'     Plots the group-specific fitted functions for the reference group and one
#'     comparison group on shared axes. For a continuous variable \eqn{x}, the
#'     plotted curves correspond to \eqn{\hat f_r(x)} for the reference group
#'     and \eqn{\hat f_j(x)} for comparison group \eqn{j}. Confidence bands for
#'     these group-specific fitted functions are controlled by
#'     \code{show_ci_fitted}. The y-axis label identifies the model scale of
#'     the covariate function, for example outcome scale, log-odds, log-mean,
#'     or log-hazard.
#'   }
#'   \item{\code{"difference"}}{
#'     Plots the pointwise fitted-function difference
#'     \eqn{\hat f_j(x) - \hat f_r(x)} for one comparison group \eqn{j} versus
#'     the reference group \eqn{r}. A horizontal zero line is added to make the
#'     absence of a fitted-function difference visually explicit. Confidence
#'     bands for the difference curve are controlled by \code{show_ci_diff}.
#'     The y-axis keeps the explicit contrast, \eqn{\hat f_j(x) - \hat f_r(x)},
#'     and appends the model scale in parentheses.
#'   }
#'   \item{\code{"both"}}{
#'     Places the fitted-function plot and the fitted-function difference plot
#'     side by side for each variable and group contrast. This option requires
#'     the \pkg{patchwork} package.
#'   }
#' }
#'
#' @section Model scale and interpretation:
#' The plotted values are estimated covariate functions, not full fitted
#' responses. For Gaussian models the fitted function is on the outcome scale.
#' For binomial models it is on the log-odds scale. For Poisson models it is on
#' the log-mean scale. For Cox models it is on the log-hazard scale.
#'
#' Difference plots show differences between estimated covariate functions on
#' the same model scale. They should not be interpreted as probabilities, odds
#' ratios, rate ratios, hazard ratios, survival probabilities, or complete
#' predicted responses unless a separate transformation is explicitly applied.
#'
#' @section Prediction source:
#' For each plotted variable, the method calls \code{predict.mfpi()} with
#' \code{type = "both"} and \code{grid = TRUE}. Standard errors are requested
#' only when at least one requested plot requires confidence bands. The plotting
#' data are taken from the long-format \code{functions} and \code{differences}
#' components returned by \code{predict.mfpi()}.
#'
#' The \code{functions} component contains one row per evaluation value and
#' group, with columns such as \code{term}, \code{x}, \code{group}, and
#' \code{fit}. The \code{differences} component contains one row per
#' evaluation value and non-reference contrast, with columns such as
#' \code{term}, \code{x}, \code{contrast}, \code{group}, \code{reference},
#' and \code{fit}. Confidence limits are read from \code{lower} and
#' \code{upper} when standard errors were requested.
#'
#' @section Group labels and reference group:
#' MFPI stores group information in two forms. Internal group codes, usually
#' \code{0}, \code{1}, ..., are used in coefficient names and prediction
#' metadata. Original group labels, such as \code{"Placebo"} and
#' \code{"Treatment"}, are used for display. If the fitted object contains a
#' \code{group_level_map} component, plotting uses that mapping preferentially.
#' Otherwise it falls back to \code{group_levels_original} and
#' \code{group_levels_new}.
#'
#' Plot list names are intentionally based on internal group codes, for example
#' \code{"grp1_vs_grp0"}. This keeps returned object names syntactically stable
#' even when original group labels contain spaces, punctuation, or non-standard
#' characters. The visual labels inside the plots use the original group labels.
#'
#' @section Term selection:
#' If \code{terms = NULL}, only selected interaction terms are plotted. Selected
#' variables are obtained from \code{x$best_interaction_model}, and prediction is
#' requested with \code{model = "best"}.
#'
#' If \code{terms} is supplied, the requested variables are plotted from the
#' full stored winner set using \code{model = "all"}. This allows inspection of
#' non-selected variables, provided the corresponding term-specific winner models
#' are stored in the fitted object.
#'
#' @section Confidence bands:
#' Confidence-band controls are separate because fitted-function plots and
#' difference plots have different inferential roles.
#'
#' \itemize{
#'   \item \code{show_ci_fitted} controls confidence bands around the
#'     group-specific fitted functions. The default is \code{FALSE} to keep the
#'     fitted-function display visually simple.
#'   \item \code{show_ci_diff} controls the confidence band around the
#'     fitted-function difference curve. The default is \code{TRUE} because the
#'     difference curve is usually the primary interaction display.
#' }
#'
#' @section Return value:
#' The method returns a nested named list invisibly. Outer names are continuous
#' variable names. Inner names are contrast labels of the form
#' \code{"grp<j>_vs_grp<r>"}, where \code{j} is the internal code of the
#' comparison group and \code{r} is the internal code of the reference group.
#' Each leaf is a \code{ggplot} object when \code{plot_type = "fitted"} or
#' \code{plot_type = "difference"}, and a \code{patchwork} object when
#' \code{plot_type = "both"} and both component plots are available.
#'
#' @param x Object of class \code{"mfpi"}.
#' @param terms Optional character vector of continuous variables to plot. If
#'   \code{NULL}, only selected interaction terms are plotted. If supplied, the
#'   requested terms are plotted regardless of selection status, provided stored
#'   term-specific winner models are available.
#' @param plot_type Character scalar. One of \code{"fitted"},
#'   \code{"difference"}, or \code{"both"}.
#' @param auto_print Logical scalar. If \code{TRUE}, each plot is printed as it
#'   is created. If \code{FALSE}, plots are returned invisibly for manual
#'   printing, modification, or saving.
#' @param show_title Logical scalar. If \code{TRUE}, plot titles and subtitles
#'   are added. Subtitles include interaction type, group contrast, and the
#'   selection metric available for the plotted variable.
#' @param show_rug Logical scalar. If \code{TRUE}, rug marks for the plotted
#'   \eqn{x} grid values are added at the bottom of each plot.
#' @param show_ci_fitted Logical scalar. If \code{TRUE}, confidence bands are
#'   shown around the group-specific fitted functions.
#' @param show_ci_diff Logical scalar. If \code{TRUE}, a confidence band is
#'   shown around the fitted-function difference curve.
#' @param legend_position Character scalar. Legend position for fitted-function
#'   plots. One of \code{"inside"}, \code{"right"}, \code{"left"},
#'   \code{"bottom"}, \code{"top"}, or \code{"none"}.
#' @param legend_inside Numeric vector of length 2 giving the legend position
#'   inside the plotting panel when \code{legend_position = "inside"}. Values
#'   are interpreted in normalized parent coordinates.
#' @param legend_justification Numeric vector of length 2 giving legend
#'   justification when \code{legend_position = "inside"}.
#' @param colour_ref Character scalar. Colour for the reference-group fitted
#'   function and its confidence band.
#' @param colour_grp Character scalar. Colour for the comparison-group fitted
#'   function and its confidence band.
#' @param line_types Optional character vector of linetypes for the two fitted
#'   functions. If unnamed, it must have length 2 and is interpreted as
#'   reference group followed by comparison group. If named, names must match the
#'   original group labels shown in the legend.
#' @param colour_diff Character scalar. Colour for the difference curve and its
#'   confidence band.
#' @param linewidth Positive numeric scalar. Line width for fitted-function and
#'   difference curves.
#' @param ribbon_alpha Numeric scalar in \eqn{[0, 1]}. Transparency of confidence
#'   bands.
#' @param rug_alpha Numeric scalar in \eqn{[0, 1]}. Transparency of rug marks.
#' @param ... Currently unused. Supplied arguments trigger a warning.
#'
#' @return Invisibly returns a nested named list of \code{ggplot} or
#'   \code{patchwork} objects.
#'
#' @seealso \code{mfpi()}, \code{predict.mfpi()}
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- mfpi(x, y, group_var = "trt", cont_vars = c("age", "bmi"))
#'
#' # Plot selected fitted functions.
#' plot(fit)
#'
#' # Return selected difference plots without printing.
#' plots <- plot(fit, plot_type = "difference", auto_print = FALSE)
#'
#' # Plot a non-selected variable if its winner model is stored.
#' plot(fit, terms = "age", plot_type = "both")
#'
#' # Add confidence bands to group-specific fitted functions.
#' plot(fit, plot_type = "fitted", show_ci_fitted = TRUE)
#'
#' # Use named line types matching original group labels.
#' plot(fit, plot_type = "fitted",
#'      line_types = c(Placebo = "solid", Treatment = "dashed"))
#'
#' # Save a returned plot.
#' ggplot2::ggsave("age_difference.pdf", plots$age$grp1_vs_grp0)
#' }
#' @method plot mfpi
#' @export
plot.mfpi <- function(x,
                      terms          = NULL,
                      plot_type      = c("fitted", "difference", "both"),
                      auto_print     = TRUE,
                      show_title     = TRUE,
                      show_rug       = FALSE,
                      show_ci_fitted = FALSE,
                      show_ci_diff   = TRUE,
                      legend_position = c("inside", "right", "left", "bottom", "top", "none"),
                      legend_inside  = c(0.02, 0.98),
                      legend_justification = c(0, 1),
                      colour_ref     = "#2166AC",
                      colour_grp     = "#D6604D",
                      line_types     = NULL,
                      colour_diff    = "#1B7837",
                      linewidth      = 1,
                      ribbon_alpha   = 0.2,
                      rug_alpha      = 0.4,
                      ...) {
  # ---------------------------------------------------------------------------
  # Basic argument validation
  # ---------------------------------------------------------------------------
  model <- x
  plot_type <- match.arg(plot_type)
  legend_position <- match.arg(legend_position)
  
  if (!inherits(model, "mfpi")) {
    stop("! `x` must be an object of class \"mfpi\".", call. = FALSE)
  }
  
  check_logical_scalar <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("! `", name, "` must be a single TRUE or FALSE value.",
           call. = FALSE)
    }
    invisible(TRUE)
  }
  
  check_logical_scalar(auto_print, "auto_print")
  check_logical_scalar(show_title, "show_title")
  check_logical_scalar(show_rug, "show_rug")
  check_logical_scalar(show_ci_fitted, "show_ci_fitted")
  check_logical_scalar(show_ci_diff, "show_ci_diff")
  
  if (!is.null(line_types) &&
      (!is.character(line_types) || length(line_types) == 0L ||
       any(is.na(line_types)))) {
    stop(
      "! `line_types` must be NULL or a non-empty character vector without NA values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(legend_inside) || length(legend_inside) != 2L ||
      any(is.na(legend_inside)) || any(legend_inside < 0) ||
      any(legend_inside > 1)) {
    stop(
      "! `legend_inside` must be a numeric vector of length 2 with values in [0, 1].",
      call. = FALSE
    )
  }
  
  if (!is.numeric(legend_justification) ||
      length(legend_justification) != 2L ||
      any(is.na(legend_justification))) {
    stop(
      "! `legend_justification` must be a numeric vector of length 2.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(linewidth) || length(linewidth) != 1L ||
      is.na(linewidth) || linewidth <= 0) {
    stop("! `linewidth` must be a single positive numeric value.",
         call. = FALSE)
  }
  
  if (!is.numeric(ribbon_alpha) || length(ribbon_alpha) != 1L ||
      is.na(ribbon_alpha) || ribbon_alpha < 0 || ribbon_alpha > 1) {
    stop("! `ribbon_alpha` must be a single numeric value in [0, 1].",
         call. = FALSE)
  }
  
  if (!is.numeric(rug_alpha) || length(rug_alpha) != 1L ||
      is.na(rug_alpha) || rug_alpha < 0 || rug_alpha > 1) {
    stop("! `rug_alpha` must be a single numeric value in [0, 1].",
         call. = FALSE)
  }
  
  dots <- list(...)
  if (length(dots) > 0L) {
    warning(
      "Unused arguments in `...`: ", paste(names(dots), collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  if (plot_type == "both" && !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "! `plot_type = \"both\"` requires the patchwork package. ",
      "Install it with install.packages(\"patchwork\").",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Resolve terms and prediction model scope
  # ---------------------------------------------------------------------------
  # Plotting data are now generated on demand by predict.mfpi().
  if (is.null(terms)) {
    selected_terms <- names(model$best_interaction_model)
    selected_terms <- selected_terms[!is.na(selected_terms) & nzchar(selected_terms)]
    
    if (length(selected_terms) == 0L) {
      message(
        "Nothing to plot: no variables in `cont_vars` were selected as\n",
        "having a significant interaction with the group variable.\n",
        "To inspect a variable's fitted curve anyway, pass its name via\n",
        "the `terms` argument, e.g. terms = \"age\"."
      )
      return(invisible(x))
    }
    
    terms <- selected_terms
    predict_model <- "best"
  } else {
    if (!is.character(terms) || length(terms) == 0L || anyNA(terms)) {
      stop("! `terms` must be a non-empty character vector or NULL.",
           call. = FALSE)
    }
    
    terms <- unique(terms)
    predict_model <- "all"
  }
  
  # ---------------------------------------------------------------------------
  # Metric/subtitle metadata
  # ---------------------------------------------------------------------------
  subtitle_lookup <- mfpi_plot_subtitle_lookup(model)
  
  # Group metadata translates internal prediction codes such as "0" and "1"
  # into the original group labels shown to users in legends and axis labels.
  group_meta <- mfpi_plot_group_metadata(model)
  group_label <- group_meta$group_var
  
  # ---------------------------------------------------------------------------
  # Helper: subtitle for one variable/contrast
  # ---------------------------------------------------------------------------
  make_subtitle <- function(var, grp_label_value, ref_label_value) {
    info <- subtitle_lookup[[var]]
    
    if (is.null(info)) {
      return(var)
    }
    
    type_str <- switch(
      info$type,
      linear = "Linear",
      fp1    = "FP1",
      fp2    = "FP2",
      info$type
    )
    
    lbl <- if (!is.null(info$label)) info$label else "metric"
    
    val_str <- if (is.null(info$metric) || is.na(info$metric)) {
      "NA"
    } else if (lbl %in% c("p", "p_raw")) {
      trimws(formatC(info$metric, format = "g", digits = 3))
    } else {
      trimws(formatC(info$metric, format = "f", digits = 2))
    }
    
    contrast_text <- sprintf(
      "%s: %s vs %s",
      group_label,
      grp_label_value,
      ref_label_value
    )
    
    metric_text <- sprintf("%s = %s", lbl, val_str)
    
    if (!is.null(info$metric_adj) && !is.na(info$metric_adj)) {
      adj_str <- trimws(formatC(info$metric_adj, format = "g", digits = 3))
      metric_text <- sprintf("%s (p_adj = %s)", metric_text, adj_str)
    }
    
    sprintf("type = %s; %s; %s", type_str, contrast_text, metric_text)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: linetype mapping for fitted-function plots
  # ---------------------------------------------------------------------------
  resolve_line_types <- function(labels) {
    if (is.null(line_types)) {
      return(stats::setNames(rep("solid", length(labels)), labels))
    }
    
    lt_names <- names(line_types)
    has_names <- !is.null(lt_names) && any(nzchar(lt_names))
    
    if (has_names) {
      missing_labels <- setdiff(labels, lt_names)
      if (length(missing_labels) > 0L) {
        stop(
          "! Named `line_types` must include entries for the original group labels used in the plot: ",
          paste(missing_labels, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      return(stats::setNames(unname(line_types[labels]), labels))
    }
    
    if (length(line_types) != 2L) {
      stop(
        "! Unnamed `line_types` must have length 2: one for the reference group and one for the comparison group.",
        call. = FALSE
      )
    }
    
    stats::setNames(line_types, labels)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: fitted-function plot for one variable and one group contrast
  # ---------------------------------------------------------------------------
  make_fitted_plot <- function(var, functions_df,
                               ref_label_value, grp_label_value) {
    if (is.null(functions_df) || !is.data.frame(functions_df) ||
        nrow(functions_df) == 0L) {
      warning(
        paste0("! Skipping fitted plot for ", var,
               ": prediction output contains no fitted functions."),
        call. = FALSE
      )
      return(NULL)
    }
    
    needed <- c("x", "group", "fit")
    if (isTRUE(show_ci_fitted)) needed <- c(needed, "lower", "upper")
    missing <- setdiff(needed, names(functions_df))
    
    if (length(missing) > 0L) {
      warning(
        paste0(
          "! Skipping fitted plot for ", var,
          ": prediction functions are missing column(s): ",
          paste(missing, collapse = ", "), "."
        ),
        call. = FALSE
      )
      return(NULL)
    }
    
    # predict.mfpi() now returns fitted functions in long format. Filter the
    # requested term defensively, then split by the user-facing group labels.
    fdat <- functions_df
    if ("term" %in% names(fdat)) {
      fdat <- fdat[fdat$term == var, , drop = FALSE]
    }
    fdat$group <- as.character(fdat$group)
    
    ref_raw <- fdat[fdat$group == ref_label_value, , drop = FALSE]
    grp_raw <- fdat[fdat$group == grp_label_value, , drop = FALSE]
    
    if (nrow(ref_raw) == 0L || nrow(grp_raw) == 0L) {
      warning(
        paste0(
          "! Skipping fitted plot for ", var,
          ": prediction functions do not contain both groups ",
          ref_label_value, " and ", grp_label_value, "."
        ),
        call. = FALSE
      )
      return(NULL)
    }
    
    ref_df <- data.frame(
      x = ref_raw$x,
      y = ref_raw$fit,
      stringsAsFactors = FALSE
    )
    
    grp_df <- data.frame(
      x = grp_raw$x,
      y = grp_raw$fit,
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(show_ci_fitted)) {
      ref_df$lo <- ref_raw$lower
      ref_df$hi <- ref_raw$upper
      grp_df$lo <- grp_raw$lower
      grp_df$hi <- grp_raw$upper
    }
    
    # Lines and ribbons must be sorted by x for visually coherent curves.
    ref_df <- ref_df[order(ref_df$x), , drop = FALSE]
    grp_df <- grp_df[order(grp_df$x), , drop = FALSE]
    
    ref_df$group <- ref_label_value
    grp_df$group <- grp_label_value
    
    line_df <- rbind(
      ref_df[, c("x", "y", "group"), drop = FALSE],
      grp_df[, c("x", "y", "group"), drop = FALSE]
    )
    
    line_df$group <- factor(line_df$group,
                            levels = c(ref_label_value, grp_label_value))
    
    group_values <- levels(line_df$group)
    colour_values <- stats::setNames(c(colour_ref, colour_grp), group_values)
    linetype_values <- resolve_line_types(group_values)
    
    p <- ggplot2::ggplot()
    
    if (isTRUE(show_ci_fitted)) {
      p <- p +
        ggplot2::geom_ribbon(
          data = ref_df,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          fill = colour_ref,
          alpha = ribbon_alpha,
          show.legend = FALSE
        ) +
        ggplot2::geom_ribbon(
          data = grp_df,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          fill = colour_grp,
          alpha = ribbon_alpha,
          show.legend = FALSE
        )
    }
    
    p <- p +
      ggplot2::geom_line(
        data = line_df,
        ggplot2::aes(x = .data$x, y = .data$y,
                     colour = .data$group, linetype = .data$group),
        linewidth = linewidth
      ) +
      ggplot2::scale_colour_manual(values = colour_values, name = group_label) +
      ggplot2::scale_linetype_manual(values = linetype_values, name = group_label) +
      ggplot2::xlab(var) +
      # Make the model scale explicit. These curves are covariate functions,
      # not full fitted responses; for non-Gaussian families they are plotted
      # on the corresponding link/linear-predictor scale.
      ggplot2::ylab(
        paste0("Covariate function: ", mfpi_plot_scale_label(model))
      ) +
      ggplot2::labs(
        title = if (show_title) "Group-specific fitted functions" else NULL,
        subtitle = if (show_title) {
          make_subtitle(var, grp_label_value, ref_label_value)
        } else {
          NULL
        }
      ) +
      ggplot2::theme_bw()
    
    if (legend_position == "inside") {
      p <- p +
        ggplot2::theme(
          legend.position = legend_inside,
          legend.justification = legend_justification,
          legend.background = ggplot2::element_rect(
            fill = grDevices::adjustcolor("white", alpha.f = 0.85),
            colour = "grey80"
          )
        )
    } else {
      p <- p + ggplot2::theme(legend.position = legend_position)
    }
    
    if (isTRUE(show_rug)) {
      p <- p +
        ggplot2::geom_rug(
          data = data.frame(x = unique(line_df$x)),
          ggplot2::aes(x = .data$x),
          sides = "b",
          alpha = rug_alpha
        )
    }
    
    p
  }
  
  # ---------------------------------------------------------------------------
  # Helper: fitted-function difference plot for one group contrast
  # ---------------------------------------------------------------------------
  make_diff_plot <- function(var, differences_df,
                             ref_label_value, grp_label_value) {
    if (is.null(differences_df) || !is.data.frame(differences_df) ||
        nrow(differences_df) == 0L) {
      warning(
        paste0("! Skipping difference plot for ", var,
               ": prediction output contains no fitted-function differences."),
        call. = FALSE
      )
      return(NULL)
    }
    
    needed <- c("x", "group", "reference", "fit")
    if (isTRUE(show_ci_diff)) needed <- c(needed, "lower", "upper")
    missing <- setdiff(needed, names(differences_df))
    
    if (length(missing) > 0L) {
      warning(
        paste0(
          "! Skipping difference plot for ", var,
          ": prediction differences are missing column(s): ",
          paste(missing, collapse = ", "), "."
        ),
        call. = FALSE
      )
      return(NULL)
    }
    
    # predict.mfpi() returns differences in long format. Filter by term and
    # display labels instead of reconstructing synthetic matrix column names.
    ddat <- differences_df
    if ("term" %in% names(ddat)) {
      ddat <- ddat[ddat$term == var, , drop = FALSE]
    }
    ddat$group <- as.character(ddat$group)
    ddat$reference <- as.character(ddat$reference)
    
    diff_raw <- ddat[
      ddat$group == grp_label_value & ddat$reference == ref_label_value,
      ,
      drop = FALSE
    ]
    
    if (nrow(diff_raw) == 0L) {
      warning(
        paste0(
          "! Skipping difference plot for ", var,
          ": prediction differences do not contain contrast ",
          grp_label_value, " vs ", ref_label_value, "."
        ),
        call. = FALSE
      )
      return(NULL)
    }
    
    plot_df <- data.frame(
      x = diff_raw$x,
      diff = diff_raw$fit,
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(show_ci_diff)) {
      plot_df$lower <- diff_raw$lower
      plot_df$upper <- diff_raw$upper
    }
    
    plot_df <- plot_df[order(plot_df$x), , drop = FALSE]
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x))
    
    if (isTRUE(show_ci_diff)) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          fill = colour_diff,
          alpha = ribbon_alpha
        )
    }
    
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$diff),
        colour = colour_diff,
        linewidth = linewidth
      ) +
      ggplot2::xlab(var) +
      # Keep the original contrast-specific label, but append the model scale so
      # the difference is not mistaken for a response-scale probability, rate,
      # ratio, survival probability, or complete prediction.
      ggplot2::ylab(sprintf(
        "f(%s = %s) - f(%s = %s) (%s)",
        group_label,
        grp_label_value,
        group_label,
        ref_label_value,
        mfpi_plot_scale_label(model)
      )) +
      ggplot2::labs(
        title = if (show_title) "Difference in fitted functions" else NULL,
        subtitle = if (show_title) {
          make_subtitle(var, grp_label_value, ref_label_value)
        } else {
          NULL
        }
      ) +
      ggplot2::theme_bw()
    
    if (isTRUE(show_rug)) {
      p <- p +
        ggplot2::geom_rug(
          data = data.frame(x = unique(plot_df$x)),
          ggplot2::aes(x = .data$x),
          sides = "b",
          alpha = rug_alpha
        )
    }
    
    p
  }
  
  # ---------------------------------------------------------------------------
  # Generate prediction objects and build plots
  # ---------------------------------------------------------------------------
  plots <- stats::setNames(vector("list", length(terms)), terms)
  
  for (var in terms) {
    # Developer note:
    # plot.mfpi() intentionally requests grid-based prediction. Observed-row
    # prediction is still available directly via predict.mfpi(grid = FALSE), but
    # smooth plotting should use the grid output.
    need_se <- (plot_type %in% c("fitted", "both") && isTRUE(show_ci_fitted)) ||
      (plot_type %in% c("difference", "both") && isTRUE(show_ci_diff))
    
    pred <- predict(
      model,
      terms = var,
      model = predict_model,
      type = "both",
      se.fit = need_se,
      grid = TRUE
    )
    
    functions_df <- pred$functions
    differences_df <- pred$differences
    
    if (plot_type %in% c("fitted", "both") &&
        (is.null(functions_df) || !is.data.frame(functions_df) ||
         nrow(functions_df) == 0L)) {
      stop(
        paste0("Prediction output for term `", var, "` does not contain fitted functions."),
        call. = FALSE
      )
    }
    
    if (plot_type %in% c("difference", "both") &&
        (is.null(differences_df) || !is.data.frame(differences_df) ||
         nrow(differences_df) == 0L)) {
      stop(
        paste0("Prediction output for term `", var, "` does not contain fitted-function differences."),
        call. = FALSE
      )
    }
    
    display_map <- pred$metadata$group_display_labels
    pred_group_codes <- names(display_map)
    
    if (is.null(pred_group_codes) || anyNA(pred_group_codes) ||
        any(!nzchar(pred_group_codes))) {
      pred_group_codes <- names(pred$metadata$coefficient_groups)
    }
    if (is.null(pred_group_codes) || anyNA(pred_group_codes) ||
        any(!nzchar(pred_group_codes))) {
      pred_group_codes <- as.character(seq_along(pred$metadata$coefficient_groups) - 1L)
    }
    
    if (is.null(display_map) || length(display_map) != length(pred_group_codes)) {
      display_map <- stats::setNames(
        vapply(
          pred_group_codes,
          function(code) mfpi_plot_original_group_label(group_meta, code),
          character(1L)
        ),
        pred_group_codes
      )
    } else {
      display_map <- stats::setNames(as.character(display_map), pred_group_codes)
    }
    
    ref_code <- as.character(pred$metadata$reference)
    if (is.null(ref_code) || length(ref_code) != 1L || is.na(ref_code) ||
        !nzchar(ref_code)) {
      ref_code <- pred_group_codes[1L]
    }
    
    comparison_codes <- setdiff(pred_group_codes, ref_code)
    
    if (length(comparison_codes) == 0L) {
      warning(
        paste0("No non-reference group contrasts are available for term `", var, "`."),
        call. = FALSE
      )
      plots[[var]] <- list()
      next
    }
    
    # Prediction metadata keeps internal group codes for stable object names and
    # user-facing labels for plot legends, subtitles, and axes. Prefer the
    # labels returned by predict.mfpi(), with fitted-object metadata as fallback.
    ref_label <- if (!is.null(pred$metadata$reference_label) &&
                     length(pred$metadata$reference_label) == 1L &&
                     !is.na(pred$metadata$reference_label) &&
                     nzchar(pred$metadata$reference_label)) {
      as.character(pred$metadata$reference_label)
    } else if (ref_code %in% names(display_map)) {
      unname(display_map[ref_code])
    } else {
      mfpi_plot_original_group_label(group_meta, ref_code)
    }
    
    comparison_labels <- vapply(
      comparison_codes,
      function(code) {
        if (code %in% names(display_map)) {
          unname(display_map[code])
        } else {
          mfpi_plot_original_group_label(group_meta, code)
        }
      },
      character(1L)
    )
    
    # Keep returned list names based on internal codes. This makes them stable
    # even when original group labels contain spaces or punctuation.
    grp_plots <- stats::setNames(
      vector("list", length(comparison_codes)),
      paste0("grp", comparison_codes, "_vs_grp", ref_code)
    )
    
    for (j in seq_along(comparison_codes)) {
      grp_code <- comparison_codes[j]
      grp_label_value <- comparison_labels[j]
      
      p_out <- switch(
        plot_type,
        fitted = make_fitted_plot(
          var = var,
          functions_df = functions_df,
          ref_label_value = ref_label,
          grp_label_value = grp_label_value
        ),
        difference = make_diff_plot(
          var = var,
          differences_df = differences_df,
          ref_label_value = ref_label,
          grp_label_value = grp_label_value
        ),
        both = {
          p_fit <- make_fitted_plot(
            var = var,
            functions_df = functions_df,
            ref_label_value = ref_label,
            grp_label_value = grp_label_value
          )
          p_diff <- make_diff_plot(
            var = var,
            differences_df = differences_df,
            ref_label_value = ref_label,
            grp_label_value = grp_label_value
          )
          
          if (!is.null(p_fit) && !is.null(p_diff)) {
            patchwork::wrap_plots(p_fit, p_diff, ncol = 2L)
          } else if (!is.null(p_fit)) {
            p_fit
          } else {
            p_diff
          }
        }
      )
      
      grp_plots[[j]] <- p_out
      
      if (isTRUE(auto_print) && !is.null(p_out)) {
        print(p_out)
      }
    }
    
    plots[[var]] <- grp_plots
  }
  
  invisible(plots)
}


# -----------------------------------------------------------------------------
# plot.mfpi() helper functions ------------------------------------------------
# -----------------------------------------------------------------------------

#' Return the Plot Scale Label for an MFPI Model
#'
#' Determines the short model-scale label used on MFPI plot y-axes. The helper
#' intentionally returns compact labels because axis titles must remain readable
#' when plots are displayed side by side. Longer interpretation details are kept
#' in the \code{plot.mfpi()} documentation.
#'
#' The returned scale describes the scale of the estimated covariate function
#' produced by MFPI prediction. It does not describe a full fitted response. In
#' particular, binomial curves are labelled as log-odds, Poisson curves as
#' log-mean, and Cox curves as log-hazard.
#'
#' @param model Object of class \code{"mfpi"}. The helper reads
#'   \code{model$family_string} when available and falls back to
#'   \code{model$family$family} for older objects.
#'
#' @return Character scalar. One of \code{"outcome scale"},
#'   \code{"log-odds"}, \code{"log-mean"}, \code{"log-hazard"}, or
#'   \code{"linear predictor"} for unrecognised families.
#'
#' @keywords internal
#' @noRd
mfpi_plot_scale_label <- function(model) {
  # Prefer the explicit family string stored by mfpi(). Use [[, exact = TRUE]]
  # rather than $ so missing or partially matching component names do not affect
  # the chosen label.
  family_string <- model[["family_string", exact = TRUE]]
  
  # Older or manually constructed objects may not contain family_string. In that
  # case, fall back to the standard family object structure used by glm-like
  # models. The fallback is also exact to avoid partial-matching surprises.
  if (is.null(family_string)) {
    family_obj <- model[["family", exact = TRUE]]
    if (is.list(family_obj)) {
      family_string <- family_obj[["family", exact = TRUE]]
    }
  }
  
  # If family information is unavailable or malformed, use the conservative
  # generic label. This keeps plotting robust for older fitted objects while
  # avoiding misleading family-specific labels.
  if (is.null(family_string) || length(family_string) == 0L ||
      is.na(family_string[1L]) || !nzchar(as.character(family_string[1L]))) {
    return("linear predictor")
  }
  
  # Normalise to a single lower-case string. This makes the switch robust to
  # accidental vector values while preserving a safe default for unknown input.
  family_string <- tolower(as.character(family_string[1L]))
  
  switch(
    family_string,
    gaussian = "outcome scale",
    binomial = "log-odds",
    poisson  = "log-mean",
    cox      = "log-hazard",
    "linear predictor"
  )
}


#' Build Subtitle Metadata for MFPI Plots
#'
#' Creates a named lookup table containing the interaction form and the
#' criterion-specific model-selection metric for each continuous variable stored
#' in an \code{"mfpi"} object.
#'
#' This helper is used only by \code{plot.mfpi()}. It reads
#' \code{model$all_model_metrics}, determines the criterion used during model
#' selection, and extracts the statistic that is most relevant for plot
#' subtitles. For p-value selection, the raw p-value is stored and the adjusted
#' p-value is stored when p-value adjustment was used. For information-criterion
#' selection, the stored AIC or BIC improvement is used.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return A named list indexed by continuous-variable name. Each element is a
#'   list with components \code{type}, \code{metric}, \code{label}, and, for
#'   adjusted p-value selection, \code{metric_adj}. If no model metrics are
#'   available, an empty list is returned.
#'
#' @keywords internal
#' @noRd
mfpi_plot_subtitle_lookup <- function(model) {
  crit_used <- if (is.null(model$criterion)) "pvalue" else model$criterion
  padj <- if (!is.null(model$p_adjust_method)) model$p_adjust_method else "none"
  adjusting <- padj != "none" && crit_used == "pvalue"
  
  out <- list()
  all_metrics <- model$all_model_metrics
  
  if (is.null(all_metrics) || nrow(all_metrics) == 0L) {
    return(out)
  }
  
  for (i in seq_len(nrow(all_metrics))) {
    v <- all_metrics$variable[i]
    form <- if ("type" %in% names(all_metrics)) all_metrics$type[i] else NA_character_
    
    if (is.na(v) || !nzchar(v)) next
    
    if (crit_used == "pvalue") {
      metric <- if ("pvalue" %in% names(all_metrics)) all_metrics$pvalue[i] else NA_real_
      metric_adj <- if (adjusting &&
                        "p_adjusted" %in% names(all_metrics) &&
                        !is.na(all_metrics$p_adjusted[i])) {
        all_metrics$p_adjusted[i]
      } else {
        NULL
      }
      out[[v]] <- list(
        type = form,
        metric = metric,
        label = "p_raw",
        metric_adj = metric_adj
      )
    } else if (crit_used == "aic") {
      metric <- if ("dAIC" %in% names(all_metrics)) {
        all_metrics$dAIC[i]
      } else if ("AIC_main_minus_int" %in% names(all_metrics)) {
        all_metrics$AIC_main_minus_int[i]
      } else {
        NA_real_
      }
      out[[v]] <- list(type = form, metric = metric, label = "dAIC")
    } else if (crit_used == "bic") {
      metric <- if ("dBIC" %in% names(all_metrics)) {
        all_metrics$dBIC[i]
      } else if ("BIC_main_minus_int" %in% names(all_metrics)) {
        all_metrics$BIC_main_minus_int[i]
      } else {
        NA_real_
      }
      out[[v]] <- list(type = form, metric = metric, label = "dBIC")
    }
  }
  
  out
}


#' Extract Group-Level Metadata for MFPI Plots
#'
#' Returns the grouping-variable name, original group labels, and internal group
#' codes used by \code{plot.mfpi()} to translate prediction metadata into
#' user-facing plot labels.
#'
#' MFPI prediction keeps internal group codes for coefficient alignment and
#' stable returned object names. Plot labels, however, should use the original
#' group labels supplied by the user, for example \code{"Placebo"} and
#' \code{"Treatment"}.
#'
#' If the fitted object contains \code{group_level_map}, this helper uses it as
#' the primary source of truth. The map is expected to contain at least columns
#' \code{original} and \code{internal}. If the map is absent, the helper falls
#' back to \code{group_levels_original} and \code{group_levels_new}. If internal
#' levels are absent or inconsistent with the original labels, consecutive
#' internal codes \code{0}, \code{1}, ... are reconstructed.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return List with elements \code{group_var}, \code{original}, and
#'   \code{internal}. All returned elements are character vectors. The
#'   \code{original} and \code{internal} vectors have the same length and define
#'   a one-to-one display mapping.
#'
#' @keywords internal
#' @noRd
mfpi_plot_group_metadata <- function(model) {
  group_var <- model$group_var
  
  if (is.null(group_var) || length(group_var) != 1L || is.na(group_var) ||
      !nzchar(group_var)) {
    stop("! The MFPI object must store a valid `group_var`.", call. = FALSE)
  }
  
  # Prefer the explicit level map when available. This is the most robust
  # representation because it keeps display labels, input values, and internal
  # modelling codes together in one object.
  level_map <- model$group_level_map
  
  if (!is.null(level_map) &&
      is.data.frame(level_map) &&
      all(c("original", "internal") %in% names(level_map))) {
    original <- level_map$original
    internal <- level_map$internal
  } else {
    original <- model$group_levels_original
    internal <- model$group_levels_new
  }
  
  if (is.null(original) || length(original) < 2L) {
    stop(
      "! `model$group_levels_original` must contain at least two group levels.",
      call. = FALSE
    )
  }
  
  # Older fitted objects may not contain stored internal levels. Reconstruct the
  # default internal coding used by MFPI so plotting remains backward-compatible.
  if (is.null(internal) || length(internal) != length(original)) {
    internal <- seq_along(original) - 1L
  }
  
  list(
    group_var = as.character(group_var),
    original  = as.character(original),
    internal  = as.character(internal)
  )
}


#' Map an Internal Group Code to the Original Group Label
#'
#' Converts an internal group code used by prediction output, such as
#' \code{"0"} or \code{"1"}, into the corresponding original group label used
#' for plot display.
#'
#' This helper is deliberately tolerant. If the internal code cannot be matched
#' to the stored group metadata, it returns the internal code unchanged. That
#' fallback avoids failing an otherwise valid plot when an older object lacks
#' complete group-label metadata.
#'
#' @param group_meta Group metadata returned by
#'   \code{mfpi_plot_group_metadata()}.
#' @param internal_code Character scalar giving an internal group code.
#'
#' @return Character scalar. The original group label when a match is available;
#'   otherwise \code{internal_code} converted to character.
#'
#' @keywords internal
#' @noRd
mfpi_plot_original_group_label <- function(group_meta, internal_code) {
  internal_code <- as.character(internal_code)[1L]
  hit <- match(internal_code, group_meta$internal)
  
  if (!is.na(hit)) {
    return(group_meta$original[hit])
  }
  
  internal_code
}
