# S3 plot method for mfpi objects
#
# plot.mfpi() produces three types of plots for each continuous variable
# in an mfpi fit (selected or explicitly requested via `terms`):
#
#   "fitted"     - group-specific fitted curves, optionally with 95% CI bands
#   "difference" - pointwise fitted-function difference f_j(x) - f_0(x),
#                  optionally with 95% CI band and reference line at zero
#   "both"       - fitted and difference plots side by side (requires patchwork)
#
# Plots can optionally include rug marks for the observed values of the
# continuous variable. Each plot includes a subtitle showing the variable
# name, functional form, and p-value or information-criterion improvement.


# -----------------------------------------------------------------------------
# plot.mfpi() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Plot Group-Specific Fitted Functions and Fitted-Function Differences from an MFPI Model
#'
#' Produces group-specific fitted-function plots, fitted-function difference
#' plots, or both for continuous variables in an \code{"mfpi"} object. By
#' default, only variables whose interaction was retained as significant are
#' plotted; any variable in \code{cont_vars} can also be plotted by passing
#' its name via \code{terms}. Fitted-function plots and difference plots have
#' separate confidence-band controls, allowing the user to suppress uncertainty
#' bands for group-specific fitted curves while retaining the confidence band
#' for the fitted-function difference curve. Plots are returned invisibly as a
#' nested named list.
#'
#' @section Plot types:
#' Three plot types are available via the \code{plot_type} argument.
#'
#' \code{"fitted"} overlays the fitted linear predictor for the reference
#' group and one non-reference group on shared axes. Confidence bands for these
#' group-specific curves are controlled by \code{show_ci_fitted}.
#'
#' \code{"difference"} plots the pointwise fitted-function difference
#' \eqn{D_j(x) = f_j(x) - f_0(x)} with a horizontal reference line at zero.
#' The confidence band for the difference curve is controlled by
#' \code{show_ci_diff}. Deviations from zero indicate effect modification.
#'
#' \code{"both"} shows fitted and difference plots side by side and requires
#' the \pkg{patchwork} package.
#'
#' @section Confidence bands:
#' The two confidence-band controls are deliberately separate because the
#' inferential roles of the plots differ.
#'
#' \itemize{
#'   \item \code{show_ci_fitted} controls confidence bands around the
#'     group-specific fitted curves. The default is \code{FALSE}, which gives
#'     cleaner and faster fitted-function plots.
#'   \item \code{show_ci_diff} controls the confidence band around the
#'     fitted-function difference curve. The default is \code{TRUE}, because
#'     the difference curve is usually the main inferential display for
#'     interaction.
#' }
#'
#' When a confidence band is suppressed, the corresponding lower and upper
#' confidence-limit columns are not required for plotting.
#'
#' @section Performance:
#' The function constructs only the plot type requested by \code{plot_type}.
#' For example, when \code{plot_type = "difference"}, fitted-function plots are
#' not built internally. This avoids unnecessary \pkg{ggplot2} object creation.
#'
#' Rug marks can be expensive for large datasets because each observation is
#' drawn as a graphical object. The default is \code{show_rug = FALSE} for
#' faster rendering. Set \code{show_rug = TRUE} to display rug marks.
#'
#' Plot printing can also be slow, especially when many variables or group
#' comparisons are plotted. The default is \code{auto_print = TRUE}, so plots
#' are printed as they are created. Set \code{auto_print = FALSE} to return
#' them invisibly for manual printing or saving.
#'
#' @section Return structure:
#' The return value is a nested named list. Outer names are variable names from
#' \code{terms}. Inner names are of the form
#' \code{"grp<j>_vs_grp<0>"} for each non-reference group level \eqn{j}.
#' Each leaf is a \code{ggplot} object when \code{plot_type = "fitted"} or
#' \code{"difference"}, or a \code{patchwork} object when
#' \code{plot_type = "both"}.
#'
#' @section Title and subtitle annotations:
#' Each plot is annotated by default with two text lines.
#'
#' The title describes the structural comparison being plotted without
#' hard-coding colour names. Group identity in fitted-function plots is shown
#' through the legend and line aesthetics.
#' \itemize{
#'   \item Fitted plots: \code{"Group-specific fitted functions"}.
#'   \item Difference plots: \code{"Difference in fitted functions"}.
#' }
#'
#' The subtitle summarises the functional form being shown, the displayed
#' group contrast, and the selection metric, for example
#' \code{"type = Linear; rx: 1 vs 0; p_raw = 0.031 (p_adj = 0.047)"}.
#' The continuous variable is not repeated in the subtitle because it is shown
#' on the x-axis.
#'
#' The metric shown depends on the criterion used to fit the model:
#' \itemize{
#'   \item \code{criterion = "pvalue"}: the subtitle shows \code{p_raw} and,
#'     when \code{p_adjust_method != "none"} and the variable was selected,
#'     also the stored adjusted p-value \code{p_adj}.
#'   \item \code{criterion = "aic"}: the subtitle shows \code{dAIC}.
#'   \item \code{criterion = "bic"}: the subtitle shows \code{dBIC}.
#' }
#'
#' To suppress title and subtitle, set \code{show_title = FALSE}. Axis labels
#' are still shown.
#'
#' @section Plotting non-selected variables:
#' By default (\code{terms = NULL}), only variables whose interaction was
#' retained as significant are plotted. If no variables were selected, an
#' informative message is printed and the function returns invisibly.
#'
#' Any variable in \code{cont_vars} can be plotted regardless of whether its
#' interaction was selected, by passing its name via \code{terms}:
#' \preformatted{plot(fit, terms = "age")}
#' This plots the fitted curve using the functional form pre-specified for that
#' variable in \code{cont_var_forms}, with the same subtitle showing the
#' interaction test metric.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param terms Character vector of variable names to plot. Default
#'   \code{NULL} plots all variables whose interaction was selected as
#'   significant. Pass one or more variable names to plot specific variables
#'   regardless of selection status.
#' @param plot_type Character string specifying the plot to produce. One of
#'   \code{"fitted"}, \code{"difference"}, or \code{"both"}. Default is
#'   \code{"fitted"}.
#' @param auto_print Logical. If \code{TRUE}, each plot is printed as it is
#'   created. Default is \code{TRUE}. With FALSE, plots are returned
#'   invisibly and can be printed, saved, or assembled manually.
#' @param show_title Logical. If \code{TRUE}, titles and subtitles are added.
#'   Default is \code{TRUE}.
#' @param show_rug Logical. If \code{TRUE}, rug marks for the observed values
#'   of the continuous variable are added at the bottom of each plot. Default is
#'   \code{FALSE} for faster rendering with large datasets.
#' @param show_ci_fitted Logical. If \code{TRUE}, 95\% confidence bands are
#'   shown around the group-specific fitted curves. Default is \code{FALSE},
#'   which gives cleaner and faster fitted-function plots.
#' @param show_ci_diff Logical. If \code{TRUE}, a 95\% confidence band is shown
#'   around the fitted-function difference curve. Default is \code{TRUE}.
#' @param colour_ref Character. Colour for the reference-group curve and
#'   confidence band. Default is \code{"#2166AC"}.
#' @param colour_grp Character. Colour for the non-reference-group curve and
#'   confidence band. Default is \code{"#D6604D"}.
#' @param line_types Optional character vector specifying linetypes for the
#'   group-specific fitted-function curves. If \code{NULL}, solid lines are
#'   used for both the reference and comparison groups. The vector may be
#'   unnamed with length 2, in which case values are matched to the reference
#'   and comparison groups in order, or named using the original group labels,
#'   for example \code{c("placebo" = "solid", "standard" = "dashed")}.
#' @param colour_diff Character. Colour for the difference curve and confidence
#'   band. Default is \code{"#1B7837"}.
#' @param linewidth Positive numeric. Line width for fitted and difference
#'   curves. Default is \code{1}.
#' @param ribbon_alpha Numeric in \eqn{[0, 1]}. Transparency of confidence
#'   bands. Default is \code{0.2}.
#' @param rug_alpha Numeric in \eqn{[0, 1]}. Transparency of rug marks. Ignored
#'   when \code{show_rug = FALSE}. Default is \code{0.4}.
#' @param legend_position Legend position for the fitted-function plot. May be
#'   one of `"inside"`, `"right"`, `"left"`, `"bottom"`, `"top"`, or `"none"`.
#'   The default is `"inside"`, which places the legend inside the plotting
#'   panel.
#' @param legend_inside Numeric vector of length 2 giving the legend position
#'   inside the plotting panel when `legend_position = "inside"`. Values are
#'   given in normalized parent coordinates, with `c(0, 0)` at the bottom-left
#'   and `c(1, 1)` at the top-right. The default is `c(0.02, 0.98)`.
#' @param legend_justification Numeric vector of length 2 giving the legend
#'   justification when `legend_position = "inside"`. The default is `c(0, 1)`,
#'   which anchors the top-left corner of the legend at `legend_inside`.
#' @param ... Currently unused.
#'
#' @return A nested named list of \code{ggplot} or \code{patchwork} objects,
#'   returned invisibly.
#'
#' @seealso \code{mfpi()}, \code{gen_fitted_values_per_group()}
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- mfpi(x, y, group_var = "trt", cont_vars = c("age", "bmi"))
#'
#' # Return selected fitted-function plots without printing
#' plots <- plot(fit)
#'
#' # Print one plot manually
#' plots$age$grp1_vs_grp0
#'
#' # Difference plot with confidence band, returned without printing
#' plots <- plot(fit, plot_type = "difference")
#'
#' # Print difference plots immediately
#' plot(fit, plot_type = "difference", auto_print = TRUE)
#'
#' # Add confidence bands to group-specific fitted curves
#' plot(fit, plot_type = "fitted",
#'      show_ci_fitted = TRUE,
#'      auto_print = TRUE)
#'
#' # Custom linetypes for the reference and comparison groups
#' plot(fit, plot_type = "fitted",
#'      line_types = c("solid", "dotdash"))
#'
#' # Suppress all confidence bands
#' plots <- plot(fit, plot_type = "both",
#'               show_ci_fitted = FALSE,
#'               show_ci_diff = FALSE)
#'
#' # Show rug marks if desired
#' plot(fit, plot_type = "difference",
#'      show_rug = TRUE,
#'      auto_print = TRUE)
#'
#' # Plot a non-selected variable by name
#' plot(fit, terms = "age", plot_type = "difference", auto_print = TRUE)
#'
#' # Save a specific plot
#' ggplot2::ggsave(
#'   "age_difference.pdf",
#'   plots$age$grp1_vs_grp0,
#'   width = 7,
#'   height = 5
#' )
#' }
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
  
  model     <- x
  plot_type <- match.arg(plot_type)
  legend_position <- match.arg(legend_position)
  
  if (!inherits(model, "mfpi")) {
    stop("! `x` must be an object of class \"mfpi\".", call. = FALSE)
  }
  
  if (!is.logical(auto_print) || length(auto_print) != 1L || is.na(auto_print)) {
    stop("! `auto_print` must be a single TRUE or FALSE value.", call. = FALSE)
  }
  
  if (!is.logical(show_title) || length(show_title) != 1L || is.na(show_title)) {
    stop("! `show_title` must be a single TRUE or FALSE value.", call. = FALSE)
  }
  
  if (!is.logical(show_rug) || length(show_rug) != 1L || is.na(show_rug)) {
    stop("! `show_rug` must be a single TRUE or FALSE value.", call. = FALSE)
  }
  
  if (!is.logical(show_ci_fitted) || length(show_ci_fitted) != 1L ||
      is.na(show_ci_fitted)) {
    stop("! `show_ci_fitted` must be a single TRUE or FALSE value.",
         call. = FALSE)
  }
  
  if (!is.logical(show_ci_diff) || length(show_ci_diff) != 1L ||
      is.na(show_ci_diff)) {
    stop("! `show_ci_diff` must be a single TRUE or FALSE value.",
         call. = FALSE)
  }
  
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
    stop("! `legend_inside` must be a numeric vector of length 2 with values in [0, 1].",
         call. = FALSE)
  }
  
  if (!is.numeric(legend_justification) || length(legend_justification) != 2L ||
      any(is.na(legend_justification))) {
    stop("! `legend_justification` must be a numeric vector of length 2.",
         call. = FALSE)
  }
  
  if (!is.numeric(linewidth) || length(linewidth) != 1L ||
      is.na(linewidth) || linewidth <= 0) {
    stop("! `linewidth` must be a single positive numeric value.", call. = FALSE)
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
  
  if (plot_type == "both" && !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "! `plot_type = \"both\"` requires the patchwork package. ",
      "Install it with install.packages(\"patchwork\").",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Extract fitted functions and metrics
  # ---------------------------------------------------------------------------
  metrics <- model$best_model_metrics
  
  # Use best_fitted_functions for selected variables (already unwrapped matrices)
  # and all_fitted_functions for non-selected variables (need to unwrap the
  # single-element named list to get the matrix directly).
  best_ff <- model$best_fitted_functions
  all_ff  <- model$all_fitted_functions
  
  if (is.null(all_ff) || length(all_ff) == 0L) {
    stop(
      "! No fitted functions found in this model object.",
      call. = FALSE
    )
  }
  
  all_fitted <- stats::setNames(
    lapply(names(all_ff), function(v) {
      # Prefer best_fitted_functions (already a matrix) for selected variables
      if (length(best_ff) > 0L && !is.null(best_ff[[v]])) {
        return(best_ff[[v]])
      }
      # For non-selected variables, unwrap the single-element named list
      inner <- all_ff[[v]]
      if (is.null(inner) || length(inner) == 0L) return(NULL)
      if (is.matrix(inner)) return(inner)          # already unwrapped
      inner[[1L]]                                   # unwrap list(form = matrix)
    }),
    names(all_ff)
  )
  
  # ---------------------------------------------------------------------------
  # Resolve terms
  # ---------------------------------------------------------------------------
  if (is.null(terms)) {
    # Default: plot only selected variables
    selected_vars <- if (!is.null(metrics) && nrow(metrics) > 0L) {
      metrics$variable
    } else {
      character(0L)
    }
    
    terms <- intersect(selected_vars, names(all_fitted))
    
    if (length(terms) == 0L) {
      message(
        "Nothing to plot: no variables in `cont_vars` were selected as\n",
        "having a significant interaction with the group variable.\n",
        "To inspect a variable's fitted curve anyway, pass its name via\n",
        "the `terms` argument, e.g. terms = \"age\"."
      )
      return(invisible(x))
    }
  } else {
    if (!is.character(terms) || length(terms) == 0L) {
      stop("! `terms` must be a non-empty character vector or NULL.",
           call. = FALSE)
    }
    
    # Error only for variables not in cont_vars at all
    truly_missing <- setdiff(terms, names(all_fitted))
    if (length(truly_missing) > 0L) {
      stop(
        "! The following `terms` are not variables in this model:\n  ",
        paste(truly_missing, collapse = ", "), ".\n",
        "  Available variables: ",
        paste(names(all_fitted), collapse = ", "), ".",
        call. = FALSE
      )
    }
    
    # Remove any terms with no fitted function (e.g. compute_fitted = FALSE)
    no_fit <- vapply(terms, function(v) {
      fv <- all_fitted[[v]]
      is.null(fv) || length(fv) == 0L
    }, logical(1L))
    
    if (any(no_fit)) {
      stop(
        "! The following variable(s) have no fitted function available:\n  ",
        paste(terms[no_fit], collapse = ", "), ".\n",
        "  Re-fit with compute_fitted = TRUE.",
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Determine metric display
  # ---------------------------------------------------------------------------
  crit_used <- if (is.null(model$criterion)) "pvalue" else model$criterion
  padj      <- if (!is.null(model$p_adjust_method)) model$p_adjust_method else "none"
  adjusting <- padj != "none" && crit_used == "pvalue"
  
  subtitle_lookup <- list()
  
  add_subtitle_info <- function(var, form, metric, label, metric_adj = NULL) {
    subtitle_lookup[[var]] <<- list(
      type       = form,
      metric     = metric,
      label      = label,
      metric_adj = metric_adj
    )
  }
  
  # Build subtitle info from all_model_metrics so all terms (selected or not)
  # get a subtitle. For selected variables, also pick up adjusted p-values
  # from best_model_metrics where available.
  all_metrics <- model$all_model_metrics
  
  if (!is.null(all_metrics) && nrow(all_metrics) > 0L) {
    for (i in seq_len(nrow(all_metrics))) {
      v    <- all_metrics$variable[i]
      form <- all_metrics$type[i]
      
      if (crit_used == "pvalue") {
        metric <- if ("pvalue" %in% names(all_metrics)) all_metrics$pvalue[i] else NA_real_
        # p_adjusted is available in all_metrics for all variables when
        # p_adjust_method != "none", regardless of selection status.
        metric_adj <- if (adjusting &&
                          "p_adjusted" %in% names(all_metrics) &&
                          !is.na(all_metrics$p_adjusted[i])) {
          all_metrics$p_adjusted[i]
        } else {
          NULL
        }
        add_subtitle_info(v, form, metric, "p_raw", metric_adj)
      } else if (crit_used == "aic") {
        metric <- if ("dAIC" %in% names(all_metrics)) all_metrics$dAIC[i]
        else if ("AIC_main_minus_int" %in% names(all_metrics)) all_metrics$AIC_main_minus_int[i]
        else NA_real_
        add_subtitle_info(v, form, metric, "dAIC", NULL)
      } else if (crit_used == "bic") {
        metric <- if ("dBIC" %in% names(all_metrics)) all_metrics$dBIC[i]
        else if ("BIC_main_minus_int" %in% names(all_metrics)) all_metrics$BIC_main_minus_int[i]
        else NA_real_
        add_subtitle_info(v, form, metric, "dBIC", NULL)
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Group-level metadata
  # ---------------------------------------------------------------------------
  group_levels_original <- model$group_levels_original
  group_levels_new      <- model$group_levels_new
  group_label           <- model$group_var
  
  if (is.null(group_levels_original) || length(group_levels_original) < 2L) {
    stop(
      "! `model$group_levels_original` must contain at least two group levels.",
      call. = FALSE
    )
  }
  
  if (is.null(group_levels_new) || length(group_levels_new) != length(group_levels_original)) {
    group_levels_new <- seq_along(group_levels_original) - 1L
  }
  
  ref_code  <- group_levels_new[1L]
  non_ref_code <- group_levels_new[-1L]
  ref_label <- as.character(group_levels_original[1L])
  non_ref_label <- as.character(group_levels_original[-1L])
  
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
  # Helper: subtitle
  # ---------------------------------------------------------------------------
  make_subtitle <- function(var, grp_label_value) {
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
    
    contrast_text <- sprintf("%s: %s vs %s", group_label, grp_label_value, ref_label)
    
    metric_text <- sprintf("%s = %s", lbl, val_str)
    
    if (!is.null(info$metric_adj) && !is.na(info$metric_adj)) {
      adj_str <- trimws(formatC(info$metric_adj, format = "g", digits = 3))
      metric_text <- sprintf("%s (p_adj = %s)", metric_text, adj_str)
    }
    
    sprintf("type = %s; %s; %s", type_str, contrast_text, metric_text)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: fitted plot for one variable / group comparison
  # ---------------------------------------------------------------------------
  make_fitted_plot <- function(var, data_mat, grp_code, grp_label_value) {
    col_ref    <- sprintf("f%s", ref_code)
    col_grp    <- sprintf("f%s", grp_code)
    col_ref_lo <- sprintf("f%s_lower", ref_code)
    col_ref_hi <- sprintf("f%s_upper", ref_code)
    col_grp_lo <- sprintf("f%s_lower", grp_code)
    col_grp_hi <- sprintf("f%s_upper", grp_code)
    
    needed <- c(var, col_ref, col_grp)
    
    if (isTRUE(show_ci_fitted)) {
      needed <- c(
        needed,
        col_ref_lo,
        col_ref_hi,
        col_grp_lo,
        col_grp_hi
      )
    }
    
    missing <- setdiff(needed, names(data_mat))
    
    if (length(missing) > 0L) {
      warning(
        paste0(
          "! Skipping fitted plot for ",
          var,
          " (group ",
          grp_label_value,
          "): missing columns: ",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
      
      return(NULL)
    }
    
    xvec <- data_mat[[var]]
    
    ref_df <- data.frame(
      x = xvec,
      y = data_mat[[col_ref]]
    )
    
    grp_df <- data.frame(
      x = xvec,
      y = data_mat[[col_grp]]
    )
    
    if (isTRUE(show_ci_fitted)) {
      ref_df$lo <- data_mat[[col_ref_lo]]
      ref_df$hi <- data_mat[[col_ref_hi]]
      grp_df$lo <- data_mat[[col_grp_lo]]
      grp_df$hi <- data_mat[[col_grp_hi]]
    }
    
    ref_df <- ref_df[order(ref_df$x), , drop = FALSE]
    grp_df <- grp_df[order(grp_df$x), , drop = FALSE]
    
    ref_df$group <- ref_label
    grp_df$group <- grp_label_value
    line_df <- rbind(ref_df[, c("x", "y", "group"), drop = FALSE],
                     grp_df[, c("x", "y", "group"), drop = FALSE])
    line_df$group <- factor(
      line_df$group,
      levels = c(ref_label, grp_label_value)
    )
    
    group_values <- levels(line_df$group)
    colour_values <- stats::setNames(c(colour_ref, colour_grp), group_values)
    linetype_values <- resolve_line_types(group_values)
    
    p <- ggplot2::ggplot()
    
    if (isTRUE(show_ci_fitted)) {
      p <- p +
        ggplot2::geom_ribbon(
          data = ref_df,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          fill  = colour_ref,
          alpha = ribbon_alpha,
          show.legend = FALSE
        ) +
        ggplot2::geom_ribbon(
          data = grp_df,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          fill  = colour_grp,
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
      ggplot2::ylab("Fitted function") +
      ggplot2::labs(
        title = if (show_title) {
          "Group-specific fitted functions"
        } else {
          NULL
        },
        subtitle = if (show_title) make_subtitle(var, grp_label_value) else NULL
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
      rug_df <- data.frame(x = xvec)
      
      p <- p +
        ggplot2::geom_rug(
          data = rug_df,
          ggplot2::aes(x = .data$x),
          sides = "b",
          alpha = rug_alpha
        )
    }
    
    p
  }
  
  # ---------------------------------------------------------------------------
  # Helper: difference plot for one variable / group comparison
  # ---------------------------------------------------------------------------
  make_diff_plot <- function(var, data_mat, grp_code, grp_label_value) {
    col_diff  <- sprintf("f%s-f%s", grp_code, ref_code)
    col_lower <- sprintf("(f%s-f%s)_lower", grp_code, ref_code)
    col_upper <- sprintf("(f%s-f%s)_upper", grp_code, ref_code)
    
    needed <- c(var, col_diff)
    
    if (isTRUE(show_ci_diff)) {
      needed <- c(needed, col_lower, col_upper)
    }
    
    missing <- setdiff(needed, names(data_mat))
    
    if (length(missing) > 0L) {
      warning(
        paste0(
          "! Skipping difference plot for ",
          var,
          " (group ",
          grp_label_value,
          "): missing columns: ",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
      
      return(NULL)
    }
    
    xvec <- data_mat[[var]]
    
    plot_df <- data.frame(
      x    = xvec,
      diff = data_mat[[col_diff]]
    )
    
    if (isTRUE(show_ci_diff)) {
      plot_df$lower <- data_mat[[col_lower]]
      plot_df$upper <- data_mat[[col_upper]]
    }
    
    plot_df <- plot_df[order(plot_df$x), , drop = FALSE]
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x))
    
    if (isTRUE(show_ci_diff)) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          fill  = colour_diff,
          alpha = ribbon_alpha
        )
    }
    
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype   = "dashed",
        colour     = "grey50"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$diff),
        colour    = colour_diff,
        linewidth = linewidth
      ) +
      ggplot2::xlab(var) +
      ggplot2::ylab(sprintf("f(%s = %s) - f(%s = %s)",
                            group_label, grp_label_value,
                            group_label, ref_label)) +
      ggplot2::labs(
        title = if (show_title) {
          "Difference in fitted functions"
        } else {
          NULL
        },
        subtitle = if (show_title) make_subtitle(var, grp_label_value) else NULL
      ) +
      ggplot2::theme_bw()
    
    if (isTRUE(show_rug)) {
      rug_df <- data.frame(x = xvec)
      
      p <- p +
        ggplot2::geom_rug(
          data = rug_df,
          ggplot2::aes(x = .data$x),
          sides = "b",
          alpha = rug_alpha
        )
    }
    
    p
  }
  
  # ---------------------------------------------------------------------------
  # Build all plots
  # ---------------------------------------------------------------------------
  plots <- stats::setNames(vector("list", length(terms)), terms)
  
  for (var in terms) {
    data_mat <- as.data.frame(all_fitted[[var]])
    
    grp_plots <- stats::setNames(
      vector("list", length(non_ref_code)),
      sprintf("grp%s_vs_grp%s", non_ref_code, ref_code)
    )
    
    for (j in seq_along(non_ref_code)) {
      grp_j <- non_ref_code[j]
      grp_label_value <- non_ref_label[j]
      
      p_out <- switch(
        plot_type,
        
        fitted = {
          make_fitted_plot(var, data_mat, grp_j, grp_label_value)
        },
        
        difference = {
          make_diff_plot(var, data_mat, grp_j, grp_label_value)
        },
        
        both = {
          p_fit  <- make_fitted_plot(var, data_mat, grp_j, grp_label_value)
          p_diff <- make_diff_plot(var, data_mat, grp_j, grp_label_value)
          
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