# -----------------------------------------------------------------------------
# plot.mfpi() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Plot MFPI Fitted Functions and Fitted-Function Differences
#'
#' Produces plots of group-specific MFPI fitted functions, fitted-function
#' differences, or both for continuous variables in an object of class
#' \code{"mfpi"}.
#'
#' Plotting data are computed on demand by \code{predict.mfpi()}. The plotting
#' method no longer reads precomputed fitted-function matrices from
#' \code{object$best_fitted_functions} or \code{object$all_fitted_functions}.
#' This means fitted functions are reconstructed only when requested, using the
#' same prediction path used by \code{predict(object, type = "both")}.
#'
#' @section Plot types:
#' Three plot types are available through \code{plot_type}.
#'
#' \describe{
#'   \item{\code{"fitted"}}{
#'     Plots the group-specific fitted functions \eqn{\hat f_j(x)} for the
#'     reference group and one comparison group on shared axes. Confidence bands
#'     for these curves are controlled by \code{show_ci_fitted}.
#'   }
#'   \item{\code{"difference"}}{
#'     Plots the pointwise fitted-function difference
#'     \eqn{\hat f_j(x) - \hat f_r(x)} with a horizontal reference line at zero.
#'     Confidence bands for the difference curve are controlled by
#'     \code{show_ci_diff}.
#'   }
#'   \item{\code{"both"}}{
#'     Places the fitted-function plot and the fitted-function difference plot
#'     side by side. This option requires the \pkg{patchwork} package.
#'   }
#' }
#'
#' @section Prediction source:
#' The method calls \code{predict.mfpi()} internally with
#' \code{type = "both"}, \code{grid = TRUE}, and \code{se.fit = TRUE}. The
#' resulting \code{wide} matrix is used for plotting. Therefore this plot method
#' depends on \code{predict.mfpi()} and the shared prediction helper
#' \code{build_group_fp_basis()}, not on the legacy fitted-function generator.
#'
#' @section Term selection:
#' If \code{terms = NULL}, only selected interaction terms are plotted. These are
#' obtained from \code{object$best_interaction_model}, using
#' \code{model = "best"} in \code{predict.mfpi()}.
#'
#' If \code{terms} is supplied, those variables are plotted from the wider
#' stored winner set using \code{model = "all"}. This allows inspection of
#' non-selected continuous variables, provided their term-specific winner models
#' are stored in the object.
#'
#' @section Confidence bands:
#' The confidence-band controls are separate because the fitted-function and
#' difference plots have different inferential roles.
#'
#' \itemize{
#'   \item \code{show_ci_fitted} controls confidence bands around the
#'     group-specific fitted functions. The default is \code{FALSE} for cleaner
#'     fitted-function plots.
#'   \item \code{show_ci_diff} controls the confidence band around the
#'     fitted-function difference curve. The default is \code{TRUE} because the
#'     difference curve is usually the primary interaction display.
#' }
#'
#' @section Return value:
#' The return value is a nested named list, returned invisibly. Outer names are
#' continuous-variable names. Inner names are contrast labels of the form
#' \code{"grp<j>_vs_grp<r>"}, where \code{j} is the comparison group and
#' \code{r} is the reference group. Each leaf is a \code{ggplot} object for
#' \code{plot_type = "fitted"} or \code{"difference"}, or a \code{patchwork}
#' object for \code{plot_type = "both"}.
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
#'   printing or saving.
#' @param show_title Logical scalar. If \code{TRUE}, plot titles and subtitles
#'   are added.
#' @param show_rug Logical scalar. If \code{TRUE}, rug marks for the plotted
#'   \eqn{x} values are added at the bottom of each plot.
#' @param show_ci_fitted Logical scalar. If \code{TRUE}, confidence bands are
#'   shown around the group-specific fitted functions.
#' @param show_ci_diff Logical scalar. If \code{TRUE}, a confidence band is
#'   shown around the fitted-function difference curve.
#' @param legend_position Character scalar. Legend position for fitted-function
#'   plots. One of \code{"inside"}, \code{"right"}, \code{"left"},
#'   \code{"bottom"}, \code{"top"}, or \code{"none"}.
#' @param legend_inside Numeric vector of length 2 giving the legend position
#'   inside the plotting panel when \code{legend_position = "inside"}.
#' @param legend_justification Numeric vector of length 2 giving legend
#'   justification when \code{legend_position = "inside"}.
#' @param colour_ref Character scalar. Colour for the reference-group fitted
#'   function and confidence band.
#' @param colour_grp Character scalar. Colour for the comparison-group fitted
#'   function and confidence band.
#' @param line_types Optional character vector of linetypes for the two fitted
#'   functions. If unnamed, it must have length 2 and is interpreted as
#'   reference then comparison. If named, names must match the original group
#'   labels shown in the legend.
#' @param colour_diff Character scalar. Colour for the difference curve and its
#'   confidence band.
#' @param linewidth Positive numeric scalar. Line width for fitted and
#'   difference curves.
#' @param ribbon_alpha Numeric scalar in \eqn{[0, 1]}. Transparency of
#'   confidence bands.
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
#' # Save a returned plot.
#' ggplot2::ggsave("age_difference.pdf", plots$age$grp1_vs_grp0)
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
  make_fitted_plot <- function(var, data_mat, ref_code, grp_code,
                               ref_label_value, grp_label_value) {
    col_ref    <- paste0("f", ref_code)
    col_grp    <- paste0("f", grp_code)
    col_ref_lo <- paste0("f", ref_code, "_lower")
    col_ref_hi <- paste0("f", ref_code, "_upper")
    col_grp_lo <- paste0("f", grp_code, "_lower")
    col_grp_hi <- paste0("f", grp_code, "_upper")
    
    needed <- c(var, col_ref, col_grp)
    
    if (isTRUE(show_ci_fitted)) {
      needed <- c(needed, col_ref_lo, col_ref_hi, col_grp_lo, col_grp_hi)
    }
    
    missing <- setdiff(needed, names(data_mat))
    
    if (length(missing) > 0L) {
      warning(
        paste0(
          "! Skipping fitted plot for ", var,
          " (group ", grp_label_value, "): missing columns: ",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
      return(NULL)
    }
    
    ref_df <- data.frame(
      x = data_mat[[var]],
      y = data_mat[[col_ref]],
      stringsAsFactors = FALSE
    )
    
    grp_df <- data.frame(
      x = data_mat[[var]],
      y = data_mat[[col_grp]],
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(show_ci_fitted)) {
      ref_df$lo <- data_mat[[col_ref_lo]]
      ref_df$hi <- data_mat[[col_ref_hi]]
      grp_df$lo <- data_mat[[col_grp_lo]]
      grp_df$hi <- data_mat[[col_grp_hi]]
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
      ggplot2::ylab("Fitted function") +
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
          data = data.frame(x = data_mat[[var]]),
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
  make_diff_plot <- function(var, data_mat, ref_code, grp_code,
                             ref_label_value, grp_label_value) {
    col_diff  <- paste0("f", grp_code, "-f", ref_code)
    col_lower <- paste0("(f", grp_code, "-f", ref_code, ")_lower")
    col_upper <- paste0("(f", grp_code, "-f", ref_code, ")_upper")
    
    needed <- c(var, col_diff)
    
    if (isTRUE(show_ci_diff)) {
      needed <- c(needed, col_lower, col_upper)
    }
    
    missing <- setdiff(needed, names(data_mat))
    
    if (length(missing) > 0L) {
      warning(
        paste0(
          "! Skipping difference plot for ", var,
          " (group ", grp_label_value, "): missing columns: ",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
      return(NULL)
    }
    
    plot_df <- data.frame(
      x = data_mat[[var]],
      diff = data_mat[[col_diff]],
      stringsAsFactors = FALSE
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
      ggplot2::ylab(sprintf(
        "f(%s = %s) - f(%s = %s)",
        group_label,
        grp_label_value,
        group_label,
        ref_label_value
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
          data = data.frame(x = data_mat[[var]]),
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
    
    data_mat <- as.data.frame(pred$wide)
    
    if (!(var %in% names(data_mat))) {
      stop(
        paste0("Prediction output for term `", var, "` does not contain an x column."),
        call. = FALSE
      )
    }
    
    pred_group_codes <- names(pred$metadata$coefficient_groups)
    if (is.null(pred_group_codes) || anyNA(pred_group_codes) ||
        any(!nzchar(pred_group_codes))) {
      pred_group_codes <- as.character(seq_along(pred$metadata$coefficient_groups) - 1L)
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
    
    ref_label <- mfpi_plot_original_group_label(group_meta, ref_code)
    comparison_labels <- vapply(
      comparison_codes,
      function(code) mfpi_plot_original_group_label(group_meta, code),
      character(1L)
    )
    
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
          data_mat = data_mat,
          ref_code = ref_code,
          grp_code = grp_code,
          ref_label_value = ref_label,
          grp_label_value = grp_label_value
        ),
        difference = make_diff_plot(
          var = var,
          data_mat = data_mat,
          ref_code = ref_code,
          grp_code = grp_code,
          ref_label_value = ref_label,
          grp_label_value = grp_label_value
        ),
        both = {
          p_fit <- make_fitted_plot(
            var = var,
            data_mat = data_mat,
            ref_code = ref_code,
            grp_code = grp_code,
            ref_label_value = ref_label,
            grp_label_value = grp_label_value
          )
          p_diff <- make_diff_plot(
            var = var,
            data_mat = data_mat,
            ref_code = ref_code,
            grp_code = grp_code,
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

#' Build Subtitle Metadata for MFPI Plots
#'
#' Creates a named lookup table containing functional-form and model-selection
#' metric information for all continuous variables stored in an \code{"mfpi"}
#' object.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return Named list indexed by continuous-variable name. Each element contains
#'   \code{type}, \code{metric}, \code{label}, and optionally
#'   \code{metric_adj}.
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
#' Reads the original and internal group labels stored in an \code{"mfpi"}
#' object and returns them in a normalized structure used by the plotting method.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return List with elements \code{group_var}, \code{original}, and
#'   \code{internal}.
#'
#' @keywords internal
#' @noRd
mfpi_plot_group_metadata <- function(model) {
  group_var <- model$group_var
  
  if (is.null(group_var) || length(group_var) != 1L || is.na(group_var) ||
      !nzchar(group_var)) {
    stop("! The MFPI object must store a valid `group_var`.", call. = FALSE)
  }
  
  original <- model$group_levels_original
  internal <- model$group_levels_new
  
  if (is.null(original) || length(original) < 2L) {
    stop(
      "! `model$group_levels_original` must contain at least two group levels.",
      call. = FALSE
    )
  }
  
  if (is.null(internal) || length(internal) != length(original)) {
    internal <- seq_along(original) - 1L
  }
  
  list(
    group_var = as.character(group_var),
    original = as.character(original),
    internal = as.character(internal)
  )
}


#' Map an Internal Group Code to the Original Group Label
#'
#' Converts the internal group code used by fitted coefficient names, for example
#' \code{"0"} or \code{"1"}, to the original group label supplied by the user
#' when the model was fitted.
#'
#' @param group_meta Group metadata returned by
#'   \code{mfpi_plot_group_metadata()}.
#' @param internal_code Character scalar giving an internal group code.
#'
#' @return Character scalar with the original group label when available;
#'   otherwise \code{internal_code}.
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
