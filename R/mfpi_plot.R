# S3 plot method for mfpi objects
#
# plot.mfpi() produces three types of plots for each significant continuous
# variable in an mfpi fit:
#
#   "fitted"      - group-specific fitted curves, optionally with 95% CI bands
#   "difference" - pointwise treatment-effect difference f_j(x) - f_0(x),
#                  optionally with 95% CI band and reference line at zero
#   "both"        - fitted and difference plots side by side (requires patchwork)
#
# Plots can optionally include rug marks for the observed values of the
# continuous variable. Each plot can include a subtitle showing the variable
# name, functional form, and p-value or information-criterion improvement.


# -----------------------------------------------------------------------------
# plot.mfpi() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Plot Fitted Functions and Treatment-Effect Differences from an MFPI Model
#'
#' Produces group-specific fitted-function plots, treatment-effect difference
#' plots, or both for each significant continuous variable in an
#' \code{"mfpi"} object. Fitted-function plots and difference plots have
#' separate confidence-band controls, allowing the user to suppress uncertainty
#' bands for group-specific fitted curves while retaining the confidence band
#' for the treatment-effect difference curve. Plots are returned invisibly as a
#' nested named list.
#'
#' @section Plot types:
#' Three plot types are available via the \code{plot_type} argument.
#'
#' \code{"fitted"} overlays the fitted linear predictor for the reference
#' group and one non-reference group on shared axes. Confidence bands for these
#' group-specific curves are controlled by \code{show_ci_fitted}.
#'
#' \code{"difference"} plots the pointwise treatment-effect function
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
#'     treatment-effect difference curve. The default is \code{TRUE}, because
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
#' Line and ribbon data are sorted by the continuous variable before plotting.
#' This avoids zig-zag paths when fitted values are stored in observation order
#' and improves rendering behaviour for \code{geom_line()} and
#' \code{geom_ribbon()}.
#'
#' Rug marks can be expensive for large datasets because each observation is
#' drawn as a graphical object. The default is \code{show_rug = FALSE} for
#' faster rendering. Set \code{show_rug = TRUE} to display rug marks.
#'
#' Plot printing can also be slow, especially when many variables or group
#' comparisons are plotted. The default is \code{auto_print = FALSE}, so plots
#' are returned invisibly and can be printed manually.
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
#' The title describes the structural comparison being plotted:
#' \itemize{
#'   \item Fitted plots:
#'     \code{"<group_var>: group <j> (red) vs group <0> (blue)"}.
#'   \item Difference plots:
#'     \code{"<group_var>: treatment effect (group <j> - group <0>)"}.
#' }
#'
#' The subtitle summarises the variable, the functional form being shown, and
#' the selection metric for that form, in the format
#' \code{"<variable> | <type> | <metric> = <value>"}.
#'
#' The metric shown depends on the criterion used to fit the model:
#' \itemize{
#'   \item \code{criterion = "pvalue"}: shows the interaction-test p-value.
#'   \item \code{criterion = "aic"}: shows \code{dAIC}, defined as
#'     \eqn{\mathrm{AIC}_{\mathrm{main}} -
#'     \mathrm{AIC}_{\mathrm{interaction}}}. Positive values favour the
#'     interaction model.
#'   \item \code{criterion = "bic"}: shows \code{dBIC}, defined analogously.
#' }
#'
#' When \code{p_adjust_method != "none"} and \code{criterion = "pvalue"}, the
#' p-value label adapts to avoid ambiguity:
#' \itemize{
#'   \item If \code{type = NULL}, selected variables are plotted and the
#'     subtitle shows both the raw and adjusted p-values.
#'   \item If \code{type} is specified, candidate curves are plotted and the
#'     subtitle shows \code{p_raw}; the adjusted p-value is not shown because
#'     it applies to the winning candidate, not necessarily to the specific
#'     candidate being displayed.
#' }
#'
#' To suppress title and subtitle, set \code{show_title = FALSE}. Axis labels
#' are still shown.
#'
#' @section Plotting non-selected variables:
#' Whether a variable can be plotted depends on \code{type}.
#'
#' Under \code{type = NULL}, only variables whose interaction was retained as
#' significant have a selected curve stored in
#' \code{model$best_fitted_functions}. If no variables were selected, the
#' function prints an informative message and returns invisibly.
#'
#' Under \code{type = "linear"}, \code{"fp1"}, or \code{"fp2"}, the function
#' plots the requested candidate curve from \code{model$all_fitted_functions},
#' provided the model was fitted with \code{compute_fitted = TRUE}. This allows
#' inspection of candidate curves even for variables that were not selected.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param terms Character vector of variable names to plot. Default
#'   \code{NULL} plots all variables with fitted values available under the
#'   current \code{type} setting.
#' @param type Character string or \code{NULL}. Controls which functional form
#'   is plotted. One of \code{"linear"}, \code{"fp1"}, \code{"fp2"}, or
#'   \code{NULL}. Default \code{NULL} plots the selected functional form for
#'   each retained variable.
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
#'   around the treatment-effect difference curve. Default is \code{TRUE}.
#' @param colour_ref Character. Colour for the reference-group curve and
#'   confidence band. Default is \code{"#2166AC"}.
#' @param colour_grp Character. Colour for the non-reference-group curve and
#'   confidence band. Default is \code{"#D6604D"}.
#' @param colour_diff Character. Colour for the difference curve and confidence
#'   band. Default is \code{"#1B7837"}.
#' @param linewidth Positive numeric. Line width for fitted and difference
#'   curves. Default is \code{1}.
#' @param ribbon_alpha Numeric in \eqn{[0, 1]}. Transparency of confidence
#'   bands. Default is \code{0.2}.
#' @param rug_alpha Numeric in \eqn{[0, 1]}. Transparency of rug marks. Ignored
#'   when \code{show_rug = FALSE}. Default is \code{0.4}.
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
#' # Inspect candidate FP2 curve for age, even if not selected
#' plot(fit, terms = "age", type = "fp2",
#'      plot_type = "difference",
#'      auto_print = TRUE)
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
                      type           = NULL,
                      plot_type      = c("fitted", "difference", "both"),
                      auto_print     = TRUE,
                      show_title     = TRUE,
                      show_rug       = FALSE,
                      show_ci_fitted = FALSE,
                      show_ci_diff   = TRUE,
                      colour_ref     = "#2166AC",
                      colour_grp     = "#D6604D",
                      colour_diff    = "#1B7837",
                      linewidth      = 1,
                      ribbon_alpha   = 0.2,
                      rug_alpha      = 0.4,
                      ...) {
  
  model     <- x
  plot_type <- match.arg(plot_type)
  
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
  
  valid_types <- c("linear", "fp1", "fp2")
  
  if (is.null(type)) {
    all_fitted <- model$best_fitted_functions
    
    if (is.null(all_fitted)) {
      stop(
        "! No fitted functions found. Re-fit with compute_fitted = TRUE.",
        call. = FALSE
      )
    }
  } else {
    type <- match.arg(type, valid_types)
    
    if (is.null(model$all_fitted_functions) ||
        length(model$all_fitted_functions) == 0L) {
      stop(
        "! No candidate fitted functions found. ",
        "Re-fit with compute_fitted = TRUE.",
        call. = FALSE
      )
    }
    
    all_fitted <- lapply(model$all_fitted_functions, function(v) v[[type]])
    all_fitted <- Filter(Negate(is.null), all_fitted)
    
    if (length(all_fitted) == 0L) {
      stop(
        paste0("! No fitted functions found for type = '", type, "'."),
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Resolve terms
  # ---------------------------------------------------------------------------
  if (is.null(terms)) {
    has_fit <- vapply(
      all_fitted,
      function(v) !is.null(v) && length(v) > 0L,
      logical(1L)
    )
    
    terms <- names(all_fitted)[has_fit]
    
    if (length(terms) == 0L) {
      if (is.null(type)) {
        message(
          "Nothing to plot: no variables in `cont_vars` were selected as\n",
          "having a significant interaction with the group variable.\n",
          "To inspect the fitted candidate curves, pass\n",
          "type = \"linear\", \"fp1\", or \"fp2\"."
        )
      } else {
        message(
          "Nothing to plot: no variables have a '", type,
          "' candidate fit.\n",
          "This can happen for flex levels that skip certain degrees,\n",
          "or when compute_fitted = FALSE."
        )
      }
      
      return(invisible(x))
    }
  } else {
    if (!is.character(terms) || length(terms) == 0L) {
      stop("! `terms` must be a non-empty character vector or NULL.",
           call. = FALSE)
    }
    
    missing_terms <- setdiff(terms, names(all_fitted))
    
    if (length(missing_terms) > 0L) {
      all_candidates <- model$all_fitted_functions
      
      not_selected <- if (!is.null(all_candidates)) {
        intersect(missing_terms, names(all_candidates))
      } else {
        character(0L)
      }
      
      truly_missing <- if (!is.null(all_candidates)) {
        setdiff(missing_terms, names(all_candidates))
      } else {
        missing_terms
      }
      
      if (length(truly_missing) > 0L) {
        available <- if (!is.null(all_candidates)) {
          paste(names(all_candidates), collapse = ", ")
        } else {
          "(none)"
        }
        
        stop(
          "! The following `terms` are not variables in this model:\n  ",
          paste(truly_missing, collapse = ", "), ".\n",
          "  Available variables: ",
          available,
          ".",
          call. = FALSE
        )
      }
      
      if (length(not_selected) > 0L) {
        message(
          "The following variable(s) were not selected as having a\n",
          "significant interaction: ",
          paste(not_selected, collapse = ", "), ".\n",
          "To plot the candidate curves anyway, pass\n",
          "type = \"linear\", \"fp1\", or \"fp2\"."
        )
        
        terms <- setdiff(terms, not_selected)
        
        if (length(terms) == 0L) {
          return(invisible(x))
        }
      }
    }
    
    no_fit <- vapply(
      terms,
      function(v) {
        fv <- all_fitted[[v]]
        is.null(fv) || length(fv) == 0L
      },
      logical(1L)
    )
    
    if (any(no_fit)) {
      bad_vars <- terms[no_fit]
      
      msg <- if (is.null(type)) {
        paste0(
          "! The following variables have no selected interaction fit:\n  ",
          paste(bad_vars, collapse = ", "), ".\n",
          "  Their interaction with the group variable did not meet the\n",
          "  selection threshold, so no winner curve was retained.\n",
          "  To plot the candidate fits anyway, pass `type = \"linear\"`,\n",
          "  `\"fp1\"`, or `\"fp2\"`."
        )
      } else {
        paste0(
          "! The following variables have no '", type, "' candidate fit:\n  ",
          paste(bad_vars, collapse = ", "), ".\n",
          "  This can happen for flex levels that skip certain degrees, or\n",
          "  when `compute_fitted = FALSE`."
        )
      }
      
      if (all(no_fit)) {
        stop(msg, call. = FALSE)
      } else {
        warning(
          paste0(
            msg,
            "\n  Continuing with the remaining variable(s): ",
            paste(terms[!no_fit], collapse = ", "),
            "."
          ),
          call. = FALSE
        )
        
        terms <- terms[!no_fit]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Determine metric display
  # ---------------------------------------------------------------------------
  crit_used <- if (is.null(model$criterion)) "pvalue" else model$criterion
  
  metric_spec <- switch(
    crit_used,
    pvalue = list(col = "pvalue",             label = "p"),
    aic    = list(col = "AIC_main_minus_int", label = "dAIC"),
    bic    = list(col = "BIC_main_minus_int", label = "dBIC"),
    list(col = "pvalue", label = "p")
  )
  
  padj      <- if (!is.null(model$p_adjust_method)) model$p_adjust_method else "none"
  adjusting <- padj != "none" && crit_used == "pvalue"
  
  subtitle_lookup <- list()
  
  if (!is.null(type)) {
    explore_label <- if (adjusting) "p_raw" else metric_spec$label
    all_metrics   <- model$all_model_metrics
    
    if (!is.null(all_metrics) && nrow(all_metrics) > 0L) {
      type_rows <- all_metrics[
        !is.na(all_metrics$type) & all_metrics$type == type,
        ,
        drop = FALSE
      ]
      
      for (i in seq_len(nrow(type_rows))) {
        v <- type_rows$variable[i]
        
        subtitle_lookup[[v]] <- list(
          type       = type,
          metric     = type_rows[[metric_spec$col]][i],
          label      = explore_label,
          metric_adj = NULL
        )
      }
    }
  } else if (!is.null(metrics) && nrow(metrics) > 0L) {
    adj_pvals <- NULL
    
    if (adjusting) {
      vw <- model$var_winners
      
      if (!is.null(vw)) {
        raw_pvals <- vapply(
          vw,
          function(w) {
            if (is.null(w$fit)) NA_real_ else w$metric$pvalue[1L]
          },
          numeric(1L)
        )
        
        adj_pvals <- stats::p.adjust(raw_pvals, method = padj)
        names(adj_pvals) <- names(vw)
      }
    }
    
    for (i in seq_len(nrow(metrics))) {
      v <- metrics$variable[i]
      
      subtitle_lookup[[v]] <- list(
        type       = metrics$type[i],
        metric     = metrics[[metric_spec$col]][i],
        label      = metric_spec$label,
        metric_adj = if (!is.null(adj_pvals) && v %in% names(adj_pvals)) {
          adj_pvals[[v]]
        } else {
          NULL
        }
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Group-level metadata
  # ---------------------------------------------------------------------------
  group_levels <- model$group_levels_original
  group_label  <- model$group_var
  
  if (is.null(group_levels) || length(group_levels) < 2L) {
    stop(
      "! `model$group_levels_original` must contain at least two group levels.",
      call. = FALSE
    )
  }
  
  ref_level <- group_levels[1L]
  non_ref   <- group_levels[-1L]
  
  # ---------------------------------------------------------------------------
  # Helper: subtitle
  # ---------------------------------------------------------------------------
  make_subtitle <- function(var) {
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
    
    lbl <- if (!is.null(info$label)) info$label else metric_spec$label
    
    val_str <- if (is.null(info$metric) || is.na(info$metric)) {
      "NA"
    } else if (lbl %in% c("p", "p_raw")) {
      formatC(info$metric, format = "g", digits = 3)
    } else {
      formatC(info$metric, format = "f", digits = 2)
    }
    
    base <- sprintf("%s  |  %s  |  %s = %s", var, type_str, lbl, val_str)
    
    if (!is.null(info$metric_adj) && !is.na(info$metric_adj)) {
      adj_str <- formatC(info$metric_adj, format = "g", digits = 3)
      base <- sprintf("%s (p_adj = %s)", base, adj_str)
    }
    
    base
  }
  
  # ---------------------------------------------------------------------------
  # Helper: fitted plot for one variable / group comparison
  # ---------------------------------------------------------------------------
  make_fitted_plot <- function(var, data_mat, grp_j) {
    col_ref    <- sprintf("f%s", ref_level)
    col_grp    <- sprintf("f%s", grp_j)
    col_ref_lo <- sprintf("f%s_lower", ref_level)
    col_ref_hi <- sprintf("f%s_upper", ref_level)
    col_grp_lo <- sprintf("f%s_lower", grp_j)
    col_grp_hi <- sprintf("f%s_upper", grp_j)
    
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
          " (grp ",
          grp_j,
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
    
    p <- ggplot2::ggplot()
    
    if (isTRUE(show_ci_fitted)) {
      p <- p +
        ggplot2::geom_ribbon(
          data = ref_df,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          fill  = colour_ref,
          alpha = ribbon_alpha
        )
    }
    
    p <- p +
      ggplot2::geom_line(
        data = ref_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        colour    = colour_ref,
        linewidth = linewidth
      )
    
    if (isTRUE(show_ci_fitted)) {
      p <- p +
        ggplot2::geom_ribbon(
          data = grp_df,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          fill  = colour_grp,
          alpha = ribbon_alpha
        )
    }
    
    p <- p +
      ggplot2::geom_line(
        data = grp_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        colour    = colour_grp,
        linewidth = linewidth
      ) +
      ggplot2::xlab(var) +
      ggplot2::ylab("Fitted linear predictor") +
      ggplot2::labs(
        title = if (show_title) {
          sprintf(
            "%s: group %s (red) vs group %s (blue)",
            group_label,
            grp_j,
            ref_level
          )
        } else {
          NULL
        },
        subtitle = if (show_title) make_subtitle(var) else NULL
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
  # Helper: difference plot for one variable / group comparison
  # ---------------------------------------------------------------------------
  make_diff_plot <- function(var, data_mat, grp_j) {
    col_diff  <- sprintf("f%s-f%s", grp_j, ref_level)
    col_lower <- sprintf("(f%s-f%s)_lower", grp_j, ref_level)
    col_upper <- sprintf("(f%s-f%s)_upper", grp_j, ref_level)
    
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
          " (grp ",
          grp_j,
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
      ggplot2::ylab(sprintf("f%s(x) - f%s(x)", grp_j, ref_level)) +
      ggplot2::labs(
        title = if (show_title) {
          sprintf(
            "%s: treatment effect (group %s - group %s)",
            group_label,
            grp_j,
            ref_level
          )
        } else {
          NULL
        },
        subtitle = if (show_title) make_subtitle(var) else NULL
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
      vector("list", length(non_ref)),
      sprintf("grp%s_vs_grp%s", non_ref, ref_level)
    )
    
    for (j in seq_along(non_ref)) {
      grp_j <- non_ref[j]
      
      p_out <- switch(
        plot_type,
        
        fitted = {
          make_fitted_plot(var, data_mat, grp_j)
        },
        
        difference = {
          make_diff_plot(var, data_mat, grp_j)
        },
        
        both = {
          p_fit  <- make_fitted_plot(var, data_mat, grp_j)
          p_diff <- make_diff_plot(var, data_mat, grp_j)
          
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