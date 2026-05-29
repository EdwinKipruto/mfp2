# S3 plot method for mfpi objects
#
# plot.mfpi() produces three types of plots for each significant continuous
# variable in an mfpi fit:
#
#   "fitted"     - group-specific fitted curves with 95% CI bands
#   "difference" - pointwise treatment-effect difference f_j(x) - f_0(x)
#                  with 95% CI band and reference line at zero
#   "both"       - fitted and difference plots side by side (requires patchwork)
#
# All plots include rug marks for the observed values of the continuous
# variable, a subtitle showing the variable name / functional form / p-value,
# and are returned as a nested named list that can be printed, saved, or
# assembled externally.


# -----------------------------------------------------------------------------
# plot.mfpi() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Plot Fitted Functions and Treatment-Effect Differences from an MFPI Model
#'
#' Produces group-specific fitted-function plots, treatment-effect difference
#' plots, or both for each significant continuous variable in an
#' \code{"mfpi"} object. Each plot includes 95\% confidence bands and rug
#' marks. All plots are printed automatically and returned invisibly as a
#' nested named list.
#'
#' @section Plot types:
#' Three plot types are available via the \code{plot_type} argument.
#' \code{"fitted"} overlays the fitted linear predictor for the reference
#' group and one non-reference group on shared axes with 95\% CI bands.
#' \code{"difference"} plots the pointwise treatment-effect function
#' \eqn{D_j(x) = f_j(x) - f_0(x)} with a 95\% CI band and a horizontal
#' reference line at zero; deviations from zero indicate effect modification.
#' \code{"both"} shows fitted and difference plots side by side (requires
#' \pkg{patchwork}).
#'
#' @section Return structure:
#' The return value is a nested named list. Outer names are variable names
#' from \code{terms}. Inner names are of the form
#' \code{"grp<j>_vs_grp<0>"} for each non-reference group level \eqn{j}.
#' Each leaf is a \code{ggplot} object (\code{plot_type = "fitted"} or
#' \code{"difference"}) or a \code{patchwork} object
#' (\code{plot_type = "both"}).
#'
#' @section Title and subtitle annotations:
#' Each plot is annotated by default with two text lines:
#'
#' The \strong{title} describes the structural comparison being plotted:
#' \itemize{
#'   \item Fitted plots: \code{"<group_var>: group <j> (red) vs group <0> (blue)"}.
#'   \item Difference plots: \code{"<group_var>: treatment effect (group <j> - group <0>)"}.
#' }
#'
#' The \strong{subtitle} summarises the variable, the functional form being
#' shown, and the selection metric for that form, in the format
#' \code{"<variable>  |  <type>  |  <metric> = <value>"}. The metric shown
#' depends on the \code{criterion} used to fit the model:
#' \itemize{
#'   \item \code{criterion = "pvalue"} (default): shows \code{p = <pvalue>},
#'     the p-value of the interaction test.
#'   \item \code{criterion = "aic"}: shows \code{dAIC = <diff>}, where
#'     \code{dAIC} is \eqn{\mathrm{AIC}_{\mathrm{main}} -
#'     \mathrm{AIC}_{\mathrm{interaction}}}. Positive values favour the
#'     interaction model.
#'   \item \code{criterion = "bic"}: shows \code{dBIC = <diff>}, defined
#'     analogously.
#' }
#'
#' When \code{type = NULL} (default), the subtitle reports the metric of
#' the \strong{selected} (winning) functional form. When \code{type} is
#' explicit, the subtitle reports the metric of the \strong{requested}
#' candidate (linear, FP1, or FP2), even if a different form was the
#' winner. This keeps the printed metric consistent with the curve being
#' shown.
#'
#' To suppress both title and subtitle for cleaner plots in reports or
#' manuscripts, set \code{show_title = FALSE}. Axis labels (\code{xlab},
#' \code{ylab}) are unaffected; they always describe what the axes mean.
#'
#' @section Plotting non-selected variables:
#' Whether a variable can be plotted depends on \code{type}:
#'
#' Under \code{type = NULL} (the default), only variables whose interaction
#' was retained as significant have a "winner" curve stored in
#' \code{best_fitted_functions}. If you call \code{plot(model)} with no
#' \code{terms}, those non-selected variables are silently dropped and you
#' see plots only for the variables that survived selection. If you name a
#' non-selected variable explicitly in \code{terms}, you get an informative
#' error pointing you to set \code{type} to a specific candidate.
#'
#' Under \code{type = "linear"}, \code{"fp1"}, or \code{"fp2"}, every
#' variable in \code{cont_vars} has a candidate curve in
#' \code{all_fitted_functions} (provided the model was fitted with
#' \code{compute_fitted = TRUE}), so you can plot the candidate fit even
#' for variables that were not selected. This is useful for diagnostic
#' inspection: you can see what each functional form looked like before
#' the selection step and judge whether the rejection was supported.
#'
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param terms Character vector of variable names to plot. Default
#'   \code{NULL} plots all variables that have a fitted curve available
#'   under the current \code{type} setting:
#'   when \code{type = NULL}, that is all variables with a significant
#'   interaction (i.e. those retained in \code{best_fitted_functions});
#'   when \code{type} is explicit, that is all variables with that
#'   candidate fit (typically all of \code{cont_vars}). If you name
#'   variables that have no fit under the current \code{type}, you get
#'   a clear message explaining why (see \emph{Plotting non-selected
#'   variables}).
#' @param type Character string or \code{NULL}. Controls which functional
#'   form to plot: \code{"linear"}, \code{"fp1"}, or \code{"fp2"}. Default
#'   \code{NULL} plots the \strong{best selected} form for each variable
#'   (from \code{model$best_fitted_functions}); variables whose interaction
#'   was not retained are silently skipped under this default. Setting
#'   \code{type} explicitly plots that specific candidate from
#'   \code{model$all_fitted_functions} for every variable in
#'   \code{cont_vars}, allowing the user to inspect a candidate curve
#'   regardless of whether it was the winner.
#' @param plot_type Character string; \code{"fitted"} (default),
#'   \code{"difference"}, or \code{"both"}. See the \emph{Plot types} section.
#' @param auto_print Logical. If \code{TRUE} (default), each plot is printed
#'   to the active graphics device as it is created.
#' @param show_title Logical. If \code{TRUE} (default), each plot is
#'   annotated with a title and subtitle as described in the
#'   \emph{Title and subtitle annotations} section. Set to \code{FALSE}
#'   for clean plots without text annotations, useful when embedding in
#'   reports or manuscripts where the caption already describes the
#'   figure.
#' @param colour_ref Character. Colour for the reference group curve and CI.
#'   Default \code{"#2166AC"} (blue).
#' @param colour_grp Character. Colour for the non-reference group curve and
#'   CI. Default \code{"#D6604D"} (red-orange).
#' @param colour_diff Character. Colour for the difference curve and CI band.
#'   Default \code{"#1B7837"} (green).
#' @param linewidth Positive numeric. Line width. Default \code{0.8}.
#' @param ribbon_alpha Numeric in \eqn{[0,1]}. Transparency of CI ribbon.
#'   Default \code{0.2}.
#' @param rug_alpha Numeric in \eqn{[0,1]}. Transparency of rug marks.
#'   Default \code{0.4}.
#' @param ... Currently unused.
#'
#' @return A nested named list of \code{ggplot} (or \code{patchwork}) objects,
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
#' fit   <- mfpi(x, y, group_var = "trt", cont_vars = c("age", "bmi"))
#' plots <- plot(fit)                        # all variables, fitted + printed
#' plots <- plot(fit, plot_type = "both")    # fitted + difference side by side
#' plots <- plot(fit, terms = "age",
#'               plot_type = "difference")   # difference only for age
#' # Save a specific plot
#' ggplot2::ggsave("age_diff.pdf", plots$age$grp1_vs_grp0, width = 7, height = 5)
#' }
plot.mfpi <- function(x,
                      terms        = NULL,
                      type         = NULL,
                      plot_type    = c("fitted", "difference", "both"),
                      auto_print   = TRUE,
                      show_title   = TRUE,
                      colour_ref   = "#2166AC",
                      colour_grp   = "#D6604D",
                      colour_diff  = "#1B7837",
                      linewidth    = 0.8,
                      ribbon_alpha = 0.2,
                      rug_alpha    = 0.4,
                      ...) {
  
  model     <- x
  plot_type <- match.arg(plot_type)
  
  if (!inherits(model, "mfpi"))
    stop("! `x` must be an object of class \"mfpi\".", call. = FALSE)
  
  if (plot_type == "both" && !requireNamespace("patchwork", quietly = TRUE))
    stop("! `plot_type = \"both\"` requires the patchwork package. ",
         "Install it with install.packages(\"patchwork\").", call. = FALSE)
  
  # ---------------------------------------------------------------------------
  # Extract fitted functions and metrics
  # ---------------------------------------------------------------------------
  metrics <- model$best_model_metrics
  
  # Resolve which fitted functions to use based on type argument:
  #   NULL (default) -> winner fitted functions from model$best_fitted_functions
  #   "linear"/"fp1"/"fp2" -> specific candidate from model$all_fitted_functions
  valid_types <- c("linear", "fp1", "fp2")
  
  if (is.null(type)) {
    # Default: use the selected (best) functional form per variable
    all_fitted <- model$best_fitted_functions
    if (is.null(all_fitted))
      stop("! No fitted functions found. Re-fit with compute_fitted = TRUE.",
           call. = FALSE)
    # If all_fitted is an empty list, no variable was selected -- this is
    # handled below (line "if (length(terms) == 0L)") with a more informative
    # message that suggests using type = "linear"/"fp1"/"fp2".
  } else {
    type <- match.arg(type, valid_types)
    if (is.null(model$all_fitted_functions) || length(model$all_fitted_functions) == 0L)
      stop("! No candidate fitted functions found. Re-fit with compute_fitted = TRUE.",
           call. = FALSE)
    # Extract fitted functions for the requested type across all variables
    all_fitted <- lapply(model$all_fitted_functions, function(v) v[[type]])
    # Drop variables where this type has no fitted values
    all_fitted <- Filter(Negate(is.null), all_fitted)
    if (length(all_fitted) == 0L)
      stop(paste0("! No fitted functions found for type = '", type, "'."),
           call. = FALSE)
  }
  
  # Resolve terms
  if (is.null(terms)) {
    # Default to all variables that have fitted values; drop those that don't.
    # When `type` is NULL this is the set of significant variables; when
    # `type` is explicit this is all variables (since all have candidates).
    has_fit <- vapply(all_fitted, function(v) !is.null(v) && length(v) > 0L,
                      logical(1L))
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
    # User asked for specific variables: validate names first.
    # When type = NULL, all_fitted is best_fitted_functions (winners only).
    # A variable may be missing because it was not selected (not significant)
    # rather than because it does not exist. Distinguish the two cases.
    missing_terms <- setdiff(terms, names(all_fitted))
    
    if (length(missing_terms) > 0L) {
      # Check if the missing terms exist in all_fitted_functions (candidates)
      all_candidates <- model$all_fitted_functions
      not_selected   <- intersect(missing_terms,  names(all_candidates))
      truly_missing  <- setdiff(missing_terms, names(all_candidates))
      
      if (length(truly_missing) > 0L) {
        stop(
          "! The following `terms` are not variables in this model:\n  ",
          paste(truly_missing, collapse = ", "), ".\n",
          "  Available variables: ",
          paste(names(all_candidates), collapse = ", "), ".",
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
        # Remove non-selected terms and continue with any remaining
        terms <- setdiff(terms, not_selected)
        if (length(terms) == 0L) return(invisible(x))
      }
    }
    
    # Then check which (if any) lack fitted values for the chosen type
    no_fit <- vapply(terms, function(v) {
      fv <- all_fitted[[v]]
      is.null(fv) || length(fv) == 0L
    }, logical(1L))
    
    if (any(no_fit)) {
      bad_vars <- terms[no_fit]
      msg <- if (is.null(type)) {
        paste0(
          "! The following variables have no selected interaction fit:\n  ",
          paste(bad_vars, collapse = ", "), ".\n",
          "  Their interaction with the group variable did not meet the\n",
          "  selection threshold, so no 'winner' curve was retained.\n",
          "  To plot the candidate fits anyway, pass `type = \"linear\"`,\n",
          "  `\"fp1\"`, or `\"fp2\"` to see what each functional form looked\n",
          "  like before the selection step."
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
            msg, "\n  Continuing with the remaining variable(s): ",
            paste(terms[!no_fit], collapse = ", "), "."
          ),
          call. = FALSE
        )
        terms <- terms[!no_fit]
      }
    }
  }
  
  # Determine which metric to display in the subtitle based on the criterion
  # used to fit the model. Stored on `model$criterion` by fit_mfpi().
  crit_used <- if (is.null(model$criterion)) "pvalue" else model$criterion
  
  # Map criterion -> (metric_col, display_label)
  metric_spec <- switch(crit_used,
                        "pvalue" = list(col = "pvalue",             label = "p"),
                        "aic"    = list(col = "AIC_main_minus_int", label = "dAIC"),
                        "bic"    = list(col = "BIC_main_minus_int", label = "dBIC"),
                        list(col = "pvalue", label = "p")   # fallback
  )
  
  # Build a lookup: variable -> (type, metric_value) for subtitle.
  # - When `type` is NULL: use the per-variable winner from best_model_metrics
  #   so the subtitle reflects what was selected.
  # - When `type` is explicit: pull from all_model_metrics matched on
  #   (variable, type), so the subtitle matches the curve being plotted.
  subtitle_lookup <- list()
  if (!is.null(type)) {
    all_metrics <- model$all_model_metrics
    if (!is.null(all_metrics) && nrow(all_metrics) > 0L) {
      type_rows <- all_metrics[!is.na(all_metrics$type) &
                                 all_metrics$type == type, , drop = FALSE]
      for (i in seq_len(nrow(type_rows))) {
        v <- type_rows$variable[i]
        subtitle_lookup[[v]] <- list(
          type     = type,
          metric   = type_rows[[metric_spec$col]][i]
        )
      }
    }
  } else if (!is.null(metrics) && nrow(metrics) > 0L) {
    for (i in seq_len(nrow(metrics))) {
      v <- metrics$variable[i]
      subtitle_lookup[[v]] <- list(
        type     = metrics$type[i],
        metric   = metrics[[metric_spec$col]][i]
      )
    }
  }
  
  # Group-level metadata
  group_levels <- model$group_levels_original
  group_label  <- model$group_var
  ref_level    <- group_levels[1L]
  non_ref      <- group_levels[-1L]
  
  # ---------------------------------------------------------------------------
  # Helper: build subtitle string
  # ---------------------------------------------------------------------------
  make_subtitle <- function(var) {
    info <- subtitle_lookup[[var]]
    if (is.null(info)) return(var)
    type_str <- switch(info$type,
                       linear = "Linear",
                       fp1    = "FP1",
                       fp2    = "FP2",
                       info$type
    )
    val_str <- if (is.null(info$metric) || is.na(info$metric))
      "NA"
    else if (metric_spec$label == "p")
      formatC(info$metric, format = "g", digits = 3)
    else
      formatC(info$metric, format = "f", digits = 2)
    sprintf("%s  |  %s  |  %s = %s",
            var, type_str, metric_spec$label, val_str)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: build fitted plot for one variable / one group comparison
  # ---------------------------------------------------------------------------
  make_fitted_plot <- function(var, data_mat, grp_j) {
    col_ref    <- sprintf("f%s",       ref_level)
    col_grp    <- sprintf("f%s",       grp_j)
    col_ref_lo <- sprintf("f%s_lower", ref_level)
    col_ref_hi <- sprintf("f%s_upper", ref_level)
    col_grp_lo <- sprintf("f%s_lower", grp_j)
    col_grp_hi <- sprintf("f%s_upper", grp_j)
    
    needed <- c(var, col_ref, col_grp, col_ref_lo, col_ref_hi,
                col_grp_lo, col_grp_hi)
    missing <- setdiff(needed, names(data_mat))
    if (length(missing) > 0L) {
      warning(paste0("! Skipping fitted plot for ", var,
                     " (grp ", grp_j, "): missing columns: ",
                     paste(missing, collapse = ", ")),
              call. = FALSE)
      return(NULL)
    }
    
    xvec      <- data_mat[[var]]
    ref_df    <- data.frame(x = xvec,
                            y = data_mat[[col_ref]],
                            lo = data_mat[[col_ref_lo]],
                            hi = data_mat[[col_ref_hi]])
    grp_df    <- data.frame(x = xvec,
                            y = data_mat[[col_grp]],
                            lo = data_mat[[col_grp_lo]],
                            hi = data_mat[[col_grp_hi]])
    rug_df    <- data.frame(x = xvec)
    
    ggplot2::ggplot() +
      # Reference group CI ribbon and curve
      ggplot2::geom_ribbon(
        data = ref_df,
        ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
        fill = colour_ref, alpha = ribbon_alpha
      ) +
      ggplot2::geom_line(
        data = ref_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        colour = colour_ref, linewidth = linewidth
      ) +
      # Non-reference group CI ribbon and curve
      ggplot2::geom_ribbon(
        data = grp_df,
        ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
        fill = colour_grp, alpha = ribbon_alpha
      ) +
      ggplot2::geom_line(
        data = grp_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        colour = colour_grp, linewidth = linewidth
      ) +
      # Rug marks
      ggplot2::geom_rug(
        data = rug_df,
        ggplot2::aes(x = .data$x),
        sides = "b", alpha = rug_alpha
      ) +
      ggplot2::xlab(var) +
      ggplot2::ylab("Fitted linear predictor") +
      ggplot2::labs(
        title    = if (show_title)
          sprintf("%s: group %s (red) vs group %s (blue)",
                  group_label, grp_j, ref_level)
        else NULL,
        subtitle = if (show_title) make_subtitle(var) else NULL
      ) +
      ggplot2::theme_bw()
  }
  
  # ---------------------------------------------------------------------------
  # Helper: build difference plot for one variable / one group comparison
  # ---------------------------------------------------------------------------
  make_diff_plot <- function(var, data_mat, grp_j) {
    col_diff  <- sprintf("f%s-f%s",         grp_j, ref_level)
    col_lower <- sprintf("(f%s-f%s)_lower", grp_j, ref_level)
    col_upper <- sprintf("(f%s-f%s)_upper", grp_j, ref_level)
    
    needed  <- c(var, col_diff, col_lower, col_upper)
    missing <- setdiff(needed, names(data_mat))
    if (length(missing) > 0L) {
      warning(paste0("! Skipping difference plot for ", var,
                     " (grp ", grp_j, "): missing columns: ",
                     paste(missing, collapse = ", ")),
              call. = FALSE)
      return(NULL)
    }
    
    plot_df <- data.frame(
      x     = data_mat[[var]],
      diff  = data_mat[[col_diff]],
      lower = data_mat[[col_lower]],
      upper = data_mat[[col_upper]]
    )
    
    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        fill = colour_diff, alpha = ribbon_alpha
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey50") +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$diff),
        colour = colour_diff, linewidth = linewidth
      ) +
      ggplot2::geom_rug(
        ggplot2::aes(x = .data$x),
        sides = "b", alpha = rug_alpha
      ) +
      ggplot2::xlab(var) +
      ggplot2::ylab(sprintf("f%s(x) - f%s(x)", grp_j, ref_level)) +
      ggplot2::labs(
        title    = if (show_title)
          sprintf("%s: treatment effect (group %s - group %s)",
                  group_label, grp_j, ref_level)
        else NULL,
        subtitle = if (show_title) make_subtitle(var) else NULL
      ) +
      ggplot2::theme_bw()
  }
  
  # ---------------------------------------------------------------------------
  # Build all plots
  # ---------------------------------------------------------------------------
  plots <- stats::setNames(vector("list", length(terms)), terms)
  
  for (var in terms) {
    data_mat  <- as.data.frame(all_fitted[[var]])
    grp_plots <- stats::setNames(
      vector("list", length(non_ref)),
      sprintf("grp%s_vs_grp%s", non_ref, ref_level)
    )
    
    for (j in seq_along(non_ref)) {
      grp_j <- non_ref[j]
      
      p_fit  <- make_fitted_plot(var, data_mat, grp_j)
      p_diff <- make_diff_plot(var, data_mat, grp_j)
      
      p_out <- switch(plot_type,
                      fitted     = p_fit,
                      difference = p_diff,
                      both       = if (!is.null(p_fit) && !is.null(p_diff))
                        patchwork::wrap_plots(p_fit, p_diff, ncol = 2L)
                      else
                        if (!is.null(p_fit)) p_fit else p_diff
      )
      
      grp_plots[[j]] <- p_out
      if (isTRUE(auto_print) && !is.null(p_out)) print(p_out)
    }
    
    plots[[var]] <- grp_plots
  }
  
  invisible(plots)
}