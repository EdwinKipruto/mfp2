# S3 plot method for mfpi objects
#
# plot.mfpi() produces three types of plots for each significant continuous
# variable in an mfpi fit:
#
#   "fitted"     — group-specific fitted curves with 95% CI bands
#   "difference" — pointwise treatment-effect difference f_j(x) - f_0(x)
#                  with 95% CI band and reference line at zero
#   "both"       — fitted and difference plots side by side (requires patchwork)
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
#' @param x An object of class \code{"mfpi"}, as returned by \code{mfpi()}.
#' @param terms Character vector of variable names to plot. Default \code{NULL}
#'   plots all variables with significant interactions.
#' @param type Character string or \code{NULL}. Controls which functional form
#'   to plot: \code{"linear"}, \code{"fp1"}, or \code{"fp2"}. Default
#'   \code{NULL} plots the \strong{best selected} form for each variable (from
#'   \code{model$best_fitted_functions}). Setting \code{type} explicitly plots
#'   that specific candidate from \code{model$all_fitted_functions} for all
#'   variables in \code{terms}, allowing direct comparison of functional forms
#'   regardless of which was selected.
#' @param plot_type Character string; \code{"fitted"} (default),
#'   \code{"difference"}, or \code{"both"}. See the \emph{Plot types} section.
#' @param auto_print Logical. If \code{TRUE} (default), each plot is printed
#'   to the active graphics device as it is created.
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
    if (is.null(all_fitted) || length(all_fitted) == 0L)
      stop("! No fitted functions found. Re-fit with compute_fitted = TRUE.",
           call. = FALSE)
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
    terms <- names(all_fitted)
  } else {
    missing_terms <- setdiff(terms, names(all_fitted))
    if (length(missing_terms) > 0L)
      stop(paste0("! Variables not found in fitted functions: ",
                  paste(missing_terms, collapse = ", "), "."),
           call. = FALSE)
  }
  
  # Build a lookup: variable -> (type, pvalue) for subtitle
  # When type is explicit, use that; otherwise use the selected type per variable
  subtitle_lookup <- list()
  if (!is.null(metrics) && nrow(metrics) > 0L) {
    for (i in seq_len(nrow(metrics))) {
      v <- metrics$variable[i]
      subtitle_lookup[[v]] <- list(
        type   = if (!is.null(type)) type else metrics$type[i],
        pvalue = metrics$pvalue[i]
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
    sprintf("%s  |  %s  |  p = %s", var, type_str, info$pvalue)
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
        title    = sprintf("%s: group %s (red) vs group %s (blue)",
                           group_label, grp_j, ref_level),
        subtitle = make_subtitle(var)
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
        title    = sprintf("%s: treatment effect (group %s - group %s)",
                           group_label, grp_j, ref_level),
        subtitle = make_subtitle(var)
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