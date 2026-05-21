# S3 plot method for mfpi objects
#
# plot.mfpi() dispatches on class "mfpi" so users call it as plot(model) or
# plot(model, terms = "age", plot_type = "difference"). It returns a named list
# of ggplot objects that can be printed, saved, or assembled with patchwork.


# -----------------------------------------------------------------------------
# plot.mfpi() -----------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Plot Fitted Functions and Treatment-Effect Differences from an MFPI Model
#'
#' Produces group-specific fitted-function plots or treatment-effect difference
#' plots for each continuous variable in an `"mfpi"` object. For each
#' continuous variable and each non-reference group, one [ggplot2::ggplot()]
#' object is created and returned in a named list.
#'
#' @section Plot types:
#' **`"fitted"`** overlays the fitted linear predictor \eqn{\hat{f}_j(x)} for
#' the reference group \eqn{j = 0} and one non-reference group \eqn{j} on a
#' shared set of axes. The two curves are distinguished by colour.
#'
#' **`"difference"`** plots the pointwise treatment-effect function
#' \deqn{
#'   \hat{D}_j(x) = \hat{f}_j(x) - \hat{f}_0(x)
#' }
#' together with a pointwise 95% confidence band
#' \deqn{
#'   \hat{D}_j(x) \;\pm\; 1.96\;\widehat{\mathrm{SE}}\!\bigl[\hat{D}_j(x)\bigr].
#' }
#' A horizontal reference line at zero is included; deviations from zero
#' indicate effect modification by the continuous covariate.
#'
#' @section Return structure:
#' The return value is a nested named list:
#' \itemize{
#'   \item Outer names: variable names from `terms`.
#'   \item Inner names: `"grp<j>_vs_grp<0>"` for each non-reference group
#'     level \eqn{j}.
#'   \item Each leaf is a [ggplot2::ggplot()] object.
#' }
#' To display a specific plot:
#' ```r
#' plots <- plot(model)
#' print(plots$age$grp1_vs_grp0)
#' ```
#'
#' @param x An object of class `"mfpi"`, as returned by [mfp2::mfpi()].
#' @param terms Character vector of variable names to plot. Default `NULL`
#'   plots all variables with significant interactions in `x$fitted_functions`.
#' @param interaction_type Character string; `"linear"`, `"fp1"`, or `"fp2"`.
#'   Selects which interaction type to visualise. The function checks that the
#'   requested type was actually fitted for each variable and warns if it was
#'   not.
#' @param plot_type Character string; `"fitted"` (default) or `"difference"`.
#'   See the *Plot types* section.
#' @param linewidth Positive numeric. Line width passed to
#'   [ggplot2::geom_line()]. Default `1`.
#' @param linetype Integer or character. Line type passed to
#'   [ggplot2::geom_line()]. Default `1` (solid).
#' @param line_colour Character string. Colour of the difference curve (used
#'   only when `plot_type = "difference"`). Default `"black"`.
#' @param rug_colour Character string. Colour of the rug marks along the
#'   x-axis. Default `"black"`.
#' @param ... Currently unused. Reserved for future extensions.
#'
#' @return A nested named list of [ggplot2::ggplot()] objects (see *Return
#'   structure*). Invisibly returned; call [print()] or
#'   [ggplot2::ggsave()] on individual elements.
#'
#' @seealso [mfp2::mfpi()], `gen_fitted_values_per_group()`
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' \dontrun{
#' fit    <- mfpi(x, y, group_var = "trt", cont_vars = c("age", "bmi"))
#' plots  <- plot(fit)
#' print(plots$age$grp1_vs_grp0)
#'
#' # Difference plots only for one variable
#' plot(fit, terms = "age", plot_type = "difference")
#' }
plot.mfpi <- function(x,
                      terms            = NULL,
                      interaction_type = c("linear", "fp1", "fp2"),
                      plot_type        = c("fitted", "difference"),
                      linewidth        = 1,
                      linetype         = 1,
                      line_colour      = "black",
                      rug_colour       = "black",
                      ...) {
  
  # Rename for clarity; `x` is the S3 convention but `model` reads better
  model <- x
  
  if (!inherits(model, "mfpi")) {
    stop("! `x` must be an object of class \"mfpi\".", call. = FALSE)
  }
  
  interaction_type <- match.arg(interaction_type)
  plot_type        <- match.arg(plot_type)
  
  # ---------------------------------------------------------------------------
  # Extract fitted functions for the requested interaction type
  # ---------------------------------------------------------------------------
  # fitted_functions is a flat named list: model$fitted_functions[[var_name]]
  # Each element is the matrix returned by gen_fitted_values_per_group().
  # The `type` column in model_evaluation_metrics records which interaction
  # type was selected for each variable; we filter to the requested type.
  # ---------------------------------------------------------------------------
  all_fitted <- model$univariable_interactions$fitted_functions
  
  # Determine which variables had the requested interaction type selected
  metrics <- model$univariable_interactions$model_evaluation_metrics
  
  if (is.null(all_fitted) || length(all_fitted) == 0L) {
    stop(
      "! No fitted functions found in the model. ",
      "Ensure `compute_fitted = TRUE` was used when calling `mfpi()`.",
      call. = FALSE
    )
  }
  
  # Filter variables by interaction type
  if (!is.null(metrics) && nrow(metrics) > 0L && "type" %in% names(metrics)) {
    vars_of_type <- metrics$variable[metrics$type == interaction_type]
    fitted       <- all_fitted[intersect(names(all_fitted), vars_of_type)]
    if (length(fitted) == 0L) {
      stop(
        paste0(
          "! No significant interactions of type '", interaction_type,
          "' were found. Available types: ",
          paste(unique(metrics$type), collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  } else {
    # Fallback: use all fitted functions (no type column available)
    fitted <- all_fitted
  }
  
  # Validate and resolve terms -------------------------------------------------
  available_vars <- names(fitted)
  
  if (is.null(terms)) {
    terms <- available_vars
  } else {
    missing_terms <- setdiff(terms, available_vars)
    if (length(missing_terms) > 0L) {
      stop(
        paste0(
          "! The following terms were not found in the fitted functions for ",
          "interaction_type = '", interaction_type, "': ",
          paste(missing_terms, collapse = ", "), ".\n",
          "i Available variables: ", paste(available_vars, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  }
  
  # Group-level metadata -------------------------------------------------------
  group_levels <- model$univariable_interactions$group_levels_original
  group_label  <- model$univariable_interactions$group_var
  ref_level    <- group_levels[1L]
  non_ref      <- group_levels[-1L]
  
  # ---------------------------------------------------------------------------
  # Build plots
  # ---------------------------------------------------------------------------
  plots <- stats::setNames(vector("list", length(terms)), terms)
  
  for (var in terms) {
    data_mat  <- as.data.frame(fitted[[var]])
    grp_plots <- stats::setNames(
      vector("list", length(non_ref)),
      sprintf("grp%s_vs_grp%s", non_ref, ref_level)
    )
    
    for (j in seq_along(non_ref)) {
      grp_j <- non_ref[j]
      
      if (plot_type == "difference") {
        # Columns: <var>, fj-f0, (fj-f0)_lower, (fj-f0)_upper
        col_diff  <- sprintf("f%s-f%s",       grp_j, ref_level)
        col_lower <- sprintf("(f%s-f%s)_lower", grp_j, ref_level)
        col_upper <- sprintf("(f%s-f%s)_upper", grp_j, ref_level)
        
        missing <- setdiff(c(var, col_diff, col_lower, col_upper),
                           names(data_mat))
        if (length(missing) > 0L) {
          warning(
            paste0("! Skipping ", var, " (grp ", grp_j, " vs ", ref_level,
                   "): missing columns: ", paste(missing, collapse = ", ")),
            call. = FALSE
          )
          next
        }
        
        plot_df <- data.frame(
          variable = data_mat[[var]],
          diff     = data_mat[[col_diff]],
          lower    = data_mat[[col_lower]],
          upper    = data_mat[[col_upper]]
        )
        
        grp_plots[[j]] <- ggplot2::ggplot(
          plot_df,
          ggplot2::aes(x = .data$variable, y = .data$diff)
        ) +
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
            fill = "blue", alpha = 0.2
          ) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                              colour = "grey50") +
          ggplot2::geom_line(linewidth = linewidth,
                             linetype = linetype,
                             colour   = line_colour) +
          ggplot2::geom_rug(ggplot2::aes(x = .data$variable),
                            sides = "b", colour = rug_colour) +
          ggplot2::xlab(var) +
          ggplot2::ylab(sprintf("f%s(x) - f%s(x)", grp_j, ref_level)) +
          ggplot2::labs(
            title = sprintf("%s: group %s vs group %s",
                            group_label, grp_j, ref_level)
          ) +
          ggplot2::theme_bw()
        
      } else {
        # plot_type == "fitted"
        # Columns: <var>, f<ref>, f<j>
        col_ref <- sprintf("f%s", ref_level)
        col_grp <- sprintf("f%s", grp_j)
        
        missing <- setdiff(c(var, col_ref, col_grp), names(data_mat))
        if (length(missing) > 0L) {
          warning(
            paste0("! Skipping ", var, " (grp ", grp_j, " vs ", ref_level,
                   "): missing columns: ", paste(missing, collapse = ", ")),
            call. = FALSE
          )
          next
        }
        
        # Reshape to long format without data.table or dplyr::recode
        plot_wide <- data.frame(
          variable = data_mat[[var]],
          ref      = data_mat[[col_ref]],
          grp      = data_mat[[col_grp]]
        )
        plot_long <- reshape(
          plot_wide,
          varying       = c("ref", "grp"),
          v.names       = "fitted",
          timevar       = "group",
          times         = c(as.character(ref_level), as.character(grp_j)),
          direction     = "long",
          new.row.names = seq_len(2L * nrow(plot_wide))
        )
        
        grp_plots[[j]] <- ggplot2::ggplot(
          plot_long,
          ggplot2::aes(
            x      = .data$variable,
            y      = .data$fitted,
            colour = factor(.data$group),
            group  = factor(.data$group)
          )
        ) +
          ggplot2::geom_line(linewidth = linewidth,
                             linetype = linetype) +
          ggplot2::geom_rug(ggplot2::aes(x = .data$variable),
                            sides = "b", colour = rug_colour) +
          ggplot2::scale_colour_discrete(
            name   = group_label,
            labels = c(as.character(ref_level), as.character(grp_j))
          ) +
          ggplot2::xlab(var) +
          ggplot2::ylab("Fitted linear predictor") +
          ggplot2::theme_bw()
      }
    }
    
    plots[[var]] <- grp_plots
  }
  
  invisible(plots)
}