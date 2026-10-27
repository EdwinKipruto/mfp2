#' Plot MFPI Fitted Functions and Group Differences
#'
#' Plot group-specific fitted functions and differences between groups for
#' variables tested for interaction by [mfpi()]. Continuous predictors are
#' shown as curves; binary and categorical predictors are shown as points
#' at each fitted level.
#'
#' @param x An object of class `"mfpi"` returned by [mfpi()].
#' @param terms An optional character vector naming variables tested for
#'   interaction. The default is `NULL`, which plots selected interactions
#'   only. Specify a variable to plot its interaction even if it was not
#'   selected.
#' @param plot_type A character string specifying the plots to display:
#'   `"fitted"` for group-specific functions, `"difference"` for the
#'   comparison-group function minus the reference-group function, or
#'   `"both"` for side-by-side plots. Default: `"both"`, which requires
#'   the optional \pkg{patchwork} package.
#' @param auto_print A logical value indicating whether to display plots
#'   automatically. Default: `TRUE`. Use `FALSE` to return the plots
#'   without displaying them.
#' @param show_title A logical value indicating whether to display titles
#'   and subtitles. Default: `TRUE`.
#' @param show_rug A logical value indicating whether to show observed
#'   continuous predictor values as rug marks along the horizontal axis.
#'   Default: `FALSE`.
#' @param show_ci_fitted A logical value indicating whether to display
#'   confidence intervals for group-specific functions. Default: `FALSE`.
#' @param show_ci_diff A logical value indicating whether to display
#'   confidence intervals for group differences. Default: `TRUE`.
#' @param show_null_line A logical value indicating whether to add a
#'   horizontal line at zero to difference plots. Default: `FALSE`.
#' @param show_maineffect_line A logical value indicating whether to add
#'   the constant group difference estimated by the model without interaction
#'   to difference plots. Default: `TRUE`. The line is omitted if that
#'   model is unavailable.
#' @param show_ci_maineffect A logical value indicating whether to display
#'   a confidence band around the main-effects reference line.
#'   Requires `show_maineffect_line = TRUE` and an available standard error.
#'   Default: `FALSE`.
#' @param legend_position A character string specifying the legend position:
#'   `"inside"` (default), `"right"`, `"left"`, `"bottom"`, `"top"`, or `"none"`.
#' @param legend_inside Two numeric values giving the horizontal and vertical
#'   position of an inside legend, between 0 and 1. The bottom-left corner
#'   is `c(0, 0)` and the top-right corner is `c(1, 1)`.
#'   Default: `c(0.02, 0.98)`.
#' @param legend_justification Two numeric values specifying which part
#'   of an inside legend is anchored at `legend_inside`.
#'   Default: `c(0, 1)`, which anchors the legend at its top-left corner.
#' @param colour_ref Colour of the reference-group curve, points, and
#'   confidence intervals. Default: `"#2166AC"`.
#' @param colour_group Colour of the comparison-group curve, points, and
#'   confidence intervals. Default: `"#D6604D"`.
#' @param line_types An optional character vector specifying curve styles.
#'   For example, `c("solid", "dashed")` gives the reference group a solid
#'   line and the comparison group a dashed line. Alternatively, assign styles
#'   to group labels, for example `c(Placebo = "solid", Treatment = "dashed")`,
#'   using the labels in your data. Default: `NULL`, which makes all lines solid.
#' @param colour_diff Colour of the difference curve, points, and confidence
#'   intervals. Default: `"#1B7837"`.
#' @param colour_null Colour of the zero reference line. Default: `"grey50"`.
#' @param colour_maineffect Colour of the main-effects reference line.
#'   Default: `"grey20"`.
#' @param linewidth A positive number controlling curve width. Default: `1`.
#' @param ribbon_alpha Confidence-band transparency, from 0 (transparent)
#'   to 1 (opaque). Default: `0.2`.
#' @param rug_alpha Rug-mark transparency, from 0 (transparent) to 1 (opaque).
#'   Default: `0.4`.
#' @param ... Not used. Supplied arguments generate a warning.
#'
#' @section Interpreting the plots:
#' Each variable is plotted using its separate fitted interaction model.
#' The reference group is the one used during fitting. With more than two
#' groups, a separate plot is produced for each comparison with the reference.
#' Multinomial models also produce separate plots for each non-reference
#' response category.
#'
#' Fitted functions include the predictor contribution, group main effects,
#' and model intercept where applicable. They exclude contributions from
#' other adjustment variables and offsets. Group differences include both
#' the group main-effect difference and differences between the predictor
#' functions (see Kipruto, Royston and Sauerbrei (2026) for details).
#'
#' Plots use the model's linear-predictor scale. For example, Gaussian
#' models with an identity link use the outcome scale, and binomial models
#' with a logit link use the log-odds scale. See [predict.mfpi()] for
#' interpretation by model family.
#'
#' Confidence intervals use a fixed confidence level of 95%. They are
#' pointwise intervals, shown as shaded bands for continuous predictors
#' and error bars for binary and categorical predictors. They are conditional
#' on the fitted model and do not account for uncertainty from selecting
#' variables, transformations, or interactions.
#'
#' @section Reference lines:
#' Difference plots can show two horizontal reference lines.
#'
#' The zero line indicates no difference between the groups.
#' Add it using `show_null_line = TRUE`.
#'
#' The main-effects line shows the group difference estimated by the model
#' without interaction. That model assumes the predictor has the same effect
#' in every group, so the group difference is constant across predictor values.
#' This line is shown by default; use `show_maineffect_line = FALSE` to hide it.
#' Use `show_ci_maineffect = TRUE` to add its 95% confidence band.
#'
#' @return
#' Invisibly returns a named list with one element per plotted interaction
#' variable. Each element contains a list of plots comparing individual
#' groups with the reference group.
#'
#' For multinomial models, each variable's results contain an additional
#' list level for the non-reference response categories.
#'
#' Individual plots are \pkg{ggplot2} objects. With `plot_type = "both"`,
#' each group comparison contains a combined \pkg{patchwork} plot.
#' Returned plots can be modified or saved using [ggplot2::ggsave()].
#' Use `auto_print = FALSE` to obtain them without displaying them.
#'
#' An empty list is returned if no interactions were selected and
#' `terms = NULL`.
#'
#' @references
#' Kipruto, E., Royston, P. and Sauerbrei, W. (2026).
#' Multivariable Fractional Polynomial Interaction Models in R.
#' *Journal of Statistical Software*.
#'
#' @examples
#' data("prostate")
#'
#' # Include the grouping variable in the formula and name it in group_var.
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   family = gaussian(),
#'   group_var = "svi",
#'   interaction_vars = c("cavol", "age"),
#'   flex = "flex1",
#'   p_interact = 0.9,
#'   include_group_var = TRUE,
#'   center = FALSE,
#'   verbose = FALSE
#' )
#'
#' # Plot an interaction even if it was not selected.
#' plot(fit, terms = "cavol", plot_type = "difference")
#'
#' # Plot group-specific functions with confidence bands.
#' plot(
#'   fit, terms = "cavol", plot_type = "fitted",
#'   show_ci_fitted = TRUE
#' )
#'
#' # Return plots without displaying them.
#' plots <- plot(
#'   fit, terms = "cavol", plot_type = "difference",
#'   auto_print = FALSE
#' )
#' print(plots$cavol[[1]])
#'
#' \donttest{
#' # Side-by-side plots require patchwork.
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   plot(fit, terms = "cavol", plot_type = "both")
#' }
#' }
#'
#' @seealso [mfpi()], [predict.mfpi()], [summary.mfpi()]
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#'
#' @method plot mfpi
#' @export
plot.mfpi <- function(x,
                      terms          = NULL,
                      plot_type      = c("both", "fitted", "difference"),
                      auto_print     = TRUE,
                      show_title     = TRUE,
                      show_rug       = FALSE,
                      show_ci_fitted = FALSE,
                      show_ci_diff   = TRUE,
                      show_null_line       = FALSE,
                      show_maineffect_line = TRUE,
                      show_ci_maineffect   = FALSE,
                      legend_position = c("inside", "right", "left", "bottom", "top", "none"),
                      legend_inside  = c(0.02, 0.98),
                      legend_justification = c(0, 1),
                      colour_ref     = "#2166AC",
                      colour_group   = "#D6604D",
                      line_types     = NULL,
                      colour_diff    = "#1B7837",
                      colour_null       = "grey50",
                      colour_maineffect = "grey20",
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

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting. Please install it.",
         call. = FALSE)
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
  check_logical_scalar(show_null_line, "show_null_line")
  check_logical_scalar(show_maineffect_line, "show_maineffect_line")
  check_logical_scalar(show_ci_maineffect, "show_ci_maineffect")

  check_colour_scalar <- function(value, name) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("! `", name, "` must be a single non-empty character string.",
           call. = FALSE)
    }
    invisible(TRUE)
  }
  check_colour_scalar(colour_null, "colour_null")
  check_colour_scalar(colour_maineffect, "colour_maineffect")

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
        "Nothing to plot: no variables in `interaction_vars` were selected as\n",
        "having a significant interaction with the group variable.\n",
        "To inspect a variable's fitted curve anyway, pass its name via\n",
        "the `terms` argument, e.g. terms = \"age\"."
      )
      return(invisible(list()))
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
  # Subtitle metadata
  # ---------------------------------------------------------------------------
  subtitle_lookup <- mfpi_plot_subtitle_lookup(model)

  # Group metadata translates internal prediction codes such as "0" and "1"
  # into the original group labels shown to users in legends and axis labels.
  group_meta <- mfpi_plot_group_metadata(model)
  group_label <- group_meta$group_var

  # For multinomial models the plotting loop iterates over the non-reference
  # logits and sets `active_logit` so subtitles can name the logit being shown.
  active_logit <- NULL
  active_reference_class <- NULL

  # ---------------------------------------------------------------------------
  # Helper: subtitle for one variable/contrast
  # ---------------------------------------------------------------------------
  make_subtitle <- function(var, grp_label_value, ref_label_value) {
    info <- subtitle_lookup[[var]]

    logit_suffix <- if (!is.null(active_logit)) {
      sprintf("; logit: %s vs %s", active_logit, active_reference_class)
    } else {
      ""
    }

    if (is.null(info)) {
      return(paste0(var, logit_suffix))
    }

    type_str <- switch(
      info$type,
      linear = "Linear",
      fp1    = "FP1",
      fp2    = "FP2",
      info$type
    )

    contrast_text <- sprintf(
      "%s: %s vs %s",
      group_label,
      grp_label_value,
      ref_label_value
    )

    sprintf("%s interaction; %s%s", type_str, contrast_text, logit_suffix)
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
    colour_values <- stats::setNames(c(colour_ref, colour_group), group_values)
    linetype_values <- resolve_line_types(group_values)

    interaction_spec <- mfpi_get_interaction_spec(model, var)
    if (isTRUE(interaction_spec$discrete)) {
      level_order <- interaction_spec$level_labels
      if (is.null(level_order)) level_order <- unique(as.character(line_df$x))
      line_df$x <- factor(as.character(line_df$x), levels = level_order)

      if (isTRUE(show_ci_fitted)) {
        interval_df <- rbind(
          transform(ref_df, group = ref_label_value),
          transform(grp_df, group = grp_label_value)
        )
        interval_df$group <- factor(
          interval_df$group,
          levels = c(ref_label_value, grp_label_value)
        )
        interval_df$x <- factor(as.character(interval_df$x), levels = level_order)
      }

      p <- ggplot2::ggplot(
        line_df,
        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)
      )
      if (isTRUE(show_ci_fitted)) {
        p <- p + ggplot2::geom_errorbar(
          data = interval_df,
          ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
          width = 0.12,
          position = ggplot2::position_dodge(width = 0.35)
        )
      }
      p <- p +
        ggplot2::geom_point(
          size = 2.2,
          position = ggplot2::position_dodge(width = 0.35)
        ) +
        ggplot2::scale_colour_manual(values = colour_values, name = group_label) +
        ggplot2::xlab(var) +
        ggplot2::ylab(mfpi_plot_fitted_ylabel(model)) +
        ggplot2::labs(
          title = if (show_title) "Group-specific fitted functions" else NULL,
          subtitle = if (show_title) {
            make_subtitle(var, grp_label_value, ref_label_value)
          } else NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())

      if (legend_position == "inside") {
        p <- p + ggplot2::theme(
          legend.position = legend_inside,
          legend.justification = legend_justification
        )
      } else {
        p <- p + ggplot2::theme(legend.position = legend_position)
      }
      return(p)
    }

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
          fill = colour_group,
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
      # These are group-specific partial linear predictors, not complete
      # fitted responses. The helper reports the family-appropriate link scale.
      ggplot2::ylab(mfpi_plot_fitted_ylabel(model)) +
      ggplot2::labs(
        title = if (show_title) "Group-specific fitted functions" else NULL,
        subtitle = if (show_title) {
          make_subtitle(var, grp_label_value, ref_label_value)
        } else {
          NULL
        }
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())

    if (legend_position == "inside") {
      p <- p +
        ggplot2::theme(
          legend.position = legend_inside,
          legend.justification = legend_justification,
          legend.background = ggplot2::element_rect(
            fill = ggplot2::alpha("white", 0.85),
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
                             ref_label_value, grp_label_value,
                             main_effect_value = NA_real_,
                             main_effect_se    = NA_real_) {
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

    interaction_spec <- mfpi_get_interaction_spec(model, var)
    if (isTRUE(interaction_spec$discrete)) {
      level_order <- interaction_spec$level_labels
      if (is.null(level_order)) level_order <- unique(as.character(plot_df$x))
      plot_df$x <- factor(as.character(plot_df$x), levels = level_order)

      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$x, y = .data$diff)
      )
      if (isTRUE(show_maineffect_line) && is.finite(main_effect_value) &&
          isTRUE(show_ci_maineffect) && is.finite(main_effect_se)) {
        crit <- stats::qnorm(0.975)
        # Train the discrete x scale before a rectangle whose endpoints are
        # infinite; otherwise ggplot2 can infer a continuous x scale first.
        p <- p + ggplot2::geom_blank() + ggplot2::annotate(
          "rect", xmin = -Inf, xmax = Inf,
          ymin = main_effect_value - crit * main_effect_se,
          ymax = main_effect_value + crit * main_effect_se,
          fill = colour_maineffect, alpha = ribbon_alpha
        )
      }
      if (isTRUE(show_ci_diff)) {
        p <- p + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.12,
          colour = colour_diff
        )
      }
      p <- p + ggplot2::geom_point(colour = colour_diff, size = 2.2)
      if (isTRUE(show_null_line)) {
        p <- p + ggplot2::geom_hline(
          yintercept = 0, linetype = "dashed", colour = colour_null
        )
      }
      if (isTRUE(show_maineffect_line) && is.finite(main_effect_value)) {
        p <- p + ggplot2::geom_hline(
          yintercept = main_effect_value,
          linetype = "dashed",
          colour = colour_maineffect
        )
      }
      return(
        p +
          ggplot2::xlab(var) +
          ggplot2::ylab(mfpi_plot_difference_ylabel(
            model = model,
            group = grp_label_value,
            reference = ref_label_value
          )) +
          ggplot2::labs(
            title = if (show_title) "Difference in fitted functions" else NULL,
            subtitle = if (show_title) {
              make_subtitle(var, grp_label_value, ref_label_value)
            } else NULL
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid = ggplot2::element_blank())
      )
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
      ggplot2::geom_line(
        ggplot2::aes(y = .data$diff),
        colour = colour_diff,
        linewidth = linewidth
      )

    # Reference lines. Drawn on top of the ribbon but below the fitted curve
    # would be preferable; drawn after the curve is acceptable and matches the
    # legacy behaviour, which placed the zero line just before the curve. The
    # `after_stat` approach is not needed here because both references are
    # constants known at plotting time.
    if (isTRUE(show_null_line)) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          colour = colour_null
        )
    }

    if (isTRUE(show_maineffect_line) &&
        is.finite(main_effect_value)) {
      # Confidence band for the main-effects contrast alpha^{(M)}.
      # Drawn before the line so the ribbon sits behind it.
      if (isTRUE(show_ci_maineffect) && is.finite(main_effect_se)) {
        me_lower <- main_effect_value - stats::qnorm(0.975) * main_effect_se
        me_upper <- main_effect_value + stats::qnorm(0.975) * main_effect_se
        x_range <- range(plot_df$x, na.rm = TRUE)
        me_band_df <- data.frame(
          x    = x_range,
          ymin = rep(me_lower, 2L),
          ymax = rep(me_upper, 2L)
        )
        p <- p +
          ggplot2::geom_ribbon(
            data = me_band_df,
            ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax),
            fill = colour_maineffect,
            alpha = ribbon_alpha,
            inherit.aes = FALSE
          )
      }
      p <- p +
        ggplot2::geom_hline(
          yintercept = main_effect_value,
          linetype = "dashed",
          colour = colour_maineffect
        )
    }

    p <- p +
      ggplot2::xlab(var) +
      # The difference includes the group main-effect contrast as well as the
      # difference between group-specific covariate functions. Label the
      # resulting family-specific contrast directly.
      ggplot2::ylab(
        mfpi_plot_difference_ylabel(
          model = model,
          group = grp_label_value,
          reference = ref_label_value
        )
      ) +
      ggplot2::labs(
        title = if (show_title) "Difference in fitted functions" else NULL,
        subtitle = if (show_title) {
          make_subtitle(var, grp_label_value, ref_label_value)
        } else {
          NULL
        }
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())

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

    # Multinomial predictions carry one fitted curve / contrast per non-reference
    # logit. Iterate over the logits, reusing the single-response plotting
    # machinery on each per-logit slice; other families run once with no slicing.
    is_multinomial_pred <- isTRUE(pred$metadata$multinomial)
    logit_codes <- if (is_multinomial_pred) pred$metadata$classes else NA_character_
    active_reference_class <- if (is_multinomial_pred) {
      pred$metadata$reference_class
    } else {
      NULL
    }
    term_plots <- list()

    for (logit_code in logit_codes) {
    active_logit <- if (is_multinomial_pred) logit_code else NULL

    functions_df <- pred$functions
    differences_df <- pred$differences
    if (is_multinomial_pred) {
      if (!is.null(functions_df) && "class" %in% names(functions_df)) {
        functions_df <- functions_df[functions_df$class == logit_code, , drop = FALSE]
      }
      if (!is.null(differences_df) && "class" %in% names(differences_df)) {
        differences_df <- differences_df[differences_df$class == logit_code, , drop = FALSE]
      }
    }

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

    # -------------------------------------------------------------------------
    # Look up the main-effects group contrast alpha^{(M)} for this variable,
    # one value per comparison group. If the main-effects model is not
    # preserved on the fitted object, `me_alpha` is NULL and the reference
    # line is silently skipped.
    # -------------------------------------------------------------------------
    me_alpha <- NULL
    if (isTRUE(show_maineffect_line) &&
        plot_type %in% c("difference", "both") &&
        length(comparison_codes) > 0L) {
      me_alpha <- mfpi_plot_maineffect_alpha(
        model            = model,
        var              = var,
        group_var_name   = group_label,
        comparison_codes = comparison_codes,
        logit            = active_logit
      )
    }

    if (length(comparison_codes) == 0L) {
      warning(
        paste0("No non-reference group contrasts are available for term `", var, "`."),
        call. = FALSE
      )
      if (is_multinomial_pred) {
        term_plots[[paste0("logit_", logit_code)]] <- list()
      } else {
        term_plots <- list()
      }
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

      me_alpha_j <- if (!is.null(me_alpha) &&
                        grp_code %in% names(me_alpha$estimate)) {
        me_alpha$estimate[[grp_code]]
      } else {
        NA_real_
      }

      me_se_j <- if (!is.null(me_alpha) &&
                     grp_code %in% names(me_alpha$se)) {
        me_alpha$se[[grp_code]]
      } else {
        NA_real_
      }

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
          grp_label_value = grp_label_value,
          main_effect_value = me_alpha_j,
          main_effect_se    = me_se_j
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
            grp_label_value = grp_label_value,
            main_effect_value = me_alpha_j,
            main_effect_se    = me_se_j
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

    if (is_multinomial_pred) {
      term_plots[[paste0("logit_", logit_code)]] <- grp_plots
    } else {
      term_plots <- grp_plots
    }
    } # end per-logit loop

    plots[[var]] <- term_plots
  }

  invisible(plots)
}


# -----------------------------------------------------------------------------
# plot.mfpi() helper functions ------------------------------------------------
# -----------------------------------------------------------------------------

#' Determine the Link Used for MFPI Plot Labels
#'
#' Extracts the fitted GLM link for use in plot-axis labels. Cox and `survreg`
#' models use dedicated values because they do not store a GLM
#' family object, and negative-binomial models use their supported log link.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return Character scalar naming the link or survival-model scale.
#'
#' @keywords internal
#' @noRd
mfpi_plot_link_name <- function(model) {
  family_string <- tolower(mfpi_family_string(model))

  if (identical(family_string, "cox")) {
    return("cox")
  }
  if (identical(family_string, "survreg")) return("survreg")
  if (identical(family_string, "negbin")) {
    return("log")
  }

  family_obj <- model[["family", exact = TRUE]]
  if (is.list(family_obj)) {
    link <- family_obj[["link", exact = TRUE]]
    if (!is.null(link) && length(link) > 0L && !is.na(link[1L]) &&
        nzchar(as.character(link[1L]))) {
      return(tolower(as.character(link[1L])))
    }
  }

  switch(
    family_string,
    gaussian = "identity",
    binomial = "logit",
    poisson  = "log",
    gamma = "inverse",
    inverse.gaussian = "1/mu^2",
    "unknown"
  )
}


#' Build the Fitted-Curve Y-Axis Label
#'
#' Returns a family- and link-specific label for the group-specific partial
#' linear predictors displayed by \code{plot.mfpi()}.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return Character scalar.
#'
#' @keywords internal
#' @noRd
mfpi_plot_fitted_ylabel <- function(model) {
  family_string <- tolower(mfpi_family_string(model))
  link <- mfpi_plot_link_name(model)

  if (identical(family_string, "cox")) {
    return("Partial linear predictor (log relative hazard)")
  }
  if (identical(family_string, "survreg")) {
    return("Partial linear predictor (survreg location scale)")
  }

  scale_label <- switch(
    link,
    identity = "outcome scale",
    logit    = "log odds",
    log      = "log mean",
    probit   = "probit scale",
    cloglog  = "complementary log-log scale",
    cauchit  = "cauchit scale",
    inverse  = "inverse-mean scale",
    `1/mu^2` = "inverse-squared-mean scale",
    sqrt     = "square-root mean scale",
    if (identical(link, "unknown")) "link scale" else paste0(link, " scale")
  )

  paste0("Partial linear predictor (", scale_label, ")")
}


#' Build the Difference-Plot Y-Axis Label
#'
#' Returns the family-specific interpretation of the comparison-group minus
#' reference-group contrast displayed by \code{plot.mfpi()}.
#'
#' @param model Object of class \code{"mfpi"}.
#' @param group Character scalar giving the comparison-group label.
#' @param reference Character scalar giving the reference-group label.
#'
#' @return Character scalar.
#'
#' @keywords internal
#' @noRd
mfpi_plot_difference_ylabel <- function(model, group, reference) {
  family_string <- tolower(mfpi_family_string(model))
  link <- mfpi_plot_link_name(model)

  measure <- if (identical(family_string, "cox")) {
    "Log hazard ratio"
  } else if (identical(family_string, "survreg")) {
    "Location-scale difference"
  } else if (identical(family_string, "gaussian") &&
             identical(link, "identity")) {
    "Mean difference"
  } else if (identical(family_string, "binomial") &&
             identical(link, "logit")) {
    "Log odds ratio"
  } else if (family_string %in% c("poisson", "negbin") &&
             identical(link, "log")) {
    "Log mean ratio"
  } else {
    "Link-scale difference"
  }

  sprintf("%s: %s vs %s", measure, group, reference)
}


#' Return the Compact Plot Scale Label for an MFPI Model
#'
#' Backward-compatible internal helper returning a compact scale label. New
#' plot axes use \code{mfpi_plot_fitted_ylabel()} and
#' \code{mfpi_plot_difference_ylabel()} so that the plotted quantity, rather
#' than only its scale, is identified.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return Character scalar.
#'
#' @keywords internal
#' @noRd
mfpi_plot_scale_label <- function(model) {
  family_string <- tolower(mfpi_family_string(model))
  link <- mfpi_plot_link_name(model)

  if (identical(family_string, "cox")) {
    return("log relative hazard")
  }
  if (identical(family_string, "survreg")) return("survreg location scale")

  switch(
    link,
    identity = "outcome scale",
    logit    = "log odds",
    log      = "log mean",
    probit   = "probit scale",
    cloglog  = "complementary log-log scale",
    cauchit  = "cauchit scale",
    inverse  = "inverse-mean scale",
    `1/mu^2` = "inverse-squared-mean scale",
    sqrt     = "square-root mean scale",
    if (identical(link, "unknown")) "linear predictor" else paste0(link, " scale")
  )
}


#' Build Subtitle Metadata for MFPI Plots
#'
#' Creates a named lookup table containing the interaction form for each
#' variable stored in an \code{"mfpi"} object's interaction metrics.
#'
#' This helper is used only by \code{plot.mfpi()}. It reads
#' \code{model$all_model_metrics} and extracts the fitted interaction form used
#' in plot subtitles. Model-selection statistics remain available through the
#' model's print and summary methods but are not repeated on the plots.
#'
#' @param model Object of class \code{"mfpi"}.
#'
#' @return A named list indexed by continuous-variable name. Each element is a
#'   list with component \code{type}. If no model metrics are available, an
#'   empty list is returned.
#'
#' @keywords internal
#' @noRd
mfpi_plot_subtitle_lookup <- function(model) {
  out <- list()
  all_metrics <- model$all_model_metrics

  if (is.null(all_metrics) || nrow(all_metrics) == 0L) {
    return(out)
  }

  for (i in seq_len(nrow(all_metrics))) {
    v <- all_metrics$variable[i]
    form <- if ("type" %in% names(all_metrics)) all_metrics$type[i] else NA_character_

    if (is.na(v) || !nzchar(v)) next
    out[[v]] <- list(type = form)
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

#' Extract the Main-Effects Group Contrast for a Variable
#'
#' Attempts to locate a preserved main-effects (no-interaction) model for the
#' given continuous variable within an \code{mfpi} object and returns a named
#' numeric vector giving the estimated group main-effect contrast
#' \eqn{\hat\alpha^{(M)}} for each internal comparison group code relative to
#' the reference group.
#'
#' The function is deliberately tolerant. It looks in several plausible
#' locations, and returns \code{NULL} without warning when no main-effects
#' model can be located. This allows \code{plot.mfpi()} to silently omit the
#' \eqn{\hat\alpha^{(M)}} reference line for older fitted objects that did not
#' preserve the main-effects model, while continuing to draw it when the model
#' is available.
#'
#' Searched locations, in order:
#' \enumerate{
#'   \item \code{model$best_main_model[[var]]$fit}
#'   \item \code{model$best_main_model[[var]]}
#'   \item \code{model$var_winners[[var]]$fit$test_results$main_model$fit}
#'   \item \code{model$var_winners[[var]]$fit$test_results$main_model}
#' }
#'
#' Group-dummy coefficients are identified by name, matching the naming
#' convention used elsewhere in \pkg{mfp2}. For a grouping variable
#' \code{group_var} coded internally as \code{0, 1, ..., K - 1}, the main
#' effect of internal code \code{k} is expected to appear as a coefficient
#' named \code{"<group_var><k>"} (for example, \code{"rx1"} for a binary
#' contrast). Non-matching coefficients are ignored.
#'
#' @param model An \code{mfpi} object.
#' @param var Character scalar naming the continuous variable of interest.
#' @param group_var_name Character scalar naming the grouping variable.
#' @param comparison_codes Character vector of internal comparison-group codes
#'   for which to extract the main-effect estimate.
#' @param logit Non-reference outcome class for a multinomial model; `NULL`
#'   for scalar-response families.
#'
#' @return A list with two named numeric vectors, each with names equal to
#'   \code{comparison_codes}:
#'   \describe{
#'     \item{\code{estimate}}{The estimated group main-effect contrast for
#'       each comparison group.}
#'     \item{\code{se}}{The corresponding standard error, obtained from the
#'       main-effects model's covariance matrix or its stored multinomial
#'       information matrix.
#'       Elements are \code{NA} when the standard error cannot be extracted.}
#'   }
#'   Returns \code{NULL} if no main-effects model can be located.
#'
#' @keywords internal
#' @noRd
mfpi_plot_maineffect_alpha <- function(model, var, group_var_name,
                                       comparison_codes, logit = NULL) {
  # ---------------------------------------------------------------------------
  # Locate the fitted main-effects model
  # ---------------------------------------------------------------------------
  candidates <- list(
    tryCatch(model$best_main_model[[var]]$fit, error = function(e) NULL),
    tryCatch(model$best_main_model[[var]],      error = function(e) NULL),
    tryCatch(model$var_winners[[var]]$fit$test_results$main_model$fit,
             error = function(e) NULL),
    tryCatch(model$var_winners[[var]]$fit$test_results$main_model,
             error = function(e) NULL)
  )

  main_fit <- NULL
  for (cand in candidates) {
    if (!is.null(cand)) {
      main_fit <- cand
      break
    }
  }
  if (is.null(main_fit)) return(NULL)

  # ---------------------------------------------------------------------------
  # Extract coefficients defensively
  # ---------------------------------------------------------------------------
  # The retained fast multinomial fit stores its named per-logit coefficients
  # here; coef.nnet() instead returns a flat vector of neural-network weights.
  cf <- if (!is.null(logit) &&
            is.matrix(main_fit$mfp2_coefficient_matrix)) {
    main_fit$mfp2_coefficient_matrix
  } else {
    tryCatch(stats::coef(main_fit), error = function(e) NULL)
  }
  if (is.null(cf) && is.list(main_fit)) {
    # Some internal fit objects wrap the model; try likely component names.
    for (comp in c("coefficients", "coef", "beta")) {
      if (!is.null(main_fit[[comp]])) {
        cf <- main_fit[[comp]]
        break
      }
    }
  }
  if (is.matrix(cf)) {
    if (is.null(logit) || !logit %in% rownames(cf)) return(NULL)
    cf <- cf[logit, , drop = TRUE]
  }
  if (is.null(cf) || length(cf) == 0L || is.null(names(cf))) return(NULL)

  # ---------------------------------------------------------------------------
  # Extract variance-covariance matrix for standard errors.
  # Some internal fit objects stored by mfpi() have the structure of a glm/lm
  # but are stored as plain lists without a class attribute. Try vcov() first;
  # if that fails, check whether the object looks like an unclassed glm/lm
  # (has qr, residuals, coefficients) and temporarily assign the class so that
  # stats::vcov.glm / vcov.lm can succeed.
  # ---------------------------------------------------------------------------
  vc <- tryCatch(stats::vcov(main_fit), error = function(e) NULL)

  # MFPI retains fast multinomial main-effects fits as nnet objects. They
  # contain the design and fitted probabilities but no vcov() method. Compute
  # the information for their freely estimated baseline-category logits using
  # the same routine as the final multinomial fit (including case weights).
  if (is.null(vc) && !is.null(logit) &&
      !is.null(main_fit$mfp2_design) &&
      !is.null(main_fit$mfp2_family$prepared) &&
      !is.null(main_fit$fitted.values)) {
    prepared <- main_fit$mfp2_family$prepared
    vc <- tryCatch({
      information <- mfp2_multinomial_information(
        x = main_fit$mfp2_design,
        probabilities = main_fit$fitted.values,
        weights = prepared$effective_weights,
        outcomes = prepared$nonreference,
        class_levels = prepared$levels
      )
      solve(information)
    }, error = function(e) NULL)
  }

  if (is.null(vc) && is.list(main_fit) && is.null(class(main_fit)) ||
      (is.null(vc) && is.list(main_fit) && identical(class(main_fit), "list"))) {
    glm_components <- c("coefficients", "residuals", "qr", "rank",
                        "family", "deviance", "df.residual")
    if (all(glm_components %in% names(main_fit))) {
      tmp <- main_fit
      class(tmp) <- c("glm", "lm")
      vc <- tryCatch(stats::vcov(tmp), error = function(e) NULL)
    } else {
      lm_components <- c("coefficients", "residuals", "qr", "rank",
                         "df.residual")
      if (all(lm_components %in% names(main_fit))) {
        tmp <- main_fit
        class(tmp) <- "lm"
        vc <- tryCatch(stats::vcov(tmp), error = function(e) NULL)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Match group-dummy coefficient names for each requested comparison code
  # ---------------------------------------------------------------------------
  comparison_codes <- as.character(comparison_codes)
  out_est <- stats::setNames(rep(NA_real_, length(comparison_codes)),
                             comparison_codes)
  out_se  <- stats::setNames(rep(NA_real_, length(comparison_codes)),
                             comparison_codes)

  for (code in comparison_codes) {
    nm <- paste0(group_var_name, code)
    if (nm %in% names(cf)) {
      value <- unname(cf[[nm]])
      if (is.numeric(value) && length(value) == 1L && !is.na(value)) {
        out_est[[code]] <- value
      }
      # Extract SE from the diagonal of vcov if available.
      vc_name <- if (is.null(logit)) nm else paste0(logit, ":", nm)
      if (!is.null(vc) && vc_name %in% rownames(vc) &&
          vc_name %in% colnames(vc)) {
        var_value <- vc[vc_name, vc_name]
        if (is.numeric(var_value) && length(var_value) == 1L &&
            is.finite(var_value) && var_value >= 0) {
          out_se[[code]] <- sqrt(var_value)
        }
      }
    }
  }

  if (all(is.na(out_est))) return(NULL)
  list(estimate = out_est, se = out_se)
}
