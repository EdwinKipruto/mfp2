#' Plot Fitted Functions from an MFPI Model
#'
#' Produces group-specific fitted-function plots, between-group difference
#' plots, or both for continuous variables evaluated by [mfpi()].
#'
#' The plotted functions are shown on the model's linear-predictor scale.
#' They represent the fitted contribution of the selected continuous variable,
#' rather than a complete predicted response.
#'
#' @section Plot types:
#' Three plot types are available through `plot_type`.
#'
#' \describe{
#'   \item{\code{"fitted"}}{
#'     Plots the fitted covariate function separately for the reference group
#'     and one comparison group on shared axes.
#'
#'     Confidence bands for the group-specific functions are controlled by
#'     `show_ci_fitted`.
#'   }
#'
#'   \item{\code{"difference"}}{
#'     Plots the pointwise difference between the comparison-group and
#'     reference-group fitted functions:
#'
#'     \deqn{\hat f_j(x) - \hat f_r(x)}
#'
#'     Two optional horizontal reference lines can be added, each answering a
#'     different question about the plot:
#'
#'     \itemize{
#'       \item A dashed line at zero, controlled by `show_null_line`, marks
#'         the value of the fitted-function difference that corresponds to a
#'         null group contrast at that covariate value. On the linear-predictor
#'         scale this means equal partial linear predictors; on the response
#'         scale it corresponds to a hazard ratio, rate ratio, or odds ratio
#'         of one for the applicable family, and to a zero mean difference for
#'         Gaussian models. This line is a treatment-effect reference, not an
#'         interaction reference: a difference curve can cross zero without
#'         implying interaction, and a flat nonzero difference implies no
#'         interaction even though it never crosses zero.
#'
#'       \item A dashed line at the main-effects contrast \eqn{\hat\alpha^{(M)}},
#'         controlled by `show_maineffect_line`, marks the constant fitted
#'         difference implied by the no-interaction (main-effects) model. This
#'         is the natural interaction reference: departure of the fitted
#'         difference curve from this horizontal line is the visual analogue
#'         of the Step 2 interaction comparison. When the main-effects and
#'         interaction models share the same FP basis (as they do under
#'         `flex1` and `flex2`), this comparison is exact. Under `flex3` and
#'         `flex4` the two fitted models may use different FP bases, in which
#'         case the line still summarises what the no-interaction model
#'         concluded but the comparison need not be strictly nested.
#'     }
#'
#'     The two lines coincide only when there is no group main effect, which
#'     is unusual in practice. Drawing `show_maineffect_line = TRUE` requires
#'     that the main-effects (no-interaction) model has been preserved on the
#'     `mfpi` object; if it is not available, the line is silently omitted.
#'
#'     The confidence band for the difference curve is controlled by
#'     `show_ci_diff`.
#'   }
#'
#'   \item{\code{"both"}}{
#'     Places the group-specific fitted-function plot and the corresponding
#'     difference plot side by side. This option requires the \pkg{patchwork}
#'     package.
#'   }
#' }
#'
#' @section Model scale and interpretation:
#' The plotted functions are expressed on the model's linear-predictor scale:
#'
#' \itemize{
#'   \item Gaussian models: outcome scale;
#'   \item binomial models: log-odds scale;
#'   \item Poisson models: log-mean scale;
#'   \item Cox models: log-hazard scale.
#' }
#'
#' A difference plot therefore displays a difference between fitted covariate
#' functions on the relevant model scale. It should not be interpreted directly
#' as a probability, odds ratio, rate ratio, hazard ratio, survival probability,
#' or complete predicted response.
#'
#' Group labels shown in the plots correspond to the original values or factor
#' levels supplied through `group_var`.
#'
#' When confidence intervals are displayed, they are pointwise intervals for
#' the fitted function or fitted-function difference. They do not account for
#' uncertainty introduced by variable selection, functional-form selection, or
#' interaction selection.
#'
#' @section Term selection:
#' If `terms = NULL`, plots are produced for interactions retained by the
#' selected MFPI criterion.
#'
#' If `terms` is supplied, the requested continuous variables are plotted
#' regardless of whether their interactions were retained, provided that the
#' corresponding fitted models are available in the `mfpi` object.
#'
#' @section Confidence bands:
#' Confidence bands for the plot elements are controlled separately:
#'
#' \itemize{
#'   \item `show_ci_fitted` controls confidence bands around the
#'     group-specific fitted functions;
#'   \item `show_ci_diff` controls the confidence band around the
#'     fitted-function difference;
#'   \item `show_ci_maineffect` controls the confidence band around the
#'     main-effects reference line \eqn{\hat\alpha^{(M)}} on difference plots.
#' }
#'
#' The default suppresses confidence bands for the group-specific functions
#' and the main-effects line, and displays them for the difference curve.
#'
#' @param x An object of class `"mfpi"`.
#'
#' @param terms Optional character vector naming the continuous variables to
#'   plot. If `NULL`, only retained interaction terms are plotted; when no
#'   interactions were retained, no plots are produced. If supplied, the
#'   requested terms are plotted when their fitted models are available in the
#'   `all_interaction_models` component of the `mfpi` object.
#'
#' @param plot_type Character scalar specifying the plot to produce. One of
#'   `"fitted"`, `"difference"`, or `"both"`.
#'
#' @param auto_print Logical scalar. If `TRUE`, each plot is printed when it is
#'   created. If `FALSE`, plots are returned invisibly for manual printing,
#'   modification, or saving.
#'
#' @param show_title Logical scalar. If `TRUE`, titles and subtitles are added.
#'   Subtitles may include the interaction type, group contrast, and available
#'   selection result.
#'
#' @param show_rug Logical scalar. If `TRUE`, rug marks are added along the
#'   horizontal axis to show the covariate evaluation values.
#'
#' @param show_ci_fitted Logical scalar. If `TRUE`, confidence bands are shown
#'   around the group-specific fitted functions.
#'
#' @param show_ci_diff Logical scalar. If `TRUE`, a confidence band is shown
#'   around the fitted-function difference.
#'
#' @param show_null_line Logical scalar. If `TRUE` (default), a dashed
#'   horizontal reference line at zero is drawn on difference plots. This line
#'   marks the null-group-contrast value on the linear-predictor scale. It is
#'   informative for reading the treatment-effect scale at a given covariate
#'   value, but it is not an interaction reference (see the plot-type
#'   documentation).
#'
#' @param show_maineffect_line Logical scalar. If `TRUE`, a dashed horizontal
#'   reference line is drawn on difference plots at the constant contrast
#'   \eqn{\hat\alpha^{(M)}} implied by the no-interaction (main-effects) model.
#'   Departure of the fitted difference curve from this line is the visual
#'   analogue of the interaction comparison. Requires that the main-effects
#'   model has been preserved on the fitted `mfpi` object; if it cannot be
#'   located, the line is silently omitted. Defaults to `TRUE`.
#'
#' @param show_ci_maineffect Logical scalar. If `TRUE`, a pointwise 95\%
#'   confidence band is drawn around the main-effects reference line
#'   \eqn{\hat\alpha^{(M)}} on difference plots, computed as
#'   \eqn{\hat\alpha^{(M)} \pm z_{0.975} \cdot \mathrm{SE}(\hat\alpha^{(M)})}.
#'   The standard error is extracted from the variance-covariance matrix of
#'   the main-effects model. If the standard error cannot be obtained, the
#'   band is silently omitted. The band is only drawn when
#'   `show_maineffect_line = TRUE` and the main-effects model is available.
#'   Defaults to `FALSE`.
#'
#' @param colour_null Character scalar specifying the colour of the null
#'   reference line at zero on difference plots.
#'
#' @param colour_maineffect Character scalar specifying the colour of the
#'   main-effects reference line at \eqn{\hat\alpha^{(M)}} on difference plots.
#'
#' @param legend_position Character scalar specifying the legend position for
#'   fitted-function plots. One of `"inside"`, `"right"`, `"left"`,
#'   `"bottom"`, `"top"`, or `"none"`.
#'
#' @param legend_inside Numeric vector of length two giving the legend position
#'   inside the plotting panel when `legend_position = "inside"`. Values are
#'   interpreted as normalized parent coordinates.
#'
#' @param legend_justification Numeric vector of length two giving the legend
#'   justification when `legend_position = "inside"`.
#'
#' @param colour_ref Character scalar specifying the colour of the
#'   reference-group fitted function and confidence band.
#'
#' @param colour_group Character scalar specifying the colour of the
#'   comparison-group fitted function and confidence band.
#'
#' @param line_types Optional character vector specifying the line types of the
#'   two group-specific functions. If unnamed, it must have length two and is
#'   interpreted as the reference group followed by the comparison group. If
#'   named, its names must match the original group labels shown in the legend.
#'
#' @param colour_diff Character scalar specifying the colour of the difference
#'   curve and its confidence band.
#'
#' @param linewidth Positive numeric scalar specifying the width of the fitted
#'   and difference curves.
#'
#' @param ribbon_alpha Numeric scalar in \eqn{[0,1]} specifying the transparency
#'   of confidence bands.
#'
#' @param rug_alpha Numeric scalar in \eqn{[0,1]} specifying the transparency
#'   of rug marks.
#'
#' @param ... Currently unused. Supplied arguments generate a warning.
#'
#' @return
#' Invisibly returns a nested named list of plots.
#'
#' The outer list is named by continuous variable. Within each variable, the
#' inner list contains one element for each displayed group contrast.
#'
#' Each element is:
#'
#' \itemize{
#'   \item a `ggplot` object when `plot_type = "fitted"` or
#'     `plot_type = "difference"`;
#'   \item a `patchwork` object when `plot_type = "both"`.
#' }
#'
#' Returned plots can be printed, modified with \pkg{ggplot2}, combined with
#' \pkg{patchwork}, or saved with [ggplot2::ggsave()].
#'
#' @seealso [mfpi()], [predict.mfpi()], [summary.mfpi()]
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#'
#' @examples
#' data("prostate")
#'
#' # Investigate whether the fitted effects of cavol and age differ by svi.
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) +
#'     fp(bph) + fp(cp),
#'   data = prostate,
#'   group_var = "svi",
#'   cont_vars = c("cavol", "age"),
#'   flex = "flex1",
#'   include_group_var = TRUE,
#'   center = FALSE,
#'   verbose = FALSE
#' )
#'
#' # Plot retained group-specific fitted functions.
#' plot(fit, plot_type = "fitted")
#'
#' # Plot retained fitted-function differences.
#' plot(fit, plot_type = "difference")
#'
#' # Return plots without printing them.
#' plots <- plot(
#'   fit,
#'   terms = "cavol",
#'   plot_type = "difference",
#'   auto_print = FALSE
#' )
#'
#' # Print the first returned contrast plot for cavol.
#' print(plots[["cavol"]][[1]])
#'
#' \donttest{
#' # Display fitted and difference plots side by side.
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   plot(fit, terms = "cavol", plot_type = "both")
#' }
#' }
#'
#' # Add confidence bands to the group-specific fitted functions.
#' plot(
#'   fit,
#'   terms = "cavol",
#'   plot_type = "fitted",
#'   show_ci_fitted = TRUE
#' )
#'
#' # Use line types named according to the original group labels.
#' plot(
#'   fit,
#'   terms = "cavol",
#'   plot_type = "fitted",
#'   line_types = c("0" = "solid", "1" = "dashed")
#' )
#'
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
    colour_values <- stats::setNames(c(colour_ref, colour_group), group_values)
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
        comparison_codes = comparison_codes
      )
    }

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
    negbin   = "log-mean",
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
#'
#' @return A list with two named numeric vectors, each with names equal to
#'   \code{comparison_codes}:
#'   \describe{
#'     \item{\code{estimate}}{The estimated group main-effect contrast for
#'       each comparison group.}
#'     \item{\code{se}}{The corresponding standard error, obtained from the
#'       diagonal of the variance-covariance matrix of the main-effects model.
#'       Elements are \code{NA} when the standard error cannot be extracted.}
#'   }
#'   Returns \code{NULL} if no main-effects model can be located.
#'
#' @keywords internal
#' @noRd
mfpi_plot_maineffect_alpha <- function(model, var, group_var_name,
                                       comparison_codes) {
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
  cf <- tryCatch(stats::coef(main_fit), error = function(e) NULL)
  if (is.null(cf) && is.list(main_fit)) {
    # Some internal fit objects wrap the model; try likely component names.
    for (comp in c("coefficients", "coef", "beta")) {
      if (!is.null(main_fit[[comp]])) {
        cf <- main_fit[[comp]]
        break
      }
    }
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
      if (!is.null(vc) && nm %in% rownames(vc) && nm %in% colnames(vc)) {
        var_value <- vc[nm, nm]
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