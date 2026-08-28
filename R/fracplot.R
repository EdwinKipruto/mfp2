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
#' @param color_line Colour used for fitted curves and for binary fitted
#'   points and confidence intervals.
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
#' Binary predictors are displayed at their two fitted levels as point estimates
#' with vertical confidence intervals. Only those two levels are labelled on the
#' x-axis; no interpolating line or confidence ribbon is drawn. For a spike-at-zero
#' term retaining only the zero indicator, the internal 0/1 coding is replaced by
#' semantic labels such as \code{x = 0} and \code{x > 0}.
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
#' @seealso [predict.mfp2()], [mfp2()]
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
                      color_line = "black",
                      color_line_spike = "black",
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

# Identify fitted terms that should be displayed as two-level effects rather
# than as continuous curves. Formula factors with two levels are binary;
# singleton numeric terms are checked on the fitted raw design values.
mfp2_plot_term_is_binary <- function(model, term) {
  factor_info <- if (!is.null(model$formula_factor_info)) {
    model$formula_factor_info[[term]]
  } else {
    NULL
  }

  if (!is.null(factor_info) && length(factor_info$levels) == 2L) {
    return(TRUE)
  }

  term_columns <- if (!is.null(model$term_to_columns) &&
                      term %in% names(model$term_to_columns)) {
    model$term_to_columns[[term]]
  } else {
    term
  }

  if (length(term_columns) != 1L ||
      is.null(model$x_original) ||
      is.null(colnames(model$x_original)) ||
      !term_columns %in% colnames(model$x_original)) {
    return(FALSE)
  }

  values <- model$x_original[, term_columns, drop = TRUE]
  values <- values[!is.na(values)]

  length(unique(values)) == 2L
}

# Formula factors are discrete effects regardless of how many fitted levels
# they contain. Keep this check separate from binary detection so multi-level
# factors are never sent through the continuous line/ribbon plotting path.
mfp2_plot_term_is_factor <- function(model, term) {
  !is.null(model$formula_factor_info) &&
    !is.null(model$formula_factor_info[[term]])
}

# Re-evaluate one discrete term on its fitted observations. This keeps the general
# predict.mfp2() API unchanged while ensuring that plot(..., terms_seq =
# "equidistant") does not display artificial intermediate binary values.
mfp2_plot_binary_prediction <- function(model, term, type, ref, alpha) {
  term_ref <- NULL
  if (identical(type, "contrasts") &&
      !is.null(ref) &&
      term %in% names(ref) &&
      !is.null(ref[[term]])) {
    term_ref <- stats::setNames(list(ref[[term]]), term)
  }

  out <- predict(
    model,
    type = type,
    terms = term,
    ref = term_ref,
    terms_seq = "data",
    alpha = alpha
  )

  out[[term]]
}

# Collapse repeated observed-value predictions to one fitted estimate per level
# and preserve fitted factor-level order for discrete axes, including factors
# with more than two levels. For a SAZ
# binary-only term, translate the internal indicator I(x == 0) into semantic
# labels on the original covariate scale: "x = 0" and "x > 0".
mfp2_plot_prepare_binary_data <- function(model, term, df, residual_df = NULL) {
  binary_df <- df[!duplicated(df$variable), , drop = FALSE]

  is_saz_binary_only <- FALSE
  if (!is.null(model$spike_dec) &&
      term %in% names(model$spike_dec) &&
      !is.na(model$spike_dec[[term]])) {
    is_saz_binary_only <- identical(
      as.integer(model$spike_dec[[term]]),
      as.integer(saz_decision_codes[["binary_only"]])
    )
  }

  if (is_saz_binary_only) {
    zero_label <- paste0(term, " = 0")
    positive_label <- paste0(term, " > 0")
    level_order <- c(zero_label, positive_label)

    saz_label <- function(x) {
      indicator <- suppressWarnings(as.numeric(as.character(x)))
      out <- rep(NA_character_, length(indicator))
      out[indicator == 1] <- zero_label
      out[indicator == 0] <- positive_label
      out
    }

    binary_df$variable <- factor(
      saz_label(binary_df$variable),
      levels = level_order
    )
    binary_df <- binary_df[order(binary_df$variable), , drop = FALSE]

    if (!is.null(residual_df)) {
      residual_df$variable <- factor(
        saz_label(residual_df$variable),
        levels = level_order
      )
    }

    return(list(
      fitted = binary_df,
      residual = residual_df,
      continuous_axis = FALSE,
      breaks = level_order
    ))
  }

  if (is.numeric(binary_df$variable)) {
    binary_df <- binary_df[order(binary_df$variable), , drop = FALSE]
    return(list(
      fitted = binary_df,
      residual = residual_df,
      continuous_axis = TRUE,
      breaks = as.numeric(binary_df$variable)
    ))
  }

  factor_info <- if (!is.null(model$formula_factor_info)) {
    model$formula_factor_info[[term]]
  } else {
    NULL
  }

  level_order <- if (!is.null(factor_info) &&
                     length(factor_info$levels) > 0L) {
    # Use the fitted factor metadata rather than observed ordering so the
    # reference and non-reference levels appear in their original order.
    as.character(factor_info$levels)
  } else {
    unique(as.character(binary_df$variable))
  }

  binary_df$variable <- factor(
    as.character(binary_df$variable),
    levels = level_order
  )
  binary_df <- binary_df[order(binary_df$variable), , drop = FALSE]

  if (!is.null(residual_df)) {
    residual_df$variable <- factor(
      as.character(residual_df$variable),
      levels = level_order
    )
  }

  list(
    fitted = binary_df,
    residual = residual_df,
    continuous_axis = FALSE,
    breaks = level_order
  )
}


# Attach fitted-model residuals to the observed-value prediction frame for
# every requested term. `predict(..., terms_seq = "data")` constructs these
# frames from `model$x_original`, so fitted-observation order is the governing
# invariant. Validate it explicitly to prevent R from recycling a shorter
# residual vector when its length happens to divide the number of rows.
mfp2_plot_attach_residuals <- function(pred_data, resid, model) {
  if (is.null(model$x_original) || is.null(nrow(model$x_original))) {
    stop(
      "Cannot construct component-plus-residual plots: the fitted object ",
      "does not contain its original predictor rows.",
      call. = FALSE
    )
  }

  if (!is.atomic(resid) || !is.null(dim(resid))) {
    stop(
      "Cannot construct component-plus-residual plots: model residuals must ",
      "be a vector.",
      call. = FALSE
    )
  }

  expected_n <- nrow(model$x_original)
  if (length(resid) != expected_n) {
    stop(
      "Cannot construct component-plus-residual plots: the model contains ",
      expected_n, " fitted observations but ", length(resid), " residuals.",
      call. = FALSE
    )
  }

  observation_ids <- rownames(model$x_original)
  residual_ids <- names(resid)

  if (!is.null(observation_ids) && !is.null(residual_ids)) {
    # When identifiers are available, use them rather than assuming the
    # residual method retained fitted-data order. Duplicate, missing, or
    # different identifiers make a one-to-one alignment impossible.
    invalid_ids <- anyNA(observation_ids) || anyNA(residual_ids) ||
      anyDuplicated(observation_ids) > 0L ||
      anyDuplicated(residual_ids) > 0L ||
      !setequal(observation_ids, residual_ids)

    if (invalid_ids) {
      stop(
        "Cannot align model residuals with the fitted observations: their ",
        "row identifiers are missing, duplicated, or different.",
        call. = FALSE
      )
    }

    resid <- resid[match(observation_ids, residual_ids)]
  }

  for (term in names(pred_data)) {
    term_data <- pred_data[[term]]
    term_n <- if (is.data.frame(term_data)) nrow(term_data) else NA_integer_

    if (is.na(term_n) || term_n != expected_n) {
      term_n_label <- if (is.na(term_n)) "a non-data-frame result" else term_n
      stop(
        "Cannot attach residuals for term `", term, "`: its prediction ",
        "frame contains ", term_n_label, " but ", expected_n,
        " fitted observations were expected.",
        call. = FALSE
      )
    }

    term_resid <- resid
    term_ids <- rownames(term_data)
    default_term_ids <- as.character(seq_len(expected_n))

    if (!is.null(observation_ids) && !is.null(term_ids) &&
        !identical(term_ids, observation_ids)) {
      if (!anyNA(term_ids) && !anyDuplicated(term_ids) &&
          setequal(term_ids, observation_ids)) {
        # Some grouped prediction frames preserve source row names. If such a
        # frame has been reordered, align residuals to that frame explicitly.
        term_resid <- resid[match(term_ids, observation_ids)]
      } else if (!identical(term_ids, default_term_ids)) {
        # Default sequential row names carry no observation identity and use
        # the documented x_original order. Non-default incompatible names,
        # however, are evidence that positional attachment is unsafe.
        stop(
          "Cannot attach residuals for term `", term, "`: its prediction ",
          "row identifiers do not match the fitted observations.",
          call. = FALSE
        )
      }
    }

    # Remove residual names after alignment. The prediction frame already
    # carries the row order, and retaining unrelated names on this column can
    # obscure that the assignment is deliberately positional at this point.
    term_data$resid <- unname(term_resid)
    pred_data[[term]] <- term_data
  }

  pred_data
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
    # Prediction frames and residuals must describe the same fitted rows.
    # The helper checks counts and uses observation names for alignment when
    # both the stored predictors and residual vector provide them.
    pred_data <- mfp2_plot_attach_residuals(pred_data, resid, model)

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

    is_binary <- (is_spike &&
                    spike_dec_v == saz_decision_codes[["binary_only"]]) ||
      mfp2_plot_term_is_binary(model, v)

    # A formula factor is discrete even when it has more than two levels.
    # Treat it like a binary effect for plotting purposes: predict only at the
    # fitted levels and display estimates with confidence intervals, rather
    # than drawing a continuous curve between category labels.
    is_factor <- mfp2_plot_term_is_factor(model, v)
    is_discrete <- is_binary || is_factor

    binary_plot_data <- NULL
    if (is_discrete) {
      binary_df <- mfp2_plot_binary_prediction(
        model = model,
        term = v,
        type = type,
        ref = ref,
        alpha = alpha
      )

      if (!is.null(binary_df)) {
        binary_plot_data <- mfp2_plot_prepare_binary_data(
          model = model,
          term = v,
          df = binary_df,
          residual_df = if (!partial_only) pred_data[[v]] else NULL
        )
        df <- binary_plot_data$fitted
        if (!partial_only) {
          pred_data[[v]] <- binary_plot_data$residual
        }
      }
    }

    power_label <- paste0(model$fp_powers[[v]], collapse = ", ")
    continuous_label <- if (is_acd) {
      sprintf("ACD FP(%s)", power_label)
    } else {
      sprintf("FP(%s)", power_label)
    }

    # Title with FP/ACD powers and explicit spike-at-zero representation.
    # The internal spike_dec code is not exposed in the plot title.
    # Factor variables get "Categorical" instead of a misleading FP label.
    plot_title <- if (is_factor) {
      "Categorical"
    } else if (is_spike) {
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
    if (is_discrete && !is.null(binary_plot_data)) {
      # Discrete terms (binary variables and formula factors of any size) use
      # one estimate and confidence interval per fitted level. Do not imply
      # values between unordered/ordered factor levels with a continuous line.
      bin_df <- binary_plot_data$fitted

      p <- p + ggplot2::geom_point(
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

      if (isTRUE(binary_plot_data$continuous_axis)) {
        p <- p + ggplot2::scale_x_continuous(
          breaks = binary_plot_data$breaks,
          labels = format(binary_plot_data$breaks, trim = TRUE),
          minor_breaks = NULL
        )
      } else {
        p <- p + ggplot2::scale_x_discrete(
          breaks = binary_plot_data$breaks,
          limits = binary_plot_data$breaks
        )
      }

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

#' Deprecated Fractional-Polynomial Plotting Interface
#'
#' `fracplot()` is retained temporarily for compatibility with code written
#' for earlier versions of `mfp2`. Use [plot.mfp2()] or `plot()` instead.
#'
#' The function forwards its arguments to [plot.mfp2()] and returns the same
#' plot objects.
#'
#' @inheritParams plot.mfp2
#' @param model A fitted object of class `"mfp2"`.
#'
#' @return
#' A named list of `ggplot2` objects, one for each plotted term, identical to
#' the result returned by [plot.mfp2()].
#'
#' @seealso [plot.mfp2()]
#'
#' @examples
#' data("prostate")
#'
#' fit <- mfp2(
#'   lpsa ~ fp(age) + fp(cavol) + fp(weight) + svi,
#'   data = prostate,
#'   verbose = FALSE
#' )
#'
#' # Recommended interface
#' plots <- plot(fit, terms = "age")
#' plots
#'
#' \donttest{
#' # Deprecated compatibility interface
#' old_plots <- fracplot(fit, terms = "age")
#' }
#'
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

  plot_args <- list(
    x = model,
    terms = terms,
    partial_only = partial_only,
    type = type,
    ref = ref,
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
    alpha_fill = alpha_fill
  )

  # Preserve plot.mfp2()'s type-dependent default when terms_seq was omitted.
  if (!terms_seq_missing) {
    plot_args$terms_seq <- terms_seq
  }

  do.call(plot.mfp2, plot_args)
}
