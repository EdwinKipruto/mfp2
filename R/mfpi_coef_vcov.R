# -----------------------------------------------------------------------------
# Coefficient and covariance accessors for MFPI models -------------------------
# -----------------------------------------------------------------------------

#' Extract Coefficients from MFPI Interaction Models
#'
#' Extracts regression coefficients from the term-specific interaction models
#' stored in an [mfpi()] fit. Because MFPI evaluates each continuous variable in
#' a separate interaction model, the method returns one coefficient vector when
#' a single interaction is requested and a named collection when several
#' interactions are requested.
#'
#' @details
#' `mfpi()` does not fit all selected interactions together in one joint final
#' model. Each variable in `cont_vars` has its own term-specific interaction
#' model. Consequently, `coef.mfpi()` reports coefficients from those stored
#' interaction models rather than from the first-stage MFP adjustment model.
#'
#' By default, only interactions retained by the MFPI selection criterion are
#' available (`model = "best"`). Set `model = "all"` to extract the fitted
#' interaction model for every evaluated variable, including interactions that
#' were not retained. Use `term` to select one or more continuous variables.
#'
#' The printed representation is designed for interpretation. For each
#' interaction it reports the interaction variable, grouping variable, FLEX
#' setting, and requested interaction form. The interaction-model FP powers are
#' then shown separately for every group level, irrespective of FLEX. This gives
#' one consistent display for FLEX1--FLEX4 and makes group-specific power choices
#' under FLEX4 explicit.
#'
#' The coefficients belonging to the group-specific continuous-variable
#' functions are displayed in a table with columns for the original group level,
#' the FP transformation, and the coefficient estimate. Coefficients for the
#' grouping variable are shown in a separate table together with the fitted
#' reference level. The table reports the coefficient attached to each
#' non-reference factor level; it is deliberately not labelled as a direct
#' comparison with the reference group because, in the presence of an
#' interaction, that coefficient is only one component of the model
#' parameterization and is not an overall group effect. For a grouping variable
#' with K levels, there are K group-specific FP blocks and K - 1 group-variable
#' coefficients. Adjustment coefficients are displayed separately, using stored
#' MFP metadata to recover readable variable or transformation labels where
#' available. A model intercept, when present, is also shown separately.
#'
#' Transformation labels identify the FP basis before centering. When the MFPI
#' fit uses `center = TRUE`, the corresponding fitted design columns are centered
#' using the stored MFPI centering constants; the coefficient estimates are
#' therefore coefficients of those centered design columns.
#'
#' The numeric values remain ordinary regression coefficients. For a single
#' requested interaction, the returned object is a numeric vector with an MFPI
#' print class and can be used in numerical calculations. If several
#' interactions are requested, a named list of such vectors is returned.
#'
#' The transformation labels reflect the interaction model actually fitted.
#' Thus, under `flex4`, different group levels can display different FP
#' transformations. Repeated FP powers use the standard repeated-power form;
#' for example, powers `c(-2, -2)` are displayed as `x^-2` and
#' `x^-2 * log(x)` (with any fitted shift included in `x`).
#'
#' @param object A fitted [mfpi()] object.
#' @param term Optional character vector naming continuous variables whose
#'   term-specific interaction models should be extracted. When `NULL`, all
#'   available terms in the requested `model` scope are returned.
#' @param model Character scalar indicating which stored interaction models are
#'   eligible. `"best"` (the default) uses interactions retained by the MFPI
#'   selection criterion; `"all"` uses every evaluated term-specific winner,
#'   whether retained or not.
#' @param ... Not used. Supplying additional arguments produces an error.
#'
#' @return
#' If one interaction model is requested, a named numeric vector of regression
#' coefficients with class `"mfpi_coef"`. If several interaction models are
#' requested, a named list with class `"mfpi_coef_list"`, with one coefficient
#' vector per interaction. The printed labels are user-facing; the original
#' fitted coefficient names are retained in the `raw_names` attribute of each
#' coefficient vector.
#'
#' @examples
#' data("prostate", package = "mfp2")
#'
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(cavol),
#'   data = prostate,
#'   group_var = "svi",
#'   cont_vars = "cavol",
#'   cont_var_forms = c(cavol = "fp2"),
#'   flex = "flex4",
#'   p_interact = 1,
#'   verbose = FALSE
#' )
#'
#' # User-facing coefficient display for the selected interaction model.
#' coef(fit)
#'
#' # Extract one term explicitly.
#' b <- coef(fit, term = "cavol")
#' is.numeric(b)
#'
#' # Inspect an evaluated interaction even when it was not retained.
#' coef(fit, term = "cavol", model = "all")
#'
#' @seealso [mfpi()], [vcov.mfpi()], [summary.mfpi()], [predict.mfpi()]
#'
#' @method coef mfpi
#' @export
coef.mfpi <- function(object,
                      term = NULL,
                      model = c("best", "all"),
                      ...) {
  mfpi_check_accessor_input(object, term, model, list(...), "coef")
  model <- match.arg(model)

  fits <- mfpi_get_accessor_fits(object, term = term, model = model)
  out <- lapply(names(fits), function(term_name) {
    info <- mfpi_coefficient_info(object, term_name, fits[[term_name]])
    beta <- info$coefficients
    class(beta) <- c("mfpi_coef", "numeric")
    attr(beta, "mfpi_info") <- info
    attr(beta, "raw_names") <- info$raw_names
    beta
  })
  names(out) <- names(fits)

  if (length(out) == 1L) {
    return(out[[1L]])
  }

  structure(out, class = c("mfpi_coef_list", "list"))
}


#' Extract Variance-Covariance Matrices from MFPI Interaction Models
#'
#' Extracts estimated variance-covariance matrices from the term-specific
#' interaction models stored in an [mfpi()] fit. The method uses the same model
#' scope, term selection, coefficient order, and user-facing labels as
#' [coef.mfpi()].
#'
#' @details
#' Each continuous variable evaluated by `mfpi()` has a separate interaction
#' model. There is therefore no single covariance matrix spanning coefficients
#' from different MFPI interactions. When one term is requested, `vcov.mfpi()`
#' returns that model's ordinary covariance matrix. When several terms are
#' requested, it returns a named list containing one covariance matrix per
#' term-specific model.
#'
#' Row and column names correspond exactly to `names(coef(object, ...))` for the
#' same interaction model. Group-specific FP coefficients are labelled using the
#' original group levels and readable FP transformations. Group-variable
#' coefficients use labels such as `treatment [Drug A]`; the reference level is
#' shown by the `coef()` print method rather than encoded as a pairwise contrast
#' in the coefficient name.
#'
#' The covariance matrices are conditional on the functional forms and model
#' selected by MFPI. They do not incorporate the additional uncertainty caused
#' by FP power selection, adjustment-model selection, or interaction selection.
#'
#' @inheritParams coef.mfpi
#'
#' @return
#' If one interaction model is requested, a numeric square matrix with row and
#' column names matching [coef.mfpi()]. If several interaction models are
#' requested, a named list with class `"mfpi_vcov_list"`, containing one matrix
#' per interaction. The original fitted coefficient names are retained in the
#' `raw_names` attribute of each matrix.
#'
#' @examples
#' data("prostate", package = "mfp2")
#'
#' fit <- mfpi(
#'   lpsa ~ fp(age) + svi + fp(cavol),
#'   data = prostate,
#'   group_var = "svi",
#'   cont_vars = "cavol",
#'   cont_var_forms = c(cavol = "fp2"),
#'   flex = "flex4",
#'   p_interact = 1,
#'   verbose = FALSE
#' )
#'
#' V <- vcov(fit, term = "cavol")
#' dim(V)
#' identical(rownames(V), names(coef(fit, term = "cavol")))
#'
#' # `model = "all"` is useful when the evaluated interaction was not retained.
#' vcov(fit, term = "cavol", model = "all")
#'
#' @seealso [coef.mfpi()], [mfpi()], [summary.mfpi()], [predict.mfpi()]
#'
#' @method vcov mfpi
#' @export
vcov.mfpi <- function(object,
                      term = NULL,
                      model = c("best", "all"),
                      ...) {
  mfpi_check_accessor_input(object, term, model, list(...), "vcov")
  model <- match.arg(model)

  fits <- mfpi_get_accessor_fits(object, term = term, model = model)
  out <- lapply(names(fits), function(term_name) {
    fit_result <- fits[[term_name]]
    info <- mfpi_coefficient_info(object, term_name, fit_result)
    fit_obj <- fit_result$test_results$interaction_model$fit
    V <- stats::vcov(fit_obj)

    raw <- info$raw_names
    missing_rows <- setdiff(raw, rownames(V))
    missing_cols <- setdiff(raw, colnames(V))
    if (length(missing_rows) > 0L || length(missing_cols) > 0L) {
      stop(
        paste0(
          "Stored coefficient names for term `", term_name,
          "` do not align with its covariance matrix."
        ),
        call. = FALSE
      )
    }

    V <- V[raw, raw, drop = FALSE]
    dimnames(V) <- list(info$display_names, info$display_names)
    attr(V, "raw_names") <- raw
    attr(V, "mfpi_term") <- term_name
    V
  })
  names(out) <- names(fits)

  if (length(out) == 1L) {
    return(out[[1L]])
  }

  structure(out, class = c("mfpi_vcov_list", "list"))
}


#' Print MFPI Coefficients
#'
#' Prints coefficients returned by [coef.mfpi()] in a user-facing form. For
#' each interaction model, the selected interaction powers are shown for every
#' group level, followed by separate tables for group-specific FP terms,
#' group-variable coefficients, and adjustment coefficients. These methods are
#' called automatically and are not normally invoked directly.
#'
#' @param x An object returned by [coef.mfpi()].
#' @param digits Number of significant digits used for coefficient estimates.
#' @param ... Not used.
#'
#' @return Invisibly returns `x`.
#'
#' @method print mfpi_coef
#' @export
print.mfpi_coef <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    mfpi_check_unused_dots(dots, "`print.mfpi_coef()`")
  }
  digits <- mfpi_validate_accessor_digits(digits)

  info <- attr(x, "mfpi_info", exact = TRUE)
  if (is.null(info)) {
    return(NextMethod("print"))
  }

  cat("MFPI coefficients\n\n")
  mfpi_print_coef_block(info, digits = digits)
  invisible(x)
}


#' @rdname print.mfpi_coef
#' @method print mfpi_coef_list
#' @export
print.mfpi_coef_list <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    mfpi_check_unused_dots(dots, "`print.mfpi_coef_list()`")
  }
  digits <- mfpi_validate_accessor_digits(digits)

  cat("MFPI coefficients\n")
  for (i in seq_along(x)) {
    cat("\n")
    info <- attr(x[[i]], "mfpi_info", exact = TRUE)
    if (is.null(info)) {
      cat("Interaction: ", names(x)[i], "\n", sep = "")
      print(unclass(x[[i]]), digits = digits)
    } else {
      mfpi_print_coef_block(info, digits = digits)
    }
  }
  invisible(x)
}


#' Print a Collection of MFPI Variance-Covariance Matrices
#'
#' Prints a compact index when [vcov.mfpi()] returns covariance matrices for
#' several term-specific interaction models. Individual matrices can be
#' displayed by requesting a single `term`.
#'
#' @param x A multi-model object returned by [vcov.mfpi()].
#' @param ... Not used.
#'
#' @return Invisibly returns `x`.
#'
#' @method print mfpi_vcov_list
#' @export
print.mfpi_vcov_list <- function(x, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    mfpi_check_unused_dots(dots, "`print.mfpi_vcov_list()`")
  }

  npar <- vapply(x, nrow, integer(1L))
  tab <- data.frame(
    Interaction = names(x),
    Parameters = npar,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cat("MFPI variance-covariance matrices\n\n")
  print(tab, row.names = FALSE)
  cat("\nUse `vcov(object, term = \"<interaction>\")` to display one matrix.\n")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Internal helpers -------------------------------------------------------------
# -----------------------------------------------------------------------------

# Validate shared coef()/vcov() inputs before extracting any stored model.
mfpi_check_accessor_input <- function(object, term, model, dots, method) {
  if (!inherits(object, "mfpi")) {
    stop("`object` must be an object of class \"mfpi\".", call. = FALSE)
  }

  if (!is.character(model) || length(model) < 1L || anyNA(model)) {
    stop("`model` must be either \"best\" or \"all\".", call. = FALSE)
  }

  if (!is.null(term) &&
      (!is.character(term) || length(term) == 0L || anyNA(term) ||
       any(!nzchar(term)))) {
    stop("`term` must be a non-empty character vector or NULL.", call. = FALSE)
  }

  if (length(dots) > 0L) {
    mfpi_check_unused_dots(dots, paste0("`", method, ".mfpi()`"))
  }

  invisible(TRUE)
}


# Resolve selected or all stored term-specific MFPI winner fits.
mfpi_get_accessor_fits <- function(object, term = NULL, model = c("best", "all")) {
  model <- match.arg(model)

  winners <- object$var_winners
  if (is.null(winners) || !is.list(winners)) {
    stop("The MFPI object does not contain stored term-specific fits.", call. = FALSE)
  }

  has_fit <- vapply(winners, function(w) {
    !is.null(w) && !is.null(w$fit) &&
      !is.null(w$fit$test_results) &&
      !is.null(w$fit$test_results$interaction_model) &&
      !is.null(w$fit$test_results$interaction_model$fit)
  }, logical(1L))

  available_all <- names(winners)[has_fit]
  available_all <- available_all[!is.na(available_all) & nzchar(available_all)]

  if (identical(model, "best")) {
    available <- names(object$best_interaction_model)
    available <- available[!is.na(available) & nzchar(available)]
    available <- intersect(available_all, available)
  } else {
    available <- available_all
  }

  if (length(available) == 0L) {
    if (identical(model, "best")) {
      stop(
        paste0(
          "No retained MFPI interaction models are available. ",
          "Use `model = \"all\"` to extract evaluated interaction models."
        ),
        call. = FALSE
      )
    }
    stop("No evaluated MFPI interaction models are available.", call. = FALSE)
  }

  if (!is.null(term)) {
    missing_terms <- setdiff(term, available)
    if (length(missing_terms) > 0L) {
      suffix <- if (identical(model, "best")) {
        " Use `model = \"all\"` to access evaluated but unselected interactions."
      } else {
        ""
      }
      stop(
        paste0(
          "The following term(s) do not have ", model,
          " MFPI interaction models: ", paste(missing_terms, collapse = ", "),
          ".", suffix
        ),
        call. = FALSE
      )
    }
    available <- unique(term)
  }

  out <- stats::setNames(vector("list", length(available)), available)
  for (term_name in available) {
    out[[term_name]] <- winners[[term_name]]$fit
  }
  out
}


# Build all coefficient labels and print metadata for one term-specific fit.
mfpi_coefficient_info <- function(object, term, fit_result) {
  interaction_model <- fit_result$test_results$interaction_model
  fit_obj <- interaction_model$fit
  beta <- interaction_model$coefficients
  if (is.null(beta)) {
    beta <- stats::coef(fit_obj)
  }

  if (!is.numeric(beta) || is.null(names(beta))) {
    stop(
      paste0("The stored interaction model for term `", term,
             "` does not expose a named numeric coefficient vector."),
      call. = FALSE
    )
  }

  raw_names <- names(beta)
  coefficient_groups <- fit_result$coefficient_groups
  if (is.null(coefficient_groups)) {
    coefficient_groups <- attr(fit_result$xinteraction, "column_groups")
  }
  mfpi_validate_coefficient_groups(coefficient_groups, term = term)

  internal_groups <- names(coefficient_groups)
  if (is.null(internal_groups) || anyNA(internal_groups) ||
      any(!nzchar(internal_groups))) {
    internal_groups <- as.character(seq_along(coefficient_groups) - 1L)
    names(coefficient_groups) <- internal_groups
  }
  display_groups <- mfpi_prediction_group_display_labels(object, internal_groups)

  powers <- fit_result$bestfp_interaction
  if (!is.list(powers) || length(powers) != length(coefficient_groups)) {
    stop(
      paste0("Term `", term,
             "` does not contain a valid group-specific interaction power map."),
      call. = FALSE
    )
  }
  if (!is.null(names(powers))) {
    if (setequal(names(powers), internal_groups)) {
      powers <- powers[internal_groups]
    } else if (length(powers) == length(internal_groups)) {
      names(powers) <- internal_groups
    }
  } else {
    names(powers) <- internal_groups
  }

  group_rows <- vector("list", length(coefficient_groups))
  interaction_raw_names <- character(0L)

  for (g in seq_along(coefficient_groups)) {
    source_names <- coefficient_groups[[g]]
    model_names <- mfpi_resolve_fitted_coefficient_names(
      interaction_model = interaction_model,
      source_names = source_names,
      coefficient_names = raw_names,
      term = term
    )

    power_g <- as.numeric(powers[[g]])
    if (length(power_g) != length(model_names) || anyNA(power_g)) {
      stop(
        paste0("Interaction powers for term `", term,
               "` do not align with its fitted group-specific coefficients."),
        call. = FALSE
      )
    }

    transformations <- mfpi_fp_transformation_labels(
      object = object,
      term = term,
      powers = power_g
    )

    estimates <- unname(beta[model_names])
    group_rows[[g]] <- data.frame(
      group_level = rep(display_groups[g], length(model_names)),
      transformation = transformations,
      estimate = estimates,
      raw_name = model_names,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    interaction_raw_names <- c(interaction_raw_names, model_names)
  }

  group_table <- do.call(rbind, group_rows)
  rownames(group_table) <- NULL

  group_var <- object$group_var
  reference_label <- display_groups[1L]
  dummy_source_names <- if (length(internal_groups) > 1L) {
    paste0(group_var, internal_groups[-1L])
  } else {
    character(0L)
  }
  dummy_model_names <- if (length(dummy_source_names) > 0L) {
    mfpi_resolve_fitted_coefficient_names(
      interaction_model = interaction_model,
      source_names = dummy_source_names,
      coefficient_names = raw_names,
      term = term
    )
  } else {
    character(0L)
  }

  dummy_labels <- stats::setNames(character(length(dummy_model_names)), dummy_model_names)
  if (length(dummy_model_names) > 0L) {
    # Do not label these as pairwise contrasts. With an interaction in the
    # model, a group-indicator coefficient is only one component of the fitted
    # parameterization and is not an overall comparison with the reference.
    dummy_labels[] <- paste0(group_var, " [", display_groups[-1L], "]")
  }

  remaining_raw_names <- raw_names[!raw_names %in% c(interaction_raw_names, dummy_model_names)]
  intercept_raw_name <- intersect("(Intercept)", remaining_raw_names)
  adjustment_raw_names <- setdiff(remaining_raw_names, intercept_raw_name)
  adjustment_labels <- mfpi_adjustment_coefficient_labels(
    object = object,
    interaction_model = interaction_model,
    coefficient_names = adjustment_raw_names
  )

  all_display <- raw_names
  all_display[match(interaction_raw_names, raw_names)] <- paste0(
    group_table$group_level, ": ", group_table$transformation
  )
  if (length(dummy_model_names) > 0L) {
    all_display[match(dummy_model_names, raw_names)] <- unname(dummy_labels[dummy_model_names])
  }
  if (length(adjustment_raw_names) > 0L) {
    all_display[match(adjustment_raw_names, raw_names)] <- adjustment_labels
  }

  all_display <- mfpi_make_unique_display_names(all_display)
  names(beta) <- all_display

  group_variable_table <- data.frame(
    level = if (length(dummy_model_names) > 0L) display_groups[-1L] else character(0L),
    estimate = if (length(dummy_model_names) > 0L) unname(beta[match(dummy_model_names, raw_names)]) else numeric(0L),
    raw_name = dummy_model_names,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  adjustment_table <- data.frame(
    term = if (length(adjustment_raw_names) > 0L) all_display[match(adjustment_raw_names, raw_names)] else character(0L),
    estimate = if (length(adjustment_raw_names) > 0L) unname(beta[match(adjustment_raw_names, raw_names)]) else numeric(0L),
    raw_name = adjustment_raw_names,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  intercept_table <- data.frame(
    term = if (length(intercept_raw_name) > 0L) "(Intercept)" else character(0L),
    estimate = if (length(intercept_raw_name) > 0L) unname(beta[match(intercept_raw_name, raw_names)]) else numeric(0L),
    raw_name = intercept_raw_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  power_table <- data.frame(
    group_level = display_groups,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  power_table$powers <- I(unname(powers))

  form <- mfpi_accessor_interaction_form(object, term, fit_result)
  shift <- mfpi_named_scalar(object$shift, term, default = 0)

  list(
    term = term,
    group_var = group_var,
    shift = shift,
    flex = toupper(as.character(object$flex)),
    form = form,
    coefficients = beta,
    raw_names = raw_names,
    display_names = all_display,
    power_table = power_table,
    group_table = group_table,
    reference_level = reference_label,
    group_variable_table = group_variable_table,
    intercept_table = intercept_table,
    adjustment_table = adjustment_table
  )
}


# Build readable labels for adjustment coefficients carried into an interaction
# model. The interaction fit stores transformed source-column names, while the
# adjustment mfp2 object stores the functional-form metadata needed to turn
# those names back into variables or FP expressions.
mfpi_adjustment_coefficient_labels <- function(object,
                                               interaction_model,
                                               coefficient_names) {
  if (length(coefficient_names) == 0L) return(character(0L))

  labels <- coefficient_names
  source_map <- interaction_model$transformed_to_model_columns
  adj <- object$adjustment_model

  if (is.null(source_map) || !is.character(source_map) ||
      is.null(names(source_map)) || is.null(adj) || !inherits(adj, "mfp2")) {
    return(labels)
  }

  source_names <- names(source_map)[match(coefficient_names, unname(source_map))]
  have_source <- !is.na(source_names) & nzchar(source_names)
  if (!any(have_source)) return(labels)

  adjustment_info <- tryCatch(
    mfp2_design_column_info(adj),
    error = function(e) NULL
  )
  if (is.null(adjustment_info) || nrow(adjustment_info) == 0L) {
    return(labels)
  }

  info_rows <- match(source_names[have_source], adjustment_info$transformed_column)
  matched <- !is.na(info_rows)
  if (any(matched)) {
    label_positions <- which(have_source)[matched]
    labels[label_positions] <- adjustment_info$basis[info_rows[matched]]
  }

  labels
}


# Map source design-column names to fitted coefficient names.
mfpi_resolve_fitted_coefficient_names <- function(interaction_model,
                                                  source_names,
                                                  coefficient_names,
                                                  term,
                                                  allow_missing = FALSE) {
  if (length(source_names) == 0L) return(character(0L))

  out <- source_names
  direct <- source_names %in% coefficient_names

  mapping <- interaction_model$transformed_to_model_columns
  if (!all(direct) && !is.null(mapping) && !is.null(names(mapping))) {
    idx <- match(source_names[!direct], names(mapping))
    mapped <- !is.na(idx)
    if (any(mapped)) {
      out[which(!direct)[mapped]] <- unname(mapping[idx[mapped]])
    }
  }

  present <- out %in% coefficient_names
  if (!all(present) && !isTRUE(allow_missing)) {
    stop(
      paste0(
        "Could not align the following stored coefficient columns for term `",
        term, "`: ", paste(source_names[!present], collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  out[present]
}


# Convert one selected FP power vector into readable transformation labels.
mfpi_fp_transformation_labels <- function(object, term, powers) {
  shift <- mfpi_named_scalar(object$shift, term, default = 0)
  mfp2_fp_basis_labels(term = term, powers = powers, shift = shift)
}


# Resolve the user-facing interaction form for one term.
mfpi_accessor_interaction_form <- function(object, term, fit_result) {
  form <- NULL
  if (!is.null(object$cont_var_forms) && term %in% names(object$cont_var_forms)) {
    form <- object$cont_var_forms[[term]]
  }
  if (is.null(form) || length(form) != 1L || is.na(form) || !nzchar(form)) {
    winner <- object$var_winners[[term]]
    if (!is.null(winner$type)) form <- winner$type
  }
  if (is.null(form) || length(form) != 1L || is.na(form) || !nzchar(form)) {
    n_basis <- length(fit_result$bestfp_interaction[[1L]])
    form <- if (n_basis <= 1L) "fp1" else "fp2"
  }

  form <- tolower(as.character(form))
  switch(
    form,
    linear = "Linear",
    fp1 = "FP1",
    fp2 = "FP2",
    toupper(form)
  )
}


# Make user-facing coefficient labels unique without changing ordinary labels.
mfpi_make_unique_display_names <- function(x) {
  if (!anyDuplicated(x)) return(x)
  out <- make.unique(x, sep = " [")
  sub("( \\[0-9]+)$", "\\1]", out)
}


# Validate display precision used by custom coefficient print methods.
mfpi_validate_accessor_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0L || digits != floor(digits)) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }
  as.integer(digits)
}


# Format one FP power vector for the user-facing coefficient display.
mfpi_format_power_vector <- function(powers) {
  vals <- vapply(
    as.numeric(powers),
    function(x) format(x, trim = TRUE, scientific = FALSE, digits = 15L),
    character(1L)
  )
  paste0("(", paste(vals, collapse = ", "), ")")
}


# Print selected interaction powers for every fitted group level.
mfpi_print_interaction_powers <- function(power_table) {
  if (is.null(power_table) || nrow(power_table) == 0L) return(invisible(NULL))

  labels <- paste0(power_table$group_level, ":")
  width <- max(nchar(labels), na.rm = TRUE)
  power_text <- vapply(power_table$powers, mfpi_format_power_vector, character(1L))

  cat("\nInteraction powers:\n")
  for (i in seq_len(nrow(power_table))) {
    padded <- format(labels[i], width = width, justify = "left")
    cat("  ", padded, " ", power_text[i], "\n", sep = "")
  }
  invisible(NULL)
}


# Print one term-specific coefficient block.
mfpi_print_coef_block <- function(info, digits) {
  cat("Interaction: ", info$term, "\n", sep = "")
  cat("Group variable: ", info$group_var, "\n", sep = "")
  cat("FLEX: ", info$flex, "\n", sep = "")
  if (identical(info$form, "Linear")) {
    cat("Form: Linear\n")
  } else {
    cat("FP form: ", info$form, "\n", sep = "")
  }

  mfpi_print_interaction_powers(info$power_table)

  tab_group <- info$group_table[, c("group_level", "transformation", "estimate"), drop = FALSE]
  names(tab_group) <- c("Group level", "Transformation", "Estimate")

  cat("\nGroup-specific FP terms:\n")
  print(tab_group, row.names = FALSE, digits = digits)

  if (!is.null(info$group_variable_table) && nrow(info$group_variable_table) > 0L) {
    tab_group_var <- info$group_variable_table[, c("level", "estimate"), drop = FALSE]
    names(tab_group_var) <- c("Level", "Estimate")
    cat("\nGroup-variable coefficients:\n")
    cat("Reference level: ", info$reference_level, "\n\n", sep = "")
    print(tab_group_var, row.names = FALSE, digits = digits)
  }

  if (!is.null(info$intercept_table) && nrow(info$intercept_table) > 0L) {
    tab_intercept <- info$intercept_table[, c("term", "estimate"), drop = FALSE]
    names(tab_intercept) <- c("Term", "Estimate")
    cat("\nModel intercept:\n")
    print(tab_intercept, row.names = FALSE, digits = digits)
  }

  if (!is.null(info$adjustment_table) && nrow(info$adjustment_table) > 0L) {
    tab_adjustment <- info$adjustment_table[, c("term", "estimate"), drop = FALSE]
    names(tab_adjustment) <- c("Term", "Estimate")
    cat("\nAdjustment coefficients:\n")
    print(tab_adjustment, row.names = FALSE, digits = digits)
  }

  invisible(NULL)
}
