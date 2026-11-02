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
#' model. Each variable in `interaction_vars` has its own term-specific interaction
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
#' Multinomial results contain the slope coefficients for each non-reference
#' logit; their class-specific intercepts are not included.
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
#'   interaction_vars = "cavol",
#'   interaction_forms = c(cavol = "fp2"),
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
#'   interaction_vars = "cavol",
#'   interaction_forms = c(cavol = "fp2"),
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

#' Validate Shared Inputs to MFPI `coef()`/`vcov()` Accessors
#'
#' Validates the arguments shared by the MFPI `coef()` and `vcov()`
#' accessors before extracting any stored fit, so that the caller's error
#' messages consistently identify the offending argument regardless of
#' which downstream helper eventually consumes it.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character scalar naming a term, or `NULL` for all.
#' @param model Character scalar, `"best"` or `"all"`.
#' @param dots List of trailing arguments captured by the accessor.
#' @param method Character scalar naming the calling accessor (`"coef"` or
#'   `"vcov"`), used in error messages.
#'
#' @return Invisibly returns `TRUE` on success; raises a helpful error
#'   otherwise.
#'
#' @keywords internal
#' @noRd
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


#' Resolve Stored MFPI Winner Fits for Accessors
#'
#' Returns either the winning term-specific fit or every stored candidate
#' fit for the requested term(s). This is the shared entry point used by
#' the `coef()`, `vcov()`, and `print()` methods for `"mfpi"` objects.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character vector of term names, or `NULL` for every stored
#'   term.
#' @param model Character scalar, one of `"best"` (winning fit only) or
#'   `"all"` (every candidate fit).
#'
#' @return A named list of stored fit lists, keyed by term name.
#'
#' @keywords internal
#' @noRd
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


#' Coefficient Labels and Print Metadata for One MFPI Term
#'
#' Assembles the complete set of coefficient labels, per-column display
#' names, group labels, and adjustment-coefficient metadata that the MFPI
#' `coef()` and `print()` methods need to render one term-specific fit.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character scalar naming the term.
#' @param fit_result Stored fit list for `term` from
#'   `mfpi_get_accessor_fits()`.
#'
#' @return A named list carrying the fitted coefficient names,
#'   user-facing display names, per-group labels, adjustment metadata, and
#'   any auxiliary tables consumed by the print helpers.
#'
#' @keywords internal
#' @noRd
mfpi_coefficient_info <- function(object, term, fit_result) {
  # Multinomial interaction models carry a per-logit coefficient matrix rather
  # than a flat named vector, so coefficient information is assembled on a
  # dedicated path that labels each coefficient by its non-reference logit. The
  # family lookup is guarded so minimal or older objects without family metadata
  # fall through to the single-response path rather than erroring.
  family_string <- tryCatch(mfpi_family_string(object),
                            error = function(e) NA_character_)
  if (identical(family_string, "multinomial")) {
    return(mfpi_multinomial_coefficient_info(object, term, fit_result))
  }

  interaction_spec <- mfpi_get_interaction_spec(object, term)
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
  all_interaction_source_names <- unlist(
    coefficient_groups,
    use.names = FALSE
  )

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

    transformations <- if (isTRUE(interaction_spec$discrete)) {
      interaction_spec$columns
    } else {
      mfpi_fp_transformation_labels(
        object = object,
        term = term,
        powers = power_g
      )
    }

    estimates <- unname(beta[model_names])
    centers <- mfpi_interaction_center_values(
      fit_result = fit_result,
      source_names = source_names,
      all_source_names = all_interaction_source_names
    )
    group_rows[[g]] <- data.frame(
      group_level = rep(display_groups[g], length(model_names)),
      transformation = transformations,
      center = centers,
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
  adjustment_metadata <- mfpi_adjustment_coefficient_metadata(
    object = object,
    interaction_model = interaction_model,
    coefficient_names = adjustment_raw_names
  )
  adjustment_labels <- adjustment_metadata$basis

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

  adjustment_table <- adjustment_metadata
  if (length(adjustment_raw_names) > 0L) {
    adjustment_table$term <- all_display[match(adjustment_raw_names, raw_names)]
    adjustment_table$estimate <- unname(
      beta[match(adjustment_raw_names, raw_names)]
    )
  } else {
    adjustment_table$estimate <- numeric(0L)
  }

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


#' Readable Labels for Shared Multinomial Interaction Columns
#'
#' Builds readable per-column labels for the shared model columns of a
#' multinomial MFPI interaction design (intercept, group-specific FP terms,
#' group-indicator dummies, and adjustment columns). Labels mirror the
#' single-response path but omit the per-logit prefix, which the caller
#' prepends when a logit-major table is assembled.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character scalar naming the interaction term.
#' @param fit_result Stored fit list for `term`.
#' @param interaction_model Fitted multinomial interaction model.
#' @param model_cols Character vector of shared model-column names.
#'
#' @return Named character vector of labels keyed by model-column name.
#'
#' @keywords internal
#' @noRd
mfpi_multinomial_column_labels <- function(object, term, fit_result,
                                           interaction_model, model_cols) {
  interaction_spec <- mfpi_get_interaction_spec(object, term)
  labels <- stats::setNames(model_cols, model_cols)

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
  if (is.list(powers)) {
    if (!is.null(names(powers)) && setequal(names(powers), internal_groups)) {
      powers <- powers[internal_groups]
    } else if (length(powers) == length(internal_groups)) {
      names(powers) <- internal_groups
    }
  }

  interaction_cols <- character(0L)
  for (g in seq_along(coefficient_groups)) {
    model_names <- mfpi_resolve_fitted_coefficient_names(
      interaction_model = interaction_model,
      source_names = coefficient_groups[[g]],
      coefficient_names = model_cols,
      term = term,
      allow_missing = TRUE
    )
    if (length(model_names) == 0L) next
    transformations <- if (isTRUE(interaction_spec$discrete)) {
      interaction_spec$columns
    } else {
      mfpi_fp_transformation_labels(
        object = object, term = term, powers = as.numeric(powers[[g]])
      )
    }
    if (length(transformations) == length(model_names)) {
      labels[model_names] <- paste0(display_groups[g], ": ", transformations)
    }
    interaction_cols <- c(interaction_cols, model_names)
  }

  # Group-indicator dummies.
  dummy_source_names <- if (length(internal_groups) > 1L) {
    paste0(object$group_var, internal_groups[-1L])
  } else {
    character(0L)
  }
  dummy_model_names <- mfpi_resolve_fitted_coefficient_names(
    interaction_model = interaction_model,
    source_names = dummy_source_names,
    coefficient_names = model_cols,
    term = term,
    allow_missing = TRUE
  )
  if (length(dummy_model_names) > 0L) {
    labels[dummy_model_names] <- paste0(
      object$group_var, " [", display_groups[-1L][seq_along(dummy_model_names)], "]"
    )
  }

  # Intercept.
  if ("(Intercept)" %in% model_cols) labels[["(Intercept)"]] <- "(Intercept)"

  # Anything left over is an adjustment coefficient.
  adjustment_cols <- setdiff(
    model_cols, c(interaction_cols, dummy_model_names, "(Intercept)")
  )
  if (length(adjustment_cols) > 0L) {
    labels[adjustment_cols] <- mfpi_adjustment_coefficient_labels(
      object = object,
      interaction_model = interaction_model,
      coefficient_names = adjustment_cols
    )
  }

  labels
}


# Coefficient information for a multinomial MFPI interaction model. Produces a
#' Multinomial Coefficient Labels and Print Metadata for One MFPI Term
#'
#' Multinomial variant of `mfpi_coefficient_info()`. Returns a flat
#' coefficient vector and matching covariance-aligned raw names in the same
#' class-major order used by [nnet::vcov.multinom()], so [coef.mfpi()] and
#' [vcov.mfpi()] can reuse the shared extraction code unchanged. Display
#' names and a per-logit table carry the non-reference class of each
#' coefficient.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character scalar naming the term.
#' @param fit_result Stored fit list for `term`.
#'
#' @return A named list carrying coefficient names, display names, per-logit
#'   labels, and any adjustment metadata consumed by the multinomial print
#'   helpers.
#'
#' @keywords internal
#' @noRd
mfpi_multinomial_coefficient_info <- function(object, term, fit_result) {
  interaction_model <- fit_result$test_results$interaction_model
  if (is.null(interaction_model)) {
    stop(paste0("No interaction model is stored for term `", term, "`."),
         call. = FALSE)
  }
  coef_mat <- interaction_model$coefficient_matrix
  if (is.null(coef_mat) || !is.matrix(coef_mat) || is.null(rownames(coef_mat)) ||
      is.null(colnames(coef_mat))) {
    stop(
      paste0("The stored multinomial interaction model for term `", term,
             "` does not expose a labelled per-logit coefficient matrix."),
      call. = FALSE
    )
  }

  classes <- rownames(coef_mat)          # non-reference logits
  # The MFPI accessors report the fitted slope functions. Multinomial
  # intercepts are class-specific nuisance constants and are deliberately
  # omitted, matching the slope-only contract used by the other families.
  model_cols <- setdiff(colnames(coef_mat), "(Intercept)")

  col_label <- mfpi_multinomial_column_labels(
    object = object, term = term, fit_result = fit_result,
    interaction_model = interaction_model, model_cols = model_cols
  )

  interaction_spec <- mfpi_get_interaction_spec(object, term)
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
  if (!is.null(names(powers)) && setequal(names(powers), internal_groups)) {
    powers <- powers[internal_groups]
  }
  all_interaction_source_names <- unlist(coefficient_groups, use.names = FALSE)
  interaction_rows <- vector("list", length(coefficient_groups))
  interaction_model_names <- character(0L)
  for (g in seq_along(coefficient_groups)) {
    source_names <- coefficient_groups[[g]]
    model_names_g <- mfpi_resolve_fitted_coefficient_names(
      interaction_model = interaction_model,
      source_names = source_names,
      coefficient_names = model_cols,
      term = term
    )
    power_g <- as.numeric(powers[[g]])
    transformations <- if (isTRUE(interaction_spec$discrete)) {
      interaction_spec$columns
    } else {
      mfpi_fp_transformation_labels(object, term, power_g)
    }
    interaction_rows[[g]] <- data.frame(
      model_column = model_names_g,
      group_level = rep(display_groups[g], length(model_names_g)),
      transformation = transformations,
      center = mfpi_interaction_center_values(
        fit_result,
        source_names,
        all_interaction_source_names
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    interaction_model_names <- c(interaction_model_names, model_names_g)
  }
  interaction_metadata <- do.call(rbind, interaction_rows)

  dummy_source_names <- if (length(internal_groups) > 1L) {
    paste0(object$group_var, internal_groups[-1L])
  } else {
    character(0L)
  }
  dummy_model_names <- mfpi_resolve_fitted_coefficient_names(
    interaction_model = interaction_model,
    source_names = dummy_source_names,
    coefficient_names = model_cols,
    term = term,
    allow_missing = TRUE
  )
  adjustment_model_names <- setdiff(
    model_cols,
    c(interaction_model_names, dummy_model_names)
  )
  adjustment_metadata <- mfpi_adjustment_coefficient_metadata(
    object = object,
    interaction_model = interaction_model,
    coefficient_names = adjustment_model_names
  )

  # Class-major flat order matches vcov.multinom(): all columns of the first
  # non-reference logit, then the next, ... Names in the covariance matrix join
  # class and column with a single colon.
  raw_names <- as.vector(t(outer(classes, model_cols,
                                 function(a, b) paste0(a, ":", b))))
  values    <- as.vector(t(coef_mat[, model_cols, drop = FALSE]))
  display_names <- as.vector(t(outer(
    classes, model_cols,
    function(a, b) paste0(a, " | ", col_label[b])
  )))
  display_names <- mfpi_make_unique_display_names(display_names)
  beta <- stats::setNames(values, display_names)

  logit_rows <- lapply(classes, function(cls) {
    data.frame(
      class = cls,
      term = unname(col_label[model_cols]),
      estimate = unname(coef_mat[cls, model_cols]),
      raw_name = paste0(cls, ":", model_cols),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  logit_table <- do.call(rbind, logit_rows)

  # `paste0()` recycles a scalar prefix even when its other argument has length
  # zero. Return an explicit empty vector so metadata tables with no group or
  # adjustment coefficients remain valid zero-row data frames.
  prefix_class <- function(cls, coefficient_names) {
    if (length(coefficient_names) == 0L) return(character(0L))
    paste0(cls, ":", coefficient_names)
  }

  group_table <- do.call(rbind, lapply(classes, function(cls) {
    data.frame(
      class = rep(cls, nrow(interaction_metadata)),
      group_level = interaction_metadata$group_level,
      transformation = interaction_metadata$transformation,
      center = interaction_metadata$center,
      raw_name = prefix_class(cls, interaction_metadata$model_column),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  group_variable_table <- do.call(rbind, lapply(classes, function(cls) {
    data.frame(
      class = rep(cls, length(dummy_model_names)),
      level = display_groups[-1L][seq_along(dummy_model_names)],
      raw_name = prefix_class(cls, dummy_model_names),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  adjustment_table <- do.call(rbind, lapply(classes, function(cls) {
    data.frame(
      class = rep(cls, nrow(adjustment_metadata)),
      variable = adjustment_metadata$variable,
      basis = adjustment_metadata$basis,
      center = adjustment_metadata$center,
      raw_name = prefix_class(cls, adjustment_metadata$raw_name),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  power_table <- data.frame(group_level = display_groups, stringsAsFactors = FALSE)
  power_table$powers <- I(unname(powers))

  list(
    term = term,
    group_var = object$group_var,
    shift = mfpi_named_scalar(object$shift, term, default = 0),
    flex = toupper(as.character(object$flex)),
    form = mfpi_accessor_interaction_form(object, term, fit_result),
    coefficients = beta,
    raw_names = raw_names,
    display_names = display_names,
    multinomial = TRUE,
    classes = classes,
    reference_class = interaction_model$reference_class,
    logit_table = logit_table,
    group_table = group_table,
    group_variable_table = group_variable_table,
    adjustment_table = adjustment_table,
    power_table = power_table,
    reference_level = display_groups[1L]
  )
}


#' Readable Labels for Adjustment Coefficients
#'
#' Builds readable labels for adjustment coefficients carried from the
#' Step 1 MFP model into one term-specific MFPI interaction model. The
#' interaction fit stores transformed source-column names, and the
#' adjustment `mfp2` object stores the functional-form metadata needed to
#' turn those names back into variable names or FP expressions.
#'
#' @param object An `"mfpi"` model object.
#' @param ... Additional named arguments carrying the fit result and the
#'   adjustment-model reference.
#'
#' @return Named character vector of adjustment-coefficient labels.
#'
#' @keywords internal
#' @noRd
mfpi_adjustment_coefficient_labels <- function(object,
                                               interaction_model,
                                               coefficient_names) {
  metadata <- mfpi_adjustment_coefficient_metadata(
    object = object,
    interaction_model = interaction_model,
    coefficient_names = coefficient_names
  )
  metadata$basis
}


#' Adjustment-Coefficient Metadata for One MFPI Interaction Model
#'
#' Builds the variable, basis, and centering metadata for adjustment
#' coefficients carried from the Step 1 MFP model into one term-specific
#' MFPI interaction model.
#'
#' @param object An `"mfpi"` model object.
#' @param ... Additional named arguments carrying the fit result and the
#'   adjustment-model reference.
#'
#' @return A data frame with columns describing each adjustment
#'   coefficient's source variable, basis, and centering constant.
#'
#' @keywords internal
#' @noRd
mfpi_adjustment_coefficient_metadata <- function(object,
                                                 interaction_model,
                                                 coefficient_names) {
  out <- data.frame(
    term = coefficient_names,
    variable = coefficient_names,
    basis = coefficient_names,
    center = rep(NA_real_, length(coefficient_names)),
    raw_name = coefficient_names,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (length(coefficient_names) == 0L) return(out)

  source_map <- interaction_model$transformed_to_model_columns
  adj <- object$adjustment_model
  if (is.null(source_map) || !is.character(source_map) ||
      is.null(names(source_map)) || is.null(adj) || !inherits(adj, "mfp2")) {
    return(out)
  }

  source_names <- names(source_map)[match(coefficient_names, unname(source_map))]
  have_source <- !is.na(source_names) & nzchar(source_names)
  if (!any(have_source)) return(out)

  adjustment_info <- tryCatch(
    mfp2_design_column_info(adj),
    error = function(e) NULL
  )
  if (is.null(adjustment_info) || nrow(adjustment_info) == 0L) return(out)

  info_rows <- match(source_names[have_source], adjustment_info$transformed_column)
  matched <- !is.na(info_rows)
  if (!any(matched)) return(out)

  positions <- which(have_source)[matched]
  rows <- info_rows[matched]
  out$variable[positions] <- adjustment_info$variable[rows]
  out$basis[positions] <- adjustment_info$basis[rows]
  out$term[positions] <- adjustment_info$basis[rows]
  if ("center" %in% names(adjustment_info)) {
    out$center[positions] <- adjustment_info$center[rows]
  }
  out
}


#' Recover Stored MFPI Centering Constants for Interaction Columns
#'
#' Returns the stored MFPI centering constants for a subset of the
#' group-specific interaction design columns. Older or uncentered fits
#' legitimately return `NA` for missing constants; callers are expected to
#' tolerate that.
#'
#' @param fit_result Stored fit list for the interaction term.
#' @param ... Additional named arguments identifying the subset of columns
#'   to look up.
#'
#' @return Named numeric vector of centering constants, with `NA` for
#'   columns whose constant was not stored.
#'
#' @keywords internal
#' @noRd
mfpi_interaction_center_values <- function(fit_result,
                                           source_names,
                                           all_source_names) {
  centers <- fit_result$center_vals
  if (is.null(centers)) return(rep(NA_real_, length(source_names)))

  centers <- as.numeric(centers)
  center_names <- names(fit_result$center_vals)
  if (!is.null(center_names) && all(source_names %in% center_names)) {
    return(unname(centers[match(source_names, center_names)]))
  }

  if (length(centers) == length(all_source_names)) {
    return(unname(centers[match(source_names, all_source_names)]))
  }

  rep(NA_real_, length(source_names))
}


#' Map Source Design-Column Names to Fitted Coefficient Names
#'
#' Resolves the mapping from source design-column names (as stored on the
#' MFPI fit) to the fitted coefficient names actually used by the
#' interaction model. Handles the naming conventions of both single-response
#' and multinomial fitters.
#'
#' @param interaction_model Fitted interaction model.
#' @param ... Additional named arguments carrying the source column names
#'   to resolve.
#'
#' @return Named character vector mapping source names to fitted coefficient
#'   names.
#'
#' @keywords internal
#' @noRd
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


#' Readable FP Transformation Labels for One MFPI Term
#'
#' Converts one selected FP power vector for a term into the readable
#' transformation labels used in the printed coefficient block.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character scalar naming the term.
#' @param powers Numeric vector of selected FP powers.
#'
#' @return Character vector of transformation labels, one per power.
#'
#' @keywords internal
#' @noRd
mfpi_fp_transformation_labels <- function(object, term, powers) {
  shift <- mfpi_named_scalar(object$shift, term, default = 0)
  mfp2_fp_basis_labels(term = term, powers = powers, shift = shift)
}


#' User-Facing Interaction Form Label for One MFPI Term
#'
#' Resolves the user-facing interaction form label for one MFPI term (for
#' example the selected FP1 or FP2 powers per group), shown in the header
#' of the printed coefficient block.
#'
#' @param object An `"mfpi"` model object.
#' @param term Character scalar naming the term.
#' @param fit_result Stored fit list for `term`.
#'
#' @return Character scalar with the form label.
#'
#' @keywords internal
#' @noRd
mfpi_accessor_interaction_form <- function(object, term, fit_result) {
  form <- NULL
  if (!is.null(object$interaction_forms) && term %in% names(object$interaction_forms)) {
    form <- object$interaction_forms[[term]]
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


#' Make User-Facing Coefficient Labels Unique
#'
#' Ensures the user-facing coefficient labels are unique without changing
#' any label that is already distinct. Duplicated labels receive a numeric
#' `.1`, `.2`, ... suffix in the order they occur, so the printed table
#' never contains ambiguous row identifiers.
#'
#' @param x Character vector of candidate display labels.
#'
#' @return Character vector of the same length as `x`, guaranteed to be
#'   unique.
#'
#' @keywords internal
#' @noRd
mfpi_make_unique_display_names <- function(x) {
  if (!anyDuplicated(x)) return(x)

  out <- make.unique(x, sep = " [")

  # make.unique() appends only the separator and number (for example, " [1").
  # Identify generated suffixes from the values it actually changed instead of
  # parsing their text with a regular expression. This remains correct for
  # multi-digit suffixes and leaves every original display label untouched.
  changed <- out != x
  out[changed] <- paste0(out[changed], "]")
  out
}


#' Validate the `digits` Argument Used by MFPI Print Methods
#'
#' Validates the display-precision argument used by the custom MFPI
#' coefficient print methods, so that a bad value stops with a clear
#' error rather than propagating a bad `format()` call.
#'
#' @param digits Value supplied by the caller.
#'
#' @return Integer scalar with the validated digit count.
#'
#' @keywords internal
#' @noRd
mfpi_validate_accessor_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0L || digits != floor(digits)) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }
  as.integer(digits)
}


#' Format One FP Power Vector for the Coefficient Display
#'
#' Formats a single selected FP power vector as a compact string for the
#' user-facing coefficient display.
#'
#' @param powers Numeric vector of selected FP powers.
#'
#' @return Character scalar such as `"(0, 3)"`.
#'
#' @keywords internal
#' @noRd
mfpi_format_power_vector <- function(powers) {
  vals <- vapply(
    as.numeric(powers),
    function(x) format(x, trim = TRUE, scientific = FALSE, digits = 15L),
    character(1L)
  )
  paste0("(", paste(vals, collapse = ", "), ")")
}


#' Print Selected MFPI Interaction Powers Per Group
#'
#' Prints the small per-group table of selected FP powers shown in the
#' header of an MFPI coefficient block. Returns early when the table is
#' empty so callers do not need to guard the call.
#'
#' @param power_table Data frame with one row per fitted group level and
#'   columns for the selected FP powers.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
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


#' Print One Term-Specific MFPI Coefficient Block
#'
#' Prints a single term-specific coefficient block for an MFPI fit,
#' dispatching to `mfpi_print_multinomial_coef_block()` when the fit is
#' multinomial and otherwise formatting the shared single-response layout.
#'
#' @param info Coefficient-info list from `mfpi_coefficient_info()` or
#'   `mfpi_multinomial_coefficient_info()`.
#' @param digits Integer scalar controlling display precision.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of
#'   printing.
#'
#' @keywords internal
#' @noRd
mfpi_print_coef_block <- function(info, digits) {
  if (isTRUE(info$multinomial)) {
    return(mfpi_print_multinomial_coef_block(info, digits = digits))
  }

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


#' Print MFPI Coefficient Block for a Multinomial Interaction Model
#'
#' Prints the coefficient block for a multinomial MFPI interaction model.
#' Coefficients are grouped by non-reference logit; each block lists that
#' logit's coefficients against the common reference class named in the
#' header.
#'
#' @param info Multinomial coefficient-info list from
#'   `mfpi_multinomial_coefficient_info()`.
#' @param digits Integer scalar controlling display precision.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of
#'   printing.
#'
#' @keywords internal
#' @noRd
mfpi_print_multinomial_coef_block <- function(info, digits) {
  cat("Interaction: ", info$term, "\n", sep = "")
  cat("Group variable: ", info$group_var, "\n", sep = "")
  cat("FLEX: ", info$flex, "\n", sep = "")
  if (identical(info$form, "Linear")) {
    cat("Form: Linear\n")
  } else {
    cat("FP form: ", info$form, "\n", sep = "")
  }
  cat("Family: multinomial (reference class: ", info$reference_class, ")\n", sep = "")

  mfpi_print_interaction_powers(info$power_table)

  for (cls in info$classes) {
    rows <- info$logit_table[info$logit_table$class == cls, , drop = FALSE]
    tab <- rows[, c("term", "estimate"), drop = FALSE]
    names(tab) <- c("Term", "Estimate")
    cat("\nLogit ", cls, " vs ", info$reference_class, ":\n", sep = "")
    print(tab, row.names = FALSE, digits = digits)
  }

  invisible(NULL)
}
