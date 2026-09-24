# =============================================================================
# summary_mfp2.R
#
# The structured summary method for `mfp2` model fits, and its dedicated
# `print` method. Provides an MFP-aware alternative to the raw
# model-specific summary output that ships from the underlying fitter,
# including `summary.survreg()` for parametric survival models.
#
# The file is organised in three sections:
#
#   1. summary.mfp2()               - the exported S3 method that returns an
#                                     object of class "summary.mfp2".
#   2. print.summary.mfp2()         - the exported S3 print method that
#                                     renders that object.
#   3. Internal helpers             - the shared computation and formatting
#                                     helpers, some of which are also used by
#                                     print.mfp2()'s Model Fit block so the
#                                     two methods stay in exact agreement.
#
# =============================================================================
#' Summarize an `mfp2` Model Fit
#'
#' Produces a structured, MFP-aware summary of a fitted [mfp2()] model. Unlike
#' the raw model-specific summary output, this
#' method separates interpretable linear terms from fractional-polynomial (FP)
#' and other nonlinear terms, whose individual coefficients are curve
#' parameters rather than per-unit effect sizes.
#'
#' For multinomial models, the structured summary reports the common selected
#' functional forms and an outcome-specific coefficient table containing the
#' estimate, standard error, z statistic, and p-value for each non-reference
#' logit. Model-fit values are reported as minus twice the log likelihood.
#'
#' @details
#' The summary is organised into the following sections:
#' \itemize{
#'   \item \strong{Selection Overview}: a compact overview of the selected
#'     functional form for each variable.
#'   \item \strong{Ordinal Intercepts}: for ordinal models, the estimated
#'     threshold intercepts with standard errors, Wald statistics, p-values,
#'     and confidence intervals.
#'   \item \strong{Linear Terms}: variables entering the model as a single
#'     linear term (including binary-only spike variables and factor levels),
#'     with the coefficient, standard error, test statistic, p-value, and
#'     confidence interval. GLM coefficients are always reported on their
#'     fitted link scale; supported survival models additionally report the
#'     usual exponentiated effect when it has a standard interpretation.
#'   \item \strong{Nonlinear Terms}: one row per variable modelled by an FP1,
#'     FP2, ACD, or spike/catzero compound form. Each row reports a joint
#'     likelihood-ratio test (LRT) of all of that variable's terms, using the
#'     package's selection-adjusted degrees of freedom. The LRT compares the
#'     final MFP model with an otherwise unchanged model that omits the variable.
#'   \item \strong{Model Fit}: full-linear and final-MFP deviances for GLMs,
#'     or minus twice the fitted likelihood for survival models (partial
#'     likelihood for Cox and Fine--Gray), with model degrees of freedom.
#' }
#'
#' The joint LRT is a diagnostic: the variable was selected by the MFP
#' procedure, not by this test. The degrees of freedom follow the
#' Royston--Sauerbrei convention used throughout the package (2 df per FP1,
#' 4 df per FP2), so the reported p-values acknowledge the FP power search
#' rather than treating the powers as fixed in advance.
#'
#' Coefficients for nonlinear terms are omitted from the default output. Set
#' \code{formulas = TRUE} to append the fitted-function formulas, or
#' \code{basis = TRUE} to append the individual fitted basis coefficients. For
#' an active ACD component, either option also prints the fitted definition of
#' \eqn{A(x)}, including the scale used during fitting. Set
#' \code{raw = TRUE} to obtain the underlying model-specific summary object
#' instead.
#'
#' The structured header names the fitted model and its link where applicable.
#' It also reports interpretation-critical nuisance parameters: the
#' distribution, fixed or estimated scale, scale strata, and censoring pattern
#' for `survreg` models; estimated or fixed dispersion for Gaussian, Gamma, and
#' inverse-Gaussian GLMs; and the estimated theta for negative-binomial models.
#' Ordinal and multinomial headers additionally report unweighted frequencies
#' of the responses used for fitting.
#'
#' @section Coefficient statistics:
#' The linear-terms table reports the coefficient statistics from the applicable
#' GLM, Cox, or `survreg` summary. For a robust Cox fit, it uses the robust
#' standard error so that the estimate, test statistic, and p-value remain
#' consistent. When these statistics cannot be read directly from a compatible
#' model summary, they are calculated from the fitted coefficients and
#' [stats::vcov()] when possible.
#'
#' @param object A fitted [mfp2()] object.
#' @param formulas Logical. If \code{TRUE}, append the fitted-function formula
#'   for each nonlinear variable. Default \code{FALSE}.
#' @param basis Logical. If \code{TRUE}, append a table of the individual fitted
#'   basis coefficients for nonlinear terms. Default \code{FALSE}.
#' @param raw Logical. If \code{TRUE}, bypass the structured summary and return
#'   the raw model-specific summary object, with
#'   the fitting call replaced by the original [mfp2()] call. Default
#'   \code{FALSE}.
#' @param notes Logical. If \code{FALSE}, suppress the explanatory model-df note
#'   when the structured summary is printed. Default \code{TRUE}.
#' @param digits Number of decimal places used when printing numeric output
#'   other than fractional-polynomial powers and degrees of freedom. Defaults
#'   to 3.
#' @param ... Further arguments. When \code{raw = TRUE}, passed to the
#'   underlying summary method; otherwise ignored.
#'
#' @return
#' When \code{raw = FALSE}, an object of class \code{"summary.mfp2"}: a list
#' with components \code{call}, \code{family}, \code{criterion},
#' \code{converged}, \code{n}, \code{nevents}, model-specific components
#' \code{distribution}, \code{scale}, \code{scale_fixed}, \code{dispersion},
#' \code{dispersion_fixed},
#' \code{scale_strata}, \code{distribution_parameters}, \code{censoring},
#' \code{link}, \code{theta}, and \code{response_frequencies}, plus \code{function_table},
#' \code{linear_terms}, \code{nonlinear_terms}, \code{basis} (or \code{NULL}),
#' \code{formulas} (or \code{NULL}), \code{acd_definitions} (or
#' \code{NULL}), \code{ordinal_intercepts} (a threshold-inference data frame
#' for ordinal models, otherwise \code{NULL}), \code{fit},
#' \code{raw_summary}, and \code{notes}. A dedicated \code{print} method
#' renders these. When \code{raw = TRUE}, the underlying summary object.
#'
#' @examples
#' data("prostate")
#'
#' fit <- mfp2(
#'   lpsa ~ fp(age) + svi,
#'   data = prostate,
#'   select = 1,
#'   alpha = 1,
#'   verbose = FALSE
#' )
#'
#' fit_summary <- summary(fit)
#'
#' # The standardized coefficient statistics for selected linear terms.
#' fit_summary$linear_terms
#'
#' # The original model-specific coefficient table remains available.
#' fit_summary$raw_summary$coefficients
#'
#' @seealso [mfp2()], [print.mfp2()], [stats::summary.glm()],
#'   [survival::summary.coxph()], [survival::summary.survreg()]
#'
#' @export
summary.mfp2 <- function(object,
                         formulas = FALSE,
                         basis = FALSE,
                         raw = FALSE,
                         notes = TRUE,
                         digits = 3L,
                         ...) {
  if (!inherits(object, "mfp2")) {
    stop("The object is not an mfp2 object.", call. = FALSE)
  }

  validate_logical_vector(notes, "notes", allowed_lengths = 1L)

  # ---------------------------------------------------------------------------
  # raw = TRUE: fall back to the underlying model's summary method.
  # ---------------------------------------------------------------------------
  if (isTRUE(raw)) {
    result <- NextMethod_summary(object, ...)
    if (!is.null(object$call_mfp)) {
      result$call <- object$call_mfp
    }
    return(result)
  }

  if (identical(object$family_string, "multinomial")) {
    base_object <- object
    class(base_object) <- setdiff(class(base_object), "mfp2")
    raw_summary <- summary(base_object)
    coefficients <- object$mfp2_coefficient_matrix
    standard_errors <- raw_summary$standard.errors
    if (is.null(dim(standard_errors))) {
      standard_errors <- matrix(
        standard_errors,
        nrow = nrow(coefficients),
        dimnames = dimnames(coefficients)
      )
    }
    z <- coefficients / standard_errors
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    coefficient_table <- do.call(rbind, lapply(seq_len(nrow(coefficients)), function(i) {
      data.frame(
        outcome = rownames(coefficients)[[i]],
        reference = object$reference_class,
        term = colnames(coefficients),
        coefficient = unname(coefficients[i, ]),
        se = unname(standard_errors[i, ]),
        z = unname(z[i, ]),
        p = unname(p[i, ]),
        stringsAsFactors = FALSE
      )
    }))
    out <- list(
      call = object$call_mfp,
      family = object$family_string,
      criterion = mfp2_summary_criterion_label(object),
      converged = isTRUE(object$convergence_mfp),
      n = mfp2_summary_nobs(object),
      nevents = NA_integer_,
      distribution = NULL,
      scale = NULL,
      scale_fixed = NULL,
      dispersion = NULL,
      dispersion_fixed = NULL,
      scale_strata = NULL,
      distribution_parameters = NULL,
      censoring = NULL,
      link = NULL,
      theta = NULL,
      multinomial = TRUE,
      class_levels = object$class_levels,
      reference_class = object$reference_class,
      n_logits = object$n_logits,
      response_frequencies = mfp2_summary_response_frequencies(object),
      function_table = mfp2_summary_classify_terms(object)$function_table,
      coefficients = coefficient_table,
      model_fit_values = mfp2_summary_model_fit_values(object),
      raw_summary = raw_summary,
      notes = notes,
      digits = digits
    )
    class(out) <- "summary.mfp2"
    return(out)
  }

  # ---------------------------------------------------------------------------
  # Gather the pieces.
  # ---------------------------------------------------------------------------
  classified <- mfp2_summary_classify_terms(object)
  # summary.orm() prints instead of constructing a reusable coefficient table.
  # Ordinal inference is built directly from coef() and the complete covariance
  # matrix, avoiding backend output while the structured summary is assembled.
  raw_summary <- if (mfp2_family_is_ordinal(object$family_string)) {
    NULL
  } else {
    mfp2_summary_raw(object)
  }
  model_metadata <- mfp2_summary_model_metadata(
    object,
    raw_summary = raw_summary
  )

  linear_terms <- mfp2_summary_linear_table(object, classified, raw_summary)
  nonlinear_terms <- mfp2_summary_nonlinear_table(object, classified)

  ordinal_intercepts <- if (mfp2_family_is_ordinal(object$family_string)) {
    mfp2_summary_ordinal_intercepts(object)
  } else {
    NULL
  }

  basis_table <- if (isTRUE(basis)) {
    mfp2_summary_basis_table(object, classified)
  } else {
    NULL
  }

  formula_strings <- if (isTRUE(formulas)) {
    mfp2_summary_formula_strings(object, classified)
  } else {
    NULL
  }

  acd_definitions <- if (isTRUE(formulas) || isTRUE(basis)) {
    mfp2_summary_acd_definitions(object, classified)
  } else {
    NULL
  }

  out <- list(
    call            = object$call_mfp,
    family          = object$family_string,
    criterion       = mfp2_summary_criterion_label(object),
    converged       = isTRUE(object$convergence_mfp),
    n               = mfp2_summary_nobs(object),
    nevents         = mfp2_summary_nevents(object),
    distribution    = model_metadata$distribution,
    scale           = model_metadata$scale,
    scale_fixed     = model_metadata$scale_fixed,
    dispersion      = model_metadata$dispersion,
    dispersion_fixed = model_metadata$dispersion_fixed,
    scale_strata    = model_metadata$scale_strata,
    distribution_parameters = model_metadata$distribution_parameters,
    censoring       = model_metadata$censoring,
    link            = model_metadata$link,
    theta           = model_metadata$theta,
    function_table  = classified$function_table,
    linear_terms    = linear_terms,
    nonlinear_terms = nonlinear_terms,
    basis           = basis_table,
    formulas        = formula_strings,
    acd_definitions = acd_definitions,
    fit                = mfp2_summary_fit_stats(object),
    raw_summary        = raw_summary,
    exp_label          = mfp2_summary_exponent_label(object),
    ordinal_intercepts = ordinal_intercepts,
    response_frequencies = mfp2_summary_response_frequencies(object),
    notes              = notes,
    digits             = digits
  )
  class(out) <- "summary.mfp2"
  out
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Number of Fitted Observations for an `mfp2` Object
#'
#' Returns the number of observations used to fit the model, robust to a
#' missing `nobs` field on the fitted object.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return Integer scalar with the number of fitted observations, or
#'   `NA_integer_` if it cannot be determined.
#'
#' @keywords internal
#' @noRd
mfp2_summary_nobs <- function(object) {
  if (!is.null(object$nobs)) return(as.integer(object$nobs))
  n <- tryCatch(NROW(object$y), error = function(e) NA_integer_)
  as.integer(n)
}

#' Unweighted Response Frequencies for Categorical `mfp2` Models
#'
#' Computes the unweighted class-frequency counts displayed for a fitted
#' ordinal or multinomial `"mfp2"` model. The display follows the original
#' response order, not the reference-first reordering used internally by the
#' multinomial fitting engine. For a grouped multinomial count response, the
#' reported frequencies are column totals of the count matrix.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A named numeric vector of class frequencies, in original response
#'   order, or `NULL` when the family is not ordinal or multinomial.
#'
#' @keywords internal
#' @noRd
mfp2_summary_response_frequencies <- function(object) {
  family_string <- object$family_string
  if (!isTRUE(family_string %in% c("ordinal", "multinomial"))) return(NULL)

  response <- object$y_original
  if (is.null(response) && identical(family_string, "ordinal")) {
    response <- object$y
  }
  if (is.null(response)) return(NULL)

  if (identical(family_string, "multinomial")) {
    prepared <- if (is.list(object$mfp2_family)) {
      object$mfp2_family$prepared
    } else {
      NULL
    }
    original_levels <- if (is.list(prepared)) prepared$original_levels else NULL

    if (is.matrix(response)) {
      response_names <- colnames(response)
      if (is.null(response_names)) {
        response_names <- original_levels
      }
      if (is.null(response_names) || length(response_names) != ncol(response)) {
        response_names <- paste0("class", seq_len(ncol(response)))
      }
      totals <- stats::setNames(as.numeric(colSums(response)), response_names)
      if (!is.null(original_levels) && all(original_levels %in% names(totals))) {
        totals <- totals[original_levels]
      }
      return(totals)
    }

    if (is.null(original_levels)) {
      original_levels <- if (is.factor(response)) {
        levels(droplevels(response))
      } else {
        levels(factor(response))
      }
    }
    counts <- table(factor(as.character(response), levels = original_levels))
    return(stats::setNames(as.numeric(counts), original_levels))
  }

  ordinal_levels <- object$mfp2_ordinal_levels
  if (is.null(ordinal_levels)) ordinal_levels <- object$yunique
  if (is.null(ordinal_levels)) return(NULL)
  ordinal_levels <- as.character(ordinal_levels)

  response_values <- as.character(response)
  if (!all(response_values %in% ordinal_levels)) {
    codes <- suppressWarnings(as.integer(response))
    if (!anyNA(codes) && all(codes >= 1L & codes <= length(ordinal_levels))) {
      response_values <- ordinal_levels[codes]
    }
  }
  counts <- table(factor(response_values, levels = ordinal_levels))
  stats::setNames(as.numeric(counts), ordinal_levels)
}

#' Number of Target Events for a Proportional-Hazards `mfp2` Model
#'
#' Returns the number of target events for Cox and Fine--Gray models. Uses
#' the stored `nevents` field when present, otherwise counts `status == 1`
#' rows of a `Surv` response.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return Integer scalar with the number of events, or `NA_integer_` for
#'   families that do not report event counts.
#'
#' @keywords internal
#' @noRd
mfp2_summary_nevents <- function(object) {
  if (!mfp2_family_uses_event_count(object$family_string)) return(NA_integer_)
  if (!is.null(object$nevents)) return(as.integer(object$nevents))
  y <- object$y
  if (inherits(y, "Surv")) {
    status_col <- ncol(y)
    return(as.integer(sum(y[, status_col] == 1, na.rm = TRUE)))
  }
  NA_integer_
}

#' Extract Model-Specific Interpretation Metadata
#'
#' Collects the family-specific metadata required to interpret parametric
#' survival, GLM, negative-binomial, and ordinal fits: distribution name,
#' scale (with strata table when applicable), fixed-vs-estimated flags,
#' dispersion, link function, negative-binomial `theta`, censoring summary,
#' and any extra distribution parameters. Centralising the extraction in one
#' helper guarantees that [print.mfp2()] and [print.summary.mfp2()] report
#' exactly the same values.
#'
#' @param object An `"mfp2"` model object.
#' @param family_string Canonical family name. Defaults to
#'   `object$family_string`.
#' @param raw_summary Optional raw summary produced by `mfp2_summary_raw()`,
#'   used when dispersion has to be pulled from the underlying-class summary.
#'
#' @return A named list carrying `distribution`, `scale`, `scale_fixed`,
#'   `dispersion`, `dispersion_fixed`, `scale_strata`,
#'   `distribution_parameters`, `censoring`, `link`, and `theta`. Fields not
#'   applicable to the family remain `NULL`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_model_metadata <- function(object,
                                        family_string = object$family_string,
                                        raw_summary = NULL) {
  out <- list(
    distribution = NULL,
    scale = NULL,
    scale_fixed = NULL,
    dispersion = NULL,
    dispersion_fixed = NULL,
    scale_strata = NULL,
    distribution_parameters = NULL,
    censoring = NULL,
    link = NULL,
    theta = NULL
  )

  if (identical(family_string, "survreg")) {
    requested_dist <- object$mfp2_survreg_distribution
    if (is.null(requested_dist) && is.list(object$family) &&
        !is.null(object$family$dist)) {
      requested_dist <- object$family$dist
    }
    if (is.null(requested_dist)) requested_dist <- object$dist

    if (is.character(requested_dist) && length(requested_dist) == 1L) {
      distributions <- survival::survreg.distributions
      registered <- distributions[[requested_dist]]
      out$distribution <- if (!is.null(registered$name)) {
        registered$name
      } else {
        requested_dist
      }
    } else if (is.list(requested_dist) &&
               is.character(requested_dist$name) &&
               length(requested_dist$name) == 1L) {
      out$distribution <- requested_dist$name
    }

    out$scale <- if (is.numeric(object$scale)) unname(object$scale) else NULL

    scale_fixed <- object$mfp2_survreg_scale_fixed
    if (is.null(scale_fixed)) {
      family_scale <- if (is.list(object$family)) object$family$scale else NULL
      fixed_by_argument <- is.numeric(family_scale) &&
        length(family_scale) == 1L && is.finite(family_scale) && family_scale > 0
      fixed_by_distribution <- FALSE
      family_dist <- if (is.list(object$family)) object$family$dist else NULL
      if (is.character(family_dist) && length(family_dist) == 1L) {
        registered <- survival::survreg.distributions[[family_dist]]
        fixed_by_distribution <- !is.null(registered$scale)
      } else if (is.list(family_dist)) {
        fixed_by_distribution <- !is.null(family_dist$scale)
      }
      scale_fixed <- fixed_by_argument || fixed_by_distribution
    }
    out$scale_fixed <- isTRUE(scale_fixed)

    if (length(out$scale) > 1L) {
      strata <- object$mfp2_strata_levels
      if (is.null(strata) || length(strata) != length(out$scale)) {
        strata <- names(object$scale)
      }
      if (is.null(strata) || length(strata) != length(out$scale) ||
          anyNA(strata) || any(!nzchar(as.character(strata)))) {
        strata <- as.character(seq_along(out$scale))
      }
      out$scale_strata <- data.frame(
        stratum = as.character(strata),
        scale = out$scale,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    parms <- object$mfp2_survreg_parms
    if (is.null(parms) && !is.null(object$parms)) parms <- object$parms
    if (is.numeric(parms) && length(parms) > 0L) {
      out$distribution_parameters <- parms
    }

    out$censoring <- mfp2_summary_survreg_censoring(object)
  }

  if (identical(family_string, "negbin")) {
    link <- if (is.list(object$family)) object$family$link else NULL
    if (is.function(link)) link <- NULL
    if (is.character(link) && length(link) == 1L && !is.na(link)) {
      out$link <- link
    } else {
      # The supported fastglm negative-binomial backend uses the log link.
      out$link <- "log"
    }
    if (is.numeric(object$theta) && length(object$theta) == 1L &&
        is.finite(object$theta) && object$theta > 0) {
      out$theta <- unname(object$theta)
    }
  }

  if (isTRUE(mfp2_family_is_glm(family_string)) &&
      !identical(family_string, "negbin")) {
    link <- if (is.list(object$family)) object$family$link else NULL
    if (is.character(link) && length(link) == 1L && !is.na(link) &&
        nzchar(link)) {
      out$link <- link
    }
  }

  # Gaussian, Gamma, and inverse-Gaussian GLMs estimate (or may explicitly
  # fix) a dispersion parameter in the underlying fitting function. Preserve
  # that native value and terminology; do not derive an SD, Gamma shape, or
  # any other transformed quantity for the custom mfp2 display.
  if (isTRUE(family_string %in% c(
    "gaussian", "Gamma", "inverse.gaussian"
  ))) {
    family_dispersion <- if (is.list(object$family)) {
      object$family$dispersion
    } else {
      NULL
    }
    fixed_dispersion <- is.numeric(family_dispersion) &&
      length(family_dispersion) == 1L && is.finite(family_dispersion) &&
      family_dispersion > 0

    if (fixed_dispersion) {
      out$dispersion <- unname(family_dispersion)
      out$dispersion_fixed <- TRUE
    } else {
      if (is.null(raw_summary)) raw_summary <- mfp2_summary_raw(object)
      fitted_dispersion <- if (!is.null(raw_summary)) {
        raw_summary$dispersion
      } else {
        NULL
      }
      if (is.numeric(fitted_dispersion) &&
          length(fitted_dispersion) == 1L &&
          is.finite(fitted_dispersion) && fitted_dispersion > 0) {
        out$dispersion <- unname(fitted_dispersion)
        out$dispersion_fixed <- FALSE
      }
    }
  }

  if (mfp2_family_is_ordinal(family_string)) {
    link <- object$mfp2_ordinal_link
    if (is.character(link) && length(link) == 1L && !is.na(link)) {
      out$link <- link
    }
  }

  out
}


#' Human-Readable Model Label for an `mfp2` Family
#'
#' Converts the normalised internal family identifier into a printable model
#' label. Link-dependent models include their fitted link so that generic
#' identifiers such as `"ordinal"` are not presented as complete statistical
#' model specifications.
#'
#' @param family_string Canonical family name.
#' @param metadata Optional metadata list produced by
#'   `mfp2_summary_model_metadata()`. When supplied, `metadata$link` is
#'   appended to the label where the family admits multiple links.
#'
#' @return Character scalar with a human-readable model label.
#'
#' @keywords internal
#' @noRd
mfp2_summary_model_label <- function(family_string, metadata = NULL) {
  if (!is.character(family_string) || length(family_string) != 1L ||
      is.na(family_string) || !nzchar(family_string)) {
    return("Unknown Model")
  }
  link <- if (is.list(metadata)) metadata$link else NULL
  if (!is.character(link) || length(link) != 1L || is.na(link) || !nzchar(link)) {
    link <- NULL
  }

  if (mfp2_family_is_ordinal(family_string)) {
    if (is.null(link)) return("Ordinal Regression")
    return(switch(
      link,
      logistic = "Logistic Ordinal Regression (proportional odds)",
      probit = "Probit Ordinal Regression",
      loglog = "Log-log Ordinal Regression",
      cloglog = "Complementary Log-log Ordinal Regression",
      cauchit = "Cauchit Ordinal Regression",
      "Ordinal Regression"
    ))
  }

  if (identical(family_string, "multinomial")) {
    return("Multinomial Logistic Regression")
  }

  if (identical(family_string, "negbin")) {
    label <- "Negative-Binomial GLM"
    if (!is.null(link)) label <- paste0(label, " (", link, " link)")
    return(label)
  }

  glm_names <- c(
    gaussian = "Gaussian",
    binomial = "Binomial",
    poisson = "Poisson",
    Gamma = "Gamma",
    inverse.gaussian = "Inverse-Gaussian"
  )
  if (isTRUE(family_string %in% names(glm_names))) {
    label <- paste0(unname(glm_names[[family_string]]), " GLM")
    if (!is.null(link)) label <- paste0(label, " (", link, " link)")
    return(label)
  }

  switch(
    family_string,
    cox = "Cox Proportional Hazards",
    finegray = "Fine--Gray Subdistribution Hazards",
    survreg = "Parametric Survival Regression",
    family_string
  )
}


#' Print the Response-Frequency Block for Categorical Outcomes
#'
#' Prints the native-style response-frequency block used for categorical
#' multinomial and ordinal outcomes. Counts are deliberately unweighted so
#' the printed frequencies match the raw response data.
#'
#' @param response_frequencies A named numeric vector of class frequencies,
#'   typically produced by `mfp2_summary_response_frequencies()`, or `NULL`.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of printing.
#'
#' @keywords internal
#' @noRd
mfp2_print_response_frequencies <- function(response_frequencies) {
  if (is.null(response_frequencies) || length(response_frequencies) == 0L) {
    return(invisible(NULL))
  }
  cat("\nFrequencies of Responses\n\n")
  print(response_frequencies)
  invisible(NULL)
}


#' Count Exact and Censored Observations for a `survreg` Model
#'
#' Counts exact events and censoring events for a fitted parametric survival
#' model according to the response's `Surv` type. Status value `0` has
#' different meanings for right- and left-censored responses, so a generic
#' `status == 1` event count is not sufficient.
#'
#' @param object An `"mfp2"` model object with a `Surv` response.
#'
#' @return A named list with `type`, `exact`, `right_censored`,
#'   `left_censored`, and `interval_censored`, or `NULL` if the response
#'   type is not supported by the reporting block.
#'
#' @keywords internal
#' @noRd
mfp2_summary_survreg_censoring <- function(object) {
  y <- object$y
  if (!inherits(y, "Surv")) return(NULL)

  type <- attr(y, "type", exact = TRUE)
  status <- y[, NCOL(y)]
  result <- list(
    type = type,
    exact = as.integer(sum(status == 1, na.rm = TRUE)),
    right_censored = 0L,
    left_censored = 0L,
    interval_censored = 0L
  )

  if (identical(type, "right")) {
    result$right_censored <- as.integer(sum(status == 0, na.rm = TRUE))
  } else if (identical(type, "left")) {
    result$left_censored <- as.integer(sum(status == 0, na.rm = TRUE))
  } else if (type %in% c("interval", "interval2")) {
    result$right_censored <- as.integer(sum(status == 0, na.rm = TRUE))
    result$left_censored <- as.integer(sum(status == 2, na.rm = TRUE))
    result$interval_censored <- as.integer(sum(status == 3, na.rm = TRUE))
  } else {
    return(NULL)
  }

  result
}


#' Print the Shared Model-Header Block
#'
#' Prints the common top-level metadata used by both the `mfp2` and
#' `summary.mfp2` print methods: model label, criterion, convergence status,
#' observation count with censoring/event breakdown, family-specific
#' parameters, and, for categorical outcomes, response frequencies.
#'
#' @param family_string Canonical family name.
#' @param criterion Character scalar identifying the selection criterion.
#' @param converged Logical indicating whether MFP selection converged, or
#'   `NA` when unknown.
#' @param n Integer number of fitted observations.
#' @param nevents Integer number of events for proportional-hazards models,
#'   or `NA_integer_`.
#' @param metadata Metadata list from `mfp2_summary_model_metadata()`.
#' @param digits Integer scalar controlling the number of digits used for
#'   numeric parameters.
#' @param response_frequencies Optional named numeric vector of response
#'   frequencies from `mfp2_summary_response_frequencies()`.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of printing.
#'
#' @keywords internal
#' @noRd
mfp2_print_model_header <- function(family_string, criterion, converged, n,
                                    nevents = NA_integer_, metadata = NULL,
                                    digits = 3L,
                                    response_frequencies = NULL) {
  model_label <- mfp2_summary_model_label(family_string, metadata)
  fields <- c(paste0("Model: ", model_label))
  if (identical(family_string, "survreg") &&
      !is.null(metadata$distribution)) {
    fields <- c(fields, paste0("Distribution: ", metadata$distribution))
  }
  fields <- c(
    fields,
    paste0("Criterion: ", criterion),
    paste0(
      "Converged: ",
      if (isTRUE(converged)) "yes" else if (identical(converged, FALSE)) "no" else "unknown"
    )
  )
  cat(paste(fields, collapse = " | "), "\n", sep = "")

  censoring <- metadata$censoring
  if (identical(family_string, "survreg") && !is.null(censoring)) {
    if (identical(censoring$type, "right")) {
      cat(sprintf(
        "Observations: %s | Events: %s | Censored: %s\n",
        n, censoring$exact, censoring$right_censored
      ))
    } else if (identical(censoring$type, "left")) {
      cat(sprintf(
        "Observations: %s | Exact: %s | Left-censored: %s\n",
        n, censoring$exact, censoring$left_censored
      ))
    } else {
      cat(sprintf(
        paste0(
          "Observations: %s | Exact: %s | Right-censored: %s | ",
          "Left-censored: %s | Interval-censored: %s\n"
        ),
        n, censoring$exact, censoring$right_censored,
        censoring$left_censored, censoring$interval_censored
      ))
    }
  } else if (!is.na(nevents)) {
    cat(sprintf("Observations: %s | Events: %s\n", n, nevents))
  } else {
    cat(sprintf("Observations: %s\n", n))
  }

  mfp2_print_model_specific_parameters(
    family_string = family_string,
    metadata = metadata,
    digits = digits
  )
  mfp2_print_response_frequencies(response_frequencies)
  invisible(NULL)
}


#' Print Family-Specific Nuisance and Model Parameters
#'
#' Prints the nuisance and model parameters (scale, distribution parameters,
#' negative-binomial `theta`, GLM dispersion) that are needed to interpret a
#' fit but do not belong in the regression coefficient table. Parameter
#' formatting follows the same fixed-decimal convention used for the
#' coefficient and fit-statistic tables.
#'
#' @param family_string Canonical family name.
#' @param metadata Metadata list from `mfp2_summary_model_metadata()`.
#' @param digits Integer scalar controlling the number of digits used when
#'   formatting numeric parameters.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of printing.
#'
#' @keywords internal
#' @noRd
mfp2_print_model_specific_parameters <- function(family_string, metadata,
                                                  digits = 3L) {
  # Model-specific numeric parameters follow the same fixed-decimal convention
  # as the coefficient and fit-statistic tables.
  fmt <- function(value) {
    format_print_decimal(value, digits)
  }

  if (identical(family_string, "survreg") && length(metadata$scale) > 0L) {
    status <- if (isTRUE(metadata$scale_fixed)) "fixed" else "estimated"
    if (length(metadata$scale) == 1L) {
      cat(sprintf("Scale: %s (%s)\n", fmt(metadata$scale), status))
    } else {
      cat(sprintf("Scale parameters (%s):\n", status))
      scale_table <- metadata$scale_strata
      names(scale_table) <- c("Stratum", "Scale")
      scale_table <- format_model_print_table(scale_table, digits)
      print.data.frame(
        scale_table,
        row.names = FALSE,
        right = FALSE
      )
    }

    parms <- metadata$distribution_parameters
    if (length(parms) > 0L) {
      parm_names <- names(parms)
      if (is.null(parm_names) || any(!nzchar(parm_names))) {
        parm_names <- paste0("parameter", seq_along(parms))
      }
      values <- paste0(parm_names, "=", vapply(parms, fmt, character(1L)))
      cat("Distribution parameters: ", paste(values, collapse = ", "), "\n", sep = "")
    }
  }

  if (identical(family_string, "negbin") && length(metadata$theta) == 1L) {
    cat(sprintf("Theta: %s (estimated)\n", fmt(metadata$theta)))
  }

  if (isTRUE(family_string %in% c(
        "gaussian", "Gamma", "inverse.gaussian"
      )) &&
      length(metadata$dispersion) == 1L) {
    status <- if (isTRUE(metadata$dispersion_fixed)) "fixed" else "estimated"
    cat(sprintf("Dispersion: %s (%s)\n", fmt(metadata$dispersion), status))
  }

  invisible(NULL)
}

#' Human-Readable Selection-Criterion Label
#'
#' Returns a printable label for the MFP selection criterion recorded on a
#' fitted model. Falls back to `"p-value"` when the criterion field is
#' absent, matching the `fp_terms` select/alpha convention.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return Character scalar such as `"p-value"`, `"AIC"`, or `"BIC"`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_criterion_label <- function(object) {
  crit <- object$criterion_mfp
  if (is.null(crit) || length(crit) != 1L || is.na(crit)) {
    # Fall back to fp_terms select/alpha convention.
    return("p-value")
  }
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(crit)))
  switch(key,
         "pvalue" = "p-value",
         "aic"    = "AIC",
         "bic"    = "BIC",
         as.character(crit)
  )
}

#' Underlying-Class Summary With `mfp2()` Call Restored
#'
#' Computes the underlying model-specific summary (for example
#' [stats::summary.glm()] or [survival::summary.coxph()]) with the fitted
#' `call` field replaced by the original [mfp2()] call. Storing this on the
#' summary result lets users access the underlying-class summary without
#' having to refit the model.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return The underlying-class summary object, or `NULL` when the
#'   dispatched summary throws an error.
#'
#' @keywords internal
#' @noRd
mfp2_summary_raw <- function(object) {
  result <- tryCatch(
    NextMethod_summary(object),
    error = function(e) NULL
  )
  if (!is.null(result) && !is.null(object$call_mfp)) {
    result$call <- object$call_mfp
  }
  result
}

#' Manual `NextMethod()` Dispatcher for `summary()`
#'
#' `NextMethod()` only works inside a method definition; this helper
#' dispatches manually to the next class after `"mfp2"` so that
#' `mfp2_summary_raw()` can build the underlying-class summary from an
#' ordinary function. For negative-binomial models fitted through `fastglm`,
#' the GLM dispersion is forced to `1` so that `summary()` uses normal/z
#' inference, matching [MASS::summary.negbin()], rather than estimating a
#' second dispersion and reporting t statistics.
#'
#' @param object An `"mfp2"` model object.
#' @param ... Additional arguments forwarded to the dispatched summary
#'   method.
#'
#' @return The result of the underlying-class `summary()` call, or `NULL`
#'   if no method exists for the remaining classes.
#'
#' @keywords internal
#' @noRd
NextMethod_summary <- function(object, ...) {
  cls <- class(object)
  next_classes <- cls[which(cls == "mfp2")[1L] + 1L]
  next_classes <- next_classes[!is.na(next_classes)]
  if (length(next_classes) == 0L) {
    return(NULL)
  }
  # Temporarily strip "mfp2" so summary() dispatches to the underlying class.
  obj2 <- object
  class(obj2) <- cls[cls != "mfp2"]

  # A negative-binomial GLM has fixed unit GLM dispersion after theta
  # has been estimated. Force dispersion = 1 so summary.fastglm() uses
  # normal/z inference, matching MASS::summary.negbin(), rather than
  # estimating a second dispersion and reporting t statistics.
  dots <- list(...)
  if (identical(object$family_string, "negbin") &&
      inherits(obj2, "fastglm")) {
    dots$dispersion <- 1
  }
  do.call(summary, c(list(object = obj2), dots))
}

# ---------------------------------------------------------------------------
# Term classification
# ---------------------------------------------------------------------------

#' Classify Selected Variables Into Linear and Nonlinear Terms
#'
#' Splits the selected variables of an `"mfp2"` model into linear and
#' nonlinear terms, records which are spike-at-zero binary-only outcomes,
#' and maps each variable to its transformed coefficient columns by name
#' prefix. Also returns the selection overview table used in the printed
#' header.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A named list with the classification vectors (`selected`,
#'   `is_linear`, `is_nonlinear`, `binary_only`, `acd`, `zero`, `catzero`,
#'   `spike`, `spike_dec`), per-variable power vectors preserving positions
#'   (`power_slots_by_var`) and with `NA` removed (`powers_by_var`), the
#'   design-column-to-variable mapping (`cols_by_var`), final degrees of
#'   freedom (`df_final`), and the printable `function_table` overview.
#'
#' @keywords internal
#' @noRd
mfp2_summary_classify_terms <- function(object) {
  fp_terms <- object$fp_terms
  variable_names <- rownames(fp_terms)

  selected <- if ("selected" %in% names(fp_terms)) {
    as.logical(fp_terms[["selected"]])
  } else {
    rep(TRUE, nrow(fp_terms))
  }
  selected[is.na(selected)] <- FALSE

  acd     <- mfp2_summary_flag(fp_terms, "acd")
  zero    <- mfp2_summary_flag(fp_terms, "zero")
  catzero <- mfp2_summary_flag(fp_terms, "catzero")
  spike   <- mfp2_summary_flag(fp_terms, "spike")

  df_final <- if ("df_final" %in% names(fp_terms)) {
    as.numeric(fp_terms[["df_final"]])
  } else {
    rep(NA_real_, nrow(fp_terms))
  }

  power_cols <- grep("^power[0-9]+$", names(fp_terms), value = TRUE)

  # Keep both representations of the selected powers:
  #
  # * power_slots_by_var preserves the original positions, including NA. This
  #   is required for ACD terms because slot 1 applies to x and slot 2 applies
  #   to A(x); c(NA, 1) is therefore different from c(1, NA).
  # * powers_by_var removes NA values for the existing selection summaries and
  #   ordinary FP classification.
  power_slots_by_var <- lapply(seq_len(nrow(fp_terms)), function(i) {
    suppressWarnings(as.numeric(unlist(
      fp_terms[i, power_cols, drop = FALSE],
      use.names = FALSE
    )))
  })
  names(power_slots_by_var) <- variable_names

  powers_by_var <- lapply(power_slots_by_var, function(p) p[!is.na(p)])

  # Map variable -> transformed coefficient column names, by stripping the
  # ".<index>" suffix from the fitted design column names. ACD component
  # columns are named A_<variable>.<index>; map those columns back to their
  # source variable so direct and ACD components are formatted together.
  design_cols <- colnames(object$x)
  if (is.null(design_cols) && !is.null(object$mfp2_design)) {
    design_cols <- setdiff(colnames(object$mfp2_design), "(Intercept)")
  }
  if (is.null(design_cols)) design_cols <- names(object$coefficients)
  base_names <- sub("\\.[0-9]+$", "", design_cols)
  source_names <- base_names

  acd_variables <- variable_names[acd]
  for (v in acd_variables) {
    source_names[base_names == paste0("A_", v)] <- v
  }

  cols_by_var <- split(design_cols, source_names)

  # SAZ decision code -> is this a "binary only" outcome?
  spike_dec <- if ("spike_dec" %in% names(fp_terms)) {
    suppressWarnings(as.integer(fp_terms[["spike_dec"]]))
  } else {
    rep(NA_integer_, nrow(fp_terms))
  }

  # Classify each selected variable.
  is_linear <- logical(length(variable_names))
  is_nonlinear <- logical(length(variable_names))
  binary_only <- logical(length(variable_names))

  for (i in seq_along(variable_names)) {
    if (!selected[i]) next

    p <- powers_by_var[[i]]
    plain_linear <- (length(p) == 1L && isTRUE(p == 1)) &&
      !acd[i] && !zero[i] && !catzero[i] && !spike[i]

    # Binary-only spike: no continuous component survived.
    saz_binary_only <- spike[i] &&
      mfp2_summary_saz_is_binary_only(spike_dec[i])

    if (plain_linear) {
      is_linear[i] <- TRUE
    } else if (saz_binary_only) {
      is_linear[i] <- TRUE
      binary_only[i] <- TRUE
    } else {
      is_nonlinear[i] <- TRUE
    }
  }

  function_table <- mfp2_summary_function_overview(
    variable_names, selected, df_final, powers_by_var,
    acd, zero, catzero, spike, spike_dec
  )

  list(
    variable_names = variable_names,
    selected       = selected,
    is_linear      = is_linear,
    is_nonlinear   = is_nonlinear,
    binary_only    = binary_only,
    acd            = acd,
    zero           = zero,
    catzero        = catzero,
    spike          = spike,
    spike_dec      = spike_dec,
    df_final          = df_final,
    powers_by_var     = powers_by_var,
    power_slots_by_var = power_slots_by_var,
    cols_by_var       = cols_by_var,
    function_table = function_table
  )
}

#' Coerce an `fp_terms` Flag Column to a Logical Vector
#'
#' Reads one boolean-valued column of the `fp_terms` data frame and coerces
#' it to a proper logical vector, tolerating logical, character, and
#' numeric encodings. Missing values are treated as `FALSE` so downstream
#' classification never has to test for `NA`.
#'
#' @param fp_terms A data frame with one row per model term (see the
#'   `fp_terms` component of an `"mfp2"` object).
#' @param name Character scalar naming the column to coerce.
#'
#' @return Logical vector of length `nrow(fp_terms)`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_flag <- function(fp_terms, name) {
  if (!name %in% names(fp_terms)) return(rep(FALSE, nrow(fp_terms)))
  v <- fp_terms[[name]]
  if (is.logical(v)) return(ifelse(is.na(v), FALSE, v))
  tolower(as.character(v)) %in% c("true", "t", "yes", "y", "1")
}

#' Is a Spike-at-Zero Decision Code the Binary-Only Outcome?
#'
#' Package spike-at-zero decisions use integer codes; `3L` means that the
#' continuous positive-component function was dropped and only the binary
#' zero indicator was retained. Matches `saz_decision_label()`.
#'
#' @param code Integer decision code, possibly `NA`.
#'
#' @return `TRUE` when `code == 3L`, `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_summary_saz_is_binary_only <- function(code) {
  if (is.na(code)) return(FALSE)
  # In mfp2, spike_dec == 3 corresponds to "binary only".
  identical(as.integer(code), 3L)
}

#' Human-Readable Functional-Form Label for One Variable
#'
#' Builds the printable functional-form label for one MFP term, combining
#' the selected FP powers with ACD, `zero`, `catzero`, and spike-at-zero
#' modifiers. Returns `"out"` when the variable was not selected and
#' `"binary indicator only"` when a `catzero`- or `spike`-handled term
#' retains no continuous component.
#'
#' @param powers Numeric vector of selected FP powers, possibly empty.
#' @param acd Logical: was the term fitted with the ACD extension?
#' @param zero Logical: does the term use exact-zero handling?
#' @param catzero Logical: does the term use `catzero` handling?
#' @param spike Logical: was the term assessed with spike-at-zero?
#' @param spike_dec Integer spike-at-zero decision code (see
#'   `mfp2_summary_saz_is_binary_only()`).
#' @param selected Logical: was the term retained in the final model?
#'
#' @return Character scalar with the functional-form label.
#'
#' @keywords internal
#' @noRd
mfp2_summary_form_label <- function(powers, acd, zero, catzero, spike,
                                    spike_dec, selected) {
  if (!isTRUE(selected)) return("out")

  has_cont <- length(powers) > 0L

  if (!has_cont) {
    if (catzero || spike) return("binary indicator only")
    return("out")
  }

  base <- if (length(powers) == 1L && isTRUE(powers == 1)) {
    "linear"
  } else {
    sprintf("FP(%s)", paste(powers, collapse = ", "))
  }
  if (acd) base <- paste0("ACD ", base)
  if (zero) base <- paste0(base, " (x > 0)")
  if (catzero) base <- paste0(base, " + binary")
  base
}

#' Build the Selection-Overview Table
#'
#' Assembles the printable variable-selection overview shown in the header
#' of `print.summary.mfp2()`. Each row reports the variable name, whether
#' it was retained, its final degrees of freedom, and its functional-form
#' label from `mfp2_summary_form_label()`.
#'
#' @param variable_names Character vector of term names.
#' @param selected Logical vector: which terms were retained.
#' @param df_final Numeric vector of final degrees of freedom per term.
#' @param powers_by_var Named list of selected FP power vectors.
#' @param acd Logical vector: ACD flags.
#' @param zero Logical vector: `zero` flags.
#' @param catzero Logical vector: `catzero` flags.
#' @param spike Logical vector: spike-at-zero flags.
#' @param spike_dec Integer vector of spike-at-zero decision codes.
#'
#' @return A data frame with columns `Variable`, `Selected`, `df`, and
#'   `Function`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_function_overview <- function(variable_names, selected, df_final,
                                           powers_by_var, acd, zero, catzero,
                                           spike, spike_dec) {
  form <- vapply(seq_along(variable_names), function(i) {
    mfp2_summary_form_label(
      powers_by_var[[i]], acd[i], zero[i], catzero[i],
      spike[i], spike_dec[i], selected[i]
    )
  }, character(1L))

  df <- ifelse(is.na(df_final), ".", format(df_final, trim = TRUE))

  data.frame(
    Variable = variable_names,
    Selected = ifelse(selected, "yes", "no"),
    df       = df,
    Function = form,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Coefficient -> display-name lookup
# ---------------------------------------------------------------------------

#' Build a Fitted-Coefficient to Display-Name Lookup
#'
#' Builds a named character vector mapping every fitted coefficient name to
#' the user-facing display name shown in the printed summary. The mapping is
#' composed from the two authoritative lookups populated by `fit_mfp()` at
#' fit time:
#'
#' - `term_to_columns`: user variable name to raw source columns.
#' - `transformed_to_model_columns`: raw source column to fitted coefficient
#'   name.
#'
#' Composing them yields an exact `coef_name -> display_name` map for every
#' selected variable, so the fitter's suffix convention is treated as an
#' implementation detail rather than something the summary depends on. The
#' caller falls back to a regex strip when either lookup is missing.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A named character vector: names are fitted coefficient names,
#'   values are display names. Empty when either authoritative lookup is
#'   missing.
#'
#' @keywords internal
#' @noRd
mfp2_summary_coef_to_display <- function(object) {
  t2c <- object$term_to_columns
  tmc <- object$transformed_to_model_columns

  if (is.null(t2c) || is.null(tmc)) {
    return(stats::setNames(character(0L), character(0L)))
  }

  out_names  <- character(0L)
  out_values <- character(0L)

  for (v in names(t2c)) {
    raw_cols <- as.character(t2c[[v]])
    if (length(raw_cols) == 0L) next

    # Rule: single raw column whose name equals v -> display is v.
    # Otherwise display is the raw column name (preserves factor levels
    # and any renaming the fitter did between user and raw column).
    display_for_raw <- if (length(raw_cols) == 1L && identical(raw_cols[[1L]], v)) {
      stats::setNames(v, raw_cols)
    } else {
      stats::setNames(raw_cols, raw_cols)
    }

    for (raw_col in raw_cols) {
      # transformed_to_model_columns is keyed by raw column; its value is the
      # final fitted coefficient name. Names and values can be identical for
      # scalar predictors, but for FP terms with multiple basis columns the
      # same raw column may map to multiple fitted names.
      hits <- which(names(tmc) == raw_col)
      if (length(hits) == 0L) next
      fitted_names <- as.character(tmc[hits])
      out_names  <- c(out_names, fitted_names)
      out_values <- c(out_values, rep(display_for_raw[[raw_col]], length(fitted_names)))
    }
  }

  stats::setNames(out_values, out_names)
}


# ---------------------------------------------------------------------------
# Linear terms table
# ---------------------------------------------------------------------------

#' Complete Ordinal Covariance Matrix in Coefficient Order
#'
#' Returns the full ordinal covariance matrix (thresholds and slopes) for a
#' fitted MFP ordinal model, ordered to match `names(object$coefficients)`.
#' New MFP ordinal fits always retain the native `orm` information matrix,
#' including for intercept-only models, so missing or misaligned covariance
#' is treated as an internal error rather than silently masked with
#' fallback `NA` values.
#'
#' @param object An `"mfp2"` ordinal model object.
#'
#' @return A numeric covariance matrix with row and column names equal to
#'   `names(object$coefficients)`.
#'
#' @keywords internal
#' @noRd
mfp2_ordinal_vcov_all <- function(object) {
  coefficient_names <- names(object$coefficients)
  if (is.null(coefficient_names) || anyNA(coefficient_names) ||
      any(!nzchar(coefficient_names))) {
    stop("Ordinal model coefficients must have non-empty names.", call. = FALSE)
  }

  covariance <- stats::vcov(object, intercepts = "all")
  if (!is.matrix(covariance) || is.null(rownames(covariance)) ||
      is.null(colnames(covariance)) ||
      !all(coefficient_names %in% rownames(covariance)) ||
      !all(coefficient_names %in% colnames(covariance))) {
    stop(
      "Ordinal covariance matrix is not aligned with model coefficients.",
      call. = FALSE
    )
  }

  covariance[coefficient_names, coefficient_names, drop = FALSE]
}


#' Ordinal Slope Covariance Matrix in Fitted-Coefficient Order
#'
#' Returns only the slope-slope block of the ordinal covariance matrix, in
#' fitted-coefficient order. [rms::vcov.orm()] has separate code paths for
#' threshold intercepts and regression slopes; requesting
#' `intercepts = "none"` avoids making slope inference depend on inversion
#' and retention of the full threshold covariance block.
#'
#' @param object An `"mfp2"` ordinal model object.
#'
#' @return A numeric covariance matrix with row and column names equal to
#'   the slope coefficient names. Returns an empty matrix when the model has
#'   no slope coefficients (intercept-only fit).
#'
#' @keywords internal
#' @noRd
mfp2_ordinal_vcov_slopes <- function(object) {
  coefficient_names <- names(object$coefficients)
  if (is.null(coefficient_names) || anyNA(coefficient_names) ||
      any(!nzchar(coefficient_names))) {
    stop("Ordinal model coefficients must have non-empty names.", call. = FALSE)
  }
  n_intercepts <- length(object$mfp2_ordinal_intercepts)
  slope_names <- if (length(coefficient_names) > n_intercepts) {
    coefficient_names[seq.int(n_intercepts + 1L, length(coefficient_names))]
  } else {
    character(0L)
  }
  if (length(slope_names) == 0L) {
    return(matrix(numeric(0L), nrow = 0L, ncol = 0L,
                  dimnames = list(character(0L), character(0L))))
  }

  covariance <- stats::vcov(object, intercepts = "none")
  if (!is.matrix(covariance) || is.null(rownames(covariance)) ||
      is.null(colnames(covariance)) ||
      !all(slope_names %in% rownames(covariance)) ||
      !all(slope_names %in% colnames(covariance))) {
    stop(
      "Ordinal slope covariance matrix is not aligned with model coefficients.",
      call. = FALSE
    )
  }

  covariance[slope_names, slope_names, drop = FALSE]
}


#' Ordinal Threshold-Intercept Summary Data Frame
#'
#' Builds a data frame of ordinal threshold intercepts with estimate,
#' standard error, z statistic, two-sided p value, and 95 % Wald confidence
#' interval, for display in [print.summary.mfp2()]. All thresholds are
#' returned, not just the middle one that [rms::vcov.orm()] defaults to.
#'
#' @param object An `"mfp2"` ordinal model object.
#'
#' @return A data frame with columns `threshold`, `estimate`, `se`, `z`,
#'   `p`, `ci_lower`, and `ci_upper`, or `NULL` when no thresholds are
#'   present.
#'
#' @keywords internal
#' @noRd
mfp2_summary_ordinal_intercepts <- function(object) {
  ord_int <- object$mfp2_ordinal_intercepts
  if (is.null(ord_int) || length(ord_int) == 0L) {
    return(NULL)
  }

  if (is.null(names(ord_int)) || anyNA(names(ord_int)) ||
      any(!nzchar(names(ord_int)))) {
    stop("Ordinal threshold intercepts must have non-empty names.", call. = FALSE)
  }

  # rms::orm's vcov() defaults to intercepts = "mid" (middle intercept only).
  # Request all intercepts explicitly so every threshold gets an SE.
  # The mfp2 object IS the orm object (class c("mfp2", "orm")), so vcov()
  # dispatches directly on it, not on a $fit slot.
  covariance_full <- mfp2_ordinal_vcov_all(object)
  variances <- diag(
    covariance_full[names(ord_int), names(ord_int), drop = FALSE]
  )
  se <- sqrt(variances)

  z    <- ord_int / se
  pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

  data.frame(
    threshold = names(ord_int),
    estimate  = unname(ord_int),
    se        = unname(se),
    z         = unname(z),
    p         = unname(pval),
    ci_lower  = unname(ord_int - 1.96 * se),
    ci_upper  = unname(ord_int + 1.96 * se),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build the Linear-Terms Table for `print.summary.mfp2()`
#'
#' Assembles the data frame of linear-term rows (plus spike-at-zero
#' binary-only rows) shown in the printed summary. Coefficients, standard
#' errors, test statistics, and p values are pulled from the raw summary
#' where possible; when the raw summary does not carry a usable coefficient
#' table the values are recomputed from `stats::vcov(object)`. Display
#' names come from `mfp2_summary_coef_to_display()` rather than a regex
#' strip so that factor levels are preserved and future fitter renaming
#' does not break the display.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#' @param raw_summary Raw underlying-class summary from
#'   `mfp2_summary_raw()`.
#'
#' @return A data frame with one row per linear-term coefficient and columns
#'   `term`, `variable`, `coef`, `se`, `statistic`, `p`, plus `exp_coef`,
#'   `ci_lower`, and `ci_upper` (on the exponentiated scale for
#'   proportional-hazards and comparable models). Returns an empty data
#'   frame when there are no linear terms.
#'
#' @keywords internal
#' @noRd
mfp2_summary_linear_table <- function(object, classified, raw_summary) {
  coefs <- object$coefficients
  if (is.null(coefs) || length(coefs) == 0L) {
    return(data.frame())
  }

  # Column names of the fitted design.
  design_cols <- names(coefs)

  # Which fitted columns belong to linear-classified variables?
  linear_vars <- classified$variable_names[classified$is_linear]

  keep_cols <- unlist(
    classified$cols_by_var[linear_vars],
    use.names = FALSE
  )
  keep_cols <- intersect(design_cols, keep_cols)
  if (length(keep_cols) == 0L) {
    return(data.frame())
  }

  # Extract coef / se / stat / p from the raw summary coefficient matrix.
  cmat <- mfp2_summary_coef_matrix(object, raw_summary)
  rows <- match(keep_cols, rownames(cmat))
  valid <- !is.na(rows)
  keep_cols <- keep_cols[valid]
  rows <- rows[valid]

  est <- cmat[rows, "estimate"]
  se  <- cmat[rows, "se"]
  stat <- cmat[rows, "statistic"]
  pval <- cmat[rows, "pvalue"]

  exponentiate <- mfp2_summary_exponentiates_coefficients(object)

  # -------------------------------------------------------------------------
  # Build a display name for each fitted coefficient.
  #
  # Rather than syntactically stripping ".N" suffixes off coefficient names,
  # we walk the two authoritative lookups populated by fit_mfp() at fit time:
  #
  #   term_to_columns             user variable name  -> raw source columns
  #   transformed_to_model_columns raw source column  -> fitted coef name
  #
  # Composing them gives an exact `fitted_coef_name -> display_name` map,
  # independent of any naming convention the internal fitter happens to use.
  # If future changes rename or renumber transformed columns, the summary
  # tracks those changes automatically because it consults the map, not the
  # regex.
  #
  # Display-name rule per user variable v:
  #
  #   - If v maps to a single raw column whose name equals v exactly, the
  #     display for every fitted coefficient descending from v is v itself.
  #     Covers plain scalar predictors, e.g. "hx" -> "hx.1" -> "hx".
  #
  #   - Otherwise (multi-level factor, or a raw column whose name differs
  #     from the user's variable name, e.g. "ekg_fb"), the display is the
  #     RAW source column name. This preserves the factor level in the
  #     printed table (e.g. "ekg_fb" rather than a bare "ekg_f" that would
  #     appear on multiple rows).
  #
  # The fallback for older or unusual fit objects that lack either lookup is
  # the previous regex strip; this preserves user-visible behaviour rather
  # than emitting cryptic internal names.
  coef_to_display <- mfp2_summary_coef_to_display(object)

  display_variable <- vapply(
    keep_cols,
    function(nm) {
      if (nm %in% names(coef_to_display)) {
        return(coef_to_display[[nm]])
      }
      # Fallback to the regex-based strip only when the lookup does not
      # cover this coefficient (e.g. a partial or corrupted fitted object).
      sub("\\.[0-9]+$", "", nm)
    },
    character(1L),
    USE.NAMES = FALSE
  )

  df <- data.frame(
    term = keep_cols,
    variable = display_variable,
    coef = est,
    se   = se,
    statistic = stat,
    p = pval,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (exponentiate) {
    df$exp_coef <- exp(est)
    df$ci_lower <- exp(est - 1.96 * se)
    df$ci_upper <- exp(est + 1.96 * se)
  } else {
    df$ci_lower <- est - 1.96 * se
    df$ci_upper <- est + 1.96 * se
  }

  attr(df, "statistic_label") <- attr(cmat, "statistic_label", exact = TRUE)

  df
}


#' Does the Summary Report Exponentiated Coefficients for This Family?
#'
#' Reports whether the printed summary should append an exponentiated
#' coefficient column (hazard ratio, subdistribution hazard ratio, time
#' ratio, or generic multiplicative effect). Exponentiation is retained for
#' Cox and Fine--Gray, and for `survreg` distributions that use log time.
#' GLM coefficients remain on their fitted link scale, matching the
#' parameter estimates returned by their fitting functions.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return `TRUE` when the summary should display exponentiated
#'   coefficients, `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_summary_exponentiates_coefficients <- function(object) {
  if (mfp2_family_is_ph(object$family_string)) return(TRUE)

  if (identical(object$family_string, "survreg")) {
    dist <- object$family$dist
    if (!is.character(dist) || length(dist) != 1L) return(FALSE)
    return(tolower(dist) %in% c(
      "weibull", "exponential", "rayleigh", "lognormal", "log normal",
      "loggaussian", "log gaussian", "loglogistic", "log logistic"
    ))
  }

  # GLM summaries report the coefficients fitted by glm()/fastglm() on their
  # native link scale. Do not add derived odds, rate, risk, or mean ratios.
  if (mfp2_family_is_glm(object$family_string)) return(FALSE)

  family <- object$family
  is.list(family) && is.character(family$link) && length(family$link) == 1L &&
    family$link %in% c("log", "logit")
}


#' Label for the Exponentiated-Coefficient Column
#'
#' Returns the printable column label used for the exponentiated-coefficient
#' column in the printed summary. Distinguishes hazard ratio (Cox),
#' subdistribution hazard ratio (Fine--Gray), time ratio (`survreg`), and
#' a generic multiplicative-effect label for other exponentiated families.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return Character scalar with the column label, or `NULL` when the
#'   family does not exponentiate.
#'
#' @keywords internal
#' @noRd
mfp2_summary_exponent_label <- function(object) {
  if (!mfp2_summary_exponentiates_coefficients(object)) return(NULL)

  if (identical(object$family_string, "cox")) return("hazard ratio")
  if (identical(object$family_string, "finegray")) {
    return("subdistribution hazard ratio")
  }
  if (identical(object$family_string, "survreg")) return("time ratio")

  "multiplicative effect"
}

#' Match One Column of a Coefficient Table by Exact Name
#'
#' Returns the position of the first candidate column label found in
#' `column_names`. Candidate order defines precedence: this matters for
#' robust Cox summaries, which can contain both `"robust se"` and
#' `"se(coef)"`. A duplicated candidate is ambiguous and deliberately
#' returns `NA` so the caller can fall back to computing the column from
#' first principles rather than pick an unspecified one.
#'
#' @param column_names Character vector of column names from a raw
#'   coefficient table.
#' @param candidates Character vector of candidate exact-match labels, in
#'   precedence order.
#'
#' @return Integer scalar column position, or `NA_integer_` when no
#'   unambiguous match exists.
#'
#' @keywords internal
#' @noRd
mfp2_summary_match_coef_column <- function(column_names, candidates) {
  if (is.null(column_names)) {
    return(NA_integer_)
  }

  for (candidate in candidates) {
    positions <- which(column_names == candidate)

    if (length(positions) > 1L) {
      return(NA_integer_)
    }
    if (length(positions) == 1L) {
      return(as.integer(positions))
    }
  }

  NA_integer_
}


#' Standardised Coefficient Matrix for the Summary Table
#'
#' Returns a standardised numeric matrix with columns `estimate`, `se`,
#' `statistic`, and `pvalue` and rows named by coefficient. Sourced from
#' the raw summary's coefficient table where possible (survreg calls it
#' `table`), falling back to `stats::vcov(object)` and Wald z or t inference
#' when a raw table is missing or in an unfamiliar layout. Also carries a
#' `"statistic_label"` attribute (`"z"`, `"t"`, or `"Statistic"`).
#'
#' @param object An `"mfp2"` model object.
#' @param raw_summary Raw underlying-class summary from
#'   `mfp2_summary_raw()`, possibly `NULL`.
#'
#' @return A numeric matrix with columns `estimate`, `se`, `statistic`, and
#'   `pvalue`, plus a `"statistic_label"` attribute.
#'
#' @keywords internal
#' @noRd
mfp2_summary_coef_matrix <- function(object, raw_summary) {
  # Try the raw summary's coefficient table first.
  cmat <- NULL
  if (!is.null(raw_summary) && !is.null(raw_summary$coefficients)) {
    cmat <- raw_summary$coefficients
  } else if (!is.null(raw_summary) && !is.null(raw_summary$table)) {
    # summary.survreg() calls this component `table` and appends one or more
    # scale rows after the regression coefficients.
    cmat <- raw_summary$table
  }

  est <- object$coefficients
  nm <- names(est)

  if (is.matrix(cmat) && !is.null(rownames(cmat)) &&
      !is.null(nm) && all(nm %in% rownames(cmat))) {
    cmat <- cmat[nm, , drop = FALSE]
  }

  if (is.matrix(cmat) && is.numeric(cmat) && nrow(cmat) == length(est)) {
    # Use exact known labels rather than grep()/partial matching. If any label
    # is absent or ambiguous, `columns` contains NA and execution continues to
    # the covariance fallback below instead of indexing cmat with NA.
    cn <- colnames(cmat)
    columns <- c(
      estimate = mfp2_summary_match_coef_column(
        cn,
        c("Estimate", "Value", "coef")
      ),
      se = mfp2_summary_match_coef_column(
        cn,
        c("robust se", "Std. Error", "se(coef)")
      ),
      statistic = mfp2_summary_match_coef_column(
        cn,
        c("z value", "t value", "z", "t")
      ),
      pvalue = mfp2_summary_match_coef_column(
        cn,
        c("Pr(>|z|)", "Pr(>|t|)", "p", "p-value", "pvalue")
      )
    )

    if (!anyNA(columns) && !anyDuplicated(columns)) {
      out <- cbind(
        estimate  = cmat[, columns[["estimate"]]],
        se        = cmat[, columns[["se"]]],
        statistic = cmat[, columns[["statistic"]]],
        pvalue    = cmat[, columns[["pvalue"]]]
      )
      rownames(out) <- rownames(cmat)
      statistic_column <- cn[columns[["statistic"]]]
      attr(out, "statistic_label") <- if (grepl("^t", statistic_column)) {
        "t"
      } else if (grepl("^z", statistic_column)) {
        "z"
      } else {
        "Statistic"
      }
      return(out)
    }
  }

  # Fallback for an absent, malformed, or unfamiliar raw coefficient table.
  # This is intentionally reached for unresolved headers rather than allowing
  # an NA column index to manufacture an all-NA result.
  is_ordinal <- mfp2_family_is_ordinal(object$family_string)
  V <- if (is_ordinal) {
    mfp2_ordinal_vcov_slopes(object)
  } else {
    tryCatch(stats::vcov(object), error = function(e) NULL)
  }
  se <- stats::setNames(rep(NA_real_, length(est)), nm)
  if (is_ordinal && is.matrix(V) && !is.null(nm) &&
      !is.null(rownames(V)) && !is.null(colnames(V))) {
    covariance_names <- intersect(
      nm,
      intersect(rownames(V), colnames(V))
    )
    if (length(covariance_names) > 0L) {
      variances <- diag(V[covariance_names, covariance_names, drop = FALSE])
      valid <- !is.na(variances) & is.finite(variances) & variances >= 0
      se[covariance_names[valid]] <- sqrt(variances[valid])
    }
  } else {
    if (is.matrix(V) && !is.null(nm) &&
        !is.null(rownames(V)) && !is.null(colnames(V)) &&
        all(nm %in% rownames(V)) && all(nm %in% colnames(V))) {
      V <- V[nm, nm, drop = FALSE]
    }
    if (is.matrix(V) && all(dim(V) == length(est))) {
      se <- sqrt(diag(V))
    }
  }
  stat <- est / se
  uses_t <- inherits(object, "glm") &&
    mfp2_glm_estimates_dispersion(object$family, object$family_string) &&
    is.numeric(object$df.residual) && length(object$df.residual) == 1L &&
    is.finite(object$df.residual)
  pval <- if (uses_t) {
    2 * stats::pt(abs(stat), df = object$df.residual, lower.tail = FALSE)
  } else {
    2 * stats::pnorm(abs(stat), lower.tail = FALSE)
  }

  out <- cbind(
    estimate  = est,
    se        = se,
    statistic = stat,
    pvalue    = pval
  )
  rownames(out) <- nm
  attr(out, "statistic_label") <- if (uses_t) {
    "t"
  } else {
    "z"
  }
  out
}

# ---------------------------------------------------------------------------
# Nonlinear terms table (LRT, one row per variable)
# ---------------------------------------------------------------------------

#' Cache the Metadata Needed to Refit Reduced Models
#'
#' Assembles the family-specific response, weights, offset, strata, control,
#' and full-model log-likelihood needed to refit reduced models when
#' computing likelihood-ratio tests for nonlinear terms. Doing this once at
#' the top of the summary avoids repeating expensive setup (in particular
#' `prepare_family_for_fit()` for `survreg` and unpacking the Fine--Gray
#' counting-process representation) for every nonlinear variable.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A named list with `valid` (flag), `full_logl`, `family`,
#'   `family_string`, `y`, `weights`, `offset`, `strata`, `method`,
#'   `fitter`, `control`, `nocenter`, and `is_finegray`. When `valid` is
#'   `FALSE`, downstream code should skip reduced-model refits and report
#'   `NA` test statistics.
#'
#' @keywords internal
#' @noRd
mfp2_summary_refit_context <- function(object) {
  family_string <- object$family_string
  is_cox <- identical(family_string, "cox")
  is_finegray <- identical(family_string, "finegray")
  is_survreg <- identical(family_string, "survreg")
  is_prepared_categorical <- family_string %in% c("ordinal", "multinomial")

  full_logl <- object$mfp_logl
  if (!is.numeric(full_logl) || length(full_logl) != 1L ||
      is.na(full_logl) || !is.finite(full_logl)) {
    full_logl <- tryCatch({
      if (is_cox || is_finegray || is_survreg) {
        unname(object$loglik[length(object$loglik)])
      } else {
        as.numeric(stats::logLik(object))
      }
    }, error = function(e) NA_real_)
  }

  fitter <- if (is.null(object$fitter)) "base" else object$fitter
  control <- object$mfp2_control
  if (is.null(control)) {
    control <- tryCatch(
      normalize_fit_control(
        control = NULL,
        family_string = family_string,
        fitter = fitter
      ),
      error = function(e) NULL
    )
  }

  nocenter_meta <- object$mfp2_nocenter
  nocenter <- if (is.list(nocenter_meta) &&
                  "value" %in% names(nocenter_meta)) {
    nocenter_meta[["value"]]
  } else if (is.numeric(nocenter_meta)) {
    nocenter_meta
  } else {
    c(-1, 0, 1)
  }

  context <- list(
    valid = is.numeric(full_logl) && length(full_logl) == 1L &&
      !is.na(full_logl) && is.finite(full_logl) && !is.null(control),
    full_logl = full_logl,
    family = if (is_prepared_categorical && !is.null(object$mfp2_family)) {
      object$mfp2_family
    } else {
      object$family
    },
    family_string = family_string,
    y = object$y,
    weights = if (is_cox || is_finegray || is_survreg) {
      object$weights
    } else {
      object$prior.weights
    },
    offset = object$offset,
    strata = if (is_cox) object$strata else NULL,
    method = if (is_cox || is_finegray) object$method else NULL,
    fitter = fitter,
    control = control,
    nocenter = nocenter,
    is_finegray = is_finegray
  )

  if (is_finegray) {
    # The retained final model already contains the expanded counting-process
    # response, design and weights. Reuse them without another finegray() call.
    row_map <- object$mfp2_finegray_row_map
    original_offset <- object$mfp2_original_offset
    if (is.null(row_map) || is.null(original_offset)) {
      context$valid <- FALSE
      return(context)
    }
    context$family <- "cox"
    context$family_string <- "cox"
    context$offset <- original_offset[row_map]
  } else if (is_survreg) {
    # Prepare the transformed response and resolved distribution once for all
    # nonlinear-term reduced refits in this summary call.
    prepared <- tryCatch(
      prepare_family_for_fit(
        family = object$family,
        family_string = "survreg",
        y = object$y_original,
        weights = object$weights,
        strata = object$mfp2_survreg_strata
      ),
      error = function(e) NULL
    )
    if (is.null(prepared)) {
      context$valid <- FALSE
      return(context)
    }
    context$family <- prepared$family
    context$y <- object$y_original
  }

  context
}


#' Build the Nonlinear-Terms Table for `print.summary.mfp2()`
#'
#' Assembles the data frame of one row per nonlinear term shown in the
#' printed summary, carrying the functional-form label, degrees of freedom
#' for the drop-variable test, likelihood-ratio chi-square, and p value.
#' Reduced models are fitted once per variable using the cached refit
#' context; families that cannot support the drop-variable refit return
#' `NA` statistics rather than aborting.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#'
#' @return A data frame with columns `variable`, `form`, `df`, `lr_chisq`,
#'   and `p`. Empty when the model has no nonlinear terms.
#'
#' @keywords internal
#' @noRd
mfp2_summary_nonlinear_table <- function(object, classified) {
  nl_vars <- classified$variable_names[classified$is_nonlinear]
  if (length(nl_vars) == 0L) {
    return(data.frame())
  }

  # Construct family-specific response/control metadata once, then reuse it for
  # every reduced model. In particular, survreg response preparation is not
  # repeated once per nonlinear variable.
  refit_context <- mfp2_summary_refit_context(object)

  rows <- lapply(nl_vars, function(v) {
    idx <- match(v, classified$variable_names)
    form <- mfp2_summary_form_label(
      classified$powers_by_var[[idx]],
      classified$acd[idx], classified$zero[idx],
      classified$catzero[idx], classified$spike[idx],
      classified$spike_dec[idx], TRUE
    )
    lrt <- mfp2_summary_lrt_drop_variable(
      object, classified, v, refit_context = refit_context
    )

    data.frame(
      variable = v,
      form     = form,
      df       = lrt$df,
      lr_chisq = lrt$lr,
      p        = lrt$p,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Drop-Variable Likelihood-Ratio Test for One Variable
#'
#' Performs the joint likelihood-ratio test for dropping every design column
#' belonging to one variable, holding all other functional forms fixed. The
#' selected model's stored log-likelihood is already on the same
#' family-specific scale, so only the reduced model is refitted. Degrees of
#' freedom use the selection-adjusted `df_final` from `fp_terms`.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#' @param v Character scalar naming the variable to drop.
#' @param refit_context Cached refit context from
#'   `mfp2_summary_refit_context()`. When `NULL`, it is recomputed.
#'
#' @return A named list `list(lr, df, p)` with the likelihood-ratio
#'   statistic, degrees of freedom, and p value. Any component may be `NA`
#'   when the reduced model cannot be fitted.
#'
#' @keywords internal
#' @noRd
mfp2_summary_lrt_drop_variable <- function(object, classified, v,
                                            refit_context = NULL) {
  cols_v <- intersect(colnames(object$x), classified$cols_by_var[[v]])
  idx <- match(v, classified$variable_names)
  df_v <- classified$df_final[idx]
  if (is.na(df_v)) df_v <- length(cols_v)

  na_result <- list(lr = NA_real_, df = df_v, p = NA_real_)
  x_full <- object$x
  if (is.null(x_full) || length(cols_v) == 0L) {
    return(na_result)
  }

  if (is.null(refit_context)) {
    refit_context <- mfp2_summary_refit_context(object)
  }
  if (!isTRUE(refit_context$valid)) return(na_result)

  keep <- setdiff(colnames(x_full), cols_v)
  x_reduced <- x_full[, keep, drop = FALSE]

  reduced_fit <- if (isTRUE(refit_context$is_finegray)) {
    # agreg.fit() is the counting-process matrix fitter. With resid = FALSE it
    # does not use row names, so avoid allocating them for summary refits too.
    tryCatch({
      fg_fit <- mfp2_with_cox_convergence_guard(
        survival::agreg.fit(
          x = x_reduced,
          y = refit_context$y,
          strata = NULL,
          offset = refit_context$offset,
          init = NULL,
          control = refit_context$control,
          weights = refit_context$weights,
          method = if (is.null(refit_context$method)) {
            "breslow"
          } else {
            refit_context$method
          },
          rownames = NULL,
          resid = FALSE,
          nocenter = refit_context$nocenter
        ),
        fast = TRUE
      )
      list(logl = unname(fg_fit$loglik[length(fg_fit$loglik)]))
    }, error = function(e) NULL)
  } else {
    x_has_intercept <- NCOL(x_reduced) > 0L &&
      !is.null(colnames(x_reduced)) &&
      identical(colnames(x_reduced)[1L], "(Intercept)")
    tryCatch(
      fit_model(
        x = x_reduced,
        y = refit_context$y,
        family = refit_context$family,
        family_string = refit_context$family_string,
        fitter = refit_context$fitter,
        weights = refit_context$weights,
        offset = refit_context$offset,
        method = refit_context$method,
        strata = refit_context$strata,
        control = refit_context$control,
        nocenter = refit_context$nocenter,
        x_has_intercept = x_has_intercept,
        keep_coefficients = FALSE,
        fast = TRUE
      ),
      error = function(e) NULL
    )
  }

  if (is.null(reduced_fit) || is.null(reduced_fit$logl) ||
      !is.finite(reduced_fit$logl)) {
    return(na_result)
  }

  lr <- 2 * (refit_context$full_logl - reduced_fit$logl)
  if (!is.finite(lr)) return(na_result)
  lr <- max(lr, 0)
  p <- stats::pchisq(lr, df = df_v, lower.tail = FALSE)

  list(lr = lr, df = df_v, p = p)
}

# ---------------------------------------------------------------------------
# Final design-column description
# ---------------------------------------------------------------------------

#' Compact Labels for a Fractional-Polynomial Basis
#'
#' Builds compact printable labels for the columns of an ordinary
#' fractional-polynomial basis, on the scale used by the final model. A
#' nonzero shift is shown explicitly (`(x + s)` or `(x - s)`), and repeated
#' powers use the standard FP `log(...)` multiplier convention.
#'
#' @param term Character scalar with the predictor name that appears in the
#'   label.
#' @param powers Finite numeric vector of FP powers.
#' @param shift Numeric scalar shift applied to `term`. Default `0`.
#'
#' @return Character vector of labels, one per element of `powers`.
#'
#' @keywords internal
#' @noRd
mfp2_fp_basis_labels <- function(term, powers, shift = 0) {
  if (!is.numeric(powers) || length(powers) == 0L || anyNA(powers) ||
      any(!is.finite(powers))) {
    stop("`powers` must be a finite non-empty numeric vector.", call. = FALSE)
  }

  if (length(shift) != 1L || is.na(shift) || !is.finite(shift)) shift <- 0

  base <- if (isTRUE(shift == 0)) {
    term
  } else if (shift > 0) {
    paste0("(", term, " + ", format(abs(shift), trim = TRUE, scientific = FALSE), ")")
  } else {
    paste0("(", term, " - ", format(abs(shift), trim = TRUE, scientific = FALSE), ")")
  }

  log_arg <- if (isTRUE(shift == 0)) {
    term
  } else if (shift > 0) {
    paste0(term, " + ", format(abs(shift), trim = TRUE, scientific = FALSE))
  } else {
    paste0(term, " - ", format(abs(shift), trim = TRUE, scientific = FALSE))
  }

  out <- character(length(powers))
  seen <- numeric(0L)

  for (i in seq_along(powers)) {
    p <- powers[[i]]
    repetition <- sum(seen == p) + 1L

    first <- if (isTRUE(p == 0)) {
      paste0("log(", log_arg, ")")
    } else if (isTRUE(p == 1)) {
      base
    } else {
      paste0(base, "^", format(p, trim = TRUE, scientific = FALSE))
    }

    if (repetition == 1L) {
      out[[i]] <- first
    } else if (isTRUE(p == 0)) {
      out[[i]] <- paste0("log(", log_arg, ")^", repetition)
    } else if (repetition == 2L) {
      out[[i]] <- paste0(first, " * log(", log_arg, ")")
    } else {
      out[[i]] <- paste0(first, " * log(", log_arg, ")^", repetition - 1L)
    }

    seen <- c(seen, p)
  }

  out
}

#' Describe One Factor Design Column
#'
#' Formats a single factor design column from its stored level-by-design
#' matrix. Treatment-coded indicator columns are shown as level indicators
#' (for example `I(x = "high")`). Other contrast columns retain their
#' model-matrix names, such as `x4.L` and `x4.Q`, because those names
#' identify the configured contrast basis directly and are the labels users
#' expect to see.
#'
#' @param object An `"mfp2"` model object with
#'   `object$formula_factor_info` populated.
#' @param term Character scalar naming the factor term.
#' @param source Character scalar naming the raw source column.
#'
#' @return A named list `list(variable, basis)`, or `NULL` when the column
#'   is not a recognised factor column.
#'
#' @keywords internal
#' @noRd
mfp2_factor_design_column_info <- function(object, term, source) {
  factor_info <- if (!is.null(object$formula_factor_info)) {
    object$formula_factor_info[[term]]
  } else {
    NULL
  }

  if (is.null(factor_info) || is.null(factor_info$design_by_level)) {
    return(NULL)
  }

  design <- as.matrix(factor_info$design_by_level)
  if (is.null(colnames(design)) || !source %in% colnames(design) ||
      is.null(rownames(design))) {
    return(NULL)
  }

  values <- suppressWarnings(as.numeric(design[, source]))
  if (length(values) != nrow(design) || anyNA(values) || any(!is.finite(values))) {
    return(NULL)
  }

  variable <- factor_info$variable
  if (is.null(variable) || length(variable) != 1L || is.na(variable) ||
      !nzchar(variable)) {
    variable <- term
  }

  tolerance <- 1e-10
  is_zero <- abs(values) <= tolerance
  is_one <- abs(values - 1) <= tolerance

  if (all(is_zero | is_one) && sum(is_one) == 1L) {
    level <- rownames(design)[which(is_one)]
    basis <- sprintf(
      "I(%s = %s)",
      variable,
      encodeString(level, quote = "\"")
    )
  } else {
    basis <- source
  }

  list(variable = variable, basis = basis)
}

#' Describe the Final Transformed Design Columns
#'
#' Builds a data frame describing every transformed design column of the
#' final fitted model, using the metadata recorded when the design matrix
#' was assembled. Source variables, component types (FP basis, ACD basis,
#' zero indicator, identity binary), zero handling, centering, and
#' model-column mappings are read directly from the fitted object; the
#' generated transformed column names are not interpreted.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A data frame with columns `variable`, `transformed_column`,
#'   `model_column`, `basis`, `center`, `centered`, `zero_handled`, and
#'   `component`, one row per transformed column. Returns an empty data
#'   frame for an intercept-only model.
#'
#' @keywords internal
#' @noRd
mfp2_design_column_info <- function(object) {
  transformed_to_model <- object$transformed_to_model_columns
  transformed_to_source <- object$transformed_column_to_source
  transformed_component <- object$transformed_column_component
  transformed_zero_handled <- object$transformed_column_zero_handled
  transformed_centered <- object$transformed_column_centered
  term_to_columns <- object$term_to_columns

  # An intercept-only final model has no transformed predictor columns. The
  # fitting backend records this as an empty named transformed-to-model map,
  # while the transformation-specific metadata are naturally NULL because no
  # design matrix was constructed. Accept that empty design only when the
  # fitted coefficient vector also contains no predictor coefficients; this
  # keeps a genuinely incomplete nonempty design from being silently hidden.
  coefficient_names <- if (is.matrix(object$coefficients)) {
    colnames(object$coefficients)
  } else {
    names(object$coefficients)
  }
  if (is.null(coefficient_names)) {
    coefficient_names <- character(0L)
  }
  predictor_coefficients <- setdiff(coefficient_names, "(Intercept)")

  if (!is.null(transformed_to_model) &&
      length(transformed_to_model) == 0L &&
      !is.null(names(transformed_to_model)) &&
      length(predictor_coefficients) == 0L) {
    return(data.frame(
      variable = character(0L),
      transformed_column = character(0L),
      model_column = character(0L),
      basis = character(0L),
      center = numeric(0L),
      centered = logical(0L),
      zero_handled = logical(0L),
      component = character(0L),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }

  metadata <- list(
    transformed_to_model = transformed_to_model,
    transformed_to_source = transformed_to_source,
    transformed_component = transformed_component,
    transformed_zero_handled = transformed_zero_handled,
    transformed_centered = transformed_centered
  )

  valid_named_vector <- function(value) {
    !is.null(value) && !is.null(names(value)) && !anyDuplicated(names(value))
  }

  if (!all(vapply(metadata, valid_named_vector, logical(1L))) ||
      is.null(term_to_columns) || is.null(names(term_to_columns))) {
    stop("Final design-column metadata is incomplete.", call. = FALSE)
  }

  transformed_columns <- names(transformed_to_model)
  aligned <- vapply(metadata[-1L], function(value) {
    setequal(transformed_columns, names(value))
  }, logical(1L))
  if (!all(aligned)) {
    stop("Final design-column metadata is not aligned.", call. = FALSE)
  }

  components_allowed <- c(
    "fp_basis", "acd_basis", "zero_indicator", "identity_binary"
  )
  unknown_components <- setdiff(
    unique(unname(transformed_component[transformed_columns])),
    components_allowed
  )
  if (length(unknown_components) > 0L) {
    stop(
      sprintf(
        "Unknown transformed-column component(s): %s.",
        paste(unknown_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Each raw source column must belong to exactly one conceptual model term.
  source_to_term <- character(0L)
  for (term in names(term_to_columns)) {
    raw_columns <- as.character(term_to_columns[[term]])
    for (raw_column in raw_columns) {
      if (raw_column %in% names(source_to_term) &&
          !identical(unname(source_to_term[[raw_column]]), term)) {
        stop(
          sprintf("Source column '%s' belongs to more than one model term.", raw_column),
          call. = FALSE
        )
      }
      source_to_term[[raw_column]] <- term
    }
  }

  sources <- unname(transformed_to_source[transformed_columns])
  missing_sources <- is.na(sources) | !nzchar(sources)
  if (any(missing_sources)) {
    stop("Final design metadata contains an unmapped source column.", call. = FALSE)
  }

  missing_terms <- setdiff(unique(sources), names(source_to_term))
  if (length(missing_terms) > 0L) {
    stop(
      sprintf(
        "Final design metadata contains unmapped source column(s): %s.",
        paste(missing_terms, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  centers <- object$centers
  if (!is.null(centers)) {
    if (is.null(names(centers)) || anyDuplicated(names(centers)) ||
        !setequal(transformed_columns, names(centers))) {
      stop("Final centering constants are not aligned with design columns.", call. = FALSE)
    }
  }

  fp_terms <- object$fp_terms
  power_columns <- if (!is.null(fp_terms)) {
    grep("^power[0-9]+$", names(fp_terms), value = TRUE)
  } else {
    character(0L)
  }

  rows <- vector("list", length(transformed_columns))

  for (i in seq_along(transformed_columns)) {
    transformed_column <- transformed_columns[[i]]
    source <- transformed_to_source[[transformed_column]]
    component <- transformed_component[[transformed_column]]
    zero_handled <- isTRUE(transformed_zero_handled[[transformed_column]])
    centered <- isTRUE(transformed_centered[[transformed_column]])
    term <- source_to_term[[source]]
    raw_columns <- as.character(term_to_columns[[term]])

    mapped_block <- term_uses_column_mapping(term, raw_columns)
    display_variable <- if (mapped_block) term else source
    basis <- source

    factor_column <- if (mapped_block) {
      mfp2_factor_design_column_info(object, term, source)
    } else {
      NULL
    }

    if (!is.null(factor_column)) {
      display_variable <- factor_column$variable
      basis <- factor_column$basis
    } else if (mapped_block) {
      # Explicit grouped design blocks retain their supplied source-column
      # labels when there is no factor level-to-design mapping.
      basis <- source
    } else if (identical(component, "zero_indicator")) {
      basis <- sprintf("I(%s = 0)", term)
    } else if (identical(component, "identity_binary")) {
      basis <- term
    } else {
      if (is.null(fp_terms) || is.null(rownames(fp_terms)) ||
          !term %in% rownames(fp_terms)) {
        stop(sprintf("Missing FP metadata for term '%s'.", term), call. = FALSE)
      }

      slots <- suppressWarnings(as.numeric(unlist(
        fp_terms[term, power_columns, drop = FALSE],
        use.names = FALSE
      )))
      term_row <- fp_terms[term, , drop = FALSE]
      is_acd <- mfp2_summary_flag(term_row, "acd")[[1L]]

      shift <- 0
      transformations <- object$transformations
      if (!is.null(transformations) && term %in% rownames(transformations) &&
          "shift" %in% colnames(transformations)) {
        shift_value <- suppressWarnings(as.numeric(transformations[term, "shift"]))
        if (length(shift_value) == 1L && !is.na(shift_value) &&
            is.finite(shift_value)) {
          shift <- shift_value
        }
      }

      if (zero_handled) {
        shift <- 0
      }

      if (identical(component, "acd_basis")) {
        power <- if (length(slots) >= 2L) slots[[2L]] else NA_real_
        basis <- sprintf("A(%s)", term)
        if (!is.na(power)) {
          basis <- mfp2_fp_basis_labels(basis, power, shift = 0)[[1L]]
        }
      } else if (identical(component, "fp_basis")) {
        if (isTRUE(is_acd)) {
          powers <- if (length(slots) >= 1L) slots[[1L]] else NA_real_
          powers <- powers[!is.na(powers)]
        } else {
          powers <- slots[!is.na(slots)]
        }

        if (length(powers) == 0L) {
          stop(sprintf("Missing FP power metadata for term '%s'.", term), call. = FALSE)
        }

        fp_columns_for_source <- transformed_columns[
          unname(transformed_to_source[transformed_columns]) == source &
            unname(transformed_component[transformed_columns]) == "fp_basis"
        ]
        k <- match(transformed_column, fp_columns_for_source)
        labels <- mfp2_fp_basis_labels(term, powers, shift = shift)

        if (is.na(k) || k > length(labels)) {
          stop(
            sprintf("FP basis metadata is not aligned for term '%s'.", term),
            call. = FALSE
          )
        }

        basis <- labels[[k]]
        if (zero_handled) {
          basis <- sprintf("I(%s > 0) * %s", term, basis)
        }
      }
    }

    rows[[i]] <- data.frame(
      variable = display_variable,
      transformed_column = transformed_column,
      model_column = unname(transformed_to_model[[transformed_column]]),
      basis = basis,
      center = if (is.null(centers)) NA_real_ else unname(centers[[transformed_column]]),
      centered = centered,
      zero_handled = zero_handled,
      component = component,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Basis coefficients table (optional)
# ---------------------------------------------------------------------------

#' Build the Optional Basis-Coefficient Table
#'
#' Builds the per-basis coefficient table shown for nonlinear terms when
#' basis details are requested in the printed summary. One row per fitted
#' basis coefficient is returned; the `variable` column is filled only for
#' the first row of each variable so the printed table stays compact.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#'
#' @return A data frame with columns `variable`, `term`, and `coef`, or
#'   `NULL` when no nonlinear terms are present or no coefficients match.
#'
#' @keywords internal
#' @noRd
mfp2_summary_basis_table <- function(object, classified) {
  nl_vars <- classified$variable_names[classified$is_nonlinear]
  if (length(nl_vars) == 0L) return(NULL)

  coefs <- object$coefficients
  rows <- list()
  for (v in nl_vars) {
    cols_v <- intersect(names(coefs), classified$cols_by_var[[v]])
    for (k in seq_along(cols_v)) {
      rows[[length(rows) + 1L]] <- data.frame(
        variable = if (k == 1L) v else "",
        term     = mfp2_summary_term_label(object, classified, v, cols_v[k]),
        coef     = unname(coefs[cols_v[k]]),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

#' Readable Transformation Label for One Fitted Column
#'
#' Constructs the printable transformation label used in the basis-details
#' table for a single fitted column. Ordinary final FP coefficients multiply
#' a shifted-but-unscaled basis, because `fit_mfp()` backscales the working
#' predictors before the final transformation. ACD component columns are
#' functions of `A(x)`; new fits use shifted, unscaled input for `A(x)`, and
#' a stored non-unit ACD scale is honoured only so that formula output
#' remains exact for legacy serialised model objects.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#' @param v Character scalar naming the variable.
#' @param col Character scalar naming the fitted design column.
#'
#' @return Character scalar with the transformation label.
#'
#' @keywords internal
#' @noRd
mfp2_summary_term_label <- function(object, classified, v, col) {
  idx <- match(v, classified$variable_names)
  cols_v <- classified$cols_by_var[[v]]
  k <- match(col, cols_v)
  acd <- classified$acd[idx]

  base_col <- sub("\\.[0-9]+$", "", col)
  is_acd_component <- isTRUE(acd) && isTRUE(base_col == paste0("A_", v))

  if (isTRUE(acd)) {
    slots <- classified$power_slots_by_var[[idx]]
    if (length(slots) < 2L) slots <- c(slots, rep(NA_real_, 2L - length(slots)))

    if (is_acd_component) {
      power <- slots[2L]
      base <- sprintf("A(%s)", v)
    } else {
      power <- slots[1L]
      ss <- mfp2_summary_shift_scale(object, v)
      base <- mfp2_summary_fp_inner_expr(v, ss$shift)
    }

    if (is.na(power)) return(base)
    return(mfp2_summary_power_expr(base, power))
  }

  powers <- classified$powers_by_var[[idx]]
  if (length(powers) == 0L) {
    return(sprintf("I(%s > 0)", v))
  }

  ss <- mfp2_summary_shift_scale(object, v)
  base <- mfp2_summary_fp_inner_expr(v, ss$shift)
  power <- if (k <= length(powers)) powers[k] else powers[length(powers)]

  # Repeated powers have the special FP2 second basis x^p * log(x). This rule
  # applies only to ordinary FP2 terms; equal ACD powers act on different inputs
  # (x and A(x)) and must not trigger the repeated-power construction.
  repeated <- length(powers) == 2L && isTRUE(powers[1L] == powers[2L])
  mfp2_summary_power_expr(base, power, repeated = repeated && k == 2L)
}

#' Apply an FP Power to a Base Expression
#'
#' Wraps a readable base expression with one FP power. Power `0` becomes
#' `log(base)`; power `1` returns the base unchanged; other powers become
#' `base^power`. When `repeated = TRUE`, the standard FP2 second-basis
#' multiplier `* log(base)` is appended.
#'
#' @param base Character scalar with the base expression.
#' @param power Numeric FP power.
#' @param repeated Logical. If `TRUE`, append the repeated-power log
#'   multiplier.
#'
#' @return Character scalar with the powered expression.
#'
#' @keywords internal
#' @noRd
mfp2_summary_power_expr <- function(base, power, repeated = FALSE) {
  wrapped <- sprintf("(%s)", base)

  out <- if (isTRUE(power == 0)) {
    sprintf("log%s", wrapped)
  } else if (isTRUE(power == 1)) {
    wrapped
  } else {
    sprintf("%s^(%s)", wrapped, format(power, trim = TRUE))
  }

  if (isTRUE(repeated)) {
    out <- sprintf("%s*log(%s)", out, base)
  }
  out
}

#' Inner Expression for the Final Ordinary FP Basis
#'
#' Builds the parenthesised inner expression used by `mfp2_summary_power_expr()`
#' for an ordinary FP term. The final fit uses shifted but unscaled
#' predictors, so preprocessing scale is intentionally absent here.
#'
#' @param v Character scalar naming the variable.
#' @param shift Numeric shift, possibly `NA`.
#'
#' @return Character scalar such as `"(x)"` or `"(x + 0.5)"`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_fp_inner_expr <- function(v, shift) {
  if (!is.na(shift) && shift != 0) {
    sprintf("(%s + %s)", v, format(shift, trim = TRUE))
  } else {
    sprintf("(%s)", v)
  }
}

#' Inner Expression for a Stored ACD Approximation
#'
#' Builds the inner expression used only by the stored ACD approximation
#' formula. New model objects store `scale = 1` because ACD variables are
#' not scaled, but a non-unit stored scale is honoured so summaries of
#' legacy serialised model objects remain exact.
#'
#' @param v Character scalar naming the variable.
#' @param shift Numeric shift, possibly `NA`.
#' @param scale Numeric scale, possibly `NA`.
#'
#' @return Character scalar with the inner expression.
#'
#' @keywords internal
#' @noRd
mfp2_summary_acd_inner_expr <- function(v, shift, scale) {
  base <- if (!is.na(shift) && shift != 0) {
    sprintf("(%s + %s)", v, format(shift, trim = TRUE))
  } else {
    sprintf("(%s)", v)
  }

  if (!is.na(scale) && scale != 1) {
    base <- sprintf("%s/%s", base, format(scale, trim = TRUE))
  }
  base
}

#' Look Up Shift and Scale for One Variable
#'
#' Reads the stored preprocessing shift and scale values for one variable
#' from the fitted model's `transformations` table.
#'
#' @param object An `"mfp2"` model object.
#' @param v Character scalar naming the variable.
#'
#' @return A named list `list(shift, scale)` of numeric scalars; either
#'   entry is `NA_real_` when the variable is absent from the table or the
#'   column is missing.
#'
#' @keywords internal
#' @noRd
mfp2_summary_shift_scale <- function(object, v) {
  tr <- object$transformations
  shift <- NA_real_
  scale <- NA_real_
  if (!is.null(tr) && v %in% rownames(tr)) {
    if ("shift" %in% colnames(tr)) shift <- suppressWarnings(as.numeric(tr[v, "shift"]))
    if ("scale" %in% colnames(tr)) scale <- suppressWarnings(as.numeric(tr[v, "scale"]))
  }
  list(shift = shift, scale = scale)
}

# ---------------------------------------------------------------------------
# Fitted-function formula strings (optional)
# ---------------------------------------------------------------------------

#' Fitted-Function Formula Strings for Nonlinear Terms
#'
#' Builds the printable fitted-function formula strings shown for nonlinear
#' terms in the printed summary. Each formula is of the form
#' `f(x) = b1 * label1 + b2 * label2 + ...`, using `signif()` rounding for
#' coefficients so the display remains readable.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#'
#' @return Character vector of formula strings, one per nonlinear term, or
#'   `NULL` when the model contains no nonlinear terms.
#'
#' @keywords internal
#' @noRd
mfp2_summary_formula_strings <- function(object, classified) {
  nl_vars <- classified$variable_names[classified$is_nonlinear]
  if (length(nl_vars) == 0L) return(NULL)

  coefs <- object$coefficients
  out <- character(length(nl_vars))
  for (i in seq_along(nl_vars)) {
    v <- nl_vars[i]
    cols_v <- intersect(names(coefs), classified$cols_by_var[[v]])
    terms <- vapply(cols_v, function(col) {
      b <- unname(coefs[col])
      lbl <- mfp2_summary_term_label(object, classified, v, col)
      sprintf("%s * %s", format(signif(b, 4), trim = TRUE), lbl)
    }, character(1L))
    body <- paste(terms, collapse = " + ")
    # Tidy "+ -" into "- ".
    body <- gsub("\\+ -", "- ", body)
    out[i] <- sprintf("f(%s) = %s", v, body)
  }
  out
}

#' Build ACD Definition Strings for the Printed Summary
#'
#' Assembles the printable definitions `A(x) = pnorm(beta0 + beta1 * ...)`
#' for every active ACD component. Both the ordinary FP term and new ACD
#' fits use shifted, unscaled predictor values; a non-unit stored ACD scale
#' is honoured only for compatibility with legacy fitted objects.
#'
#' @param object An `"mfp2"` model object.
#' @param classified Classification result from
#'   `mfp2_summary_classify_terms()`.
#'
#' @return Character vector of ACD definition strings, or `NULL` when no
#'   ACD components are active.
#'
#' @keywords internal
#' @noRd
mfp2_summary_acd_definitions <- function(object, classified) {
  acd_vars <- classified$variable_names[
    classified$is_nonlinear & classified$acd
  ]
  if (length(acd_vars) == 0L) return(NULL)

  definitions <- vapply(acd_vars, function(v) {
    cols_v <- classified$cols_by_var[[v]]
    base_cols <- sub("\\.[0-9]+$", "", cols_v)
    if (!any(base_cols == paste0("A_", v))) return(NA_character_)
    mfp2_summary_acd_definition(object, v)
  }, character(1L))

  definitions <- definitions[!is.na(definitions) & nzchar(definitions)]
  if (length(definitions) == 0L) NULL else unname(definitions)
}

#' ACD Definition String for One Variable
#'
#' Builds the printable ACD definition `A(v) = pnorm(beta0 + beta1 * ...)`
#' for one variable, from the stored ACD parameter list. Returns
#' `NA_character_` when any required piece is missing.
#'
#' @param object An `"mfp2"` model object.
#' @param v Character scalar naming the variable.
#'
#' @return Character scalar with the definition, or `NA_character_`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_acd_definition <- function(object, v) {
  par <- object$acd_parameter[[v]]
  if (is.null(par)) return(NA_character_)

  required <- c("beta0", "beta1", "power", "shift", "scale")
  if (!all(required %in% names(par))) return(NA_character_)

  ss <- mfp2_summary_shift_scale(object, v)
  preprocessing_shift <- if (is.na(ss$shift)) 0 else ss$shift
  acd_shift <- suppressWarnings(as.numeric(par$shift)[1L])
  acd_scale <- suppressWarnings(as.numeric(par$scale)[1L])
  acd_power <- suppressWarnings(as.numeric(par$power)[1L])
  beta0 <- suppressWarnings(as.numeric(par$beta0)[1L])
  beta1 <- suppressWarnings(as.numeric(par$beta1)[1L])

  if (anyNA(c(acd_shift, acd_scale, acd_power, beta0, beta1))) {
    return(NA_character_)
  }

  total_shift <- preprocessing_shift + acd_shift
  inner <- mfp2_summary_acd_inner_expr(v, total_shift, acd_scale)
  acd_basis <- mfp2_summary_power_expr(inner, acd_power)

  sprintf(
    "A(%s) = pnorm(%s + %s * %s)",
    v,
    format(signif(beta0, 4), trim = TRUE),
    format(signif(beta1, 4), trim = TRUE),
    acd_basis
  )
}

# ---------------------------------------------------------------------------
# Model-fit statistics
# ---------------------------------------------------------------------------

#' Model-Fit Statistics Wrapper
#'
#' Returns the values used by the shared Model Fit block renderer.
#' Everything else previously in this list (LR test statistic, R-squared,
#' etc.) has been removed: the block is now purely descriptive, and
#' inferential statements are made only by the per-variable nonlinear LRTs.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A named list with a single element `model_fit` from
#'   `mfp2_summary_model_fit_values()`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_fit_stats <- function(object) {
  list(
    model_fit = mfp2_summary_model_fit_values(object)
  )
}

#' Coerce a Scalar Value to a Finite Numeric
#'
#' Convenience helper that coerces a scalar value to `numeric`, returning
#' `NA_real_` for `NULL`, length mismatches, or coercion failures.
#'
#' @param x Value to coerce.
#'
#' @return Numeric scalar, possibly `NA_real_`.
#'
#' @keywords internal
#' @noRd
mfp2_summary_num <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NA_real_)
  suppressWarnings(as.numeric(x))
}

#' FP-Adjusted Model Degrees of Freedom
#'
#' Returns the total FP-adjusted model degrees of freedom: the sum of
#' `df_final` over the selected variables. Each variable is charged its
#' selection-adjusted df (1 linear, 2 for FP1, 4 for FP2, etc.) rather than
#' its raw coefficient count, so the overall LR test is consistent with the
#' per-variable nonlinear LRTs and with the MFP df accounting used
#' elsewhere. The intercept is not counted, because it is already present
#' in the null model.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return Integer scalar with the FP-adjusted model degrees of freedom.
#'   When `fp_terms$df_final` is unavailable, a fallback count based on
#'   `coef()` (excluding an intercept, when present) is returned.
#'
#' @keywords internal
#' @noRd
mfp2_summary_model_df <- function(object) {
  fp_terms <- object$fp_terms
  if (is.null(fp_terms) || !"df_final" %in% names(fp_terms)) {
    # Fallback: coefficient count excluding an intercept, if present.
    nm <- names(object$coefficients)
    n <- length(nm)
    if (!is.null(nm) && "(Intercept)" %in% nm) n <- n - 1L
    return(as.integer(n))
  }
  selected <- if ("selected" %in% names(fp_terms)) {
    as.logical(fp_terms[["selected"]])
  } else {
    rep(TRUE, nrow(fp_terms))
  }
  selected[is.na(selected)] <- FALSE
  df_final <- suppressWarnings(as.numeric(fp_terms[["df_final"]]))
  total <- sum(df_final[selected], na.rm = TRUE)
  as.integer(total)
}

#' Assemble the Model Fit Block Values
#'
#' Assembles the numeric values for the printed Model Fit block from a
#' fitted `"mfp2"` object. Used by both [print.mfp2()] and
#' [print.summary.mfp2()] so the two methods display identical
#' family-specific fit statistics.
#'
#' The `df` column follows a single convention: the number of regression
#' coefficients, excluding the intercept. `linear_df` as stored comes from
#' `fit_model()$df`, which counts the intercept for models that have one
#' and includes any estimated GLM dispersion, negative-binomial `theta`, or
#' `survreg` scale parameter; those adjustments are stripped here so the
#' number matches the promise made in the printed note.
#'
#' @param object An `"mfp2"` model object.
#'
#' @return A data frame with columns `label`, `fit_statistic`, and `df` and
#'   rows `"Full linear model"` and `"MFP model"`. Carries a
#'   `"statistic_label"` attribute (`"Deviance"` for GLMs, `"-2 log L"` for
#'   survival and multinomial models), and, where applicable, `"n_logits"`
#'   or `"n_ordinal_intercepts"` attributes used to render the df note.
#'
#' @keywords internal
#' @noRd
mfp2_summary_model_fit_values <- function(object) {
  family_string <- object$family_string
  fit_statistic <- c(
    mfp2_summary_num(object$linear_deviance),
    mfp2_summary_num(object$mfp_deviance)
  )

  # The df column follows a single convention: number of regression
  # coefficients, EXCLUDING the intercept. `linear_df` as stored comes from
  # fit_model()$df, which counts the intercept for models that have one and
  # includes any estimated GLM dispersion, negative-binomial theta, or `survreg`
  # scale parameter. These adjustments are stripped here so the number the
  # user sees matches the promise made in the note.
  linear_df <- if (!is.null(object$linear_df)) {
    d <- as.integer(object$linear_df)
    if (identical(family_string, "negbin")) {
      # Strip one for the intercept and one for theta.
      d <- d - 2L
    } else if (identical(family_string, "survreg")) {
      # survreg$idf counts the intercept plus estimated scale parameter(s).
      d <- d - as.integer(if (is.null(object$idf)) 1L else object$idf)
    } else if (identical(family_string, "multinomial")) {
      d <- d - as.integer(if (is.null(object$n_logits)) 1L else object$n_logits)
    } else if (mfp2_family_is_ordinal(family_string)) {
      d <- d - length(object$mfp2_ordinal_intercepts)
    } else if (mfp2_family_has_intercept(family_string)) {
      # Strip the intercept and, when applicable, estimated GLM dispersion.
      d <- d - 1L - as.integer(mfp2_glm_estimates_dispersion(
        family = object$family,
        family_string = family_string
      ))
    }
    d
  } else {
    NA_integer_
  }
  if (!is.na(linear_df) && linear_df < 0L) linear_df <- 0L

  # The MFP row uses the FP-adjusted convention (sum of df_final over selected
  # variables), which is already intercept- and sigma-free by construction.
  mfp_df <- mfp2_summary_model_df(object)

  values <- data.frame(
    label         = c("Full linear model", "MFP model"),
    fit_statistic = fit_statistic,
    df            = c(linear_df, mfp_df),
    stringsAsFactors = FALSE
  )
  attr(values, "statistic_label") <- if (
    mfp2_family_is_survival(family_string) ||
    identical(family_string, "multinomial")
  ) {
    "-2 log L"
  } else {
    "Deviance"
  }
  if (identical(family_string, "multinomial")) {
    attr(values, "n_logits") <- as.integer(object$n_logits)
  } else if (mfp2_family_is_ordinal(family_string)) {
    attr(values, "n_ordinal_intercepts") <- length(object$mfp2_ordinal_intercepts)
  }
  values
}

#' Print the Model Fit Block
#'
#' Prints the shared Model Fit block: heading, two-row table (Full linear /
#' MFP), and family-specific df note. Shared by [print.mfp2()] and
#' [print.summary.mfp2()] so the two methods produce identical output.
#' Callers pass a pre-computed values data frame to avoid recomputing it.
#'
#' `heading_printer` is a function of one string that draws the section
#' heading in the caller's own style, so [print.mfp2()] can reuse its
#' boxed-rule helper without leaking those internals into the summary path.
#' The summary printer supplies its own rule-based heading.
#'
#' @param values Values data frame from
#'   `mfp2_summary_model_fit_values()`.
#' @param digits Integer scalar controlling the number of decimal places
#'   used when formatting the fit statistic.
#' @param heading_printer Function of one string that prints the section
#'   heading.
#' @param notes Logical. If `TRUE`, print the explanatory df note.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of
#'   printing.
#'
#' @keywords internal
#' @noRd
mfp2_format_model_fit_block <- function(values, digits, heading_printer,
                                        notes = TRUE) {
  heading_printer("Model Fit")

  # Model-fit statistics use a fixed number of decimal places. The `digits`
  # argument controls decimal places here, rather than significant digits, so
  # values in the Deviance / -2 log L column remain vertically consistent.
  decimal_places <- suppressWarnings(as.integer(digits[1L]))
  if (length(decimal_places) != 1L ||
      is.na(decimal_places) ||
      decimal_places < 0L) {
    decimal_places <- 3L
  }

  # Right-align numeric columns; left-align the label. Column widths are
  # chosen from the widest formatted value so the note below reads under the
  # correct table width regardless of magnitude.
  statistic <- vapply(values$fit_statistic, function(v) {
    if (is.na(v) || !is.finite(v)) return("NA")
    formatC(v, format = "f", digits = decimal_places)
  }, character(1L))
  df_fmt <- vapply(values$df, function(v) {
    if (is.na(v)) return("NA")
    format(v, trim = TRUE)
  }, character(1L))

  statistic_label <- attr(values, "statistic_label", exact = TRUE)
  if (is.null(statistic_label) || length(statistic_label) != 1L) {
    statistic_label <- "Deviance"
  }

  label_width <- max(nchar(values$label))
  statistic_width <- max(nchar(statistic_label), max(nchar(statistic)))
  df_width <- max(nchar("df"), max(nchar(df_fmt)))

  # Column separators: two spaces after the label, four spaces before df.
  header_fmt <- sprintf(
    "%%-%ds  %%%ds    %%%ds\n",
    label_width, statistic_width, df_width
  )
  row_fmt <- header_fmt

  cat(sprintf(header_fmt, "", statistic_label, "df"))
  for (i in seq_len(nrow(values))) {
    cat(sprintf(row_fmt, values$label[i], statistic[i], df_fmt[i]))
  }
  cat("\n")

  # Note about the df column. Kept identical between print.mfp2() and
  # print.summary.mfp2() so users see one consistent explanation.
  if (isTRUE(notes)) {
    n_logits <- attr(values, "n_logits", exact = TRUE)
    n_ordinal_intercepts <- attr(
      values, "n_ordinal_intercepts", exact = TRUE
    )
    if (!is.null(n_logits)) {
      q <- n_logits
      note_lines <- c(
        sprintf("df counts coefficients across %d non-reference logits, excluding their intercepts,", q),
        "plus 1 df for each estimated FP power shared across logits. A linear term",
        sprintf("adds %d df, FP1 adds %d df, and FP2 adds %d df.", q, q + 1L, 2L * q + 2L)
      )
    } else if (!is.null(n_ordinal_intercepts)) {
      note_lines <- c(
        sprintf(
          "df excludes the %d ordinal threshold intercept%s and counts fitted regression",
          n_ordinal_intercepts,
          if (identical(n_ordinal_intercepts, 1L)) "" else "s"
        ),
        "coefficients, plus 1 df for each estimated FP power (FP1 = 2 df, FP2 = 4 df).",
        "A retained catzero or spike-at-zero binary indicator adds 1 df."
      )
    } else {
      note_lines <- c(
        "df counts fitted regression coefficients, excluding the intercept, plus 1 df",
        "for each estimated FP power (FP1 = 2 df, FP2 = 4 df). A retained catzero",
        "or spike-at-zero binary indicator adds 1 df; binary-only SAZ uses 1 df."
      )
    }
    for (ln in note_lines) cat(ln, "\n", sep = "")
  }
}

#' Print a Summary of an `mfp2` Model Fit
#'
#' Renders the structured summary produced by [summary.mfp2()].
#'
#' @param x An object of class \code{"summary.mfp2"}.
#' @param notes Logical. If \code{FALSE}, suppress the explanatory model-df
#'   note. Defaults to the value stored by [summary.mfp2()].
#' @param ... Not used.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso [summary.mfp2()], [print.mfp2()]
#'
#' @export
print.summary.mfp2 <- function(x, notes = x$notes, ...) {
  validate_logical_vector(notes, "notes", allowed_lengths = 1L)

  if (isTRUE(x$multinomial)) {
    digits <- if (!is.null(x$digits)) x$digits else 3L
    rule <- paste(rep("=", 78L), collapse = "")
    dash <- paste(rep("-", 78L), collapse = "")
    cat(rule, "\nMFP Model Summary\n", rule, "\n\n", sep = "")
    if (!is.null(x$call)) {
      cat("Call:\n")
      print(x$call)
      cat("\n")
    }
    mfp2_print_model_header(
      family_string = "multinomial",
      criterion = x$criterion,
      converged = x$converged,
      n = x$n,
      metadata = list(link = "logit"),
      digits = digits,
      response_frequencies = x$response_frequencies
    )
    cat(sprintf(
      "\nReference class: %s | Non-reference logits: %d\n",
      x$reference_class, x$n_logits
    ))
    cat("FP powers: common across logits\n\n")
    cat(dash, "\nSelection Overview\n", dash, "\n", sep = "")
    print.data.frame(x$function_table, row.names = FALSE, right = FALSE)
    cat("\n", dash, "\nCoefficient Tests\n", dash, "\n", sep = "")
    display <- x$coefficients
    display$logit <- paste0(display$outcome, " vs ", display$reference)
    display <- display[, c("logit", "term", "coefficient", "se", "z", "p")]
    display <- format_model_print_table(display, digits)
    print.data.frame(display, row.names = FALSE)
    cat("\n")
    mfp2_format_model_fit_block(
      x$model_fit_values,
      digits = digits,
      heading_printer = function(title) cat(dash, "\n", title, "\n", dash, "\n", sep = ""),
      notes = notes
    )
    cat("\n", rule, "\n", sep = "")
    return(invisible(x))
  }

  digits <- if (!is.null(x$digits)) x$digits else 3L
  width <- 78L
  rule_eq <- paste(rep("=", width), collapse = "")
  rule_dash <- paste(rep("-", width), collapse = "")

  section <- function(title) {
    cat(rule_dash, "\n", title, "\n", rule_dash, "\n", sep = "")
  }

  # --- Banner --------------------------------------------------------------
  cat(rule_eq, "\n", "MFP Model Summary", "\n", rule_eq, "\n\n", sep = "")

  # --- Call ----------------------------------------------------------------
  #
  # The summary opens with a plain "Call:" label rather than the boxed
  # "Model Call" section heading used by print.mfp2(). print()'s style is
  # consistent internally -- every one of its sections uses the same dashed
  # rule -- but for summary() the opening call reads better as a light
  # preamble than as a heavy first section.
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  # --- Model metadata ------------------------------------------------------
  model_metadata <- list(
    distribution = x$distribution,
    scale = x$scale,
    scale_fixed = x$scale_fixed,
    dispersion = x$dispersion,
    dispersion_fixed = x$dispersion_fixed,
    scale_strata = x$scale_strata,
    distribution_parameters = x$distribution_parameters,
    censoring = x$censoring,
    link = x$link,
    theta = x$theta
  )
  mfp2_print_model_header(
    family_string = x$family,
    criterion = x$criterion,
    converged = x$converged,
    n = x$n,
    nevents = x$nevents,
    metadata = model_metadata,
    digits = digits,
    response_frequencies = x$response_frequencies
  )
  cat("\n")

  # --- Selection overview --------------------------------------------------
  section("Selection Overview")
  ft <- x$function_table
  # Selected first, then excluded.
  ft <- ft[order(ft$Selected != "yes"), , drop = FALSE]
  print.data.frame(ft, row.names = FALSE, right = FALSE)
  cat(sprintf("\nVariables selected: %d of %d\n\n",
              sum(ft$Selected == "yes"), nrow(ft)))

  # --- Ordinal intercepts --------------------------------------------------
  # Printed as a dedicated section before linear/nonlinear terms so that
  # threshold parameters are not confused with MFP-selected predictor terms.
  if (!is.null(x$ordinal_intercepts) && nrow(x$ordinal_intercepts) > 0L) {
    section("Ordinal Intercepts")
    oi  <- x$ordinal_intercepts
    fmt_ordinal <- function(v) format_print_decimal(v, digits)
    disp_oi <- data.frame(
      Threshold       = oi$threshold,
      Estimate        = fmt_ordinal(oi$estimate),
      `Std. Error`    = fmt_ordinal(oi$se),
      z               = fmt_ordinal(oi$z),
      p               = mfp2_summary_format_p(oi$p, digits),
      `[95% CI]`      = sprintf(
        "[%s, %s]",
        fmt_ordinal(oi$ci_lower),
        fmt_ordinal(oi$ci_upper)
      ),
      check.names     = FALSE,
      stringsAsFactors = FALSE
    )
    print.data.frame(disp_oi, row.names = FALSE, right = FALSE)
    link_label <- switch(
      x$link,
      logistic = "logistic (proportional-odds)",
      probit = "probit",
      loglog = "log-log",
      cloglog = "complementary log-log",
      cauchit = "cauchit",
      "ordinal"
    )
    cat(sprintf(
      "\nThresholds for the %s cumulative-link model.\n\n",
      link_label
    ))
  }

  # --- Linear terms --------------------------------------------------------
  section("Linear Terms")
  if (nrow(x$linear_terms) == 0L) {
    cat("(none)\n\n")
  } else {
    lt <- x$linear_terms
    fmt <- function(v, d = digits) format_print_decimal(v, d)

    disp <- data.frame(
      # Show the user-facing variable name, not the internal fitted-column
      # name. `lt$variable` is populated by mfp2_summary_linear_table() with
      # the ".N" suffix stripped from `lt$term`; the stored `term` remains
      # available for programmatic mapping back to coef(fit).
      Term = lt$variable,
      coef = fmt(lt$coef),
      `se(coef)` = fmt(lt$se),
      stat = fmt(lt$statistic),
      p = mfp2_summary_format_p(lt$p, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    stat_name <- attr(lt, "statistic_label", exact = TRUE)
    if (is.null(stat_name) || !stat_name %in% c("t", "z")) {
      stat_name <- "Statistic"
    }
    names(disp)[names(disp) == "stat"] <- stat_name

    if (!is.null(lt$exp_coef)) {
      disp[["exp(coef)"]] <- fmt(lt$exp_coef)
      disp[["[95% CI]"]] <- sprintf("[%s, %s]", fmt(lt$ci_lower), fmt(lt$ci_upper))
    } else {
      disp[["[95% CI]"]] <- sprintf("[%s, %s]", fmt(lt$ci_lower), fmt(lt$ci_upper))
    }

    print.data.frame(disp, row.names = FALSE, right = FALSE)
    cat("\n")

    if (!is.null(lt$exp_coef)) {
      cat(sprintf(
        "exp(coef) is the %s.\n\n",
        if (is.null(x$exp_label)) "multiplicative effect" else x$exp_label
      ))
    } else {
      cat("Coefficients are on the link scale.\n\n")
    }
  }

  # --- Nonlinear terms -----------------------------------------------------
  section("Nonlinear Terms")
  if (nrow(x$nonlinear_terms) == 0L) {
    cat("(none)\n\n")
  } else {
    nt <- x$nonlinear_terms
    disp <- data.frame(
      Variable = nt$variable,
      Function = nt$form,
      df = nt$df,
      `LR chi-sq` = format_print_decimal(nt$lr_chisq, digits),
      p = mfp2_summary_format_p(nt$p, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    print.data.frame(disp, row.names = FALSE, right = FALSE)
    cat("\n")
    cat(
      "Joint likelihood-ratio tests for each variable in the final MFP model,\n",
      "with all other selected functional forms held fixed. df follow the MFP\n",
      "convention (FP1 = 2, FP2 = 4).\n\n",
      sep = ""
    )

    # Optional: fitted-function formulas.
    if (!is.null(x$formulas)) {
      cat("Fitted functions:\n")
      for (f in x$formulas) cat("  ", f, "\n", sep = "")
      cat("\n")
    }

    # Optional: raw basis coefficients.
    if (!is.null(x$basis)) {
      cat("Basis coefficients:\n")
      bt <- x$basis
      disp_b <- data.frame(
        Variable = bt$variable,
        Term = bt$term,
        coef = format_print_decimal(bt$coef, digits),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      print.data.frame(disp_b, row.names = FALSE, right = FALSE)
      cat("\n")
    }

    # ACD component labels use the compact A(x) notation above. Print the
    # stored transformation once so the internal training scale remains
    # explicit without incorrectly scaling the ordinary FP basis.
    if (!is.null(x$acd_definitions)) {
      cat("ACD definitions:\n")
      for (definition in x$acd_definitions) {
        cat("  ", definition, "\n", sep = "")
      }
      cat("\n")
    }
  }

  # --- Model fit -----------------------------------------------------------
  #
  # The Model Fit block is rendered by the shared helper
  # mfp2_format_model_fit_block(), which is called identically by
  # print.mfp2(). This ensures the two methods display the same family-specific
  # fit statistic and df convention, together with the same optional note. The
  # `section()` helper defined above draws the dash-rule heading used elsewhere
  # in the summary output.
  mfp2_format_model_fit_block(
    values          = x$fit$model_fit,
    digits          = digits,
    heading_printer = section,
    notes           = notes
  )

  cat("\n", rule_eq, "\n", sep = "")
  invisible(x)
}

#' Format a P-Value Vector for Table Display
#'
#' Thin wrapper around `format_print_pvalue()` used by the summary tables so
#' that all p-value formatting flows through a single named helper.
#'
#' @param p Numeric vector of p values.
#' @param digits Integer scalar controlling the number of significant
#'   digits.
#'
#' @return Character vector of formatted p values.
#'
#' @keywords internal
#' @noRd
mfp2_summary_format_p <- function(p, digits = 3L) {
  format_print_pvalue(p, digits)
}
