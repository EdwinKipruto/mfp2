#' Survival Family Specifications
#'
#' `survreg_family()` specifies a parametric accelerated failure-time or
#' location-scale model. `finegray_family()` specifies a proportional
#' subdistribution-hazards model for competing risks.
#'
#' For parametric survival models, repeated candidate models are fitted with
#' [survival::survreg.fit()], and the selected model is refitted with
#' [survival::survreg()].
#' For log-time distributions, an interval-censored observation with lower
#' boundary zero is converted to the mathematically equivalent left-censored
#' observation at its finite upper boundary before applying the logarithmic
#' time transformation. This includes the log-logistic distribution.
#'
#' For Fine--Gray models, [survival::finegray()] expands the multi-state
#' response once. Repeated candidate models then use
#' [survival::agreg.fit()], while the selected model is returned as a weighted
#' robust `coxph` fit. The public model interface remains
#' `family = finegray_family(etype = ...)`. An `id` argument is required by
#' [mfp2()] only for a start--stop multi-state response. Fine--Gray models are
#' not supported by [mfpi()].
#'
#' @param dist A distribution accepted by [survival::survreg()].
#' @param scale Fixed scale. The default, zero, estimates the scale.
#' @param parms Optional distribution parameters accepted by
#'   [survival::survreg()].
#' @param etype Endpoint of interest, using the same interpretation as the
#'   `etype` argument of [survival::finegray()]. `NULL` uses the first event
#'   state, matching `finegray()`.
#' @param timefix Logical; pass event times through survival's round-off check.
#'
#' @return An `mfp2_survreg_family` specification for the `family` argument of
#'   [mfp2()] or [mfpi()], or an `mfp2_finegray_family` specification for
#'   [mfp2()].
#' @examples
#' survreg_spec <- survreg_family(dist = "weibull")
#' finegray_spec <- finegray_family(etype = "cause1")
#'
#' \donttest{
#' # Parametric Weibull survival model.
#' data("gbsg")
#' fit_survreg <- mfp2(
#'   survival::Surv(rectime, censrec) ~ age + nodes,
#'   data = gbsg,
#'   family = survreg_spec,
#'   df = 1,
#'   select = 1,
#'   cycles = 1,
#'   verbose = FALSE
#' )
#'
#' # Fine--Gray model for the first of two competing event types.
#' set.seed(42)
#' n <- 160
#' x <- stats::rnorm(n)
#' cause1_time <- stats::rexp(n, rate = exp(0.3 * x) / 8)
#' cause2_time <- stats::rexp(n, rate = exp(-0.2 * x) / 10)
#' censor_time <- stats::rexp(n, rate = 1 / 15)
#' event <- ifelse(
#'   censor_time < pmin(cause1_time, cause2_time), 0L,
#'   ifelse(cause1_time <= cause2_time, 1L, 2L)
#' )
#' competing_data <- data.frame(
#'   time = pmin(cause1_time, cause2_time, censor_time),
#'   event = factor(
#'     event,
#'     levels = 0:2,
#'     labels = c("censor", "cause1", "cause2")
#'   ),
#'   x = x
#' )
#' fit_finegray <- mfp2(
#'   survival::Surv(time, event) ~ x,
#'   data = competing_data,
#'   family = finegray_spec,
#'   df = 1,
#'   select = 1,
#'   cycles = 1,
#'   verbose = FALSE
#' )
#' }
#' @seealso [survival::survreg()], [survival::finegray()], [mfp2()], [mfpi()]
#' @rdname survival_families
#' @export
survreg_family <- function(dist = "weibull", scale = 0, parms = NULL) {
  scale_supplied <- !missing(scale)
  if (!(is.character(dist) && length(dist) == 1L && !is.na(dist) && nzchar(dist)) &&
      !is.list(dist)) {
    stop(
      "! `dist` must be one distribution name or a distribution list accepted by `survival::survreg()`.",
      call. = FALSE
    )
  }

  # Fail fast on an unrecognized built-in distribution name rather than deferring
  # the error to fit time. A `dist` supplied as a list is treated as a custom
  # distribution object and passed through unchecked; custom distributions
  # registered into `survival::survreg.distributions` at runtime are likewise
  # not validated here.
  if (is.character(dist)) {
    valid_dists <- names(survival::survreg.distributions)
    if (!dist %in% valid_dists) {
      stop(
        "! `dist` = \"", dist, "\" is not a recognized `survreg` distribution. ",
        "Built-in choices are: ", paste(valid_dists, collapse = ", "), ". ",
        "A custom distribution can be supplied as a list instead.",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(scale) || length(scale) != 1L || anyNA(scale) ||
      !is.finite(scale) || scale < 0) {
    stop("! `scale` must be one finite non-negative number.", call. = FALSE)
  }

  if (!is.null(parms)) {
    parms_flat <- unlist(parms, use.names = TRUE)
    if (!is.numeric(parms_flat) || length(parms_flat) < 1L ||
        anyNA(parms_flat) || any(!is.finite(parms_flat))) {
      stop("! `parms` must be `NULL` or a finite numeric value/vector.", call. = FALSE)
    }
  }

  structure(
    list(
      family = "survreg",
      dist = dist,
      scale = unname(scale),
      scale_supplied = scale_supplied,
      parms = parms,
      prepared = NULL
    ),
    class = c("mfp2_survreg_family", "mfp2_family")
  )
}


#' @param strata_action A character string controlling how `strata()` terms in
#'   the model formula are used in Fine--Gray models. Three options are
#'   available, following Zhou et al. (2011, *Biometrics*):
#'   \describe{
#'     \item{`"both"`}{(default) Strata stratify **both** the censoring
#'       distribution (IPCW weights estimated within each stratum) **and** the
#'       baseline subdistribution hazard (each stratum gets its own baseline
#'       hazard in the weighted Cox fit). This is the standard stratified
#'       Fine--Gray model of Zhou et al. (2011).}
#'     \item{`"censoring"`}{Strata stratify the censoring distribution only.
#'       The baseline subdistribution hazard is common across all strata.}
#'     \item{`"baseline"`}{Strata stratify the baseline subdistribution hazard
#'       only. The censoring distribution is estimated pooling across strata.
#'       Equivalent to `ctype = 2` in Zhou et al. (2011).}
#'   }
#'
#' @references
#' Zhou B, Fine J, Laird G (2011). Competing risks regression for stratified
#' data. *Biometrics*, **67**(2), 661--670.
#' \doi{10.1111/j.1541-0420.2010.01493.x}
#'
#' @rdname survival_families
#' @export
finegray_family <- function(etype = NULL, timefix = TRUE,
                            strata_action = c("both", "censoring", "baseline")) {
  if (!is.null(etype)) {
    if (!is.atomic(etype) || is.matrix(etype) || !is.null(dim(etype)) ||
        length(etype) < 1L || anyNA(etype)) {
      stop("! `etype` must be `NULL` or a non-missing event-state value.", call. = FALSE)
    }
    if (length(etype) > 1L) {
      warning("Only the first value of `etype` is used.", call. = FALSE)
      etype <- etype[1L]
    }
  }
  if (!is.logical(timefix) || length(timefix) != 1L || is.na(timefix)) {
    stop("! `timefix` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  strata_action <- match.arg(strata_action)

  structure(
    list(
      family = "finegray",
      etype = etype,
      timefix = timefix,
      strata_action = strata_action,
      prepared = NULL
    ),
    class = c("mfp2_finegray_family", "mfp2_family")
  )
}


#' Multinomial Logistic Family Specification
#'
#' Specifies an unpenalized baseline-category multinomial logistic model for
#' [mfp2()] or [mfpi()]. Fractional-polynomial powers are selected once per
#' predictor and shared by all non-reference logits; regression coefficients
#' remain outcome specific. Candidate fits use [nnet::nnet.default()] and the
#' retained fit uses [nnet::multinom()]. The response normalization, effective
#' case weights, and any class-offset contrasts are prepared once and reused by
#' every candidate fit.
#'
#' @param reference Optional response class or count-matrix column used as the
#'   reference outcome. `NULL` uses the first factor level, the first level
#'   obtained by converting a class-label vector to a factor, or the first
#'   count-matrix column. For a count matrix, a numeric value may instead be
#'   the one-based column index.
#'
#' @details
#' The response may be a factor or character/numeric class-label vector with at
#' least three classes, or a numeric matrix with at least three columns of
#' nonnegative integer class counts. Non-factor label vectors are converted to
#' factors. Every count-response row must contain at least one trial, and every
#' response class must have positive weighted support.
#'
#' If there are \eqn{C} classes and \eqn{Q = C - 1} non-reference logits, a
#' linear term contributes \eqn{Q} regression degrees of freedom. An FP of
#' degree \eqn{m} contributes \eqn{Qm + m}: \eqn{Qm} logit-specific
#' coefficients and \eqn{m} shared-power search degrees of freedom.
#'
#' @return An `mfp2_multinomial_family` specification.
#' @examples
#' set.seed(43)
#' n <- 150
#' multinomial_data <- data.frame(
#'   x1 = stats::rnorm(n),
#'   x2 = stats::runif(n, -1, 1)
#' )
#' eta_b <- with(multinomial_data, 0.4 * x1 - 0.2 * x2)
#' eta_c <- with(multinomial_data, -0.3 * x1 + 0.5 * x2)
#' denominator <- 1 + exp(eta_b) + exp(eta_c)
#' probabilities <- cbind(
#'   A = 1 / denominator,
#'   B = exp(eta_b) / denominator,
#'   C = exp(eta_c) / denominator
#' )
#' multinomial_data$y <- factor(vapply(
#'   seq_len(n),
#'   function(i) sample(c("A", "B", "C"), 1, prob = probabilities[i, ]),
#'   character(1)
#' ))
#'
#' fit_multinomial <- mfp2(
#'   y ~ x1 + x2,
#'   data = multinomial_data,
#'   family = multinomial_family(reference = "A"),
#'   df = 1,
#'   select = 1,
#'   cycles = 1,
#'   verbose = FALSE
#' )
#' predict(
#'   fit_multinomial,
#'   newdata = multinomial_data[1:3, ],
#'   type = "response",
#'   se.fit = FALSE
#' )
#' @seealso [nnet::multinom()], [mfp2()], [mfpi()]
#' @export
multinomial_family <- function(reference = NULL) {
  if (!is.null(reference) &&
      (length(reference) != 1L || is.na(reference) ||
       !(is.character(reference) || is.numeric(reference)))) {
    stop(
      "! `reference` must be `NULL` or one non-missing response level.",
      call. = FALSE
    )
  }

  structure(
    list(family = "multinomial", reference = reference, prepared = NULL),
    class = c("mfp2_multinomial_family", "mfp2_family")
  )
}


#' Ordinal (Proportional-Odds) Family Specification
#'
#' Specifies a proportional-odds ordinal regression model for [mfp2()] or
#' [mfpi()], fitted with the ordinal regression model engine from the `rms`
#' package. Candidate fits during the FP search call [rms::orm.fit()] directly
#' on the model matrix and integer-coded response for speed; the retained model
#' is a native `rms::orm` object, including when variable selection removes all
#' predictors. `rms` is an optional dependency and must be installed to use
#' this family.
#'
#' @param link Character string selecting the cumulative-link function passed to
#'   [rms::orm.fit()]. One of `"logistic"` (default, the proportional-odds
#'   model), `"probit"`, `"loglog"`, `"cloglog"`, or `"cauchit"`.
#'
#' @details
#' The model has one regression coefficient per predictor (a common slope across
#' all cut-points, the proportional-odds assumption) and \eqn{k - 1} intercepts
#' for a response with \eqn{k} distinct ordered values. A linear term therefore
#' contributes one regression degree of freedom and an FP of degree \eqn{m}
#' contributes \eqn{m}; the \eqn{k - 1} intercepts are common to every candidate
#' and are counted in the model degrees of freedom used for AIC and BIC.
#'
#' The response may be an ordered factor (its level order is used), or a numeric,
#' integer, character, or unordered factor. In the latter cases the category
#' order is taken as `sort(unique(y))` -- ascending for numeric and alphabetical
#' for character or unordered factors -- matching the behaviour of
#' [rms::orm()]. Because an alphabetical order may not be the intended clinical
#' order, [mfp2()] emits an informational message stating the assumed order when
#' it is inferred from an unordered categorical response; supply an ordered
#' factor or numeric codes to control the direction. The model is parameterized
#' on the \eqn{P(Y \ge j)} scale, so a positive coefficient means larger
#' predictor values are associated with higher response categories.
#' Case weights are not currently supported for ordinal models; omitted weights
#' or an explicit all-ones vector are accepted.
#'
#' @return An `mfp2_ordinal_family` specification.
#' @examples
#' ordinal_family()
#' ordinal_family(link = "probit")
#'
#' \donttest{
#' if (requireNamespace("rms", quietly = TRUE)) {
#'   set.seed(44)
#'   n <- 180
#'   ordinal_data <- data.frame(
#'     x1 = stats::rnorm(n),
#'     x2 = stats::rbinom(n, 1, 0.5)
#'   )
#'   latent <- with(ordinal_data, 0.6 * x1 + 0.5 * x2 + stats::rlogis(n))
#'   ordinal_data$y <- ordered(
#'     cut(
#'       latent,
#'       breaks = stats::quantile(latent, probs = seq(0, 1, length.out = 5)),
#'       include.lowest = TRUE,
#'       labels = c("low", "medium", "high", "very high")
#'     )
#'   )
#'
#'   fit_ordinal <- mfp2(
#'     y ~ x1 + x2,
#'     data = ordinal_data,
#'     family = ordinal_family(link = "logistic"),
#'     df = 1,
#'     select = 1,
#'     cycles = 1,
#'     verbose = FALSE
#'   )
#'   predict(
#'     fit_ordinal,
#'     newdata = ordinal_data[1:3, ],
#'     type = "response"
#'   )
#' }
#' }
#' @seealso [rms::orm()], [mfp2()], [mfpi()]
#' @export
ordinal_family <- function(link = c("logistic", "probit", "loglog",
                                    "cloglog", "cauchit")) {
  link <- match.arg(link)

  structure(
    list(family = "ordinal", link = link, prepared = NULL),
    class = c("mfp2_ordinal_family", "mfp2_family")
  )
}


#' Does This Family Fit an Ordinary Intercept?
#'
#' Predicate used by the MFP/MFPI engine to decide whether the fitted-model
#' design matrix should carry an intercept column. Cox and Fine--Gray models
#' condition on the baseline hazard and therefore have no ordinary intercept;
#' every other supported family does.
#'
#' @param family_string Character scalar giving the canonical family name, one
#'   of `"gaussian"`, `"binomial"`, `"poisson"`, `"Gamma"`,
#'   `"inverse.gaussian"`, `"negbin"`, `"multinomial"`, `"ordinal"`, `"cox"`,
#'   `"survreg"`, or `"finegray"`.
#'
#' @return `TRUE` when the family fits an ordinary intercept and `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_has_intercept <- function(family_string) {
  !family_string %in% c("cox", "finegray")
}

#' Does the Candidate-Model Design Matrix Need an Explicit Intercept Column?
#'
#' The repeated candidate-model fits used by MFP need an explicit all-ones
#' intercept column only when the low-level fitter would otherwise estimate an
#' ordinary intercept from `x`. Ordinal models handled by `rms::orm.fit()` are
#' the exception: they estimate threshold intercepts internally, so no explicit
#' column is required in the caller-supplied design matrix.
#'
#' @param family_string Canonical family name, as documented in
#'   `mfp2_family_has_intercept()`.
#'
#' @return `TRUE` if the candidate-model design matrix requires an explicit
#'   intercept column, `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_candidate_matrix_has_intercept <- function(family_string) {
  mfp2_family_has_intercept(family_string) &&
    !identical(family_string, "ordinal")
}

#' Does This Family Report an Event Count Rather Than a Continuous Deviance?
#'
#' Cox and Fine--Gray models are compared by event count when constructing
#' summary tables, since they lack an ordinary log-likelihood based
#' observation count.
#'
#' @param family_string Canonical family name.
#'
#' @return `TRUE` for Cox and Fine--Gray models, `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_uses_event_count <- function(family_string) {
  family_string %in% c("cox", "finegray")
}

#' Is This Family a Proportional-Hazards Model?
#'
#' Cox and Fine--Gray subdistribution hazards are both fitted through the
#' partial likelihood of a proportional-hazards structure.
#'
#' @param family_string Canonical family name.
#'
#' @return `TRUE` for `"cox"` and `"finegray"`; `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_is_ph <- function(family_string) {
  family_string %in% c("cox", "finegray")
}

#' Is This a Survival Family?
#'
#' Survival families are Cox proportional hazards, parametric accelerated
#' failure time via `survival::survreg()`, and Fine--Gray subdistribution
#' hazards.
#'
#' @param family_string Canonical family name.
#'
#' @return `TRUE` for `"cox"`, `"survreg"`, and `"finegray"`; `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_is_survival <- function(family_string) {
  family_string %in% c("cox", "survreg", "finegray")
}

#' Is This a Likelihood GLM Family?
#'
#' Predicate distinguishing the six standard likelihood GLM families accepted
#' by `mfp2()` from multinomial, ordinal, and survival families that require
#' specialised fitters.
#'
#' @param family_string Canonical family name.
#'
#' @return `TRUE` when `family_string` is one of `"gaussian"`, `"binomial"`,
#'   `"poisson"`, `"Gamma"`, `"inverse.gaussian"`, or `"negbin"`; `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_is_glm <- function(family_string) {
  family_string %in% c(
    "gaussian", "binomial", "poisson", "Gamma",
    "inverse.gaussian", "negbin"
  )
}

#' Is This the Multinomial Family?
#'
#' @param family_string Canonical family name.
#'
#' @return `TRUE` when `family_string` is exactly `"multinomial"`, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_is_multinomial <- function(family_string) {
  identical(family_string, "multinomial")
}

#' Is This the Ordinal (Proportional-Odds) Family?
#'
#' @param family_string Canonical family name.
#'
#' @return `TRUE` when `family_string` is exactly `"ordinal"`, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @noRd
mfp2_family_is_ordinal <- function(family_string) {
  identical(family_string, "ordinal")
}

#' Number of Logits Fitted by the Family
#'
#' Reports the number of logit predictors estimated by the family. Every
#' non-multinomial family fits a single linear predictor. Multinomial models
#' fit one predictor per non-reference class; that number is stored in the
#' prepared-family metadata during model setup.
#'
#' @param family Prepared family object.
#' @param family_string Canonical family name.
#'
#' @return Integer scalar giving the number of fitted logit predictors.
#'
#' @keywords internal
#' @noRd
mfp2_family_n_logits <- function(family, family_string) {
  if (!mfp2_family_is_multinomial(family_string)) return(1L)
  if (is.null(family$prepared$n_logits)) {
    stop("Internal error: multinomial family has not been prepared.", call. = FALSE)
  }
  as.integer(family$prepared$n_logits)
}

#' Number of Response Classes in a Multinomial Outcome
#'
#' Determines the number of response classes for a candidate multinomial
#' outcome. Class-count matrices are counted directly by the number of
#' columns; label vectors are counted by their number of unique levels after
#' coercion to a factor. Values that do not match any accepted response shape
#' return `NA_integer_`.
#'
#' @param y Response object: a class-count matrix, a factor, or a vector of
#'   character, numeric, or logical class labels.
#'
#' @return Integer scalar with the number of response classes, or
#'   `NA_integer_` for unsupported shapes.
#'
#' @keywords internal
#' @noRd
mfp2_multinomial_n_classes <- function(y) {
  if (is.matrix(y)) return(ncol(y))
  if (is.factor(y)) return(nlevels(y))
  if (is.null(dim(y)) &&
      (is.character(y) || is.numeric(y) || is.logical(y))) {
    return(nlevels(factor(y)))
  }
  NA_integer_
}

#' Does the Family Estimate a Dispersion Parameter?
#'
#' Matches the nuisance-parameter convention of [stats::logLik()] so that
#' MFP information-criterion comparisons treat estimated dispersion as an
#' additional parameter. Recent versions of R let a family declare a fixed
#' dispersion explicitly; older family objects omit that field, in which case
#' the three traditional dispersion families (Gaussian, Gamma, inverse
#' Gaussian) estimate it from the fitted model.
#'
#' @param family Family object.
#' @param family_string Canonical family name.
#'
#' @return `TRUE` when the family estimates a dispersion parameter from the
#'   fitted model, `FALSE` otherwise. Non-GLM and negative-binomial families
#'   always return `FALSE`.
#'
#' @keywords internal
#' @noRd
mfp2_glm_estimates_dispersion <- function(family, family_string) {
  if (!mfp2_family_is_glm(family_string) ||
      identical(family_string, "negbin")) {
    return(FALSE)
  }

  dispersion <- if (is.list(family)) family$dispersion else NULL
  if (!is.null(dispersion)) {
    return(length(dispersion) == 1L && is.na(dispersion))
  }

  family_string %in% c("gaussian", "Gamma", "inverse.gaussian")
}

#' Strip Prepared-Family Metadata Before Returning the Fitted Model
#'
#' Removes the internal caches attached to a prepared family object so that
#' the returned `mfp2` result does not carry MFP-specific fitting metadata.
#' The `mfp2_fit_flags` attribute and any `prepared` list slot are cleared;
#' the remainder of the family is left unchanged.
#'
#' @param family Prepared family object, possibly inheriting from
#'   `"mfp2_family"`.
#'
#' @return The same family object with the `mfp2_fit_flags` attribute and
#'   `prepared` slot removed.
#'
#' @keywords internal
#' @noRd
mfp2_strip_prepared_family <- function(family) {
  attr(family, "mfp2_fit_flags") <- NULL
  if (inherits(family, "mfp2_family")) family$prepared <- NULL
  family
}


#' Prepare a `survreg`-Compatible Family for Repeated MFP Candidate Fits
#'
#' Builds the transformed response expected by [survival::survreg.fit()] once
#' per MFP analysis so it is not recomputed for every candidate model. This
#' function reproduces the non-formula portion of [survival::survreg()],
#' including distribution look-up, response validation, and reparameterisation
#' for censored, left-truncated, and interval-censored data.
#'
#' @param family Family object created by [survreg_family()]; must inherit
#'   from `"mfp2_survreg_family"`.
#' @param y Response, typically a [survival::Surv()] object of one of the
#'   censoring types supported by [survival::survreg()].
#' @param weights Optional numeric vector of strictly positive observation
#'   weights. If `NULL`, all weights are set to `1`.
#' @param strata Optional stratification vector or factor used by
#'   `survival::survreg()` to fit distinct scale parameters per stratum. May be
#'   `NULL` for an unstratified analysis.
#'
#' @return A prepared family object carrying the transformed response,
#'   validated weights, distribution list, and strata definition in a
#'   `prepared` list slot, ready for repeated calls to
#'   [survival::survreg.fit()] during MFP selection.
#'
#' @keywords internal
#' @noRd
prepare_survreg_family <- function(family, y, weights, strata = NULL) {
  if (!inherits(family, "mfp2_survreg_family")) {
    stop("Internal error: invalid survreg family specification.", call. = FALSE)
  }

  n <- NROW(y)
  if (is.null(weights)) weights <- rep.int(1, n)
  if (!is.numeric(weights) || length(weights) != n || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("! `weights` must contain one finite positive value per observation.", call. = FALSE)
  }

  distributions <- survival::survreg.distributions
  requested_dist <- family$dist
  if (is.character(requested_dist)) {
    dist_name <- tryCatch(
      match.arg(requested_dist, names(distributions)),
      error = function(e) {
        stop("Invalid `survreg` distribution: ", conditionMessage(e), call. = FALSE)
      }
    )
    dlist <- distributions[[dist_name]]
  } else {
    dist_name <- requested_dist
    dlist <- requested_dist
  }

  distribution_test <- utils::getFromNamespace("survregDtest", "survival")
  if (!is.list(dlist) || !isTRUE(distribution_test(dlist))) {
    stop("! Invalid distribution object supplied to `survreg_family()`.", call. = FALSE)
  }

  y_fit <- y
  type <- attr(y_fit, "type", exact = TRUE)
  logcorrect <- 0

  # For log-time distributions, an interval (0, upper] is mathematically left
  # censored at its finite upper boundary. Convert it before applying log() so
  # log(0) never enters the transformed interval. Keep this explicit local
  # allow-list because upstream survival releases have historically misspelled
  # the Log logistic distribution name in the corresponding boundary branch.
  log_time_names <- c(
    "Weibull", "Exponential", "Rayleigh", "Log Normal",
    "Log logistic"
  )
  if (identical(type, "interval") && dlist$name %in% log_time_names) {
    fix <- y_fit[, 1L] == 0 & y_fit[, 3L] == 3
    if (any(fix)) {
      y_fit[fix, ] <- cbind(y_fit[fix, 2L], 1, 2)
    }
  }

  if (!is.null(dlist$trans)) {
    exact <- y_fit[, NCOL(y_fit)] == 1
    if (any(exact)) {
      logcorrect <- sum(
        weights[exact] * log(dlist$dtrans(y_fit[exact, 1L]))
      )
    }

    if (identical(type, "interval")) {
      if (any(y_fit[, 3L] == 3)) {
        y_fit <- cbind(dlist$trans(y_fit[, 1:2, drop = FALSE]), y_fit[, 3L])
      } else {
        y_fit <- cbind(dlist$trans(y_fit[, 1L]), y_fit[, 3L])
      }
    } else if (identical(type, "left")) {
      y_fit <- cbind(dlist$trans(y_fit[, 1L]), 2 - y_fit[, 2L])
    } else {
      y_fit <- cbind(dlist$trans(y_fit[, 1L]), y_fit[, 2L])
    }

    if (!all(is.finite(y_fit))) {
      stop("! Invalid survival times for the requested `survreg` distribution.", call. = FALSE)
    }
  } else {
    if (identical(type, "left")) {
      y_fit[, 2L] <- 2 - y_fit[, 2L]
    } else if (identical(type, "interval") && all(y_fit[, 3L] < 3)) {
      y_fit <- y_fit[, c(1L, 3L), drop = FALSE]
    }
  }

  distribution_has_fixed_scale <- !is.null(dlist$scale)
  fit_scale <- family$scale
  if (!is.null(dlist$scale)) {
    if (isTRUE(family$scale_supplied)) {
      warning(
        dlist$name, " has a fixed scale; the value supplied to `survreg_family()` is ignored.",
        call. = FALSE
      )
    }
    fit_scale <- dlist$scale
  }

  if (!is.null(dlist$dist)) {
    dlist <- if (is.atomic(dlist$dist)) {
      distributions[[dlist$dist]]
    } else {
      dlist$dist
    }
  }

  default_parms <- dlist$parms
  requested_parms <- family$parms
  if (is.null(default_parms)) {
    if (!is.null(requested_parms)) {
      stop("! The requested `survreg` distribution has no optional parameters.", call. = FALSE)
    }
    fit_parms <- NULL
  } else {
    if (!is.numeric(default_parms)) {
      stop("Internal error: default survreg distribution parameters are not numeric.", call. = FALSE)
    }
    fit_parms <- default_parms
    if (!is.null(requested_parms)) {
      supplied <- unlist(requested_parms, use.names = TRUE)
      if (is.null(names(supplied)) || any(!nzchar(names(supplied)))) {
        if (length(supplied) != length(fit_parms)) {
          stop(
            "! Unnamed `parms` must have the same length as the distribution defaults.",
            call. = FALSE
          )
        }
        fit_parms[] <- supplied
      } else {
        invalid <- setdiff(names(supplied), names(fit_parms))
        if (length(invalid) > 0L) {
          stop(
            "! Invalid `survreg` parameter name(s): ",
            paste(invalid, collapse = ", "), ".",
            call. = FALSE
          )
        }
        fit_parms[names(supplied)] <- supplied
      }
    }
  }

  if (is.null(strata)) {
    strata_factor <- factor(rep.int(1L, NROW(y_fit)))
  } else {
    strata_factor <- droplevels(normalize_cox_strata(strata, nobs = NROW(y_fit)))
  }
  nstrat <- nlevels(strata_factor)
  if (fit_scale > 0 && nstrat > 1L) {
    stop("! A fixed `survreg` scale cannot be combined with multiple scale strata.", call. = FALSE)
  }

  family$prepared <- list(
    y = y_fit,
    y_original = y,
    dist = dlist,
    dist_requested = requested_dist,
    scale = fit_scale,
    distribution_has_fixed_scale = distribution_has_fixed_scale,
    parms = fit_parms,
    logcorrect = logcorrect,
    strata = as.integer(strata_factor),
    strata_factor = strata_factor,
    strata_supplied = !is.null(strata),
    nstrat = nstrat,
    n_original = NROW(y)
  )
  family
}


#' Prepare a Fine--Gray Subdistribution Family for MFP Candidate Fits
#'
#' Runs [survival::finegray()] once at the beginning of an MFP analysis and
#' retains the input-row mapping needed to expand each dynamically generated
#' fractional-polynomial candidate design matrix to the expanded weighted
#' data set. Repeating this expansion for every candidate model would be both
#' expensive and unnecessary, so the mapping is cached in the prepared
#' family's `prepared` slot.
#'
#' @param family Family object created by [finegray_family()]; must inherit
#'   from `"mfp2_finegray_family"`.
#' @param y Multi-state right-censored response, typically a
#'   [survival::Surv()] object with a factor-valued event column identifying
#'   the competing causes.
#' @param weights Optional numeric vector of strictly positive observation
#'   weights aligned with `y`. If `NULL`, unit weights are used.
#' @param strata Optional stratification vector or factor used by
#'   [survival::finegray()] to stratify estimation of the censoring
#'   distribution.
#' @param id Optional subject identifier vector. Required for a start--stop
#'   multi-state response and optional for ordinary one-row-per-subject data.
#' @param offset Optional numeric offset vector aligned with `y`.
#'
#' @return A prepared Fine--Gray family object carrying the expanded response,
#'   expansion weights, row-index mapping, event-of-interest indicator, and
#'   strata definition in its `prepared` slot.
#'
#' @keywords internal
#' @noRd
prepare_finegray_family <- function(family, y, weights, strata = NULL, id = NULL,
                                    offset = NULL) {
  if (!inherits(family, "mfp2_finegray_family")) {
    stop("Internal error: invalid Fine--Gray family specification.", call. = FALSE)
  }

  n <- NROW(y)
  if (is.null(weights)) weights <- rep.int(1, n)
  if (!is.numeric(weights) || length(weights) != n || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("! `weights` must contain one finite positive value per observation.", call. = FALSE)
  }
  if (is.null(offset)) offset <- rep.int(0, n)
  if (!is.numeric(offset) || length(offset) != n || anyNA(offset) ||
      any(!is.finite(offset))) {
    stop("! `offset` must contain one finite numeric value per observation.", call. = FALSE)
  }
  if (!is.null(id)) {
    if (!is.atomic(id) || is.matrix(id) || !is.null(dim(id))) {
      stop("! `id` must be an atomic vector or factor.", call. = FALSE)
    }
    if (length(id) != n || anyNA(id)) {
      stop("! `id` must contain one non-missing subject identifier per observation.", call. = FALSE)
    }
  }

  type <- attr(y, "type", exact = TRUE)
  if (identical(type, "mcounting") && is.null(id)) {
    stop("! Start--stop Fine--Gray data require the `id` argument.", call. = FALSE)
  }

  original_row <- seq_len(n)
  subject_id <- if (is.null(id)) original_row else id
  d <- data.frame(
    .mfp2_row_id = original_row,
    .mfp2_subject_id = subject_id,
    check.names = FALSE
  )
  d[[".mfp2_response"]] <- y

  strata_action <- if (!is.null(family$strata_action)) family$strata_action else "both"
  strata_factor <- if (is.null(strata)) {
    NULL
  } else {
    normalize_cox_strata(strata, nobs = n)
  }

  rhs <- c(".mfp2_row_id", ".mfp2_subject_id")
  # Pass strata to finegray() for stratified censoring weights when
  # strata_action is "both" or "censoring". When "baseline", the censoring
  # distribution is estimated pooling across strata (Zhou et al. ctype = 2).
  stratify_censoring <- strata_action %in% c("both", "censoring")
  if (!is.null(strata_factor) && stratify_censoring) {
    d[[".mfp2_censor_strata"]] <- strata_factor
    rhs <- c(rhs, "strata(.mfp2_censor_strata)")
  }
  fg_formula <- stats::as.formula(
    paste(".mfp2_response ~", paste(rhs, collapse = " + ")),
    env = environment()
  )

  fg_args <- list(
    formula = fg_formula,
    data = d,
    weights = weights,
    prefix = ".mfp2_fg",
    count = ".mfp2_fg_count",
    timefix = family$timefix
  )
  if (!is.null(family$etype)) fg_args$etype <- family$etype
  if (identical(type, "mcounting")) fg_args$id <- subject_id

  fg <- tryCatch(
    do.call(survival::finegray, fg_args),
    error = function(e) {
      stop("Fine--Gray data preparation failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  required <- c(
    ".mfp2_row_id", ".mfp2_subject_id", ".mfp2_fgstart",
    ".mfp2_fgstop", ".mfp2_fgstatus", ".mfp2_fgwt"
  )
  missing_required <- setdiff(required, names(fg))
  if (length(missing_required) > 0L) {
    stop(
      "Internal error: `finegray()` output is missing: ",
      paste(missing_required, collapse = ", "), ".",
      call. = FALSE
    )
  }

  row_map <- as.integer(fg[[".mfp2_row_id"]])
  y_fg <- survival::Surv(
    fg[[".mfp2_fgstart"]],
    fg[[".mfp2_fgstop"]],
    fg[[".mfp2_fgstatus"]]
  )

  # When strata_action is "both" or "baseline", strata must be passed through
  # to the weighted Cox fit for baseline subdistribution hazard stratification.
  # Expand the original-row strata along the row_map so it aligns with the
  # Fine--Gray pseudo-observations.
  stratify_baseline <- strata_action %in% c("both", "baseline")
  strata_expanded <- if (!is.null(strata_factor) && stratify_baseline) {
    strata_factor[row_map]
  } else {
    NULL
  }

  family$prepared <- list(
    y = y_fg,
    weights = as.numeric(fg[[".mfp2_fgwt"]]),
    row_map = row_map,
    # The offset is invariant across all FP candidates. Expand it once along
    # with the response and Fine--Gray weights instead of indexing the original
    # vector for every weighted Cox fit.
    offset_expanded = as.numeric(offset[row_map]),
    subject_id = fg[[".mfp2_subject_id"]],
    subject_id_original = subject_id,
    # Expanded strata for baseline subdistribution hazard stratification.
    # NULL when strata_action = "censoring" or when no strata are present.
    strata_expanded = strata_expanded,
    strata_expanded_codes = if (!is.null(strata_expanded)) {
      as.integer(strata_expanded)
    } else {
      NULL
    },
    strata_original = if (!is.null(strata_factor) && stratify_baseline) {
      strata_factor
    } else {
      NULL
    },
    event = attr(fg, "event", exact = TRUE),
    nevents = sum(fg[[".mfp2_fgstatus"]] > 0),
    n_original = n,
    n_expanded = nrow(fg)
  )
  family
}


#' Prepare a Multinomial Family for MFP Candidate Fits
#'
#' Normalises the multinomial response, weights, and optional offset to the
#' compact representation expected by the direct baseline-category candidate
#' fitter and by the retained [nnet::multinom()] fit. This includes coercing
#' the response to a factor with the requested reference level, storing the
#' class order, and caching the number of response classes so that later
#' calls do not re-derive it.
#'
#' @param family Family object created by [multinomial_family()]; must
#'   inherit from `"mfp2_multinomial_family"`.
#' @param y Response, one of: a factor, a character or numeric class-label
#'   vector with at least three classes, or a numeric matrix with at least
#'   three nonnegative integer class-count columns.
#' @param weights Optional numeric vector of strictly positive observation
#'   weights. Unit weights are used when `NULL`.
#' @param offset Optional class-offset matrix (`n` by `C`) or
#'   reference-logit offset matrix (`n` by `C - 1`).
#' @param has_offset Logical flag indicating whether an offset was supplied.
#'   Defaults to `!is.null(offset)` and is exposed as an argument so callers
#'   can override the detection.
#'
#' @return A prepared family object carrying the normalised response, class
#'   labels, reference class, weights, offset, and class count in its
#'   `prepared` slot.
#'
#' @keywords internal
#' @noRd
prepare_multinomial_family <- function(family, y, weights, offset = NULL,
                                       has_offset = !is.null(offset)) {
  if (!inherits(family, "mfp2_multinomial_family")) {
    stop("Internal error: invalid multinomial family specification.", call. = FALSE)
  }

  n <- NROW(y)
  if (is.null(weights)) weights <- rep.int(1, n)

  if (!is.matrix(y)) {
    if (!is.factor(y)) {
      response_names <- names(y)
      y <- factor(y)
      names(y) <- response_names
    }
    levels_original <- levels(y)
    reference <- family$reference
    if (is.null(reference)) reference <- levels_original[[1L]]
    reference <- as.character(reference)
    if (!reference %in% levels_original) {
      stop("! Multinomial `reference` is not a response level.", call. = FALSE)
    }
    levels_internal <- c(reference, setdiff(levels_original, reference))
    y_native <- factor(y, levels = levels_internal)
    y_matrix <- nnet::class.ind(y_native)
    colnames(y_matrix) <- levels_internal
    row_totals <- rep.int(1, n)
  } else {
    y_matrix <- unclass(y)
    levels_original <- colnames(y_matrix)
    if (is.null(levels_original) || any(!nzchar(levels_original))) {
      levels_original <- paste0("class", seq_len(ncol(y_matrix)))
      colnames(y_matrix) <- levels_original
    }
    if (anyDuplicated(levels_original)) {
      stop("! Multinomial count-response column names must be unique.", call. = FALSE)
    }
    reference <- family$reference
    if (is.null(reference)) {
      reference_index <- 1L
    } else if (is.numeric(reference)) {
      reference_index <- as.integer(reference)
      if (reference != reference_index || reference_index < 1L ||
          reference_index > ncol(y_matrix)) {
        stop("! Numeric multinomial `reference` is outside the response columns.", call. = FALSE)
      }
    } else {
      reference_index <- match(as.character(reference), levels_original)
      if (is.na(reference_index)) {
        stop("! Multinomial `reference` is not a response column.", call. = FALSE)
      }
    }
    reference <- levels_original[[reference_index]]
    levels_internal <- c(reference, setdiff(levels_original, reference))
    y_matrix <- y_matrix[, levels_internal, drop = FALSE]
    y_native <- y_matrix
    row_totals <- rowSums(y_matrix)
  }

  class_totals <- colSums(y_matrix * as.numeric(weights))
  empty_classes <- names(class_totals)[class_totals <= 0]
  if (length(empty_classes) > 0L) {
    stop(
      "! Every multinomial class must have positive weighted support; empty class(es): ",
      paste(empty_classes, collapse = ", "), ".",
      call. = FALSE
    )
  }

  row_totals <- as.numeric(row_totals)
  case_weights <- as.numeric(weights)

  family$reference <- reference
  family$prepared <- list(
    y = y_native,
    y_matrix = y_matrix,
    # nnet::multinom() performs these two transformations internally. Cache
    # them here because the direct nnet.default() candidate path reuses the
    # same response and case weights for every transformed design.
    y_fit = y_matrix / row_totals,
    case_weights = case_weights,
    effective_weights = case_weights * row_totals,
    levels = levels_internal,
    original_levels = levels_original,
    reference = reference,
    nonreference = levels_internal[-1L],
    n_classes = ncol(y_matrix),
    n_logits = ncol(y_matrix) - 1L,
    row_totals = row_totals,
    effective_n = sum(case_weights * row_totals),
    has_offset = isTRUE(has_offset)
  )
  family$prepared$offset_matrix <- mfp2_prepare_multinomial_offset(
    offset = offset,
    prepared = family$prepared,
    nobs = n,
    has_offset = has_offset
  )
  family
}


#' Normalise a Multinomial Offset for Fitting
#'
#' Normalises a user-supplied multinomial offset once at the public fitting
#' boundary. Both the direct candidate fitter and the retained
#' [nnet::multinom()] fit then share the same class order and
#' reference-contrast representation, so the offset does not have to be
#' re-derived for every candidate model.
#'
#' @param offset Optional class-offset matrix (`n` by `C`) or
#'   reference-logit offset matrix (`n` by `C - 1`).
#' @param prepared Prepared-multinomial metadata list, typically the
#'   `prepared` slot of a family object.
#' @param nobs Integer number of observations.
#' @param has_offset Logical flag indicating whether an offset was supplied.
#'
#' @return An `n` by `C - 1` reference-logit offset matrix, or `NULL` when
#'   `has_offset` is `FALSE`.
#'
#' @keywords internal
#' @noRd
mfp2_prepare_multinomial_offset <- function(offset, prepared, nobs,
                                            has_offset) {
  if (!isTRUE(has_offset)) return(NULL)

  c_classes <- prepared$n_classes
  q_logits <- prepared$n_logits

  if (!is.numeric(offset) || anyNA(offset) || any(!is.finite(offset))) {
    stop("! Multinomial offsets must be finite and numeric.", call. = FALSE)
  }
  if (is.null(dim(offset))) {
    if (q_logits != 1L || length(offset) != nobs) {
      stop(
        "! Multinomial offsets must be an n x C class matrix or an n x (C - 1) reference-logit matrix.",
        call. = FALSE
      )
    }
    offset <- matrix(offset, ncol = 1L)
  }
  if (!is.matrix(offset) || nrow(offset) != nobs ||
      !ncol(offset) %in% c(q_logits, c_classes)) {
    stop(
      "! Multinomial offsets must be an n x C class matrix or an n x (C - 1) reference-logit matrix.",
      call. = FALSE
    )
  }

  if (ncol(offset) == q_logits) {
    logit_offset <- unclass(offset)
    if (!is.null(colnames(logit_offset)) &&
        all(prepared$nonreference %in% colnames(logit_offset))) {
      logit_offset <- logit_offset[, prepared$nonreference, drop = FALSE]
    }
    result <- cbind(0, logit_offset)
    colnames(result) <- prepared$levels
    return(result)
  }

  result <- unclass(offset)
  if (!is.null(colnames(result)) &&
      all(prepared$levels %in% colnames(result))) {
    result <- result[, prepared$levels, drop = FALSE]
  } else {
    original_index <- match(prepared$levels, prepared$original_levels)
    result <- result[, original_index, drop = FALSE]
  }
  result <- sweep(result, 1L, result[, 1L], "-")
  colnames(result) <- prepared$levels
  result
}


#' Multinomial Offset for Prediction
#'
#' Returns a class-logit offset matrix suitable for prediction. Prediction
#' paths need an explicit zero matrix so that class-logit arithmetic can
#' remain branch-free even when no offset was supplied at fitting time.
#' Fitting preparation instead stores `NULL` for no offset, so this helper
#' materialises the zero matrix on demand.
#'
#' @param offset Optional class-offset matrix supplied at prediction, or
#'   `NULL`.
#' @param family Prepared multinomial family object.
#' @param nobs Integer number of prediction rows.
#' @param has_offset Logical flag indicating whether a prediction-time offset
#'   was supplied.
#'
#' @return An `n` by `C - 1` reference-logit offset matrix. When `has_offset`
#'   is `FALSE`, the matrix is filled with zeros.
#'
#' @keywords internal
#' @noRd
mfp2_multinomial_offset <- function(offset, family, nobs, has_offset) {
  prepared <- family$prepared
  if (!isTRUE(has_offset)) {
    return(matrix(
      0,
      nrow = nobs,
      ncol = prepared$n_classes,
      dimnames = list(NULL, prepared$levels)
    ))
  }
  mfp2_prepare_multinomial_offset(
    offset = offset,
    prepared = prepared,
    nobs = nobs,
    has_offset = TRUE
  )
}


# Encode the ordinal response into integer codes 1..k once, so every candidate
# fit reuses the same compact representation. The category order follows
# rms::orm(): an ordered factor uses its level order; a numeric response is
# ordered ascending; a character or unordered factor is ordered by
# sort(unique()), i.e. alphabetically. When the order is inferred from an
# unordered categorical, an informational message names the assumed order so a
# wrong (alphabetical) ordering is catchable.
#' Prepare an Ordinal (Proportional-Odds) Family for MFP Candidate Fits
#'
#' Coerces the ordinal response to the compact representation used by both
#' the MFP candidate-model path and the retained `rms::orm()` fit. The
#' category order follows [rms::orm()]: ordered factors keep their level
#' order; numeric responses are ordered ascending; character or unordered
#' factor responses are ordered by `sort(unique())` (alphabetically). When
#' the order has to be inferred from an unordered categorical response, an
#' informational message names the assumed order so that an unintended
#' alphabetical ordering can be spotted and corrected.
#'
#' @param family Family object created by [ordinal_family()]; must inherit
#'   from `"mfp2_ordinal_family"`.
#' @param y Ordinal response. An ordered factor is recommended so that the
#'   category order is explicit.
#' @param weights Optional numeric vector of observation weights. Ordinal
#'   models currently accept only omitted weights or an explicit all-ones
#'   vector.
#' @param offset Optional numeric offset vector.
#' @param has_offset Logical flag indicating whether an offset was supplied.
#'
#' @return A prepared family object carrying the ordered response, category
#'   labels, weights, and offset in its `prepared` slot.
#'
#' @keywords internal
#' @noRd
prepare_ordinal_family <- function(family, y, weights, offset = NULL,
                                   has_offset = !is.null(offset)) {
  if (!inherits(family, "mfp2_ordinal_family")) {
    stop("Internal error: invalid ordinal family specification.", call. = FALSE)
  }
  if (!requireNamespace("rms", quietly = TRUE)) {
    stop(
      "! `family = ordinal_family()` requires the `rms` package. ",
      "Install it with install.packages(\"rms\").",
      call. = FALSE
    )
  }

  # Ordinal likelihoods are currently unweighted. Validate this once while the
  # response is prepared instead of scanning the invariant weight vector for
  # every FP candidate fit.
  if (!is.null(weights) && any(as.numeric(weights) != 1)) {
    stop(
      "! `family = ordinal_family()` does not support case weights, ",
      "because `rms::orm()` fits an unweighted ordinal model.",
      call. = FALSE
    )
  }

  if (!is.null(dim(y))) {
    stop("! The ordinal response must be a vector, not a matrix.", call. = FALSE)
  }

  inferred_order <- FALSE
  if (is.ordered(y)) {
    # Matrix-interface callers can supply a factor that still carries levels
    # removed by an earlier subset. orm.fit() bases its intercept count on the
    # observed response values, so retain the same invariant here and recode the
    # response contiguously to 1..k.
    y <- droplevels(y)
    levels_ordered <- levels(y)
    codes <- as.integer(y)
  } else if (is.factor(y)) {
    # Unordered factor: rms orders by sort(unique()) = sorted level labels.
    y <- droplevels(y)
    levels_ordered <- sort(unique(as.character(y)))
    codes <- match(as.character(y), levels_ordered)
    inferred_order <- TRUE
  } else if (is.numeric(y)) {
    levels_ordered <- sort(unique(y))
    codes <- match(y, levels_ordered)
  } else if (is.character(y) || is.logical(y)) {
    levels_ordered <- sort(unique(as.character(y)))
    codes <- match(as.character(y), levels_ordered)
    inferred_order <- TRUE
  } else {
    stop(
      "! The ordinal response must be an ordered factor, factor, numeric, ",
      "integer, character, or logical vector.",
      call. = FALSE
    )
  }

  levels_labels <- as.character(levels_ordered)
  k <- length(levels_labels)

  if (inferred_order) {
    message(
      "Ordinal response order inferred as: ",
      paste(levels_labels, collapse = " < "),
      ". Supply an ordered factor or numeric codes to set a different order."
    )
  }

  family$prepared <- list(
    y = as.integer(codes),          # integer codes 1..k in ascending order
    levels = levels_labels,
    n_classes = k,
    n_intercepts = k - 1L,
    link = family$link,
    inferred_order = inferred_order,
    n_original = length(codes),
    has_offset = isTRUE(has_offset),
    offset = if (isTRUE(has_offset)) as.numeric(offset) else NULL
  )
  family
}


#' Dispatch Family-Specific Preparation Before MFP Selection
#'
#' Top-level dispatcher that routes each supported family to its dedicated
#' `prepare_*` helper. This is called once at the start of an MFP analysis so
#' that the family-specific response transformation, weight handling, offset
#' setup, and any cached metadata are constructed a single time and reused
#' across every candidate model fit.
#'
#' @param family Family object supplied by the user. Its concrete class
#'   determines which specialised preparation routine is invoked.
#' @param family_string Canonical family name used by the MFP engine.
#' @param y Response object.
#' @param weights Optional numeric vector of observation weights.
#' @param strata Optional stratification vector or matrix used by survival
#'   families.
#' @param id Optional subject identifier used by Fine--Gray models.
#' @param offset Optional numeric offset (or class-offset matrix for
#'   multinomial models).
#' @param has_offset Logical flag indicating whether an offset was supplied.
#'   Defaults to `!is.null(offset)`.
#'
#' @return The prepared family object, with any family-specific caches
#'   populated in its `prepared` slot. Families that require no special
#'   preparation are returned unchanged.
#'
#' @keywords internal
#' @noRd
prepare_family_for_fit <- function(family, family_string, y, weights,
                                   strata = NULL, id = NULL, offset = NULL,
                                   has_offset = !is.null(offset)) {
  if (identical(family_string, "survreg")) {
    family <- prepare_survreg_family(family, y, weights, strata)
    return(list(family = family, strata = strata))
  }
  if (identical(family_string, "finegray")) {
    family <- prepare_finegray_family(
      family, y, weights, strata, id, offset = offset
    )
    # Strata handling depends on strata_action (Zhou et al. 2011):
    #   "both"     — strata used for censoring (in finegray) AND baseline hazard
    #   "censoring"— strata used for censoring only; baseline hazard is common
    #   "baseline" — strata used for baseline hazard only; censoring is pooled
    # The expanded strata for the Cox fit are stored in family$prepared by
    # prepare_finegray_family() and consumed by fit_finegray(); return NULL
    # here because fit_finegray reads them from prepared$strata_expanded
    # (not from the strata argument passed through fit_model).
    return(list(family = family, strata = NULL))
  }
  if (identical(family_string, "multinomial")) {
    family <- prepare_multinomial_family(
      family, y, weights, offset = offset, has_offset = has_offset
    )
    return(list(family = family, strata = strata))
  }
  if (identical(family_string, "cox")) {
    if (!is.null(strata)) {
      # Keep the normalized factor for the retained formula fit and cache the
      # low-level integer representation used by every coxph.fit() candidate.
      attr(strata, "mfp2_integer_codes") <- as.integer(strata)
    }
    return(list(family = family, strata = strata))
  }
  if (identical(family_string, "ordinal")) {
    family <- prepare_ordinal_family(
      family, y, weights, offset = offset, has_offset = has_offset
    )
    return(list(family = family, strata = strata))
  }
  if (mfp2_family_is_glm(family_string)) {
    # These family properties are invariant across every candidate fit. Cache
    # them on the resolved family object so fit_glm() does not rediscover them
    # inside the repeated FP search.
    attr(family, "mfp2_fit_flags") <- list(
      is_gaussian = identical(family_string, "gaussian"),
      is_negbin = identical(family_string, "negbin"),
      estimates_dispersion = mfp2_glm_estimates_dispersion(
        family = family,
        family_string = family_string
      )
    )
  }
  list(family = family, strata = strata)
}


#' Normalize and validate a model family
#'
#' @param family Character family name, GLM family function/object, or an mfp2
#'   survival-family specification.
#' @param family_arg Character label used in error messages for function inputs.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{family}}{A resolved fitting-family object or character label.}
#'   \item{\code{family_string}}{The normalized character family name.}
#' }
#'
#' @keywords internal
#' @noRd
normalize_family_argument <- function(family, family_arg = deparse(substitute(family))) {
  allowed_glm_families <- c(
    "gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian"
  )
  allowed_families <- c(
    allowed_glm_families, "negbin", "multinomial", "cox", "survreg",
    "finegray", "ordinal"
  )
  family_arg <- paste(family_arg, collapse = " ")

  if (inherits(family, "mfp2_family")) {
    family_string <- family$family
    if (!is.character(family_string) || length(family_string) != 1L ||
        !family_string %in% c("survreg", "finegray", "multinomial", "ordinal")) {
      stop("! Invalid mfp2 family specification.", call. = FALSE)
    }
    return(list(family = family, family_string = family_string))
  }

  if (is.character(family)) {
    if (length(family) != 1L) {
      stop(
        sprintf(
          "! `family` must be a single character string; got %d values: %s.",
          length(family),
          paste(family, collapse = ", ")
        ),
        "\ni Supported character families are: gaussian, binomial, poisson, Gamma, inverse.gaussian, negbin, multinomial, cox, survreg, finegray.",
        call. = FALSE
      )
    }

    family_key <- tolower(family)
    family <- switch(
      family_key,
      gaussian = "gaussian",
      binomial = "binomial",
      poisson = "poisson",
      gamma = "Gamma",
      inverse.gaussian = "inverse.gaussian",
      negbin = "negbin",
      multinomial = "multinomial",
      cox = "cox",
      survreg = "survreg",
      finegray = "finegray",
      ordinal = "ordinal",
      family
    )

    if (grepl("^quasi", family_key)) {
      stop(
        sprintf("! Family '%s' has no likelihood and is not supported by likelihood-based MFP selection.", family_key),
        call. = FALSE
      )
    }

    if (!family %in% allowed_families) {
      stop(
        sprintf("! Invalid family: '%s'.", family),
        sprintf(
          "\ni Supported character families are: %s.",
          paste(allowed_families, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    family_obj <- if (family %in% c("cox", "negbin")) {
      family
    } else if (identical(family, "multinomial")) {
      multinomial_family()
    } else if (identical(family, "survreg")) {
      survreg_family()
    } else if (identical(family, "finegray")) {
      finegray_family()
    } else if (identical(family, "ordinal")) {
      ordinal_family()
    } else {
      switch(
        family,
        gaussian = stats::gaussian(),
        binomial = stats::binomial(),
        poisson  = stats::poisson(),
        Gamma = stats::Gamma(),
        inverse.gaussian = stats::inverse.gaussian()
      )
    }

    return(list(
      family = family_obj,
      family_string = family
    ))
  }

  if (is.function(family)) {
    family_obj <- tryCatch(
      family(),
      error = function(e) {
        stop(
          sprintf(
            "! Could not create a family object from the provided function `%s`. Error: %s",
            family_arg,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    if (!inherits(family_obj, "family")) {
      stop(
        sprintf(
          "! The provided function `%s` did not return a valid GLM family object.",
          family_arg
        ),
        call. = FALSE
      )
    }

    family_string <- family_obj$family

    if (identical(family_string, "cox")) {
      stop(
        "! `cox` must be specified as a character string, not as a function.",
        call. = FALSE
      )
    }

    if (grepl("^quasi", tolower(family_string))) {
      stop(
        sprintf(
          "! Family returned by `%s` has no likelihood and is not supported by likelihood-based MFP selection.",
          family_arg
        ),
        call. = FALSE
      )
    }

    if (!family_string %in% allowed_glm_families) {
      stop(
        sprintf(
          "! Invalid family returned by `%s`: '%s'. Supported families are: %s.",
          family_arg,
          family_string,
          paste(allowed_glm_families, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(list(
      family = family_obj,
      family_string = family_string
    ))
  }

  if (inherits(family, "family")) {
    family_string <- family$family

    if (identical(family_string, "cox")) {
      stop(
        "! `cox` must be specified as a character string, not as a family object.",
        call. = FALSE
      )
    }

    if (grepl("^quasi", tolower(family_string))) {
      stop(
        "! Quasi families have no likelihood and are not supported by likelihood-based MFP selection.",
        call. = FALSE
      )
    }

    if (!family_string %in% allowed_glm_families) {
      stop(
        sprintf(
          "! Invalid family: '%s'. Supported families are: %s.",
          family_string,
          paste(allowed_glm_families, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(list(
      family = family,
      family_string = family_string
    ))
  }

  stop(
    "! `family` must be a character string, a GLM family function/object, or an mfp2 family specification.",
    call. = FALSE
  )
}

#' Extract family name for the formula interface
#'
#' Internal helper used by `mfp2.formula()` to determine the family name before
#' calling `mfp2.default()`. This is needed because the formula interface must
#' know whether Cox-specific formula terms such as `strata()` are allowed, but
#' family validation and normalization are still handled by `mfp2.default()`.
#'
#' @param family A character family name, a GLM family function, or a GLM family
#'   object.
#'
#' @return A single character string giving the family name.
#'
#' @keywords internal
#' @noRd
get_family_string_formula <- function(family, family_arg = deparse(substitute(family))) {
  normalize_family_argument(
    family,
    family_arg = family_arg
  )$family_string
}


#' Resolve a model family once for repeated internal fits
#'
#' @param family Character family name, family function, family object, or an
#'   mfp2 survival-family specification.
#' @return A resolved GLM family object, character Cox/negbin label, or an mfp2
#'   survival-family specification.
#' @keywords internal
#' @noRd
resolve_fit_model_family <- function(family) {
  normalize_family_argument(family)$family
}

#' Validate response against model family
#'
#' Validates the response object `y` after the family has been normalized.
#' This function checks only response shape and family-specific admissibility;
#' it does not modify `y`.
#'
#' @param y Response vector, matrix, factor, or `survival::Surv()` object.
#' @param family_string Normalized supported family name.
#' @param nobs Expected number of observations.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @keywords internal
#' @noRd
validate_family_response <- function(y, family_string, nobs) {

  if (!is.character(family_string) || length(family_string) != 1L) {
    stop("! `family_string` must be a single character string.", call. = FALSE)
  }

  if (!is.numeric(nobs) || length(nobs) != 1L ||
      anyNA(nobs) || !is.finite(nobs) || nobs < 1L) {
    stop("! `nobs` must be a positive finite scalar.", call. = FALSE)
  }

  if (family_string == "cox") {
    if (!survival::is.Surv(y)) {
      stop(
        "! For `family = 'cox'`, `y` must be a `survival::Surv()` object.",
        call. = FALSE
      )
    }

    if (NROW(y) != nobs) {
      stop(
        paste0(
          "! `y` has ", NROW(y), " rows but `x` has ", nobs,
          " rows; they must match."
        ),
        call. = FALSE
      )
    }

    type <- attr(y, "type", exact = TRUE)
    if (!identical(type, "right")) {
      stop(
        paste0(
          "! Only right-censored survival data are currently supported; ",
          "`y` has censoring type '", type, "'."
        ),
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(as.matrix(y)))) {
      stop(
        "! For `family = 'cox'`, `y` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "survreg") {
    if (!survival::is.Surv(y)) {
      stop(
        "! For `family = survreg_family()`, `y` must be a `survival::Surv()` object.",
        call. = FALSE
      )
    }
    if (NROW(y) != nobs) {
      stop("! `y` and `x` must contain the same number of observations.", call. = FALSE)
    }
    type <- attr(y, "type", exact = TRUE)
    if (!type %in% c("right", "left", "interval", "interval2")) {
      stop(
        sprintf("! `survreg` does not support survival response type '%s'.", type),
        call. = FALSE
      )
    }
    if (anyNA(y) || any(!is.finite(as.matrix(y)))) {
      stop("! The parametric survival response must contain only finite, non-missing values.", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (family_string == "finegray") {
    if (!survival::is.Surv(y)) {
      stop(
        "! For `family = finegray_family()`, `y` must be a multi-state `survival::Surv()` object.",
        call. = FALSE
      )
    }
    if (NROW(y) != nobs) {
      stop("! `y` and `x` must contain the same number of observations.", call. = FALSE)
    }
    type <- attr(y, "type", exact = TRUE)
    if (!type %in% c("mright", "mcounting")) {
      stop(
        sprintf("! Fine--Gray requires a multi-state response; `y` has type '%s'.", type),
        call. = FALSE
      )
    }
    states <- attr(y, "states", exact = TRUE)
    if (length(states) < 2L) {
      stop("! Fine--Gray requires at least two event states.", call. = FALSE)
    }
    if (anyNA(y) || any(!is.finite(as.matrix(y)))) {
      stop("! The Fine--Gray response must contain only finite, non-missing values.", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (survival::is.Surv(y)) {
    stop(
      paste0(
        "! Response is a `survival::Surv()` object but family = '",
        family_string, "'. Select `family = 'cox'`, `survreg_family()`, or ",
        "`finegray_family()` as appropriate."
      ),
      call. = FALSE
    )
  }

  if (NROW(y) != nobs) {
    stop(
      paste0(
        "! `y` has ", NROW(y), " observations but `x` has ", nobs,
        " rows; they must match."
      ),
      call. = FALSE
    )
  }

  if (is.data.frame(y)) {
    stop(
      "! `y` must not be a data frame.",
      call. = FALSE
    )
  }

  if (family_string == "multinomial") {
    is_label_vector <- is.factor(y) ||
      (is.null(dim(y)) &&
       (is.character(y) || is.numeric(y) || is.logical(y)))
    if (is_label_vector) {
      if (anyNA(y)) {
        stop("! Multinomial class-label responses must not contain missing values.", call. = FALSE)
      }
      if (is.numeric(y) && any(!is.finite(y))) {
        stop("! Numeric multinomial class labels must be finite.", call. = FALSE)
      }
      if (nlevels(if (is.factor(y)) y else factor(y)) < 3L) {
        stop(
          "! Multinomial class-label responses must contain at least three classes.",
          call. = FALSE
        )
      }
      return(invisible(TRUE))
    }
    if (!is.matrix(y) || !is.numeric(y) || ncol(y) < 3L) {
      stop(
        "! For `family = \"multinomial\"`, `y` must be a factor or character/numeric class-label vector with at least three classes, or a numeric count matrix with at least three columns.",
        call. = FALSE
      )
    }
    if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) {
      stop("! Multinomial response counts must be finite non-negative integers.", call. = FALSE)
    }
    if (any(rowSums(y) <= 0)) {
      stop("! Every multinomial count-response row must contain at least one trial.", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (family_string == "ordinal") {
    if (!is.null(dim(y))) {
      stop("! For `family = ordinal_family()`, `y` must be a vector, not a matrix.", call. = FALSE)
    }
    if (!(is.factor(y) || is.numeric(y) || is.character(y) || is.logical(y))) {
      stop(
        "! For `family = ordinal_family()`, `y` must be an ordered factor, factor, numeric, integer, character, or logical vector.",
        call. = FALSE
      )
    }
    if (anyNA(y)) {
      stop("! The ordinal response must not contain missing values.", call. = FALSE)
    }
    if (is.numeric(y) && any(!is.finite(y))) {
      stop("! A numeric ordinal response must contain only finite values.", call. = FALSE)
    }
    n_levels <- if (is.factor(y)) nlevels(droplevels(y)) else length(unique(y))
    if (n_levels < 3L) {
      stop(
        "! An ordinal response must have at least three distinct ordered categories; ",
        "use `family = 'binomial'` for a two-level response.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }

  if (is.matrix(y)) {
    if (family_string != "binomial") {
      stop(
        "! Matrix responses are only supported for `family = 'binomial'`.",
        "i For grouped binomial counts, use `y = cbind(successes, failures)`.",
        call. = FALSE
      )
    }

    if (!is.numeric(y) || ncol(y) != 2L) {
      stop(
        "! For `family = 'binomial'`, matrix `y` must be a numeric two-column matrix.",
        "i Use `y = cbind(successes, failures)`.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y < 0) ||
        any(y != floor(y))) {
      stop(
        paste0(
          "! For `family = 'binomial'`, matrix `y` counts must be finite, ",
          "non-missing, non-negative integers."
        ),
        call. = FALSE
      )
    }

    if (any(rowSums(y) <= 0)) {
      stop(
        "! For `family = 'binomial'`, each row of matrix `y` must contain at least one trial.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "gaussian") {
    if (!is.numeric(y)) {
      stop(
        "! For `family = 'gaussian'`, `y` must be a numeric vector.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y))) {
      stop(
        "! For `family = 'gaussian'`, `y` must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string %in% c("Gamma", "inverse.gaussian")) {
    if (!is.numeric(y)) {
      stop(
        sprintf("! For `family = '%s'`, `y` must be a numeric vector.", family_string),
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y <= 0)) {
      stop(
        sprintf(
          "! For `family = '%s'`, `y` must contain only finite, non-missing, strictly positive values.",
          family_string
        ),
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "negbin") {
    if (!is.numeric(y)) {
      stop(
        "! For `family = 'negbin'`, `y` must be a numeric vector of counts.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) {
      stop(
        "! For `family = 'negbin'`, `y` must contain only finite, non-missing, non-negative integer counts.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "poisson") {
    if (!is.numeric(y)) {
      stop(
        "! For `family = 'poisson'`, `y` must be a numeric vector.",
        call. = FALSE
      )
    }

    if (anyNA(y) || any(!is.finite(y)) || any(y < 0) ||
        any(y != floor(y))) {
      stop(
        paste0(
          "! For `family = 'poisson'`, `y` must contain only finite, ",
          "non-missing, non-negative integer counts."
        ),
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  if (family_string == "binomial") {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop(
          "! For `family = 'binomial'`, factor `y` must have exactly two levels.",
          call. = FALSE
        )
      }

      if (anyNA(y)) {
        stop(
          "! For `family = 'binomial'`, factor `y` must not contain missing values.",
          call. = FALSE
        )
      }

      return(invisible(TRUE))
    }

    if (is.numeric(y)) {
      if (anyNA(y) || any(!is.finite(y))) {
        stop(
          "! For `family = 'binomial'`, numeric `y` must contain only finite, non-missing values.",
          call. = FALSE
        )
      }

      if (any(y < 0 | y > 1)) {
        stop(
          "! For `family = 'binomial'`, numeric `y` must contain values in [0, 1].",
          call. = FALSE
        )
      }

      return(invisible(TRUE))
    }

    stop(
      "! For `family = 'binomial'`, `y` must be a numeric vector, two-level factor, or numeric two-column matrix.",
      call. = FALSE
    )
  }

  stop(
    paste0("! Unsupported family: '", family_string, "'."),
    call. = FALSE
  )
}
