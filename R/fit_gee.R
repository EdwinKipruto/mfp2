# =============================================================================
# Generalized Estimating Equations (GEE) support for mfp2
# -----------------------------------------------------------------------------
# GEE candidate and final fits are delegated to the geepack package:
#   * geepack::geese.fit() is used for the many lightweight candidate fits in
#     the fractional-polynomial search (fast = TRUE).
#   * geepack::geeglm()   is used for the single retained final fit
#     (fast = FALSE), so the returned object is a native `geeglm` object that is
#     numerically identical to a direct geeglm() call.
#
# GEE is not a full-likelihood model, so mfp2's three selection criteria are
# mapped onto quasi-likelihood quantities:
#
#   * criterion = "pvalue": differences of negative overall robust Wald
#     chi-square values from separate fits, matching the `mfp: xtgee`
#     selection convention. This is not a single-fit block Wald test.
#   * criterion = "aic": MFP-adjusted QICu = -2Q/phi + 2p, where p includes the
#     searched-power degrees of freedom (based on Pan's QICu form, 2001).
#   * criterion = "bic": QBIC = -2Q/phi + log(n) p with the same adjusted p,
#     where n is the selected cluster or observation count and phi is the
#     dispersion from a shared full reference model.
#
# The quasi-likelihood Q is computed in R directly from the fitted means and
# effective prior weights.  The weights are essential for ordinary weighted
# GEE fits and for grouped-binomial responses, whose binomial totals are stored
# as effective weights.
# =============================================================================

#' Supported GEE response families and their canonical geepack variance names
#'
#' @keywords internal
#' @noRd
mfp2_gee_supported_families <- function() {
  c("gaussian", "binomial", "poisson", "Gamma")
}

#' Is this the GEE meta-family?
#'
#' @param family_string Canonical family name.
#' @return `TRUE` when `family_string` is exactly `"gee"`, `FALSE` otherwise.
#' @keywords internal
#' @noRd
mfp2_family_is_gee <- function(family_string) {
  identical(family_string, "gee")
}

# Keep clusters together without changing the observation order within a
# cluster. The result indexes the fitted sample in its original row order.
# geepack requires this grouping for both geese.fit() and geeglm().
mfp2_gee_cluster_order <- function(id) {
  cluster <- match(id, unique(id))
  order(cluster, seq_along(cluster), method = "radix")
}

#' Robust Wald test for a coefficient block
#'
#' Requires full rank in the tested sandwich covariance block. Rank-deficient
#' blocks have no valid joint Wald statistic and return missing test values.
#' @param compute_pvalue Whether to calculate the chi-square tail probability;
#'   the model-selection score only uses the statistic.
#'
#' @keywords internal
#' @noRd
gee_robust_wald_test <- function(coefficients, robust_vcov, indices,
                                 compute_pvalue = TRUE) {
  indices <- unique(as.integer(indices))
  indices <- indices[is.finite(indices) & indices >= 1L &
                       indices <= length(coefficients)]

  if (length(indices) == 0L) {
    return(list(statistic = 0, pvalue = 1, df = 0L, dev_diff = 0))
  }

  covariance_dim <- dim(robust_vcov)
  if (is.null(robust_vcov) || length(covariance_dim) != 2L ||
      any(covariance_dim < max(indices))) {
    return(list(statistic = NA_real_, pvalue = NA_real_, df = NA_integer_,
                dev_diff = NA_real_))
  }

  beta <- as.numeric(coefficients[indices])
  covariance <- as.matrix(robust_vcov)[indices, indices, drop = FALSE]
  if (anyNA(beta) || any(!is.finite(beta)) || anyNA(covariance) ||
      any(!is.finite(covariance))) {
    return(list(statistic = NA_real_, pvalue = NA_real_, df = NA_integer_,
                dev_diff = NA_real_))
  }

  covariance <- (covariance + t(covariance)) / 2
  eig <- tryCatch(eigen(covariance, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || length(eig$values) == 0L) {
    return(list(statistic = NA_real_, pvalue = NA_real_, df = NA_integer_,
                dev_diff = NA_real_))
  }

  scale <- max(abs(eig$values))
  tolerance <- max(1, length(beta)) * sqrt(.Machine$double.eps) * scale
  estimable <- eig$values > tolerance
  df <- sum(estimable)
  if (df != length(indices)) {
    return(list(statistic = NA_real_, pvalue = NA_real_, df = df,
                dev_diff = NA_real_))
  }

  projected <- crossprod(eig$vectors[, estimable, drop = FALSE], beta)
  statistic <- sum(as.numeric(projected)^2 / eig$values[estimable])
  statistic <- max(0, unname(statistic))

  list(
    statistic = statistic,
    pvalue = if (compute_pvalue) {
      stats::pchisq(statistic, df = df, lower.tail = FALSE)
    } else {
      NA_real_
    },
    df = unname(df),
    dev_diff = statistic
  )
}

#' Apply Stata's finite-cluster correction to a sandwich covariance
#'
#' Stata's clustered robust covariance multiplies the sandwich middle by
#' \eqn{K/(K-1)}, where \eqn{K} is the number of independent clusters. geepack
#' reports the uncorrected asymptotic sandwich covariance.
#'
#' @keywords internal
#' @noRd
gee_robust_vcov <- function(robust_vcov, n_clusters = NULL) {
  if (is.null(n_clusters)) return(robust_vcov)

  n_clusters <- as.numeric(n_clusters)
  if (length(n_clusters) != 1L || !is.finite(n_clusters) || n_clusters <= 1) {
    corrected <- as.matrix(robust_vcov)
    corrected[,] <- NA_real_
    return(corrected)
  }

  robust_vcov * (n_clusters / (n_clusters - 1))
}

#' Select the covariance estimator requested for a GEE fit
#'
#' The component names follow `geepack::vcov.geeglm()`. The ordinary
#' sandwich estimator receives the package's documented finite-cluster
#' correction; jackknife estimators are returned unchanged.
#'
#' @keywords internal
#' @noRd
gee_selected_vcov <- function(geese, std.err = "san.se", n_clusters = NULL) {
  covariance <- switch(
    std.err,
    "jack" = geese$vbeta.ajs,
    "j1s" = geese$vbeta.j1s,
    "fij" = geese$vbeta.fij,
    "san.se" = geese$vbeta,
    stop("Internal error: unsupported GEE covariance estimator.", call. = FALSE)
  )

  if (is.null(covariance)) {
    stop(
      sprintf("Internal error: GEE covariance estimator '%s' was not computed.",
              std.err),
      call. = FALSE
    )
  }

  covariance <- as.matrix(covariance)
  if (identical(std.err, "san.se")) {
    covariance <- gee_robust_vcov(covariance, n_clusters)
  }
  covariance
}

#' Configure geese.fit() to compute the requested covariance estimator
#'
#' This mirrors the flag mapping performed internally by
#' `geepack::geeglm()`.
#'
#' @keywords internal
#' @noRd
gee_control_for_std_err <- function(control, std.err) {
  control$jack <- as.integer(identical(std.err, "jack"))
  control$j1s <- as.integer(identical(std.err, "j1s"))
  control$fij <- as.integer(identical(std.err, "fij"))
  control
}

# Scale the invariant GEE adjustment columns and the shared FP/ACD basis once
# per focal search. Each basis column represents the same n observations in
# every candidate that uses it, so its maximum magnitude is invariant too.
# The candidate loop replaces only the focal columns and their scale factors.
# This working design is private: the C++ candidate copier modifies it in place.
prepare_gee_candidate_scaling <- function(design, basis) {
  intercept <- which(colnames(design) == "(Intercept)")
  result <- prepare_gee_candidate_scaling_cpp(
    design, basis,
    intercept_col = if (length(intercept) == 1L) as.integer(intercept) else 0L
  )
  names(result$column_scales) <- colnames(design)
  names(result$basis_scales) <- colnames(basis)
  result
}

gee_design_column_scales <- function(x) {
  apply(x, 2L, function(col) {
    s <- max(abs(col))
    if (!is.finite(s) || s <= 0) 1 else s
  })
}

#' Refit a Retained geeglm Object with Exact Integer Wave Distances
#'
#' `geeglm()` coerces waves to a factor, compressing absent integer levels (for
#' example, visits 1 and 3 become codes 1 and 2). Construct the native return
#' object with `geeglm()`, then replace its GEE core when exact integer gaps
#' require the lower-level `geese.fit()` call.
#'
#' @keywords internal
#' @noRd
gee_refit_exact_waves <- function(fit, y, id, waves, offset, weights,
                                  control, response_family, corstr,
                                  scale.fix, scale.value) {
  x <- fit$geese$X
  ans <- geepack::geese.fit(
    x = x,
    y = y,
    id = id,
    offset = offset,
    weights = weights,
    waves = waves,
    control = control,
    b = fit$coefficients,
    family = response_family,
    corstr = corstr,
    scale.fix = scale.fix,
    scale.value = scale.value
  )

  ans <- c(
    ans,
    list(call = fit$geese$call, formula = fit$geese$formula)
  )
  class(ans) <- "geese"
  ans$X <- x
  ans$id <- id
  ans$waves <- waves
  ans$weights <- weights

  fit$geese <- ans
  fit$weights <- weights
  fit$coefficients <- ans$beta
  fit$offset <- offset
  fit$linear.predictors <- as.numeric(x %*% ans$beta) + offset
  fit$fitted.values <- response_family$linkinv(fit$linear.predictors)
  fit$modelInfo <- ans$model
  fit$id <- id
  fit$corstr <- ans$model$corstr
  fit$cor.link <- ans$model$cor.link
  fit$control <- ans$control
  fit
}

#' Negative overall robust Wald chi-square used by Stata mfp/xtgee
#'
#' @keywords internal
#' @noRd
gee_selection_deviance <- function(coefficients,
                                         robust_vcov,
                                         n_clusters = NULL) {
  intercept <- match("(Intercept)", names(coefficients), nomatch = 0L)
  tested <- setdiff(seq_along(coefficients), intercept)
  test <- gee_robust_wald_test(
    coefficients,
    gee_robust_vcov(robust_vcov, n_clusters),
    tested,
    compute_pvalue = FALSE
  )
  if (is.finite(test$statistic)) -test$statistic else NA_real_
}

#' Quasi-likelihood for a GEE fit
#'
#' Computes the quasi-likelihood \eqn{Q(\mu; y)} used by the QIC family of
#' criteria.  Contributions are multiplied by the effective prior weights;
#' this also handles grouped-binomial responses, for which the response is a
#' success proportion and the total number of trials is part of the effective
#' weight.  The Gamma expression is the quasi-likelihood obtained by integrating
#' \eqn{(y-\mu)/\mu^2}, namely \eqn{-y/\mu-\log(\mu)}.
#'
#' @param y Numeric response vector. For grouped binomial data this is the
#'   observed proportion of successes, matching the value stored by
#'   \code{geepack::geeglm()}.
#' @param mu Numeric vector of fitted means.
#' @param weights Numeric vector of effective prior weights.
#' @param family_string Inner response family: one of \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}, or \code{"Gamma"}.
#'
#' @return A single numeric value, the quasi-likelihood.
#'
#' @references
#' Pan, W. (2001). Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, 57, 120--125.
#'
#' @keywords internal
#' @noRd
gee_quasi_likelihood <- function(y, mu, family_string, weights = NULL) {
  if (is.null(weights)) weights <- rep.int(1, length(y))
  weights <- as.numeric(weights)

  if (length(y) != length(mu) || length(y) != length(weights)) {
    stop("Internal error: GEE quasi-likelihood inputs must have equal length.",
         call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Internal error: GEE quasi-likelihood weights must be finite and non-negative.",
         call. = FALSE)
  }

  # Zero-weight observations do not contribute. Subset before evaluating the
  # family expression so an otherwise irrelevant boundary value cannot produce
  # 0 * Inf = NaN.
  positive_weight <- weights > 0
  y <- y[positive_weight]
  mu <- mu[positive_weight]
  weights <- weights[positive_weight]

  contribution <- switch(
    family_string,
    gaussian = ((y - mu)^2) / -2,
    poisson = {
      # Use the limiting value 0 * log(0) = 0 for zero counts. This also
      # preserves -Inf for an impossible positive count with fitted mean zero.
      log_term <- numeric(length(y))
      positive_y <- y > 0
      log_term[positive_y] <- y[positive_y] * log(mu[positive_y])
      log_term - mu
    },
    binomial = {
      # y*log(mu) + (1-y)*log(1-mu), evaluated without 0 * Inf at the
      # legitimate y/mu boundaries zero and one.
      success_term <- numeric(length(y))
      failure_term <- numeric(length(y))
      has_success <- y > 0
      has_failure <- y < 1
      success_term[has_success] <- y[has_success] * log(mu[has_success])
      failure_term[has_failure] <-
        (1 - y[has_failure]) * log1p(-mu[has_failure])
      success_term + failure_term
    },
    Gamma = {
      gamma_term <- -y / mu - log(mu)
      gamma_term[mu == 0 & y > 0] <- -Inf
      gamma_term
    },
    stop(
      sprintf(
        "Internal error: unsupported GEE response family '%s'.",
        family_string
      ),
      call. = FALSE
    )
  )
  sum(weights * contribution)
}

#' Normalise a control list for the GEE fitter
#'
#' Resolves the user-supplied control list to a \code{geepack::geese.control()}
#' object. GLM-style names are accepted as aliases so that the same
#' \code{control = list(maxit = ..., epsilon = ...)} works across families:
#' \code{maxit} maps to the maximum number of Fisher-scoring iterations and
#' \code{epsilon} to the convergence tolerance.
#'
#' @param control User-supplied control list, or \code{NULL} for GEE defaults.
#'
#' @return A validated list returned by \code{geepack::geese.control()}.
#' @keywords internal
#' @noRd
normalize_gee_control <- function(control = NULL) {
  if (is.null(control)) {
    return(geepack::geese.control())
  }

  if (!is.list(control)) {
    stop(
      "For GEE models, `control` must be `NULL` or a list accepted by ",
      "`geepack::geese.control()` (e.g. list(maxit = 50, epsilon = 1e-6)).",
      call. = FALSE
    )
  }

  # geese.control() uses the same `maxit`/`epsilon`/`trace` names as
  # glm.control(); pass recognised entries straight through. Unknown entries
  # are rejected by geese.control() itself.
  tryCatch(
    do.call(geepack::geese.control, control),
    error = function(e) {
      stop(
        "Invalid GEE `control`: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

#' Prepare a GEE family once for repeated MFP candidate fits
#'
#' Validates the clustering structure and response, resolves the inner GLM
#' response family, converts grouped-binomial responses to the proportion /
#' cluster-weight representation used by \code{geepack::geeglm()}, and caches
#' every invariant GEE input in \code{family$prepared}. Candidate fits then read
#' the cached objects and never re-validate clusters, re-parse the family, or
#' rebuild the control list.
#'
#' @param family A family object created by \code{gee_family()}; must inherit
#'   from \code{"mfp2_gee_family"}.
#' @param y Response accepted by the inner GLM family. For binomial models this
#'   may be a two-column \code{cbind(successes, failures)} matrix.
#' @param id Cluster identifier, one value per observation. Observations from
#'   the same cluster must be contiguous.
#' @param waves Optional integer-like visit identifier, one value per
#'   observation, used to order repeated measurements and to represent gaps for
#'   AR(1) working correlation.
#' @param weights Optional numeric vector of prior weights.
#' @param offset Optional numeric offset vector.
#' @param control Resolved \code{geepack::geese.control()} list.
#'
#' @return The family object with a populated \code{prepared} slot.
#' @keywords internal
#' @noRd
prepare_gee_family <- function(family, y, id, waves = NULL, weights = NULL,
                               offset = NULL, control = NULL) {
  if (!inherits(family, "mfp2_gee_family")) {
    stop("Internal error: `family` is not a GEE family object.", call. = FALSE)
  }

  if (!requireNamespace("geepack", quietly = TRUE)) {
    stop(
      "`gee_family()` requires the 'geepack' package. ",
      "Install it with install.packages(\"geepack\").",
      call. = FALSE
    )
  }

  response_family <- family$response_family
  inner_string <- response_family$family
  if (!inner_string %in% mfp2_gee_supported_families()) {
    stop(
      sprintf(
        "! GEE currently supports the %s response families; got '%s'.",
        paste(mfp2_gee_supported_families(), collapse = ", "),
        inner_string
      ),
      call. = FALSE
    )
  }

  nobs <- NROW(y)

  # ---- Clustering identifier ------------------------------------------------
  if (is.null(id)) {
    stop(
      "! `id` is required for GEE models: supply a cluster identifier with one ",
      "value per observation.",
      call. = FALSE
    )
  }
  if (!is.atomic(id) || is.matrix(id) || !is.null(dim(id)) ||
      length(id) != nobs || anyNA(id)) {
    stop(
      "! `id` must be an atomic vector or factor with one non-missing value ",
      "per observation.",
      call. = FALSE
    )
  }

  # mfp2.default() groups clusters before preparing this family. Preserve this
  # invariant for direct calls to this internal helper as well.
  id_codes <- if (is.factor(id)) as.integer(id) else as.integer(factor(id, levels = unique(id)))
  run_starts <- c(TRUE, id_codes[-1L] != id_codes[-nobs])
  first_seen <- id_codes[run_starts]
  if (anyDuplicated(first_seen)) {
    stop(
      "Internal error: GEE observations must be grouped by `id` before ",
      "preparing the family.",
      call. = FALSE
    )
  }

  clus_run <- c(which(diff(id_codes) != 0), nobs)
  clusz <- c(clus_run[1L], diff(clus_run))
  n_clusters <- length(clusz)
  if (n_clusters < 2L) {
    stop(
      "! GEE models require at least two independent `id` clusters.",
      call. = FALSE
    )
  }

  # ---- Waves ----------------------------------------------------------------
  waves_int <- NULL
  if (!is.null(waves)) {
    if (!is.atomic(waves) || is.matrix(waves) || !is.null(dim(waves)) ||
        length(waves) != nobs || anyNA(waves)) {
      stop(
        "! `waves` must be an atomic vector with one non-missing value per ",
        "observation.",
        call. = FALSE
      )
    }
    waves_num <- if (is.factor(waves)) {
      # Preserve numeric visit labels and therefore real gaps. Factor codes
      # would turn levels such as "1" and "3" into consecutive values 1 and 2.
      suppressWarnings(as.numeric(as.character(waves)))
    } else {
      suppressWarnings(as.numeric(waves))
    }
    if (any(!is.finite(waves_num)) || any(waves_num != round(waves_num)) ||
        any(waves_num < 1) || any(waves_num > .Machine$integer.max)) {
      stop(
        paste0(
          "! `waves` must contain positive integer-like visit identifiers ",
          "within R's integer range."
        ),
        call. = FALSE
      )
    }
    waves_int <- as.integer(waves_num)

    # Waves must be unique within each cluster; duplicates make the working
    # correlation ill-defined.
    split_dup <- tapply(waves_int, id_codes, function(w) anyDuplicated(w) > 0L)
    if (any(unlist(split_dup))) {
      stop(
        "! `waves` must be unique within each `id` cluster.",
        call. = FALSE
      )
    }
  }

  # ---- Response, weights and grouped-binomial handling ----------------------
  if (is.null(weights)) weights <- rep.int(1, nobs)
  weights <- as.numeric(weights)
  if (length(weights) != nobs || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop(
      paste0(
        "! GEE `weights` must be finite and non-negative with one value per ",
        "observation."
      ),
      call. = FALSE
    )
  }

  if (is.null(offset)) offset <- rep.int(0, nobs)
  offset <- as.numeric(offset)
  if (length(offset) != nobs || anyNA(offset) || any(!is.finite(offset))) {
    stop(
      "! GEE `offset` must be finite with one value per observation.",
      call. = FALSE
    )
  }

  y_model <- y
  w_eff <- weights
  if (identical(inner_string, "binomial") && is.matrix(y) && ncol(y) == 2L) {
    # Match geepack::geeglm(): collapse cbind(successes, failures) to the
    # success proportion with cluster totals as weights.
    totals <- rowSums(y)
    if (any(totals <= 0)) {
      stop("! Grouped binomial rows must have a positive total count.", call. = FALSE)
    }
    y_model <- y[, 1L] / totals
    w_eff <- weights * totals
  } else if (is.matrix(y)) {
    stop(
      "! GEE response must be a numeric vector (grouped binomial may use a ",
      "two-column success/failure matrix).",
      call. = FALSE
    )
  } else {
    y_model <- as.numeric(y)
  }

  family$prepared <- list(
    response_family = response_family,
    inner_string    = inner_string,
    y_model         = y_model,
    weights         = w_eff,
    offset          = offset,
    id              = id,
    id_codes        = id_codes,
    waves           = waves_int,
    clusz           = clusz,
    n_clusters      = n_clusters,
    corstr          = family$corstr,
    std.err         = family$std.err,
    scale.fix       = family$scale.fix,
    scale.value     = family$scale.value,
    control         = normalize_gee_control(control),
    nobs            = nobs
  )

  family
}

#' Fit a GEE model for one MFP candidate or the final retained model
#'
#' @param x Design matrix (without intercept) for the current candidate.
#' @param y Response (ignored for the modelling response, which is taken from
#'   the prepared family so grouped-binomial and weight handling stay
#'   consistent across candidates; retained only for interface symmetry).
#' @param family Prepared GEE family (\code{"mfp2_gee_family"} with a populated
#'   \code{prepared} slot).
#' @param weights,offset Passed through by \code{fit_model()}; the authoritative
#'   values live in the prepared family.
#' @param control Resolved control (also cached in the prepared family).
#' @param fast Logical: \code{TRUE} uses \code{geepack::geese.fit()} for a
#'   lightweight candidate fit; \code{FALSE} uses \code{geepack::geeglm()} for
#'   the retained native fit.
#' @param calculate_fit_statistics Logical: return null/model deviance.
#' @param keep_fit Logical: retain the underlying fitted object.
#' @param keep_fitted_values Logical: retain only the fitted-value vector.
#' @param has_offset Logical: whether a user offset was supplied.
#' @param x_has_intercept Logical: whether \code{x} already carries an
#'   intercept column.
#' @param alpha_start Optional starting working-correlation parameter for a
#'   candidate fit. The correlation is still re-estimated for this model.
#' @param column_scales Optional per-column factors for a GEE FP/ACD candidate
#'   whose design is already scaled. Other fits scale their own design.
#' @param selection_criterion AIC/BIC candidate fits use quasi-likelihood only;
#'   p-value and retained fits still calculate their covariance and Wald score.
#' @param allow_failed_candidate Whether an unsuccessful candidate can be
#'   returned with `converged = FALSE` for exclusion by an FP power search.
#' @param reserved_names Reserved user-facing names for the final formula fit.
#'
#' @return A list following the \code{fit_model()} contract, with \code{logl}
#'   holding the quasi-likelihood surrogate and \code{is_gee = TRUE}.
#' @keywords internal
#' @noRd
fit_gee <- function(x,
                    y,
                    family,
                    weights = NULL,
                    offset = NULL,
                    control = NULL,
                    fast = TRUE,
                    calculate_fit_statistics = FALSE,
                    keep_fit = !fast,
                    keep_fitted_values = FALSE,
                    has_offset = FALSE,
                    x_has_intercept = FALSE,
                    alpha_start = NULL,
                    column_scales = NULL,
                    selection_criterion = NULL,
                    allow_failed_candidate = FALSE,
                    reserved_names = character()) {

  prepared <- family$prepared
  if (is.null(prepared)) {
    stop("Internal error: GEE family has not been prepared.", call. = FALSE)
  }

  response_family <- prepared$response_family
  inner_string    <- prepared$inner_string
  y_model         <- prepared$y_model
  w_eff           <- prepared$weights
  gee_offset      <- prepared$offset
  id_codes        <- prepared$id_codes
  waves_int       <- prepared$waves
  corstr          <- prepared$corstr
  std.err         <- prepared$std.err
  scale.fix       <- prepared$scale.fix
  scale.value     <- prepared$scale.value
  gee_control     <- prepared$control
  nobs            <- prepared$nobs
  reference_dispersion <- prepared$reference_dispersion
  qbic_n <- prepared$qbic_n

  # AIC/BIC select on quasi-likelihood alone. Their numerous candidate fits
  # need no jackknife, sandwich covariance extraction, or Wald eigendecomposition.
  # The retained fit still honors the requested standard-error method.
  ic_candidate <- isTRUE(fast) && !isTRUE(keep_fit) &&
    isTRUE(selection_criterion %in% c("aic", "bic"))
  gee_control <- gee_control_for_std_err(
    gee_control, if (ic_candidate) "san.se" else std.err
  )

  has_predictors <- !is.null(x) && NCOL(x) > 0L

  if (fast) {
    # ---- Candidate path: geese.fit() on the raw design matrix ----------------
    if (isTRUE(x_has_intercept)) {
      xx <- x
    } else {
      xx <- assemble_design_matrix(
        blocks = list(x),
        nobs = nobs,
        intercept = TRUE
      )
    }

    # Fractional-polynomial candidate bases can span many orders of magnitude
    # (e.g. FP2 terms with large powers), which makes geese.fit()'s internal
    # glm.fit() start values overflow and abort. Column-scaling the design
    # leaves the linear predictor, fitted means, and quasi-likelihood exactly
    # invariant (it is a pure reparametrisation), so we scale each non-intercept
    # column to unit maximum magnitude for the fit and unscale beta afterwards.
    pre_scaled <- !is.null(column_scales)
    if (pre_scaled) {
      if (!is.numeric(column_scales) || length(column_scales) != ncol(xx) ||
          any(!is.finite(column_scales)) || any(column_scales <= 0)) {
        stop("Internal error: invalid pre-scaled GEE column factors.", call. = FALSE)
      }
      col_scale <- column_scales
      xx_scaled <- xx
    } else {
      col_scale <- gee_design_column_scales(xx)
      intercept_col <- which(colnames(xx) == "(Intercept)")
      if (length(intercept_col) == 1L) col_scale[intercept_col] <- 1
      xx_scaled <- sweep(xx, 2L, col_scale, "/")
    }

    rank <- ncol(xx)

    # A previous candidate's alpha is only an initial value. If that start
    # fails for this mean model, retry with geepack's ordinary initialization.
    # This also protects the selection search from a bad warm start.
    if (identical(corstr, "independence") ||
        is.null(alpha_start) || !is.numeric(alpha_start) ||
        length(alpha_start) != 1L || any(!is.finite(alpha_start)) ||
        abs(alpha_start) >= 1) {
      alpha_start <- NULL
    }
    fit_with_alpha <- function(start) {
      tryCatch(
        geepack::geese.fit(
          x = xx_scaled,
          y = y_model,
          id = id_codes,
          offset = gee_offset,
          weights = w_eff,
          waves = waves_int,
          control = gee_control,
          family = response_family,
          corstr = corstr,
          scale.fix = scale.fix,
          scale.value = scale.value,
          alpha = start
        ),
        error = function(e) e
      )
    }

    # A power candidate may fail numerically. Only the FP power search is
    # allowed to skip it; required null/linear and retained models must fail.
    converged <- FALSE
    beta <- rep(NA_real_, rank)
    mu <- rep(NA_real_, nobs)
    quasi <- NA_real_

    fam_dev <- NA_real_
    robust_vcov <- matrix(NA_real_, nrow = rank, ncol = rank)
    selection_deviance <- NA_real_

    for (start in if (is.null(alpha_start)) list(NULL) else
                  list(alpha_start, NULL)) {
      ans <- fit_with_alpha(start)
      converged <- !inherits(ans, "error") &&
        (is.null(ans$error) || isTRUE(ans$error == 0))
      if (converged) {
        beta <- ans$beta / col_scale
        names(beta) <- colnames(xx)
        if (!ic_candidate) {
          inverse_scale <- 1 / col_scale
          raw_vcov <- tryCatch(
            gee_selected_vcov(ans, std.err, n_clusters = NULL),
            error = function(e) NULL
          )
          if (is.null(raw_vcov) ||
              !identical(dim(raw_vcov), c(rank, rank)) ||
              any(!is.finite(raw_vcov))) {
            converged <- FALSE
          } else {
            raw_vcov <- raw_vcov * outer(inverse_scale, inverse_scale)
            robust_vcov <- if (identical(std.err, "san.se")) {
              gee_robust_vcov(raw_vcov, prepared$n_clusters)
            } else {
              raw_vcov
            }
            dimnames(robust_vcov) <- list(colnames(xx), colnames(xx))
            selection_deviance <- gee_selection_deviance(
              beta,
              robust_vcov
            )
          }
        }
        eta <- if (pre_scaled) {
          as.numeric(xx_scaled %*% ans$beta) + gee_offset
        } else {
          as.numeric(xx %*% beta) + gee_offset
        }
        mu <- response_family$linkinv(eta)
        q <- gee_quasi_likelihood(y_model, mu, inner_string, w_eff)
        if (is.finite(q)) {
          quasi <- q
          if (isTRUE(calculate_fit_statistics) || isTRUE(keep_fit)) {
            fam_dev <- gee_family_deviance(y_model, mu, w_eff, response_family)
          }
        } else {
          converged <- FALSE
          selection_deviance <- NA_real_
          robust_vcov[,] <- NA_real_
        }
        if (converged && !ic_candidate &&
            !is.finite(selection_deviance)) {
          converged <- FALSE
        }
      }
      if (converged) break
    }

    if (!converged) {
      if (!isTRUE(allow_failed_candidate)) {
        validate_mfp_fit_result(
          logl = quasi, df = rank, converged = FALSE,
          family_string = "gee", fast = fast
        )
      }
      return(list(
        logl = NA_real_, selection_deviance = NA_real_,
        coefficients = rep(NA_real_, rank), rank = unname(rank),
        df = unname(rank), is_gee = TRUE, converged = FALSE
      ))
    }

    result <- list(
      logl = quasi,
      # geese.fit() estimates gamma even for lightweight IC candidates. Keep
      # only this scalar so the full reference can use the same successful
      # candidate path without requesting a Wald covariance or the full fit.
      gee_scale = if (!is.null(ans$gamma) && length(ans$gamma) > 0L)
        unname(ans$gamma[[1L]]) else NA_real_,
      reference_dispersion = reference_dispersion,
      qbic_n = qbic_n,
      family_deviance = fam_dev,
      selection_deviance = selection_deviance,
      robust_vcov = robust_vcov,
      coefficients = beta,
      rank = unname(rank),
      df = unname(rank),
      is_gee = TRUE,
      converged = converged
    )
    # Retain only the small successful correlation estimate, not the backend
    # fit, so the next FP candidate in this search can start from it.
    if (converged && !identical(corstr, "independence") &&
        is.numeric(ans$alpha) && length(ans$alpha) == 1L &&
        is.finite(ans$alpha) && abs(ans$alpha) < 1) {
      result$gee_alpha <- unname(ans$alpha)
    }

    # Only successful candidates enter the model-selection metrics.
    validate_mfp_fit_result(
      logl = result$logl,
      df = result$df,
      converged = converged,
      family_string = "gee",
      fast = fast
    )

    if (isTRUE(keep_fitted_values)) {
      result$fitted_values <- mu
    }
    if (isTRUE(keep_fit)) {
      result$fit <- ans
    }
    if (isTRUE(calculate_fit_statistics)) {
      result$null_deviance <- NA_real_
      result$model_deviance <- if (converged) {
        gee_family_deviance(y_model, mu, w_eff, response_family)
      } else {
        NA_real_
      }
    }

    return(result)
  }

  # ---- Final path: geepack::geeglm() for a native retained object -----------
  x_formula <- if (isTRUE(x_has_intercept)) {
    x[, -1L, drop = FALSE]
  } else {
    x
  }
  has_formula_predictors <- !is.null(x_formula) && NCOL(x_formula) > 0L

  if (has_formula_predictors) {
    if (is.null(colnames(x_formula)) || any(colnames(x_formula) == "")) {
      stop("Internal error: x must have non-empty column names.", call. = FALSE)
    }
    data <- data.frame(x_formula, check.names = FALSE)
    rhs <- paste(sprintf("`%s`", colnames(x_formula)), collapse = " + ")
  } else {
    data <- data.frame(row.names = seq_len(nobs))
    rhs <- "1"
  }

  used_names <- unique(c(names(data), as.character(reserved_names)))

  response_col <- mfp2_internal_name("response", used_names, preferred = "..mfp2_y")
  used_names <- c(used_names, response_col)
  data[[response_col]] <- y_model
  lhs <- response_col

  id_col <- mfp2_internal_name("id", used_names, preferred = "..mfp2_id")
  used_names <- c(used_names, id_col)
  data[[id_col]] <- id_codes

  waves_arg <- NULL
  if (!is.null(waves_int)) {
    waves_col <- mfp2_internal_name("waves", used_names, preferred = "..mfp2_waves")
    used_names <- c(used_names, waves_col)
    data[[waves_col]] <- waves_int
    waves_arg <- waves_col
  }

  weights_col <- mfp2_internal_name("weights", used_names, preferred = "..mfp2_w")
  used_names <- c(used_names, weights_col)
  data[[weights_col]] <- w_eff

  internal_names <- list(response = response_col, offset = NULL)
  if (isTRUE(has_offset)) {
    offset_col <- mfp2_internal_name("offset", used_names, preferred = "offset_")
    used_names <- c(used_names, offset_col)
    data[[offset_col]] <- gee_offset
    internal_names$offset <- offset_col
    rhs <- if (identical(rhs, "1")) {
      paste0("offset(", offset_col, ")")
    } else {
      paste0(rhs, " + offset(", offset_col, ")")
    }
  }

  formula <- stats::as.formula(paste(lhs, "~", rhs))

  # geepack::geeglm() evaluates `id`/`waves` in `data`; build the call so those
  # symbols resolve to the helper columns created above.
  cl <- list(
    quote(geepack::geeglm),
    formula = formula,
    family = response_family,
    data = quote(data),
    id = as.name(id_col),
    weights = as.name(weights_col),
    corstr = corstr,
    std.err = std.err,
    scale.fix = scale.fix,
    control = gee_control
  )
  if (!is.null(waves_arg)) {
    cl$waves <- as.name(waves_arg)
  }

  eval_env <- new.env(parent = environment())
  eval_env$data <- data
  fit <- eval(as.call(cl), envir = eval_env)

  unique_waves <- if (is.null(waves_int)) integer() else sort(unique(waves_int))
  needs_exact_waves <- identical(corstr, "ar1") &&
    length(unique_waves) > 1L && any(diff(unique_waves) > 1L)
  if (needs_exact_waves) {
    fit <- gee_refit_exact_waves(
      fit = fit,
      y = y_model,
      id = id_codes,
      waves = waves_int,
      offset = gee_offset,
      weights = w_eff,
      control = gee_control,
      response_family = response_family,
      corstr = corstr,
      scale.fix = scale.fix,
      scale.value = scale.value
    )
  }

  fit$mfp2_internal_names <- internal_names

  beta <- fit$coefficients
  mu <- fit$fitted.values
  quasi <- gee_quasi_likelihood(y_model, mu, inner_string, w_eff)
  rank <- length(beta)

  fam_dev <- gee_family_deviance(y_model, mu, w_eff, response_family)
  raw_vcov <- gee_selected_vcov(fit$geese, std.err, n_clusters = NULL)
  robust_vcov <- if (identical(std.err, "san.se")) {
    gee_robust_vcov(raw_vcov, prepared$n_clusters)
  } else {
    raw_vcov
  }
  dimnames(robust_vcov) <- list(names(beta), names(beta))

  # vcov.geeglm() and summary.geeglm() read the selected field directly. Store
  # the corrected sandwich covariance on the retained native object so every
  # downstream method (including mfp2 prediction and confidence intervals)
  # observes the same estimator used during selection.
  if (identical(std.err, "san.se")) {
    fit$geese$vbeta <- robust_vcov
  }
  # geeglm() does not retain its mean-model design by default. Summary and
  # downstream joint-term inference need the exact fitted columns, including
  # the intercept, so store the already constructed native GEE design rather
  # than reconstructing columns from coefficient-name patterns.
  fit$x <- fit$geese$X
  selection_deviance <- gee_selection_deviance(
    beta,
    robust_vcov
  )

  result <- list(
    logl = quasi,
    reference_dispersion = reference_dispersion,
    qbic_n = qbic_n,
    family_deviance = fam_dev,
    selection_deviance = selection_deviance,
    robust_vcov = robust_vcov,
    coefficients = beta,
    rank = unname(rank),
    df = unname(rank),
    is_gee = TRUE,
    converged = is.null(fit$geese$error) || isTRUE(fit$geese$error == 0)
  )

  validate_mfp_fit_result(
    logl = result$logl,
    df = result$df,
    converged = result$converged,
    family_string = "gee",
    fast = fast
  )

  if (isTRUE(keep_fitted_values)) {
    result$fitted_values <- mu
  }
  if (isTRUE(keep_fit)) {
    result$fit <- fit
  }
  if (isTRUE(calculate_fit_statistics)) {
    result$model_deviance <- fam_dev
    # Null model: intercept (+ offset) only, same GEE settings.
    null_mu <- gee_null_fitted_mean(
      y_model, w_eff, gee_offset, id_codes, waves_int,
      response_family, corstr, scale.fix, scale.value, gee_control
    )
    result$null_deviance <- gee_family_deviance(y_model, null_mu, w_eff, response_family)
  }

  result
}

#' Family deviance for a GEE fit
#'
#' The GLM family deviance evaluated at the GEE fitted means. This equals the
#' deviance Stata's `xtgee` reports as `e(deviance)` and is retained for model-
#' fit reporting. It is distinct from Stata `mfp`'s negative-Wald selection
#' quantity and from the default robust block-Wald tests.
#'
#' @keywords internal
#' @noRd
gee_family_deviance <- function(y, mu, weights, response_family) {
  dev <- sum(response_family$dev.resids(y, mu, weights))
  unname(dev)
}

#' Intercept-only GEE fitted mean, for the null-model deviance
#'
#' @keywords internal
#' @noRd
gee_null_fitted_mean <- function(y, weights, offset, id_codes, waves_int,
                                 response_family, corstr, scale.fix,
                                 scale.value, control) {
  xx <- matrix(1, nrow = length(y), ncol = 1L,
               dimnames = list(NULL, "(Intercept)"))
  ans <- geepack::geese.fit(
    x = xx, y = y, id = id_codes,
    offset = offset, weights = weights, waves = waves_int,
    control = control, family = response_family, corstr = corstr,
    scale.fix = scale.fix, scale.value = scale.value
  )
  eta <- as.numeric(xx %*% ans$beta) + offset
  response_family$linkinv(eta)
}

#' Model metrics for GEE candidate and final fits (parallel to
#' `calculate_model_metrics()`)
#'
#' Produces the same named vector as \code{calculate_model_metrics()} so the
#' shared MFP selection engine can use GEE semantics: \code{logl} holds the
#' quasi-likelihood \eqn{Q} used for information-criterion selection.
#' \code{deviance_rs} contains the negative overall robust Wald chi-square
#' used by \code{mfp: xtgee} for p-value selection; information-criterion
#' candidate fits omit this Wald calculation. \code{aic} uses the QICu form
#' \eqn{-2Q/\phi + 2p} of Pan (2001) with the package's extra FP-power df charge.
#' \code{bic} is the package's BIC-like quasi-likelihood score
#' \eqn{-2Q/\phi + \log(n)\,p}, where \eqn{\phi} is the shared reference
#' dispersion and \eqn{n} is the configured number of clusters or observations.
#'
#' The p-value statistics are kept separate from \eqn{Q}; information criteria
#' always use \eqn{-2Q/\phi}.
#'
#' @param obj A GEE fit-result list from \code{fit_gee()} (\code{is_gee = TRUE}).
#' @param n_obs Number of independent clusters \eqn{K}, used for residual
#'   degrees of freedom and as the fallback QBIC penalty count for legacy fits.
#' @param df_additional Extra degrees of freedom for FP/ACD power searches,
#'   added once to the parameter count.
#'
#' @return A named numeric vector with entries \code{logl}, \code{df},
#'   \code{deviance_rs}, \code{deviance_gaussian}, \code{aic}, \code{bic}, and
#'   \code{df_resid}. \code{deviance_gaussian} is \code{NA}: GEE uses either
#'   the Stata-compatible negative-Wald calculation rather than the Gaussian
#'   F-test.
#'
#' @references
#' Pan, W. (2001). Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, 57, 120--125.
#' @keywords internal
#' @noRd
calculate_gee_metrics <- function(obj, n_obs, df_additional = 0) {
  quasi <- obj$logl
  df <- obj$df + df_additional
  dispersion <- if (is.null(obj$reference_dispersion)) 1 else
    obj$reference_dispersion
  qbic_n <- if (is.null(obj$qbic_n)) n_obs else obj$qbic_n

  fit_rank <- obj$rank
  regression_df <- if (is.numeric(fit_rank) && length(fit_rank) == 1L &&
                       !is.na(fit_rank) && is.finite(fit_rank)) {
    unname(fit_rank)
  } else {
    sum(!is.na(obj$coefficients))
  }
  # GEE QICu counts only mean-model parameters; there are no nuisance df.
  df_for_resid <- df - max(0, obj$df - regression_df)

  # GEE p-value selection uses the negative overall robust Wald chi-square.
  # IC candidates deliberately omit the Wald covariance, so their descriptive
  # deviance slot uses -2Q; failed fits are excluded before reaching here.
  stata_dev <- obj$selection_deviance
  dev_rs <- if (!is.null(stata_dev) && is.finite(stata_dev)) {
    stata_dev
  } else {
    -2 * quasi
  }

  c(
    logl = quasi,
    df = df,
    deviance_rs = dev_rs,
    deviance_gaussian = NA_real_,
    # QICu and QBIC use reference-scaled Q, not the Wald deviance D.
    aic = -2 * quasi / dispersion + 2 * df,
    bic = -2 * quasi / dispersion + log(qbic_n) * df,
    df_resid = n_obs - df_for_resid
  )
}

#' Final Selection Score for a Retained GEE Model
#'
#' Returns the same criterion-scale value used by the MFP selection engine:
#' the method-specific p-value statistic, QICu, or QBIC. Keeping this value on
#' the returned object makes the selected model auditable without recomputing
#' or guessing which GEE comparison scale was active.
#'
#' @param obj The retained GEE fit-result list from `fit_gee()`.
#' @param criterion The active selection criterion: `"pvalue"`, `"aic"`, or
#'   `"bic"`.
#' @param n_clusters Number of independent clusters used for the QBIC penalty.
#' @param df_additional Additional degrees of freedom for retained FP power
#'   searches. This does not affect the p-value selection statistic.
#'
#' @return The final score on the scale of `criterion`.
#' @keywords internal
#' @noRd
gee_final_selection_score <- function(obj, criterion, n_clusters,
                                      df_additional = 0) {
  criterion <- match.arg(criterion, c("pvalue", "aic", "bic"))
  metrics <- calculate_gee_metrics(
    obj, n_obs = n_clusters, df_additional = df_additional
  )
  metric_name <- switch(
    criterion,
    pvalue = "deviance_rs",
    aic = "aic",
    bic = "bic"
  )
  score <- unname(metrics[[metric_name]])

  if (!is.numeric(score) || length(score) != 1L || !is.finite(score)) {
    stop(
      "Internal error: the retained GEE model has no finite selection score.",
      call. = FALSE
    )
  }
  score
}
