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

#' Group GEE Observations by Cluster
#'
#' Returns a row permutation that places observations from each cluster
#' together, as required by [geepack::geese.fit()] and [geepack::geeglm()].
#' Clusters follow their order of first appearance in `id`; observations
#' retain their original order within each cluster.
#'
#' @param id Vector of cluster identifiers aligned with the fitted sample.
#'
#' @return Integer vector of row indices for reordering the fitted sample.
#'
#' @keywords internal
#' @noRd
mfp2_gee_cluster_order <- function(id) {
  cluster <- match(id, unique(id))
  order(cluster, seq_along(cluster), method = "radix")
}

#' Robust Wald Test for a Coefficient Block
#'
#' Tests the joint null hypothesis that the coefficients selected by `indices`
#' are zero, using their supplied sandwich covariance matrix. The statistic
#' is compared with a chi-square distribution with one degree of freedom per
#' tested coefficient. A rank-deficient covariance block cannot support this
#' joint test and produces missing test values.
#'
#' @param coefficients Named numeric vector of fitted coefficients.
#' @param robust_vcov Sandwich covariance matrix in the same coefficient order.
#' @param indices Integer positions of the coefficients to test.
#' @param compute_pvalue Logical. If `FALSE`, omit the chi-square tail
#'   probability; model-selection scoring requires only the statistic.
#'
#' @return A list containing `statistic`, `pvalue`, `df`, and `dev_diff`.
#'   `dev_diff` is an alias for the Wald statistic used by downstream
#'   comparison code. For an empty block, the statistic and degrees of
#'   freedom are zero and the p-value is one. Invalid or rank-deficient
#'   blocks return missing test values.
#'
#' @keywords internal
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

#' Apply a Finite-Cluster Correction to a Robust Covariance
#'
#' Multiplies a sandwich covariance matrix by \eqn{K/(K-1)}, where \eqn{K}
#' is the number of independent clusters. This matches the finite-cluster
#' multiplier used by Stata's `xtgee` with `vce(robust)`.
#'
#' @param robust_vcov Sandwich covariance matrix.
#' @param n_clusters Number of independent clusters, or `NULL` to leave the
#'   covariance unchanged.
#'
#' @return The corrected covariance matrix. If `n_clusters` is `NULL`,
#'   returns `robust_vcov` unchanged. If it is invalid or no greater than
#'   one, returns a matrix of missing values.
#'
#' @keywords internal
#' @noRd
gee_robust_vcov <- function(robust_vcov, n_clusters = NULL) {
  if (is.null(n_clusters)) {
     return(robust_vcov)
  }

  # The K/(K - 1) correction requires one finite cluster count greater than 1.
  # Return a missing covariance matrix when it cannot be calculated. That
  # prevents a Wald test from using an invalid correction.
  n_clusters <- as.numeric(n_clusters)
  if (length(n_clusters) != 1L || !is.finite(n_clusters) || n_clusters <= 1) {
    corrected <- as.matrix(robust_vcov)
    corrected[,] <- NA_real_
    return(corrected)
  }

  robust_vcov * (n_clusters / (n_clusters - 1))
}

#' Select the Covariance Matrix for a GEE Fit
#'
#' Extracts the covariance estimator requested by `std.err` from a
#' [geepack::geese.fit()] result. The sandwich estimator (`"san.se"`)
#' receives the finite-cluster correction when `n_clusters` is supplied;
#' jackknife estimators are returned unchanged.
#'
#' @param geese Fitted `geese` object containing the covariance matrices.
#' @param std.err Estimator to select: `"san.se"`, `"jack"`, `"j1s"`, or
#'   `"fij"`.
#' @param n_clusters Number of independent clusters used to correct
#'   `"san.se"`, or `NULL` to leave it uncorrected.
#'
#' @return The selected covariance matrix. Raises an error if the requested
#'   estimator was not computed.
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

#' Configure GEE Covariance Estimation
#'
#' Sets the `jack`, `j1s`, and `fij` control flags for
#' [geepack::geese.fit()] to match the requested standard-error method.
#' For `"san.se"`, all three flags are zero and the default sandwich
#' covariance is used. This follows the mapping in [geepack::geeglm()].
#'
#' @param control GEE fitting control list.
#' @param std.err Requested method: `"san.se"`, `"jack"`, `"j1s"`, or
#'   `"fij"`.
#'
#' @return The control list with the three covariance flags set.
#'
#' @keywords internal
#' @noRd
gee_control_for_std_err <- function(control, std.err) {
  control$jack <- as.integer(identical(std.err, "jack"))
  control$j1s <- as.integer(identical(std.err, "j1s"))
  control$fij <- as.integer(identical(std.err, "fij"))
  control
}

#' Prepare Scaled Matrices for GEE Candidate Fits
#'
#' Scales the candidate design template and shared FP/ACD basis once per
#' focal search. Each non-intercept column is divided by its maximum absolute
#' value; zero or non-finite maxima use a scale factor of one. Candidate fits
#' reuse the scaled basis and update only the focal columns of a private
#' design copy.
#'
#' @param design Numeric candidate design template, including any intercept,
#'   focal placeholders, and adjustment columns.
#' @param basis Numeric matrix of shared focal transformations, with the same
#'   number and order of rows as `design`.
#'
#' @return A list containing private scaled matrices `design` and `basis`,
#'   and named vectors `column_scales` and `basis_scales`. The input matrices
#'   are not modified.
#'
#' @keywords internal
#' @noRd
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

#' Calculate GEE Design Column Scales
#'
#' Returns the maximum absolute value in each design column. Columns whose
#' maximum is zero or non-finite receive a scale factor of one. This function
#' calculates the factors; it does not modify the matrix.
#'
#' @param x Numeric design matrix.
#'
#' @return Numeric vector of positive scale factors, named by the column
#'   names of `x` when present.
#'
#' @keywords internal
#' @noRd
gee_design_column_scales <- function(x) {
  scales <- vapply(seq_len(ncol(x)), function(j) {
    s <- max(abs(x[, j]))
    if (!is.finite(s) || s <= 0) 1 else s
  }, numeric(1L))
  names(scales) <- colnames(x)
  scales
}

#' Refit a Final GEE Model Using Exact Wave Distances
#'
#' `geeglm()` converts waves to factor codes, which can collapse gaps between
#' observed integer waves. For an AR(1) model with such gaps, this function
#' refits the GEE core with [geepack::geese.fit()] using the original wave
#' numbers. It then updates the corresponding fields of the retained
#' `geeglm` object. The caller determines whether a refit is needed.
#'
#' @param fit Final model fitted by [geepack::geeglm()].
#' @param y Response vector aligned with the fitted rows.
#' @param id Cluster identifiers aligned with `y`.
#' @param waves Integer wave numbers retaining the intended gaps.
#' @param offset,weights Offset and case-weight vectors aligned with `y`.
#' @param control GEE fitting control list.
#' @param response_family GLM family used for the response.
#' @param corstr Working correlation structure.
#' @param scale.fix,scale.value Scale settings passed to
#'   [geepack::geese.fit()].
#'
#' @return The updated `geeglm` object, with its GEE core and fitted-value
#'   fields calculated using the supplied wave numbers.
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

#' Calculate the GEE Selection Deviance
#'
#' Returns the negative overall Wald chi-square for all coefficients except
#' the intercept. This follows the `-e(chi2)` scoring convention used by
#' Stata's `mfp: xtgee`; it is not a likelihood deviance. The Wald statistic
#' uses the supplied covariance estimator, so results depend on whether
#' sandwich or jackknife covariance was selected.
#'
#' @param coefficients Named vector of fitted coefficients.
#' @param robust_vcov Covariance matrix in coefficient order.
#' @param n_clusters Optional number of clusters. When supplied, applies
#'   the finite-cluster correction to `robust_vcov`.
#'
#' @return Numeric scalar: the negative Wald chi-square, zero when there
#'   are no non-intercept coefficients, or `NA` if the test is invalid.
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

#' Calculate the Quasi-Likelihood for a GEE Fit
#'
#' Sums the family-specific quasi-likelihood contributions used for GEE
#' information-criterion selection. Contributions are multiplied by effective
#' prior weights; zero-weight observations are excluded before evaluation.
#' For grouped binomial data, `y` is the success proportion and `weights`
#' includes the number of trials.
#'
#' The Gamma contribution is \eqn{-y/\mu-\log(\mu)}, obtained by integrating
#' \eqn{(y-\mu)/\mu^2}.
#'
#' @param y Numeric response vector; a success proportion for grouped
#'   binomial data.
#' @param mu Numeric vector of fitted means.
#' @param family_string Response family: `"gaussian"`, `"binomial"`,
#'   `"poisson"`, or `"Gamma"`.
#' @param weights Optional vector of finite, nonnegative effective prior
#'   weights. Defaults to one per observation.
#'
#' @return Numeric scalar containing the weighted quasi-likelihood.
#'
#' @references
#' Pan, W. (2001). Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, 57, 120--125.
#'
#' @keywords internal
#' @noRd
gee_quasi_likelihood <- function(y, mu, family_string, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep.int(1, length(y))
  }

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

#' Normalise Controls for the GEE Fitter
#'
#' Passes user-supplied settings to [geepack::geese.control()] and fills
#' unspecified fields with its defaults. Common settings include `maxit`
#' for the iteration limit, `epsilon` for the convergence tolerance, and
#' `trace` for iteration output. Invalid settings raise an error.
#'
#' @param control Named list of settings accepted by
#'   [geepack::geese.control()], or `NULL` to use its defaults.
#'
#' @return A validated control list from [geepack::geese.control()].
#'
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

#' Prepare a GEE Family for Repeated MFP Fits
#'
#' Validates the GEE family, grouped cluster identifiers, optional waves,
#' weights, and offsets. For a two-column binomial response, converts counts
#' to success proportions and multiplies prior weights by the trial totals.
#' It then caches the inputs shared by candidate fits in `family$prepared`,
#' including cluster sizes and the normalised GEE control list.
#'
#' @param family Object created by [gee_family()].
#' @param y Response vector, or a two-column matrix of successes and
#'   failures for grouped binomial data.
#' @param id Cluster identifier with one non-missing value per observation.
#'   Observations must already be contiguous within each cluster.
#' @param waves Optional positive integer-like visit identifiers, unique
#'   within each cluster.
#' @param weights Optional finite, nonnegative prior weights.
#' @param offset Optional finite numeric offset vector.
#' @param control Optional list of settings accepted by
#'   [geepack::geese.control()], or `NULL` for its defaults.
#'
#' @return `family` with a populated `prepared` component.
#'
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

#' Fit a GEE Candidate or Retained Final Model
#'
#' Uses [geepack::geese.fit()] for matrix-based MFP candidate fits and
#' [geepack::geeglm()] for the retained final model. Response, cluster,
#' weight, offset, and control inputs are reused from `family$prepared`.
#' The `logl` result contains quasi-likelihood, not a likelihood.
#'
#' @param x Numeric design matrix. It may already contain an intercept when
#'   `x_has_intercept` is `TRUE`.
#' @param y Response passed for interface consistency; the fitted response
#'   is read from the prepared family.
#' @param family Prepared object created by [gee_family()].
#' @param weights,offset,control Arguments passed through by `fit_model()`.
#'   The values cached in the prepared family are used for fitting.
#' @param fast Logical. Use `geese.fit()` for a candidate when `TRUE`, or
#'   `geeglm()` for the retained final model when `FALSE`.
#' @param calculate_fit_statistics Logical. Return the family deviance of
#'   the fitted model. The final path also fits an intercept-only GEE to
#'   obtain null deviance; the candidate path returns `NA` for null deviance.
#' @param keep_fit Logical. Include the underlying fitted object in `fit`.
#'   Defaults to `!fast`.
#' @param keep_fitted_values Logical. Include fitted means in
#'   `fitted_values`.
#' @param has_offset Logical. Include the prepared offset in the formula
#'   used for the final `geeglm()` fit.
#' @param x_has_intercept Logical. Whether `x` already has an intercept
#'   column.
#' @param alpha_start Optional starting working-correlation estimate for
#'   a candidate fit. The correlation is re-estimated; a failed start is
#'   retried with the fitter's default initialisation.
#' @param column_scales Optional factors for an already scaled candidate
#'   design. Otherwise, candidate columns are scaled before fitting.
#' @param selection_criterion Selection criterion. AIC and BIC candidates
#'   omit covariance and Wald-score calculations.
#' @param allow_failed_candidate Logical. If `TRUE`, a failed candidate
#'   can return `converged = FALSE` for exclusion from the search.
#' @param reserved_names Names to avoid when creating columns for the
#'   final formula-based fit.
#'
#' @return A list following the `fit_model()` contract, including
#'   quasi-likelihood in `logl`, coefficients, degrees of freedom,
#'   convergence status, and `is_gee = TRUE`. Optional components depend
#'   on the reporting and retention arguments.
#'
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

#' Calculate Family Deviance for a GEE Fit
#'
#' Sums the GLM family's deviance residuals at the GEE fitted means for
#' model-fit reporting. This is separate from the negative overall Wald
#' statistic used for GEE p-value selection.
#'
#' @param y Response vector used in the GEE fit.
#' @param mu Vector of fitted means.
#' @param weights Effective prior weights aligned with `y`.
#' @param response_family GLM family object providing `dev.resids`.
#'
#' @return Unnamed numeric scalar containing the family deviance.
#'
#' @keywords internal
#' @noRd
gee_family_deviance <- function(y, mu, weights, response_family) {
  dev <- sum(response_family$dev.resids(y, mu, weights))
  unname(dev)
}

#' Fit an Intercept-Only GEE for Null Deviance
#'
#' Fits an intercept-only GEE with the supplied offset and the same cluster,
#' wave, weight, family, and working-correlation settings as the final model.
#' Returns fitted means for calculating the null-model family deviance.
#'
#' @param y Response vector used in the GEE fit.
#' @param weights Effective prior weights aligned with `y`.
#' @param offset Numeric offset vector aligned with `y`.
#' @param id_codes Grouped integer cluster identifiers.
#' @param waves_int Optional integer wave identifiers, or `NULL`.
#' @param response_family GLM response family.
#' @param corstr Working correlation structure.
#' @param scale.fix,scale.value Scale settings passed to
#'   [geepack::geese.fit()].
#' @param control GEE fitting control list.
#'
#' @return Numeric vector of fitted means in observation order.
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

#' Calculate Selection Metrics for a GEE Fit
#'
#' Returns the metric vector used by the shared MFP selection engine.
#' `logl` stores quasi-likelihood \eqn{Q}, not a log-likelihood. When an
#' overall Wald statistic was calculated, `deviance_rs` stores its negative
#' for p-value selection; otherwise it stores \eqn{-2Q}.
#'
#' The `aic` entry is the package's QICu-style score
#' \eqn{-2Q/\phi + 2p}. The `bic` entry is the BIC-like score
#' \eqn{-2Q/\phi + \log(n)p}. Here \eqn{p} includes any additional FP/ACD
#' selection degrees of freedom, \eqn{\phi} is the shared reference
#' dispersion, and \eqn{n} is the configured cluster or observation count.
#'
#' @param obj GEE fit-result list from `fit_gee()`.
#' @param n_obs Number of independent clusters in the GEE fitting path.
#'   Also used as the QBIC penalty count if `obj$qbic_n` is absent.
#' @param df_additional Additional degrees of freedom charged for the
#'   FP/ACD selection search.
#'
#' @return Named numeric vector containing `logl`, `df`, `deviance_rs`,
#'   `deviance_gaussian`, `aic`, `bic`, and `df_resid`.
#'   `deviance_gaussian` is `NA` because the Gaussian F-test does not
#'   apply to this GEE selection path.
#'
#' @references
#' Pan, W. (2001). Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, 57, 120--125.
#'
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

#' Get the Final GEE Selection Score
#'
#' Calculates the retained model's score on the active selection scale:
#' negative overall Wald chi-square for `"pvalue"`, QICu-style score for
#' `"aic"`, or the BIC-like quasi-likelihood score for `"bic"`. The value
#' is stored on the returned MFP model for inspection.
#'
#' @param obj Retained GEE fit-result list from `fit_gee()`.
#' @param criterion Selection criterion: `"pvalue"`, `"aic"`, or `"bic"`.
#' @param n_clusters Number of independent clusters. Used as the QBIC
#'   penalty count when `obj$qbic_n` is absent.
#' @param df_additional Additional FP selection degrees of freedom for
#'   the information-criterion penalties. Does not change the Wald score.
#'
#' @return A finite, unnamed numeric selection score. Raises an error if
#'   no finite score is available.
#'
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
