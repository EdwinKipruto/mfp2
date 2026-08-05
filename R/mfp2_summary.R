# =============================================================================
# summary_mfp2.R
#
# The structured summary method for `mfp2` model fits, and its dedicated
# `print` method. Provides an MFP-aware alternative to the raw
# `summary.glm()` / `summary.coxph()` output that ships from the underlying
# fitter.
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
#' the raw [stats::summary.glm()] or [survival::summary.coxph()] output, this
#' method separates interpretable linear terms from fractional-polynomial (FP)
#' and other nonlinear terms, whose individual coefficients are curve
#' parameters rather than per-unit effect sizes.
#'
#' @details
#' The summary is organised into the following sections:
#' \itemize{
#'   \item \strong{Selection Overview}: the same selection overview shown by
#'     [print.mfp2()].
#'   \item \strong{Linear Terms}: variables entering the model as a single
#'     linear term (including binary-only spike variables and factor levels),
#'     with the coefficient, standard error, test statistic, p-value, and --
#'     for non-Gaussian families -- the exponentiated coefficient and its
#'     confidence interval.
#'   \item \strong{Nonlinear Terms}: one row per variable modelled by an FP1,
#'     FP2, ACD, or spike/catzero compound form. Each row reports a joint
#'     likelihood-ratio test (LRT) of all of that variable's terms, using the
#'     selection-adjusted degrees of freedom stored in \code{fp_terms}. The
#'     LRT compares the final MFP model with the model obtained by dropping
#'     that variable's columns while holding all other functional forms fixed.
#'   \item \strong{Model Fit}: full-linear and final-MFP deviances for GLMs,
#'     or minus twice the partial log-likelihood for Cox models, with model
#'     degrees of freedom.
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
#' \code{basis = TRUE} to append the raw basis coefficients. For an active
#' ACD component, either option also prints the stored definition of
#' \eqn{A(x)}, including its internal training scale. Set
#' \code{raw = TRUE} to obtain the underlying [stats::summary.glm()] or
#' [survival::summary.coxph()] object instead.
#'
#' @param object A fitted [mfp2()] object.
#' @param formulas Logical. If \code{TRUE}, append the fitted-function formula
#'   for each nonlinear variable. Default \code{FALSE}.
#' @param basis Logical. If \code{TRUE}, append a table of raw basis
#'   coefficients for the nonlinear terms. Default \code{FALSE}.
#' @param raw Logical. If \code{TRUE}, bypass the structured summary and return
#'   the raw [stats::summary.glm()] or [survival::summary.coxph()] object, with
#'   the fitting call replaced by the original [mfp2()] call. Default
#'   \code{FALSE}.
#' @param digits Number of significant digits used when printing. Defaults to
#'   \code{max(3L, getOption("digits") - 3L)}.
#' @param ... Further arguments. When \code{raw = TRUE}, passed to the
#'   underlying summary method; otherwise ignored.
#'
#' @return
#' When \code{raw = FALSE}, an object of class \code{"summary.mfp2"}: a list
#' with components \code{call}, \code{family}, \code{criterion},
#' \code{converged}, \code{n}, \code{nevents}, \code{function_table},
#' \code{linear_terms}, \code{nonlinear_terms}, \code{basis} (or \code{NULL}),
#' \code{formulas} (or \code{NULL}), \code{acd_definitions} (or
#' \code{NULL}), \code{fit}, and \code{raw_summary}. A
#' dedicated \code{print} method renders these. When \code{raw = TRUE}, the
#' underlying summary object.
#'
#' @seealso [mfp2()], [print.mfp2()], [stats::summary.glm()],
#'   [survival::summary.coxph()]
#'
#' @export
summary.mfp2 <- function(object,
                         formulas = FALSE,
                         basis = FALSE,
                         raw = FALSE,
                         digits = max(3L, getOption("digits") - 3L),
                         ...) {
  if (!inherits(object, "mfp2")) {
    stop("The object is not an mfp2 object.", call. = FALSE)
  }

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

  # ---------------------------------------------------------------------------
  # Gather the pieces.
  # ---------------------------------------------------------------------------
  classified <- mfp2_summary_classify_terms(object)
  raw_summary <- mfp2_summary_raw(object)

  linear_terms <- mfp2_summary_linear_table(object, classified, raw_summary)
  nonlinear_terms <- mfp2_summary_nonlinear_table(object, classified)

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
    function_table  = classified$function_table,
    linear_terms    = linear_terms,
    nonlinear_terms = nonlinear_terms,
    basis           = basis_table,
    formulas        = formula_strings,
    acd_definitions = acd_definitions,
    fit             = mfp2_summary_fit_stats(object),
    raw_summary     = raw_summary,
    digits          = digits
  )
  class(out) <- "summary.mfp2"
  out
}


# =============================================================================
# Internal helpers
# =============================================================================

# Number of observations, robust to missing nobs field.
mfp2_summary_nobs <- function(object) {
  if (!is.null(object$nobs)) return(as.integer(object$nobs))
  n <- tryCatch(NROW(object$y), error = function(e) NA_integer_)
  as.integer(n)
}

# Number of events for Cox models; NA otherwise.
mfp2_summary_nevents <- function(object) {
  if (!identical(object$family_string, "cox")) return(NA_integer_)
  if (!is.null(object$nevents)) return(as.integer(object$nevents))
  y <- object$y
  if (inherits(y, "Surv")) {
    status_col <- ncol(y)
    return(as.integer(sum(y[, status_col] == 1, na.rm = TRUE)))
  }
  NA_integer_
}

# Human-readable selection criterion.
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

# Compute (once) the underlying summary.glm / summary.coxph object, with the
# call replaced by the original mfp2() call. Stored on the result so users can
# access it without refitting.
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

# NextMethod() only works inside a method; this dispatches manually to the
# next class after "mfp2" so the raw summary can be built from a helper.
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

# Split selected variables into linear vs nonlinear, and map each variable to
# its transformed coefficient columns (by name prefix). Also returns the
# selection overview table used in the printed header.
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

mfp2_summary_flag <- function(fp_terms, name) {
  if (!name %in% names(fp_terms)) return(rep(FALSE, nrow(fp_terms)))
  v <- fp_terms[[name]]
  if (is.logical(v)) return(ifelse(is.na(v), FALSE, v))
  tolower(as.character(v)) %in% c("true", "t", "yes", "y", "1")
}

# SAZ decision codes: package uses integer codes; "binary only" is the code
# meaning the continuous component was dropped. Matches saz_decision_label().
mfp2_summary_saz_is_binary_only <- function(code) {
  if (is.na(code)) return(FALSE)
  # In mfp2, spike_dec == 3 corresponds to "binary only".
  identical(as.integer(code), 3L)
}

# Build the human-readable functional-form label for one variable.
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

# Selection overview table (Variable, Selected, df, Function).
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

# Build a named character vector mapping each fitted coefficient name to a
# user-facing display name, drawing on the two authoritative lookups stored on
# the fitted mfp2 object at fit time:
#
#   term_to_columns             user variable name  -> raw source columns
#   transformed_to_model_columns raw source column  -> fitted coef name
#
# Composing the two gives an exact `coef_name -> display_name` map for every
# selected variable, with the internal fitter's suffix convention treated as
# an implementation detail rather than something the summary depends on.
#
# Returns a named character vector: names are fitted coefficient names, values
# are display names. Returns an empty named vector if either authoritative
# lookup is missing (in which case the caller falls back to the regex strip).
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

  is_gaussian <- identical(object$family_string, "gaussian")

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

  if (!is_gaussian) {
    df$exp_coef <- exp(est)
    df$ci_lower <- exp(est - 1.96 * se)
    df$ci_upper <- exp(est + 1.96 * se)
  } else {
    df$ci_lower <- est - 1.96 * se
    df$ci_upper <- est + 1.96 * se
  }

  df
}

# Standardised coefficient matrix (estimate, se, statistic, pvalue) with row
# names equal to coefficient names, sourced from the raw summary where
# possible and falling back to vcov().
mfp2_summary_coef_matrix <- function(object, raw_summary) {
  # Try the raw summary's coefficient table first.
  cmat <- NULL
  if (!is.null(raw_summary) && !is.null(raw_summary$coefficients)) {
    cmat <- raw_summary$coefficients
  }

  est <- object$coefficients
  nm <- names(est)

  if (!is.null(cmat) && nrow(cmat) == length(est)) {
    # Cox: columns are coef, exp(coef), se(coef), z, Pr(>|z|)
    # GLM: columns are Estimate, Std. Error, z/t value, Pr(>|.|)
    cn <- colnames(cmat)
    se_col <- grep("se|Std", cn, ignore.case = TRUE)[1]
    stat_col <- grep("^z$|value|^t$", cn, ignore.case = TRUE)[1]
    p_col <- grep("Pr|p.?value", cn, ignore.case = TRUE)[1]
    est_col <- grep("coef|Estimate", cn, ignore.case = TRUE)[1]

    out <- cbind(
      estimate  = cmat[, est_col],
      se        = cmat[, se_col],
      statistic = cmat[, stat_col],
      pvalue    = cmat[, p_col]
    )
    rownames(out) <- rownames(cmat)
    return(out)
  }

  # Fallback: compute from vcov().
  V <- tryCatch(stats::vcov(object), error = function(e) NULL)
  se <- if (!is.null(V)) sqrt(diag(V)) else rep(NA_real_, length(est))
  stat <- est / se
  pval <- 2 * stats::pnorm(abs(stat), lower.tail = FALSE)

  out <- cbind(
    estimate  = est,
    se        = se,
    statistic = stat,
    pvalue    = pval
  )
  rownames(out) <- nm
  out
}

# ---------------------------------------------------------------------------
# Nonlinear terms table (LRT, one row per variable)
# ---------------------------------------------------------------------------

mfp2_summary_nonlinear_table <- function(object, classified) {
  nl_vars <- classified$variable_names[classified$is_nonlinear]
  if (length(nl_vars) == 0L) {
    return(data.frame())
  }

  rows <- lapply(nl_vars, function(v) {
    idx <- match(v, classified$variable_names)
    form <- mfp2_summary_form_label(
      classified$powers_by_var[[idx]],
      classified$acd[idx], classified$zero[idx],
      classified$catzero[idx], classified$spike[idx],
      classified$spike_dec[idx], TRUE
    )
    lrt <- mfp2_summary_lrt_drop_variable(object, classified, v)

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

# Joint likelihood-ratio test for dropping one variable's columns, holding all
# other functional forms fixed. The LR statistic is computed by refitting the
# reduced model with the same deviance definition as the full model; the df is
# the selection-adjusted df_final from fp_terms.
mfp2_summary_lrt_drop_variable <- function(object, classified, v) {
  cols_v <- intersect(colnames(object$x), classified$cols_by_var[[v]])
  idx <- match(v, classified$variable_names)
  df_v <- classified$df_final[idx]
  if (is.na(df_v)) df_v <- length(cols_v)

  na_result <- list(lr = NA_real_, df = df_v, p = NA_real_)

  x_full <- object$x
  if (is.null(x_full) || length(cols_v) == 0L) {
    return(na_result)
  }

  keep <- setdiff(colnames(x_full), cols_v)
  x_reduced <- x_full[, keep, drop = FALSE]

  refit <- function(x) {
    tryCatch(
      fit_model(
        x = x,
        y = object$y,
        family = object$family,
        family_string = object$family_string,
        fitter = if (is.null(object$fitter)) "base" else object$fitter,
        weights = object$prior.weights,
        offset = object$offset,
        fast = TRUE
      ),
      error = function(e) NULL
    )
  }

  # Refit BOTH models with the same routine so the deviance scale matches.
  fit_full <- refit(x_full)
  fit_red  <- if (ncol(x_reduced) == 0L) {
    refit(x_reduced)
  } else {
    refit(x_reduced)
  }

  if (is.null(fit_full) || is.null(fit_red) ||
      is.null(fit_full$logl) || is.null(fit_red$logl)) {
    return(na_result)
  }

  lr <- 2 * (fit_full$logl - fit_red$logl)
  if (!is.finite(lr) || lr < 0) lr <- max(lr, 0)
  p <- stats::pchisq(lr, df = df_v, lower.tail = FALSE)

  list(lr = lr, df = df_v, p = p)
}

# ---------------------------------------------------------------------------
# Basis coefficients table (optional)
# ---------------------------------------------------------------------------

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

# Construct a readable transformation label for a single fitted column.
#
# Ordinary final FP coefficients multiply a shifted-but-unscaled basis, because
# fit_mfp() backscales the working predictors before the final transformation.
# ACD component columns are functions of A(x). New fits also use shifted,
# unscaled input for A(x); a stored non-unit ACD scale is honored only so formula
# output remains correct for legacy serialized model objects.
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

# Apply one FP power to a readable base expression.
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

# Inner expression for the final ordinary FP basis. The final fit uses shifted
# but unscaled predictors, so preprocessing scale is intentionally absent here.
mfp2_summary_fp_inner_expr <- function(v, shift) {
  if (!is.na(shift) && shift != 0) {
    sprintf("(%s + %s)", v, format(shift, trim = TRUE))
  } else {
    sprintf("(%s)", v)
  }
}

# Inner expression used only by the stored ACD approximation. New model objects
# store scale = 1 because ACD variables are not scaled. Retain support for a
# non-unit stored scale so summaries of legacy serialized objects remain exact.
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

# Look up shift and scale for a variable from the transformations table.
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

# Build explicit definitions for active ACD component columns. Both the ordinary
# FP term and new ACD fits use shifted, unscaled predictor values. A non-unit
# stored ACD scale is included only for compatibility with legacy fitted objects.
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

mfp2_summary_fit_stats <- function(object) {
  # Values used by the shared Model Fit block renderer. Everything else that
  # was previously in this list (LR test statistic, R-squared, and so on) has
  # been removed: the block is now purely descriptive, and inferential
  # statements are made only by the per-variable Nonlinear LRTs above.
  list(
    model_fit = mfp2_summary_model_fit_values(object)
  )
}

mfp2_summary_num <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NA_real_)
  suppressWarnings(as.numeric(x))
}

# Total FP-adjusted model degrees of freedom: the sum of df_final over the
# selected variables. Each variable is charged its selection-adjusted df
# (1 linear, 2 FP1, 4 FP2, etc.) rather than its raw coefficient count, so the
# overall LR test is consistent with the per-variable Nonlinear LRTs and with
# the MFP df accounting used elsewhere in the package. The intercept, which is
# already present in the null model, is not counted.
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

# Assemble the Model Fit rows (values only, no formatting) from a fitted mfp2
# object. Used by both print.mfp2() and print.summary.mfp2() so the two
# methods display exactly the same family-specific fit statistics.
#
# Returns a data frame with two rows (Full linear / MFP), and columns
# `fit_statistic` (numeric) and `df` (integer). The `statistic_label` attribute
# is "Deviance" for GLMs and "-2 log L" for Cox models. The Model Fit block
# renderer, mfp2_format_model_fit_block(), formats and prints it.
mfp2_summary_model_fit_values <- function(object) {
  family_string <- object$family_string
  fit_statistic <- c(
    mfp2_summary_num(object$linear_deviance),
    mfp2_summary_num(object$mfp_deviance)
  )

  # The df column follows a single convention: number of regression
  # coefficients, EXCLUDING the intercept. `linear_df` as stored comes from
  # fit_model()$df, which follows logLik.glm()'s convention -- it counts the
  # intercept for GLMs (Cox has none) and adds +1 for either the
  # Gaussian residual variance sigma^2 or negative-binomial theta. These
  # adjustments are stripped here so the number the
  # user sees matches the promise made in the note.
  linear_df <- if (!is.null(object$linear_df)) {
    d <- as.integer(object$linear_df)
    if (family_string %in% c("gaussian", "negbin")) {
      # Strip 1 for the intercept AND 1 for sigma^2 or theta.
      d <- d - 2L
    } else if (!identical(family_string, "cox")) {
      # Non-Gaussian GLM: strip the intercept only.
      d <- d - 1L
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
  attr(values, "statistic_label") <- if (identical(family_string, "cox")) {
    "-2 log L"
  } else {
    "Deviance"
  }
  values
}

# Print the Model Fit block: heading, table, and note. Shared by print.mfp2()
# and print.summary.mfp2(). Callers pass an already-computed values frame to
# avoid recomputing it.
#
# The `heading_printer` argument accepts a function of one string that draws
# the section heading in the caller's own style, so print.mfp2() can reuse
# its existing print_section_heading() helper (which draws the boxed rules
# used elsewhere in its output) without leaking that helper's internals into
# the summary path. summary()'s printer supplies its own rule-based heading.
mfp2_format_model_fit_block <- function(values, digits, heading_printer) {
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
  note_lines <- c(
    "df counts regression coefficients (excluding the intercept). For the MFP",
    "row it also includes 1 df for each estimated FP power: an FP1 term is",
    "charged 2 df total (1 coefficient + 1 power), an FP2 term 4 df (2",
    "coefficients + 2 powers)."
  )
  for (ln in note_lines) cat(ln, "\n", sep = "")
}

#' Print a Summary of an `mfp2` Model Fit
#'
#' Renders the structured summary produced by [summary.mfp2()].
#'
#' @param x An object of class \code{"summary.mfp2"}.
#' @param ... Not used.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso [summary.mfp2()], [print.mfp2()]
#'
#' @export
print.summary.mfp2 <- function(x, ...) {
  digits <- if (!is.null(x$digits)) x$digits else max(3L, getOption("digits") - 3L)
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

  # --- One-line meta -------------------------------------------------------
  meta <- sprintf(
    "Family: %s | Criterion: %s | Converged: %s",
    x$family, x$criterion, if (isTRUE(x$converged)) "yes" else "no"
  )
  cat(meta, "\n")
  if (!is.na(x$nevents)) {
    cat(sprintf("Observations: %s | Events: %s\n", x$n, x$nevents))
  } else {
    cat(sprintf("Observations: %s\n", x$n))
  }
  cat("\n")

  # --- Selection overview --------------------------------------------------
  section("Selection Overview")
  ft <- x$function_table
  # Selected first, then excluded.
  ft <- ft[order(ft$Selected != "yes"), , drop = FALSE]
  print.data.frame(ft, row.names = FALSE, right = FALSE)
  cat("\n")

  # --- Linear terms --------------------------------------------------------
  section("Linear Terms")
  if (nrow(x$linear_terms) == 0L) {
    cat("(none)\n\n")
  } else {
    is_gaussian <- identical(x$family, "gaussian")
    lt <- x$linear_terms
    fmt <- function(v, d = digits) formatC(v, format = "g", digits = d)

    disp <- data.frame(
      # Show the user-facing variable name, not the internal fitted-column
      # name. `lt$variable` is populated by mfp2_summary_linear_table() with
      # the ".N" suffix stripped from `lt$term`; the stored `term` remains
      # available for programmatic mapping back to coef(fit).
      Term = lt$variable,
      coef = fmt(lt$coef),
      `se(coef)` = fmt(lt$se),
      stat = fmt(lt$statistic),
      p = mfp2_summary_format_p(lt$p),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    stat_name <- if (is_gaussian) "t" else "z"
    names(disp)[names(disp) == "stat"] <- stat_name

    if (!is_gaussian && !is.null(lt$exp_coef)) {
      disp[["exp(coef)"]] <- fmt(lt$exp_coef)
      disp[["[95% CI]"]] <- sprintf("[%s, %s]", fmt(lt$ci_lower), fmt(lt$ci_upper))
    } else {
      disp[["[95% CI]"]] <- sprintf("[%s, %s]", fmt(lt$ci_lower), fmt(lt$ci_upper))
    }

    print.data.frame(disp, row.names = FALSE, right = FALSE)
    cat("\n")

    if (!is_gaussian) {
      exp_name <- switch(x$family,
                         "binomial" = "odds ratio",
                         "poisson"  = "rate ratio",
                         "cox"      = "hazard ratio",
                         "exp(coef)"
      )
      cat(sprintf(
        "exp(coef) is the %s (per unit for continuous terms; versus\nreference for factor levels).\n\n",
        exp_name
      ))
    } else {
      cat("Coefficients are on the outcome scale (mean difference).\n\n")
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
      `LR chi-sq` = formatC(nt$lr_chisq, format = "f", digits = 2),
      p = mfp2_summary_format_p(nt$p),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    print.data.frame(disp, row.names = FALSE, right = FALSE)
    cat("\n")
    cat(
      "Joint likelihood-ratio test per variable (H0: variable has no\n",
      "effect): the final MFP model versus the same model with that\n",
      "variable removed, all other functional forms held fixed. df uses\n",
      "the selection-adjusted convention (2 df per FP1, 4 per FP2).\n",
      "Diagnostic only; selection was performed by the MFP procedure.\n",
      "Use plot() to interpret the fitted function.\n\n",
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
        coef = formatC(bt$coef, format = "g", digits = digits),
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
  # fit statistic and df convention, together with the same note. The
  # `section()` helper defined above draws the dash-rule heading used elsewhere
  # in the summary output.
  mfp2_format_model_fit_block(
    values          = x$fit$model_fit,
    digits          = digits,
    heading_printer = section
  )

  cat("\n", rule_eq, "\n", sep = "")
  invisible(x)
}

# Format a p-value vector for table display.
mfp2_summary_format_p <- function(p) {
  vapply(p, function(v) {
    if (is.na(v)) return("NA")
    if (v < 0.0001) return("<0.001")
    formatC(v, format = "f", digits = 4)
  }, character(1L))
}