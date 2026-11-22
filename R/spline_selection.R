#' Select a nested spline term and its optional exact-zero indicator
#'
#' Candidate fits reuse the precomputed basis in `x`. RCS prefixes use one
#' restricted cubic basis; B-spline prefixes combine columns from nested
#' fixed-knot B-spline spaces. The mapped block controls the selected prefix
#' and whether an exact-zero indicator accompanies it.
#' Under AIC/BIC, failed optional candidates are skipped with a warning.
#' Null fits and all closed-testing fits remain required and fail strictly.
#'
#' @param args Arguments prepared by `find_best_fp_step()`.
#' @param max_df Maximum number of spline basis columns.
#' @param candidate_df Increasing spline column counts to compare, including
#'   the linear and full states. B-splines omit the knotless cubic state.
#' @return Selected spline df (NA for removal, zero for indicator only), the
#'   spike decision, model metrics, and adjustment cache metadata.
#' @keywords internal
#' @noRd
select_spline_step <- function(args, max_df, candidate_df = seq_len(max_df)) {
  xi <- args$xi
  full_columns <- args$term_to_columns[[xi]]
  spec <- attr(args$term_to_columns, "spline_specs", exact = TRUE)[[xi]]
  has_indicator <- isTRUE(spec$indicator)
  spike <- isTRUE(spec$spike)
  criterion <- tolower(args$criterion)
  is_ic <- criterion %in% c("aic", "bic")
  fits <- list()
  failed_candidates <- character(0L)

  # The adjustment model is fixed throughout this focal step. Reuse both
  # its design and null fit across all continuous/indicator candidates.
  adjustment_args <- args[c(
    "x", "xi", "powers_current", "powers", "acdx", "zero", "catzero",
    "spike", "spike_decision", "acd_parameter", "prev_adj_params",
    "transform_cache", "term_to_columns"
  )]
  precomputed_adj <- do.call(build_adjustment_step, adjustment_args)
  fit_args <- args[setdiff(names(args), c(
    "keep", "degree", "criterion", "ftest", "select", "alpha",
    "force_max_fp", "transform_cache"
  ))]
  fit_args$precomputed_adj <- precomputed_adj
  if (mfp2_family_is_gee(args$family_string)) {
    fit_args$gee_selection_criterion <- criterion
  }
  fit_null <- do.call(fit_null_step, fit_args)

  fit_candidate <- function(df, indicator) {
    key <- if (df == 0L) "binary" else paste0(df, if (indicator) "B" else "C")
    columns <- if (df == 0L) character(0L) else full_columns[seq_len(df)]
    if (indicator) columns <- c(columns, full_columns[[max_df + 1L]])
    candidate <- fit_args
    candidate$term_to_columns[[xi]] <- columns
    attr(candidate$term_to_columns, "spline_specs") <-
      attr(args$term_to_columns, "spline_specs", exact = TRUE)
    # The grouped linear fitter scores an arbitrary fixed design block.
    # An IC candidate is optional, including the maximum form. Closed tests
    # require their specified reference and reduced fits, so do not skip them.
    fitted <- if (is_ic) {
      tryCatch(
        do.call(fit_linear_step, candidate),
        mfp2_candidate_failure = function(e) {
          failed_candidates[key] <<- conditionMessage(e)
          NULL
        }
      )
    } else {
      do.call(fit_linear_step, candidate)
    }
    if (is.null(fitted)) return(key)
    fits[[key]] <<- list(
      metrics = rbind(fit_null$metrics, fitted$metrics),
      current_adj_params = fitted$current_adj_params
    )
    key
  }

  with_indicator <- has_indicator
  full_key <- fit_candidate(max_df, with_indicator)
  for (df in setdiff(candidate_df, max_df)) fit_candidate(df, with_indicator)
  if (spike) {
    for (df in candidate_df) fit_candidate(df, FALSE)
    fit_candidate(0L, TRUE)
  }

  if (length(failed_candidates)) {
    warning(structure(
      list(
        message = sprintf(
          "i %d spline candidate model(s) for '%s' could not be fitted and were skipped (%s).\ni First problem: %s",
          length(failed_candidates), xi,
          paste(names(failed_candidates), collapse = ", "),
          failed_candidates[[1L]]
        ),
        call = NULL
      ),
      class = c("mfp2_skipped_candidates", "warning", "condition")
    ))
  }

  metric <- function(key) fits[[key]]$metrics[2L, , drop = FALSE]
  null <- fit_null$metrics[1L, , drop = FALSE]
  compare <- function(reduced, larger) {
    pair <- rbind(reduced, larger)
    if (isTRUE(args$ftest)) {
      calculate_f_test(pair[, "deviance_gaussian"], pair[, "df_resid"], args$n_obs)
    } else if (mfp2_family_is_gee(args$family_string)) {
      calculate_deviance_test(pair[, "deviance_rs"], pair[, "df"])
    } else {
      calculate_lr_test(pair[, "logl"], pair[, "df"])
    }
  }
  pvalues <- statistics <- numeric(0L)
  test <- function(name, reduced, larger) {
    result <- compare(reduced, larger)
    pvalues[name] <<- result$pvalue
    statistics[name] <<- result$statistic
    result$pvalue
  }

  selected <- max_df
  chosen_key <- full_key
  stage1_key <- full_key
  decision <- saz_decision_codes[[if (spike) "cont_binary" else "continuous_only"]]
  if (criterion == "pvalue") {
    full <- metric(full_key)
    p <- test("maximum vs null", null, full)
    if (mfp_pvalue_exceeds(p, args$select) && !xi %in% args$keep) {
      selected <- NA_integer_
      chosen_key <- "null"
    } else {
      label <- function(df) paste0(df, if (with_indicator) "B" else "C")
      p <- test("maximum vs linear", metric(label(1L)), full)
      if (mfp_pvalue_exceeds(p, args$alpha)) {
        selected <- 1L
      } else if (length(candidate_df) > 2L) {
        for (df in setdiff(candidate_df, c(1L, max_df))) {
          p <- test(paste0("maximum vs df ", df), metric(label(df)), full)
          if (mfp_pvalue_exceeds(p, args$alpha)) {
            selected <- df
            break
          }
        }
      }
      chosen_key <- label(selected)
      stage1_key <- chosen_key
      if (spike) {
        both <- metric(chosen_key)
        continuous_key <- paste0(selected, "C")
        p_binary <- test("drop indicator", metric(continuous_key), both)
        p_continuous <- test("drop spline", metric("binary"), both)
        if (mfp_pvalue_exceeds(p_binary, args$alpha) &&
            !mfp_pvalue_exceeds(p_continuous, args$alpha)) {
          decision <- saz_decision_codes[["continuous_only"]]
          chosen_key <- continuous_key
        } else if (!mfp_pvalue_exceeds(p_binary, args$alpha) &&
                   mfp_pvalue_exceeds(p_continuous, args$alpha)) {
          decision <- saz_decision_codes[["binary_only"]]
          chosen_key <- "binary"
          selected <- 0L
        } else if (mfp_pvalue_exceeds(p_binary, args$alpha) &&
                   mfp_pvalue_exceeds(p_continuous, args$alpha)) {
          if (metric(continuous_key)[1L, "bic"] < metric("binary")[1L, "bic"]) {
            decision <- saz_decision_codes[["continuous_only"]]
            chosen_key <- continuous_key
          } else {
            decision <- saz_decision_codes[["binary_only"]]
            chosen_key <- "binary"
            selected <- 0L
          }
        }
      }
    }
  } else {
    keys <- names(fits)
    scores <- vapply(keys, function(key) metric(key)[1L, criterion], numeric(1L))
    scores <- c(null = unname(null[1L, criterion]), scores)
    if (xi %in% args$keep) scores[["null"]] <- Inf
    if (!any(is.finite(scores))) {
      stop("No finite ", toupper(criterion), " spline candidates are available for '",
           xi, "' under the requested keep constraint.", call. = FALSE)
    }
    scores[!is.finite(scores)] <- Inf
    chosen_key <- names(scores)[which.min(scores)]
    if (chosen_key == "null") selected <- NA_integer_
    else if (chosen_key == "binary") selected <- 0L
    else selected <- as.integer(sub("[BC]$", "", chosen_key))
    if (spike) decision <- if (chosen_key == "binary") {
      saz_decision_codes[["binary_only"]]
    } else if (grepl("B$", chosen_key)) {
      saz_decision_codes[["cont_binary"]]
    } else {
      saz_decision_codes[["continuous_only"]]
    }
  }
  spike_decision <- args$spike_decision
  spike_decision[[xi]] <- decision
  list(
    power_best = selected, searched_fp = FALSE,
    spike_decision = spike_decision,
    current_adj_params = if (chosen_key == "null") fit_null$current_adj_params else
      fits[[chosen_key]]$current_adj_params,
    transform_cache = precomputed_adj$transform_cache,
    metrics = do.call(rbind, c(list(null = null), lapply(fits, function(f) f$metrics[2L, , drop = FALSE]))),
    model_best = chosen_key, stage1_key = stage1_key,
    candidate_keys = c("null", names(fits)),
    pvalue = pvalues, statistic = statistics
  )
}

