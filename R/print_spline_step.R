#' Print one spline selection step during verbose fitting
#'
#' Uses the same total-model DF as the FP progress tables. For p-value
#' selection, each reduced form is shown against the maximum spline. A test
#' that was not reached by the selection procedure has no printed p-value.
#' @param xi Name of the predictor being evaluated.
#' @param criterion Selection criterion: p-value, AIC, or BIC.
#' @param fit Result returned by `select_spline_step()`.
#' @param spec Fitted RCS or cubic B-spline specification.
#' @param family_string Family identifier, used for GEE progress labels.
#' @return Invisibly returns `NULL` after printing the candidate table.
#' @keywords internal
#' @noRd
print_spline_step <- function(xi, criterion, fit, spec, family_string) {
  criterion <- tolower(criterion)
  is_gee <- mfp2_family_is_gee(family_string)
  form <- if (identical(spec$type, "bspline")) "bs" else "rcs"
  keys <- fit$candidate_keys
  metric <- fit$metrics
  stopifnot(length(keys) == nrow(metric))

  # Sort by the term's complexity, while DF below remains the total fitted
  # model DF. Indicator-only and combined forms get their own rows.
  complexity <- vapply(keys, function(key) {
    if (key == "null") return(-1)
    if (key == "binary") return(0)
    as.numeric(sub("[BC]$", "", key)) + if (endsWith(key, "B")) 0.1 else 0
  }, numeric(1L))
  positions <- order(complexity, decreasing = TRUE)
  keys <- keys[positions]
  metric <- metric[positions, , drop = FALSE]

  label <- function(key) {
    if (key == "null") return("null")
    if (key == "binary") return("zero indicator only")
    df <- as.integer(sub("[BC]$", "", key))
    name <- if (df == 1L) "linear" else sprintf("%s(df=%d)", form, df)
    if (endsWith(key, "B")) paste0(name, " + zero indicator") else name
  }
  format_metric <- function(value) {
    ifelse(is.finite(value), sprintf("%.3f", value), ".")
  }
  format_pvalue <- function(value) {
    if (!length(value) || !is.finite(value)) return(".")
    if (value < 0.001) "<0.001" else sprintf("%.3f", value)
  }

  output <- data.frame(
    DF = sprintf("%g", metric[, "df"]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  deviance <- metric[, "deviance_rs"]
  if (is_gee && criterion == "pvalue") {
    output[["Wald chi-sq"]] <- format_metric(ifelse(deviance <= 0, -deviance, NA_real_))
  } else {
    output[["Deviance"]] <- format_metric(deviance)
  }

  if (criterion == "pvalue") {
    full_key <- paste0(spec$df, if (isTRUE(spec$indicator)) "B" else "C")
    comparison <- difference <- pvalue <- rep(".", length(keys))
    for (i in seq_along(keys)) {
      key <- keys[[i]]
      if (key == full_key) next
      reference_key <- full_key
      test_name <- if (key == "null") {
        "maximum vs null"
      } else if (key == "binary") {
        "drop spline"
      } else {
        df <- as.integer(sub("[BC]$", "", key))
        matching_indicator <- endsWith(key, "B") == isTRUE(spec$indicator)
        if (!matching_indicator) NA_character_ else if (df == 1L) {
          "maximum vs linear"
        } else sprintf("maximum vs df %d", df)
      }
      if (key == "binary" && "drop spline" %in% names(fit$pvalue)) {
        reference_key <- fit$stage1_key
      } else if (endsWith(fit$stage1_key, "B") &&
                 key == sub("B$", "C", fit$stage1_key) &&
                 "drop indicator" %in% names(fit$pvalue)) {
        reference_key <- fit$stage1_key
        test_name <- "drop indicator"
      }
      ref <- match(reference_key, keys)
      comparison[[i]] <- label(reference_key)
      # deviance_rs = -W for the GEE p-value selector, so the difference
      # still follows the comparison's full-minus-reduced Wald convention.
      difference[[i]] <- format_metric(deviance[[i]] - deviance[[ref]])
      if (!is.na(test_name) && test_name %in% names(fit$pvalue)) {
        pvalue[[i]] <- format_pvalue(fit$pvalue[[test_name]])
      }
    }
    output[["Versus"]] <- comparison
    output[[if (is_gee) "Wald diff." else "Deviance diff."]] <- difference
    output[["P-value"]] <- pvalue
  } else {
    output[[mfp_progress_criterion_label(list(is_gee = is_gee), criterion)]] <-
      format_metric(metric[, criterion])
  }

  cat(sprintf("\nVariable: %s\n", xi))
  # Keep model labels in the unheaded first column and align all fields.
  columns <- c(list(c("", vapply(keys, label, character(1L)))),
               lapply(names(output), function(name) c(name, output[[name]])))
  widths <- vapply(columns, function(column) max(nchar(column, type = "width")),
                    integer(1L))
  widths[[1L]] <- max(16L, widths[[1L]])
  for (i in seq_len(nrow(output) + 1L)) {
    fields <- vapply(seq_along(columns), function(j) {
      sprintf("%-*s", widths[[j]], columns[[j]][[i]])
    }, character(1L))
    cat(sub(" +$", "", paste(fields, collapse = "  ")), "\n", sep = "")
  }
  cat(sprintf("Selected: %s\n", label(fit$model_best)))
  invisible(NULL)
}
