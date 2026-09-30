#' Print intermediate MFP selection results
#'
#' Internal helper used when verbose output is requested during model
#' selection.
#'
#' @param xi Name of the variable currently being evaluated.
#' @param criterion Model-selection criterion.
#' @param fit Intermediate model-selection results.
#' @param stage2 Whether to print Stage 2 of the spike-at-zero procedure.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
print_mfp_step <- function(xi, criterion, fit, stage2 = FALSE) {
  
  print_mat_fct <- switch(
    criterion, 
    "pvalue" = print_mfp_pvalue_step, 
    print_mfp_ic_step
  ) 
  
  # matrix for printing
  # use whitespace in column names to try to make printing "fixed width"
  # by making the column names longer than its entries
  # longest name is linear(., A(x)) -> 15 symbols
  mat_print <- cbind(
    # remove NAs from printed powers
    "Powers   " = apply(fit$powers, 1, 
                        function(row) {
                          if (!fit$acd) {
                            row_prep = na.omit(row)
                            if (length(row_prep) == 0)
                              return("NA")  
                          } else row_prep = row
                          
                          paste0(row_prep, collapse = ", ")
                        }), 
    "DF   " = sprintf("%d", fit$metrics[, "df"]), 
    print_mat_fct(xi, fit, criterion)
  )

  # ensure fixed width model names
  rownames(mat_print) <- sprintf("%-16s", rownames(fit$metrics))
  
  # Joint AIC/BIC SAZ selection is a single comparison, not Stage 1 followed
  # by component removal. Keep its heading distinct so verbose output mirrors
  # the statistical procedure actually used.
  joint_saz_ic <- identical(fit$selection_mode, "joint_ic")

  # Stage 1 (p-value SAZ) or the only step (ordinary and joint-IC selection).
  if (!stage2) {
    active_flags <- c(
      if (isTRUE(fit$keep)) "keep = TRUE",
      if (isTRUE(fit$spike)) "spike = TRUE",
      if (isTRUE(fit$acd)) "acd = TRUE"
    )
    flag_suffix <- if (length(active_flags)) {
      sprintf(" (%s)", paste(active_flags, collapse = ", "))
    } else {
      ""
    }
    cat(sprintf("\nVariable: %s%s\n", xi, flag_suffix))
    if (joint_saz_ic) {
      cat("Joint Spike at Zero Selection:\n")
    } else if (fit$spike) {
      cat("Stage 1 of Spike at Zero Algorithm:\n")
    }
    print(mat_print, quote = FALSE, na.print = ".", print.gap = 1)
    selected_model <- rownames(fit$metrics)[fit$model_best]
    cat(sprintf("Selected: %s\n", selected_model))
  }
  
  # Print only stage 2 of SAZ if the variable was selected
  if (fit$spike && stage2) {
    selected_model <- rownames(fit$metrics)[fit$model_best]
    if (!selected_model %in% "null") {
    cat("Stage 2 of Spike at Zero Algorithm: \n")
    # Extract the dynamic row name
    model_name <- trimws(rownames(fit$metrics)[fit$model_best])
    
    # Ensure mat_print row names have no extra spaces
    rownames(mat_print) <- trimws(rownames(mat_print))
    
    # Find the row index dynamically
    row_index <- match(model_name, rownames(mat_print))
    mat_print <- cbind(
      # remove NAs from printed powers
      "Powers   " = c(mat_print[row_index, 1], mat_print[row_index, 1], 1), 
      "DF   " = fit$spike_metrics$metrics[,"df"], 
      print_mat_fct(xi, fit, criterion, spike = TRUE)
    )
    print(mat_print, quote = FALSE, na.print = ".", print.gap = 1)
    
    stage2_model_names <- rownames(fit$spike_metrics$metrics)
    
    selected_stage2_model <- saz_decision_label(
      fit$spike_metrics$spike_decision,
      style = "stage2_model",
      full_label = stage2_model_names[1],
      continuous_label = stage2_model_names[2],
      binary_label = stage2_model_names[3]
    )
    
    selected_stage2_row <- match(selected_stage2_model, rownames(mat_print))
    
    cat(sprintf("Selected: %s\n", rownames(mat_print)[selected_stage2_row]))
    } 
  }
}

#' Build the Printed Table for a P-Value MFP Selection Step
#'
#' Constructs the character matrix printed after one selection step under
#' `criterion = "pvalue"`. Rows correspond to the FP candidate forms
#' evaluated at that step (or the SAZ candidate models when `spike = TRUE`)
#' and columns hold the deviance and p-value comparison against the current
#' reference model.
#'
#' @param xi Character scalar naming the predictor being evaluated.
#' @param fit Step-selection result from the MFP engine.
#' @param criterion Character scalar; the active selection criterion,
#'   passed for consistent labelling.
#' @param spike Logical. If `TRUE`, format the output for a spike-at-zero
#'   selection step; otherwise format for a standard MFP step.
#'
#' @return A character matrix containing the printed selection table.
#'
#' @keywords internal
#' @noRd
print_mfp_pvalue_step <- function(xi, fit, criterion, spike = FALSE) {

  is_gee <- isTRUE(fit$is_gee)

  # Format to 3 decimals, keeping NA as a true NA so the caller's
  # `na.print = "."` renders a dot rather than the string "NA".
  fmt3 <- function(x) {
    x <- unname(x)
    ifelse(is.na(x), NA_character_, sprintf("%.3f", x))
  }

  # GEE p-value selection scores models by the negative overall robust Wald
  # chi-square, stored internally as deviance_rs = -W. That negative sign is an
  # internal convention for the shared "smaller-is-better" selection machinery;
  # for display we report the Wald chi-square W itself. A genuine -W is always
  # <= 0, so a candidate with deviance_rs > 0 (a divergent fit that fell back to
  # a losing sentinel) has no valid Wald and prints as NA (".").
  gee_wald <- function(dev) {
    ifelse(is.finite(dev) & dev <= 0, -dev, NA_real_)
  }

  if (spike) {
    m <- fit$spike_metrics$metrics
    fpmax_row <- rownames(m)[1]
    dev <- m[, "deviance_rs"]

    if (is_gee) {
      wald <- gee_wald(dev)
      # Statistic = W(full) - W(reduced); NA propagates when either lacks a Wald.
      wdiff <- c(wald[1] - wald[2], wald[1] - wald[3])
      mat <- cbind(
        "Wald chi-sq" = fmt3(wald),
        "Versus          " = c(NA, rep(fpmax_row, nrow(m) - 1)),
        "Wald diff." = c(NA, fmt3(wdiff)),
        "P-value" = c(NA, sprintf("%.3f", fit$spike_metrics$pvalue))
      )
      return(mat)
    }

    mat <- cbind(
      "Deviance   " = sprintf("%.3f", dev),
      "Versus          " = c(NA, rep(fpmax_row, nrow(m) - 1)),
      "Deviance diff." = c(NA, sprintf("%.3f", dev[2] - dev[1]),
                           sprintf("%.3f", dev[3] - dev[1])),
      "P-value" = c(NA, sprintf("%.3f", fit$spike_metrics$pvalue))
    )
    return(mat)
  }

  fpmax <- rownames(fit$metrics)[1]
  dev <- fit$metrics[, "deviance_rs"]

  if (is_gee) {
    wald <- gee_wald(dev)
    # Statistic for each reduced row vs fpmax = W(full) - W(reduced), matching
    # the sign convention of the deviance-difference test; NA propagates when
    # either model lacks a valid Wald.
    wdiff <- if (identical(fpmax, "null")) {
      wald[-1] - wald[1]
    } else {
      wald[1] - wald[-1]
    }
    mat <- cbind(
      "Wald chi-sq" = fmt3(wald),
      "Versus          " = c(NA, rep(fpmax, nrow(fit$metrics) - 1)),
      "Wald diff." = c(NA, fmt3(wdiff)),
      "P-value" = c(NA, sprintf("%.3f", fit$pvalue))
    )
    return(mat)
  }

  # p-value specific matrix for printing (likelihood families)
  mat <- cbind(
    "Deviance   " = sprintf("%.3f", dev),
    "Versus          " = c(NA, rep(fpmax, nrow(fit$metrics) - 1)),
    "Deviance diff." = c(NA, switch(fpmax,
                                    "null" = sprintf("%.3f",
                                                     dev[fpmax] - dev[-1]),
                                    sprintf("%.3f",
                                            dev[-1] - dev[fpmax]))),
    "P-value" = c(NA, sprintf("%.3f", fit$pvalue))
  )
  return(mat)
}

#' Build the Printed Table for an AIC/BIC MFP Selection Step
#'
#' Constructs the character matrix printed after one selection step under
#' `criterion = "aic"` or `"bic"`. Rows correspond to the FP candidate
#' forms evaluated at that step (or the SAZ candidate models when
#' `spike = TRUE`) and the sole numeric column holds the requested
#' information-criterion value.
#'
#' @param xi Character scalar naming the predictor being evaluated.
#' @param fit Step-selection result from the MFP engine.
#' @param criterion Character scalar; the active information criterion,
#'   either `"aic"` or `"bic"`. Used both to select the metric column and to
#'   label it in the printed output.
#' @param spike Logical. If `TRUE`, format the output for a spike-at-zero
#'   selection step; otherwise format for a standard MFP step.
#'
#' @return A character matrix containing the printed selection table.
#'
#' @keywords internal
#' @noRd
print_mfp_ic_step <- function(xi, fit, criterion, spike = FALSE) {
  
  if (spike){
    mat_print <- cbind(
      sprintf("%.3f", fit$spike_metrics$metrics[, tolower(criterion)])
    )
    colnames(mat_print) <- mfp_progress_criterion_label(fit, criterion)
    return(mat_print)
  }
  
  # IC specific matrix for printing
  mat_print <- cbind(
    sprintf("%.3f", fit$metrics[, tolower(criterion)])
  )
  colnames(mat_print) <- mfp_progress_criterion_label(fit, criterion)
  
  mat_print
}

#' Label an Information Criterion in Progress Output
#'
#' GEE models are scored with quasi-likelihood information criteria, so their
#' progress-table labels use QAIC and QBIC rather than the likelihood-based
#' AIC and BIC names used for other model families.
#'
#' @param fit Intermediate model-selection results.
#' @param criterion Character scalar naming the selection criterion.
#'
#' @return An uppercase character label for the progress output.
#'
#' @keywords internal
#' @noRd
mfp_progress_criterion_label <- function(fit, criterion) {
  label <- toupper(criterion)
  if (isTRUE(fit$is_gee) && label %in% c("AIC", "BIC")) {
    label <- paste0("Q", label)
  }
  label
}
