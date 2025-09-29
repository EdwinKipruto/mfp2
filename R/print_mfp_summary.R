#' Verbose printing of the function selection procedure (FSP)
#'
#' This function prints intermediate results from the fractional polynomial
#' function selection process, including details of the spike-at-zero algorithm
#' when relevant. It is primarily called from within `mfp_step` during model
#' selection. 
#'
#' @param xi Character string. The name of the covariate currently being tested.
#' @param criterion Character string. The selection criterion, either `"pvalue"`
#'   or `"aic"`/`"bic"` (depending on implementation in
#'   `print_mfp_ic_step()`).
#' @param fit A list containing the intermediate fit for the current variable.
#'   Typically produced by `mfp_step`.
#' @param stage2 Logical. When `FALSE` (default), only Stage 1 results of the 
#'   spike-at-zero (SAZ) procedure are printed, or the usual FSP output for 
#'   `xi` when `fit$spike[xi] = FALSE`. When `TRUE`, Stage 1 results are printed, 
#'   and if `fit$spike[xi] = TRUE` and the Stage 1 selected model is not `"null"`, 
#'   the Stage 2 spike-at-zero results are printed immediately after Stage 1.
#'   
#' @details
#' Printing is split into two parts:
#' - **Stage 1**: All candidate FP/ACD models for `xi` are listed with their powers,
#'   degrees of freedom, and criterion-specific statistics. If `fit$spike[xi] = TRUE`,
#'   this corresponds to the first stage of the spike-at-zero (SAZ) algorithm.
#' - **Stage 2**: If `stage2 = TRUE` and the Stage 1 selected model for `xi` is
#'   not `"null"`, the results of the spike-at-zero are printed immediately
#'   after Stage 1. The chosen model is highlighted at the end of the table.
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
  
  # Stage 1: only if not in Stage 2 mode
  if (!stage2) {
    cat(sprintf(
      "\nVariable: %s (keep = %s, spike = %s)\n", 
      xi, fit$keep, fit$spike
    ))
    if (fit$spike) {
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
    cat(sprintf("Selected: %s\n", rownames(mat_print)[fit$spike_metrics$spike_decision]))
    } 
  }
}

#' @param spike Logical. Indicates whether `xi` was treated as a spike-at-zero 
#' variable. Default is `FALSE`.
#' @describeIn print_mfp_step Helper for verbose printing based on p-value.
print_mfp_pvalue_step <- function(xi, fit, criterion, spike = FALSE) {
  
  fpmax <- rownames(fit$metrics)[1]
  
  if (spike) {
    mat <- cbind(
      "Deviance   " = sprintf("%.1f", fit$spike_metrics$metrics[, "deviance_rs"]),
      #"DF" = fit$spike_metrics$metrics[, "df"],
      "Versus          " = c(NA, rep(rownames(fit$spike_metrics$metrics)[1], nrow(fit$spike_metrics$metrics) - 1)), 
      "Deviance diff." = c(NA, sprintf("%.1f", 
                                       fit$spike_metrics$metrics[2, "deviance_rs"] - 
                                         fit$spike_metrics$metrics[1, "deviance_rs"]),
                           sprintf("%.1f", 
                                   fit$spike_metrics$metrics[3, "deviance_rs"] - 
                                     fit$spike_metrics$metrics[1, "deviance_rs"]))
      , 
      "P-value" = c(NA, sprintf("%.4f", fit$spike_metrics$pvalue))
    )
    return(mat)
  }
  
  # p-value specific matrix for printing
  mat <- cbind(
    "Deviance   " = sprintf("%.1f", fit$metrics[, "deviance_rs"]), 
    "Versus          " = c(NA, rep(fpmax, nrow(fit$metrics) - 1)), 
    "Deviance diff." = c(NA, switch(fpmax,
                                    "null" = sprintf("%.1f",
                                                     fit$metrics[fpmax, "deviance_rs"]-
                                                     fit$metrics[-1, "deviance_rs"] 
                                                       ),
                                    sprintf("%.1f", 
                                            fit$metrics[-1, "deviance_rs"] - 
                                              fit$metrics[fpmax, "deviance_rs"]))
                                    ), 
    "P-value" = c(NA, sprintf("%.4f", fit$pvalue))
  )
  return(mat)
}

#' @describeIn print_mfp_step Helper for verbose printing based on information criterion.
print_mfp_ic_step <- function(xi, fit, criterion, spike = FALSE) {
  
  if (spike){
    mat_print <- cbind(
      sprintf("%.1f", fit$spike_metrics$metrics[, tolower(criterion)])
    )
    colnames(mat_print) <- toupper(criterion) 
    return(mat_print)
  }
  
  # IC specific matrix for printing
  mat_print <- cbind(
    sprintf("%.1f", fit$metrics[, tolower(criterion)])
  )
  colnames(mat_print) <- toupper(criterion)
  
  mat_print
}
