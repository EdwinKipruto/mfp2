#' Function for verbose printing of function selection procedure (FSP)
#' 
#' @param fit intermediary model fit in `mfp_step`.
#' @inheritParams find_best_fp_step
print_mfp_step <- function(xi, criterion, fit) {
  
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
  
  # print
  #cat(sprintf("\nVariable: %s (keep = %s)\n", xi, fit$keep))
  cat(sprintf(
  "\nVariable: %s (keep = %s, spike = %s)\n", 
  xi, fit$keep, fit$spike
  ))
  
  if (fit$spike) {
    cat("Stage 1 of Spike at Zero Algorithm: \n")
  }

  print(mat_print, quote = FALSE, na.print = ".", print.gap = 1)
  cat(sprintf("Selected: %s\n", rownames(fit$metrics)[fit$model_best]))
  
  if (fit$spike) {
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
