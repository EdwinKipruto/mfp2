print_mfp_step <- function(xi, criterion, fit) {
  
  print_fct <- switch(
    criterion, 
    "pvalue" = print_mfp_pvalue_step, 
    print_mfp_ic_step
  ) 
  
  cat(sprintf("\nVariable: %s (keep = %s)\n", xi, fit$keep))
  print_fct(xi, fit, criterion)
  cat(sprintf("Selected: %s\n", rownames(fit$metrics)[fit$model_best]))
}

#' Function for verbose printing of FSP based on p-value
#' 
#' @details 
#' Prints out details on the variabe
print_mfp_pvalue_step <- function(xi, fit, ...) {
  
  fpmax <- rownames(fit$metrics)[1]
  
  # matrix for printing
  # use whitespace in column names to try to make printing "fixed width"
  # by making the column names longer than its entries
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
    "Deviance   " = sprintf("%.1f", fit$metrics[, "deviance_rs"]), 
    "Versus        " = c(NA, rep(fpmax, nrow(fit$metrics) - 1)), 
    "Deviance diff." = c(NA, 
                         sprintf("%.1f", 
                                 fit$metrics[-1, "deviance_rs"] - 
                                   fit$metrics[fpmax, "deviance_rs"])), 
    "P-value" = c(NA, sprintf("%.4f", fit$pvalue))
  )
  # ensure fixed width model names
  # longest name is FP1(x, A(x)) -> 12 symbols
  rownames(mat_print) <- sprintf("%-14s", rownames(fit$metrics))
  
  print(mat_print, quote = FALSE, na.print = ".", print.gap = 1)
}

print_mfp_ic_step <- function(xi, fit, criterion, ...) {
  
}
