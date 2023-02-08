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

#' Helper function to print cycle header for verbose printing
print_mfp_summary <- function(criterion, ftest) {
  # the header for criterion = p value
  if (criterion == "pvalue") {
    if (ftest) {
      cat("\n", rep("-", 174), sep = "")
      pos <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21)
      print_structured("Variable", "Model", "(vs.)", "Deviance", "df.resid", "Dev.diff", "F", "df.num", "df.den", "P value", "Powers", pos = pos)
      cat("\n", rep("-", 174), sep = "")
    } else {
      cat("\n", rep("-", 126), sep = "")
      pos <- c(1, 3, 5, 7, 9, 11, 13, 15)
      print_structured("Variable", "Model", "(vs.)", "Deviance", "Dev diff.", "DF", "P value", "Powers", pos = pos)
      cat("\n", rep("-", 126), sep = "")
    }
    # the header for criterion = aic/bic
  } else {
    cat("\n", rep("-", 76), sep = "")
    pos <- c(1, 3, 5, 7) # output formatting
    print_structured("Variable", "Model", ifelse(criterion == "AIC", "AIC", "BIC"), "Powers", pos = pos)
    cat("\n", rep("-", 76), sep = "")
  }
}

# for printing when criterion is pvalue and test is chi-square
print_mfp_summary_1 <- function(namex, dev.all, dev.diff, pvalues, best.function, index.bestmodel, acd) {
  if (acd) {
    pos <- c(1, 3, 5, 7, 9, 11, 13, 15)
    print_structured(paste0("(A)", namex), " ", " ", pos = pos)
    mnames <- c(
      "NULL", paste0("Lin(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ")"),
      paste0("FP1(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ")"),
      paste0("FP1(a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"),
      paste0("FP1(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ",a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"),
      paste0("Lin(a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"), "Final"
    )
    # degrees of freedom
    df <- c(4, 3, 2, 2, 1)
    for (i in seq_along(dev.all)) {
      print_structured("        ", mnames[i],
        if (i == 1) {
          paste0(mnames[5])
        } else {
          if (i == 6) {
            mnames[4]
          } else {
            " "
          }
        },
        format(round(dev.all[i], 3), nsmall = 3),
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(dev.diff[5], 3), nsmall = 3)
          } else {
            format(round(dev.diff[i], 3), nsmall = 3)
          }
        },
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(df[5], 3), nsmall = 3)
          } else {
            format(round(df[i], 3), nsmall = 3)
          }
        },
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(pvalues[5], 3), nsmall = 3)
          } else {
            format(round(pvalues[i], 3), nsmall = 3)
          }
        },
        if (i == 1) {
          "."
        } else {
          best.function[[i]]
        },
        pos = pos
      )
    }
    print_structured(" ", paste0("Final:", mnames[index.bestmodel]), "", format(round(dev.all[index.bestmodel], 3), nsmall = 3), ".", ".", ".", if (index.bestmodel == 1) {
      "."
    } else {
      best.function[[index.bestmodel]]
    }, pos = pos)
  } else {
    pos <- c(1, 3, 5, 7, 9, 11, 13, 15)
    print_structured(namex, " ", " ", pos = pos)
    kk <- length(dev.all)
    m <- kk - 2 # remove null and lin to get the max permitted degree of FP
    mnames <- c("NULL", "Lin.", if (m != 0) {
      paste0("FP", seq_len(m))
    }, "Final")
    # degrees of freedom
    if (m != 0) {
      d <- seq(m - 1, 1)
      df <- c(2 * m, 2 * m - 1, if (m >= 2) {
        2 * (d)
      }) # same as 2*(m-d) where d = seq(1,m-1)
    } else {
      df <- c(1) # df of NULL vs linear
    }
    for (i in seq_along(dev.all)) {
      print_structured("        ", mnames[i],
        if (i == 1) {
          paste0(mnames[kk])
        } else {
          " "
        },
        format(round(dev.all[i], 3), nsmall = 3),
        if (i == kk) {
          "."
        } else {
          format(round(dev.diff[i], 3), nsmall = 3)
        },
        if (i == kk) {
          "."
        } else {
          format(round(df[i], 3), nsmall = 3)
        },
        if (i == kk) {
          "."
        } else {
          format(round(pvalues[i], 3), nsmall = 3)
        },
        if (i == 1) {
          "."
        } else {
          best.function[[i - 1]]
        },
        pos = pos
      )
    }
    print_structured("        ", paste0("Final:", mnames[index.bestmodel]), " ", format(round(dev.all[index.bestmodel], 3), nsmall = 3), ".", ".", ".", if (index.bestmodel == 1) {
      "."
    } else {
      best.function[[index.bestmodel - 1]]
    }, pos = pos)
  }
  cat("\n")
}
# for printing when criterion is pvalue and test is f test. best.function is only a list while others are vectors
print_mfp_summary_2 <- function(namex, dev.all, df.res, dev.diff, f, df.den, pvalues, best.function, index.bestmodel, pos = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21), acd) {
  if (acd) {
    print_structured(paste0("(A)", namex), " ", " ", pos = pos)
    mnames <- c(
      "NULL", paste0("Lin(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ")"),
      paste0("FP1(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ")"),
      paste0("FP1(a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"),
      paste0("FP1(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ",a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"),
      paste0("Lin(a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"), "Final"
    )
    # degrees of freedom
    df <- c(4, 3, 2, 2, 1)
    for (i in seq_along(dev.all)) {
      print_structured("        ", mnames[i],
        if (i == 1) {
          paste0(mnames[5])
        } else {
          if (i == 6) {
            mnames[4]
          } else {
            " "
          }
        },
        format(round(dev.all[i], 3), nsmall = 3),
        format(round(df.res[i], 3), nsmall = 3),
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(dev.diff[5], 3), nsmall = 3)
          } else {
            format(round(dev.diff[i], 3), nsmall = 3)
          }
        },
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(f[5], 3), nsmall = 3)
          } else {
            format(round(f[i], 3), nsmall = 3)
          }
        },
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(df[5], 3), nsmall = 3)
          } else {
            format(round(df[i], 3), nsmall = 3)
          }
        },
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(df.den[4], 3), nsmall = 3)
          } else {
            format(round(df.den[5], 3), nsmall = 3)
          }
        },
        if (i == 5) {
          "."
        } else {
          if (i == 6) {
            format(round(pvalues[5], 3), nsmall = 3)
          } else {
            format(round(pvalues[i], 3), nsmall = 3)
          }
        },
        if (i == 1) {
          "."
        } else {
          best.function[[i]]
        },
        pos = pos
      )
    }
    print_structured("        ", paste0("Final:", mnames[index.bestmodel]), " ", format(round(dev.all[index.bestmodel], 3), nsmall = 3), ".", ".", ".", ".", ".", ".",
      if (index.bestmodel == 1) {
        "."
      } else {
        best.function[[index.bestmodel]]
      },
      pos = pos
    )
  } else {
    print_structured(namex, " ", " ", pos = pos)
    kk <- length(dev.all)
    m <- kk - 2 # remove null and lin to get the max permitted degree of FP
    mnames <- c("NULL", "Lin.", if (m != 0) {
      paste0("FP", seq_len(m))
    }, "Final")
    # degrees of freedom
    if (m != 0) {
      d <- seq(m - 1, 1)
      df <- c(2 * m, 2 * m - 1, if (m >= 2) {
        2 * (d)
      }) # same as 2*(m-d) where d = seq(1,m-1)
    } else {
      df <- c(1) # df of NULL vs linear
    }

    for (i in seq_along(dev.all)) {
      print_structured("        ", mnames[i],
        if (i == 1) {
          paste0(mnames[kk])
        } else {
          " "
        },
        format(round(dev.all[i], 3), nsmall = 3),
        format(round(df.res[i], 3), nsmall = 3),
        if (i == kk) {
          "."
        } else {
          format(round(dev.diff[i], 3), nsmall = 3)
        },
        if (i == kk) {
          "."
        } else {
          format(round(f[i], 3), nsmall = 3)
        },
        if (i == kk) {
          "."
        } else {
          format(round(df[i], 3), nsmall = 3)
        },
        if (i == kk) {
          "."
        } else {
          format(round(df.den[kk], 3), nsmall = 3)
        },
        if (i == kk) {
          "."
        } else {
          format(round(pvalues[i], 3), nsmall = 3)
        }, if (i == 1) {
          "."
        } else {
          best.function[[i - 1]]
        },
        pos = pos
      )
    }
    print_structured("        ", paste0("Final:", mnames[index.bestmodel]), " ", format(round(dev.all[index.bestmodel], 3), nsmall = 3), ".", ".", ".", ".", ".", ".",
      if (index.bestmodel == 1) {
        "."
      } else {
        best.function[[index.bestmodel - 1]]
      },
      pos = pos
    )
  }
  cat("\n")
}

# for printing aic and bic
print_mfp_summary_3 <- function(namex, gic, best.function, keep, pos = c(1, 3, 5, 7), acd) {
  min.index <- which.min(unlist(gic))
  if (namex %in% keep) {
    min.index <- which.min(unlist(gic)[-1]) + 1
  }
  if (acd) {
    print_structured(paste0("(A)", namex), " ", " ", pos = pos)

    mnames <- c(
      "NULL", paste0("Lin(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ")"),
      paste0("FP1(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ")"),
      paste0("FP1(a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"),
      paste0("FP1(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), ",a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"),
      paste0("Lin(a(", ifelse(nchar(namex) > 3, paste0(substring(namex, 1, 3)), namex), "))"), "Final"
    )
    for (i in seq_along(gic)) {
      print_structured("        ", mnames[i], format(round(gic[i], 3)), if (i == 1) {
        "."
      } else {
        best.function[[i]]
      }, pos = pos)
    }
    print_structured("        ", paste0("Final:", mnames[min.index]), format(round(gic[min.index], 3)), if (min.index == 1) {
      "."
    } else {
      best.function[[min.index]]
    }, pos = pos)
  } else {
    print_structured(namex, " ", " ", pos = pos)
    nn <- length(gic) - 2 # remove null and lin
    mnames <- c("NULL", "Lin.", if (nn != 0) {
      paste0("FP", seq_len(nn))
    }, "Final")
    for (i in seq_along(gic)) {
      print_structured("        ", mnames[i], format(round(gic[i], 3), nsmall = 3), if (i == 1) {
        "."
      } else {
        best.function[[i - 1]]
      }, pos = pos)
    }
    print_structured("        ", paste0("Final:", mnames[min.index]), format(round(gic[min.index], 3), nsmall = 3), if (min.index == 1) {
      "."
    } else {
      best.function[[min.index - 1]]
    }, pos = pos)
  }
  cat("\n")
}

#' Helper to print information in a structured form
#' 
#' To be used in the `print_mfp_summary` functions.
print_structured <- function(..., pos) {

  x <- list(...)
  nx <- length(x)
  tabs <- 0
  cat("\n")
  for (i in 1:nx) {
    wordi <- x[[i]]
    if (pos[i] > tabs) {
      cat(rep("\t", (pos[i] - tabs)))
      tabs <- pos[i]
    }
    wordi <- wordi[!is.na(wordi)]
    # deal with powers
    if (length(wordi) > 0) {
      cat(wordi)
      if (length(wordi) > 1) {
        nwordi <- 1 + sum(nchar(wordi))
      } else {
        nwordi <- nchar(wordi)
      }
      if (nwordi > 7) {
        tabs <- tabs + 1
      }
    } else {
      cat("\t")
    }
  }
  
  invisible()
}
