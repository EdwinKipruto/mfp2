
#' @method print mfpa
#' @export
print.mfpa <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = FALSE, ...) {
  # shift and scaling factors with centering values
  cat("Shifting, Scaling and Centering of covariates", "\n")
  ss <- as.data.frame(x$transformations)
  ss[is.na(ss)] <- "."
  print.data.frame(ss)
  cat("\n")
  # Final MFP Powers
  cat("Final Multivariable Fractional Polynomial for y", "\n")
  dd <- as.data.frame(x$fp_terms)
  dd[is.na(dd)] <- "."
  print.data.frame(dd)
  cat("\n")
  family <- ifelse(is.character(x$family), x$family, x$family$family)
  if (family != "cox") {
    cat("\nCall:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = ""
    )
    if (length(coef(x))) {
      cat("Coefficients")
      if (is.character(co <- x$contrasts)) {
        cat(
          "  [contrasts: ",
          apply(cbind(names(co), co), 1L, paste, collapse = "="), "]"
        )
      }
      cat(":\n")
      print.default(format(x$coefficients, digits = digits),
        print.gap = 2, quote = FALSE
      )
    } else {
      cat("No coefficients\n\n")
    }
    cat(
      "\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
      x$df.residual, "Residual\n"
    )
    if (nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat(
      "Null Deviance:	   ", format(signif(x$null.deviance, digits)),
      "\nResidual Deviance:", format(signif(x$deviance, digits)),
      "\tAIC:", format(signif(x$aic, digits))
    )
    cat("\n")
    invisible(x)
    # cox model
  } else {
    if (!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }
    if (!is.null(x$fail)) {
      cat(" Coxph failed.", x$fail, "\n")
      return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- x$coefficients
    se <- sqrt(diag(x$var))
    if (is.null(coef) | is.null(se)) {
      stop("Input is not valid")
    }

    if (is.null(x$naive.var)) {
      tmp <- cbind(
        coef, exp(coef), se, coef / se,
        pchisq((coef / se)^2, 1, lower.tail = FALSE)
      )
      dimnames(tmp) <- list(names(coef), c(
        "coef", "exp(coef)",
        "se(coef)", "z", "p"
      ))
    } else {
      nse <- sqrt(diag(x$naive.var))
      tmp <- cbind(
        coef, exp(coef), nse, se, coef / se,
        pchisq((coef / se)^2, 1, lower.tail = FALSE)
      )
      dimnames(tmp) <- list(names(coef), c(
        "coef", "exp(coef)",
        "se(coef)", "robust se", "z", "p"
      ))
    }

    if (inherits(x, "coxphms")) {
      # print it group by group
      # lazy: I don't want to type x$cmap many times
      #  remove transitions with no covariates
      cmap <- x$cmap[, colSums(x$cmap) > 0, drop = FALSE]
      cname <- colnames(cmap)
      printed <- rep(FALSE, length(cname))
      for (i in 1:length(cname)) {
        # if multiple colums of tmat are identical, only print that
        #  set of coefficients once
        if (!printed[i]) {
          j <- apply(cmap, 2, function(x) all(x == cmap[, i]))
          printed[j] <- TRUE

          tmp2 <- tmp[cmap[, i], , drop = FALSE]
          names(dimnames(tmp2)) <- c(paste(cname[j], collapse = ", "), "")
          # restore character row names
          rownames(tmp2) <- rownames(cmap)[cmap[, i] > 0]
          printCoefmat(tmp2,
            digits = digits, P.values = TRUE,
            has.Pvalue = TRUE,
            signif.stars = signif.stars, ...
          )
          cat("\n")
        }
      }

      cat(" States: ", paste(paste(seq(along.with = x$states), x$states, sep = "= "),
        collapse = ", "
      ), "\n")
      # cat(" States: ", paste(x$states, collapse=", "), '\n')
      if (FALSE) { # alternate forms, still deciding which I like
        stemp <- x$states
        names(stemp) <- 1:length(stemp)
        print(stemp, quote = FALSE)
      }
    } else {
      printCoefmat(tmp,
        digits = digits, P.values = TRUE, has.Pvalue = TRUE,
        signif.stars = signif.stars, ...
      )
    }

    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) {
      df <- sum(!is.na(coef))
    } else {
      df <- round(sum(x$df), 2)
    }
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
      df, " df,", " p=",
      format.pval(pchisq(logtest, df, lower.tail = FALSE), digits = digits),
      "\n",
      sep = ""
    )
    omit <- x$na.action
    cat("n=", x$n)
    if (!is.null(x$nevent)) {
      cat(", number of events=", x$nevent, "\n")
    } else {
      cat("\n")
    }
    if (length(omit)) {
      cat("\   (", naprint(omit), ")\n", sep = "")
    }
    invisible(x)
  }
}
fpout <- function(x) print(x)
