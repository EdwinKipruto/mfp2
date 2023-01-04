#' Function to extract and transform adjustment variables
#' 
#' @details
#' After extracting the adjustment variables this function, using their
#' corresponding FP powers stored in allpowers, transforms them. This function is
#' critical because, when evaluating x of interest, we must account for other
#' variables, which can be transformed or untransformed, depending on the
#' individual powers. It's worth noting that some powers can be NA, indicating
#' that the variable has been left out of the adjustment variables. It also
#' returns the FP data, which is dependent on the degrees of freedom. For example,
#' df = 2 is equivalent to FP degree one, resulting in the generation of 8 
#' variables. If acdx is set to TRUE, however, 64 variables are generated.
#' For more information, see Royston (2016).
#' 
#' @param x a matrix of predictors that includes the xi variable of interest. 
#' It is assumed that continuous variables have already been shifted and scaled.
#' @param xi name of the continuous predictor for which the FP function will be
#' estimated. There are no binary or two-level variables allowed.
#' @param allpowers a named list of FP powers of all variables of interest, 
#' including xi. Note that these powers are updated during backfitting or MFP 
#' cycles.
#' @param df vector of degrees of freedom in which 1 indicates linear, 
#' 2 indicates FP1, 4 indicates FP2, and so on.
#' @param powers a set of FP powers.
#' 
#' @return 
#' A list containing the following elements: FP powers ("powers"), 
#' FP data ("fpdata"), adjustment data ("adjdata"), and adjustment FP powers 
#' (adjustpowers).
extract_adjustment_data <- function(x,
                                    xi,
                                    allpowers,
                                    df, 
                                    powers,
                                    acdx) {
  # number of observations
  N <- dim(x)[1L]
  # Sort x based on the names of the powers
  pwrs.namex <- names(allpowers)
  x <- x[, pwrs.namex, drop = F]
  # Extract the matrix of adjustment variables
  adjvars <- pwrs.namex[!(pwrs.namex %in% xi)]
  xadjust <- x[, adjvars, drop = F]
  # FP Powers of adjustment variables
  adjustpowers <- allpowers[adjvars]
  # acd of adjustment variables. acdx is a named vector of TRUE or FALSE
  acdxa <- unname(acdx[adjvars])
  # check whether all adjustment powers = NA. Meaning all adjustment variables
  # were eliminated in MFP backfitting process
  if (!all(is.na(unlist(adjustpowers, use.names = F)))) {
    nv <- ncol(xadjust)
    xadjust.T <- vector(mode = "list", length = nv)
    for (i in 1:nv) {
      # we distinguish between acd and non-acd variables. To avoid treating acd
      # variables as FP2 variables. we use a special function for acd
      if (acdxa[i]) {
        # scale = 1 and shift = 0 means no scaling and shifting
        xadjust.T[[i]] <- generate_acd_fp(
          x = xadjust[, i, drop = TRUE], power = adjustpowers[[i]],
          scale = 1, shift = 0, s = powers
        )
      } else {
        xadjust.T[[i]] <- generate_fp(
          x = xadjust[, i, drop = TRUE], power = adjustpowers[[i]],
          scale = 1, shift = 0
        )
      }
    }
    # combine the transformed adjustment varaibles
    xadjust.T <- do.call(cbind, xadjust.T)
    # assign arbitrary names to adjustment matrix. After all we are not interested
    # with the names of adjustment variables at this stage
    colnames(xadjust.T) <- paste0("adjvar", 1:ncol(xadjust.T))
  } else {
    xadjust.T <- NULL
    adjustpowers <- NULL
  }
  # =============================================================================
  # Generate FP data:
  # =============================================================================
  xinterest <- x[, xi, drop = T]
  nxx <- length(unique(xinterest))
  # if a variable is binary or has less than 4 levels we do not generate FP data.
  # Similarly when df = 1 we do not generate FP data because we assume linearity
  if (nxx <= 3) {
    fpdata <- xinterest
  } else {
    if (df == 1) {
      fpdata <- xinterest
    } else {
      # Generate fpdata if the variable has more than two levels or df>1
      # If a variable is acd we generate 64 pairs of variables
      if (acdx[xi]) {
        fpd <- generate_acd_fp_data(xinterest, powers = powers)
        fpdata <- fpd$data
        powers <- fpd$powers
      } else {
        # Generate 8 variables if degree = 1, 36 if degree = 2 and so on
        fpd <- generate_fp_data(xinterest, degree = df / 2, powers = powers)
        fpdata <- fpd$data
        powers <- fpd$powers
      }
    }
  }
  return(list(powers = powers, fpdata = fpdata, adjdata = xadjust.T, adjustpowers = adjustpowers))
}

#' Helper to generate FPa data for variable of interest
#' 
#' Utilizes acd(x) transformation. To be used in [extract_adjustment_data()].
#' 
#' @param x a vector of a predictor for which FPs is required. x is assumed to 
#' be shifted and scaled appropriately. 
#' @param powers vector of FP powers.
generate_acd_fp_data <- function(x, powers) {
  # Possible combination of powers given degree = 2 i.e p1,p2
  # we set degree = 2 because we are interested in two powers for x and ax = acd(x)
  combs <- generate_fp_powers(degree = 2, powers = powers)
  # interchange the columns of combs, select unrepeated powers and rbind
  combs1 <- combs[, 2:1]
  keep <- apply(combs1, 1, function(x) length(unique(x[!is.na(x)])) != 1)
  # Final pairs of powers.if default FP set is used then its length is 64
  combs2 <- rbind(combs, combs1[keep, ])
  nfp <- dim(combs2)[1L]
  # Save FP transformed data as list
  fpdt <- vector(mode = "list", length = nfp)
  for (i in seq_len(nfp)) {
    fpdt[[i]] <- generate_acd_fp(
      x = x, power = combs2[i, ], scale = 1,
      shift = 0, center = F
    )
  }
  
  list(data = fpdt, powers = combs2)
}
