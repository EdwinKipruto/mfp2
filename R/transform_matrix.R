#' Function to transform each column of matrix using final FP powers or acd
#' 
#' @param x a matrix with all continuous variables shifted and scaled.
#' @param power.list a named list of FP powers. 
#' @param center a logical vector specifying whether the columns in `x` should 
#' be centered.
#' 
#' @return 
#' If all elements of `power.list` are `NA` then this function returns `NULL`.
#' Otherwise a matrix is returned with transformed data, possibly with more 
#' columns than the input matrix due to higher order FP transformations.
transform_matrix <- function(x,
                             power.list, 
                             center, 
                             acdx) {
  # check whether all power.list are equal to NA-all variables were eliminated
  if (all(is.na(unlist(power.list)))) {
    # return NULL
    xtr.out <- NULL
  } else {
    namx <- names(power.list)
    x <- x[, namx, drop = F]
    # subset scale, center and shift using names of powers selected-
    # scale <- scale[namx] # we assume x has been shifted and scaled
    center <- center[namx]
    acdx <- acdx[namx]
    # transform acd variables if any
    if (any(acdx)) {
      # select FP powers for acd variables
      pow.acd <- power.list[acdx]
      # check whether the acd variables were selected
      if (all(is.na(unlist(pow.acd)))) {
        xtr1 <- NULL
      } else {
        xnames1 <- names(pow.acd)
        x.acdx <- x[, xnames1, drop = F]
        xtrans.acd <- vector(mode = "list", length = length(pow.acd))
        for (i in seq_along(xnames1)) {
          xtrans.acd[[i]] <- transform_acd_vector(
            x = x.acdx[, i, drop = T],
            power = pow.acd[[i]], scale = 1, # scale = scale[i],
            center = center[xnames1][i], shift = NULL
          ) # shift = shift[i])
        }
        # cbind acd transformed variables
        xtr1 <- do.call(cbind, xtrans.acd)
        # rename acd transformed variables. they start with letter "A"
        snames <- unlist(sapply(xnames1, function(x) paste0(c(x, paste0("(A)", x)), ".", c(1, 1)),
          simplify = F, USE.NAMES = F
        ))
        colnames(xtr1) <- snames[which(!is.na(unlist(pow.acd)))]
      }
      # Variables without acd transformations
      pow <- power.list[!acdx]
      if (all(is.na(unlist(pow)))) {
        xtr2 <- NULL
      } else {
        # get rid of unselected variables
        # transform_fp_vector returns NULL when power = NA
        pow2 <- pow[!is.na(pow)]
        xnames2 <- names(pow2)
        xx <- x[, xnames2, drop = F]
        xtrans <- vector(mode = "list", length = length(pow2))
        for (i in seq_along(xnames2)) {
          # x is already scaled and shifted so transform_fp_vector
          # shifting is irrelevant
          xtrans[[i]] <- transform_fp_vector(
            x = xx[, i, drop = T], power = pow2[[i]], scale = 1, # scale = scale[i],
            center = center[xnames2][i], shift = NULL
          ) # shift = shift[i])
        }
        # cbind non-acd transformed variables-
        xtr2 <- do.call(cbind, xtrans)
        # rename xtr such that if x1 is fp2 then we have x1.1, x1.2 whereas fp1 is x1.1
        # we can as well use make.names(rep(names(fpp), lapply(fpp, length)), sep = ".")
        # but has undesirable names like x1, x1.1 instead of x1.1, x1.2
        colnames(xtr2) <- unlist(sapply(xnames2, function(x) paste0(x, ".", seq_along(pow2[[x]])),
          simplify = F, USE.NAMES = F
        ))
      }
      # combine xtr1 and xtr2
      xtr.out <- cbind(xtr1, xtr2)
      # Usual MFP without acd variables
    } else {
      # # Get rid of unselected variables denoted by NA in power.list
      fpp <- power.list[!is.na(power.list)]
      namxx <- names(fpp)
      # subset x
      x <- x[, namxx, drop = F]
      xtransx <- vector(mode = "list", length = length(fpp))
      for (i in seq_along(namxx)) {
        # x is already scaled and shifted so transform_fp_vector shifting is irrelevant
        xtransx[[i]] <- transform_fp_vector(
          x = x[, i, drop = T], power = fpp[[i]], scale = 1, # scale = scale[i],
          center = center[namxx][i], shift = 0
        ) # shift = shift[i])
      }
      xtr.out <- do.call(cbind, xtransx)
      colnames(xtr.out) <- unlist(sapply(namxx, function(x) paste0(x, ".", seq_along(fpp[[x]])),
        simplify = F, USE.NAMES = F
      ))
    }
  }
  
  xtr.out
}
