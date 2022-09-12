# a function that transforms each column of matrix x using final fp powers
# x  = matrix with all continuous variables shifted and scaled and its column
# names are equal to the names of power.list.
# power.list =  is a named list of fp powers e.g list(x1 = c(0,1), x2 = 0, x3 = 1)
# if all elements of power.list are NA then we return NULL otherwise we return
# a matrix
# center = a vector of TRUE or FALSE specifying whether x should be centered
xtransform <- function(x, power.list, center, acdx){
  # check whether all power.list are equal to NA-all variables were eliminated
  if(all(is.na(unlist(power.list)))){
    # return NULL
    xtr.out <- NULL
  }else{
      namx <- names(power.list)
      x <- x[,namx, drop = F]
      # subset scale, center and shift using names of powers selected-
      # scale <- scale[namx] # we assume x has been shifted and scaled
      center <- center[namx]
      acdx <- acdx[namx]
      # transform acd variables if any
      if(any(acdx)){
        # select FP powers for acd variables
        pow.acd <- power.list[acdx]
        # check whether the acd variables were selected
        if(all(is.na(unlist(pow.acd)))){
          xtr1 <- NULL
        }else{
        xnames1 <- names(pow.acd)
        x.acdx <- x[,xnames1, drop = F]
        xtrans.acd <- vector(mode = "list", length = length(pow.acd))
        for(i in seq_along(xnames1)){
          xtrans.acd[[i]] <- afracgen(x = x.acdx[,i,drop = T],
                                         power = pow.acd[[i]], scale = 1,#scale = scale[i],
                                     center = center[xnames1][i],shift = NULL) #shift = shift[i])
        }
        # cbind acd transformed variables
        xtr1 <- do.call(cbind, xtrans.acd)
        # rename acd transformed variables. they start with letter "A"
        snames = unlist(sapply(xnames1, function(x) paste0(c(x,paste0("(A)",x)), ".", c(1,1)),
                      simplify = F, USE.NAMES = F))
        colnames(xtr1) <- snames[which(!is.na(unlist(pow.acd)))]
        }
        # Variables without acd transformations
        pow <- power.list[!acdx]
        if(all(is.na(unlist(pow)))){
          xtr2 <- NULL
        }else{
          # get rid of unselected variables
          pow2 <-pow[!is.na(pow)] # fracgen returns NULL when power = NA
        xnames2 <- names(pow2)
        xx <- x[,xnames2,drop = F]
        xtrans <- vector(mode = "list", length = length(pow2))
        for(i in seq_along(xnames2)){
          # x is already scaled and shifted so fracgen shifting is irrelevant
          xtrans[[i]] <- fracgen(x = xx[,i,drop = T], power = pow2[[i]], scale = 1,#scale = scale[i],
                                 center = center[xnames2][i],shift = NULL) #shift = shift[i])
        }
        # cbind non-acd transformed variables-
        xtr2 <- do.call(cbind, xtrans)
        # rename xtr such that if x1 is fp2 then we have x1.1, x1.2 whereas fp1 is x1.1
        # we can as well use make.names(rep(names(fpp), lapply(fpp, length)), sep = ".")
        # but has undesirable names like x1, x1.1 instead of x1.1, x1.2
        colnames(xtr2) <- unlist(sapply(xnames2, function(x) paste0(x, ".", seq_along(pow2[[x]])),
                                       simplify = F, USE.NAMES = F))
        }
        # combine xtr1 and xtr2
        xtr.out <- cbind(xtr1, xtr2)
        # Usual MFP without acd variables
      }else{
        # # Get rid of unselected variables denoted by NA in power.list
        fpp <-power.list[!is.na(power.list)]
        namxx <- names(fpp)
        # subset x
        x <- x[,namxx,drop = F]
        xtransx <- vector(mode = "list", length = length(fpp))
        for(i in seq_along(namxx)){
          # x is already scaled and shifted so fracgen shifting is irrelevant
          xtransx[[i]] <- fracgen(x = x[,i,drop = T], power = fpp[[i]], scale = 1,#scale = scale[i],
                                 center = center[namxx][i],shift = 0) #shift = shift[i])
        }
        xtr.out <- do.call(cbind, xtransx)
        colnames(xtr.out) <- unlist(sapply(namxx, function(x) paste0(x, ".", seq_along(fpp[[x]])),
                                        simplify = F, USE.NAMES = F))
      }

      }
  return(xtr.out)
}
# acdx = rep(F,ncol(x))
# acdx = setNames(acdx, colnames(x))
# acdx = replace(acdx, c(3,7), rep(T,2))
# pds = list(age = NA, bph = 1, cavol = c(1,1),cp = NA,pgg45 = NA,svi = NA,weight = c(1,1))
# xtransform(x = x, power.list = PWWW, center = center, acdx = acdx)
# xtransform(x = x, power.list = pds, center = center, acdx = acdx)
#
#
# ps = pw[c("cavol","weight")]
# ps = list(cavol = c(1,NA), weight = c(NA,1))
# unlist(sapply(namx, function(x) paste0(x, ".", seq_along(pw[[x]])),
#               simplify = F, USE.NAMES = F))
# unlist(sapply(names(ps), function(x) paste0(c(x,paste0(x,"a")), ".", c(1,1)),
#               simplify = F, USE.NAMES = F))
#
#
#
# fracgen.acd <- function(x, power = c(1,1),shift = NULL, s = NULL, scale = NULL,center = F){
#   # the length of powers must be equal to 2
#   np <- length(power)
#   if(np!=2) stop("The length of powers supplied is ",np, ". Two powers are needed")
#   # Approximate cumulative distribution (ACD)
#   xa <- acd(x, power = NULL, shift = NULL, s = NULL, scale = NULL)$acd
#   # Transform x and acdx using the supplied powers
#   x<- fracgen(x = x, power = power[1], scale = scale, shift = shift, center = F)
#   Xa <- fracgen(x = xa, power = power[2], scale = scale, shift = shift, center = F)
#   # return a matrix of two variables
#   xx <-cbind(x, xa)
#   colnames(xx) <- c("x","xa")
#   return(xx)
# }
#
