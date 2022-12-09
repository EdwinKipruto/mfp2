# Function that reset acdx of variable with less than 5 unique values to false
reset.acd <- function(x, acdx) {
  # calculate the length of unique values of each column of x
  nu.acdx <- apply(x, 2, function(x) length(unique(x)))
  nu.acdx <- nu.acdx[names(acdx)]
  # check variables with unique length<5
  index <- which(nu.acdx <= 4)
  if (length(index) == 0) {
    acdx <- acdx
  } else { # reset acdx of variables with less than 5 unique values to FALSE
    # check whether the acdx for variables less than 5 unique values are all false
    dd <- which(acdx[index])
    if (length(dd) == 0) { # meaning that all acdx for this specific variables are all false. no need of replacement
      acdx <- acdx
    } else {
      warning(paste0(c("Variable", names(acdx[index])[dd], " has fewer than 5 unique values, cannot perform acd transformation"), collapse = " "), call. = F)
      acdx <- replace(acdx, index, rep(F, length(index)))
    }
  }
  return(acdx)
}
