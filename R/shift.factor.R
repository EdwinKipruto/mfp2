# estimates factors for shifting continuous variables to positive values. If x
# is a two level variable or all values are positive then a factor of 0 is return
# hence no shifting will be conducted
find_shift_factor <- function(x) {
  nxx <- length(unique(x))
  if (all(x > 0) || nxx <= 2) {
    adj.factors <- 0
  } else {
    # see section 4.7 in Royston and Sauerbrei (2008) book for the formula
    difx <- diff(sort(x))
    # get rid of 0 differences and find the minimum
    eps <- min(difx[difx != 0])
    # adjustment factor
    adj.factors <- eps - min(x, na.rm = T)
  }
  return(adj.factors)
}
