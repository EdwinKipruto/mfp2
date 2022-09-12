# converts a nested list with same or different length into a matrix
powermat <- function(power.list){
  # Check the maximum number of powers i.e  FP2 has 2 while FP1 has 1
  psize <- sapply(power.list, length)
  maxp <- max(psize)
  # Create a new nested list of same length. This means that if FP1 was choosen
  # for x then the second power should be NA
  new.list.powers <- vector(mode = "list", length = length(power.list))
  for(i in 1:maxp){
    new.list.powers[[i]]<- sapply(power.list,function(x) x[i])
  }
  # combine the powers and rename.
  matp <- do.call(cbind, new.list.powers)
  colnames(matp) <- paste0("power", 1:maxp)
  return(matp)
}
