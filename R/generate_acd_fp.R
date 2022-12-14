# similar to generate_fp but for acd transformation.
generate_acd_fp <- function(x, 
                                               power = c(1, 1),
                                               shift = NULL, 
                                               s = NULL, 
                                               scale = NULL, 
                                               center = F) {
  # the length of powers must be equal to 2
  np <- length(power)
  if (np != 2) stop("The length of powers supplied is ", np, ". Two powers are needed")
  if (all(is.na(power))) {
    xx <- NULL
  } else {
    # Transform x and acdx using the supplied powers
    # Approximate cumulative distribution (ACD)
    xa <- acd(x, power = NULL, shift = shift, s = s, scale = scale)$acd
    xa <- generate_fp(x = xa, power = power[2],
                                         scale = scale, shift = shift, 
                                         center = center)
    x2 <- generate_fp(x = x, power = power[1],
                                         scale = scale, shift = shift,
                                         center = center)
    # return a matrix of two variables
    xx <- cbind(x2, xa)
  }

  return(xx)
}
