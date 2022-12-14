# calculates the final degrees of freedom for the selected model
# x is the power(s) so if x = c(1,2) then df = 4 but if x = NA then df = 0
calculate_df <- function(x) {
  if (all(is.na(x))) {
    # df for unselected variable
    df <- 0
  } else {
    # Remove NAs in powers (2,NA) or (NA,1) etc. Happens because of acd
    x <- x[!is.na(x)]
    # df of linear function
    if (identical(x, 1)) {
      df <- 1
      # df of fpm. Note that df = 2m where m is the degree. if length = 1 then
      # degree = 1 and df = 2 etc
    } else {
      df <- 2 * length(x)
    }
  }
  return(df)
}
