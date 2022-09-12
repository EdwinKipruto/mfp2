# A function that determines the number of unique values in a variable. Variables
# with fewer than or equal to three unique values, for example, will be assigned
# df = 1. df = 2 will be assigned to variables with 4-5 unique values, and
# df = 4 will be assigned to variables with unique values greater than or equal to 6.
# x = matrix of predictors
# df.default = default df which is 4
df.assign <- function(x, df.default){
  # Set default degrees of freedom for each variable
  df <- rep(df.default, ncol(x))
  # calculate the length of unique values of each column of x
  nu <- apply(x, 2,function(x) length(unique(x)))
  # Find the indices of variables with <=3 uniques values and set their df = 1
  # again find the variables with uniques values between 4 and 5 inclusive and set their df = 2
  index1 <- which(nu<=3)
  index2 <- which(nu>=4&nu<=5)
  if(length(index1)!=0){df <- replace(df, index1, rep(1,length(index1)))}
  if(length(index2)!=0){df <- replace(df, index1, rep(1,length(index1)))}
  return(df)
}
