#This script generates the samples for global sensitivity sobol2007 function
genSamples <- function(n, l){  #n=number of samples and l=named list with parameter value ranges
  df <- as.data.frame(l)
  len <- length(df)
  X <- data.frame(matrix(ncol=len, nrow=n))
  names(X) <- names(df)
  for(i in 1:len){
    r <- runif(n, df[1,i], df[2,i])
    X[,i] <- r
  }
  return(X)
}