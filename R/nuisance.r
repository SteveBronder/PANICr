#'@title Estimate the nuisance parameters of the error term
#'
#'@description This function estimates the short-run, long-run, and one sided
#' variance of the error term
#'
#'@usage nuisance(res,k)
#'
#'@param res A matrix consisting of the residuals from a factor model.
#'
#'@param k If fixk is 0, then automatic bandwidth selection is performed.
#' Otherwise, the integer placed here will be the selected bandwidth.
#'
#'@return A data frame containing the following columns
#'
#'@return sig2 The vector of short run variances
#'
#'@return omega2 The vector of long run variances
#'
#'@return half The vector of one-sided variances
#'
#'@export

nuisance <- function(res, k) {
  
  # convert to data frame so that vectors are treated as vectors instead of lists of numbers
  # TODO: Have an if test, if vector, do not use lapply
  nw_output<- lapply(as.data.frame(res),function(x){
    
    NW <- nw(x, k)
    
    Tn <- length(x)
    K <- NW$k
    
    Kw <- NW$w
    
    sig2 = crossprod(x)/Tn
    
    omega2 <- sig2
    
    half = 0
    
    for (j in 1:K) {
      temp <- t(x[1:I(Tn - j - 1)]) %*% x[I(j + 1):I(Tn - 1)]/Tn
      
      omega2 <- omega2 + 2 * Kw[j] * temp
      
      half <- half + Kw[j] * temp
    }
    c(sig2,omega2,half)
  })
  
  output <- data.frame(do.call(rbind,nw_output))
  
  colnames(output) <- c("sig2","omega2","half")
  
  return(output)
} 