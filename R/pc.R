#'@title Principle Component Analysis for PANIC (2004)
#'
#'@description This function performs the principle component analysis in order to determine
#' the Common and Idiosyncratic Components of the factor model.
#'
#'@usage pc(y,nfac)
#'
#'
#'@param y An NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@return ehat A matrix with the Idiosyncratic component of the factor model
#'
#'@return fhat A matrix with the factors of the approximate factor model
#'
#'@return lambda A matrix with the factor loadings of the approximate factor model
pc <- function(y, nfac) {
  
  bigt <- dim(y)[1]
  
  bign <- dim(y)[2]
  
  eig <- svd(crossprod(y))
  
  Fhat0 <- eig$u
  
  eigval <- as.matrix(eig$d)
  
  Fhat1 <- eig$v
  
  lambda <- Fhat0[, 1:nfac] * sqrt(bign)
  
  fhat <- y %*% (lambda/bign)
  
  ehat <- y - tcrossprod(fhat, lambda)
  
  return(list(ehat = ehat, fhat = fhat, lambda = lambda))
} 
