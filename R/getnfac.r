#'@title Determining The Number of Factors In Approximate Factor Model
#'
#'@description This function approximates the number of factors in an approximate factor model
#'  for large N by T matrices using the methods and criteria found in Bai and Ng (2002)
#'
#'@usage getnfac(x,kmax,criteria)
#'
#'@param x A matrix containing the data.
#'
#'@param kmax An integer with the maximum number of common factors to search over. This methodology
#' is weak to underestimation of the number of common factors so setting this value higher is preferred.
#'
#' @param criteria a character vector of length one with values of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'
#' @details This function approximates the number of factors in an approximate
#' factor model. Amongst the penalty functions BIC(3) has been found to be
#' strict against cross-sectional dependence and is recommended for panels with greater than 18 series.
#'  IC(1) is most commonly used. BIC(1) is not recommended for small N relative to T. AIC(3) and BIC(3) take into
#' account the panel structure of the data. AIC(3) performs consistently
#' across configurations of the data while BIC(3) performs better on
#' large N data sets.
#'
#' @return ic Integer of the approximate number of factors based off of the chosen
#' penalty function
#'
#' @return lambda A matrix of the estimated factor loadings associated with common factors.
#'
#' @return Fhat A matrix of the estimated common components
#'
#'
#' @references Jushan Bai and Serena Ng. 'Determining the Number of Factors in
#' Approximate Factor Models.' Econometrica 70.1 (2002): 191-221. Print.
#'
#'@export
getnfac <- function(x, kmax = NULL, criteria = NULL) {
  
  
  #######
  ## Begin: Tests
  #######
  all(sapply(x, is.numeric) == TRUE)  || stop("All columns must be numeric")
  is.xts(x) || stop("x must be an xts object so lags and differences are taken properly")
  if (!(kmax %% 1 == 0)){
    stop(" k1 must be an integer.")
  }
  if (is.null(criteria)){
    warning("criteria is NULL, setting criteria to BIC3")
    criteria <- "BIC3"
  }
  ########
  
    Tn <- dim(x)[1]
    
    N <- dim(x)[2]
    
    NT <- N * Tn
    
    NT1 <- N + Tn
    
    CT <- matrix(0, 1, kmax)
    
    ii <- seq(1:kmax)
    
    GCT <- min(N, Tn)
    
    
    if(is.null(criteria)){
      warning("Criteria is NULL, setting to BIC3")
      criteria <- "BIC3"
    }
    if (criteria == "IC1") {
        CT[1, ] <- log(NT/NT1) * ii * NT1/NT
    }
    
    if (criteria == "IC2") {
        CT[1, ] <- (NT1/NT) * log(GCT) * ii
    }
    
    if (criteria == "IC3") {
        CT[1, ] <- ii * log(GCT)/GCT
    }
    
    if (criteria == "AIC1") {
        CT[1, ] <- 2 * ii/Tn
    }
    
    if (criteria == "BIC1") {
        CT[1, ] <- log(T) * ii/Tn
    }
    
    if (criteria == "AIC3") {
        CT[1, ] <- 2 * ii * NT1/NT
    }
    
    if (criteria == "BIC3") {
        CT[1, ] <- log(NT) * ii * NT1/NT
    }
    
      
    
    IC1 <- matrix(0, dim(CT)[1], I(kmax + 1))
    
    Sigma <- matrix(0, 1, I(kmax + 1))
    
    XX <- tcrossprod(x)
    
    eig <- svd(t(XX))
    
    Fhat0 <- eig$u
    
    eigval <- as.matrix(eig$d)
    
    Fhat1 <- eig$v
    
    sumeigval <- apply(eigval, 2, cumsum)/sum(eigval)
    
    if (criteria != "eigen") {
        for (i in kmax:1) {
            
            Fhat <- Fhat0[, 1:i]
            
            lambda <- crossprod(Fhat, x)
            
            chat <- Fhat %*% lambda
            
            ehat = x - chat
            
            Sigma[i] <- mean(sum(ehat * ehat/Tn))
            
            IC1[, i] <- log(Sigma[i]) + CT[, i]
        }
        
        Sigma[kmax + 1] <- mean(sum(x * x/Tn))
        
        IC1[, kmax + 1] <- log(Sigma[kmax + 1])
        
        ic1 <- minindc(t(IC1))
        
        ic1 <- ifelse(ic1 <= kmax, ic1 * 1, ic1 * 0)
        
    }
    
    if (criteria == "eigen") {
        
        for (j in 1:I(nrow(sumeigval))) {
            
            if (sumeigval[j] >= 0.5) {
                ic1 <- j
                break
            }
        }
        
    }
    
    if (ic1 == 0) {
        
        Fhat = matrix(0, T, N)
        
        lambda = matrix(0, N, T)
        
        chat = matrix(0, T, N)
    } else {
        
        Fhat <- Fhat0[, 1:ic1]
        
        lambda <- crossprod(x, Fhat)
        
        chat <- Fhat %*% t(lambda)
    }
    
    
    output <- list(ic = ic1, lambda = lambda, Fhat = Fhat)
    
    return(output)
} 
