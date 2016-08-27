#'@title Pooling Function for Cointegration test PANIC (2004)
#'
#'@description This function find the P values for the pooled cointegration test in PANIC (2010)
#'
#'@usage poolcoint(a,x,r)
#'
#'@param a A matrix containing the p values
#'
#'@param x A matrix containing the adf test to be pooled
#'
#'@param r An integer for the number of factors determined by getnfac()
#'
#'@return pvala a numeric vector of the fisher sum of the p-values for the cointegration test
#'
#'@return pvalb a numeric vector containing the critical value of the cointegration test
#'
#'@export
poolcoint <- function(a, x, r) {
    
    x <- as.matrix(x)
    
    N <- ncol(x)
    
    pval <- matrix(0, N, 1)
    
    r <- ifelse(r > 4, 4, r)
    
    for (i in 1:N) {
        
        aa <- abs(a[, r] - x[, i])
        
        j1 <- min(aa)
        
        j2 <- which.min(aa)
        
        pval[i] <- a[j2, 4]
    }
    
    pvala <- -2 * sum(log(pval))
    
    pvalb <- (pvala - 2 * N)/sqrt(4 * N)
    
    results <- list(pvala = pvala, pvalb = pvalb)
    
    return(results)
    
} 
