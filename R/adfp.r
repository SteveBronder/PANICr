#'@title Generalized Least Squares Modified Dickey-Fuller t test
#'
#'@description This function performs a modified Dickey-Fuller t test for a
#' unit root in which the series has been modified by a generalized least
#' squares regression.
#'
#'@usage adfp(y, penalty, kmax, kmin, p)
#'
#'@param y A matrix of data
#'
#'@param penalty An integer value of either 0 or 1. 0 uses the MAIC, a penalty on
#' k that accounts for the bias in the sum of the autoregressive coefficient.
#'  1 uses the more general form MIC.
#'
#'@param kmax An integer of the maximum number of lags for the vector autoregressions. An
#' upper bound of (12*(T/100)^.25)^8 is suggested
#' in Schwert (1989)
#'
#'@param kmin An integer of the minimum number of lags for the vector autoregression. k = 0
#' is a reasonable point.
#'
#'@param p An integer with value of either 0 or -1. a value of -1 will modify the series
#' with a generalized least squares regression.
#'
#'@return adf A numeric vector of t tests for the dfgls of each column. Will have to
#' find rejection levels
#'
#'@return kstar A numeric vector of the lags for each column's vector autoregression.
#'
#'@export
adfp <- function(y, penalty, kmax, kmin, p) {
    ####### adfg??
    y <- as.matrix(y)
    
    if (p == -1) {
        
        yt = y
    } else {
        
        yt = glsd(y, p)
    }
    
    reg <- lagn(yt, 1)
    
    dyt <- mydiff(yt, 1)
    
    kstar <- s2ar(yt, penalty, kmax, kmin)
    
    for (i in 1:kstar) {
        
        reg <- cbind(reg, lagn(dyt, i))
    }
    
    dyt <- as.matrix(trimr(dyt, I(kstar + 1), 0))
    
    reg <- as.matrix(trimr(reg, I(kstar + 1), 0))
    
    rho <- myols(reg, dyt)
    
    e <- dyt - reg %*% rho
    
    nef <- dim(dyt)[1]
    
    s2e <- t(e) %*% e/nef
    
    xx <- solve(crossprod(reg))
    
    sre <- xx[1, 1] * s2e
    
    adf <- rho[1, 1]/sqrt(sre)
    
    rho1 <- rho[1, 1] + 1
    
    if (kstar > 0) {
        
        sumb = sum(rho[2:I(kstar + 1), ])
    } else {
        
        sumb = 0
    }
    
    s2vec <- s2e/((1 - sumb)^2)
    
    results <- list(adf = adf, kstar = kstar)
    
    return(results)
} 
