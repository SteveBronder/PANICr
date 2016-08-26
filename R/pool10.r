#'@title Pooling Function for PANIC (2010)
#'
#'@description This function find the P values for the pooled test in PANIC (2010)
#'
#'@usage pool(p_values,test_values)
#'
#'
#'@param p_values a data frame containing the p values
#'
#'@param test_values a matrix of adf test values to be pooled
#'
#'@return adf31a a numeric vector of the fisher sum of the p-values
#'
#'@return adf31b a numeric vector containing the critical value of the pooling test
#'
#'@export

pool <- function(p_values, test_values) {
    
    p_values <- as.matrix(p_values)
    
    test_values <- as.matrix(test_values)
    
    N <- ncol(test_values)
    
    pval <- matrix(0, N, 1)
    
    for (i in 1:N) {
        
        aa <- abs(p_values[, 1] - test_values[, i])
        
        j1 <- min(aa)
        
        j2 <- which.min(aa)
        
        pval[i, ] <- p_values[j2, 2]
        
    }
    
    pvala <- -2 * sum(log(pval))
    
    pvalb <- (pvala - 2 * N)/sqrt(4 * N)
    
    output <- list(adf31a = pvala, adf31b = pvalb)
    
    return(output)
} 
