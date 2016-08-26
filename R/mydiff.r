#'@title Difference a matrix
#'
#'@description This function differences an entire matrix by some defined
#' constant k.
#'
#'@usage mydiff(x,k)
#'
#'@param x The matrix to be differenced.
#'
#'@param k An integer of the length of difference to perform on each vector of the NxT matrix.
#'
#'@return xx A matrix that has been differenced k times.
#'
#'@export

mydiff <- function(x, k) {
    
    nt <- dim(x)[1]
    
    nc <- dim(x)[2]
    
    if (k == 0) {
        
        xx <- x
    } else {
        
        x1 <- as.matrix(trimr(x, k, 0))
        
        x2 <- as.matrix(trimr(lagn(x, k), k, 0))
        
        zero <- matrix(0, k, nc)
        
        xx <- rbind(zero, I(x1 - x2))
    }
} 
