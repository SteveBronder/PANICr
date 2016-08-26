#'@title General Least Squares Detrending
#'
#'@description This function detrends the data by general least squares.
#'
#'@usage glsd(y,p)
#'
#'@param y A matrix containing the data to be detrended
#'
#'@param p An integer with value of either 0 or 1 which decides the value of the penalty term, chat.
#' This is either either -7 or -13.5, respectively
#'
#'@return yt A matrix of the detrended data
#'
#'@export
glsd <- function(y, p) {
    
    y <- as.matrix(y)
    
    if (p == 0) {
        cbar = -7
    }
    
    if (p == 1) {
        
        cbar = -13.5
    }
    
    nt <- nrow(as.matrix(y))
    
    z <- matrix(1, nt, 1)
    
    if (p == 1) {
        
        z <- cbind(z, t(seq(1:nt)))
    }
    
    abar <- 1 + cbar/nt
    
    ya <- matrix(0, nt, 1)
    
    za <- matrix(0, nt, ncol(z))
    
    ya[1, ] <- y[1]
    
    za[1, ] <- z[1, ]
    
    ya[2:nt, 1] <- y[2:nt, 1] - abar * y[1:I(nt - 1), 1]
    
    za[2:nt, ] <- z[2:nt, ] - abar * z[1:I(nt - 1), ]
    
    bhat <- qr.solve(za, ya)
    
    yt <- y - z %*% bhat
    
    
    return(yt)
    
    
} 
