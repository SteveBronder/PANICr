#'@title Beta Coefficients for standard OLS
#'
#'@description Returns the Beta values for an Ordinary Least Squares.
#'
#'@usage myols(x,y)
#'
#'@param x A matrix of dependent variables.
#'
#'@param y A vector that is used as the independent variable.
#'
#'@export
myols <- function(x, y) {
    
    bhat <- qr.solve(crossprod(x)) %*% t(x) %*% y
    return(bhat)
} 
