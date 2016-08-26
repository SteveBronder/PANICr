#'@title Bandwidth Selection
#'
#'@description This function at fixk 0 finds the bandwidth, or lag length, of the
#' short, long, and run variance of a residual.
#'
#'@usage nw(v,fixk)
#'
#'@param v A cector of error terms from a factor model
#'
#'@param fixk If fixk is 0, then this function will perform automatic bandwidth
#' selection. Otherwise, the integer placed here will be the selected bandwidth.
#'
#'@return k An integer that is the bandwidth chosen by \code{nw()}
#'
#'@return w The vector of long run variance
#'
#'@references Moon, R. & B. Perron (2004) Testing for a unit root in panels with
#' dynamic factors. Journal of Econometrics 122, 81-126.
#' 
#'@export
#'

nw <- function(v, fixk) {
  
  Tn <- length(v)
  
  rho <- NULL
  
  sigma <- NULL
  
  if (fixk == 0) {
    # auto bandwidth selection
    
    bot <- 0
    
    top <- 0
    
 
    rho <- qr.solve(v[1:I(Tn - 1)], v[2:Tn])
      
    e <- v[2:Tn] - rho[1] * v[1:I(Tn - 1)]
      
    sigma <- crossprod(e)/(Tn - 1)
      
    top <- top + 4 * (rho^2) * (sigma^2)/(((1 - rho)^6) * (1 + rho)^2)
      
    bot <- bot + (sigma^2)/((1 - rho)^4)
    
    
    alpha <- top/bot
    
    k <- ceiling(1.1447 * (alpha * Tn)^(1/3))
    # Trying Something
    if (k > I(Tn)) {
      k = Tn-2
    }
  } else {
    
    k <- fixk
  }
  
  w <- matrix(0, k, 1)
  
  for (i in 1:k) {
    
    x <- i/k
    
    w[i] <- 1 - i/(k + 1)
  }
  
  
  output <- list(k = k, w = w)
  return(output)
} 
