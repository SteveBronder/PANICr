#'@title PANIC (2004) MCMC Non-Stationarity Tests on Common and Idiosyncratic Components
#'
#'@description This function performs an MCMC over PANIC (2010) Model C, PAC, and
#' PMSB tests. PAC estimates the pooled autoregressive coefficient, PMSB uses a sample
#' moment, and Model C performs the MP test while projecting on intercept and trend.
#'  The sample moments test is based off of the modified Sargan-Bhargava test (PMSB).
#'
#'@usage MCMCpanic04(x, nfac, k1, jj,burn = 1000, mcmc = 10000, thin = 10, verbose = 0,
#'seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
#'a0 = 0.001, b0 = 0.001, std.var = TRUE,...)
#'
#'
#'@param x A NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 The maximum lag allowed in the ADF test.
#'
#'@param jj an Integer 1 through 8. Choices 1 through 7 are respectively, IC(1),
#' IC(2), IC(3), AIC(1), BIC(1), AIC(3), and BIC(3), respectively. Choosing 8
#' makes the number of factors equal to the number of columns whose sum of
#' eigenvalues is less than  or equal to .5.
#'
#'@param burn The number of burn in iterators for the sampler
#'
#'@param mcmc The number of iterations in the sampler
#'
#'@param thin The thinning interval used in the simulation. mcmc must be divisible by this value.
#'
#'@param verbose A positive integer which determines whether or not the progress of the
#' sampler is printed to the screen. If verbose is greater than 0 the iteration 
#' number and the factor loadings and uniquenesses are printed to the screen 
#' every verboseth iteration.
#'
#'@param seed The seed for the random number generator.
#'
#'@param lambda.start Starting values for the factor loading matrix Lambda.
#'
#'@param psi.start Starting values for the uniqueness
#'
#'@param l0 The means of the independent Normal prior on the factor loadings
#'
#'@param L0 A scalar or a matrix with the same dimensions as lambda. The precision (inverse variances)
#' of the independent Normal prior on the factor loadings.
#' 
#'@param a0 scalar or a k-vector. Controls the shape of the inverse Gamma prior on the uniqueness.
#'
#'@param b0 Controls the scale of the inverse Gamma prior on the uniqueness.
#'
#'@param std.var if TRUE the variables are rescaled to have zero mean and unit variance.
#' Otherwise, the variables are rescaled to have zero mean, but retain their observed variances
#' 
#'@param ... extra parameters to pass to MCMCfactanal
#' 
#'@return adf.mcmc A list of the MCMC samples of the test statistics. Returns the test statistics
#'  Pooled Cointegration a,  Pooled Cointegration b, Pooled Idiosyncratic a, 
#'  Pooled Idiosyncratic b, Pooled Demeaned test, and tests on Common components. The critical values for the Pooled Cointegration test
#'  can be found on this packages vignette or in Bai and Ng (2004). The pooled idiosyncratic test has a critical
#'  value of 1.64. The Pooled Demeaned test has a critical value of 2.87. The common components have a critical
#'  value of -2.86.
#'
#'@return factor_mcmc The MCMC results from MCMCfactanal()
#' 
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#' 
#' @references ndrew D. Martin, Kevin M. Quinn, Jong Hee Park (2011). MCMCpack: Markov Chain Monte Carlo
#' in R. Journal of Statistical Software. 42(9): 1-21. URL http://www.jstatsoft.org/v42/i09/.
#'
#'@export
MCMCpanic04 <- function(x, nfac, k1, jj, burn = 1000, mcmc = 10000, thin = 10, verbose = 0, seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, a0 = 0.001, 
                        b0 = 0.001, std.var = TRUE,...) {
  
  ###  
  # Begin TESTS
  ###
  is.xts(x) || stop("x must be an xts object so lags and differences are taken properly")
  
  
  Tn <- dim(x)[1]
  N <- dim(x)[2]
  
  nfac < N || stop(" nfac must be less than the number of series.")
  if (is.null(nfac)){
    warning("nfac is NULL, setting the maximum number of factors equal to the number of columns")
    nfac <- dim(x)[2]
  }
  
  if (is.null(k1)){
    warning("k1 is NULL, setting k1 equal to  k1 4 * ceiling((T/100)^(1/4))")
    k1 <- 4 * ceiling((dim(x)[1]/100)^(1/4))
  }
  if (!(k1 %% 1 == 0)){
    stop(" k1 must be an integer.")
  }
  if (is.null(criteria)){
    warning("criteria is NULL, setting criteria to BIC3")
    criteria <- "BIC3"
  }
  ####
  ## End Tests
  ####
  
  # center, trim, and difference x
  x_dm <- scale(x,center = TRUE,scale = FALSE)
  x_trim <- x_dm[2:Tn, ]
  x_diff <- diff(x, 1)[2:Tn,]
  
  # Value for scaling?
  scaler <- sqrt(N) * Tn
  
  # approximate factor model
  factors <- getnfac(x_diff, nfac, criteria)
  
  
  ic <- factors$ic
  
  if (is.null(ic)) {
    ic <- 1
  }
  
  PC <- pc(x_diff, ic)
  
  if (ic == 1) {
    p = 0
  } else if (ic == 0) {
    p = -1
  } else if (ic > 1) {
    p = 1
  }
  
  lamhat <- PC$lambda
  
  dfhat <- PC$fhat
  
  fac.test <- MCMCfactanal(~.,
                           factors = ic,
                           data = as.data.frame(x_diff),
                           burnin = burn,
                           mcmc = mcmc,
                           thin = thin,
                           verbose = verbose,
                           seed = seed,
                           lambda.start = lambda.start, 
                           psi.start = psi.start,
                           l0 = l0,
                           L0 = L0,
                           a0 = a0,
                           b0 = b0,
                           store.scores = TRUE,
                           std.var = std.var)
  
  
  
  for (j in 1:I(mcmc/thin)) {
    
    lamhat <- matrix(fac.test[j, 1:I(N * ic)], N, ic)
    
    dfhat <- matrix(fac.test[j, I(N * ic + N + 1):I((Tn - 1) * ic + N * ic + N)], I(Tn - 1), ic, byrow = TRUE)
    
    dehat <- x_diff - tcrossprod(PC$fhat, PC$lambda)
    
    ehat0 <- apply(x_diff - tcrossprod(dfhat, lamhat), 2, cumsum)
    
    fhat0 <- apply(dfhat, 2, cumsum)
    
    reg <- cbind(matrix(1, I(Tn - 1), 1), fhat0[[i]])
    
    ehat1 <- matrix(0, I(Tn - 1), N)
    
    beta1 <- matrix(0, I(ic + 1), N)
    
    for (i in 1:N) {
      
      beta1[, i] <- qr.solve(reg, x_trim[, i])
      
      ehat1[, i] <- x_trim[, i] - reg %*% beta1[, i]
    }
    
    
    adf20 <- matrix(0, 1, ic)
    
    for (i in 1:ic) {
      
      adf20[, i] <- adf04(fhat0[, i], k1, p)
    }
    
    
    adf30 <- adf04(ehat0, k1, -1)  # test ehat0
    
    adf50 <- adf04(ehat1, k1, -1)  # test ehat1
    
    # now do the pooled test
    
    
    padf30 <- pool(adfnc, adf30)
    
    adf30a <- padf30$adf31a
    
    adf30b <- padf30$adf31b
    
    padf50 <- poolcoint(coint0, adf50, ic)
    
    adf50a <- padf50$pvala
    
    adf50b <- padf50$pvalb
    
  }
  
  adf20ab <- as.data.frame(matrix(unlist(adf20), mcmc, ic, byrow = TRUE))
  
  for (i in 1:ic) {
    colnames(adf20ab)[i] <- paste0("Common", i)
  }
  
  adf.tests <- cbind(adf50a, adf50b, adf30a, adf30b, adf20ab)
  
  
  colnames(c("Pooled Cointegration a", "Pooled Cointegration b", "Pooled Idiosyncratic a", "Pooled Idiosyncratic b", "Pooled Demeaned"))
  results <- list(adf.mcmc = adf.tests, factor_MCMC = fac.test)
} 