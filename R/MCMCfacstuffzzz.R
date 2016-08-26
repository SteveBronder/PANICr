#'@title PANIC (2004) MCMC Non-Stationarity Tests on Common and Idiosyncratic Components
#'
#'@description This function performs an MCMC over the 
#' tests on the idiosyncratic and common component from PANIC (2004).
#'
#'@usage MCMCpanic04(x, nfac, k1, criteria = NULL,burn = 1000,
#' mcmc = 10000, thin = 10, verbose = 0, seed = NA,
#' lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
#' a0 = 0.001, b0 = 0.001, std.var = TRUE,...)
#'
#'
#'@param x An object of class xts with each column being a time series
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 an integer that is the maximum lag allowed in the ADF test.
#'
#'@param criteria a character vector of length one with a value of either
#' IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen. Choosing eigen makes
#' the number of factors equal to the number of columns whose sum of
#' eigenvalues is less than  or equal to .5.
#'
#'@param burn Integer of the number of burn in iterators for the sampler
#'
#'@param mcmc Integer of the number of iterations in the sampler
#'
#'@param thin Integer of the thinning interval used in the simulation. mcmc must be divisible by this value.
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
#'@return mcmc_tests An mcmc object containing the resampled tests on the common components as well as
#' the test on the idiosyncratic component.
#' 
#'@return factor_mcmc The results from MCMCfactanal()
#' 
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#' 
#' @references Andrew D. Martin, Kevin M. Quinn, Jong Hee Park (2011). MCMCpack: Markov Chain Monte Carlo
#' in R. Journal of Statistical Software. 42(9): 1-21. URL http://www.jstatsoft.org/v42/i09/.
#'
#'@export
#'
#'@import coda
#'@import MCMCpack
#'@import xts
MCMCpanic04 <- function(x, nfac, k1, criteria = NULL, burn = 1000, mcmc = 10000, thin = 10, verbose = 0, seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, a0 = 0.001, 
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
  
  
  
  mcmc_stuff <- lapply(1:I(mcmc/thin), function(j){
    
    lamhat <- matrix(fac.test[j, 1:I(N * ic)], N, ic)
    
    dfhat <- matrix(fac.test[j, I(N * ic + N + 1):I((Tn - 1) * ic + N * ic + N)], I(Tn - 1), ic, byrow = TRUE)
    fhat0 <- xts::reclass(dfhat, match.to = x_diff)
    fhat0 <- cumsum(fhat0)
    
    dehat <- x_diff - tcrossprod(dfhat, lamhat)
    ehat0 <- cumsum(dehat)
    
    # Set up matrices for 
    # 1. Idiosyncratic components
    # 2. Regression of PC on x
    # 3. coefficients of regression
    reg <- cbind(matrix(1, I(Tn - 1), 1), fhat0)

    # Do regressions for each idiosync component
    
    beta1 <- qr.solve(reg,x_trim)
    ehat1 <- x_trim - reg %*% beta1
    
    adf_fhat0 <- matrix(0, factors$ic, 1)
    adf_ehat0 <- matrix(0, N, 1)
    adf40 <- matrix(0, N, 1)
    adf_ehat1 <- matrix(0, N, 1)

    # test fhat0 for a unit root
    for (i in 1:factors$ic) {
      
      adf_fhat0[i, ] <- adf04(fhat0[, i], k1, p)  
    }
    
    # test ehat0 and ehat1 for a unit root
    adf_ehat0 <- adf04(ehat0, k1, -1)  
    adf_ehat1 <- adf04(ehat1, k1, -1)  
    
    padf_ehat0 <- pool(adfnc, adf_ehat0)
    adf_ehat0a <- padf_ehat0$adf31a
    adf_ehat0b <- padf_ehat0$adf31b
    
    # MQ test
    # coint0 is the Phillips-Ouliaris Zt test with a constant.
    padf_ehat1 <- poolcoint(coint0, adf_ehat1, factors$ic)
    adf_ehat1a <- padf_ehat1$pvala
    adf_ehat1b <- padf_ehat1$pvalb
    
    output <- as.data.frame(matrix(c(adf_fhat0,adf_ehat0b),1,(length(adf_fhat0) + 1)))
    colnames(output) <- c(paste0("Common",1:I(length(adf_fhat0))),"idiosyncratic")
    return(output)
  })
  
  mcmc_stuff <- do.call(rbind,mcmc_stuff)
  
  mcmc_stuff <- coda::as.mcmc(mcmc_stuff)
  attributes(mcmc_stuff)$mcpar[3] <- thin
 
  results <- list(mcmc_tests = mcmc_stuff, factor_MCMC = fac.test)
} 