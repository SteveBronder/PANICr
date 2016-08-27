#'@title MCMC PANIC (2010) Sample Moment and PAC tests for Idiosyncratic Component
#'
#'@description This function performs the tests of PANIC (2010) with a Monte Carlo
#' Markov chain based on a Gibbs sampler. One test estimates the pooled autoregressive
#'  coefficient, and one uses a sample moment. The sample moments test is based off of the modified
#' Sargan-Bhargava test (PMSB) while the pooled autoregressive component is based on the
#' Moon and Perron test as well a biased corrected pooled coefficient from PANIC (2004).
#'
#'@usage MCMCpanic10(x = NULL, nfac = NULL, k1 = NULL, criteria = NULL,
#' demean = FALSE, burn = 100, mcmc = 100, thin = 10,
#' verbose = 0, seed = NA, lambda.start = NA, psi.start = NA,
#'  l0 = 0, L0 = 0, a0 = 0.001, b0 = 0.001, std.var = TRUE,...)
#'
#'
#'@param x An object of class xts holding the time series data
#'
#'@param nfac An integer speciyfing the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 an Integer that is the maximum lag allowed in the ADF test.
#' 
#'@param criteria a character vector with values of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'
#'@param demean logical argument. If TRUE, function performs tests on demeaned
#' data. If FALSE, uses non-demeanded data generating process.
#'
#'@param burn The number of burn in iterators for the sampler
#'
#'@param mcmc The number of iterations in the sampler
#'
#'@param thin The thinning interval used in the simulation. mcmc must be divisible by this value.
#'
#'@param verbose A positive integer which determines whether or not the progress of the
#' sampler is printed to the screen. If verbose is greater than 0 the iteration 
#' number and the factor loadings and uniqueness are printed to the screen 
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
#'@param ... extra parameters to be passed to MCMCfactanal
#'
#'@return mcmc_tests An mcmc object containing the resamples of the test statistics.
#' When demeaned, the results will be for model P, PMSB, Model C, and rho1. 
#' When not demeaned, the results will be for model A, model B, PMSB, rho1, and
#' the pooled values on the idiosyncratic component of PANIC (2004).
#'
#'@references Bai, Jushan, and Serena Ng.
#''Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation.'
#' Econometric Theory 26.04 (2010): 1088-1114. Print.
#' 
#'@references Andrew D. Martin, Kevin M. Quinn, Jong Hee Park (2011). MCMCpack: Markov Chain Monte Carlo
#' in R. Journal of Statistical Software. 42(9): 1-21. URL http://www.jstatsoft.org/v42/i09/.
#'
#'@export

MCMCpanic10 <- function(x = NULL,
                        nfac = NULL,
                        k1 = NULL,
                        criteria = NULL,
                        demean = FALSE,
                        burn = 100,
                        mcmc = 100,
                        thin = 10,
                        verbose = 0,
                        seed = NA,
                        lambda.start = NA,
                        psi.start = NA,
                        l0 = 0,
                        L0 = 0, 
                        a0 = 0.001,
                        b0 = 0.001,
                        std.var = TRUE,...) {
  
  #######
  # BEGIN: TESTS
  #######
  all(sapply(x, is.numeric) == TRUE)  || stop("All columns must be numeric")
  xts::is.xts(x) || stop("x must be an xts object so lags and differences are taken properly")
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
  if (is.null(demean)){
    warning("demean is NULL, setting to TRUE")
    demean = TRUE
  }
  #########
  ## END: TESTS
  #########
  
  
  if (demean == FALSE) {
    x_diff <- diff(x, 1)[2:nrow(x),]
    
    Tn <- dim(x_diff)[1]
    N <- dim(x_diff)[2]
    
    # test nfac < N
    nfac < N || stop(" nfac must be less than the number of series.")
    
    scaler <- sqrt(N) * Tn
    
    factors <- getnfac(x_diff, nfac, criteria)
    
    ic <- factors$ic
    
    fac.test <- MCMCpack::MCMCfactanal(~.,
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
    A1 <- 2
    B1 <- 1
    U1 <- (1/2)
    V1 <- (1/3)
    one <- matrix(1, I(Tn - 1), 1)
    
    Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
    
    
    A2 <- 3
    B2 <- 2
    
    output_list <- lapply(1:I(mcmc/thin), function(i){
      
      lamhat <- matrix(fac.test[i, 1:I(N * ic)], N, ic)
      
      dfhat <- matrix(fac.test[i, I(N * ic + N + 1):I((Tn) * ic + N * ic + N)], I(Tn), ic, byrow = TRUE)
      
      dehat <- x_diff - tcrossprod(dfhat, lamhat)
      
      ehat0 <- cumsum(dehat)
      lagehat0 <- xts::lag.xts(ehat0,1)[2:nrow(ehat0),]
      trim_ehat0 <- ehat0[2:nrow(ehat0),]
      
      # Do old panic
      adf30 <- adf(trim_ehat0, k1, -1)
      Pool <- pool(adfnc, t(adf30))
      adf30a <- Pool$adf31a
      adf30b <- Pool$adf31b
      
      # Model A compute rho0 (no demeaning)
      top0 <- sum(sum(lagehat0 * trim_ehat0))
      bottom0 <- sum(sum(lagehat0 * lagehat0))
      rho0 <- top0/bottom0
      res0 <- trim_ehat0 - lagehat0 * rho0
      
      # setup nuisance
      Nuisance <- nuisance(res0, 0)
      sig2 <- Nuisance$sig2
      omega2 <- Nuisance$omega2
      half <- Nuisance$half
      OMEGA2 <- mean(omega2)
      PHI4 <- mean(omega2 * omega2)
      SIG2 <- mean(sig2)
      HALF <- mean(half)
      
      
      # tests using rho- (do not project on deterministic trends)
      
      # setup constants
      ADJ <- N * Tn * HALF
      A1 <- 2
      B1 <- 1
      U1 <- (1/2)
      V1 <- (1/3)
      
      
      rho1 <- (top0 - ADJ)/bottom0
      
      # P = 0, -1 MP Tests Model A
      t_a <- scaler * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
      t_b <- scaler * (rho1 - 1) * sqrt(bottom0/(scaler^2)) * sqrt(B1 * OMEGA2/PHI4)
      
      # P = 0 , -1 PMSB test
      t_c <- sum(diag(tcrossprod(trim_ehat0)))
      t_c <- t_c / (N * Tn^2) - U1 * OMEGA2
      t_c <- sqrt(N) *  t_c /sqrt(V1 * PHI4)
      
      # Model B tests that project on constant
      
      one <- matrix(1, I(Tn - 1), 1)
      Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
      ehat <- Q_T %*% trim_ehat0
      lagehat <- Q_T %*% lagehat0
      
      top <- sum(sum(lagehat * ehat))
      bottom <- sum(sum((lagehat * lagehat)[2:nrow(lagehat0),]))
      rho1 <- top/bottom
      res1 <- ehat - lagehat * rho1
      
      # setup nuisance parameters
      Nuisance <- nuisance(res1, 0)
      sig2 <- Nuisance$sig2
      omega2 <- Nuisance$omega2
      half <- Nuisance$half
      OMEGA2 <- mean(omega2)
      PHI4 <- mean(omega2 * omega2)
      SIG2 <- mean(sig2)
      HALF <- mean(half)
      
      # setup constants
      A1 <- 3
      B1 <- 2
      ADJ <- -N * Tn * SIG2/2
      rho1 <- (top - ADJ)/bottom
      
      # Model B for P = 0
      t_a1 <- scaler * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
      t_a2 <- scaler * (rho1 - 1) * sqrt(bottom/(scaler^2)) * sqrt(B1 * OMEGA2/PHI4)
      
      c(t_a,t_b,t_a1,t_a2,t_c,rho1,adf30b)
    })
    
    mcmc_tests <- do.call(rbind,output_list)
    colnames(mcmc_tests) <- c("model_a_ta", "model_a_tb", "model_b_ta", "model_b_tb", "pmsb", "rho1","pool_adf")
    mcmc_tests <- coda::as.mcmc(mcmc_tests)

    
    
    return(mcmc_tests)
    ######################################################### 
  } else {
    
    
    # difference and demean
    x_diff <- diff(x, 1)[2:nrow(x),]
    x_diff <- scale(x_diff,center = TRUE,scale = FALSE)
    
    Tn <- dim(x_diff)[1]
    N <- dim(x_diff)[2]
    
    # test nfac < N
    nfac < N || stop(" nfac must be less than the number of series.")
    scaler <- sqrt(N) * Tn
    
    scale <- sqrt(N) * Tn
    
    factors <- getnfac(x_diff, nfac, criteria)
    
    ic <- factors$ic
    
    fac.test <- MCMCfactanal(~., factors = ic,
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
    
    
    A1 <- 36/5
    B1 <- 5/6
    U1 <- 1/6
    V1 <- 1/45
    A2 <- 15/4
    B2 <- 4
    one <- cbind(matrix(1, I(Tn - 1), 1), as.matrix(seq(1, I(Tn - 1))))
    
    Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
    
    
    output_list <- lapply(1:(mcmc/thin), function(i) {
      
      lamhat <- matrix(fac.test[i, 1:I(N * ic)], N, ic)
      
      dfhat <- matrix(fac.test[i, I(N * ic + N + 1):I((Tn) * ic + N * ic + N)], I(Tn), ic, byrow = TRUE)
      
      dehat <- x_diff - tcrossprod(dfhat, lamhat)
      
      ehat0 <- cumsum(dehat)
      lagehat0 <- xts::lag.xts(ehat0,1)[2:nrow(ehat0),]
      trim_ehat0 <- ehat0[2:nrow(ehat0),]
      
      
      # Do old panic
      
      adf31 <- adf(trim_ehat0, k1, -1)
      Pool <- pool(lm1, t(adf31))
      adf31a <- Pool$adf31a
      adf31b <- Pool$adf31b
      
      # set up for nuisance parameters
      bottom0 <- sum(sum(lagehat0 * lagehat0))
      top0 <- sum(sum(lagehat0 * trim_ehat0))
      rho0 <- top0/bottom0
      res0 <- trim_ehat0 - lagehat0 * rho0
      
      # Find nuisance parameters
      Nuisance <- nuisance(res0, 0)
      sig2 <- Nuisance$sig2
      omega2 <- Nuisance$omega2
      half <- Nuisance$half
      OMEGA2 <- mean(omega2)
      PHI4 <- mean(omega2 * omega2)
      SIG2 <- mean(sig2)
      HALF <- mean(half)
      
      
      # No longer do detrending
      
      ADJ <- SIG2/OMEGA2
      A1 <- 36/5
      B1 <- 5/6
      U1 <- 1/6
      V1 <- 1/45
      
      
      # P = 1 for Pa and Pb
      t_a <- scaler * (rho0 - 1 + ADJ * 3/Tn)
      t_a <- t_a /sqrt(A1 * PHI4 * SIG2^2/(OMEGA2^4))
      
      t_b <- scaler * (rho0 - 1 + ADJ * 3/Tn) 
      t_b <- t_b * sqrt(bottom0/(scaler^2)) 
      t_b <- t_b * sqrt(B1 * (OMEGA2^3)/(PHI4 * (SIG2^2)))
      
      
      # P = 1 PMSB
      t_c <- sum(diag(tcrossprod(trim_ehat0)))
      t_c <- t_c / (N * Tn^2)
      t_c <- t_c - U1 * OMEGA2
      t_c <- sqrt(N) * t_c / sqrt(V1 * PHI4)
      
      # Tests that project on intercept and trends
      one <- matrix(data = c(rep(1,I(Tn-1)),
                             seq(1,I(Tn-1))),
                    nrow = (Tn-1),
                    ncol = 2)
      
      Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
      
      ehat <- Q_T %*% trim_ehat0
      lagehat <- Q_T %*% lagehat0
      
      # setup for nuisance controls
      top <- sum(sum(lagehat * ehat))
      bottom <- sum(sum(lagehat * lagehat))
      rho1 <- (top)/bottom
      res1 <- ehat - lagehat * rho1
      Nuisance <- nuisance(res1, 0)
      
      sig2 <- Nuisance$sig2
      omega2 <- Nuisance$omega2
      half <- Nuisance$half
      OMEGA2 <- mean(omega2)
      PHI4 <- mean(omega2 * omega2)
      SIG2 <- mean(sig2)
      HALF <- mean(half)
      
      # setup constants
      A1 <- 15/4
      B1 <- 4
      ADJ <- -N * Tn * SIG2/2
      rho1 <- (top - ADJ)/bottom
      
      # Model C
      t_a1 <- scaler * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
      
      t_a2 <- scaler * (rho1 - 1) * sqrt(bottom/(scaler^2)) * sqrt(B1 * OMEGA2/PHI4)
      
      c(t_a, t_b, t_c, t_a1, t_a2, rho1)
    })
    
    mcmc_tests <- do.call(rbind,output_list)
    
    colnames(mcmc_tests) <- c("Pa", "Pb","PMSB", "Model C ta", "Model C tb",  "rho1")
    
    mcmc_tests <- coda::as.mcmc(mcmc_tests)
    
    
    return(mcmc_tests)
    
  }
} 