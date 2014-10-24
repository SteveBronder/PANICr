#'@title MCMC PANIC (2010) Sample Moment and PAC tests for Idiosyncratic Component
#'
#'@description This function performs the tests of PANIC (2010) with a Monte Carlo
#' Markov chain based on a Gibbs sampler. One test estimates the pooled autoregressive
#'  coefficient, and one uses a sample moment. The sample moments test is based off of the modified
#' Sargan-Bhargava test (PMSB) while the pooled autoregressive component is based on the
#' Moon and Perron test as well a biased corrected pooled coefficient from PANIC (2004).
#'
#'@usage MCMCpanic10(x, nfac, k1, jj, demean = FALSE, burn = 1000, mcmc = 10000, thin = 10,
#' verbose = 0, seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0,
#'   a0 = 0.001, b0 = 0.001, std.var = TRUE)
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
#'@return adf.mcmc A list of the MCMC samples of the test statistics. If demeaned is set to TRUE, adf.mcmc
#' will have the tests Pa, Pb, Model C, PMSB, and rho1. If FALSE, adf.mcmc will have Model A, Model B, 
#' PMSB, and rho. Pa, Pb, and the MP tests have a critical value of 1.96. PMSB is a degenerating critical value.
#' The critical values can be found in this packages vignette or from Stock (1990). 
#'
#'@references Bai, Jushan, and Serena Ng.
#'"Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation."
#' Econometric Theory 26.04 (2010): 1088-1114. Print.
#' 
#'@references Andrew D. Martin, Kevin M. Quinn, Jong Hee Park (2011). MCMCpack: Markov Chain Monte Carlo
#' in R. Journal of Statistical Software. 42(9): 1-21. URL http://www.jstatsoft.org/v42/i09/.
#'


MCMCpanic10<- function(x, nfac, k1, jj, demean=FALSE, burn = 1000, mcmc = 10000, thin = 10, verbose = 0,
                   seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
                   a0 = 0.001, b0 = 0.001, std.var = TRUE){
  

if (demean == FALSE){
  x  <- as.matrix(x)
  

  
  dx <- trimr(mydiff(x, 1), 1, 0)
  
  Tn <- dim(dx)[1]
  
  N  <- dim(dx)[2]
  
  scale <- sqrt(N) * Tn
  
  factors <- getnfac(dx, nfac, jj)
  
  ic <- factors$ic
  
  fac.test<- MCMCfactanal(~., factors = ic, data = as.data.frame(dx), burnin = burn,
                          mcmc = mcmc, thin = thin, verbose = verbose, seed = seed, 
                          lambda.start = lambda.start, psi.start = psi.start, l0 = l0
                          , L0 = L0,  a0 = a0, b0 = b0, store.scores = TRUE,
                          std.var = std.var)
  lamhat <- NULL
  dfhat <- NULL
  dehat <- NULL
  
  for (i in 1:I(mcmc/thin)){
    lamhat[[i]] <- matrix(fac.test[i,1:I(N*ic)],N,ic)
    
    
    
    dfhat[[i]] <- matrix(fac.test[i,I(N*ic+N + 1):I((Tn)* ic + N*ic + N)],I(Tn),ic,byrow=TRUE)
    
  }
  
  
  lagehat0 <- NULL
  ehat1 <- NULL
  beta1 <- NULL
  fhat <- NULL
  ehat0a <- NULL
  ehat0 <- NULL

  #Model A
  #compute rho0 (no demeaning)
  
  
  top0 <- NULL
  bottom0 <- NULL
  rho0 <- NULL
  rho0 <- NULL
  res0 <- NULL
  Nuisance <- NULL
  sig2 <- NULL
  omega2 <- NULL
  half <- NULL
  OMEGA2 <- NULL
  PHI4 <- NULL
  SIG2 <- NULL
  HALF <- NULL
  
  for (i in 1:I(mcmc/thin)){
  top0[[i]]     <- sum(sum(trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * trimr(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)))
  
  bottom0[[i]]  <- sum(sum(trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0)))
  
  rho0[[i]]     <- top0[[i]]/bottom0[[i]]
  
  res0[[i]]    <-  trimr(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0) - trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * rho0[[i]]
  
  Nuisance[[i]] <- nuisance(res0[[i]], 0)
  

  
  sig2[[i]]     <- Nuisance[[i]]$sig2
  
  omega2[[i]]   <- Nuisance[[i]]$omega2
  
  half[[i]]     <- Nuisance[[i]]$half
  
  OMEGA2[[i]]   <- mean(omega2[[i]])
  
  PHI4[[i]]     <- mean(omega2[[i]] * omega2[[i]])
  
  SIG2[[i]]     <- mean(sig2[[i]])
  
  HALF[[i]]     <- mean(half[[i]])
  }
  
  # tests using rho- (do not project on deterministic trends)
  
  A1  <- 2
  
  B1  <- 1
  
  U1  <- (1/2)
  
  V1  <- (1/3)
  
  res0 <- NULL
  ADJ <- NULL
  rho1 <- NULL
  t_a <- NULL
  t_b <- NULL
  t_c <- NULL
  ehat <- NULL
  
  for (i in 1:I(mcmc/thin)){
  ADJ[[i]] <- N * Tn * HALF[[i]]

  rho1[[i]] <- (top0[[i]] - ADJ[[i]]) / bottom0[[i]]
  # P = 0, -1 MP Tests Model A
  t_a[[i]] <- scale * (rho1[[i]] - 1) / sqrt(A1 * PHI4[[i]] / (OMEGA2[[i]] * OMEGA2[[i]]))
  
  t_b[[i]] <- scale * (rho1[[i]] - 1) * sqrt( bottom0[[i]] / (scale^2)) * sqrt(B1 * OMEGA2[[i]] / PHI4[[i]])
  # P = 0 , -1 PMSB test
  t_c[[i]] <- sqrt(N) * (sum(diag( trimr(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0) %*% t( trimr(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)))) / (N * Tn^2) - U1 * OMEGA2[[i]]) / sqrt(V1 * PHI4[[i]])
  }
  #Model B
  # tests that project on constant
  
  one <- NULL
  Q_T <- NULL
  
  lagehat <- NULL
  top <- NULL
  bottom <- NULL
  rho1 <- NULL
  res1 <- NULL
  Nuisance <- NULL
  sig2 <- NULL
  omega2 <- NULL
  half <- NULL
  OMEGA2 <- NULL
  PHI4 <- NULL
  SIG2 <- NULL
  HALF <- NULL
  
  one     <- matrix(1, I(Tn-1), 1)
  
  Q_T     <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
  
  for (i in 1:I(mcmc/thin)){
  
  
  
  lagehat[[i]] <- Q_T %*% trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0)
  
  top[[i]]     <- sum(sum(trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * (Q_T %*%  trimr(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0))))
  
  bottom[[i]]  <- sum(sum(trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * trimr(lagn(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0)))
  
  rho1[[i]]    <- top[[i]] / bottom[[i]]
  
  res1[[i]]    <- (Q_T %*%  trimr(apply( dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)) - lagehat[[i]] * rho1[[i]]
  
  Nuisance[[i]] <- nuisance(res1[[i]], 0)
  
  sig2[[i]]   <- Nuisance[[i]]$sig2
  
  omega2[[i]] <- Nuisance[[i]]$omega2
  
  half[[i]]   <- Nuisance[[i]]$half
  
  OMEGA2[[i]] <- mean(omega2[[i]])
  
  PHI4[[i]]   <- mean(omega2[[i]] * omega2[[i]])
  
  SIG2[[i]]   <- mean(sig2[[i]])
  
  HALF[[i]]   <- mean(half[[i]])
  
  }
  A1     <- 3
  
  B1     <- 2
  
  res1<- NULL
  ADJ <- NULL
  rho1 <- NULL
  t_a1 <- NULL
  t_a2 <- NULL
  lagehat <- NULL
  lagehat0 <- NULL
  
  for (i in 1:I(mcmc/thin)){
  ADJ[[i]]    <- -N * Tn * SIG2[[i]] / 2
  
  rho1[[i]]   <- (top[[i]] - ADJ[[i]]) / bottom[[i]]
  
  #Model B for P = 0
  t_a1[[i]] <- scale * (rho1[[i]] - 1) / sqrt(A1 * PHI4[[i]] / (OMEGA2[[i]] * OMEGA2[[i]]))
  
  t_a2[[i]] <- scale * (rho1[[i]] - 1) * sqrt(bottom[[i]] / (scale^2)) * sqrt(B1 * OMEGA2[[i]] / PHI4[[i]])
  }
  

  
  output <- cbind(t_a,t_b,t_a1,t_a2,t_c,rho1)
  
  colnames(output) <- c("model A ta","model A tb", "model B ta","model B tb", "PMSB","rho" )
 
  return(output)
  #########################################################
  #########################################################
}else{
  
  
  x  <-as.matrix(x)
  
  dx  <- trimr(mydiff(x, 1), 1, 0)
  
  intdX <- as.matrix(t(apply(dx, 2, mean)))
  
  repmat<-intdX[rep(seq_len(nrow(intdX)), each=I(nrow(x) - 1)),]
  
  Dx <- dx - repmat
  
  Tn <- dim(Dx)[1]
  
  N  <- dim(Dx)[2]
  
  scale <- sqrt(N) * Tn
  
  factors <- getnfac(dx, nfac, jj)
  
  ic <- factors$ic
  
  fac.test<- MCMCfactanal(~., factors = ic, data = as.data.frame(Dx), burnin = burn,
                          mcmc = mcmc, thin = thin, verbose = verbose, seed = seed, 
                          lambda.start = lambda.start, psi.start = psi.start, l0 = l0
                          , L0 = L0,  a0 = a0, b0 = b0, store.scores = TRUE,
                          std.var = std.var)
  lamhat <- NULL
  dfhat <- NULL
  dehat <- NULL
  
  for (i in 1:I(mcmc/thin)){
    lamhat[[i]] <- matrix(fac.test[i,1:I(N*ic)],N,ic)
    
    
    
    dfhat[[i]] <- matrix(fac.test[i,I(N*ic+N + 1):I((Tn)* ic + N*ic + N)],I(Tn),ic,byrow=TRUE)
  }
  
  lagehat0 <- NULL
  ehat1 <- NULL
  beta1 <- NULL
  fhat <- NULL
  ehat0a <- NULL
  ehat0 <- NULL
  for (i in 1:I(mcmc/thin)){
    
    lagehat0[[i]] <- trimr(lagn(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0)
    
    ehat0[[i]]  <- trimr(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)
    
  }
  
  ###
  
  
  one <- NULL
  Q_T <- NULL
  lagehat <- NULL
  top0 <- NULL
  bottom0 <- NULL
  rho0 <- NULL
  res0 <- NULL
  Nuisance <- NULL
  sig2 <- NULL
  omega2 <- NULL
  half <- NULL
  OMEGA2 <- NULL
  PHI4 <- NULL
  SIG2 <- NULL
  HALF <- NULL
  #########
  
  for (i in 1:I(mcmc/thin)){
  bottom0[[i]] <- sum(sum(trimr(lagn(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * trimr(lagn(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0)))
  
  top0[[i]]    <- sum(sum(trimr(lagn(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * trimr(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)))
  
  rho0[[i]]    <- top0[[i]] / bottom0[[i]]
  
  res0[[i]]    <- trimr(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0) - trimr(lagn(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0) * rho0[[i]]
  
  Nuisance[[i]]<- nuisance(res0[[i]], 0)
  
  sig2[[i]]    <- Nuisance[[i]]$sig2
  
  omega2[[i]]  <- Nuisance[[i]]$omega2
  
  half[[i]]   <- Nuisance[[i]]$half
  
  OMEGA2[[i]]  <- mean(omega2[[i]])
  
  PHI4[[i]]    <- mean(omega2[[i]] * omega2[[i]])
  
  SIG2[[i]]    <- mean(sig2[[i]])
  
  HALF[[i]]    <- mean(half[[i]])
  
  }
  # No longer do detrending
  
  A1  <- 36/5
  
  B1  <- 5/6
  
  U1  <- 1/6
  
  V1  <- 1/45
  
  ADJ <- NULL
  t_a <- NULL
  t_b <- NULL
  t_c <- NULL
  for (i in 1:I(mcmc/thin)){
  ADJ[[i]] <- SIG2[[i]]/OMEGA2[[i]]
  
  # P = 1 for Pa and Pb
  t_a[[i]] <- scale * (rho0[[i]] - 1 + ADJ[[i]]*3/Tn) / sqrt(A1 * PHI4[[i]] * SIG2[[i]]^2 / (OMEGA2[[i]]^4))
  
  t_b[[i]] <- scale * (rho0[[i]] - 1 + ADJ[[i]]*3/Tn) * sqrt(bottom0[[i]] / (scale^2)) * sqrt(B1 * (OMEGA2[[i]]^3) / (PHI4[[i]] * (SIG2[[i]]^2)))
  # P = 1 PMSB
  t_c[[i]] <- sqrt(N)*(sum(diag(trimr(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0) %*% t(trimr(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)))) / (N*Tn^2) - U1 * OMEGA2[[i]]) / sqrt(V1 * PHI4[[i]])
  }
  # Tests that project on intercept and trends
  
  
  ehat <- NULL
  lagehat <- NULL
  top <- NULL
  bottom <- NULL
  rho1 <- NULL
  res1 <- NULL
  Nuisance <- NULL
  sig2 <- NULL
  omega2 <- NULL
  half <- NULL
  OMEGA2 <- NULL
  PHI4 <- NULL
  SIG2 <- NULL
  HALF <- NULL
  
  one     <-  cbind(matrix(1,I(Tn-1),1),as.matrix(seq(1,I(Tn-1))))
  
  Q_T     <- diag(I(Tn-1)) - one %*% solve(crossprod(one)) %*% t(one)
  
  
  for (i in 1:I(mcmc/thin)){
  ehat[[i]]    <- Q_T %*% trimr(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1, 0)
  
  lagehat[[i]] <- Q_T %*% trimr(lagn(apply(dx- tcrossprod(dfhat[[i]],lamhat[[i]]),2,cumsum), 1), 1, 0)
  
  top[[i]]     <- sum(sum(lagehat[[i]] * ehat[[i]]))
  
  bottom[[i]]  <- sum(sum(lagehat[[i]] * lagehat[[i]]))
  
  rho1[[i]]    <- (top[[i]]) / bottom[[i]]
  
  res1[[i]]    <- ehat[[i]] - lagehat[[i]] * rho1[[i]]
  
  Nuisance[[i]]<- nuisance(res1[[i]], 0)
  
  res1[[i]]    <- NULL
  
  sig2[[i]]    <- Nuisance[[i]]$sig2
  
  omega2[[i]]  <- Nuisance[[i]]$omega2
  
  half[[i]]    <- Nuisance[[i]]$half
  
  OMEGA2[[i]]  <- mean(omega2[[i]])
  
  PHI4[[i]]   <- mean(omega2[[i]] * omega2[[i]])
  
  SIG2[[i]]    <- mean(sig2[[i]])
  
  HALF[[i]]    <- mean(half[[i]])
  }
  A1     <- 15/4
  
  B1     <- 4
  
  ADJ <- NULL
  rho1 <- NULL
  t_a1 <- NULL
  t_a2 <- NULL
  
  for (i in 1:I(mcmc/thin)){
  ADJ[[i]]    <- -N * Tn * SIG2[[i]] / 2
  
  rho1[[i]]   <- (top[[i]] - ADJ[[i]]) / bottom[[i]]
  #Model C
  t_a1[[i]]   <- scale * (rho1[[i]] - 1) / sqrt(A1 * PHI4[[i]] / (OMEGA2[[i]] * OMEGA2[[i]]))
  
  t_a2[[i]]   <- scale * (rho1[[i]] - 1) * sqrt(bottom[[i]] / (scale^2)) * sqrt(B1 * OMEGA2[[i]] / PHI4[[i]])
  
  }
  

  adf.mcmc <- cbind(t_a,t_b,t_a1,t_a2,t_c,rho1)

  colnames(adf.mcmc) <- c("Pa","Pb","Model C ta","Model C tb","PMSB","rho1")
 
  output <- list(adf.mcmc = adf.mcmc  )
  
  
  return(output)
  
}
}