#'@title MCMC PANIC (2010) Sample Moment and PAC tests for Idiosyncratic Component
#'
#'@description This function performs the tests of PANIC (2010) with a Monte Carlo
#' Markov chain based on a gibbs sampler. One test estimates the pooled autoregressive
#'  coefficient, and one uses a sample moment. The sample moments test is based off of the modified
#' Sargan-Bhargava test (PMSB).
#'
#'@usage panic10(x, nfac, k1, jj, demean, burn = 1000, mcmc = 10000, thin = 10,
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
#'@return rho1 Estimation of the Pooled Autoregressive Coefficient.
#'
#'@return Model This function shows MP models A, B, and C. A assumes no deterministic component.
#' B assumes a constant and allows for a fixed effect model. C allows a constant
#' and trend.
#'
#'@return test.A.B A matrix containing t_a, t_b, t_a1, and t_a2.
#'
#'@return test.C.P A matrix containing the Pooled tests as well as Model C.
#'
#'@return extra.test a matrix containing t_c, rho1, adf30b.
#' Or a matrix containing t_c, rho1, and the LM test for PANIC (2004)'s
#' Idiosyncratic Components.
#'
#'
#'@return pa-pb Pooled test from PANIC (2010). Null of nonstationarity.
#' If both reject conclude stationarity. However, if only one rejects the panel
#' is nonstationary.
#'
#'
#'@return PMSB Unit root test tends to zero. The unit root hypothesis is rejected
#'in favor of stationarity when the PMSB test goes below a critical value.
#'
#'@references Bai, Jushan, and Serena Ng.
#'"Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation."
#' Econometric Theory 26.04 (2010): 1088-1114. Print.


MCMCpanic10<- function(x, nfac, k1, jj, demean, burn = 1000, mcmc = 10000, thin = 10, verbose = 0,
                   seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
                   a0 = 0.001, b0 = 0.001, std.var = TRUE){
  

if (demean == FALSE){
  x  <- as.matrix(x)
  
  Tn <- dim(x)[1]
  
  N  <- dim(x)[2]
  
  dx <- trimr(mydiff(x, 1), 1, 0)
  
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
    
    
    
    dfhat[[i]] <- matrix(fac.test[i,I(N*ic+N + 1):I((Tn-1)* ic + N*ic + N)],I(Tn-1),ic,byrow=TRUE)
    
    
    dehat[[i]] <- dx- tcrossprod(dfhat[[i]],lamhat[[i]])
  }
  
  lagehat0 <- NULL
  ehat1 <- NULL
  beta1 <- NULL
  fhat <- NULL
  ehat0a <- NULL
  ehat0 <- NULL
  for (j in 1:I(mcmc/thin)){
    ehat0a[[j]] <- apply(dehat[[j]],2,cumsum)
    
    fhat[[j]] <- apply(dfhat[[j]],2,cumsum)
  
  
  
  
  lagehat0[[j]] <- trimr(lagn(ehat0a[[j]], 1), 1, 0)
  
  ehat0[[j]]  <- trimr(ehat0a[[j]], 1, 0)
  
  }
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
  top0[[i]]     <- sum(sum(lagehat0[[i]] * ehat0[[i]]))
  
  bottom0[[i]]  <- sum(sum(lagehat0[[i]] * lagehat0[[i]]))
  
  rho0[[i]]     <- top0[[i]]/bottom0[[i]]
  
  res0[[i]]    <- ehat0[[i]] - lagehat0[[i]] * rho0[[i]]
  
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
  
  ADJ <- NULL
  rho1 <- NULL
  t_a <- NULL
  t_b <- NULL
  t_c <- NULL
  
  for (i in 1:I(mcmc/thin)){
  ADJ[[i]] <- N * Tn * HALF[[i]]

  rho1[[i]] <- (top0[[i]] - ADJ[[i]]) / bottom0[[i]]
  # P = 0, -1 MP Tests Model A
  t_a[[i]] <- scale * (rho1[[i]] - 1) / sqrt(A1 * PHI4[[i]] / (OMEGA2[[i]] * OMEGA2[[i]]))
  
  t_b[[i]] <- scale * (rho1[[i]] - 1) * sqrt( bottom0[[i]] / (scale^2)) * sqrt(B1 * OMEGA2[[i]] / PHI4[[i]])
  # P = 0 , -1 PMSB test
  t_c[[i]] <- sqrt(N) * (sum(diag(ehat0[[i]] %*% t(ehat0[[i]]))) / (N * Tn^2) - U1 * OMEGA2[[i]]) / sqrt(V1 * PHI4[[i]])
  }
  #Model B
  # tests that project on constant
  
  one <- NULL
  Q_T <- NULL
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
  
  one     <- matrix(1, I(Tn-2), 1)
  
  Q_T     <- diag(I(Tn - 2)) - one %*% solve(crossprod(one)) %*% t(one)
  
  for (i in 1:I(mcmc/thin)){
  
  
  ehat[[i]]    <- Q_T %*% ehat0[[i]]
  
  lagehat[[i]] <- Q_T %*% lagehat0[[i]]
  
  top[[i]]     <- sum(sum(lagehat[[i]] * ehat[[i]]))
  
  bottom[[i]]  <- sum(sum(lagehat[[i]] * lagehat[[i]]))
  
  rho1[[i]]    <- top[[i]] / bottom[[i]]
  
  res1[[i]]    <- ehat[[i]] - lagehat[[i]] * rho1[[i]]
  
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
  
  ADJ <- NULL
  rho1 <- NULL
  t_a1 <- NULL
  t_a2 <- NULL
  
  for (i in 1:I(mcmc/thin)){
  ADJ[[i]]    <- -N * Tn * SIG2[[i]] / 2
  
  rho1[[i]]   <- (top[[i]] - ADJ[[i]]) / bottom[[i]]
  
  #Model B for P = 0
  t_a1[[i]] <- scale * (rho1[[i]] - 1) / sqrt(A1 * PHI4[[i]] / (OMEGA2[[i]] * OMEGA2[[i]]))
  
  t_a2[[i]] <- scale * (rho1[[i]] - 1) * sqrt(bottom[[i]] / (scale^2)) * sqrt(B1 * OMEGA2[[i]] / PHI4[[i]])
  }
  
  output <- cbind(t_a,t_a1,t_a2,t_b,t_c,rho1)
  
 
  return(output)
}else{
  
  
  x  <-as.matrix(x)
  
  dX  <- trimr(mydiff(x, 1), 1, 0)
  
  intdX <- as.matrix(t(apply(dX, 2, mean)))
  
  repmat<-intdX[rep(seq_len(nrow(intdX)), each=I(nrow(x) - 1)),]
  
  Dx <- dX - repmat
  
  Tn <- dim(x)[1]
  
  N  <- dim(x)[2]
  
  scale <- sqrt(N) * Tn
  
  factors <- getnfac(dx, nfac, jj)
  
  ic <- factors$ic1
  
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
    
    
    
    dfhat[[i]] <- matrix(fac.test[i,I(N*ic+N + 1):I((Tn-1)* ic + N*ic + N)],I(Tn-1),ic,byrow=TRUE)
    
    
    dehat[[i]] <- dx- tcrossprod(dfhat[[i]],lamhat[[i]])
  }
  
  lagehat0 <- NULL
  ehat1 <- NULL
  beta1 <- NULL
  fhat <- NULL
  ehat0a <- NULL
  ehat0 <- NULL
  for (j in 1:I(mcmc/thin)){
    ehat0a[[j]] <- apply(dehat[[j]],2,cumsum)
    
    fhat[[j]] <- apply(dfhat[[j]],2,cumsum)
    
    
    
    
    lagehat0[[j]] <- trimr(lagn(ehat0a[[j]], 1), 1, 0)
    
    ehat0[[j]]  <- trimr(ehat0a[[j]], 1, 0)
    
  }
  ###
  
  one <- NULL
  Q_T <- NULL
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
  #########
  
  for (i in 1:I(mcmc/thin)){
  bottom0[[i]] <- sum(sum(lagehat0[[i]] * lagehat0[[i]]))
  
  top0[[i]]    <- sum(sum(lagehat0[[i]] * ehat0[[i]]))
  
  rho0[[i]]    <- top0[[i]] / bottom0[[i]]
  
  res0[[i]]    <- ehat0[[i]] - lagehat0[[i]] * rho0[[i]]
  
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
  t_c[[i]] <- sqrt(N)*(sum(diag(ehat0[[i]] %*% t(ehat0[[i]]))) / (N*Tn^2) - U1 * OMEGA2[[i]]) / sqrt(V1 * PHI4[[i]])
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
  
  one     <-  cbind(matrix(1,I(Tn-2),1),as.matrix(seq(1,I(Tn-2))))
  
  Q_T     <- diag(I(Tn-2)) - one %*% solve(crossprod(one)) %*% t(one)
  
  
  for (i in 1:I(mcmc/thin)){
  ehat[[i]]    <- Q_T %*% ehat0[[i]]
  
  lagehat[[i]] <- Q_T %*% lagehat0[[i]]
  
  top[[i]]     <- sum(sum(lagehat[[i]] * ehat[[i]]))
  
  bottom[[i]]  <- sum(sum(lagehat[[i]] * lagehat[[i]]))
  
  rho1[[i]]    <- (top[[i]]) / bottom[[i]]
  
  res1[[i]]    <- ehat[[i]] - lagehat[[i]] * rho1[[i]]
  
  Nuisance[[i]]<- nuisance(res1[[i]], 0)
  
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
  t_a1[[i]]   <- scale * (rho1 - 1) / sqrt(A1 * PHI4[[i]] / (OMEGA2[[i]] * OMEGA2[[i]]))
  
  t_a2[[i]]   <- scale * (rho1[[i]] - 1) * sqrt(bottom[[i]] / (scale^2)) * sqrt(B1 * OMEGA2[[i]] / PHI4[[i]])
  
  }
  

  adf.test <- cbind(t_a,t_b,t_a1,t_a2,t_c,rho1)

 
  output <- list(adf.tests = adf.test  )
  
  
  return(output)
  
}
}