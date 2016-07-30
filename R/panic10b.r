#'@title PANIC (2010) Sample Moment and PAC tests for Idiosyncratic Component
#'
#'@description This function performs the tests of PANIC (2010) and models A, B, and C from Moon and Perron (2004).
#' One PMSB test estimates the pooled autoregressive coefficient, and the other uses a sample
#' moment. The sample moments test is based off of the modified
#' Sargan-Bhargava test (MSB). Each test rejects after
#' the test statistic goes below the critical value of -1.64.
#'
#'@usage panic10(x, nfac, k1, criteria, demean)
#'
#'
#'@param x A NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 The maximum lag allowed in the ADF test.
#'
#'@param criteria a character vector with values of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
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
#''Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation.'
#' Econometric Theory 26.04 (2010): 1088-1114. Print.

panic10 <- function(x, nfac = NULL, k1 = NULL, criteria = NULL, demean = NULL) {
  
  # checks
  all(sapply(x, is.numeric) == TRUE)  || stop("All columns must be numeric")
  is.xts(x) || stop("x must be an xts object so lags and differences are taken properly")
  if (is.null(nfac)){
    warning("nfac is NULL, setting the maximum number of factors equal to the number of columns")
    nfac <- dim(x)[2]
  }
  if (is.null(k1)){
    warning("k1 is NULL, setting k1 equal to  k1 4 * ceiling((T/100)^(1/4))")
    k1 <- 4 * ceiling((dim(x)[1]/100)^(1/4))
  }
  if (is.null(criteria)){
    warning("criteria is NULL, setting criteria to BIC3")
    criteria <- "BIC3"
  }
  if (is.null(demean)){
    warning("demean is NULL, setting to TRUE")
    demean = TRUE
  }
    
  
    if (demean == FALSE) {
        
      #x_diff  
      x_diff <- diff(x, 1)[2:nrow(x),]
      
      Tn <- dim(x_diff)[1]
       N <- dim(x_diff)[2]
        
      scaler <- sqrt(N) * Tn
        
      factors <- getnfac(x_diff, nfac, criteria)
        
      ic <- factors$ic
        
      lamhat <- factors$lambda
        
      dfhat <- as.matrix(factors$Fhat)
        
      fhat <- apply(dfhat, 2, cumsum)
        
        if (sum(lamhat) == 0) {
            dehat <- x_diff
        } else {
            dehat <- x_diff - dfhat %*% t(lamhat)
        }
        ehat0 <- apply(dehat, 2, cumsum)
        
        lagehat0 <- trimr(lagn(ehat0, 1), 1, 0)
        
        ehat0 <- trimr(ehat0, 1, 0)
        
        # Do old panic
        
        adf30 <- adf(ehat0, k1, -1)
        
        Pool <- pool(adfnc, t(adf30))
        
        adf30a <- Pool$adf31a
        
        adf30b <- Pool$adf31b
        
        # Model A compute rho0 (no demeaning)
        
        top0 <- sum(sum(lagehat0 * ehat0))
        
        bottom0 <- sum(sum(lagehat0 * lagehat0))
        
        rho0 <- top0/bottom0
        
        res0 <- ehat0 - lagehat0 * rho0
        
        Nuisance <- nuisance(res0, 0)
        
        sig2 <- Nuisance$sig2
        
        omega2 <- Nuisance$omega2
        
        half <- Nuisance$half
        
        OMEGA2 <- mean(omega2)
        
        PHI4 <- mean(omega2 * omega2)
        
        SIG2 <- mean(sig2)
        
        HALF <- mean(half)
        
        
        # tests using rho- (do not project on deterministic trends)
        
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
        t_c <- sqrt(N) * (sum(diag(ehat0 %*% t(ehat0)))/(N * Tn^2) - U1 * OMEGA2)/sqrt(V1 * PHI4)
        
        # Model B tests that project on constant
        
        
        one <- matrix(1, I(Tn - 1), 1)
        
        Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
        
        ehat <- Q_T %*% ehat0
        
        lagehat <- Q_T %*% lagehat0
        
        top <- sum(sum(lagehat * ehat))
        
        bottom <- sum(sum(lagehat * lagehat))
        
        rho1 <- top/bottom
        
        res1 <- ehat - lagehat * rho1
        
        Nuisance <- nuisance(res1, 0)
        
        sig2 <- Nuisance$sig2
        
        omega2 <- Nuisance$omega2
        
        half <- Nuisance$half
        
        OMEGA2 <- mean(omega2)
        
        PHI4 <- mean(omega2 * omega2)
        
        SIG2 <- mean(sig2)
        
        HALF <- mean(half)
        
        A1 <- 3
        
        B1 <- 2
        
        ADJ <- -N * Tn * SIG2/2
        
        rho1 <- (top - ADJ)/bottom
        
        # Model B for P = 0
        t_a1 <- scaler * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
        
        t_a2 <- scaler * (rho1 - 1) * sqrt(bottom/(scaler^2)) * sqrt(B1 * OMEGA2/PHI4)
      
        test_A_B <- data.frame( mp_test = c("ta","tb"), model_a = c(t_a, t_b), model_b = c(t_a1,t_a2))
        
        PAC_test <- data.frame(PMSB = t_c, rho1 = rho1, pool_adf = adf30b)
        
        output <- list(MP.tests = test_A_B, PAC.tests = PAC_test)
        
        return(output)
    } else {
        
        
      #x_diff  
      x_diff <- diff(x, 1)[2:nrow(x),]
        
      x_diff <- scale(x_diff,center = TRUE,scale = FALSE)
        
      Tn <- dim(x_diff)[1]
      N <- dim(x_diff)[2]
        
      scaler <- sqrt(N) * Tn
        
      factors <- getnfac(x_diff, nfac, criteria)
        
      ic <- factors$ic1
        
      lamhat <- factors$lambda
        
      dfhat <- as.matrix(factors$Fhat)
        
      fhat <- apply(dfhat, 2, cumsum)
        
      if (sum(sum(lamhat)) == 0) {
            
          dehat <- x_diff
      } else {
            
          dehat <- x_diff - dfhat %*% t(lamhat)
      }
        
        ehat0 <- apply(dehat, 2, cumsum)
        
        lagehat0 <- trimr(lagn(ehat0, 1), 1, 0)
        
        ehat0 <- trimr(ehat0, 1, 0)
        
        # Do old panic
        
        adf31 <- adf(ehat0, k1, -1)
        
        Pool <- pool(lm1, t(adf31))
        
        adf31a <- Pool$adf31a
        
        adf31b <- Pool$adf31b
        
        bottom0 <- sum(sum(lagehat0 * lagehat0))
        
        top0 <- sum(sum(lagehat0 * ehat0))
        
        rho0 <- top0/bottom0
        
        res0 <- ehat0 - lagehat0 * rho0
        
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
        t_a <- scaler * (rho0 - 1 + ADJ * 3/Tn)/sqrt(A1 * PHI4 * SIG2^2/(OMEGA2^4))
        
        t_b <- scaler * (rho0 - 1 + ADJ * 3/Tn) * sqrt(bottom0/(scaler^2)) * sqrt(B1 * (OMEGA2^3)/(PHI4 * (SIG2^2)))
        # P = 1 PMSB
        t_c <- sqrt(N) * (sum(diag(ehat0 %*% t(ehat0)))/(N * Tn^2) - U1 * OMEGA2)/sqrt(V1 * PHI4)
        
        # Tests that project on intercept and trends
        
        one <- cbind(matrix(1, I(Tn - 1), 1), as.matrix(seq(1, I(Tn - 1))))
        
        Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
        
        ehat <- Q_T %*% ehat0
        
        lagehat <- Q_T %*% lagehat0
        
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
        
        A1 <- 15/4
        
        B1 <- 4
        
        ADJ <- -N * Tn * SIG2/2
        
        rho1 <- (top - ADJ)/bottom
        # Model C
        t_a1 <- scaler * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
        
        t_a2 <- scaler * (rho1 - 1) * sqrt(bottom/(scaler^2)) * sqrt(B1 * OMEGA2/PHI4)
        
        test.C.P <- data.frame(pool_test = c("Pa","Pb"),P = c(t_a,t_b), mp_test = c("ta","tb"),model_c =c(t_a1, t_a2))
        
        extra.test <- data.frame(PMSB = t_c,rho1 = rho1, pool_lm_04 = adf31b)
        
        output <- list(MP.tests = test.C.P, PAC.tests = extra.test)
        
        return(output)
        
    }
} 
