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
#'
#'@details This function gives results for Moon and Perron tests with models A, B, and C as well as
#' the pooled tests from PANIC (2010) and the panel. A assumes no deterministic component.
#' B assumes a constant and allows for a fixed effect model. C allows a constant
#' and trend. pa-pb Pooled test from PANIC (2010). Null of nonstationarity.
#' If both reject conclude stationarity. However, if only one rejects the panel
#' is nonstationary.
#' 
#'@return rho1 Estimation of the Pooled Autoregressive Coefficient.
#'
#'@return MP.tests A data frame containing either the test statistics for the Moon and Perron model A and B tests if demeaned = FALSE
#' or a data frame containing the test statistics for the pooled tests as well as Moon and Perron's model C
#'
#'@return PMSB.tests a data frame containing the test statistic for the PMSB test, the rho coefficient, and either the 
#' LM test from Bai and Ng (2004) or the ADF test from Bai and Ng (2004)
#'
#'@return nfac An integer speciyfing the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@return k1 an integer that is the maximum lag allowed in the ADF test.
#' 
#'@return criteria a character vector with values of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'
#'@return func a character vector representing which function was run
#'
#'@return ic a numeric vector containing the number of components that were estimated
#'
#'
#'@return PMSB Unit root test tends to zero. The unit root hypothesis is rejected
#'in favor of stationarity when the PMSB test goes below a critical value.
#'
#'@references Bai, Jushan, and Serena Ng.
#''Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation.'
#' Econometric Theory 26.04 (2010): 1088-1114. Print.
#' 
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#' 
#' @export

panic10 <- function(x, nfac = NULL, k1 = NULL, criteria = NULL, demean = NULL) {
  
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
    
    # Get number of common components and PC
    factors <- getnfac(x_diff, nfac, criteria)
    
    ic <- factors$ic
    lamhat <- factors$lambda
    dfhat <- as.matrix(factors$Fhat)
    
    # Pg. 1133 Bai and Ng (2004)
    # Take cumulative sum of common factors by column
    fhat <- apply(dfhat, 2, cumsum)
    
    # If factor score d.n.e., use differenced data
    if (sum(lamhat) == 0) {
      dehat <- x_diff
    } else {
      dehat <- x_diff - dfhat %*% t(lamhat)
    }
    
    # take cumsum, lags of cumsum, and remove NA of error from PC
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
    
    # results
    test_A_B <- data.frame( mp_test = c("ta","tb"),
                            model_a = c(t_a, t_b),
                            model_b = c(t_a1,t_a2))
    
    PAC_test <- data.frame(PMSB = t_c,
                           rho1 = rho1,
                           pool_adf = adf30b)
    
    
    results <- list(MP.tests = test_A_B,
                    PMSB.tests = PAC_test,
                    k1 = k1,
                    nfac = nfac,
                    ic = factors$ic,
                    criteria = criteria,
                    func = "panic10nm")
    
    output <- structure(results,class = 'panic')
    return(output)
  } else {
    
    # difference and demean
    x_diff <- diff(x, 1)[2:nrow(x),]
    x_diff <- scale(x_diff,center = TRUE,scale = FALSE)
    
    Tn <- dim(x_diff)[1]
    N <- dim(x_diff)[2]
    
    # test nfac < N
    nfac < N || stop(" nfac must be less than the number of series.")
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
    
    # take cumsum, lags of cumsum, and remove NA of error from PC
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
    
    # results
    test.C.P <- data.frame(pool_test = c("Pa","Pb"),
                           P = c(t_a,t_b),
                           mp_test = c("ta","tb"),model_c =c(t_a1, t_a2))
    
    extra.test <- data.frame(PMSB = t_c,
                             rho1 = rho1,
                             pool_lm_04 = adf31b)
    
    results <- list(MP.tests = test.C.P,
                    PMSB.tests = extra.test,
                    k1 = k1,
                    nfac = nfac,
                    ic = factors$ic,
                    criteria = criteria,
                    func = "panic10m")
    
    output <- structure(results,class = 'panic')
    return(output)
    
  }
} 