#'@title PANIC (2004) Non-Stationarity Tests on Common and Idiosyncratic Components
#'
#'@description This function performs the pooled tests on the idiosyncratic and common component from PANIC (2004).
#'
#'
#'@usage panic04(x, nfac, k1, criteria)
#'
#'
#'@param x An object of class xts holding the time series data
#'
#'@param nfac An integer speciyfing the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 The maximum lag allowed in the ADF test.
#' 
#'@param criteria a character vector with values of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'
#'@return adff A data frame containing pooled demeaned critical values, demeaned error term critical values ,demeaned
#'  and detrended critical values ,R squared for principle component
#'  , and the significance of the error components.
#'
#'@return adf.ind A matrix containing the critical values for the pooled Demeaned ADF test, the
#' pooled ADF test on the common components, the pooled demeaned ADF test on the
#' Idiosyncratic component, and the pooled first differenced and demeaned ADF test on the
#' Idiosyncratic component.
#' 
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#'
panic04 <- function(x, nfac = NULL, k1 = NULL, criteria = NULL) {
    
  # checks
    is.xts(x) || stop("x must be an xts object so lags and differences are taken properly")

    Tn <- dim(x)[1]
    N <- dim(x)[2]
    
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
    
    #x_dm
    x_dm <- scale(x,center = TRUE,scale = FALSE)
    
    #x_trim
    x_trim <- x_dm[2:Tn, ]
    
    #dx
    x_diff <- diff(x, 1)[2:Tn,]

    scaler <- sqrt(N) * Tn

    factors <- getnfac(x_diff, nfac, criteria)
    
    
    ic <- factors$ic
    
    if (ic == 0) {
        ic <- 1
    }
    
    PC <- pc(x_diff, ic)

    if (sum(sum(PC$lambda)) == 0) {
        
        dehat <- x_diff
    } else {
        
        dehat <- x_diff - tcrossprod(PC$fhat, PC$lambda)
    }
    
    
    
    fhat0 <- apply(PC$fhat, 2, cumsum)
    
    ehat0 <- apply(dehat, 2, cumsum)
    
    ehat1 <- matrix(0, I(Tn - 1), N)
    
    reg <- cbind(matrix(1, I(Tn - 1), 1), fhat0)
    
    beta1 <- matrix(0, I(ic + 1), N)
    
    for (i in 1:N) {
        
        beta1[, i] <- qr.solve(reg, x_trim[, i])
        
        ehat1[, i] <- x_trim[, i] - reg %*% beta1[, i]
    }
    
    
    # some diagnostics to see the importance of the factors
    
    R21 <- matrix(0, N, 1)
    
    R22 <- matrix(0, N, 1)
    
    fit <- matrix(0, dim(fhat0)[1], N)
    
    fit1 <- matrix(0, I(Tn - 1), N)
    
    fit2 <- matrix(0, I(Tn - 1), N)
    
    PC$lambda <- as.matrix(PC$lambda)
    for (i in 1:N) {
        
        fit1[, i] <- fhat0 %*% PC$lambda[i, ]
        
        fit2[, i] <- PC$fhat %*% PC$lambda[i, ]
        
        R21[i, ] <- sd(dehat[, i])^2/sd(x_diff[, i])^2
        
        R22[i, ] <- sd(fit1[, i])/sd(ehat0[, i])
    }
    
    
    
    
    adf10 <- matrix(0, N, 1)
    
    adf20 <- matrix(0, ic, 1)
    
    adf30 <- matrix(0, N, 1)
    
    adf40 <- matrix(0, N, 1)
    
    adf50 <- matrix(0, N, 1)
    
    if (ic == 1) {
        p = 0
    }
    if (ic == 0) {
        p = -1
    }
    if (ic > 1) {
        p = 1
    }
    
    adf10 <- adf04(x_dm, k1, p)
    
    for (i in 1:ic) {
        
        adf20[i, ] <- adf04(fhat0[, i], k1, p)  # test fhat0 for a unit root
    }
    
    adf30 <- adf04(ehat0, k1, -1)  # test ehat0
    
    adf50 <- adf04(ehat1, k1, -1)  # test ehat1
    
    # now do the pooled test
    
    padf10 <- pool(adfc2, adf10)
    
    adf10a <- padf10$adf31a
    
    adf10b <- padf10$adf31b
    
    
    padf30 <- pool(adfnc, adf30)
    
    adf30a <- padf30$adf31a
    
    adf30b <- padf30$adf31b
    
    
    padf50 <- poolcoint(coint0, adf50, ic)
    
    adf50a <- padf50$pvala
    
    adf50b <- padf50$pvalb
    
    
    
    adfr <- data.frame(series = seq(1:N),adf = adf10,ehat = adf30, ehat1 = adf50,R2 = R21,sifF_sige= R22)

    Common <- data.frame(common_test = adf20)
    
    adf_ind <- data.frame(Pooled_demeaned = c(adf10a, adf10b), pooled_idiosyncratic = c(adf30a, adf30b), pooled_cointegration_test = c(adf50a, adf50b))
    
    results <- list(adff = adfr, pooladf = adf_ind, Common = Common)
    
    return(results)
    
} 

