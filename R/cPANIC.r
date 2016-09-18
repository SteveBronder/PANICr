#'@title PANIC (2004) Non-Stationarity Tests on Common and Idiosyncratic Components
#'
#'@description Performs the tests on the idiosyncratic and common component from PANIC (2004).
#'
#'
#'@usage panic04(x, nfac, k1, criteria)
#'
#'
#'@param x An object of class xts with each column being a time series
#'
#'@param nfac An integer speciyfing the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 an integer that is the maximum lag allowed in the ADF test.
#' 
#'@param criteria a character vector of length one with a value of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'  
#'
#'@return pooladf A data frame containing the pooled tests for the demeaned data,
#' idiosyncratic component, and the cointegration test. The first row is Fisher's method applied to
#' the p values of the respective test. The second row is the correction from PANIC (2004) applied to
#' the first row.
#' 
#' 
#'@return Common A data frame of the test results on the common component
#'
#'@return adff A data frame containing pooled demeaned critical values,
#' demeaned error term critical values, demeaned
#' and detrended critical values, R squared for principle component,
#' and the significance of the error components.
#'
#'@return nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@return k1 An integer that is the maximum lag allowed in the ADF test.
#' 
#'@return criteria A character vector with a value of either IC1, IC2, IC3, AIC1, BIC1, AIC3, BIC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'
#'@return func A character vector representing which function was run
#'
#'@return ic A numeric vector containing the number of components that were estimated
#'
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#' 
#'@export
#'
panic04 <- function(x, nfac = NULL, k1 = NULL, criteria = NULL) {
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
    x_diff <- diff(x_dm, 1)[2:Tn,]

    # Value for scaling?
    scaler <- sqrt(N) * Tn

    # approximate factor model
    factors <- getnfac(x_diff, nfac, criteria)
    
    # If no common factors, enforce one
    if (factors$ic == 0) {
      warning(" There were no estimated common factors, defaulting to 1")
        factors$ic <- 1
    }
    
    # take principle components
    PC <- pc(x_diff, factors$ic)

    # If the factor loadings sum to zero, assume independence and use x_diff
    if (sum(sum(PC$lambda)) == 0) {
        warning(" The factor loadings sum to zero, assuming independence and using x")
        dehat <- x_diff
    } else {
        
        dehat <- x_diff - tcrossprod(PC$fhat, PC$lambda)
    }
    
    # Pg. 1133 Bai and Ng (2004)
    # Take cumulative sum of common factors by column
    fhat0 <- reclass(PC$fhat, match.to = x_diff)
    fhat0 <- cumsum(fhat0)
    
    #fhat0 <- apply(PC$fhat, 2, cumsum)
    # cum sum of errors of PC by column
    ehat0 <- cumsum(dehat)
    #ehat0 <- apply(dehat, 2, cumsum)
    
    # Set up matrices for 
    # 1. Idiosyncratic components
    # 2. Regression of PC on x
    # 3. coefficients of regression
    ehat1 <- matrix(0, I(Tn - 1), N)
    reg <- cbind(matrix(1, I(Tn - 1), 1), fhat0)
    beta1 <- matrix(0, I(factors$ic + 1), N)
    
    # Do regressions for each idiosync component

    beta1 <- qr.solve(reg,x_trim)
    ehat1 <- x_trim - reg %*% beta1

    
    
    # Allocate mem
    R21 <- matrix(0, N, 1)
    R22 <- matrix(0, N, 1)
    fit <- matrix(0, dim(fhat0)[1], N)
    fit1 <- matrix(0, I(Tn - 1), N)
    fit2 <- matrix(0, I(Tn - 1), N)
    
    # the factor loadings need to be converted to matrix
    # TODO: is this necessary?
    PC$lambda <- as.matrix(PC$lambda)
    
    # common components wrt factor scores and cumsum factor scores
    # TODO: vectorize, but it's not much faster
    for (i in 1:N) {
        
        fit1[, i] <- fhat0 %*% PC$lambda[i, ]
        fit2[, i] <- PC$fhat %*% PC$lambda[i, ]
        R21[i, ] <- stats::sd(dehat[, i])^2/stats::sd(x_diff[, i])^2
        R22[i, ] <- stats::sd(fit1[, i])/stats::sd(ehat0[, i])
    }
    
    
    
    # allocate mem
    adf_dm <- matrix(0, N, 1)
    adf_fhat0 <- matrix(0, factors$ic, 1)
    adf_ehat0 <- matrix(0, N, 1)
    adf40 <- matrix(0, N, 1)
    adf_ehat1 <- matrix(0, N, 1)
    
    # p rule for approximate factor model
    if (factors$ic == 1) {
        p = 0
    }
    if (factors$ic == 0) {
        p = -1
    }
    if (factors$ic > 1) {
        p = 1
    }
    
    # test demeaned data for a unit root
    adf_dm <- adf04(x_dm, k1, p)
    
    # test fhat0 for a unit root
    for (i in 1:factors$ic) {
        
      adf_fhat0[i, ] <- adf04(fhat0[, i], k1, p)  
    }

    # test ehat0 and ehat1 for a unit root
    adf_ehat0 <- adf04(ehat0, k1, -1)  
    adf_ehat1 <- adf04(ehat1, k1, -1)  
    
    # pooled test with constant
    # adfc2 is the critical vals / p-vals with constant
    padf_dm <- pool(adfc2, adf_dm)
    adf_dma <- padf_dm$adf31a
    adf_dmb <- padf_dm$adf31b
    
    # pooled DF test with no constant
    # adfnc is the critical vals / p-vals with no constant
    padf_ehat0 <- pool(adfnc, adf_ehat0)
    adf_ehat0a <- padf_ehat0$adf31a
    adf_ehat0b <- padf_ehat0$adf31b
    
    # MQ test
    # coint0 is the Phillips-Ouliaris Zt test with a constant.
    padf_ehat1 <- poolcoint(coint0, adf_ehat1, factors$ic)
    adf_ehat1a <- padf_ehat1$pvala
    adf_ehat1b <- padf_ehat1$pvalb
    
    adfr <- data.frame(series = seq(1:N),
                       adf = adf_dm,
                       ehat = adf_ehat0,
                       ehat1 = adf_ehat1,
                       R2 = R21,
                       sifF_sige= R22)

    Common <- data.frame(Common = adf_fhat0)
    
    # NOTES:
    # pooled_demeaned: The pooled test on the demeaned data
    #  adf_dma is the fisher sum
    #  adf_dmb is the result of Thm 4 pg. 1140 of Bai and Ng (2004)
    # crit vals in common_pvals
    ###
    # pooled_idiosyncratic: The pooled test on cumsum of factor error
    #  adf_ehat0a: fisher sum
    #  adf_ehat0b: result of Thm 4 pg. 1140 of Bai and Ng (2004)
    # crit vals in ehat_pvals
    ###
    # pooled_cointegration_test: MQ test of Bai and Ng (2004)
    #  adf_ehat1a: fisher sum
    #  adf_ehat1b: result of Thm 4 pg. 1140 of Bai and Ng (2004)
    adf_ind <- data.frame(Demeaned = c(adf_dma, adf_dmb),
                          Idiosyncratic = c(adf_ehat0a, adf_ehat0b),
                          Cointegration = c(adf_ehat1a, adf_ehat1b))
    
    results <- list(pooladf = adf_ind,
                    Common = Common,
                    adff = adfr,
                    k1 = k1,
                    nfac = nfac,
                    ic = factors$ic,
                    criteria = criteria,
                    func = "panic04")
    
    output <- structure(results, class = "panic")
    return(output)
    
} 

