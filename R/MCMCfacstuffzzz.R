#'@title PANIC (2004) MCMC Non-Stationarity Tests on Common and Idiosyncratic Components
#'
#'@description This function performs PANIC (2010) Model C, PAC, and PMSB tests.
#' PAC estimates the pooled autoregressive coefficient, PMSB uses a sample
#' moment, and Model C performs the MP test while projecting on intercept and trend.
#'  The sample moments test is based off of the modified Sargan-Bhargava test (PMSB).
#'
#'@usage panic04(x, nfac, k1, jj)
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
#'@return adff A data frame containing ADF p value, Demeaned Error term p value,Demeaned
#'  and detrended p value,R squared for principle component
#'  , And the significance of the Error components.
#'
#'@return adf.ind A matrix containing the pooled Demeaned ADF test on the data, the
#' pooled ADF test on the common components, the pooled demeaned ADF test on the
#' Idiosyncratic component, and the pooled first differenced and demeand ADF test on the
#' Idiosyncratic component.
#' 
#'@references Bai, Jushan, and Serena Ng. 
#'"A PANIC Attack on Unit Roots and Cointegration."
#' Econometrica 72.4 (2004): 1127-1177. Print.
#'
MCMCpanic04 <- function(x, nfac, k1, jj, burn = 1000, mcmc = 10000, thin = 10, verbose = 0,
                        seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
                        a0 = 0.001, b0 = 0.001, std.var = TRUE){
x<-as.matrix(x)

Tn <- dim(x)[1]

N  <- dim(x)[2]

intx<-as.matrix(t(apply(x, 2, mean)))

repmat<-intx[rep(seq_len(nrow(intx)), each=I(nrow(x))),]

x1 <- x - repmat

x2 <- x1[2:Tn,]

dx <- trimr(mydiff(x, 1), 1, 0)

scale <- sqrt(N) * Tn

k1 <- 4 * ceiling((Tn / 100)^(1/4))


factors <- getnfac(dx, nfac, jj)


ic <- factors$ic

if (is.null(ic)){
  ic <- 1
}

PC <- pc(dx,ic)

lamhat <- PC$lambda

dfhat <- PC$fhat

fac.test<- MCMCfactanal(~., factors = ic, data = as.data.frame(dx), burnin = burn,
                        mcmc = mcmc, thin = thin, verbose = verbose, seed = seed, 
                        lambda.start = lambda.start, psi.start = psi.start, l0 = l0
                        , L0 = L0,  a0 = a0, b0 = b0, store.scores = TRUE,
                        std.var = std.var)


dehat <- NULL
lamhat<- NULL
dfhat <- NULL
ehat0 <- NULL
fhat0 <- NULL

for (i in 1:I(mcmc/thin)){
  lamhat[[i]] <- matrix(fac.test[i,1:I(N*ic)],N,ic)



  dfhat[[i]] <- matrix(fac.test[i,I(N*ic+N + 1):I((Tn-1)* ic + N*ic + N)],I(Tn-1),ic,byrow=TRUE)


  dehat[[i]] <- dx- tcrossprod(dfhat[[i]],lamhat[[i]])
}

reg   <- NULL
ehat1 <- NULL
beta1 <- NULL
fhat0 <- NULL
for (j in 1:I(mcmc/thin)){
  ehat0[[j]] <- apply(dehat[[j]],2,cumsum)

  fhat0[[j]] <- apply(dfhat[[j]],2,cumsum)
}

for (j in 1:I(mcmc/thin)){
  
reg[[j]] <- cbind( matrix(1, I(Tn-1), 1), fhat0[[i]])

ehat1[[j]] <- matrix(0, I(Tn-1), N)

beta1[[j]] <- matrix(0, I(ic+1), N)

for (i in 1:N){
  
  beta1[[j]][,i] <- qr.solve(reg[[j]], x2[,i])
  
  ehat1[[j]][,i] <- x2[,i] - reg[[j]] %*% beta1[[j]][,i]
}

}



if (ic == 1){
  p = 0
}
if (ic == 0){
  p = -1
}
if (ic > 1){
  p = 1
}
adf30 <- NULL
adf20 <- NULL
adf50 <- NULL
adf10a <- NULL
adf10b <- NULL
padf30 <- NULL
adf30a <- NULL
adf30b <- NULL
padf50 <- NULL
adf50a <- NULL
adf50b <- NULL
padf10 <- NULL

adf10 <- adf04(x1, k1, p)

padf10 <- pool(adfc2, adf10)

adf10a <- padf10$adf31a

adf10b <- padf10$adf31b

adf20 <- NULL
for (j in 1:I(mcmc/thin)){
 adf20[[j]] <- matrix(0,1,2)

for (i in 1:ic){
  
  adf20[[j]][,i] <- adf04(fhat0[[j]][,i], k1, p)
}


adf30[[j]] <- adf04(ehat0[[j]], k1, -1)     # test ehat0

adf50[[j]] <- adf04(ehat1[[j]], k1, -1)     # test ehat1

# now do the pooled test


padf30[[j]] <- pool(adfnc, adf30[[j]])

adf30a[[j]] <- padf30[[j]]$adf31a

adf30b[[j]] <- padf30[[j]]$adf31b

padf50[[j]] <- poolcoint(coint0, adf50[[j]], ic)

adf50a[[j]] <- padf50[[j]]$pvala

adf50b[[j]] <- padf50[[j]]$pvalb

}

adf20ab<- matrix(unlist(adf20), mcmc, ic, byrow=TRUE)

adf.tests <- cbind(adf50a,adf50b,adf30a,adf30b,adf20ab)

colnames(c("Demeaned Error", ))
results<- list(adf.mcmc = adf.tests)
}