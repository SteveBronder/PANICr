#'@title PANIC (2010) Demeanded tests for Idiosyncratic Component
#'
#'@description This function performs PANIC (2010) Model C, PAC, and PMSB tests.
#' PAC estimates the pooled autoregressive coefficient, PMSB uses a sample
#' moment, and Model C performs the MP test while projecting on intercept and trend.
#'  The sample moments test is based off of the modified Sargan-Bhargava test (PMSB).
#'
#'@usage panic10lm(x, nfac, p, k1, jj)
#'
#'
#'@param x A NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 The maximum lag allowed in the ADF test.
#'
#'@param p A binary selection for 0 or 1. p is the order of the determinisitic
#' function in the regression. 0 is for constant only and 1 is to include a trend.
#'
#'@param jj an Integer 1 through 8. Choices 1 through 7 are respectively, IC(1),
#' IC(2), IC(3), AIC(1), BIC(1), AIC(3), and BIC(3), respectively. Choosing 8
#' makes the number of factors equal to the number of columns whose sum of
#' eigenvalues is less than  or equal to .5.
#'
#'@return rho1 Estimation of the Pooled Autoregressive Coefficient.
#'
#'@return Model This function shows MP model C. Model C assumes a constant and trend
#'
#'@return test.C.P A matrix containing the Pooled tests as well as Model C.
#'
#'@return extra.test a matrix containing t_c, rho1, and the LM test for PANIC (2004)'s
#' Idiosyncratic Components.
#'
#'@references pa-pb The pooled tests on the idiosyncratic component
#'
#'@return t_a-t_b Pooled test from PANIC (2010). Null of nonstationarity.
#' If both reject conclude stationarity. However, if only one rejects the panel
#' is nonstationary.
#'
#'@return t_a1-t_a2 Modified Moon and Perron test with constant and trend (Model C).
#' Null of nonstationarity. If both reject conclude stationarity. However,
#' if only one rejects the panel is nonstationary.
#'
#'@return PMSB Unit root test tends to zero. The unit root hypothesis is rejected
#' in favor of stationarity when the PMSB test goes below a critical value
#'
#'@references Bai, Jushan, and Serena Ng. 
#'"Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation."
#' Econometric Theory 26.04 (2010): 1088-1114. Print.

panic10lm <- function(x, nfac,p,k1,jj){


     x  <-as.matrix(x)

    dX  <- trimr(mydiff(x, 1), 1, 0)

  intdX <- as.matrix(t(apply(dX, 2, mean)))

  repmat<-intdX[rep(seq_len(nrow(intdX)), each=I(nrow(x) - 1)),]

     dx <- dX - repmat

      Tn <- dim(x)[1]

      N  <- dim(x)[2]

  scale <- sqrt(N) * Tn

factors <- getnfac(dx, nfac, jj)

     ic <- factors$ic1

 lamhat <- factors$lambda

  dfhat <- factors$Fhat

   fhat <- apply(dfhat, 2, cumsum)

    if (sum(sum(lamhat)) == 0){

        dehat <- dx
      }else{

        dehat <- dx - dfhat %*% t(lamhat)
      }

  ehat0    <- apply(dehat, 2, cumsum)

  lagehat0 <- trimr(lagn(ehat0, 1), 1, 0)

  ehat0    <- trimr(ehat0, 1, 0)

  # Do old panic

  adf31   <- adf(ehat0, k1, -1)

  Pool    <- pool(lm1, t(adf31))

  adf31a  <- Pool$adf31a

  adf31b  <- Pool$adf31b

  bottom0 <- sum(sum(lagehat0 * lagehat0))

  top0    <- sum(sum(lagehat0 * ehat0))

  rho0    <- top0 / bottom0

  res0    <- ehat0 - lagehat0 * rho0

  Nuisance<- nuisance(res0, 0)

  sig2    <- Nuisance$sig2

  omega2  <- Nuisance$omega2

  half    <- Nuisance$half

  OMEGA2  <- mean(omega2)

  PHI4    <- mean(omega2 * omega2)

  SIG2    <- mean(sig2)

  HALF    <- mean(half)


  # No longer do detrending



  ADJ <- SIG2/OMEGA2

  A1  <- 36/5

  B1  <- 5/6

  U1  <- 1/6

  V1  <- 1/45


# P = 1 for Pa and Pb
  t_a <- scale * (rho0 - 1 + ADJ*3/Tn) / sqrt(A1 * PHI4 * SIG2^2 / (OMEGA2^4))

  t_b <- scale * (rho0 - 1 + ADJ*3/Tn) * sqrt(bottom0 / (scale^2)) * sqrt(B1 * (OMEGA2^3) / (PHI4 * (SIG2^2)))
# P = 1 PMSB
  t_c <- sqrt(N)*(sum(diag(ehat0 %*% t(ehat0))) / (N*Tn^2) - U1 * OMEGA2) / sqrt(V1 * PHI4)

  # Tests that project on intercept and trends

  one     <-  cbind(matrix(1,I(Tn-2),1),as.matrix(seq(1,I(Tn-2))))

  Q_T     <- diag(I(Tn-2)) - one %*% solve(crossprod(one)) %*% t(one)

  ehat    <- Q_T %*% ehat0

  lagehat <- Q_T %*% lagehat0

  top     <- sum(sum(lagehat * ehat))

  bottom  <- sum(sum(lagehat * lagehat))

  rho1    <- (top) / bottom

  res1    <- ehat - lagehat * rho1

  Nuisance<- nuisance(res1, 0)

  sig2    <- Nuisance$sig2

  omega2  <- Nuisance$omega2

  half    <- Nuisance$half

  OMEGA2  <- mean(omega2)

  PHI4    <- mean(omega2 * omega2)

  SIG2    <- mean(sig2)

  HALF    <- mean(half)

  A1     <- 15/4

  B1     <- 4

  ADJ    <- -N * Tn * SIG2 / 2

  rho1   <- (top - ADJ) / bottom
#Model C
  t_a1   <- scale * (rho1 - 1) / sqrt(A1 * PHI4 / (OMEGA2 * OMEGA2))

  t_a2   <- scale * (rho1 - 1) * sqrt(bottom / (scale^2)) * sqrt(B1 * OMEGA2 / PHI4)



  test.C.P<- as.data.frame(matrix(c( "Pa",        t_a,     "ta",     t_a1,
                       "Pb",        t_b,     "tb",     t_a2),2,4,byrow=TRUE))

  colnames(test.C.P) <- c("Pool Test","P","MP Test", "Model C")

  extra.test <-as.data.frame(t(matrix(c(t_c, rho1, adf31b),byrow=TRUE)))

  colnames(extra.test) <- c("PMSB","rho1","04 Pool LM")

  output <- list(MP.tests = test.C.P, PAC.tests = extra.test  )


  return(output)

}
