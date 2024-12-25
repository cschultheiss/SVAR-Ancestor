boot_sd <- function(Data, cons, Ahat, Bhat, u_res, ord, p, nboot, verbose = FALSE) {

  # Bootstrap to get standard errors of coefficients of reduced and structural
  # VAR using VECM with cointegration rank 3, and Cholesky
  # 
  # INPUT
  # Data: original Data 
  # cons: estimated constant term of reduced VAR
  # Ahat: list of estimated coefficients of reduced VAR
  #       Ahat[[1]] = A_1, ... Ahat[[p]] = A_p
  # Bhat: list of estimated coefficients of structural VAR
  #       Bhat[[1]] = B_0, ... Bhat[[p+1]] = B_p
  # u_res: residuals from reduced VAR
  # ord: causal order (f.ex. found by LiNGAM)
  # p: time lags
  # nboot: number of bootstrap samples
  # 
  # OUTPUT
  # prints coefficients, standard errors, p-values and an indicator of
  # significance at level 0.01 for matrices A_1, ... A_p, and B_0, B_1, ... B_p

  dims <- dim(u_res)
  t <- dims[1] # number of samples
  k <- dims[2] # number of variables

  kurt <- array(0,dim=c(nboot,p)) # to save kurtosis of new residuals
  c <- array(0,dim=c(nboot,p))

  # to save bootstrap results
  MM <- list()
  BB <- list()
  for (i in 1:p) {
    MM[[i]] <- array(0,dim=c(nboot,k^2))
    BB[[i]] <- array(0,dim=c(nboot,k^2))
  }
  BB0 <- array(0,dim=c(nboot,k^2))

  # to get exactly the same results as in paper
  set.seed(1)

  start <- proc.time()
  for (i in 1:nboot) {
    
    if (verbose && (i%%10==0)) {
      cat("bootstrap run", i, "out of", nboot, "\n")
    }

    # generating the artificial data
    ind <- sample(t, t, replace=TRUE)
    unew <- u_res[ind,] # randomly sample residuals with replacement
    Ynew <- matrix(0,nrow=(t+p),ncol=k)
    Ynew[1:p,] <- as.matrix(Data[1:p,]) # initial time points
    for(ii in (p+1):(t+p)) {
      for(j in 1:p) {
	      Y <- Ahat[[j]] %*% t(Data[ii-j,]) # generate new Data
	      Ynew[ii,] <- Ynew[ii,] + Y 
      }
      Ynew[ii,] <- cons + Ynew[ii,] + unew[ii-p,]
    }

    Ynew <- as.data.frame(Ynew)
    Data_can_new <- tsdata2canonicalform(Ynew, p)

    # estimate reduced form VAR using a vecm
    res <- VAR_estim(Data_can_new, "ols", regstats=FALSE, corank=3)

    #kurt[i,] <- kurtosis(res$resid) # kurtosis of reduced form VAR
    #c[i,] <- res$const

    # write results columnwise in rows of MM[[j]]
    for (j in 1:p) {
      MM[[j]][i,] <- res$Mhat[[j]]
    }
    
    # calculate instantaneous effects matrix using Cholesky
    B0ch <- choleski(res$resid,ord) # using causal order of original data
    
    # write results columnwise in rows of BB[[j]]
    BB0[i,] <- B0ch
    Gamma0 <- diag(k) - B0ch
    for (j in 1:p) {
      BB[[j]][i,] <- Gamma0 %*% res$Mhat[[j]]
    }
  }
  end <- proc.time()
  print(end-start)


  # standard errors of bootstrap samples

  sdMM <- list()
  sdBB <- list()
  for (i in 1:p) {
    #sdMM[[i]] <- sd(MM[[i]])
    sdMM[[i]] <- apply(MM[[i]], 2, sd)
    dim(sdMM[[i]]) <- c(k,k)
    #sdBB[[i]] <- sd(BB[[i]])
    sdBB[[i]] <- apply(BB[[i]], 2, sd)
    dim(sdBB[[i]]) <- c(k,k)
  }
  #sdBB0 <- sd(BB0)
  sdBB0 <- apply(BB0, 2, sd)
  dim(sdBB0) <- c(k,k)


  # get t-statistics
  alpha <- 0.01
  n <- dim(Data)[1]
  df <- k^2*p + k*(k-1)/2 + k # number of parameters:
    # p kxk matrices of lagged effects (full matrices estimated)
    # one kxk matrix of instantaneous effects (only lower triangle estimated)
    # one kx1 vector of constant terms
  tvalueBs <- qt(1-alpha/2, n-df-1)
  df <- k^2*p + k # number of parameters:
    # p kxk matrices of lagged effects (full matrices estimated)
    # one kx1 vector of constant terms
  tvalueMs <- qt(1-alpha/2, n-df-1)

  tstatsMM <- list()
  tstatsBB <- list()
  for (i in 1:p) {
    tstatsMM[[i]] <- Ahat[[i]]/sdMM[[i]]
    tstatsBB[[i]] <- Bhat[[i+1]]/sdBB[[i]]
    if (verbose && (i==1 | i==2)) {
      cat('\ninformation for A', i, ': coeffs, sd, pvalue, significant\n')
      print(round(Ahat[[i]],4))
      print(round(sdMM[[i]],4))
      print(2*pt(abs(tstatsMM[[i]]),n,lower.tail=FALSE))
      print(abs(tstatsMM[[i]]) > tvalueMs)
      cat('\ninformation for B', i, ': coeffs, sd, pvalue, significant\n')
      print(round(Bhat[[i+1]],4))
      print(round(sdBB[[i]],4))
      print(2*pt(abs(tstatsBB[[i]]),n,lower.tail=FALSE))
      print(abs(tstatsBB[[i]]) > tvalueBs)
    }
  }
  tstatsBB0 <- Bhat[[1]]/sdBB0
  if(verbose){
    cat('\ninformation for B0', ': coeffs, sd, pvalue, significant\n')
    print(round(Bhat[[1]],4))
    print(round(sdBB0,4))
    print(2*pt(abs(tstatsBB0),n,lower.tail=FALSE))
    print(abs(tstatsBB0) > tvalueBs)
  }
  tstatB <- cbind(tstatsBB0, do.call(cbind, tstatsBB)) 
  Bp <- 2*(1 - pt(abs(tstatB), df = n-df-1))
  Bp[is.na(Bp) & do.call(cbind, Bhat) == 0] <- 1
  Bp[is.na(Bp)] <- 0
  
  return(Bp)
}



boot_dist <- function(Data, Bhat, nboot, verbose = FALSE) {
  n <- nrow(Data) # number of samples
  p <- ncol(Data) # number of variables
  degree <- length(Bhat) - 1
  
  S0.boot <- array(NA, dim = c(p, p, nboot))
  Slag.boot <- array(NA, dim = c(p, p, nboot))

  for (i in 1:nboot) {
    # generating the artificial data
    X_new <- apply(Data, 2, function(xj) sample(xj, n, replace = TRUE))
    X_can_new <- tsdata2canonicalform(X_new, degree)
    
    res <- VARLiNGAM(X_can_new, pruning = FALSE, ntests = FALSE)
    
    # measure of instantaneous effects
    S0.boot[, , i] <- S0_f(res$Bhat[[1]], X_new)
    Slag.boot[, , i] <- Slag_f(res$Bhat, X_new)
  }
  
  return(list(S0 = S0_f(Bhat[[1]], Data), S0.boot = S0.boot, 
              Slag = Slag_f(Bhat, Data), Slag.boot = Slag.boot))
}


# measure of instantaneous effects
S0_f <- function(B0, X){
  x.var <- apply(X, 2, var)
  X.var.R <- sapply(x.var, function(v) x.var/v) 
  B0**2 * X.var.R
}

#measures how strong the total lagged causal influence
Slag_f <- function(B, X){
  degree <- length(B)
  sapply(1:ncol(X), function(j){
    tauS <- lapply(1:degree, function(tau){
      sapply(B[[tau]][, j], function(b) b * X[-c(degree-tau, nrow(X)-(tau-1):0), j])
    })
    apply(Reduce("+", tauS), 2, var) / apply(X, 2, var)
  })
}

