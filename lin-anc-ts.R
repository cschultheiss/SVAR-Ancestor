require(tsutils)
source("helpers.R")

lin.anc.ts <- function(x, degree, targets = colnames(x),  f  = function(x) x^3){
  # function to perform ancestor regression for SVAR
  # Input
  # x (numeric, matrix): the observational data
  # degree (integer): order of the SVAR process to be considered
  # targets (character, vector): variables whose ancestors should be estimated, all by default
  # f (function): non-linearity used for ancestor regression
  # Output
  # z.val (numeric, matrix): test statistics
  # p.val (numeric, matrix): p-values
  
  cols <- colnames(x) # all variables
  ind <- which(cols %in% targets) # considered columns
  p <- ncol(x) # number of variables
  xt <- matrix(NA, nrow = nrow(x), ncol = p * (degree +1)) # matrix to store data with lagged predictors
  for(j in 1:p){
    xtj <- lagmatrix(x[,j], 0:degree) # all lags for given variable
    xt[,j + p * (0:degree)] <- xtj # store to combined matrix
  }
  xt <- xt[complete.cases(xt),] # ignore incomplete rows
  colnames(xt) <- paste(rep(cols, degree + 1), ".", rep(0:degree, each = p), sep ="") # column names with degree information

  if(degree > 0){
    # residuals from VAR process
    u <- xt[,1:p] - xt[,-c(1:p)] %*% solve(crossprod(xt[,-c(1:p)])) %*% crossprod(xt[,-c(1:p)], xt[,1:p])
  } else {
    # original data if no time series considered
    u <- xt
  }
  n <- nrow(u) # number of available full observations
  # matrix to store test-statistics
  z.val <- matrix(NA, nrow = length(targets), ncol = p * (degree + 1), dimnames = list(targets, colnames(xt)))
  
  for (s in 0:degree){
    # loop over lags
    if(degree > 0){
      # project out data s + 1 to s + degree steps before
      us <- xt[(s+1):n,1:p] - xt[1:(n-s),-c(1:p)] %*% solve(crossprod(xt[1:(n-s),-c(1:p)]))  %*% 
        crossprod(xt[1:(n-s),-c(1:p)], xt[(s+1):n,1:p])
    } else {
      # original data if no time series considered
      us <- xt
    }
    
    colnames(us) <- colnames(x)
    wiz <- 1:p + s*p # position to store the z-value
    
    for (i in 1:length(targets)){
      # loop over target
      tari <- f(us[,targets[i]]) # apply non-linearity to target
      su <- summary(lm(tari ~u[1: (n - s), ])) # fit linear model
      z.val[i, wiz] <- su$coefficients[-1,3] # store z-values
    }
  }

  p.val <- pnorm(abs(z.val), lower.tail = FALSE)*2 # p-values
  return(list(z.val = z.val, p.val = p.val))
}


lingam.anc.ts <- function(x, degree, targets = colnames(x), n_boot = 100){
  source("VARLiNGAM/sourcedir.R")
  source("VARLiNGAM/main1.R")
  sourceDir("VARLiNGAM/", FALSE)
  sourceDir("VARLiNGAM/lingam/code", FALSE)
  # code from https://sites.google.com/site/dorisentner/publications/VARLiNGAM
  # A. Moneta, D. Entner, P.O. Hoyer, and A. Coad; Causal Inference by Independent 
  # Component Analysis: Theory and Applications (OBES 2013)
  
  # Function to perform LiNGAM for SVAR with boosting to estimate causal effect
  # as in section 6 in 
  # - A. Hyvärinen, K. Zhang, S. Shimizu, P.O. Hoyer (JMLR-2010). Estimation of
  #   a Structural Vector Autoregression Model Using Non-Gaussianity
  
  # Input
  # x (numeric, matrix): the observational data
  # degree (integer): order of the SVAR process to be considered
  # targets (character, vector): variables whose ancestors should be estimated, all by default
  # n_boot (integer): number of bootstrap samples
  # Output
  # S0: initial significance statistic of unpruned instant effects
  # S0.boot: bootstrap distribution of significance statistic of unpruned instant effects under H0
  # Slag: initial significance statistic of unpruned laged effects
  # Slag.boot: bootstrap distribution of significance statistic of unpruned instant effects under H0
  dat <- tsdata2canonicalform(x,degree)
  res <- VARLiNGAM(dat, pruning = FALSE, ntests = FALSE)
  boot_res <- boot_dist(Data = as.data.frame(x), Bhat = res$Bhat, n_boot)
  return(boot_res)
}

lingam2.anc.ts <- function(x, degree, targets = colnames(x), n_boot = 100){
  # code from https://sites.google.com/site/dorisentner/publications/VARLiNGAM
  # A. Moneta, D. Entner, P.O. Hoyer, and A. Coad; Causal Inference by Independent 
  # Component Analysis: Theory and Applications (OBES 2013)
  
  # Function to perform LiNGAM for SVAR for instant effects only
  
  # Input
  # x (numeric, matrix): the observational data
  # degree (integer): order of the SVAR process to be considered
  # targets (character, vector): variables whose ancestors should be estimated, all by default
  # n_boot (integer): number of bootstrap samples
  # Output
  # pruned instant effects
  source("VARLiNGAM/sourcedir.R")
  source("VARLiNGAM/main1.R")
  sourceDir("VARLiNGAM/", FALSE)
  dat <- tsdata2canonicalform(x,degree)
  VARres <- VAR_estim(dat, "ols", FALSE, NA, FALSE)
  return(pcalg::lingam(VARres$residuals)$Bpruned)
}

lingam3.anc.ts <- function(x, degree, targets = colnames(x), n_boot = 100){
  source("VARLiNGAM/sourcedir.R")
  source("VARLiNGAM/main1.R")
  sourceDir("VARLiNGAM/", FALSE)
  sourceDir("VARLiNGAM/lingam/code", FALSE)
  # code from https://sites.google.com/site/dorisentner/publications/VARLiNGAM
  # A. Moneta, D. Entner, P.O. Hoyer, and A. Coad; Causal Inference by Independent 
  # Component Analysis: Theory and Applications (OBES 2013)
  
  # Function to perform LiNGAM for SVAR with boosting to estimate p-values for causal effects
  
  # Input
  # x (numeric, matrix): the observational data
  # degree (integer): order of the SVAR process to be considered
  # targets (character, vector): variables whose ancestors should be estimated, all by default
  # n_boot (integer): number of bootstrap samples
  # Output
  # p-values from the bootstrapping t-distribution
  dat <- tsdata2canonicalform(x,degree)
  res <- VARLiNGAM(dat, pruning = FALSE, ntests = FALSE)
  boot_res <- boot_sd(Data = as.data.frame(x), cons = res$const, Ahat = res$Mhat, 
                      Bhat = res$Bhat, u_res = res$resid, ord = res$var_order, 
                      p = degree, nboot = n_boot)
  return(boot_res)
}
