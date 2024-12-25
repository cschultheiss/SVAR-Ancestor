rm(list = ls(all = TRUE))
rand.dist <- function(n, type){
  switch (type,
          rt(n, 7) / sqrt(1.4), 
          runif(n, -sqrt(3), sqrt(3)), 
          rt(n, 7) / sqrt(1.4), 
          rexp(n) * (2 * rbinom(n, 1, 0.5) - 1) / sqrt(2), 
          rnorm(n), 
          runif(n, -sqrt(3), sqrt(3))
  )
}

sample.rand.dist <- function(n, type.vec = 1:6){
  sapply(type.vec, rand.dist, n = n)
}


rand_simulation <- function(nsim = 1000, n.vec = 10^(2:6), n.init = 10^4, p = 6, 
                            mc.cores = 16, h = 0){
  # when called executes the simulation for Section 4
  # Input
  # nsim (integer): number of simulation runs per sample size
  # n.vec (integer vector): sample sizes to be considered
  # n.init (integer): number of burn-in samples
  # Output (string): save location of the results folder
  

  require(tictoc)
  require(doRNG)
  require(doSNOW)
  require(parallel)
  require(MultiRNG)
  require(git2r)
  require(pcalg)
  
  source('lin-anc-ts.R', local = TRUE)
  
  commit <- revparse_single(revision = "HEAD")
  print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))
  
  save <- TRUE
  # create save location, adjust depending on folder structure
  if (save) {
    newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
    dir.create(paste("results/", newdir, sep="")) 
  }
  
  progress <- function(n, tag) {
    mod <- 16
    if (n %% mod == 0 ) {
      cat(sprintf('tasks completed: %d; tag: %d\n', n, tag))
    }
    if (n %% mod == 0 ) {
      toc()
      tic()
    }
  }
  
  opts <- list(progress = progress)
  
   # number of variables
  nlag <- 1 # maximum lag considered
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42) # make it reproducible
  seed.vec <- sample(1:10000, length(n.vec))
  print(seed.vec) # 3588 3052 2252 5257 8307
  
  As <- B1s <- array(0, c(p, p, nsim)) # to store effect matrices
  for (s in 1:nsim){
    rd <- randomDAG(p, 0.2, lB = 0.5, uB = 1) # create a random DAG
    
    B <- matrix(0, p, p) # represent DAG as matrix
    for (i in 2:p){
      for(j in 1:(i-1)){
        B[i,j] <- max(0, rd@edgeData@data[[paste(j,"|",i, sep="")]]$weight) # store edge weights
      }
    }
    
    A <- solve(diag(p) - B) # total effects
    for (j in 2:p){
      varj <- sum(A[j,]^2) - 1 
      if(!isTRUE(all.equal(varj, 0))) {
        B[j,] <- B[j,] / sqrt(varj) * runif(1, sqrt(1/2), sqrt(2)) # standardize
        A <- solve(diag(p) - B) # total effects
      }
    }
    As[, , s] <- A # store to array
    
    # random matrix for lagged effects
    B1 <- matrix(runif(p^2, 0.2, 0.8) * rbinom(p^2, 1, 0.1) * sample(c(-1, 1), p^2, TRUE), nrow = p)
    Btild <- A %*% B1 # VAR matrix
    # shrink if (almost) not stable
    eig.max <- max(abs(eigen(Btild)$values))
    if(eig.max > 0.95){
      B1 <- B1 / eig.max * 0.95
    }
    B1s[, , s] <- B1
  }
  
  if(h != 0){
    hidden <- matrix(replicate(nsim, sample(1:p, h)), ncol = h, byrow = T)
    hidden_lagged <- t(apply(hidden, 1, function(hid) c(sapply(0:nlag, function(lag) hid + p * lag))))
  }else{
    hidden <- matrix(p+1, nrow = nsim)
    hidden_lagged <- hidden
  }
  setup <- list(As = As, B1s = B1s, hidden = hidden, hidden_lagged = hidden_lagged)
  
  resname <- paste0("setup ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(setup, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  seed.n <- 0
  
  
  for (n in n.vec) {
    print(n)
    n <- n + n.init # total sample size including burn-in
    seed.n <- seed.n + 1
    set.seed(seed.vec[seed.n])
    
    cl<-makeSOCKcluster(mc.cores) 
    registerDoSNOW(cl)
    clusterExport(cl = cl, c('sample.rand.dist', 'rand.dist'))
    tic()
    res<-foreach(gu = 1:nsim, .combine = rbind,
                 .packages = c("tsutils"), .options.snow = opts) %dorng%{
                   psi <- sample.rand.dist(n, sample(1:6, p, replace = T)) # get noise
                   
                   x <- psi
                   # generate data from SVAR
                   x[1,] <- As[, , gu] %*% psi[1,]
                   for(i in 2:n){
                     x[i,] <- As[, , gu] %*% (B1s[, , gu] %*% x[i - 1, ] + psi[i, ])
                   }
                   colnames(x) <- paste("x", 1:p, sep = "")
                   x <- x[-(1:n.init),] # discard burn-in
                    
                   laa <- lin.anc.ts(x[, -hidden[gu, ]], degree = nlag) # apply ancestor regression
                   Lin_laa <- lingam.anc.ts(x[, -hidden[gu, ]], degree = nlag, n_boot = 100)
                   Lin <- lingam2.anc.ts(x[, -hidden[gu, ]], degree = nlag)
                   Lin3 <- lingam3.anc.ts(x[, -hidden[gu, ]], degree = nlag, n_boot = 100)
                   
                   if(h == 0){
                     z.val <- laa[[1]]
                     b.val <- Lin$Bhat[[1]]
                     b2.val <- Lin3
                   }else{
                    z.val <- matrix(0, nrow = p, ncol = p * (1 + nlag))
                    z.val[-hidden[gu, ], -hidden_lagged[gu, ]] <- laa[[1]]
                    colnames(z.val) <- paste0('x', 1:p, '.', rep(0:nlag, each = p))
                    row.names(z.val) <- paste0('x', 1:p)
                    
                    b.val <- matrix(0, nrow = p, ncol = p)
                    b.val[-hidden[gu, ], -hidden[gu, ]] <- Lin$Bhat[[1]]
                    colnames(b.val) <- paste0('x', 1:p)
                    row.names(b.val) <- paste0('x', 1:p)
                    
                    b2.val <- matrix(1, nrow = p, ncol = p * (1 + nlag))
                    b2.val[-hidden[gu, ], -hidden_lagged[gu, ]] <- Lin3
                    colnames(b2.val) <- paste0('x', 1:p, '.', rep(0:nlag, each = p))
                    row.names(b2.val) <- paste0('x', 1:p)
                   }
                   
                   outmat <- z.val # store test statistics
                   
                   out <- list()
                   out$res <- outmat
                   out$gu <- gu
                   out$b_res <- b.val
                   out$lingam <- Lin_laa
                   out$df <- Lin_laa$df
                   out$b2_res <- b2.val
                   out                           
                 } 
    toc()
    stopCluster(cl)
    # store output list to matrix
    res.mat <- array(unlist(res[,"res"]), dim = c(p, (nlag + 1) * p, nsim), 
                     dimnames = list(rownames(res[1,"res"][[1]]),
                                     colnames(res[1,"res"][[1]]), NULL))
    
    b2_res.mat <- array(unlist(res[,"b2_res"]), dim = c(p, (nlag + 1) * p, nsim), 
                     dimnames = list(rownames(res[1,"b2_res"][[1]]),
                                     colnames(res[1,"b2_res"][[1]]), NULL))
     
    b_res.mat <- array(unlist(res[,"b_res"]), dim = c(p, p, nsim), dimnames = list(rownames(res[1,"b_res"][[1]]),
                                                                                          colnames(res[1,"b_res"][[1]]), NULL))
    ind <- unlist(res[,"gu"])
    names(ind) <- NULL
    
    #df <- unlist(res[, "df"])
    #names(df) <- NULL
    
    n <- n - n.init # adjust n
    # store results
    simulation <- list(res = res.mat, b_res = b_res.mat, b2_res = b2_res.mat, lingam = res[, "lingam"], n = n, ind = ind,
                       r.seed = attr(res, "rng"), "commit" = commit)
    # create unique filename based on sample size and time
    resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
    # save the file to the folder
    if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
    
    As1 <- As2 <- array(NA, dim(As))
    As1[abs(As) > 1e-5] <- 1
    As2[abs(As) < 1e-5] <- 1
    # print average z-statistics for instantaneous effects for ancestors and non-ancestors
    print(apply(res.mat[,(1:p), ] * As1[, , ind], 1:2, mean, na.rm = TRUE))
    print(apply(res.mat[,(1:p), ] * As2[, , ind], 1:2, mean, na.rm = TRUE))
  }
  return(paste("results/", newdir, sep = ""))
}

rand_simulation(mc.cores = 125, p = 10, h = 1)
