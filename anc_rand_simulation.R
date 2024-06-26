rm(list = ls(all = TRUE))
rand_simulation <- function(nsim = 1000, n.vec = 10^(2:6), n.init = 10^4){
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
  
  p <- 6 # number of variables
  nlag <- 1 # maximum lag considered
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42) # make it reproducible
  seed.vec <- sample(1:10000, length(n.vec))
  print(seed.vec) # 3588 3052 2252 5257 8307
  
  As <- B1s <- array(0, c(p, p, nsim)) # to store effect matrices
  pers <- matrix(NA, p, nsim) # to store permutations of the error distributions
  for (s in 1:nsim){
    rd <- randomDAG(p, 0.2, lB = 0.5, uB = 1) # create a random DAG
    
    B <- matrix(0, p , p) # represent DAG as matrix
    for (i in 2:p){
      for(j in 1:(i-1)){
        B[i,j] <- max(0, rd@edgeData@data[[paste(j,"|",i, sep="")]]$weight) # store edge weights
      }
    }
    
    pers[, s] <- sample.int(p) # random permutation
    
    A <- solve(diag(p) - B) # total effects
    for (j in 2:p){
      varj <- sum(A[j,]^2) - 1
      if(varj != 0) {
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
    B1s[,, s] <- B1
  }
  setup <- list(As = As, B1s = B1s, pers = pers)
  
  resname <- paste0("setup ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  if (save) save(setup, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  seed.n <- 0
  
  
  for (n in n.vec) {
    print(n)
    n <- n + n.init # total sample size including burnin
    seed.n <- seed.n + 1
    set.seed(seed.vec[seed.n])
    
    cl<-makeSOCKcluster(16) 
    registerDoSNOW(cl)
    tic()
    res<-foreach(gu = 1:nsim, .combine = rbind,
                 .packages = c("tsutils"), .options.snow = opts) %dorng%{
                   
                   
                   psi <- cbind(rt(n, 7) / sqrt(1.4), runif(n, -sqrt(3), sqrt(3)), rt(n, 7) / sqrt(1.4),
                                rexp(n) * (2 * rbinom(n, 1, 0.5) - 1) / sqrt(2), rnorm(n),
                                runif(n, -sqrt(3), sqrt(3))) # get noise
                   
                   # permute distribution
                   psi <- psi[, pers[, gu]]
                   
                   x <- psi
                   # generate data from SVAR
                   x[1,] <- As[, , gu] %*% psi[1,]
                   for(i in 2:n){
                     x[i,] <- As[, , gu] %*% (B1s[, , gu] %*% x[i - 1, ] + psi[i, ])
                   }
                   colnames(x) <- paste("x", 1:p, sep = "")
                   x <- x[-(1:n.init),] # discard burn-in
                   
                   laa <- lin.anc.ts(x, degree = 1) # apply ancestor regression
                   outmat <- laa[[1]] # store test statistics
                   
                   out <- list()
                   out$res <- outmat
                   out$gu <- gu
                   out                           
                 } 
    toc()
    stopCluster(cl)
    # store output list to matrix
    res.mat <- array(unlist(res[,"res"]), dim = c(p, (nlag + 1) * p, nsim), dimnames = list(rownames(res[1,"res"][[1]]),
                                                                                            colnames(res[1,"res"][[1]]), NULL))
    ind <- unlist(res[,"gu"])
    names(ind) <- NULL
    
    
    n <- n - n.init # adjust n
    # store results
    simulation <- list(res = res.mat, n = n, ind = ind,
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


