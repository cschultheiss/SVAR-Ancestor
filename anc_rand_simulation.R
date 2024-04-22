rm(list = ls(all = TRUE))
rand_simulation <- function(nsim = 1000, n.vec = 10^(2:6), n.init = 10^4){
  # when called executes the simulation for Section 4
  # Input
  # nsim (integer): number of simulation runs per sample size
  # n.vec (integer vector): sample sizes to be considered
  # n.init (integer): number of burn-in samples
  # Output (string): save location of the results folder
  
  require(MASS)
  require(hdi)
  require(Matrix)
  require(tictoc)
  require(doRNG)
  require(doSNOW)
  require(parallel)
  require(MultiRNG)
  require(git2r)
  require(expm)
  require(pcalg)
  
  source('lin-anc-ts.R')
  
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
  
  p <- 6
  nlag <- 1
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42)
  seed.vec <- sample(1:10000, length(n.vec))
  print(seed.vec) # 3588 3052 2252 5257 8307
  
  As <- B1s <- array(0, c(p, p, nsim))
  pers <- matrix(NA, p, nsim)
  for (s in 1:nsim){
    rd <- randomDAG(p, 0.2, lB = 0.5, uB = 1)
    
    B <- matrix(0, p , p)
    for (i in 2:p){
      for(j in 1:(i-1)){
        B[i,j] <- max(0, rd@edgeData@data[[paste(j,"|",i, sep="")]]$weight)
      }
    }
    
    pers[, s] <- sample.int(p)
    
    A <- solve(diag(p) - B)
    for (j in 2:p){
      varj <- sum(A[j,]^2) - 1
      if(varj != 0) {
        B[j,] <- B[j,] / sqrt(varj) * runif(1, sqrt(1/2), sqrt(2))
        A <- solve(diag(p) - B)
      }
    }
    As[, , s] <- A
    
    B1 <- matrix(runif(p^2, 0.2, 0.8) * rbinom(p^2, 1, 0.1) * sample(c(-1, 1), p^2, TRUE), nrow = p)
    Btild <- A %*% B1
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
    n <- n + n.init
    seed.n <- seed.n + 1
    set.seed(seed.vec[seed.n])
    
    cl<-makeSOCKcluster(16) 
    registerDoSNOW(cl)
    tic()
    res<-foreach(gu = 1:nsim, .combine = rbind,
                 .packages = c("MASS", "Matrix", "hdi", "MultiRNG", "tictoc", "pcalg", "tsutils"), .options.snow = opts) %dorng%{
                   
                   
                   psi <- cbind(rt(n, 7) / sqrt(1.4), runif(n, -sqrt(3), sqrt(3)), rt(n, 7) / sqrt(1.4),
                                rexp(n) * (2 * rbinom(n, 1, 0.5) - 1) / sqrt(2), rnorm(n),
                                runif(n, -sqrt(3), sqrt(3)))
                   # permute distribution
                   psi <- psi[, pers[, gu]]
                   
                   x <- psi
                   laa <- list()
                   x[1,] <- As[, , gu] %*% psi[1,]
                   for(i in 2:n){
                     x[i,] <- As[, , gu] %*% (B1s[, , gu] %*% x[i - 1, ] + psi[i, ])
                   }
                   colnames(x) <- paste("x", 1:p, sep = "")
                   x <- x[-(1:n.init),]
                   
                   laa <- lin.anc.ts(x, degree = 1)
                   outmat <- laa[[1]]
                   
                   
                   
                   out <- list()
                   out$res <- outmat
                   out$gu <- gu
                   out                           
                 } 
    toc()
    stopCluster(cl)
    res.mat <- array(unlist(res[,"res"]), dim = c(p, (nlag + 1) * p, nsim), dimnames = list(rownames(res[1,"res"][[1]]),
                                                                                            colnames(res[1,"res"][[1]]), NULL))
    ind <- unlist(res[,"gu"])
    names(ind) <- NULL
    
    
    n <- n - n.init
    simulation <- list(res = res.mat, n = n, ind = ind,
                       r.seed = attr(res, "rng"), "commit" = commit)
    resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
    if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
    
    As1 <- As2 <- array(NA, dim(As))
    As1[abs(As) > 1e-5] <- 1
    As2[abs(As) < 1e-5] <- 1
    for (r in 0:0) {
      print(apply(res.mat[,(1:p) + (r * p), ] * As1[, , ind], 1:2, mean, na.rm = TRUE))
      print(apply(res.mat[,(1:p) + (r * p), ] * As2[, , ind], 1:2, mean, na.rm = TRUE))
    }
  }
  return(paste("results/", newdir, sep = ""))
}



