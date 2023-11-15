source("lin-anc-ts.R")
nrep <- 1e3
n <- 1e2
p <- 2
deg <- 1
ptot <- p * (deg + 1)
zv <- pv <- array(NA, c(nrep, p, ptot))
lg <- array(NA, c(nrep, p, p))
for (r in 1:nrep){
  eps1 <- runif(n, -1, 1)
  eps2 <- runif(n, -1, 1)
  x1 <- x2  <- numeric(n)
  x1[1] <- eps1[1]
  x2[1] <- eps2[1]
  
  for(i in 2:n){
    x2[i] <- 0.5 * x2[i-1] + eps2[i]
    x1[i] <- 0.5 * x1[i-1] + 0.7 * x2[i] +  eps1[i]
  }
  lat <- lin.anc.ts(cbind(x1, x2), deg)
  pv[r,,] <- lat$p.val
  zv[r,,] <- lat$z.val
  X_can <- tsdata2canonicalform(cbind(x1, x2),deg) # put data into canonical form
  result <- VARLiNGAM(X_can, est_meth="ols", ntests=FALSE, pruning=TRUE)
  lg[r,,] <- result$Bhat[[1]]
}
dimnames(pv) <- dimnames(zv) <- list(1:nrep, rownames(lat$p.val), colnames(lat$p.val))

par(mfrow = c(2,2))
for(j in 1:p){
  for(k in 1:ptot){
    plot.ecdf(pv[,j,k], main = paste(dimnames(pv)[[3]][k], "on", dimnames(pv)[[2]][j]))
  }
}
apply(lg != 0, 2:3, mean)
