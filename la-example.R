source("lin-anc-ts.R")
source("source-LiNGAM.R")
nrep <- 1e3
n <- 1e3
p <- 3
deg <- 1
ptot <- p * (deg + 1)
zv <- pv <- array(NA, c(nrep, p, ptot))
# lg <- array(NA, c(nrep, p, p))
for (r in 1:nrep){
  eps1 <- rnorm(n)
  eps1 <- sign(eps1) * abs(eps1)^1.0
  eps2 <- rnorm(n)
  eps2 <- sign(eps2) * abs(eps2)^1.5
  eps3 <- rnorm(n)
  eps3 <- sign(eps3) * abs(eps3)^1.0
  x1 <- x2  <- x3 <- numeric(n)
  x1[1] <- eps1[1]
  x2[1] <- eps2[1]
  x3[1] <- eps3[1]
  
  for(i in 2:n){
    x1[i] <- 0.5 * x1[i-1] + eps1[i]
    x2[i] <- 0.5 * x2[i-1] + 0.35 * x1[i] + eps2[i]
    x3[i] <- 0.5 * x3[i-1] + 0.7 * x2[i-1] +  eps3[i]
  }
  xx <- cbind(x1, x2, x3)
  lat <- lin.anc.ts(xx, deg)
  pv[r,,] <- lat$p.val
  zv[r,,] <- lat$z.val
  # X_can <- tsdata2canonicalform(xx,deg) # put data into canonical form
  # result <- VARLiNGAM(X_can, est_meth="ols", ntests=FALSE, pruning=TRUE)
  # lg[r,,] <- result$Bhat[[1]]
}
dimnames(pv) <- dimnames(zv) <- list(1:nrep, rownames(lat$p.val), colnames(lat$p.val))
pvp <- aperm(pv, c(2, 3, 1))

spv <- array(NA, c(dim(summary.p.val(pvp[,,1])), nrep))
dimnames(spv) <- append(dimnames(summary.p.val(pvp[,,1])), list(1:nrep))
ig <- sg <- spv
spv[] <- apply(pvp, 3, summary.p.val)
sg[] <- apply(pvp, 3, summary.graph)
ig[] <- apply(pvp, 3, function(p) instant.graph(p)$rec.ancs)
# par(mfrow = c(p, deg + 1))
# for(j in 1:p){
#   for(k in 1:ptot){
#     plot.ecdf(pv[,j,k], main = paste(dimnames(pv)[[3]][k], "on", dimnames(pv)[[2]][j]))
#   }
# }

apply(pv, 2:3, median)
apply(spv, 1:2, median)
apply(sg, 1:2, mean)
apply(ig, 1:2, mean)

# apply(lg != 0, 2:3, mean)
