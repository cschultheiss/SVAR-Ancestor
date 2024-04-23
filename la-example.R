source("lin-anc-ts.R")

nrep <- 1e3 # number of simulations
n <- 1e3 # sample size
n.init <- 1e3 # burn-in
p <- 3 # number of variables
deg <- 1 # maximum lag
ptot <- p * (deg + 1)
zv <- pv <- array(NA, c(nrep, p, ptot)) # initiate array to store information

n <- n + n.init
for (r in 1:nrep){
  eps1 <- rnorm(n)
  eps1 <- sign(eps1) * abs(eps1)^1.0 # use potentices to control non-Gaussianity
  eps2 <- rnorm(n)
  eps2 <- sign(eps2) * abs(eps2)^1.5
  eps3 <- rnorm(n)
  eps3 <- sign(eps3) * abs(eps3)^1.0
  x1 <- x2  <- x3 <- numeric(n)
  # initiate time series
  x1[1] <- eps1[1]
  x2[1] <- eps2[1]
  x3[1] <- eps3[1]
  
  # simulate time series from SVAR
  for(i in 2:n){
    x1[i] <- 0.5 * x1[i-1] + eps1[i]
    x2[i] <- 0.5 * x2[i-1] + 0.35 * x1[i] + eps2[i]
    x3[i] <- 0.5 * x3[i-1] + 0.7 * x2[i-1] +  eps3[i]
  }
  xx <- cbind(x1, x2, x3)[-(1:n.init),] # discard burn-in
  lat <- lin.anc.ts(xx, deg) # run ancestor regression
  pv[r,,] <- lat$p.val # store p-values
  zv[r,,] <- lat$z.val # store z-values
}
dimnames(pv) <- dimnames(zv) <- list(1:nrep, rownames(lat$p.val), colnames(lat$p.val))
pvp <- aperm(pv, c(2, 3, 1)) # different permutation to run apply

spv <- array(NA, c(dim(summary.p.val(pvp[,,1])), nrep)) # array to store summary p-values
dimnames(spv) <- append(dimnames(summary.p.val(pvp[,,1])), list(1:nrep))
ig <- sg <- ipv <- spv # arrays to store summary and instantaneous quantities

spv[] <- apply(pvp, 3, summary.p.val) # summary p-values
ipv[] <- apply(pvp, 3, instant.p.val) # instantaneous p-values
sg[] <- apply(pvp, 3, summary.graph) # summary graph
ig[] <- apply(pvp, 3, function(p) instant.graph(p)$rec.ancs) # instantaneous effects

cat("\n Mean absolute z-statistics \n")
print(apply(abs(zv), 2:3, mean))
cat("\n Median p-values \n")
print(apply(pv, 2:3, median))
cat("\n Median summary p-values \n")
print(apply(spv, 1:2, median))
cat("\n Median instantaneous p-values \n")
print(apply(ipv, 1:2, median))
cat("\n Detection rate for summary effects \n")
print(apply(sg, 1:2, mean))
cat("\n Detection rate for instantaneous effects \n")
print(apply(ig, 1:2, mean))

