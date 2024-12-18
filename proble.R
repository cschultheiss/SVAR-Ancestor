p <- 6 # number of variables
nlag <- 0 # maximum lag considered

RNGkind("L'Ecuyer-CMRG")
nsim <- 1

As <- B1s <- array(0, c(p, p, nsim)) # to store effect matrices
pers <- matrix(NA, p, nsim) # to store permutations of the error distributions

rd <- randomDAG(p, 0.2, lB = 0.5, uB = 1) # create a random DAG
plot(rd)

B <- matrix(0, p, p) # represent DAG as matrix
for (i in 2:p){
  for(j in 1:(i-1)){
    B[i,j] <- max(0, rd@edgeData@data[[paste(j,"|",i, sep="")]]$weight) # store edge weights
  }
}
  
pers[, 1] <- sample.int(p) # random permutation
  
A <- solve(diag(p) - B) # total effects
for (j in 2:p){
  varj <- sum(A[j,]^2) - 1
  if(varj != 0) {
    B[j,] <- B[j,] / sqrt(varj) * runif(1, sqrt(1/2), sqrt(2)) # standardize
    A <- solve(diag(p) - B) # total effects
  }
}
As[, , 1] <- A # store to array
  
# random matrix for lagged effects
B1 <- matrix(runif(p^2, 0.2, 0.8) * rbinom(p^2, 1, 0.1) * sample(c(-1, 1), p^2, TRUE), 
             nrow = p)
Btild <- A %*% B1 # VAR matrix
# shrink if (almost) not stable
eig.max <- max(abs(eigen(Btild)$values))
if(eig.max > 0.95){
  B1 <- B1 / eig.max * 0.95
}
B1s[, , 1] <- B1

n <- 1000
n <- n + n.init # total sample size including burnin

psi <- cbind(rt(n, 7) / sqrt(1.4), runif(n, -sqrt(3), sqrt(3)), rt(n, 7) / sqrt(1.4),
                            rexp(n) * (2 * rbinom(n, 1, 0.5) - 1) / sqrt(2), rnorm(n),
                            runif(n, -sqrt(3), sqrt(3))) # get noise
               
               # permute distribution
gu <- 1
psi <- psi[, pers[, 1]]
               
x <- psi
               # generate data from SVAR
x[1,] <- As[, , gu] %*% psi[1,]
for(i in 2:n){
  x[i,] <- As[, , gu] %*% (B1s[, , gu] %*% x[i - 1, ] + psi[i, ])
}
colnames(x) <- paste("x", 1:p, sep = "")
x <- x[-(1:n.init),] # discard burn-in
               
laa <- lin.anc.ts(x, degree = 1) # apply ancestor regression

pfit <- laa$p.val < 0.0001
pfit[pfit] <- 1
pfit

dim(x)

library(pcalg)
fit <- pcalg::lingam(x)
fit
B
B1

list.files('VARLiNGAM')

source("VARLiNGAM/sourcedir.R")
source("VARLiNGAM/main1.R")
sourceDir("VARLiNGAM/", FALSE)
sourceDir("VARLiNGAM/lingam/code", FALSE)
dat <- tsdata2canonicalform(x,1)
res <- VARLiNGAM(dat, pruning = FALSE, ntests = FALSE)
res$Mhat
res$Bhat
boot_res <- boot_sd(as.data.frame(x), res$const, res$Mhat, res$Bhat, res$resid, res$var_order, 1, 100)
boot_res
res$Bhat
laa$p.val < 0.01

Data <-as.data.frame(x)
cons <- res$const
Ahat <- res$Mhat
Bhat <- res$Bhat
u_res <- res$resid 
ord <- res$var_order
p <- 1
nboot <- 10

boot_sd_firmgrowth(Growth_can, nboot, result$Mhat, result$Bhat)

# random DAGS for simulation
set.seed(1234)

library(pcalg)
p <- 6 #number of nodes
DAG <- randomDAG(p, prob = 0.5)
plot(DAG)

B <- matrix(0, p, p) # represent DAG as matrix
for (i in 2:p){
  for(j in 1:(i-1)){
    B[i,j] <- max(0, DAG@edgeData@data[[paste(j,"|",i, sep="")]]$weight) # store edge weights
  }
}
B

library(MASS)
library(mvtnorm)

# solution in terms of noise
Bprime <- ginv(diag(p) - B)

n <- 10000
N <- matrix(rnorm(n *p), ncol = p)
X <- t(Bprime %*% t(N))
X <- data.frame(X)

alpha <- 0.05

B
lingam(X)$Bpruned
lin.anc.ts(X, degree = 0)$p.val < alpha
res <- fci(list(C = cor(X), n = n), indepTest=gaussCItest, alpha = alpha, p = p)
plot(res)
summary(res)
