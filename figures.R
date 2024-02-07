require(latex2exp)
source("helpers-figures.R")

# figures for randomized graph
folder <- "results/31-Jan-2024 11.53"
savefolder <- "Figures/rand-graph-Gauss"
flz <- list.files(folder)
load(paste(folder, "/", flz[1], sep = ""))
lf <- length(flz)

inst.col <- dimnames(simulation$res)[[2]][which(grepl("\\.0", dimnames(simulation$res)[[2]]))]
lag.col <- dimnames(simulation$res)[[2]][which(grepl("\\.1", dimnames(simulation$res)[[2]]))]

p <- length(inst.col)
nsim <- dim(simulation$res)[3]
j <- 4

if (any(sapply(flz, function(str) grepl("setup", str)))) {
  load(paste(folder, "/", flz[lf], sep = ""))
  As <- setup$As
  B1s <- setup$B1s
  lf <- lf - 1
}

# reference value for test
alpha <- 0.05
# which ancestors to consider, all, inst (instantaneous), lag (lagged),
# lag.dir (lagged with direct edge), lag.indir (lagged without direct edge)
mode <- "all"
# correct for not-considered as well?
all.cor <- TRUE

As[abs(As) > 1e-5] <- 1
As[abs(As) < 1e-5] <- 0

B1s[abs(B1s) > 1e-5] <- 1
B1s[abs(B1s) < 1e-5] <- 0

M1s <- B1s
Btot <- B1s
for(i in 1:(dim(B1s)[3])){
  M1s[,,i] <- As[,,i] %*% B1s[,,i]
  Btot[,,i] <- M1s[,,i] %*% As[,,i]
}

for (l in 1:p){
  As[l , l, ] <- NA
}

Asj <- t(As[j, ,])[, -j]
Btotj <- t(Btot[j, ,])
B1j <- t(B1s[j, ,])
M1j <- t(M1s[j, ,])

TAR.p <- matrix(NA, nsim + 1, lf)
mean.z <- matrix(NA, lf, 5)

alpha.perf.p <- matrix(NA, 2, lf)
i <- 0
for (file in flz[1:lf]) {
  i <- i + 1
  load(paste(folder, "/", file, sep = ""))
  z.inst <- t(simulation$res[j, inst.col,])
  z.lag <- t(simulation$res[j, lag.col,])
  
  mean.z[i, 1] <- simulation$n
  # instantaneous 
  mean.z[i, 2] <- mean(abs(z.inst[,-j])[!!Asj])
  # lagged with direct edge
  mean.z[i, 3] <- mean(abs(z.lag)[!!M1j])
  # lagged without direct edge
  mean.z[i, 4] <- mean(abs(z.lag)[!M1j & !!Btotj])
  # non-ancestor
  mean.z[i, 5] <- mean(c(abs(z.inst[,-j])[!Asj], abs(z.lag)[!Btotj]))
  
  p.inst <- 2 * pnorm(-abs(z.inst[,-j]))
  p.lag <- 2 * pnorm(-abs(z.lag))
  all.p <- cbind(p.inst, p.lag)
  all.p.adj <- t(apply(all.p, 1, holm.uncut))
  non.anc <- cbind(!Asj, !Btotj)
  {switch(mode,
         all = {
           anc <- cbind(!!Asj, !!Btotj)
           p.sub <- all.p.adj
         },
         inst = {
           if(!all.cor){
             all.p.adj <- t(apply(p.inst, 1, holm.uncut))
             non.anc <- !Asj
           }
           anc <- !!Asj
           p.sub <- all.p.adj[,colnames(all.p.adj) %in% inst.col]
         },
         lag = {
           if(!all.cor){
             all.p.adj <- t(apply(p.lag, 1, holm.uncut))
             non.anc <- !Btotj
           }
           anc <- !!Btotj
           p.sub <- all.p.adj[,colnames(all.p.adj) %in% lag.col]
         },
         lag.dir = {
           if(!all.cor){
             stop("Inversion of considered has no null guarantees")
           }
           anc <- !!B1j
           p.sub <- all.p.adj[,colnames(all.p.adj) %in% lag.col]
         },
         lag.indir = {
           if(!all.cor){
             stop("Inversion of considered has no null guarantees")
           }
           anc <- !B1j & !!Btotj
           p.sub <- all.p.adj[,colnames(all.p.adj) %in% lag.col]
         })}
  non.anc[!non.anc] <- NA
  p.min <- apply(all.p.adj * non.anc, 1, min, na.rm = TRUE)
  lims.p <- c(0, sort(p.min))
  TAR.p[,i] <- sapply(lims.p, function(lim) mean(p.sub[anc] < lim, na.rm = TRUE))
  alpha.perf.p[,i] <- c(mean(p.min < alpha), mean(p.sub[anc] < alpha, na.rm = TRUE))
}

var.ind <- c(1:p)[-j]
pp <- length(var.ind)
var.ind.label <- var.ind
labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", var.ind.label, "}$')", sep = "", collapse = ","), ")")))
labels.roc <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:6, "$')", sep = "", collapse = ","), ")")))
ord <- matrix(1:pp, nrow = 2, ncol = 3, byrow = T)
plotfac <- 4
pointfrac <- 0.8
cx <- 0.75


# png(paste(savefolder, "/z+ROC-noleg.png", sep = ""), width = 600 * plotfac,
    # height = 300 * plotfac, res = 75 * plotfac)
par(mfrow = c(1,2))
matplot(mean.z[, 1], mean.z[, -1], log ="xy", xlab = "n",
        ylab = "Average absolute z-statistics",
        pch = 1:pp, col = (1:(lf + 1))[-5], lwd = 2, las = 1)
# legend("topleft", ncol = 3, legend = labels[ord][1:pp],
# pch = (1:pp)[ord], col = (1:(lf + 1))[-5][ord], pt.lwd = 2)
which.line <- c(1:3)
for (wl in which.line){
  lines(mean.z[, 1], sqrt(mean.z[, 1]) * mean.z[4, wl + 1] / sqrt(mean.z[4, 1]), lty = 2, col = "grey")
}
abline(h = sqrt(2 / pi), lty = 2, col =" grey")

matplot((0:nsim)/nsim, TAR.p[,], type = "s", xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
        col = (1:(lf + 1))[-5], las = 1)
points(alpha.perf.p[1, ], alpha.perf.p[2, ], col = (1:(lf + 1))[-5], pch = 3)
lines(c(alpha, alpha), c(0, 1), col = "gray", lty = 2)

# legend('bottomright', col = (1:(lf + 1))[-5], ncol = 1, lwd = 2, legend = labels.roc[-lf], lty = (1:(lf + 1)))
# dev.off()

all.anc <- pmax(As, Btot) > 0
for(j in 1:p) all.anc[j, j, ] <- FALSE
dimnames(all.anc)[[1]] <- dimnames(all.anc)[[2]] <- dimnames(simulation$res)[[1]]
all.anc[] <- apply(all.anc, 3, p.to.anc)
non.anc <- !all.anc
non.anc[all.anc] <- NA
all.anc[!all.anc] <- NA

TAR <- matrix(NA, nsim + 2, 2 * lf)
alpha.ind <- integer(lf)
i <- 0
for (file in flz[1:lf]) {
  i <- i + 1
  cat("\n", i, "\n")
  load(paste(folder, "/", file, sep = ""))
  z <- simulation$res
  pv <- 2 * pnorm(-abs(z))
  sum.pv <- pv[,inst.col,]
  dimnames(sum.pv)[[2]] <- dimnames(sum.pv)[[1]]
  sum.pv[] <- apply(pv, 3, summary.p.val)
  pv.adj <- sum.pv
  pv.adj[] <- apply(sum.pv, 3, function(pv) holm.corr(pv, cut = TRUE))
  # pv.nonanc <- t(apply(pv.adj, 3, function(pv) pv[which(!ancmat)]))
  
  p.min <- pmin(apply(pv.adj * non.anc, 3, min, na.rm = TRUE))
  lims <- sort(unique(c(0, alpha, p.min)))
  alpha.ind[i] <- which(lims == alpha)
  
  stru <- stru.anc <- stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], length(lims), nsim))
  stru[] <- apply(pv.adj, 3, find.structures, lims = lims)
  for (k in 1:length(lims)){
    stru.anc[,,k,] <- stru[,,k,] * all.anc
    stru.nonanc[,,k,] <- stru[,,k,] * non.anc
  }
  pwr <- apply(stru.anc, 3, mean, na.rm = TRUE)
  FWER <- apply(apply(stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
  
  TAR[1:length(lims),c(i, lf + i)] <- c(FWER, pwr)
}

par(mfrow = c(1,1))
matplot(TAR[,1:lf], TAR[,lf + (1:lf)], type = "s",
        xlim = c(0,1), ylim = c(0,1), xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
        col = (1:p)[-5], las = 1)
points(diag(TAR[alpha.ind,1:lf]), diag(TAR[alpha.ind,lf + (1:lf)]),
       col = (1:p)[-5], pch = 3)
lines(c(0.05, 0.05), c(0, 1), col = "gray", lty = 2)

