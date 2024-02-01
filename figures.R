require(latex2exp)
source("lin-anc-ts.R")

# figures for randomized graph
folder <- "results/31-Jan-2024 11.53"
savefolder <- "Figures/rand-graph-Gauss"
flz <- list.files(folder)
load(paste(folder, "/", flz[1], sep = ""))
lf <- length(flz)

inst.col <- which(grepl("\\.0", dimnames(simulation$res)[[2]]))
lag.col <- which(grepl("\\.1", dimnames(simulation$res)[[2]]))

p <- length(inst.col)
nsim <- dim(simulation$res)[3]
j <- 4

if (any(sapply(flz, function(str) grepl("setup", str)))) {
  load(paste(folder, "/", flz[lf], sep = ""))
  As <- setup$As
  B1s <- setup$B1s
  lf <- lf - 1
}

alpha <- 0.05
mode <- "all"

As[abs(As) > 1e-5] <- 1
As[abs(As) < 1e-5] <- 0

B1s[abs(B1s) > 1e-5] <- 1
B1s[abs(B1s) < 1e-5] <- 0

Btot <- B1s
for(i in 1:(dim(B1s)[3])){
  Btot[,,i] <- As[,,i] %*% B1s[,,i] %*% As[,,i]
}

for (l in 1:p){
  As[l , l, ] <- NA
}

Asj <- t(As[j, ,])[, -j]
Btotj <- t(Btot[j, ,])

TAR.p <- matrix(NA, nsim + 1, lf)
mean.z <- matrix(NA, lf, 4)

alpha.perf.p <- matrix(NA, 2, lf)
i <- 0
for (file in flz[1:lf]) {
  i <- i + 1
  load(paste(folder, "/", file, sep = ""))
  z.inst <- t(simulation$res[j, inst.col,])
  z.lag <- t(simulation$res[j, lag.col,])
  
  mean.z[i, 1] <- simulation$n
  mean.z[i, 2] <- mean(abs(z.inst[,-j])[!!Asj])
  mean.z[i, 3] <- mean(abs(z.lag)[!!Btotj])
  mean.z[i, 4] <- mean(c(abs(z.inst[,-j])[!Asj], abs(z.lag)[!Btotj]))
  
  p.inst <- 2 * pnorm(-abs(z.inst[,-j]))
  p.lag <- 2 * pnorm(-abs(z.lag))
  all.p <- cbind(p.inst, p.lag)
  switch(mode,
         all = {
           all.p.adj <- t(apply(all.p, 1, holm.uncut))
           non.anc <- cbind(!Asj, !Btotj)
           anc <- cbind(!!Asj, !!Btotj)
         },
         inst = {
           all.p.adj <- t(apply(p.inst, 1, holm.uncut))
           non.anc <- !Asj
           anc <- !!Asj
         },
         lag = {
           all.p.adj <- t(apply(p.lag, 1, holm.uncut))
           non.anc <- !Btotj
           anc <- !!Btotj
         })
  non.anc[!non.anc] <- NA
  p.min <- apply(all.p.adj * non.anc, 1, min, na.rm = TRUE)
  lims.p <- c(0, sort(p.min))
  TAR.p[,i] <- sapply(lims.p, function(lim) mean(all.p.adj[anc] < lim, na.rm = TRUE))
  alpha.perf.p[,i] <- c(mean(p.min < alpha), mean(all.p.adj[anc] < alpha, na.rm = TRUE))
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
which.line <- c(1:2)
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