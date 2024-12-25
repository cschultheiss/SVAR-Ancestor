require(latex2exp)
source("helpers-figures.R")

# figures for randomized graph

folder <- "results/25-Dez-2024 15.34"
savefolder <- "Figures/test"

one_target <- function(folder, j = 4, alpha = 0.05, mode = "all", all.cor = TRUE){
  # wrapper function to generate plots from Section 2
  # Input
  # folder (string): location where the result files are stored
  # alpha (numeric): reference value for test
  # mode (string): which ancestors to consider, all, inst (instantaneous), lag (lagged),
  # lag.dir (lagged with direct edge), lag.indir (lagged without direct edge)
  # all.cor (boolean): correct for not-considered as well?
  flz <- list.files(folder)
  load(paste(folder, "/", flz[1], sep = ""))
  lf <- length(flz)
  
  # columns representing instantaneous effects
  inst.col <- dimnames(simulation$res)[[2]][which(grepl("\\.0", dimnames(simulation$res)[[2]]))]
  # columns representing lagged effects
  lag.col <- dimnames(simulation$res)[[2]][which(grepl("\\.1", dimnames(simulation$res)[[2]]))]
  
  p <- length(inst.col) # number of covariates
  nsim <- dim(simulation$res)[3] # number of simulations
  
  if (any(sapply(flz, function(str) grepl("setup", str)))) {
    load(paste(folder, "/", flz[lf], sep = ""))
    # read of causal connections
    As <- setup$As
    B1s <- setup$B1s
    lf <- lf - 1
  }
  
  # encode as binary
  As[abs(As) > 1e-5] <- 1
  As[abs(As) < 1e-5] <- 0
  
  B1s[abs(B1s) > 1e-5] <- 1
  B1s[abs(B1s) < 1e-5] <- 0
  
  M1s <- B1s # lagged effect starting with instantaneous edge
  Btot <- B1s # every way of lagged effect
  for(i in 1:(dim(B1s)[3])){
    M1s[,,i] <- As[,,i] %*% B1s[,,i]
    Btot[,,i] <- M1s[,,i] %*% As[,,i]
  }
  
  # no self-effects
  for (l in 1:p){
    As[l , l, ] <- NA
  }
  
  # effects relevant for target
  Asj <- t(As[j, ,])[, -j]
  Btotj <- t(Btot[j, ,])
  B1j <- t(B1s[j, ,])
  M1j <- t(M1s[j, ,])
  
  TAR.p <- matrix(NA, nsim + 1, 2 * lf) # matrix to store detection rate
  mean.z <- matrix(NA, lf, 5) # matrix to store average z-statistics
  
  alpha.perf.p <- matrix(NA, 2, lf) # matrix to store performance at level alpha
  i <- 0
  for (file in flz[1:lf]) {
    # loop over sample sizes
    i <- i + 1
    load(paste(folder, "/", file, sep = ""))
    # z-statistics
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
    
    # p-values
    p.inst <- 2 * pnorm(-abs(z.inst[,-j]))
    p.lag <- 2 * pnorm(-abs(z.lag))
    all.p <- cbind(p.inst, p.lag)
    all.p.adj <- t(apply(all.p, 1, holm.corr))
    non.anc <- cbind(!Asj, !Btotj)
    {switch(mode,
            all = {
              # consider all ancestors
              anc <- cbind(!!Asj, !!Btotj)
              p.sub <- all.p.adj
            },
            inst = {
              if(!all.cor){
                # multiplicity correction only over instantaneous
                all.p.adj <- t(apply(p.inst, 1, holm.corr))
                non.anc <- !Asj
              }
              anc <- !!Asj
              p.sub <- all.p.adj[,colnames(all.p.adj) %in% inst.col]
            },
            lag = {
              if(!all.cor){
                # multiplicity correction only over lagged
                all.p.adj <- t(apply(p.lag, 1, holm.corr))
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
    # minimum p-value leading to wrong rejection
    p.min <- apply(all.p.adj * non.anc, 1, min, na.rm = TRUE)
    lims.p <- c(0, sort(p.min))
    # power and FWER at each limit
    pwr <- sapply(lims.p, function(lim) mean(p.sub[anc] < lim, na.rm = TRUE))
    FWER <- sapply(lims.p, function(lim) mean(lims.p < lim))
    TAR.p[,c(i , lf + i)] <- c(FWER, pwr) # ROC
    # power and FWER at alpha
    alpha.perf.p[,i] <- c(mean(p.min < alpha), mean(p.sub[anc] < alpha, na.rm = TRUE))
  }
  
  var.ind <- c(1:p)[-j]
  pp <- length(var.ind)
  var.ind.label <- var.ind
  labels <- eval(parse(text = paste("c(", paste("TeX('$X_{", var.ind.label, "}$')", sep = "", collapse = ","), ")")))
  labels.roc <- eval(parse(text = paste("c(", paste("TeX('$n=10^", 2:6, "$')", sep = "", collapse = ","), ")")))
  ord <- matrix(1:pp, nrow = 2, ncol = 3, byrow = T)
  
  par(mfrow = c(1,2))
  # plot z-statistics
  matplot(mean.z[, 1], mean.z[, -1], log ="xy", xlab = "T",
          ylab = "Average absolute z-statistics",
          pch = 1:pp, col = (1:(lf + 1))[-5], lwd = 2, las = 1)

  which.line <- c(1:3)
  for (wl in which.line){
    # plot root-n growth if applicable.
    lines(mean.z[, 1], sqrt(mean.z[, 1]) * mean.z[4, wl + 1] / sqrt(mean.z[4, 1]), lty = 2, col = "grey")
  }
  abline(h = sqrt(2 / pi), lty = 2, col =" grey") # absolute moment of standard normal
  
  # plot ROC
  matplot(TAR.p[, 1:lf], TAR.p[,lf + (1:lf)], type = "s", xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
          col = (1:(lf + 1))[-5], las = 1, xlim = c(0, max(TAR.p[,1:lf])), ylim = c(0, 1))
  # add performance at alpha
  points(alpha.perf.p[1, ], alpha.perf.p[2, ], col = (1:(lf + 1))[-5], pch = 3)
  # plot target alpha
  lines(c(alpha, alpha), c(0, 1), col = "gray", lty = 2)
}

plotfac <- 4
# png(paste(savefolder, "/z+ROC-noleg.png", sep = ""), width = 600 * plotfac,
#     height = 300 * plotfac, res = 75 * plotfac)
# one_target(folder)
# dev.off()

network <- function(folder, alpha = 0.05){
  # wrapper function to generate plots from Section 2
  # Input
  # folder (string): location where the result files are stored
  # alpha (numeric): reference value for test
  flz <- list.files(folder)
  load(paste(folder, "/", flz[1], sep = ""))
  lf <- length(flz)
  
  # columns representing instantaneous effects
  inst.col <- dimnames(simulation$res)[[2]][which(grepl("\\.0", dimnames(simulation$res)[[2]]))]
  # columns representing lagged effects
  lag.col <- dimnames(simulation$res)[[2]][which(grepl("\\.1", dimnames(simulation$res)[[2]]))]
  
  p <- length(inst.col) # number of covariates
  nsim <- dim(simulation$res)[3] # number of simulations
  
  if (any(sapply(flz, function(str) grepl("setup", str)))) {
    load(paste(folder, "/", flz[lf], sep = ""))
    # read of causal connections
    As <- setup$As
    B1s <- setup$B1s
    lf <- lf - 1
    hidden <- setup$hidden
    print(hidden[1])
    hidden_lagged <- setup$hidden_lagged
  }
  
  # encode as binary
  As[abs(As) > 1e-5] <- 1
  As[abs(As) < 1e-5] <- 0
  
  B1s[abs(B1s) > 1e-5] <- 1
  B1s[abs(B1s) < 1e-5] <- 0
  
  M1s <- B1s # lagged effect starting with instantaneous edge
  Btot <- B1s # every way of lagged effect
  for(i in 1:(dim(B1s)[3])){
    M1s[,,i] <- As[,,i] %*% B1s[,,i]
    Btot[,,i] <- M1s[,,i] %*% As[,,i]
  }
  
  # no self-effects
  for (l in 1:p){
    As[l , l, ] <- NA
  }
  
  # combine instantaneous and lagged ancestors
  all.anc <- pmax(As, Btot, na.rm = TRUE) > 0
  dimnames(all.anc)[[1]] <- dimnames(all.anc)[[2]] <- dimnames(simulation$res)[[1]]
  # recursively get all
  all.anc[] <- apply(all.anc, 3, p.to.anc)
  non.anc <- !all.anc
  non.anc[all.anc] <- NA
  all.anc[!all.anc] <- NA
  
  for (j in 1:p){
    # self-effects not considered
    all.anc[j, j,] <- non.anc[j, j,] <- NA
  }
  
  # consider only instantaneous
  inst.anc <- As
  inst.anc[!inst.anc] <- NA
  non.inst.anc <- !As
  non.inst.anc[!non.inst.anc] <- NA
  
  if(hidden[1, 1] <= p){
    # to and from hidden not considered
    for(i in 1:nrow(hidden)){
      idx <- hidden[i, ]
      all.anc[idx, , i] <- NA
      all.anc[, idx, i] <- NA
      
      non.anc[idx, , i] <- NA
      non.anc[, idx, i] <- NA
      
      inst.anc[idx, , i] <- NA
      inst.anc[, idx, i] <- NA
      
      non.inst.anc[idx, , i] <- NA
      non.inst.anc[, idx, i] <- NA      
    }
  }
  
  # list for the different analysis
  TARs <- list()
  alpha.inds <- list()
  LINGAM_perf <- list()
  LINGAM_perf2 <- list()
  LINGAM_bl_perf <- list()
  lin.alpha.inds <- list()
  lin2.alpha.inds <- list()
  
  for (s in 1:2){
    if(s == 1){
      cat("Analysing instantaneous effects")
    } else {
      cat("Analysing summary graphs")
    }
    TAR <- matrix(NA, nsim + 2, 2 * lf) # matrix to store ROC
    TAR_lingam <- matrix(NA, nsim + 2, 2 * lf) # matrix to store ROC
    TAR_lingam2 <- matrix(NA, nsim + 2, 2 * lf) # matrix to store ROC
    TAR_bl <- matrix(NA, lf, 2) # matrix to store ROC
    alpha.ind <- integer(lf) # vector to find alpha performance
    lin.alpha.ind <- integer(lf) # vector to find alpha performance
    lin2.alpha.ind <- integer(lf) # vector to find alpha performance
    i <- 0
    for (file in flz[1:lf]) {
      # loop over files
      i <- i + 1
      load(paste(folder, "/", file, sep = ""))
      cat("\n", "T: ", simulation$n, "\n")
      z <- simulation$res # z-statistics
      b <- simulation$lingam # lingam results
      b2.pv <- simulation$b2_res
      dimnames(b2.pv) <- list(paste0('x', 1:p), paste0('x', 1:p, '.', rep(0:1, each = p)))
      
      pv <- 2 * pnorm(-abs(z)) # p-values
      b.pv <- array(0, dim = dim(simulation$b_res), 
                    dimnames = list(paste0('x', 1:p), paste0('x', 1:p)))
      b.pv[simulation$b_res == 0] <- 1
      
      
      if (s == 1){
        inst.pv <- pv[,inst.col,] # instantaneous effects
        dimnames(inst.pv)[[2]] <- dimnames(inst.pv)[[1]]
        for (j in 1:p){
          inst.pv[j, j,] <- 1
        }
        # multiplicity correction
        pv.adj <- inst.pv
        pv.adj[-hidden, -hidden_lagged, ] <- apply(inst.pv[-hidden, -hidden_lagged, ], 
                                                   3, function(pv) holm.corr(pv, cut = TRUE))
        
        # lowest p-value for null
        p.min <- pmin(apply(pv.adj * non.inst.anc, 3, min, na.rm = TRUE))
        lims <- sort(unique(c(0, alpha, p.min)))

        alpha.ind[i] <- which(lims == alpha)
        
        # find output structures at different alphas
        stru <- stru.anc <- stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], 
                                                             length(lims), nsim))
        stru[] <- apply(pv.adj, 3, find.instant.structures, lims = lims)
        
        for (k in 1:length(lims)){
          stru.anc[,,k,] <- stru[,,k,] * inst.anc # non-ancestors not considered
          stru.nonanc[,,k,] <- stru[,,k,] * non.inst.anc # ancestors not considered
        }
        
        
        lin_pv.adj <- array(1, dim(pv.adj), 
                            dimnames = list(paste0('x', 1:p), paste0('x', 1:p)))

        for(run in 1:nsim){
          S0 <- b[[run]][['S0']]
          S0.boot <- b[[run]][['S0.boot']]          
          
          n.boot <- dim(S0.boot)[3]
            
          lin.pv <- matrix(rowMeans(sapply(1:n.boot, function(boot){
            S0.boot[, , boot] >= S0
          })), ncol = ncol(S0))
          lin_pv.adj[-hidden[run, ], -hidden_lagged[run, ], run] <- lin.pv
          
        }
        for (j in 1:p){
          lin_pv.adj[j, j,] <- 1
        }
        
        lin_pv.adj[-hidden, -hidden_lagged, ] <- apply(lin_pv.adj[-hidden, -hidden_lagged, ], 
                                                   3, function(pv) holm.corr(pv, cut = TRUE))
        # lowest p-value for null
        p.min <- pmin(apply(lin_pv.adj * non.inst.anc, 3, min, na.rm = TRUE))
        lev <- sort(unique(c(0, alpha, p.min)))
        
        lin.alpha.ind[i] <- which(lev == alpha)
        
        b_stru <- b_stru.anc <- b_stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], 
                                                                 length(lev), nsim))
        b_stru[] <- apply(lin_pv.adj, 3, find.instant.structures, lims = lev)
        for (k in 1:length(lev)){
          b_stru.anc[,,k,] <- b_stru[,,k,] * inst.anc # non-ancestors not considered
          b_stru.nonanc[,,k,] <- b_stru[,,k,] * non.inst.anc # ancestors not considered
        }
        
        bl_stru <- bl_stru.anc <- bl_stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], 
                                                                   1, nsim))
        bl_stru[] <- apply(b.pv, 3, find.instant.structures, lims = 0.5)
        bl_stru.anc[,,1,] <- bl_stru[,,1,] * inst.anc # non-ancestors not considered
        bl_stru.nonanc[,,1,] <- bl_stru[,,1,] * non.inst.anc # ancestors not considered

        b2.inst.pv <- b2.pv[,inst.col,] # instantaneous effects
        dimnames(b2.inst.pv)[[2]] <- dimnames(b2.inst.pv)[[1]]
        for (j in 1:p){
          b2.inst.pv[j, j,] <- 1
        }
        # multiplicity correction
        b2.pv.adj <- b2.inst.pv
        b2.pv.adj[-hidden, -hidden_lagged, ] <- apply(b2.inst.pv[-hidden, -hidden_lagged, ], 
                                                   3, function(pv) holm.corr(pv, cut = TRUE))
        
        # lowest p-value for null
        p.min <- pmin(apply(b2.pv.adj * non.inst.anc, 3, min, na.rm = TRUE))
        b2.lims <- sort(unique(c(0, alpha, p.min)))
        
        lin2.alpha.ind[i] <- which(b2.lims == alpha)
        
        # find output structures at different alphas
        b2_stru <- b2_stru.anc <- b2_stru.nonanc <- array(NA, dim = c(dim(b2.pv.adj)[1:2], 
                                                             length(b2.lims), nsim))
        b2_stru[] <- apply(b2.pv.adj, 3, find.instant.structures, lims = b2.lims)
        
        for (k in 1:length(b2.lims)){
          b2_stru.anc[,,k,] <- b2_stru[,,k,] * inst.anc # non-ancestors not considered
          b2_stru.nonanc[,,k,] <- b2_stru[,,k,] * non.inst.anc # ancestors not considered
        }
        
        # power and FWER
        pwr <- apply(stru.anc, 3, mean, na.rm = TRUE)
        FWER <- apply(apply(stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        b_pwr <- apply(b_stru.anc, 3, mean, na.rm = TRUE)
        b_FWER <- apply(apply(b_stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        bl_pwr <- apply(bl_stru.anc, 3, mean, na.rm = TRUE)
        bl_FWER <- apply(apply(bl_stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        b2_pwr <- apply(b2_stru.anc, 3, mean, na.rm = TRUE)
        b2_FWER <- apply(apply(b2_stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        TAR_bl[i, ] <- c(bl_pwr, bl_FWER)
        TAR_lingam[1:length(lev),c(i, lf + i)] <- c(b_pwr, b_FWER)
        TAR_lingam2[1:length(b2.lims),c(i, lf + i)] <- c(b2_pwr, b2_FWER)
        TAR[1:length(lims),c(i, lf + i)] <- c(FWER, pwr)
      } else if (s == 2){
        # summary effects
        b2.sum.pv <- b2.pv[,inst.col,]
        dimnames(b2.sum.pv)[[2]] <- dimnames(b2.sum.pv)[[1]]

        # get summary p-values
        b2.sum.pv[] <- apply(b2.pv, 3, summary.p.val)
        # multiplicity correction
        b2.pv.adj <- b2.sum.pv
        b2.pv.adj[-hidden, -hidden_lagged, ] <- apply(b2.sum.pv[-hidden, -hidden_lagged, ], 
                                                   3, function(pv) holm.corr(pv, cut = TRUE))
        
        # lowest p-value for null        
        p.min <- pmin(apply(b2.pv.adj * non.anc, 3, min, na.rm = TRUE))
        b2.lims <- sort(unique(c(0, alpha, p.min)))
        b2.lims <- b2.lims[is.finite(b2.lims)]

        lin2.alpha.ind[i] <- which(b2.lims == alpha)
        
        # find output structures at different alphas
        b2_stru <- b2_stru.anc <- b2_stru.nonanc <- array(NA, dim = c(dim(b2.pv.adj)[1:2], 
                                                             length(b2.lims), nsim))
        b2_stru[] <- apply(b2.pv.adj, 3, find.structures, lims = b2.lims)

        for (k in 1:length(b2.lims)){
          b2_stru.anc[,,k,] <- b2_stru[,,k,] * all.anc # non-ancestors not considered
          b2_stru.nonanc[,,k,] <- b2_stru[,,k,] * non.anc # ancestors not considered
        }
        
        #lin2
        # summary effects
        sum.pv <- pv[,inst.col,]
        dimnames(sum.pv)[[2]] <- dimnames(sum.pv)[[1]]
        
        # get summary p-values
        sum.pv[] <- apply(pv, 3, summary.p.val)
        # multiplicity correction
        pv.adj <- sum.pv
        pv.adj[-hidden, -hidden_lagged, ] <- apply(sum.pv[-hidden, -hidden_lagged, ], 
                                                   3, function(pv) holm.corr(pv, cut = TRUE))
        
        # lowest p-value for null        
        p.min <- pmin(apply(pv.adj * non.anc, 3, min, na.rm = TRUE))
        lims <- sort(unique(c(0, alpha, p.min)))
        lims <- lims[is.finite(lims)]
        
        alpha.ind[i] <- which(lims == alpha)
        
        # find output structures at different alphas
        stru <- stru.anc <- stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], 
                                                             length(lims), nsim))
        stru[] <- apply(pv.adj, 3, find.structures, lims = lims)
        
        for (k in 1:length(lims)){
          stru.anc[,,k,] <- stru[,,k,] * all.anc # non-ancestors not considered
          stru.nonanc[,,k,] <- stru[,,k,] * non.anc # ancestors not considered
        }
        
        
        lin_pv.adj <- array(1, dim(pv.adj), 
                            dimnames = list(paste0('x', 1:p), paste0('x', 1:p)))
        
        for(run in 1:nsim){
          Slag <- b[[run]][['Slag']]
          Slag.boot <- b[[run]][['Slag.boot']]          
          
          n.boot <- dim(S0.boot)[3]
          
          lin.pv <- matrix(rowMeans(sapply(1:n.boot, function(boot){
            Slag.boot[, , boot] >= Slag
          })), ncol = ncol(Slag))
          lin_pv.adj[-hidden[run, ], -hidden_lagged[run, ], run] <- lin.pv
        }
        for (j in 1:p){
          lin_pv.adj[j, j,] <- 1
        }
        
        lin_pv.adj[-hidden, -hidden_lagged, ] <- apply(lin_pv.adj[-hidden, -hidden_lagged, ], 
                                                       3, function(pv) holm.corr(pv, cut = TRUE))
        # lowest p-value for null
        p.min <- pmin(apply(lin_pv.adj * non.inst.anc, 3, min, na.rm = TRUE))
        lev <- sort(unique(c(0, alpha, p.min)))
        
        lin.alpha.ind[i] <- which(lev == alpha)
        
        b_stru <- b_stru.anc <- b_stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], 
                                                                   length(lev), nsim))
        b_stru[] <- apply(lin_pv.adj, 3, find.structures, lims = lev)
        for (k in 1:length(lev)){
          b_stru.anc[,,k,] <- b_stru[,,k,] * all.anc
          b_stru.nonanc[,,k,] <- b_stru[,,k,] * non.anc
        }
        
        bl_stru <- bl_stru.anc <- bl_stru.nonanc <- array(NA, dim = c(dim(pv.adj)[1:2], 
                                                                   1, nsim))
        bl_stru[] <- apply(b.pv, 3, find.structures, lims = 0.5)
        bl_stru.anc[,,1,] <- bl_stru[,,1,] * all.anc
        bl_stru.nonanc[,,1,] <- bl_stru[,,1,] * non.anc
        
        # power and FWER
        pwr <- apply(stru.anc, 3, mean, na.rm = TRUE)
        FWER <- apply(apply(stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)

        b_pwr <- apply(b_stru.anc, 3, mean, na.rm = TRUE)
        b_FWER <- apply(apply(b_stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        bl_pwr <- apply(bl_stru.anc, 3, mean, na.rm = TRUE)
        bl_FWER <- apply(apply(bl_stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        b2_pwr <- apply(b2_stru.anc, 3, mean, na.rm = TRUE)
        b2_FWER <- apply(apply(b2_stru.nonanc, 3:4, max, na.rm = TRUE) == 1, 1, mean)
        
        TAR_bl[i, ] <- c(bl_pwr, bl_FWER)
        TAR_lingam[1:length(lev),c(i, lf + i)] <- c(b_pwr, b_FWER)
        TAR_lingam2[1:length(b2.lims),c(i, lf + i)] <- c(b2_pwr, b2_FWER)
        TAR[1:length(lims),c(i, lf + i)] <- c(FWER, pwr)
      } else {
        stop("Wrong type")
      }
    }
    # save for convenience
    LINGAM_perf[[s]] <- TAR_lingam
    LINGAM_perf2[[s]] <- TAR_lingam2
    LINGAM_bl_perf[[s]] <- TAR_bl
    TARs[[s]] <- TAR
    alpha.inds[[s]] <- alpha.ind
    lin.alpha.inds[[s]] <- lin.alpha.ind
    lin2.alpha.inds[[s]] <- lin2.alpha.ind
  }
  
  par(mfrow = c(3,2))
  for (s in 1:2){
    # read off from lists
    TAR <- TARs[[s]]
    alpha.ind <- alpha.inds[[s]]
    # plot ROC
    matplot(TAR[-1,1:lf], TAR[-1,lf + (1:lf)], type = "s",
            xlim = c(0, max(c(TAR[,1:lf]), na.rm = TRUE)), ylim = c(0,1), xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
            col = (1:p)[-5], las = 1, main = "Ancestor")
    # add performance at alpha
    points(diag(TAR[alpha.ind,1:lf]), diag(TAR[alpha.ind,lf + (1:lf)]),
           col = (1:p)[-5], pch = 3)
    # plot target alpha
    lines(c(alpha, alpha), c(0, 1), col = "gray", lty = 2)
    if(s == 1){
      # add performance of simple pruned lingam
      points(x = LINGAM_bl_perf[[s]][, 2], y = LINGAM_bl_perf[[s]][, 1],
             col = (1:p)[-5], pch = 2)
    }
  }
  for (s in 1:2){
    # read off from lists
    TAR <- LINGAM_perf[[s]]
    alpha.ind <- lin.alpha.inds[[s]]
    # plot ROC
    matplot(TAR[-1,1:lf], TAR[-1,lf + (1:lf)], type = "s",
            xlim = c(0, max(c(TAR[,1:lf]), na.rm = TRUE)), ylim = c(0,1), xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
            col = (1:p)[-5], las = 1, main = 'LINGAM')
    # add performance at alpha
    points(diag(TAR[alpha.ind,1:lf]), diag(TAR[alpha.ind,lf + (1:lf)]),
           col = (1:p)[-5], pch = 3) 
    # plot target alpha
    lines(c(alpha, alpha), c(0, 1), col = "gray", lty = 2)
    
    if(s == 1){
    # add performance of simple pruned lingam
    points(x = LINGAM_bl_perf[[s]][, 2], y = LINGAM_bl_perf[[s]][, 1],
           col = (1:p)[-5], pch = 2)
    }
  }
  for (s in 1:2){
    # read off from lists
    TAR <- LINGAM_perf2[[s]]
    alpha.ind <- lin2.alpha.inds[[s]]
    # plot ROC
    matplot(TAR[-1,1:lf], TAR[-1,lf + (1:lf)], type = "s",
            xlim = c(0, max(c(TAR[,1:lf]), na.rm = TRUE)), ylim = c(0,1), xlab = "Type I FWER", ylab ="Fraction of detected ancestors",
            col = (1:p)[-5], las = 1, main = 'LINGAM')
    # add performance at alpha
    points(diag(TAR[alpha.ind,1:lf]), diag(TAR[alpha.ind,lf + (1:lf)]),
           col = (1:p)[-5], pch = 3) 
    # plot target alpha
    lines(c(alpha, alpha), c(0, 1), col = "gray", lty = 2)
    
    if(s == 1){
      # add performance of simple pruned lingam
      points(x = LINGAM_bl_perf[[s]][, 2], y = LINGAM_bl_perf[[s]][, 1],
            col = (1:p)[-5], pch = 2)
    }
  }
}

# png(paste(savefolder, "/ROC-graph-noleg.png", sep = ""), width = 600 * plotfac,
# height = 300 * plotfac, res = 75 * plotfac)
# network(folder)
# dev.off()


network('results/test')

