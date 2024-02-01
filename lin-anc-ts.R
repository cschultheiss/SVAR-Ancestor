require(tsutils)
lin.anc.ts <- function(x, degree, targets = colnames(x),  f  = function(x) x^3){
  cols <- colnames(x)
  ind <- which(cols %in% targets)
  p <- ncol(x)
  xt <- matrix(NA, nrow = nrow(x), ncol = p * (degree +1))
  for(j in 1:p){
    xtj <- lagmatrix(x[,j], 0:degree)
    
    xt[,j + p * (0:degree)] <- xtj
  }
  xt <- xt[complete.cases(xt),]
  colnames(xt) <- paste(rep(cols, degree + 1), ".", rep(0:degree, each = p), sep ="")

  if(degree > 0){
    u <- xt[,1:p] - xt[,-c(1:p)] %*% solve(crossprod(xt[,-c(1:p)])) %*% crossprod(xt[,-c(1:p)], xt[,1:p])
  } else {
    u <- xt
  }
  n <- nrow(u)
  z.val <- matrix(NA, nrow = length(targets), ncol = p * (degree + 1), dimnames = list(targets, colnames(xt)))
  
  for (s in 0:degree){
    if(degree > 0){
      us <- xt[(s+1):n,1:p] - xt[1:(n-s),-c(1:p)] %*% solve(crossprod(xt[1:(n-s),-c(1:p)]))  %*% 
        crossprod(xt[1:(n-s),-c(1:p)], xt[(s+1):n,1:p])
    } else {
      us <- xt
    }
    
    colnames(us) <- colnames(x)
    wiz <- 1:p + s*p
    
    for (i in 1:length(targets)){
      tari <- f(us[,targets[i]])
      su <- summary(lm(tari ~u[1: (n - s), ]))
      z.val[i, wiz] <- su$coefficients[-1,3]
    }
  }

  p.val <- pnorm(abs(z.val), lower.tail = FALSE)*2
  return(list(z.val = z.val, p.val = p.val))
}

summary.p.val<- function(lin.anc){
  if(is.list(lin.anc)) {
    pv <- lin.anc$p.val
  } else {
    pv <- lin.anc
  }
  
  preds <- colnames(pv)
  targets <- rownames(pv)
  p <- sum(grepl("\\.0", preds))
  npv <- ncol(pv)/p
  spv <- pv[,1:p]
  colnames(spv) <- gsub("\\.0.*","",colnames(spv))
  spv[] <- 1
  for(i in 1:length(targets)){
    for (j in 1:p){
      if(targets[i] == colnames(spv)[j]) next()
      pij <- pv[i, seq(j, j + p*(npv - 1), length.out = npv)]
      spv[i,j] <- min(min(sort(pij) * npv /(1:npv)) * sum(1/(1:npv)), 1)
    }
  }
  return(spv)
}

summary.graph <- function(lin.anc, alpha = 0.05) {
  spv <- summary.p.val(lin.anc)
  pv.corr <- holm.corr(spv)
  pmat <- pv.corr < alpha
  return(p.to.anc(pmat))
}

p.to.anc <- function(pmat) {
  preds <- colnames(pmat)
  targets <- rownames(pmat)
  ancmat <- pmat
  for (i in 1:length(targets)){
    tested <- setdiff(preds, targets)
    an <- names(which(pmat[i,]))
    while (length(setdiff(an, tested)) > 0) {
      for (k in setdiff(an, tested)) {
        j <- which(targets == k)
        an <- unique(c(an, names(which(pmat[j,]))))
        tested <- c(tested, k)
      }
    }
    ancmat[i, which(preds %in% an)] <- TRUE
  }
  return(ancmat)
}

instant.graph <- function(lin.anc, alpha = 0.05, verbose = FALSE, corr = TRUE){
  if(is.list(lin.anc)) {
    pv <- lin.anc$p.val
  } else {
    pv <- lin.anc
  }
  preds <- colnames(pv)
  targets <- rownames(pv)
  if (corr){
    p <- sum(grepl("\\.0", preds))
    pv <- pv[,1:p]
    colnames(pv) <- gsub("\\.0.*","",colnames(pv))
    for(ta in targets) pv[ta, ta] <- 1
    pv.corr <- holm.corr(pv)
  } else {
    pv.corr <- pv
  }
  anc1 <- pv.corr < alpha
  anc <- p.to.anc(anc1)

  if(sum(sapply(targets, function(ta) anc[ta , ta])) == 0){
    pv.corr[] <- anc
    return(list(rec.ancs = pv.corr, alpha = alpha))
  } else {
    loop.vars <- loop.vars <- names(which(sapply(targets, function(ta) anc[ta , ta])))
    pvs.mat <- pv.corr[loop.vars, loop.vars]
    pvs.sub <- pvs.mat[pvs.mat < alpha]
    new.alpha <- max(pvs.sub)
    if(verbose) print(paste("Try decreasing alpha from", round(alpha, 3), "to", round(new.alpha, 3)))
    
    out <- instant.graph(pvs.mat, new.alpha, verbose = verbose, corr = FALSE)
    anc1[loop.vars, loop.vars] <- out$rec.ancs == 1
    return(list(rec.ancs = p.to.anc(anc1), alpha = out$alpha))
  }
}

holm.corr <- function(pv, cut = TRUE){
  np <- length(pv)
  ntest <- np - sum(rownames(pv) %in% colnames(pv))
  i <- c(seq_len(ntest), rep(ntest, np - ntest))
  o <- order(pv)
  ro <- order(o)
  if (cut) {
    pv[] <- pmin(1, cummax((ntest + 1L - i) * pv[o]))[ro]
  }
  else {
    pv[] <- cummax((ntest + 1L - i) * pv[o])[ro]
  }
  pv
}

holm.uncut <- function(pv){
  p <- length(pv)
  i <- seq_len(p)
  o <- order(pv)
  ro <- order(o)
  cummax((p + 1L - i) * pv[o])[ro]
}