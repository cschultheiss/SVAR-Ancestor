source("helpers.R")
holm.uncut <- function(pv){
  p <- length(pv)
  i <- seq_len(p)
  o <- order(pv)
  ro <- order(o)
  cummax((p + 1L - i) * pv[o])[ro]
}

find.structures <- function(pvs, lims){
  su <- sapply(lims, function(lim) sum(c(pvs) < lim))
  di <- diff(su)
  wi <- c(0, which(di > 0)) + 1
  nl <- length(lims)
  nw <- length(wi)
  wi <- c(wi, nl)
  out <- array(NA, dim = c(dim(pvs), nl))
  dimnames(out)[1:2] <- dimnames(pvs)
  dimnames(out)[[3]] <- lims
  for (i in 1:nw){
    
    out[,,wi[i] : wi[i+1]] <- p.to.anc(pvs < lims[wi[i]])
  }
  out
}