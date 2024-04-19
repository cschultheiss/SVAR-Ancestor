source("helpers.R")

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

find.instant.structures <- function(pvs, lims, verbose.final = FALSE, verbose = FALSE){
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
    stru <- instant.graph(pvs, alpha = lims[wi[i]], verbose = verbose, corr = FALSE)
    out[,,wi[i] : wi[i+1]] <- stru[[1]]
  }
  if(verbose.final) print(paste("Used alpha = ", round(stru[[2]], 3), sep = ""))
  out
}