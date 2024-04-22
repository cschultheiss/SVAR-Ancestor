source("helpers.R")

find.structures <- function(pvs, lims){
  # function to find ancestral structures at different significance levels, it allows for cycles
  # Input
  # pvs (numeric, matrix): p-values
  # lims (numeric, vector): significance levels to be considered, sorted from low to high
  # Output
  # boolean array: constructed ancestral structure at each level
  
  su <- sapply(lims, function(lim) sum(c(pvs) < lim)) # number of detected connections at each level
  # find jumps, only consider jump points
  di <- diff(su)
  wi <- c(0, which(di > 0)) + 1
  nl <- length(lims) # number of significance levels considered
  nw <- length(wi) # number of jumps
  wi <- c(wi, nl) # include the last one to store
  out <- array(NA, dim = c(dim(pvs), nl)) # array to store structures
  dimnames(out)[1:2] <- dimnames(pvs)
  dimnames(out)[[3]] <- lims
  for (i in 1:nw){
    # loop over jump points
    out[,,wi[i] : wi[i+1]] <- p.to.anc(pvs < lims[wi[i]]) # store structure to all entries from jump to jump
  }
  out
}

find.instant.structures <- function(pvs, lims, verbose.final = FALSE, verbose = FALSE){
  # function to find ancestral structures at different significance levels, it does not allow for cycles
  # Input
  # pvs (numeric, matrix): p-values
  # lims (numeric, vector): significance levels to be considered, sorted from low to high
  # verbose.final (boolean): print significance level needed in the final step, if 1 is used as level, this a model p-value
  # verbose (boolean): print statements from (instant.graph)
  # Output
  # boolean array: constructed ancestral structure at each level
  
  su <- sapply(lims, function(lim) sum(c(pvs) < lim)) # number of detected connections at each level
  # find jumps, only consider jump points
  di <- diff(su)
  wi <- c(0, which(di > 0)) + 1
  nl <- length(lims) # number of significance levels considered
  nw <- length(wi) # number of jumps
  wi <- c(wi, nl) # include the last one to store
  out <- array(NA, dim = c(dim(pvs), nl)) # array to store structures
  dimnames(out)[1:2] <- dimnames(pvs)
  dimnames(out)[[3]] <- lims
  for (i in 1:nw){
    # loop over jump points
    # find structure for given level
    stru <- instant.graph(pvs, alpha = lims[wi[i]], verbose = verbose, corr = FALSE)
    out[,,wi[i] : wi[i+1]] <- stru[[1]] # store structure to all entries from jump to jump
  }
  if(verbose.final) print(paste("Used alpha = ", round(stru[[2]], 3), sep = ""))
  out
}