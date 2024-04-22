summary.p.val<- function(lin.anc){
  # function to provide summary p-values
  # Input
  # lin.anc: full output from lin.anc.ts() or only p-values
  # Ouptut
  # summary p-values, numeric matrix 
  
  # check input format
  if(is.list(lin.anc)) {
    pv <- lin.anc$p.val
  } else {
    pv <- lin.anc
  }
  
  # get variables names
  preds <- colnames(pv)
  targets <- rownames(pv)
  p <- sum(grepl("\\.0", preds)) # number of predictors per time step
  npv <- ncol(pv)/p # number of time steps
  spv <- pv[,1:p] # matrix shape for summary p-value
  colnames(spv) <- gsub("\\.0.*","",colnames(spv)) # no time indicator
  spv[] <- 1
  for(i in 1:length(targets)){
    # loop through targets
    for (j in 1:p){
      # loop through predictors
      if(targets[i] == colnames(spv)[j]) next() # omit self-effects
      # get p-values at all time steps for that predictor
      pij <- pv[i, seq(j, j + p*(npv - 1), length.out = npv)]
      # combine p-values
      spv[i,j] <- min(min(sort(pij) * npv /(1:npv)) * sum(1/(1:npv)), 1)
    }
  }
  return(spv)
}

summary.graph <- function(lin.anc, alpha = 0.05) {
  # function to get the summary graph
  # Input
  # lin.anc: full output from lin.anc.ts() or only p-values
  # alpha (numeric): significance level
  # Ouptut
  # boolean matrix: Indicator whether one variable affects another
  
  spv <- summary.p.val(lin.anc) # get summary p-values
  pv.corr <- holm.corr(spv) # multiplicity correction
  pmat <- pv.corr < alpha # convert to boolean
  return(p.to.anc(pmat)) # construct recursively
}

p.to.anc <- function(pmat) {
  # function to recursively construct all ancestral connections
  # Input
  # pmat (boolean, matrix): Indicator whether connection was detected
  # Ouput
  # boolean matrix: Indicator whether connection was detected or constructed
  
  # get variable names
  preds <- colnames(pmat)
  targets <- rownames(pmat)
  ancmat <- pmat # matrix shape for output
  for (i in 1:length(targets)){
    # loop through targets
    # which predictors have been checked, if for some no ancestors where search, do not use these
    tested <- setdiff(preds, targets)  
    an <- names(which(pmat[i,])) # names of detected ancestors
    while (length(setdiff(an, tested)) > 0) {
      # if not all ancestors have been checked
      for (k in setdiff(an, tested)) {
        # loop through unchecked
        j <- which(targets == k) # find position of ancestor
        an <- unique(c(an, names(which(pmat[j,])))) # add ancestors of ancestors
        tested <- c(tested, k) # declare it as checked
      }
    }
    ancmat[i, which(preds %in% an)] <- TRUE # add additional ancestors to output matrix
  }
  return(ancmat)
}

instant.graph <- function(lin.anc, alpha = 0.05, verbose = FALSE, corr = TRUE){
  # function to get the instant graph
  # Input
  # lin.anc: full output from lin.anc.ts() or only p-values
  # alpha (numeric): significance level
  # verbose (boolean): should information be printed?
  # corr (boolean): is multiplicity correction required?
  # Ouptut
  # boolean matrix:
  # rec.ancs (boolean, matrix): indicator whether one variable affects another instantaneously
  # alpha (numeric): significance level to avoid cycles
  
  # check input format
  if(is.list(lin.anc)) {
    pv <- lin.anc$p.val
  } else {
    pv <- lin.anc
  }
  
  # get variable names
  preds <- colnames(pv)
  targets <- rownames(pv)

  p <- sum(grepl("\\.0", preds)) # number of instant predictors
  if (p > 0){
    pv <- pv[,1:p] # instant p-values
    colnames(pv) <- gsub("\\.0.*","",colnames(pv)) # shorten column names
  }

  for(ta in targets) pv[ta, ta] <- 1 # no self-effects
  # apply Bonferroni-Holm if needed
  if (corr){
    pv.corr <- holm.corr(pv)
  } else {
    pv.corr <- pv
  }
  
  anc1 <- pv.corr < alpha # detected instantaneous ancestors
  anc <- p.to.anc(anc1) # reconstruct further
  if(sum(sapply(targets, function(ta) anc[ta , ta])) == 0){
    # if there are no cycles
    pv.corr[] <- anc # format for output matrix
    return(list(rec.ancs = pv.corr > 0, alpha = alpha))
  } else {
    loop.vars <- names(which(sapply(targets, function(ta) anc[ta , ta]))) # find variables that cause loops
    pvs.mat <- pv.corr[loop.vars, loop.vars] # get submatrix of p-values
    pvs.sub <- pvs.mat[pvs.mat < alpha] # get significant p-values
    new.alpha <- max(pvs.sub) # new significance level
    if(verbose) print(paste("Try decreasing alpha from", round(alpha, 3), "to", round(new.alpha, 3)))
    
    # create graph with submatrix and new significance level
    out <- instant.graph(pvs.mat, new.alpha, verbose = verbose, corr = FALSE) 
    anc1[loop.vars, loop.vars] <- out$rec.ancs == 1 # adapt detected ancestors
    # construct recursively again
    return(list(rec.ancs = p.to.anc(anc1), alpha = out$alpha))
  }
}

holm.corr <- function(pv, cut = TRUE){
  # function to perform Bonferroni-Holm correction on matrix with possibility to ignore some entries
  # Input
  # pv (numeric, matrix): p-values, unused p-values should be 1
  # cut (boolean): should correction cut p-values at 1
  # Output
  # numeric matrix: corrected p-values
  
  np <- length(pv) # total number of p-values
  ntest <- np - sum(rownames(pv) %in% colnames(pv)) # total number of tests, i.e., without self-effects
  i <- c(seq_len(ntest), rep(ntest, np - ntest)) # used for correction factor
  o <- order(pv) # ranking of p-values
  ro <- order(o) # ranking of rankings to sort back
  if (cut) {
    # correct by factor according to ordering, cut at 1, sort back
    pv[] <- pmin(1, cummax((ntest + 1L - i) * pv[o]))[ro]
  }
  else {
    # correct by factor according to ordering, sort back
    pv[] <- cummax((ntest + 1L - i) * pv[o])[ro]
  }
  pv
}