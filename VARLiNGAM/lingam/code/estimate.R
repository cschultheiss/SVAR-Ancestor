# ESTIMATE Estimate LiNGAM model from data.
#
# Estimates a LiNGAM model from data. The returned weight matrix B
# is strictly lower triangular, but the weights are not pruned. To
# prune the matrix, call PRUNE.
#
# SYNTAX:
# res <- estimate(X)
#
# INPUT:
# X     - Data matrix: each row is an observed variable, each
#         column one observed sample. The number of columns should
#         be far larger than the number of rows.
#
# OUTPUT:
# B     - Matrix of estimated connection strenghts
# stde  - Standard deviations of disturbance variables
# ci    - Constants
# k     - An estimated causal ordering
# W     - The demixing matrix of ICs, in estimated row ordering.
#
# (Note that B, stde, ci, and optionally W are all ordered
#  according to the same variable ordering as X.)
#
# See also PRUNE, LINGAM.

estimate <- function( X, verbose = FALSE ) {

  # Using the fastICA R package so make sure it is loaded
  library('fastICA')
  
  # Call the fastICA algorithm
  icares <- fastICA( t(X), nrow(X), tol=1e-14 )
  W <- t((icares$K) %*% (icares$W)) # un-mixing matrix
  #A <- solve(W) # mixing matrix
  ICs <- icares$S
  
  if(verbose) print(W)

  if(verbose) print('Test for Gaussianity of the independent components')
  if(verbose) print(Gauss_Tests(ICs))
  
  # [Here, we really should perform some tests to see if the 'icasig' 
  # really are independent. If they are not very independent, we should
  # issue a warning. This is not yet implemented.]

  # Try to permute the rows of W so that sum(1./abs(diag(W))) is minimized
  if(verbose) cat('Performing row permutation...\n');
  dims <- nrow(X)
  if (dims <= 8) {  
    if(verbose) cat('(Small dimensionality, using brute-force method.)\n')
    temp <- nzdiagbruteforce( W )
    Wp <- temp$Wopt
    rowp <- temp$rowp
  }
  else {
    if(verbose) cat('(Using the Hungarian algorithm.)\n')
    stop('Not implemented yet!')
  }
  if(verbose) cat('Done!\n')

  # Divide each row of Wp by the diagonal element
  estdisturbancestd <- 1/diag(abs(Wp))
  Wp <- Wp/diag(Wp)

  # Compute corresponding B
  Best <- diag(dims)-Wp

#  if(verbose) print(Best)

  # Estimate the constants c
  m <- rowMeans(X)
  dim(m) <- c(dims,1)
  cest <- Wp %*% m

  # Next, identically permute the rows and columns of B so as to get an
  # approximately strictly lower triangular matrix
  if(verbose) cat('Performing permutation for causal order...\n');

  if (dims <= 8) {  
    if(verbose) cat('(Small dimensionality, using brute-force method.)\n');
    temp <- sltbruteforce( Best )
    Bestcausal <- temp$Bopt
    causalperm <- temp$optperm
  }
  else {
    if(verbose) cat('(Using pruning algorithm.)\n')
    stop('Not implemented yet!') 
  }
  if(verbose) cat('Done!\n');

  if(verbose) print(Bestcausal)
  
  # Here, we report how lower triangular the result was, and in 
  # particular we issue a warning if it was not so good!
  percentinupper <- sltscore(Bestcausal)/sum(Bestcausal^2)
  if (percentinupper>0.2) 
    if(verbose) cat('WARNING: Causal B not really triangular at all!!\n')
  else if (percentinupper>0.05)
    if(verbose) cat('WARNING: Causal B only somewhat triangular!\n')
  else
    if(verbose) cat('Causal B nicely triangular. No problems to report here.\n')

  # Set the upper triangular to zero
  Bestcausal[upper.tri(Bestcausal,diag=FALSE)] <- 0

  # Finally, permute 'Bestcausal' back to the original variable
  # ordering and rename all the variables to the way we defined them
  # in the function definition
  icausal <- iperm(causalperm);
  res <- list()
  res$B <- Bestcausal[icausal, icausal];
  res$stde <- estdisturbancestd
  res$ci <- cest
  res$k <- causalperm
  res$W <- W
  res$ICs <- ICs
  res$pinup <- percentinupper # added for firmgrowth conintegration analysis

  # Return the result
  res

}
