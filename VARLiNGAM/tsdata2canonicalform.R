tsdata2canonicalform <- function(X,nlags=1, verbose = FALSE) {

  # INPUT
  # X: timeseries data matrix with:
  #    no missing values
  #    more observations than variables
  #    observations for every time-step given
  #    contains only the observations of the varibales, no time indices,
  #      variable names or similar
  # nlags: number of time lags used, nlags >= 1
  #
  # OUTPUT
  # X in canonical form
  #
  ## canonical format - columns:
  # 1st col: 'id number' for panel data, each panel gets an id number
  #           here, this will just be a columns of ones
  # 2nd col: 'time-point t of current value, starting from 1'
  # following: 'past values', nvar*nlags colums:
  #            nvar variables at time t-1
  #            nvar variables at time t-2
  #            nvar variables at time t-nlags
  # following: 'current values', nvar columns:
  #            nvar variables at time t
  # where nvar = number of variables

  # get variables in columns, observations in rows
  dims <- dim(X)
  if (dims[1]>dims[2]) {
    nobs <- dims[1] # number of observations
    nvar <- dims[2] # number of variables
  }
  else {
    X <- t(X)
    nobs <- dims[2] # number of observations
    nvar <- dims[1] # number of variables
  }

  p <- ncol(X)
  xt <- matrix(NA, nrow = nrow(X), ncol = p * (nlags +1))
  for(j in 1:p){
    xtj <- lagmatrix(X[,j], c(1:nlags, 0))

    xt[,j + p * (0:nlags)] <- xtj
  }
  xt <- xt[complete.cases(xt),]
  xt <- cbind(rep(1, nrow(xt)), 1:nrow(xt), xt)

  # get labels for colnames
  # curvali = i-the variable (at current time t)
  # pastvalij = i-th variable at time t-j
  temp1 <- array(0,dim=c(1,nvar));
  temp2 <- NULL
  for (i in 1:nlags) {
    temp1 <- paste(sep="",'pastval', 1:nvar, i)
    temp2 <- c(temp2,temp1)
  }
  temp1 <- paste(sep="",'curval', 1:nvar)
  temp2 <- c(temp2,temp1)
  
  colnames(xt) <- c("id","curtime",temp2)

  xt

}
