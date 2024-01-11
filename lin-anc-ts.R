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
  u <- xt[,1:p] - xt[,-c(1:p)] %*% solve(crossprod(xt[,-c(1:p)])) %*% crossprod(xt[,-c(1:p)], xt[,1:p])
  n <- nrow(u)
  z.val <- matrix(NA, nrow = length(targets), ncol = p * (degree + 1), dimnames = list(targets, colnames(xt)))
  
  for (s in 0:degree){
    us <- xt[(s+1):n,1:p] - xt[1:(n-s),-c(1:p)] %*% solve(crossprod(xt[1:(n-s),-c(1:p)]))  %*% 
      crossprod(xt[1:(n-s),-c(1:p)], xt[(s+1):n,1:p])
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