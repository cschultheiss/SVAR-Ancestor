require(tsutils)
lin.anc.ts <- function(x, degree, targets = colnames(x),  f  = function(x) x^3){
  cols <- colnames(x)
  ind <- which(cols %in% targets)
  p <- ncol(x)
  xt <- matrix(NA, nrow = nrow(x), ncol = p * (degree +1))
  for(j in 1:p){
    xtj <- lagmatrix(x[,j], 0:degree)
    
    xt[,1:(degree + 1) + (j -1) * (degree +1)] <- xtj
  }
  xt <- xt[complete.cases(xt),]
  colnames(xt) <- paste(rep(cols, each = degree + 1), ".", 0:degree, sep ="")
  z.val <- matrix(NA, nrow = length(targets), ncol = p * (degree + 1), dimnames = list(targets, colnames(xt)))
  targets <- paste(targets, ".", 0, sep ="")
  n <- nrow(xt)
  
  for (i in 1:length(targets)){
    tari <- f(xt[,targets[i]])
    for(s in 0:degree){
      wico <- seq(2, (p - 1)*(degree + 1) + 2, degree + 1)
      wiz <- seq(1, (p - 1)*(degree + 1) + 1, degree + 1) + s
      su <- summary(lm(tari[(s + 1):n] ~xt[1: (n - s), ]))
      z.val[i, wiz] <- su$coefficients[wico,3]
    }
  }
  p.val <- pnorm(abs(z.val), lower.tail = FALSE)*2
  return(list(z.val = z.val, p.val = p.val))
}