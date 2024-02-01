source("helpers.R")
holm.uncut <- function(pv){
  p <- length(pv)
  i <- seq_len(p)
  o <- order(pv)
  ro <- order(o)
  cummax((p + 1L - i) * pv[o])[ro]
}