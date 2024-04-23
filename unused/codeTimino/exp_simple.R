# Copyright (c) 2010-2012  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
set.seed(1)
n <- 1000
w <- rep(0,n)
x <- rep(0,n)
y <- rep(0,n)
epsw <- rnorm(n)^3
epsx <- rnorm(n)^3
epsy <- rnorm(n)^3

for(i in 3:n)
{
    x[i] <- 0.3*x[i-1]+0.5*epsx[i]
    y[i] <- 0.8*y[i-1]+0.8*x[i-1]+0.5*epsy[i]
    w[i] <- -0.6*w[i-1]+0.8*x[i-2]+0.5*epsw[i]
}
d <- timino_dag(cbind(x,y,w), alpha = 0.05, max_lag = 2, model = traints_linear, indtest = indtestts_crosscov, output = TRUE)

show("====")
show("DONE")
show("====")

show("true summary time graph:")
show(cbind(c(0,0,0),c(1,0,0),c(1,0,0)))

show("estimated summary time graph:")
show(d)

