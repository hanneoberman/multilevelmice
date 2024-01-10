# 3Rosen.R JN 220409
# To show three variants of an extended Rosenbrock function

fnR <- function (x) 
{ # adagio fnRosenbrock
  n <- length(x)
  x1 <- x[2:n]
  x2 <- x[1:(n - 1)]
  sum(100 * (x1 - x2^2)^2 + (1 - x2)^2)
}

genrose.f<- function(x, gs=NULL){ # From Nash/Walker-Smith 1987 
  ## One generalization of the Rosenbrock banana valley function (n parameters)
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
  return(fval)
}

fnexr = function(par) { # funconstrain::ex_rosen
  n <- length(par)
  if (n %% 2 != 0) {
    stop("Extended Rosenbrock: n must be even")
  }
  
  fsum <- 0
  for (i in 1:(n / 2)) {
    p2 <- 2 * i
    p1 <- p2 - 1
    
    f_p1 <- 10 * (par[p2] - par[p1] ^ 2)
    f_p2 <- 1 - par[p1]
    fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
  }
  
  fsum
}

xx<-rep(1.1,6)
cat("  fnR      fnexr(xx)    genrose.f\n")
for (i in 1:6) {
  xx[i] <- .2
  cat(fnR(xx)," ",fnexr(xx)," ",genrose.f(xx)," at: ")
  print(xx)
  xx[i]<- 1.1
}
