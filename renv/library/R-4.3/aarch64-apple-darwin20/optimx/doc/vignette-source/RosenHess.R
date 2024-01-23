# RosenHess.R
require(optimx)
f <- function(x){ #Rosenbrock banana valley function
  return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
attr(f,"fname")<-"RosenHess"
gr <- function(x){ #gradient
  return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}
h <- function(x) { #Hessian
  a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
  return(matrix(c(a11, a21, a21, 200), 2, 2))
}
x0 <- c(-1.2, 1)
t1 <- snewton(x0, fn=f, gr=gr, hess=h, control=list(trace=0))
proptimr(t1)

# we can also use nlm and nlminb
fght <- function(x){ ## combine f, g and h into single function for nlm
  ff <- f(x)
  gg <- gr(x)
  hh <- h(x)
  attr(ff, "gradient") <- gg
  attr(ff, "hessian") <- hh
  ff
}

t1nlmo <- optimr(x0, f, gr, hess=h, method="nlm", control=list(trace=0))
proptimr(t1nlmo)

t1so <- optimr(x0, f, gr, hess=h, method="snewton", control=list(trace=0))
proptimr(t1so)

t1smo <-  optimr(x0, f, gr, hess=h, method="snewtonm", control=list(trace=0))
proptimr(t1smo)


## nlminb 
tst <- try(t1nlminbo <- optimr(x0, f, gr, hess=h, method="nlminb", 
                               control=list(trace=0)))
if (class(tst) == "try-error"){
  cat("try-error on attempt to run nlminb in optimr()\n")
} else { proptimr(t1nlminbo) }

tstnoh <- try(t1nlminbnoho <- optimr(x0, f, gr, hess=NULL, method="nlminb", 
                               control=list(trace=0)))
if (class(tstnoh) == "try-error"){
  cat("try-error on attempt to run nlminb in optimr()\n")
} else { proptimr(t1nlminbnoho) }

