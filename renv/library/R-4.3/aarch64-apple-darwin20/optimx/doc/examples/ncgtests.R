# ncgtests.R
##  author: John C. Nash
rm(list=ls()) # comment out this line if you do not want the workspace cleared
require(optimx)
# source("optimx/R/Rcgminb.R")
sessionInfo()

# simplefun.R
##  author: John C. Nash
simfun.f = function(x) { 
  fun <- sum(x^2 )
  #	print(c(x = x, fun = fun))
  fun
}
simfun.g = function(x) { 
  grad<-2.0*x
  grad
}
simfun.h = function(x) { 
  n<-length(x)
  t<-rep(2.0,n)
  hess<-diag(t)
}

n<-4
lo<-rep(0,n)
up<-lo # to get arrays set
bdmsk<-rep(1,n)
for (i in 1:n) { 
  lo[i]<-1.0*(i-1)*(n-1)/n
  up[i]<-1.0*i*(n+1)/n
}
x0<-0.5*(lo+up)
cat("Now force a mask upper=lower for parameter 3 and see what happens\n")
lo[3] <- up[3]
x0[3] <- lo[3] # MUST reset parameter also
ncgbdm <- optimr(x0, simfun.f, simfun.g, lower=lo, upper=up, method="ncg", 
                 control=list(trace=4, watch=TRUE))
proptimr(ncgbdm)

# reset
for (i in 1:n) { 
  lo[i]<-1.0*(i-1)*(n-1)/n
  up[i]<-1.0*i*(n+1)/n
}
x0<-0.5*(lo+up)
sncgb <- optimr(x0, fn=simfun.f, gr=simfun.g, lower=lo, upper=up, method="ncg", 
            control=list(trace=4, maxit=600))
proptimr(sncgb)

sbvm <- optimr(x0, fn=simfun.f, gr=simfun.g,lower=lo,
               upper=up, method="Rvmmin", control=list(trace=1))
proptimr(sbvm)
# conv code 2 means point with small gradient found

# Extended Rosenbrock Function ex_rosen.R from 
# https://github.com/jlmelville/funconstrain by jlmelville
# Test function 21 from the More', Garbow and Hillstrom paper.
#
# The objective function is the sum of \code{m} functions, each of \code{n}
# parameters.
#
xrosn.f = function(par) {
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

xrosn.g = function(par) {
  n <- length(par)
  if (n %% 2 != 0) {
    stop("Extended Rosenbrock: n must be even")
  }
  grad <- rep(0, n)
  for (i in 1:(n / 2)) {
    p2 <- 2 * i
    p1 <- p2 - 1
    xx <- par[p1] * par[p1]
    
    yx <- par[p2] - xx
    f_p1 <- 10 * yx
    f_p2 <- 1 - par[p1]
    grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
    grad[p2] <- grad[p2] + 200 * yx
  }
  grad
}


n<-4
lo<-rep(-2,n)
up<- -lo # to get arrays set
bdmsk<-rep(1,n)
for (i in 1:n) { 
  lo[i]<-1.2*(i-1)*(n-1)/n
  up[i]<-3*i*(n+1)/n
}
xx<-0.5*(lo+up)
cat("lower:"); print(lo)
cat("upper:"); print(up)
cat("start:"); print(xx)

xncg <- optimr(xx, fn=xrosn.f, gr=xrosn.g, lower=lo, upper=up, method="ncg", control=list(trace=1, maxit=1000))
proptimr(xncg)
xbvm <- optimr(xx, fn=xrosn.f, gr=xrosn.g,lower=lo,
               upper=up, method="Rvmmin", control=list(trace=1))
proptimr(xbvm)

xbvm0 <- optimr(xx, fn=xrosn.f, gr=xrosn.g, method="Rvmmin", control=list(trace=1))
proptimr(xbvm0)


xall <- opm(xx, fn=xrosn.f, gr=xrosn.g, lower=lo, upper=up, method="ALL")
summary(xall, order=value)


# woodfn.R
##  author: John C. Nash
#Example: Wood function
#
wood.f <- function(x){
  res <- 100*(x[1]^2-x[2])^2+(1-x[1])^2+90*(x[3]^2-x[4])^2+(1-x[3])^2+
    10.1*((1-x[2])^2+(1-x[4])^2)+19.8*(1-x[2])*(1-x[4])
  return(res)
}
#gradient:
wood.g <- function(x){
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  return(c(g1,g2,g3,g4))
}
#hessian:
wood.h <- function(x){
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2; h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  H <- matrix(c(h11,h12,h13,h14,h12,h22,h23,h24,
                h13,h23,h33,h34,h14,h24,h34,h44),ncol=4)
  return(H)
}

#################################################
x0 <- c(-3,-1,-3,-1) # Wood standard start
lo <- c(-5, -5, -5, -5)
up <- c(0, 10, 10, 10)

wncg <- optimr(x0, fn=wood.f, gr=wood.g, method="ncg", lower=lo, upper=up, control=list(trace=1, stepredn=.2, maxit=600))
proptimr(wncg)

wbvm <- optimr(x0, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo,
               upper=up, method="Rvmmin", control=list(trace=1))
proptimr(wbvm)

