## Trig function trig1507.R
##  author: John C. Nash
# ref  More' Garbow and Hillstrom, 1981, problem 26, from
#  Spedicato (their ref. 25)
#
# This script runs optimx routines on a trigonometric function
# that has multiple minima. Even with some bounds constraints
# on the parameters, different solutions are suggested. 
#
# If interested in exploring and documenting this function, 
# please contact nashjc _at_ uottawa.ca
#
#
rm(list=ls())
require(optimx)

trig.f <- function(x){
  res <- trig.res(x)
  f <- sum(res*res)
#  cat("FV=",f," at ")
#  print(x)
  f
}

trig.res <- function(x){
   n <- length(x)
   i <- 1:n
   res <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
   return(res)
}

trig.jac <- function(x) { # not vectorized. Can it be?
## stop("Not defined")
   n <- length(x)
   J<-matrix(0,n,n)
   for (i in 1:n) {
      for (j in 1:n) {
         J[i,j]<-sin(x[j]) # we overwrite J[i,i]
      }
      J[i,i] <- (1+i) * sin(x[i])  - cos(x[i])
   }
   return(J)
}

trig.g <- function(x) { # unvectorized
  n<-length(x)
  res<-trig.res(x)
  J<-trig.jac(x)
  g<- as.vector(2.0 * ( t(J) %*% res ))
  return(g)
}


x<-rep(2,2)
cat("optim BFGS and optimr Rvmmin from (2,2)\n")
opt2<-optim(x, trig.f, trig.g, method="BFGS")
proptimr(opt2)
cat("optimr Rvmmin from (2,2)\n")
opt2r<-optimr(x, trig.f, trig.g, method="Rvmmin")
proptimr(opt2r)
cat("====================")
x<-rep(2,4)
cat("optim(BFGS) vs optimr(BFGS) from rep(2,4)\n")
opt4<-optim(x, trig.f, trig.g, method="BFGS")
proptimr(opt4)
opt4r<-optimr(x, trig.f, trig.g, method="Rvmmin")
proptimr(opt4r)
cat("====================")
x<-rep(2,8)
cat("optim BFGS vs optrimr Rvmmin from rep(2,8)\n")
opt8<-optim(x, trig.f, trig.g, method="BFGS")
proptimr(opt8)
opt8r<-optimr(x, trig.f, trig.g, method="Rvmmin")
proptimr(opt8r)

cat("opm ALL from rep(2,2) -- several solutions\n")
ttrig2<-opm(rep(2,2), trig.f, trig.g, method="ALL")
summary(ttrig2, order=value)
cat("reorder by p1")
summary(ttrig2, order=p1)

cat("opm ALL from rep(2,4)\n")
ttrig4<-opm(rep(2,4), trig.f, trig.g, method="ALL")
summary(ttrig4, order=value)
cat("reorder by p1")
summary(ttrig4, order=p1)

cat("opm ALL from (2,8)\n")
ttrig8<-opm(rep(2,8), trig.f, trig.g, method="ALL")
summary(ttrig8, order=p1)

cat("opm ALL from rep(2,2) with bounds (.1,pi) on each \n")
ttrig2b <- opm(rep(2,2), trig.f, trig.g, lower=rep(.1,2), upper=rep(pi,2), method="ALL")
summary(ttrig2b, order=p1)

cat("opm ALL from rep(2,4) with bounds (.1,pi) on each\n")
ttrig4b <- opm(rep(2,4), trig.f, trig.g, lower=rep(.1,4), upper=rep(pi,4), method="ALL")
summary(ttrig4b, order=p1)

cat("opm ALL from rep(2,8) with bounds (.1,pi) on each\n")
ttrig8b <- opm(rep(2,8), trig.f, trig.g, lower=rep(.1,8), upper=rep(pi,8), method="ALL")
summary(ttrig8b, order=p1)
