# hessian-used.R -- J C Nash 2021-12-16
library(optimx)
# Try to use hessian in optimizations
## Methods
meth <- c("snewton", "snewtonm", "nlm", "nlminb")

#Rosenbrock banana valley function
rb.f <- function(x){
  return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
#gradient
rb.g <- function(x){
  return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}
#Hessian
rb.h <- function(x) {
  a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
  return(matrix(c(a11, a21, a21, 200), 2, 2))
}

x0 <- c(-1.2, 1)

sr <- snewton(x0, fn=rb.f, gr=rb.g, hess=rb.h, control=list(trace=1))
proptimr(sr)
srm <- snewtonm(x0, fn=rb.f, gr=rb.g, hess=rb.h, control=list(trace=1))
proptimr(srm)

meth <- c("snewton", "snewtonm", "nlm", "nlminb")
trb <- opm(x0, rb.f, rb.g, hess=rb.h, method=meth)
summary(trb, order=value)

trbg <- opm(x0, rb.f, rb.g, method=meth)
summary(trbg, order=value)

trbf <- opm(x0, rb.f, method=meth)
summary(trbf, order=value)


#Example 2: Wood function
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
w0 <- c(-3, -1, -3, -1)

wd <- snewton(w0, fn=wood.f, gr=wood.g, hess=wood.h, control=list(trace=1))
proptimr(wd)

wdm <- snewtonm(w0, fn=wood.f, gr=wood.g, hess=wood.h, control=list(trace=1))
proptimr(wdm)


twood <- opm(w0, wood.f, wood.g, hess=wood.h, method=meth)
summary(twood, order=value)

twoodg <- opm(w0, wood.f, wood.g, method=meth)
summary(twoodg, order=value)

twoodf <- opm(w0, wood.f, method=meth)
summary(twoodf, order=value)
