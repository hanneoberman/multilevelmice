# woodtest.R
##  author: John C. Nash
rm(list=ls())
require(optimx)
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

# wood.fgh <- function(x){ # For info only
#   fval <- wood.f(x)
#   gval <- wood.g(x)
#   hval <- wood.h(x)
#   attr(fval,"gradient") <- gval
#   attr(fval,"hessian")<- hval
#   fval
# }

#################################################
sessionInfo()
cat("\n\n with bounds\n")

x0 <- c(-3,-1,-3,-1) # Wood standard start
lo <- c(-5, -5, -5, -5)
up <- c(0, 10, 10, 10)
up2 <- c(-1, 10, 10, 10)

mlst<-c("snewtonm", "nlminb", "Rvmmin", "Rcgmin", "ncg", "nvm")

bdtst <- opm(x0, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo, upper=up, method=mlst)
summary(bdtst, order=value)

x1 <- rep(-.5,4)

bdtst1 <- opm(x1, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo, upper=up, method=mlst)
summary(bdtst1, order=value)

x2 <- rep(0,4)

bdtst2 <- opm(x2, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo, upper=up, method=mlst)
summary(bdtst2, order=value)

x3 <- rep(-1e6,4)

bdtst3 <- opm(x3, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo, upper=up, method=mlst)
summary(bdtst3, order=value)


mth <- c("snewtonm", "Rvmmin", "L-BFGS-B", "nlminb", "Rcgmin")
wdo2<-opm(x0, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo,
          upper=up2, method=mth, control=list(trace=0))
summary(wdo2, order=value)

mth <- c("ncg", "Rcgmin", "CG")
wdo <- opm(x0, fn=wood.f, gr=wood.g, hess=wood.h, method=mth)
summary(wdo, order=value)

wdob <- opm(x0, fn=wood.f, gr=wood.g, hess=wood.h, method=mth, lower=lo, upper=up2)
summary(wdob, order=value)
