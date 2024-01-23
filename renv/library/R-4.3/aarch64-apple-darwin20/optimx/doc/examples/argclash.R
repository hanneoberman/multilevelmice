# argclash.R -- dotargs have a name that clashes with internals 
# In particular numDeriv grad has "x" as an argument.
# Might be that pracma works, since x0 is used instead of x??
rm(list=ls())
sqmod<-function(z, x){
   nn<-length(z)
   yy<-x^(1:nn)
   f<-sum((yy-z)^2)
#   cat("Fv=",f," at ")
#   print(z)
   f
}
sqmod.g <- function(z, x){
   nn<-length(z)
   yy<-x^(1:nn)
   gg<- 2*(z - yy)
}

library(optimx)
sessionInfo()
nn<-2
st<-rep(0.5, nn)
# No gradient in Nelder-Mead, so it would work
# t2n <- optimx(st, fn=sqmod, method="Nelder-Mead", control=list(trace=4), x=2)
#t2n
# Modify the functions to use variable info
x <- 2
sqmod1 <- function(z){ sqmod(z, x=x) }
dots <- list(x=2)
str(dots)
sqmod2 <- function(z){ sqmod(z, unlist(dots)) }
# simple test
library(numDeriv)
xpar <- c(2,4)
tryg1 <- grad(sqmod1, xpar )
tryg1 # OK
# Compare to actual gradient
cat("sqmod2 for x=2:", sqmod2(xpar), "\n")
print(sqmod.g(xpar, x=2))
tryg2 <- grad(sqmod2, xpar )
tryg2 # OK


t2a <- optimr(st, fn=sqmod, gr=sqmod.g, method="nvm", control=list(trace=4), x=2)
proptimr(t2a)
t2n <- try(optimr(st, fn=sqmod, gr="grnd", method="nvm", control=list(trace=4), x=2))
proptimr(t2n)
t2p <- try(optimr(st, fn=sqmod, gr="grpracma", method="nvm", control=list(trace=4), x=2))
proptimr(t2p)

t2n1 <- try(optimr(st, fn=sqmod1, gr="grnd", method="nvm", control=list(trace=4)))
proptimr(t2n1)
t2p1 <- try(optimr(st, fn=sqmod1, gr="grpracma", method="nvm", control=list(trace=4)))
proptimr(t2p1)

sqq<-function(z, q){
  nn<-length(z)
  yy<-q^(1:nn)
  f<-sum((yy-z)^2)
  #   cat("Fv=",f," at ")
  #   print(z)
  f
}
sqq.g <- function(z, q){
  nn<-length(z)
  yy<-q^(1:nn)
  gg<- 2*(z - yy)
}

q2a <- optimr(st, fn=sqq, gr=sqq.g, method="nvm", control=list(trace=4), q=2)
proptimr(q2a)
q2n <- try(optimr(st, fn=sqq, gr="grnd", method="nvm", control=list(trace=4), q=2))
proptimr(q2n)
q2p <- try(optimr(st, fn=sqq, gr="grpracma", method="nvm", control=list(trace=4), q=2))
proptimr(q2p)
