# optimrgrapprox.R
##  author: John C. Nash
rm(list=ls())
require(optimx)
sessionInfo()
# test to see if grapprox safeguards work

jones<-function(xx){
   x<-xx[1]
   y<-xx[2]
   ff<-sin(x*x/2 - y*y/4)*cos(2*x-exp(y))
   ff<- -ff
}

# Note: No hessian, so newton methods will fail
xx<-0.5*c(pi,pi)
ans<-optimr(xx, jones, "grcentral", method="Rvmmin", control=list(trace=0))
proptimr(ans)

# non-existent method. Should fail
ans<-try(optimr(xx, jones, "grrubbish", method="Rvmmin", control=list(trace=0)))
proptimr(ans)

