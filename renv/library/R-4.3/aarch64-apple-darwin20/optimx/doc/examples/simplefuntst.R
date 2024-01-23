# simplefuntst.R
##  author: John C. Nash
rm(list=ls())
library(optimx)
sessionInfo()

# Simple Test Function 1:
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

strt <- c(1,2,3)
ansfgh <- optimr(strt, simfun.f, simfun.g, simfun.h, method="nlm", 
     hessian=TRUE, control=list(trace=2))
proptimr(ansfgh)

ansall <- opm(strt, simfun.f, simfun.g, simfun.h, method="ALL")
summary(ansall, order=value)
