# trystarttests.R -- see if we can run start tests
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

library(optimx)
sessionInfo()
meth<-c("Rvmmin", "Rtnmin", "Rcgmin", "Rvmmin")
strt<-c(1,2,3,4)
# sink(file="~/temp/st1.txt", split=TRUE)
t1 <- try(optimx(strt, simfun.f, simfun.g, method=meth, control=list(trace=1, starttests=TRUE)))
# sink()
# sink(file="~/temp/st2.txt", split=TRUE)
t2 <- try(optimx(strt, simfun.f, simfun.g, method=meth, control=list(trace=1, starttests=FALSE))) 
# sink(
