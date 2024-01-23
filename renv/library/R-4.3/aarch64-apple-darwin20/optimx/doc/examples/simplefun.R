# simplefun.R
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
