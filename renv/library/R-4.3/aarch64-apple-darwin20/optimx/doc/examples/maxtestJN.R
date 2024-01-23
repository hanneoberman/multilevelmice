# maxtestJN.R built on maxtestj.R
# check that optimx runs maximize for all methods

rm(list=ls())
library(optimx)

maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}


negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

x0<-rep(pi,4)
print(maxfn(x0))
print(negmaxfn(x0))
x00<-c(1,2,3,4) # solution
print(maxfn(x00))
print(negmaxfn(x00))
amx<-opm(x0,maxfn,method="MOST", control=list(maximize=TRUE,save.failures=TRUE,trace=0))
print(summary(amx, order=-value))

ans.mxn<-opm(x0,negmaxfn,method="MOST", control=list(save.failures=TRUE))
print(summary(ans.mxn, order = value))

# Test if things work when we provide the solution!
ans.mx0<-opm(x00,maxfn,method="MOST", control=list(maximize=TRUE,save.failures=TRUE))
print(summary(ans.mx0, order = value))

ans.mx0n<-opm(x00,negmaxfn,method="MOST", control=list(save.failures=TRUE))
print(summary(ans.mx0n, order = value))



