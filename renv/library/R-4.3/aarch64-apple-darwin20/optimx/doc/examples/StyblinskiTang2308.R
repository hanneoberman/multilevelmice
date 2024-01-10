# StyblinskiTang2308.R
library(optimx)
fnStyblinskiTang <- function(x) 0.5 * sum(x^4 - 16*x^2 + 5*x)
grST <- function(x) {
   2*x^3 - 16*x + 2.5
}
n <- 4
cat("Globalmin ",n," = ",(-39.16599)*length(x0),"\n")
methlist<-c("nlnm", "hjn", "nvm", "ncg", "newuoa", "Nelder-Mead")
tst0n<-opm(x0, fnStyblinskiTang, gr="grcentral", method=methlist, control=list(trace=0))
summary(tst0n, order=value)
set.seed(4321)
x0<-rep(0, n) # start at lower corner
x1<-runif(n)
x2<-runif(n,0,5)
x3<-runif(n,-5,5)
x4<-rep(pi,n)
x5<-rep(n,n)
x6<-rep(-5,n)
stmat<-rbind(x0, x1, x2, x3, x4, x5, x6)
tst0<-opm(x1, fnStyblinskiTang, gr=grST, method=methlist, control=list(trace=0))
summary(tst0, order=value)
tst1<-opm(x1, fnStyblinskiTang, gr="grcentral", method=methlist, control=list(trace=0))
summary(tst1, order=value)
tst2<-opm(x2, fnStyblinskiTang, gr="grcentral", method=methlist, control=list(trace=0))
summary(tst2, order=value)
tst3<-opm(x3, fnStyblinskiTang, gr="grcentral", method=methlist, control=list(trace=0))
summary(tst3, order=value)
tst4<-opm(x4, fnStyblinskiTang, gr="grcentral", method=methlist, control=list(trace=0))
summary(tst4, order=value)
tst5<-opm(x5, fnStyblinskiTang, gr="grcentral", method=methlist, control=list(trace=0))
summary(tst5, order=value)
for (mstr in methlist){
  cat("Method ",mstr,":\n")
  mtst <- multistart(stmat, fnStyblinskiTang, grST, method=mstr)
  print(mtst)
}
