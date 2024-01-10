# snewtonbtest.R -- see if we can run bounds on snewtonmb
rm(list=ls())
# sink(file="~/temp/snbt.txt", split=TRUE)
library(optimx)
# source("optimx/tests/simplefun.R")
# simplefun.R
##  author: John C. Nash

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


sessionInfo()
strt<-c(1,2,3,4)
lo<-c(.1, .2, .3, .4)
up<-c(2, 3, 4, 5)

# Should FAIL -- wrong start
t1snmallmd <- try(snewtonm(strt, simfun.f, simfun.g, simfun.h, lower=lo, 
                           upper=lo, control=list(trace=3)))

# But this works, since optimr allows shift to bound
t1snmallm <- try(optimr(strt, simfun.f, simfun.g, simfun.h, lower=lo, 
                       upper=lo, method="snewtonm", control=list(trace=3)))
proptimr(t1snmallm)

# All masked -- but output does not show M on all params with direct call to snewtonm
t1snmallmds <- try(snewtonm(par=lo, simfun.f, simfun.g, simfun.h, lower=lo, 
                          upper=lo, control=list(trace=1)))
proptimr(t1snmallmds)

t1snmallms <- try(optimr(strt, simfun.f, simfun.g, simfun.h, lower=lo, 
                        upper=lo, method="snewtonm", control=list(trace=3, shift2bound=FALSE)))
proptimr(t1snmallms)




t1vm<-optimr(strt, simfun.f, simfun.g, lower=lo, 
             upper=up, method="Rvmmin", control=list(trace=1))
proptimr(t1vm) # Note L rather than M, since bounds not equal

x1vm<-t1vm$par
bc1vm<-bmchk(x1vm, lower=lo, upper=up)
bc1vm$bchar
attr(t1vm$par,"status") <- bc1vm$bchar
proptimr(t1vm)


t1snmb <- try(snewtonm(strt, simfun.f, simfun.g, simfun.h, lower=lo, 
             upper=up, control=list(trace=1)))
t1snmb$par<-bmchk(t1snmb$par, lower=lo, upper=up)$bvec
proptimr(t1snmb)

cat("All parameters masked test\n")
lo<-c(.1, 3, .3, .4)
# strt<-c(1, 3, 3, 4)
strt<-lo

t1snmallms <- try(snewtonm(par=strt, simfun.f, simfun.g, simfun.h, lower=lo, 
                          upper=lo, control=list(trace=1)))


proptimr(t1snmallms)

t1vmm<-Rvmmin(strt, simfun.f, simfun.g, lower=lo, upper=up, 
              method="Rvmmin", control=list(trace=1))
proptimr(t1vmm)


t1snmm <- try(snewtonm(strt, simfun.f, simfun.g, simfun.h, lower=lo, 
                       upper=up, control=list(trace=1)))
proptimr(t1snmm)

# sink()
