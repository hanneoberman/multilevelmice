rm(list=ls()) # comment out this line if you do not want the workspace cleared
##  author: John C. Nash
require(optimx)
sessionInfo()
## Optimization test function HOBBS
## Nash and Walker-Smith (1987, 1989) ...

hobbs.f<- function(x){ # # Hobbs weeds problem -- function
    if (abs(12*x[3]) > 500) { # check computability
       fbad<-.Machine$double.xmax
       return(fbad)
    }
    res<-hobbs.res(x)
    f<-sum(res*res)
}


hobbs.res<-function(x){ # Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    t<-1:12
    if(abs(12*x[3])>50) {
       res<-rep(Inf,12)
    } else {
       res<-x[1]/(1+x[2]*exp(-x[3]*t)) - y
    }
}

hobbs.jac<-function(x){ # Jacobian of Hobbs weeds problem
   jj<-matrix(0.0, 12, 3)
   t<-1:12
    yy<-exp(-x[3]*t)
    zz<-1.0/(1+x[2]*yy)
     jj[t,1] <- zz
     jj[t,2] <- -x[1]*zz*zz*yy
     jj[t,3] <- x[1]*zz*zz*yy*x[2]*t
   return(jj)
}

hobbs.g<-function(x){ # gradient of Hobbs weeds problem
    # NOT EFFICIENT TO CALL AGAIN
    jj<-hobbs.jac(x)
    res<-hobbs.res(x)
    gg<-as.vector(2.*t(jj) %*% res)
    return(gg)
}


hobbs.rsd<-function(x) { # Jacobian second derivative
    rsd<-array(0.0, c(12,3,3))
    t<-1:12
    yy<-exp(-x[3]*t)
    zz<-1.0/(1+x[2]*yy)
    rsd[t,1,1]<- 0.0
    rsd[t,2,1]<- -yy*zz*zz
    rsd[t,1,2]<- -yy*zz*zz
    rsd[t,2,2]<- 2.0*x[1]*yy*yy*zz*zz*zz
    rsd[t,3,1]<- t*x[2]*yy*zz*zz
    rsd[t,1,3]<- t*x[2]*yy*zz*zz
    rsd[t,3,2]<- t*x[1]*yy*zz*zz*(1-2*x[2]*yy*zz)
    rsd[t,2,3]<- t*x[1]*yy*zz*zz*(1-2*x[2]*yy*zz)
##    rsd[t,3,3]<- 2*t*t*x[1]*x[2]*x[2]*yy*yy*zz*zz*zz
    rsd[t,3,3]<- -t*t*x[1]*x[2]*yy*zz*zz*(1-2*yy*zz*x[2])
    return(rsd)
}


hobbs.h <- function(x) { ## compute Hessian
#   cat("Hessian not yet available\n")
#   return(NULL)
    H<-matrix(0,3,3)
    res<-hobbs.res(x)
    jj<-hobbs.jac(x)
    rsd<-hobbs.rsd(x)
##    H<-2.0*(t(res) %*% rsd + t(jj) %*% jj)
    for (j in 1:3) {
       for (k in 1:3) {
          for (i in 1:12) {
             H[j,k]<-H[j,k]+res[i]*rsd[i,j,k]
          }
       }
    }
    H<-2*(H + t(jj) %*% jj)
    return(H)
}


allm <- c("BFGS", "CG", "Nelder-Mead",  "nlm", "nlminb", 
          "lbfgsb3c", "Rcgmin", "Rtnmin", "Rvmmin",
          "spg", "ucminf", "bobyqa", "hjkb", "hjn", 
          "subplex")
# Dropped "L-BFGS-B", "newuoa", "nmkb", "snewton", "snewtonm","lbfgs",  as they give trouble
badm <- c("L-BFGS-B", "newuoa", "nmkb", "snewton", "snewtonm","lbfgs")

x0 <- c(200, 50, .3)
# This start seems to be OK for all methods, and most do a reasonable job
cat("Start for Hobbs:")
print(x0)
cat("Initial value of hobbs.f = ",hobbs.f(x0),"\n")

ahobb0 <- opm(x0, hobbs.f, hobbs.g, hess=hobbs.h, method=allm)
print(summary(ahobb0, order=value))
# conv code 1 means method does not consider it has returned a solution
#           3 implies line search failure (may not be serious)

# Several methods fail because f or g becomes Inf.
# Try the "BAD" methods -- which may still work sometimes
x1 <- c(1, 1, 1)
cat("Start for Hobbs:")
print(x1)
cat("Initial value of hobbs.f = ",hobbs.f(x1),"\n")
badhobb0 <- opm(x1, hobbs.f, hobbs.g, hess=hobbs.h, method=badm)
print(summary(badhobb0, order=value))
# Following shows method tries to evaluate function as inadmissible point
LBFGSBhobb0 <- try(ptim(x1, hobbs.f, hobbs.g, method="L-BFGS-B", control=list(trace=2)))

ahobb1 <- opm(x1, hobbs.f, hobbs.g, hess=hobbs.h, method=allm)
print(summary(ahobb1, order=value))

badhobb1 <- opm(x1, hobbs.f, hobbs.g, hess=hobbs.h, method=badm)
print(summary(badhobb1, order=value))

# ahobb1lbfgsb<- optim(x1, hobbs.f, hobbs.g, method="L-BFGS-B", control=list(trace=3))
# Note that optim alone fails in the above

x1s <- c(100, 10, .1)
# L-BFGS-B and lbfgb3 both fail because f or g becomes Inf.
cat("Start for Hobbs:")
print(x1s)
cat("Initial value of hobbs.f = ",hobbs.f(x1s),"\n")
ahobb1s <- opm(x1s, hobbs.f, hobbs.g, hess=hobbs.h, method=allm)
print(summary(ahobb1s, order=value))
bobyqahobb1s <- optimr(x1s, hobbs.f, hobbs.g, method="bobyqa", control=list(maxfeval=20000))
proptimr(bobyqahobb1s)
lbgfsb3cahobb1s <- optimr(x1s, hobbs.f, hobbs.g, method="lbfgsb3c", control=list(maxit=20000))
proptimr(lbgfsb3cahobb1s)
require(lbfgsb3c)
lbgfsb3cD1s <- lbfgsb3c(par=x1s, fn=hobbs.f, gr=hobbs.g, control=list(maxit=2000))
proptimr(lbgfsb3cD1s)
# lbgfsb3cD1sn <- lbfgsb3c(par=x1s, fn=hobbs.f, control=list(maxit=500000))
# Still fails to get a good answer!
# lbgfsb3cD1sn

# Note the following examples run with optim(), not optimr()
tscale <- try(optim(par=c(1,1,1), fn=hobbs.f, method="Nelder-Mead", control=list(trace=1, parscale=c(100, 10, 0.1))))
tscale1 <- try(optim(par=c(1,1,1), fn=hobbs.f, method="Nelder-Mead", control=list(trace=1, parscale=c(1, 1, 1))))
tscalex <- try(optim(par=c(1,1,1), fn=hobbs.f, method="Nelder-Mead", control=list(trace=1, parscale=c(.01, .1, 10))))

btscale <- try(optim(par=c(1,1,1), fn=hobbs.f, method="BFGS", control=list(trace=1, parscale=c(100, 10, 0.1))))
btscale1 <- try(optim(par=c(1,1,1), fn=hobbs.f, method="BFGS", control=list(trace=1, parscale=c(1, 1, 1))))
btscalex <- try(optim(par=c(1,1,1), fn=hobbs.f, method="BFGS", control=list(trace=1, parscale=c(.01, .1, 10))))
btscalexr <- try(optimr(par=c(1,1,1), fn=hobbs.f, gr="grcentral", method="BFGS", control=list(trace=1, parscale=c(.01, .1, 10))))

otscale <- optimr(par=c(1,1,1), fn=hobbs.f, gr=hobbs.g, method="Rvmmin", control=list(trace=1, parscale=c(100, 10, 0.1)))
otscale1 <- optimr(par=c(1,1,1), fn=hobbs.f, gr=hobbs.g, method="Rvmmin", control=list(trace=1, parscale=c(1, 1, 1)))
proptimr(otscale)
proptimr(otscale1)

