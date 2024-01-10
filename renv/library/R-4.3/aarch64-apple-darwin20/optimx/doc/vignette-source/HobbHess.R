# HobbHess.R
## Optimization test function HOBBS
## Nash and Walker-Smith (1987, 1989) ...
require(optimx)

hobbs.f<- function(x){ # # Hobbs weeds problem -- function
  if (abs(12*x[3]) > 500) { # check computability
    fbad<-.Machine$double.xmax
    return(fbad)
  }
  res<-hobbs.res(x)
  f<-sum(res*res)
}
attr(hobbs.f, "fname")<- "Hobbs"

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

x0 <- c(200, 50, .3)
cat("Good start for Hobbs:")
print(x0)
solx0 <- optimr(x0, hobbs.f, hobbs.g, hobbs.h, method="snewton", hessian=TRUE)
## Note that we exceed count limit, but have answer
proptimr(solx0)
print(eigen(solx0$hessian)$values)
## Note that we exceed count limit, but have answer

## Setting relative check offset larger gets quicker convergence
solx0a <- optimr(x0, hobbs.f, hobbs.g, hobbs.h, method="snewton", 
                  control=list(offset=1000.))
proptimr(solx0a)


x1s <- c(100, 10, .1)
cat("Scaled start for Hobbs:")
print(x1s)
solx1s <- optimr(x1s, hobbs.f, hobbs.g, hobbs.h, method="snewton", hessian=TRUE , control=list(trace=0))
proptimr(solx1s)
print(eigen(solx1s$hessian)$values)
solx1m <- optimr(x1s, hobbs.f, hobbs.g, hobbs.h, method="snewtonm", hessian=TRUE , control=list(trace=0))
proptimr(solx1m)
print(eigen(solx1m$hessian)$values)

cat("Following test fails in snewton with ERROR \n
     -- Not run as function infinite.\n")
x3 <- c(1, 1, 1)
# solx3 <- try(optimr(x3, hobbs.f, hobbs.g, hobbs.h, method="snewton", control=list(trace=4)))
# if ((solx3$convergence != 0) || class(solx3) != "try-error") {
#   proptimr(solx3)
#   print(eigen(solx3$hessian)$values)
# }
# dirx3 <- try(snewton(x3, hobbs.f, hobbs.g, hobbs.h, control=list(trace=4)))
# if ((dirx3$convergence != 0) || class(dirx3) != "try-error") {
#   proptimr(dirx3)
#   print(eigen(dirx3$Hess)$values)
# }
cat("But Marquardt variant succeeds\n")
solx3m <- optimr(x3, hobbs.f, hobbs.g, hobbs.h, method="snewtonm", 
                  hessian=TRUE, control=list(trace=0))
proptimr(solx3m)
print(eigen(solx3m$hessian)$values)
# we could also use nlm and nlminb and call them from optimr

solx3 <- try(optimr(x3, hobbs.f, hobbs.g, hobbs.h, method="snewton", control=list(trace=0)))
if ((class(solx3) != "try-error") && (solx3$convergence == 0)) {
  proptimr(solx3)
  print(eigen(solx3$hessian)$values)
} else cat("solx3 failed!\n")
dirx3 <- try(snewton(x3, hobbs.f, hobbs.g, hobbs.h, control=list(trace=0)))
if ((class(dirx3) != "try-error") && (dirx3$convcode == 0)) {
  proptimr(dirx3)
  print(eigen(dirx3$Hess)$values)
} else cat("dirx3 failed!\n")

