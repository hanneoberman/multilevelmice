# GenRoseHess.R
# genrosa function code -- attempts to match the rosenbrock at gs=100 and x=c(-1.2,1)
genrosa.f<- function(x, gs=NULL){ # objective function
  ## One generalization of the Rosenbrock banana valley function (n parameters)
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  # Note do not at 1.0 so min at 0
  fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}
attr(genrosa.f, "fname")<-"genrosa"

genrosa.g <- function(x, gs=NULL){
  # vectorized gradient for genrose.f
  # Ravi Varadhan 2009-04-03
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  gg <- as.vector(rep(0, n))
  tn <- 2:n
  tn1 <- tn - 1
  z1 <- x[tn] - x[tn1]^2
  z2 <- 1 - x[tn1]
  # f = gs*z1*z1 + z2*z2
  gg[tn] <- 2 * (gs * z1)
  gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1 - 2 *z2 
  return(gg)
}

genrosa.h <- function(x, gs=NULL) { ## compute Hessian
  if(is.null(gs)) { gs=100.0 }
  n <- length(x)
  hh<-matrix(rep(0, n*n),n,n)
  for (i in 2:n) {
    z1<-x[i]-x[i-1]*x[i-1]
    #		z2<-1.0 - x[i-1]
    hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
    hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
    hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
    hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
  }
  return(hh)
}

require(optimx)
cat("Generalized Rosenbrock tests\n")

cat("original n and x0")

x0 <- c(-1.2, 1)
# solorigs <- snewton(x0, genrosa.f, genrosa.g, genrosa.h) # WORKS OK if optimx loaded
solorig <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewton", hessian=TRUE)

proptimr(solorig)
print(eigen(solorig$hessian)$values)
solorigm <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewtonm", hessian=TRUE)
proptimr(solorigm)
print(eigen(solorigm$hessian)$values)

# Start with 50 values of pi and scale factor 10
x0 <- rep(pi, 50)
sol50pi <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewton", 
                  hessian=TRUE, gs=10)
proptimr(sol50pi)
print(eigen(sol50pi$hessian)$values)
hhi <- genrosa.h(sol50pi$par, gs=10)
print(eigen(hhi)$values)
sol50pim <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewtonm", 
                   hessian=TRUE, gs=10)
proptimr(sol50pim)
hhm <- genrosa.h(sol50pim$par, gs=10)
print(eigen(hhm)$values)

# Bounds constraints

lo<-rep(3,50)
up<-rep(4,50)
sol50pimb <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, lower=lo, upper=up, method="snewtonm", 
                     hessian=TRUE, gs=10)
proptimr(sol50pimb)

# approximate hessian
solom01 <- optimr(x0, genrosa.f, gr=NULL, hess="approx", method="snewtonm", hessian=TRUE)
proptimr(solom01)
print(eigen(solom01$hessian)$values)
solomg1 <- optimr(x0, genrosa.f, genrosa.g, hess="approx", method="snewtonm", hessian=TRUE)
proptimr(solomg1)
print(eigen(solomg1$hessian)$values)
# Following should fail
solomrr <- try(optimr(x0, genrosa.f, gr=NULL, hess="rubbish", method="snewtonm", hessian=TRUE))

