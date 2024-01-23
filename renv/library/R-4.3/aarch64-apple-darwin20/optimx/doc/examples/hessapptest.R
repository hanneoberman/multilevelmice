# hessapptest.R
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

x0 <- c(-1.2, 1, -2.2, -1)

optim(x0, genrosa.f, genrosa.g, method="BFGS")
optimr(x0, genrosa.f, genrosa.g, method="BFGS")

# mm <- "nlm"
mm <- "snewtonm"
tan <- optimr(x0, fn=genrosa.f, gr=genrosa.g, hess=genrosa.h, method=mm)
proptimr(tan)
tnd <- optimr(x0, fn=genrosa.f, gr=genrosa.g, hess="approx", method=mm)
proptimr(tnd)
tpr <- optimr(x0, fn=genrosa.f, gr=genrosa.g, hess="approx", method=mm, control=list(hesspkg="pracma"))
proptimr(tpr)

