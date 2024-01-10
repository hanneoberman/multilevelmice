
options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test broydt-x.f ...\n")

# THIS VERSION IS NOT AS PER funconstrain
# btf <- function(x) {
# n <- length(x)
# f <- rep(NA, n)
# f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
# tnm1 <- 2:(n-1)
# f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
# f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
# sum(f*f)
# }

broydt.f <- function(par) {
      n <- length(par)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * par) * par + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * par[2:n]
      fi[2:n] <- fi[2:n] - par[1:(n - 1)]
      sum(fi * fi)
    }

broydt.g <- function(x) {
      n <- length(x)
      if (n < 1) {
        stop("Broyden Tridiagonal: n must be positive")
      }

      fi <- (3 - 2 * x) * x + 1
      fi[1:(n - 1)] <- fi[1:(n - 1)] - 2 * x[2:n]
      fi[2:n] <- fi[2:n] - x[1:(n - 1)]

      grad <- 2 * fi * (3 - 4 * x)
      grad[1:(n - 1)] <- grad[1:(n - 1)] - 2 * fi[2:n]
      grad[2:n] <- grad[2:n] - 4 * fi[1:(n - 1)]

      grad
    }

broydt.h <- function(x) { 
       n <- length(x)
       h <- matrix(0.0, nrow=n, ncol=n)
#       ! For i <- 1
       t  <- ( 3.0 - 2.0*x[1] )*x[1] - 2.0*x[2] + 1.0
       t1 <- 3.0 - 4.0*x[1]
       h[1,1] <-   2.0*( t1 ^ 2 - 4.0*t )
       h[1,2] <- - 4.0*t1
       h[2,2] <-   8.0

       for (i in 2:(n-1)) {
          t  <- ( 3.0 - 2.0*x[i] )*x[i] - x[i-1] - 2.0*x[i+1] + 1.0
          t1 <- 3.0 - 4.0*x[i]
          h[i-1,i-1] <- h[i-1,i-1] + 2.0
          h[i-1,i  ] <- h[i-1,i  ] - 2.0*t1
          h[i,  i  ] <- h[i,  i  ] + 2.0*( t1 ^ 2 - 4.0*t )
          h[i-1,i+1] <- h[i-1,i+1] + 4.0
          h[i  ,i+1] <- h[i  ,i+1] - 4.0*t1
          h[i+1,i+1] <- h[i+1,i+1] + 8.0
       }

#       ! For i <- n
       t  <- ( 3.0 - 2.0*x[n] )*x[n] - x[n-1] + 1.0
       t1 <- 3.0 - 4.0*x[n]
       h[n-1,n-1] <- h[n-1, n-1] + 2.0
       h[n-1,  n] <- h[n-1, n  ] - 2.0*t1
       h[  n,  n] <- h[ n,  n  ] + 2.0*( t1 ^ 2 - 4.0*t )

       for (j in 1:(n-1)) { # symmetrize
         for (k in (j+1):n) {
           h[k,j] <- h[j,k]        
         }
       }
       h
}
x0 <- rep(-1,5)
# f0<-btf(x0)
# f0

ff<-broydt.f(x0)
ff
gg<-broydt.g(x0)
hh<-broydt.h(x0)
gg
hh
library(pracma)
ggn <- grad(broydt.f, x0)
ggn
max(abs(gg-ggn))
hhes<-hessian(broydt.f, x0)
hhes
max(abs(hh-hhes))

hjac<-jacobian(broydt.g, x0)
hjac
max(abs(hh-hjac))
hhes <- hessian(broydt.f, x0)
hhes
library(numDeriv)
ggn <- grad(broydt.f, x0)
ggn
max(abs(gg-ggn))
hhesn <- hessian(broydt.f, x0)
hhesn
max(abs(hh-hhesn))
hjacn<-jacobian(broydt.g, x0)
hjacn
max(abs(hh-hjacn))

# p0 <- rnorm(10, sd=1)
# system.time(ans.optx <- optimx(par=p0, fn=broydt.f,control=list(maxit=25000,save.failures=TRUE,all.methods=TRUE)))[1]
# print(ans.optx)

topm <- opm(x0, broydt.f, broydt.g, broydt.h, method="MOST")
summary(topm, order=value)
