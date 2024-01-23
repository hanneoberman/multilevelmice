
#######################################################################################################
 rosbkext <- function(x){
# Extended Rosenbrock function
 n <- length(x)
 sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
 }

np <- 10
set.seed(123)

p.0 <- rnorm(np)
xm1 <- nmk(fn=rosbkext, par=p.0) # maximum `fevals' is not sufficient to find correct minimum
xm2 <- nmk(fn=rosbkext, par=p.0, control=list(maxfeval=5000)) # finds the correct minimum 
ans.optim <- optim(fn=rosbkext, par=p.0, method="Nelder-Mead", control=list(maxit=5000))   # terminates with inferior estimates
ans.hj <- hjk(fn=rosbkext, par=p.0)   # Hooke-Jeeves algorithm
xmb <- nmkb(fn=rosbkext, par=p.0, lower=-2, upper=2)
 
#######################################################################################################
### A non-smooth problem
nsf <- function(x) {
	f1 <- x[1]^2 + x[2]^2
	f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
	f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
	max(f1, f2, f3)
}

p0 <- rnorm(3)
xm3 <- nmk(fn=nsf, par=p0)
xm3.hj <- hjk(fn=nsf, par=p0)

ans3 <- mads(p0, nsf, control=list(trace=FALSE))
c(xm3$value, xm3.hj$value, ans3$val)

#######################################################################################################
### Another non-smooth problem
rosen <- function(x) {
# Rosen JB & Suzuki S (1965), Construction of non-linear programming test problems, Comm. ACM, 8, p. 113
	f1 <- x[1]^2 + x[2]^2 + 2*x[3]^2 + x[4]^2 - 5*x[1] - 5*x[2] - 21*x[3] + 7*x[4]
	f2 <- f1 + 10 * (sum(x^2) + x[1] - x[2] + x[3] - x[4] - 8)
	f3 <- f1 + 10 * (sum(x^2) + x[2]^2 + x[4]^2 - x[1] - x[4] - 10)
	f4 <- f1 + 10 * (sum(x^2) + x[1]^2 - x[4]^2 + 2*x[1] - x[2] - x[4] - 5)
	max(f1, f2, f3, f4)
}
# Global minimum value is -44 @ (0, 1, 2, -1)

p0 <- rnorm(4)
xm4 <- nmk(fn=rosen, par=p0)
xm4.hj <- hjk(fn=rosen, par=p0)
xm4b <- nmkb(fn=rosen, par=p0, lower=-2, upper=3)
ans3 <- mads(p0, rosen, control=list(trace=FALSE))

#######################################################################################################
### Non-smooth problem #3
hald <- function(x) {
#Hald J & Madsen K (1981), Combined LP and quasi-Newton methods for minimax optimization, Mathematical Programming, 20, p.42-62.
	i <- 1:21
	t <- -1 + (i - 1)/10
	f <- (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
	max(abs(f))
	}
# Correct solution:  x* = (
# Minimum value = 0.002

p0 <- runif(5)
xm5 <- nmk(fn=hald, par=p0)
xm5.hj <- hjk(fn=hald, par=p0)
xm5b <- nmkb(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)
ans3 <- mads(p0, hald, control=list(trace=FALSE))

#################################
## Rosenbrock Banana function
#
fr <- function(x) {   
  n <- length(x)
  x1 <- x[2:n]
  x2 <- x[1:(n-1)]
  sum(100 * (x2 - x1 * x1)^2 + (1 - x1)^2)
}

n <- 10
p0 <- runif(n, 0, 2)

ans1 <- nmk(p0, fr, control=list(maxfeval=20000))
ans2 <- hjk(p0, fr, control=list(maxfeval=20000))
ans3 <- mads(p0, fr, control=list(trace=FALSE))
c(ans1$value, ans2$value, ans3$val)

################################################
# EVD52 
evd52 <- function(x){
  f <- rep(NA, 6)
  f[1] <- sum(x[1:3]^2) - 1
  f[2] <- sum(x[1:2]^2) + (x[3] - 2)^2
  f[3] <- sum(x[1:3]) - 1
  f[4] <- x[1] + x[2] - x[3] + 1
  f[5] <- 2*x[1]^3 + 6*x[2]^2 + 2*(5*x[3] - x[1] + 1)^2
  f[6] <- x[1]^2 - 9*x[3]
  return(max(f))
}

# True mimimum = 3.5997193
p0 <- runif(6)

ans1 <- nmk(p0, evd52, control=list(maxfeval=20000))
ans2 <- hjk(p0, evd52, control=list(maxfeval=20000))
ans3 <- mads(p0, evd52, control=list(trace=FALSE))
c(ans1$value, ans2$value, ans3$val)

###################################################
hs78 <- function(x){
  f <- rep(NA, 3)
  f[1] <- sum(x^2) - 10
  f[2] <- x[2]*x[3] - 5*x[4]*x[5]
  f[3] <- x[1]^3 + x[2]^3 + 1
  F <- prod(x) + 10*sum(abs(f))
  return(F)
}

# True mimimum = -2.9197004

p0 <- c(-2,1.5,2,-1,-1) + runif(5)

ans1 <- nmk(p0, hs78, control=list(maxfeval=20000))
ans2 <- hjk(p0, hs78, control=list(maxfeval=20000))
ans3 <- mads(p0, hs78, control=list(trace=FALSE))
c(ans1$value, ans2$value, ans3$val)

###################################################
elattar <- function(x){
  i <- 1:51
  ti <- 0.1*(i-1)
  yi <- 0.5*exp(-ti) - exp(-2*ti) + 0.5*exp(-3*ti) + 1.5*exp(-1.5*ti)*sin(7*ti) + exp(-2.5*ti)*sin(5*ti)
  F <- sum(abs(x[1]*exp(-x[2]*ti)*cos(x[3]*ti + x[4]) + x[5]*exp(-x[6]*ti) - yi))
    return(F)
}

# True mimimum = 0.5598131

p0 <- c(2,2,7,0,-2,1) + runif(6)

ans1 <- nmk(p0, elattar, control=list(maxfeval=20000, regsimp=TRUE))
ans2 <- hjk(p0, elattar, control=list(maxfeval=20000))
ans3 <- mads(p0, elattar, control=list(trace=FALSE))
c(ans1$value, ans2$value, ans3$val)
