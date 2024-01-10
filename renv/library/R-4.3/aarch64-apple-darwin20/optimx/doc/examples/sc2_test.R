options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

#########################################################################################
cat("optimx test sc2 ...\n")

sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / 10
}

neg.sc2.f <- function(x){
n <- length(x)
vec <- 1:n
-sum(vec * (exp(x) - x)) / 10
}

neg.sc2.g <- function(x){
n <- length(x)
vec <- 1:n
-vec * (exp(x) - 1) / 10
}

p0 <- runif(50,min=-1, max=1)
system.time(ans.optxf <- opm(par=p0, fn=sc2.f, method="MOST"))[1]
print(summary(ans.optxf, order=value))

system.time(neg.ans.optxf <- optimx(par=p0, fn=neg.sc2.f, 
              control=list(maximize=TRUE)))[1]
print(summary(neg.ans.optxf, order=-value))

system.time(ans.optxg <- opm(par=p0, fn=sc2.f, gr=sc2.g,
   method="MOST"))[1]
print(summary(ans.optxg, order=value))

system.time(neg.ans.optxg <- opm(par=p0, fn=neg.sc2.f, gr=neg.sc2.g,
   method="MOST", control=list(maximize=TRUE)))[1]
print(summary(neg.ans.optxg, order=-value))



