options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test chen-x.f ...\n")

chen.f <- function(x) {
v <- log(x) + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}

p0 <- rexp(10)
system.time(chenopf <- opm(par=p0, fn=chen.f, lower=0, upper=Inf, method="MOST"))[1]
print(summary(chenopf, order=value))
system.time(chenopfg <- opm(par=p0, fn=chen.f, gr="grpracma", lower=0, upper=Inf, method="MOST"))[1]
print(summary(chenopfg, order=value))
