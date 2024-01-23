# J C Nash 2023-6-15
options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx::opm test brown-x.f ...\n")
cat("Some methods fail on this problem.")

brown.f <- function(x) {
p <- x
n <- length(p)
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)
sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))
}

npar<-24 # Down from 500
p0 <- rnorm(npar,sd=2)
system.time(ans.opm <- opm(par=p0, fn=brown.f, method="MOST", control=list(save.failures=TRUE, maxit=2500)))[1]
print(summary(ans.opm, order=value))
# aa <- optimr(par=p0, fn=brown.f, gr="grpracma", method="nlm", control=list(trace=1))
# aa <- optimr(par=p0, fn=brown.f, method="nlm", control=list(trace=1))
# aa
