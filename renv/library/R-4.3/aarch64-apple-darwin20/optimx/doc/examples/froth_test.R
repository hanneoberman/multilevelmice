options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


#########################################
cat("optimx test froth-x ...\n")

froth.f <- function(p){
# Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
f <- rep(NA,length(p))
f[1] <- -13 + p[1] + (p[2]*(5 - p[2]) - 2) * p[2]
f[2] <- -29 + p[1] + (p[2]*(1 + p[2]) - 14) * p[2]
sum (f * f)
}

p0 <- rpois(2,10)
system.time(optmf <- opm(par=p0, fn=froth.f, method="MOST"))[1]
print(summary(optmf, order=value))

system.time(optmfgn <- opm(par=p0, fn=froth.f, gr="grpracma", method="MOST"))[1]
print(summary(optmfgn, order=value))

