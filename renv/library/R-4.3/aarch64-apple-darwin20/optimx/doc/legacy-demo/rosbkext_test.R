fuzz <- 1e-1 #1e-3 # 1e-10
options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test rosbkext.f ...\n")

rosbkext.f <- function(x){
p <- x
n <- length(p)
sum (100*(p[1:(n-1)]^2 - p[2:n])^2 + (p[1:(n-1)] - 1)^2)
}

p0 <- rnorm(50,sd=2)
system.time(ans.optx <- optimx(par=p0, fn=rosbkext.f, control=list(maxit=2500,save.failures=TRUE,all.methods=TRUE)))[1]

print(ans.optx)

