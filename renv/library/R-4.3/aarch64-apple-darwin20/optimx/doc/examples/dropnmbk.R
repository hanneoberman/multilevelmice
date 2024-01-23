# dropnmbk.R -- Show case where nmbk fails by trying to start on a bound
rm(list=ls())
meth <- c("L-BFGS-B", "nmkb", "Rcgmin")
up <- c(1, 1)
lo <- c(0, 0)
start <- c(0.9999999999, .5)
test <- function(x){sum(x^2)}
library(optimx)
trybd <- opm(start, test, gr="grfwd", method=meth, lower=lo, upper=up, control=list(trace=1))
print(summary(trybd, order=value))
cat("But watch out if starting on bound\n\n")
start <- c(1, .5)
trybd1 <- opm(start, test, gr="grfwd", method=meth, lower=lo, upper=up, control=list(trace=1))
print(trybd1)
cat("And make sure at least 1 method requested\n")
meth<-c() # should cause a failure -- no method
trybd0 <- try(opm(start, test, gr="grfwd", method=meth, lower=lo, upper=up, control=list(trace=1)))
