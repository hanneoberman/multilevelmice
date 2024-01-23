## Optimization test function GENROSE
rm(list = ls())
library(optimx)

genrose.f <- function(x, gs = NULL) {
    # objective function
    ## One generalization of the Rosenbrock banana valley
    #   function (n parameters)
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    fval <- 1 + sum(gs * (x[1:(n - 1)]^2 - x[2:n])^2 + (x[2:n] - 
        1)^2)
    return(fval)
}

genrose.h <- function(x, gs = NULL) {
    ## compute Hessian
    if (is.null(gs)) {
        gs = 100
    }
    n <- length(x)
    hh <- matrix(rep(0, n * n), n, n)
    for (i in 2:n) {
        z1 <- x[i] - x[i - 1] * x[i - 1]
        z2 <- 1 - x[i]
        hh[i, i] <- hh[i, i] + 2 * (gs + 1)
        hh[i - 1, i - 1] <- hh[i - 1, i - 1] - 4 * gs * z1 - 
            4 * gs * x[i - 1] * (-2 * x[i - 1])
        hh[i, i - 1] <- hh[i, i - 1] - 4 * gs * x[i - 1]
        hh[i - 1, i] <- hh[i - 1, i] - 4 * gs * x[i - 1]
    }
    return(hh)
}

genrose.g <- function(x, gs = NULL) {
    # vectorized gradient for genrose.f
    # Ravi Varadhan 2009-04-03
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    gg <- as.vector(rep(0, n))
    tn <- 2:n
    tn1 <- tn - 1
    z1 <- x[tn] - x[tn1]^2
    z2 <- 1 - x[tn]
    gg[tn] <- 2 * (gs * z1 - z2)
    gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
    gg
}

xx <- rep(2, 6)
g6o <- opm(xx, genrose.f, genrose.g, method="MOST", gs = 100)
print(summary(g6o, order=value))

