## cyq_test.R -- Fletcher's Chebyquad problem on optimx::opm()
# bounds seem to cause trouble

# rm(list = ls())
library(optimx)

# Ref: Fletcher, R. (1965) Function minimization without
#   calculating derivatives -- a review,
#         Computer J., 8, 33-41.

# Note we do not have all components here e.g., .jsd, .h

cyq.f = function(par) {
  n <- length(par)
  if (n < 1) {
    stop("Chebyquad: n must be positive")
  }
  
  # y is the shifted x
  y <- 2 * par - 1
  
  t0 <- rep(1, n)
  t1 <- y
  ti <- t1
  fsum <- 0
  for (i in 1:n) {
    if (i > 1) {
      ti <- 2 * y * t1 - t0
      t0 <- t1
      t1 <- ti
    }
    
    fi <- sum(ti) / n
    if (i %% 2 == 0) {
      fi <- fi + 1 / (i * i - 1)
    }
    fsum <- fsum + fi * fi
  }
  
  fsum
}

cyq.g = function(par) {
  n <- length(par)
  if (n < 1) {
    stop("Chebyquad: n must be positive")
  }
  
  y <- 2 * par - 1
  
  g0 <- rep(0, n)
  # grad of T1 wrt to shift is 2, but we'll account for that later
  g1 <- rep(1, n)
  
  t0 <- rep(1, n)
  t1 <- y
  grad <- rep(0, n)
  
  gi <- g1
  ti <- t1
  for (i in 1:n) {
    if (i > 1) {
      # Gn = 2T_old + 2xG_old - G_old_old
      gi <- 2 * t1 + 2 * y * g1 - g0
      g0 <- g1
      g1 <- gi
      
      # Update T
      ti <- 2 * y * t1 - t0
      t0 <- t1
      t1 <- ti
    }
    
    # Correct for the fact we want want deriv wrt x, not y
    # G wrt x is 2 * G wrt shifted (y)
    gi <- 2 * gi / n
    
    fi <- sum(ti) / n
    if (i %% 2 == 0) {
      fi <- fi + 1 / (i * i - 1)
    }
    
    grad <- grad + 2 * fi * gi
  }
  grad
}


cyq.h = function(x) { 
  n <- length(x)
  h <- matrix(0.0, nrow=n, ncol=n)
  fvec <- rep(0.0,n)
  gvec <- rep(0.0,n)
  for (j in 1:n) {
    t1 <- 1.0
    t2 <- 2.0*x[j] - 1.0
    t <- 2.0*t2
    for (i in 1:n){ 
      fvec[i] <- fvec[i] + t2
      th <- t*t2 - t1
      t1 <- t2
      t2 <- th
    }
  }
  d1 <- 1.0/n
  for (i in 1:n) {
    fvec[i] <- d1*fvec[i]
    if ( (i %% 2) == 0 ) {
      fvec[i] <- fvec[i] + 1.0/( i^2 - 1.0 )
    }
  }
  d2 <- 2.0*d1
  for (j in 1:n) { 
    h[j,j] <- 4.0*d1
    t1 <- 1.0
    t2 <- 2.0*x[j] - 1.0
    t <- 2.0*t2
    s1 <- 0.0
    s2 <- 2.0
    p1 <- 0.0
    p2 <- 0.0
    gvec[1] <- s2
    for (i in 2:n) {
      th <- 4.0*t2 + t*s2 - s1
      s1 <- s2
      s2 <- th
      th <- t*t2 - t1
      t1 <- t2
      t2 <- th
      th <- 8.0*s1 + t*p2 - p1
      p1 <- p2
      p2 <- th
      gvec[i] <- s2
      h[j,j] <- h[j,j] + fvec[i]*th + d1*s2 ^ 2
    }
    h[j,j] <- d2*h[j,j]
    if (j > 1) { # needed for R
      for (k in 1:(j-1)) {
        h[k,j] <- 0.0
        tt1 <- 1.0
        tt2 <- 2.0*x[k] - 1.0
        tt <- 2.0*tt2
        ss1 <- 0.0
        ss2 <- 2.0
        for (i in 1:n) {
          h[k,j] <- h[k,j] + ss2*gvec[i]
          tth <- 4.0*tt2 + tt*ss2 - ss1
          ss1 <- ss2
          ss2 <- tth
          tth <- tt*tt2 - tt1
          tt1 <- tt2
          tt2 <- tth
        }
        h[k,j] <- d2*d1*h[k,j]
      }
    } # end if
  }
  
  for (j in 1:(n-1)) { # symmetrize
    for (k in (j+1):n) {
      h[k,j] <- h[j,k]        
    }
  }
  h
}

cat("Fletcher chebyquad function in file cyq.R\n")

# nn <- c(2, 3, 5, 8, 10, 15, 20)  # added 15, dropped 30 20100525 to reduce testing time

nn <- c(2, 5)

for (n in nn) {
    cat("Chebyquad in ", n, " parameters\n")
    afname <- paste("acyq", n, "G.txt", sep = "")
    lower <- rep(-20, n)
    upper <- rep(20, n)
    bdmsk <- rep(1, n)  # free all parameters
    x0 <- 1:n
    x0 <- x0/(n + 1)  # Initial value suggested by Fletcher
    ut <- system.time(ans <- opm(x0, cyq.f, cyq.g, cyq.h, method="MOST",
                lower=lower, upper=upper))[1]
    print(summary(ans, order=value))
    cat("time = ", ut, "\n")
}
