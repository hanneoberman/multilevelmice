# SimpHess.R
x0<-c(1,2,3,4)
lo <- x0-0.5
up <- x0+1.0
fnt <- function(x, fscale=10){
  yy <- length(x):1
  val <- sum((yy*x)^2)*fscale
  val
}
attr(fnt,"fname")<-"SimpHess"
grt <- function(x, fscale=10){
  nn <- length(x)
  yy <- nn:1
  #    gg <- rep(NA,nn)
  gg <- 2*(yy^2)*x*fscale
  gg
}
hesst <- function(x, fscale=10){
  nn <- length(x)
  yy <- nn:1
  hh <- diag(2*yy^2*fscale)
  hh
}

require(optimx)
t1 <- optimr(x0, fnt, grt, hesst, method="snewton", control=list(trace=0, usexxxmeth=TRUE), fscale=3.0)
proptimr(t1)
t1m <- optimr(x0, fnt, grt, hesst, method="snewtonm", control=list(trace=0), fscale=3.0)
proptimr(t1m)
# Check alternate name works OK
t1mm <- optimr(x0, fnt, grt, hesst, method="snewtm", control=list(trace=0), fscale=3.0)
proptimr(t1mm)

t1nlmo <- optimr(x0, fnt, grt, hess=hesst, method="nlm", fscale=3.0, 
                 control=list(trace=0))
proptimr(t1nlmo)

## BUT ... nlminb may not be using a true Newton-type method
tst <- try(t1nlminbo <- optimr(x0, fnt, grt, hess=hesst, method="nlminb", 
                               fscale=3.0, control=list(trace=0)))

# Bounded
tstb <- try(t1nlminbb <- optimr(x0, fnt, grt, hess=hesst, method="nlminb", 
                lower=lo, upper=up, fscale=3.0, control=list(trace=0)))
proptimr(t1nlminbb) 

t1smb <-  optimr(x0, fnt, grt, hess=hesst, method="snewtonm", fscale=3.0, 
                 lower=lo, upper=up, control=list(trace=0))
proptimr(t1smb)

# Masked
lo[1]<-x0[1]
up[1]<-x0[1]
lo[4]<-x0[4]
up[4]<-x0[4]
tstm <- try(t1nlminbm <- optimr(x0, fnt, grt, hess=hesst, method="nlminb", 
                                 lower=lo, upper=up, fscale=3.0, control=list(trace=0)))
proptimr(t1nlminbm) 

t1smm <-  optimr(x0, fnt, grt, hess=hesst, method="snewtonm", fscale=3.0, 
                 lower=lo, upper=up, control=list(trace=0))
proptimr(t1smm)
