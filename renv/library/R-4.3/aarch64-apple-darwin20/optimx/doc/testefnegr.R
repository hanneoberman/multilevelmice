# testefnegr.R -- A test of creation of internal function
#  J C Nash -- 2021-12-06
rm(list=ls())

maxfn<-function(x) {# fn to be MAXIMIZED
  # max = 10 at 1:n
  n<-length(x)
  ss<-4^seq(1,n)
  f<-10-sum((x-ss)^2)
  f
}

maxfn.g <- function(x) { # gradient
   n <- length(x)
   ss<-4^seq(1,n)
   gg <- -2*(x-ss)
   gg
}


negmaxfn<-function(x) {# explicit negative of maxfn
  f<-(-1)*maxfn(x)
  return(f)
}
negmaxfn.g<-function(x) {# explicit negative of maxfn
  gg<-(-1)*maxfn.g(x)
  gg
}

optsp <- new.env()
optsp$deps <- 1e-8


grfwd <- function(par, userfn, fbase=NULL, env=optsp, ...) {
  # Forward different gradient approximation
  eps<-env$deps
  if (is.null(fbase)) fbase <- userfn(par, ...)  # ensure we have function value at par
  df <- rep(NA, length(par))
  teps <- eps * (abs(par) + eps)
  for (i in 1:length(par)) {
    dx <- par
    dx[i] <- dx[i] + teps[i]
    df[i] <- (userfn(dx, ...) - fbase)/teps[i]
  }
  df
}


tfg <- function(par, fn, gr, control=list(), ...) {
  control$trace <- 1
  if (is.null(control$fnscale)) control$fnscale <- 1
  npar <- length(par)
  if (is.null(control$parscale)) { 
        pscale <- rep(1,npar)
        if(control$trace > 0) { cat("Unit parameter scaling\n") }
  } else { 
        pscale <- control$parscale 
        if(control$trace > 0) {
          cat("Parameter scaling:")
          print(pscale)
        }
  }
  spar <- par/pscale # scaled parameters
  fnscale <- 1.0 # default to ensure defined and MINIMIZING
  if (! is.null(control$maximize)){ 
      if ( control$maximize ) {fnscale <- -1.0} 
  }
  else { # control$maximize is NULL, so control$fnscale defines behaviour
      fnscale <- control$fnscale # default is 1.0
      if (fnscale < 0) control$maximize<-TRUE # reset maximize if it was null and maximizing
  } # control$maximize has precedence over control$fnscale
  control$fnscale <- fnscale # to ensure set again

  efn <- function(spar, ...) {
      # rely on pscale being defined in this enclosing environment
      par <- spar*pscale
      val <- fn(par, ...) * fnscale
  }

  if (is.character(gr)) { # approximation to gradient
     egr <- function(spar, ...){
       if (control$trace>0) cat("Using numerical approximation '",gr,"' to gradient in optimr()\n")
       if (control$trace > 1) {
         cat("fnscale =",fnscale,"  pscale=")
         print(pscale)
         cat("gr:")
         print(gr)
         cat("par:")
         print(par)
       }
       par <- spar*pscale
       result <- do.call(gr, list(par, userfn=fn, ...)) * pscale * fnscale
    }
  } else { 
    if (is.null(gr)) {egr <- NULL}
    else {
       egr <- function(spar, ...) {
         par <- spar*pscale
         result <- gr(par, ...) * pscale * fnscale
       }
    }
  } # end egr definition

  cat("par =")
  print(par)
  cat("fn = ", fn(par,...),"\n")
  cat("gr =")
  if (is.character(gr))
    print(do.call(gr, list(par, userfn=fn, ...)))
  else
    print(gr(par, ...))
  cat("spar = ")
  print(spar)
  cat("efn = ",efn(spar, ...),"\n")
  cat("egr:")
  print(egr(spar, ...))
 
  return("Done")
}

## Examples/tests

par <- rep(1,4)
cat("Unscaled  - maxfn with fnscale=-1\n")
tfg(par, maxfn, maxfn.g, control=list(parscale=rep(1,length(par)), fnscale=-1))
cat("Unscaled  - negmaxfn\n")
tfg(par, negmaxfn, negmaxfn.g, control=list())
cat("Unscaled  - maxfn with fnscale=-1 and grfwd\n")
tfg(par, maxfn, "grfwd", control=list(parscale=rep(1,length(par)), fnscale=-1))
cat("Unscaled  - negmaxfn and grfwd\n")
tfg(par, negmaxfn, "grfwd", control=list())

par<-4^seq(1,3)
cat("Unscaled  - maxfn with fnscale=-1\n")
tfg(par, maxfn, maxfn.g, control=list(parscale=rep(1,length(par)), fnscale=-1))
cat("Unscaled  - negmaxfn\n")
tfg(par, negmaxfn, negmaxfn.g, control=list())
cat("Unscaled  - maxfn with fnscale=-1 and grfwd\n")
tfg(par, maxfn, "grfwd", control=list(parscale=rep(1,length(par)), fnscale=-1))
cat("Unscaled  - negmaxfn and grfwd\n")
tfg(par, negmaxfn, "grfwd", control=list())



par <- rep(1,4)
cat("Scaled  - maxfn with fnscale=-1\n")
pscale<-4^seq(1,length(par))
tfg(par, maxfn, maxfn.g, control=list(parscale=pscale, fnscale=-1))
cat("Scaled  - negmaxfn\n")
tfg(par, negmaxfn, negmaxfn.g, control=list(parscale=pscale))
cat("Scaled  - maxfn with fnscale=-1 and grfwd\n")
tfg(par, maxfn, "grfwd", control=list(parscale=pscale, fnscale=-1))
cat("Scaled  - negmaxfn and grfwd\n")
tfg(par, negmaxfn, "grfwd", control=list(parscale=pscale))


par <- 4^seq(1,3)
cat("Scaled  - maxfn with fnscale=-1\n")
pscale<-4^seq(1,length(par))
tfg(par, maxfn, maxfn.g, control=list(parscale=pscale, fnscale=-1))
cat("Scaled  - negmaxfn\n")
tfg(par, negmaxfn, negmaxfn.g, control=list(parscale=pscale))
cat("Scaled  - maxfn with fnscale=-1 and grfwd\n")
tfg(par, maxfn, "grfwd", control=list(parscale=pscale, fnscale=-1))
cat("Scaled  - negmaxfn and grfwd\n")
tfg(par, negmaxfn, "grfwd", control=list(parscale=pscale))

