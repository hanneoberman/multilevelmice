# 220610 NOT stopping when keepgoing is FALSE ??
snewtonm<-function(par,fn,gr,hess,lower=NULL, upper=NULL, 
           control=list(trace=0, maxit=500),...) {
## Safeguarded Newton minimizer with Marquardt stabilization
##
##Input
##       - par is the initial set of parameter values in vector
##       - fn is the function we wish to minimize
##       - gr is gradient function (required)
##       - hess is hessian function (required)
##       - control is control list (see ctrldefault.R)
##       - ... (dotargs) is exogenous data used in the function fn
##Output (list) -- Note that this does not perfectly match optim() output!!
##       - xs is the parameter vector at the minimum
##       - fv is the fn evaluated at xs
##       - grd is the gradient value (vector)
##       - Hess is the Hessian matrix (note capitalization of this list element name)
##       - niter is the number of interations needed (gradient and Hessian evals).
##       - counts a list of work measures 
##             niter = number of "iterations" i.e., Newton steps
##             nfn   = number of function evaluations
##             ngr   = number of gradient evaluations
##             nhess = number of hessian evaluations
##       - convcode is a number indicating status at termination (0 means OK)
##       - message is a text string returned to indicate status on termination 

npar <- length(par) # par are current best parameters
ctrl <- ctrldefault(npar) # get all the defaults
ncontrol <- names(control) # get the names of controls passed to routine
nctrld <- names(ctrl) # and those in dhe defaults
for (onename in nctrld) { # check each name, copy values NOT supplied
  if (! (onename %in% ncontrol)) {
    control[onename]<-ctrl[onename]
  }
}

trace <- control$trace # convenience to simplify code
if (trace > 0) { dispdefault(control) } # display all the controls
reltest<- control$reltest 
ceps <- .Machine$double.eps * reltest # used for testing closeness to bounds of parameters

#################################################################
  bm <- bmchk(par, lower, upper, trace, bdmsk=NULL) # check bounds
  if (trace > 1) { cat("snewtonm bm:"); print(bm) }
  bounds <- bm$bounds # TRUE if we have bounds
  if (bounds) { # if have bounds, then set some variables
    if (bm$parchanged) stop("Parameters violate bounds")
    nolower <- bm$nolower
    noupper <- bm$noupper
    lower<-bm$lower # ensures fully expanded
    upper<-bm$upper
    # if (trace > 1) { ## NOTE: optimr does the notification
    #    cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, " bounds = ", bounds, "\n") 
    # }
    bdmsk <- bm$bdmsk
    # if (trace > 2) { cat("bdmsk:"); print(bdmsk) }
    if (bm$parchanged) {
        par <- bm$bvec # change parameters
        warning("Parameter(s) changed to nearest bound")
    }
  } # end section for bounds setup
## top of algorithm
  convcode <- 0 # return with good min if 0
  nfn <- 0 # counters for function, gradient, hessian and "iterations" (Newton eqn solves)
  ngr <- 0
  nhess <- 0
  itn <- 0
  eps0<-.Machine$double.eps
  eps <- 10*eps0
  fval <- fn(par, ...) # evaluate function. ?? try()?
  nfn <- nfn + 1
  gtest <- eps*(abs(fval) + 10.0) # tolerance for for gradient test?
  if (trace > 0) { cat("  f0=",fval,"  at  "); print(par)  }
  lambdamin<-(eps0^(1/4)) ## Could do better?!! ?? ctrldefault
  laminc <- control$stepinc
  lamdec <- control$stepredn
  lambda <- lambdamin/lamdec ##  Could do better?!!
  if(trace > 2) cat("laminc, lamdec, lambda:",laminc," ",lamdec," ",lambda,"\n")
  keepgoing <- TRUE # loop control
  while ( (itn <= control$maxit) && keepgoing) { ## main loop  20220608: need to check fval < fbest
    fbest <- fval
    lambda <- lambda * lamdec # decrease Marquardt parameter (initial value will be lamdamin)
    if (trace>0) {cat("ifgh:",itn," ",nfn," ",ngr," ",nhess," fbest=",fbest,"\n")}
    if (trace>1) {cat("par:"); print(par) }
    grd<-gr(par,...) # gradient. ?? use try()?
    ngr <- ngr + 1
    if (trace > 1) { cat("gradient:"); print(grd) }
    if (bounds) {
      for (j in 1:npar) { 
        if ((bdmsk[j] == 0)) { grd[j] <- 0 } # masked; gradient component is zero
        else { # bounds not masks
          if (bdmsk[j] == 1) {
            if (trace > 1) { cat("Parameter ", j, " is free\n") }
          } else {
            if ((bdmsk[j] + 2) * grd[j] < 0) {
               # test for -ve gradient at upper bound, +ve at lower bound
               grd[j] <- 0  # so active constraint and zero gradient component
               if (trace > 2) { cat("fixing parameter ", j, "\n") }
            }
            else {
              bdmsk[j] <- 1  # freeing parameter j
              if (trace > 2) { cat("freeing parameter ", j, "\n") }
            }
          } # end else not free
        } # else bounds not masks
      }  # end masking loop on j
    } # if bounds
    if (trace > 1) { cat("constrained gradient:"); print(grd) }
    if (max(abs(grd)) <= gtest) {
       if (trace > 0) cat("Small or zero (constrained?) gradient\n")
       keepgoing <- FALSE # convergence; small gradient
       xn <- par # to ensure we have output parameters
       H <- matrix(0, nrow=npar, ncol=npar)
       break
    }
    H<-hess(par,...) # ?? do we need to adjust for bounds / masks?
    nhess <- nhess + 1
    if (trace > 2) {cat(nhess,"  H:"); print(H)}
    badstep<-TRUE # keep this until we have a good step
#    cat("badstep, keepgoing:",badstep, keepgoing,"\n")
    while (badstep && keepgoing){ # INNER loop
     # we exit inner loop EITHER with good step (new lower point) OR termination
       itn <- itn+1 # count solves of Newton equations; this is <= nhess
       Haug<-H + diag(npar)*lambda # To avoid singularity. UNSCALED HERE
       solveOK <- TRUE
       stp<-try(solve(Haug, -grd)) # ?? put in try()??
       if ( inherits( stp, "try-error") ) {
          solveOK<-FALSE
          keepgoing<-FALSE # ?? should we stop in this situation?
          warning("Solve failure of Newton equations")
#          cat("H:"); print(H)
#          cat("lambda=",lambda,"  for npar=",npar,"\n")
#          tmp <- readline("Cont.")
       }
       stp[which(grd == 0)] <- 0 # Constrained step 230622
       # Here we need to project step to deal with constraints
       changed <- NA
# 	cat("solveOK=",solveOK," badstep=",badstep,"  keepgoing=",keepgoing,"\n") 
       if (solveOK && keepgoing) { # the solution is OK, but ...
         if (trace > 2) { cat("stp:"); print(stp) }
         # apply mask and current bounds constraint ?? do we need?
         gradproj <- sum(stp*grd) # gradient projection
         if (trace > 1) { cat("gradproj=",gradproj,"\n") }
         if (gradproj < 0) { # don't stop("bad gradproj")
           slen <- 1 # Start with unit steplength
           if (bounds) { # Box constraint -- adjust step length BEFORE trying 
             for (i in 1:npar) { # loop on parameters -- vectorize?
               if ((bdmsk[i] == 1) && (stp[i] != 0)) { # only free parameters / stp != 0
                 if (stp[i] < 0) { # going down. Look at lower bound
                   trystep <- (lower[i] - par[i])/stp[i]  # stp[i] < 0 so this is positive
                 }
                 else { # going up, check upper bound
                   trystep <- (upper[i] - par[i])/stp[i]  # stp[i] > 0 so this is positive
                 }
                 slen <- min(slen, trystep)  # reduce as necessary
               }  # end slen reduction
             }  # end loop on i to reduce step length
           } # end bounds
           if ( (trace>1) && (slen < 1) ) { cat("Step length reduced to ",slen,"\n") }
           xn <- par + stp*slen # try step step
           changed <- ! identical(par,xn) 
           if (changed) { # have gradproj < 0 and changed at this point, so get fn value
             fval <- fn(xn, ...) # ?? try()
             nfn <- nfn + 1
             if (trace > 1) {cat(" lambda =", lambda,"  fval=", fval,"\n")}
             if (fval < fbest) badstep <- FALSE # we have a success, so exit
           } else {
              fval <- fbest # just in case, ensure fval set to fbest
              keepgoing<-FALSE # done! solveOK is TRUE, increasing lambda will not change
           }
         } # end if gradproj < 0 -- now decide if lambda to increase or not
       } # solveOK
       if (trace > 2) { cat("lambda=",lambda,"  solveOK=",solveOK," badstep=",badstep,
                 "  keepgoing=",keepgoing,"  changed=",changed,"  gradproj=",gradproj,"\n") }
       if(keepgoing) { 
         if (badstep) {# INCREASE lambda
           lambda <- max(lambdamin,lambda)*laminc # increase lambda
          if (trace > 1) cat("Increased lambda to ",lambda,"\n")
         }
         else { # success -- must have good step
           if (trace > 0) cat(" Success at lambda =", lambda,"  fval=", fval,"\n")
           par<-xn # save parameters, but fval --> fbest at top of cycle
           fbest <- fval # just in case
           if (bounds) { ## Reactivate constraints? 
             for (i in 1:npar) {
               if (bdmsk[i] == 1) { # only interested in free parameters
                 if (is.finite(lower[i])) { # JN091020 -- use abs in case bounds negative
                   if ((par[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) {# near / < lower bd
                     if (trace > 2) cat("(re)activate lower bd ", i, " at ", lower[i], "\n")
                     bdmsk[i] <- -3
                   }  # end lower bd reactivate
                 }
                 if (is.finite(upper[i])) {# JN091020 -- use abs in case bounds negative
                   if ((upper[i] - par[i]) < ceps * (abs(upper[i]) + 1)) {# near / > upper bd
                     if (trace > 2) cat("(re)activate upper bd ", i, " at ", upper[i], "\n")
                     bdmsk[i] <- -1
                   } # end upper bd reactivate
                 }   
               }  # end test on free params
             }  # end loop reactivate constraints
           }  # end if bounds
         } # end good step
       } # end keepgoing
     } # end while badstep (and keepgoing) <====== !!
    if (control$watch) { tmp <- readline("end iteration") }
  } # end main while loop
  if (itn >= control$maxit) {
     msg <- "snewtonm: Too many iterations!"
     if(trace > 0) cat(msg,"\n")
     convcode <- 1
  } else { msg <- "snewtonm: Normal exit" }
  out<-NULL
  out$par<-xn
  out$value<-fn(xn,...)
  out$grad<-grd
  out$hessian<-H
  out$counts <- list(niter=itn, nfn=nfn, ngr=ngr, nhess=nhess)
  out$convcode <- convcode
  out$message <- msg
  out 
} # end snewtonm
