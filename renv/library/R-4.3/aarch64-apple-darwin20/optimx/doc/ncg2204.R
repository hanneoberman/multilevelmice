ncg <- function(par, fn, gr, lower=NULL, upper=NULL, bdmsk = NULL, control = list(), ...) {

## Feb 15 -- working but inefficient for exrosen
    ## An R version of the conjugate gradient minimization using the Dai-Yuan ideas
    # ncg.R 20220212 JN
    # control defaults -- idea from spg
    badbd <- function(x, lo, up){
         val <- any(x + control$offset < lo + control$offset) | 
                any(x + control$offset > up + control$offset)
         if (val) {
           cat("BAD! x: "); print(x)
           cat("lower : "); print(lo)
           cat("upper : "); print(up)
         }
         val          
    }
#    print(control)
#    stepredn <- 0.15  # Step reduction in line search
#    stepredn <- 0.25  # Step reduction in line search
      # Step reduction in line search
#    ctrl <- list(maxit = 500, maximize = FALSE, trace = 0, eps = 1e-07, 
#        stepredn = 0.2, dowarn = TRUE, tol=0)
    npar<-length(par)
    ctrl <- ctrldefault(npar)
#    print(ctrl)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$tol == 0) tol <- npar * (npar * .Machine$double.eps)  
       # for gradient test.  Note -- integer overflow if n*n*d.eps
    else tol<-ctrl$tol
    maxit <- ctrl$maxit  # limit on function evaluations
#    cat("maxit = ",maxit,"\n")
    maximize <- ctrl$maximize  # TRUE to maximize the function
    trace <- ctrl$trace  # 0 for no output, >0 for output (bigger => more output)
    stepredn <- ctrl$stepredn
#    cat("stepredn set to ",stepredn,"\n")
    if (trace > 2) cat("trace = ", trace, "\n")
    eps <- ctrl$eps
    fargs <- list(...)  # the ... arguments that are extra function / gradient data
    grNULL <- is.null(gr)
    dowarn <- ctrl$dowarn  #
    #############################################
    if (maximize) stop("rcg does NOT maximize")
    if (is.null(gr)) stop("A gradient calculation (analytic or numerical) MUST be provided for ncg") 
    if ( is.character(gr) ) {
       # Convert string to function call, assuming it is a numerical gradient function
       mygr<-function(par=par, userfn=fn, ...){
           do.call(gr, list(par, userfn, ...))
       }
    } else { mygr<-gr }
   ############# end test gr ####################
    ## Set working parameters (See CNM Alg 22)
    if (trace > 0) cat("ncg -- J C Nash 2022 - bounds constraint version of CG\n")
    bvec <- par  # copy the parameter vector
    if (badbd(bvec, lower, upper)) stop("Initial par out of bounds")
    n <- length(bvec)  # number of elements in par vector
    maxfeval <- round(sqrt(n + 1) * maxit)  # change 091219
    ig <- 0  # count gradient evaluations
    ifn <- 1  # count function evaluations (we always make 1 try below)
    acctol <- ctrl$acctol  # acceptable point tolerance
    reltest <- ctrl$reltest  # relative equality test 
#    cat("reltest=", reltest,"\n")
    ceps <- .Machine$double.eps * reltest
    pceps <- ceps*max(abs(bvec))
    accpoint <- as.logical(FALSE)  # so far do not have an acceptable point
    cyclimit <- min(2.5 * n, 10 + sqrt(n))  #!! upper bound on when we restart CG cycle
    fargs <- list(...)  # function arguments
    if (trace > 2) {cat("Extra function arguments:");  print(fargs)}
    # set default masks if not defined
    if (is.null(bdmsk)) { bdmsk <- rep(1, n) }
    if (trace > 2) { cat("bdmsk:"); print(bdmsk)}
    # Still do checks to get nolower, noupper, bounds
    if (is.null(lower) || !any(is.finite(lower))) nolower = TRUE
    else nolower = FALSE
    if (is.null(upper) || !any(is.finite(upper))) noupper = TRUE
    else noupper = FALSE
    if (nolower && noupper && all(bdmsk == 1))  bounds = FALSE
    else bounds = TRUE
    if (trace > 2) cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, " bounds = ", bounds, "\n")
    if (nolower) lower <- rep(-Inf, n)
    if (noupper) upper <- rep(Inf, n)
    ######## check bounds and masks #############
    ## NOTE: do this inline to avoid call to external routine
    if (bounds) {
      # Make sure to expand lower and upper
      if (!nolower & (length(lower) < n)) { # have lower bounds, length must be 1 or n 
        if (length(lower) == 1) {lower <- rep(lower, n) }
        else { stop("NOT 1<length(lower)<n") }
      }  # else lower OK
      if (!noupper & (length(upper) < n)) {
        if (length(upper) == 1) { upper <- rep(upper, n) }
        else { stop("NOT 1<length(upper)<n") }
      }  # else upper OK
      # At this point, we have full bounds in play
      # This implementation as a loop, but try later to vectorize
      for (i in 1:n) { #       cat('i = ',i,'\n')
        if (bdmsk[i] == 0) { # NOTE: we do not change masked parameters, even if out of bounds
          if (!nolower) {
            if (bvec[i] < lower[i]) {
              wmsg <- paste(bvec[i], " = MASKED x[", i, "] < lower bound = ", lower[i], sep = "")
              if (dowarn) warning(wmsg)
            }
          }
          if (!noupper) {
            if (bvec[i] > upper[i]) {
              wmsg <- paste(bvec[i], " = MASKED x[", i, "] > upper bound = ", upper[i], sep = "")
              if (dowarn) warning(wmsg)
            }
          }
        }
        else { # not masked, so must be free or active constraint
          if (!nolower) {
            if (bvec[i] <= lower[i]) { # changed 090814 to ensure bdmsk is set
              wmsg <- paste("x[", i, "], set ", bvec[i], " to lower bound = ", lower[i], sep = "")
              if (dowarn && (bvec[i] != lower[i])) warning(wmsg)
              bvec[i] <- lower[i]
              bdmsk[i] <- -3  # active lower bound
            }
          }
          if (!noupper) {
            if (bvec[i] >= upper[i]) {# changed 090814 to ensure bdmsk is set
              wmsg <- paste("x[", i, "], set ", bvec[i], " to upper bound = ", upper[i], sep = "")
              if (dowarn && (bvec[i] != upper[i])) warning(wmsg)
              bvec[i] <- upper[i]
              bdmsk[i] <- -1  # active upper bound
            }
          }
        }  # end not masked
      }  # end loop for bound/mask check
    } 
    ############## end bounds check #############
    # Initial function value -- may NOT be at initial point
    #   specified by user.
    if (trace > 2) {cat("Try function at initial point:"); print(bvec) }
    f <- try(fn(bvec, ...), silent = TRUE)  # Compute the function at initial point.
    if (trace > 0) {cat("Initial function value=", f, "\n") }
    if (inherits(f,"try-error")) {
      msg <- "Initial point is infeasible."
      if (trace > 0) cat(msg, "\n")
      ans <- list(par, NA, c(ifn, 0), 2, msg, bdmsk)
      names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
      return(ans)
    }
    fmin <- f
    if (trace > 0) cat("Initial fn=", f, "\n")
    if (trace > 1) print(bvec)
    # Start the minimization process
    keepgoing <- TRUE
    msg <- "not finished"  # in case we exit somehow
    oldstep <- 0.75  #!! 2/3 #!!? WHY? formerly 0.8
    ####################################################################
    fdiff <- NA  # initially no decrease
    cycle <- 0  # !! cycle loop counter
    haveg <- FALSE
    while (keepgoing) {
      # main loop -- must remember to break out of it!!
      t <- as.vector(rep(0, n))  # zero step vector
      c <- t  # zero 'last' gradient
      while (keepgoing && (cycle < cyclimit)) {
        ## cycle loop
        cycle <- cycle + 1
        if (trace > 0) cat(ifn, " ", ig, " ", cycle, " ", fmin, "  last decrease=", fdiff, "\n")
        if (trace > 1) { print(bvec) }
        if (ifn > maxfeval) {
          msg <- paste("Too many function evaluations (> ", maxfeval, ") ", sep = "")
          if (trace > 0) cat(msg, "\n")
          ans <- list(par, fmin, c(ifn, ig), 1, msg, bdmsk)  # 1 indicates not converged in function limit
          names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
          return(ans)
        }
        if (! haveg) {
          par <- bvec  # save best parameters
          ig <- ig + 1
          if (ig > maxit) {
            msg <- paste("Too many gradient evaluations (> ", maxit, ") ", sep = "")
            if (trace > 0) cat(msg, "\n")
            ans <- list(par, fmin, c(ifn, ig), 1, msg, bdmsk)  
            # 1 indicates not converged in function or gradient limit
            names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
            return(ans)
          }
          g <- mygr(bvec, ...) # gradient
          haveg <- TRUE
        }
        if (bounds) { ## Bounds and masks adjustment of gradient
          ## first try with looping -- later try to vectorize
          if (trace > 2) { cat("bdmsk:"); print(bdmsk) }
          for (i in 1:n) {
            if ((bdmsk[i] == 0)) { g[i] <- 0 } # masked, so gradient component is zero
            else {
              if (bdmsk[i] != 1) {
                if ((bdmsk[i] + 2) * g[i] < 0) { # test for -ve gradient at upper bound, +ve at lower bound
                  g[i] <- 0  # active mask or constraint and zero gradient component
                }
                else {
                  bdmsk[i] <- 1  # freeing parameter i
                  if (trace > 1) cat("freeing parameter ", i, "\n")
                }
              }
            }
          }  # end masking loop on i
          if (trace > 2) {cat("bdmsk adj:\n"); print(bdmsk); cat("proj-g:\n"); print(g) }
        }  # end if bounds
        ## end bounds and masks adjustment of gradient
        g1 <- sum(g * (g - c))  # gradient * grad-difference
        g2 <- sum(t * (g - c))  # oldsearch * grad-difference
        gradsqr <- sum(g * g)
        if (trace > 1) { cat("Gradsqr = ", gradsqr, " g1, g2 ", g1, " ", g2, " fmin=", fmin, "\n") }
        c <- g  # save last gradient
        g3 <- 1  # !! Default to 1 to ensure it is defined -- t==0 on first cycle
#        cat(gradsqr," ",tol," ",fmin," ",reltest,"\n")
        if (gradsqr > tol * (abs(fmin) + reltest) ) {
          if (g2 > 0) {
            betaDY <- gradsqr/g2
            betaHS <- g1/g2
            g3 <- max(0, min(betaHS, betaDY))  # g3 is our new 'beta' !! Dai/Yuan 2001, (4.2)
          }
        }
        else {
          msg <- paste("Very small gradient -- gradsqr =", gradsqr, sep = " ")
          if (trace > 0) cat(msg, "\n")
          keepgoing <- FALSE  # done loops -- should we break?
          break  # to leave inner loop
        }
        if (trace > 2) cat("Betak = g3 = ", g3, "\n")
        if (g3 == 0 || cycle >= cyclimit) { # we are resetting to gradient in this case
          if (trace > 0) {
            if (cycle < cyclimit) cat("Yuan/Dai cycle reset\n")
            else cat("Cycle limit reached -- reset\n")
          }
          fdiff <- NA
          cycle <- 0 # but haveg == TRUE
          oldstep<-0.75
          break  #!!
        }
        else { # drop through if not Yuan/Dai cycle reset
          t <- t * g3 - g  # t starts at zero, later is step vector
          gradproj <- sum(t * g)  # gradient projection
          if (trace > 1) cat("Gradproj =", gradproj, "\n")
          if (bounds) { ## Adjust search direction for masks
            if (trace > 2) { cat("t:\n"); print(t) }
            t[which(bdmsk <= 0)] <- 0  # apply mask constraint
            if (trace > 2) { cat("adj-t:\n"); print(t) }
            ## end adjust search direction for masks
          }  # end if bounds
          # Why do we not check gradproj size??
          ########################################################
          ####                  Line search                   ####
          ####  2022-2-12 Use backtrack only
##          OKpoint <- FALSE # 
          accpoint <- FALSE
          f1 <- fmin # to ensure it is defined
          f <- fmin # and as large as fmin
          if (trace > 2) cat("Start linesearch with oldstep=", oldstep, "\n")
          steplength <- oldstep * 1.5  #!! try a bit bigger
          stepstrt<-steplength
          if (bounds) { # Box constraint -- adjust step length
            for (i in 1:n) { # loop on parameters -- vectorize?
              if ((bdmsk[i] == 1) && (abs(t[i]) > pceps)) { # only free params and search != 0
                if (t[i] < 0) { # going downhill. Look at lower bound
                  trystep <- (lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
                }
                else { # going uphill, check upper bound
                  trystep <- (upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
                }
                # if (trace > 3) cat("steplength, trystep:", steplength, trystep, "\n")
                steplength <- min(steplength, trystep)  # reduce as necessary
              }  # end steplength reduction
              else {t[i] <- 0} # small dirn gets set to zero just in case 
            }  # end loop on i to reduce step length
            bdlim <- (steplength < stepstrt) # TRUE if bound limits step
            stepstrt <- steplength # save this value
            if (trace > 1) cat("reset steplength (",bdlim,") = ", steplength, "\n")
             # end box constraint max step length
          }  # end if bounds
          changed <- TRUE  # Need to set so loop will start
          accpoint <- FALSE
          while ((f >= fmin) && changed && ! accpoint) {
          # may need f test to stop cycling
#          while (changed && ! accpoint) {
            bvec <- par + steplength * t
            # ?? see if we can remove this test ?? Should be able to avoid
            if (badbd(bvec, lower, upper)) stop("Top of backtrack")
#            changed <- (! identical((bvec + reltest), (par + reltest)))
            changed <- any(par+reltest != bvec+reltest)
            if (trace > 1 && changed) cat("CHANGED - ")
            if (trace > 1) {cat("trial bvec:"); print(bvec)}
            if (changed) { 
              f <- fn(bvec, ...)  # Because we need the value for linesearch, don't use try()??
              # instead preferring to fail out, which will hopefully be unlikely.
              ifn <- ifn + 1
              if (is.na(f) || (!is.finite(f))) {
                warning("ncg - undefined function")
                f <- .Machine$double.xmax
              }
              savestep<-steplength
              if (f < fmin) { f1 <- f } # Hold onto value (not needed for backtrack only??)
              else {
                steplength <- steplength * stepredn # reduce step size
                if (steplength >= savestep) changed<-FALSE
                if (trace > 0) cat("*")
              }
              accpoint <- (f1 <= fmin + gradproj * steplength * acctol) # changed MUST be TRUE or f1>...
##              OKpoint <- (f1 < fmin)
            } # changed 
          }  # end while. At this point f == f1 < fmin if changed TRUE
          if (accpoint) {
            fdiff <- fmin - f
            fmin <- f
            par <- bvec
            haveg <- FALSE # want new cycle
            if (trace > 2) { cat("new fmin=",fmin,"\n"); print(bvec)}
            oldstep <- steplength
          }
          else {
            msg <- "No acceptable point -- exit loop"
            if (trace > 0) cat("\n", msg, "\n")
            if (cycle == 1) {
              msg <- " Converged -- no progress on new CG cycle"
              if (trace > 0) cat("\n", msg, "\n")
            }
            keekpgoing <- FALSE
            break  #!! ?? where does this lead?
          }
        }  # end of test on Yuan/Dai condition
        #### End line search ####
        if (bounds) { ## Reactivate constraints? -- should check for infinite bounds!!?
          for (i in 1:n) {
            if (bdmsk[i] == 1) { # only interested in free parameters
              if (is.finite(lower[i])) {# JN091020 -- need to use abs in case bounds negative
                if ((bvec[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) { # are we near or lower than lower bd
                  if (trace > 2) cat("(re)activate lower bd ", i, " at ", lower[i], "\n")
                  bdmsk[i] <- -3
                }  # end lower bd reactivate
              }
              if (is.finite(upper[i])) { # JN091020 -- need to use abs in case bounds negative
                if ((upper[i] - bvec[i]) < ceps * (abs(upper[i]) + 1)) { # are we near or above upper bd
                  if (trace > 2) cat("(re)activate upper bd ", i, " at ", upper[i], "\n")
                  bdmsk[i] <- -1
                }  # end upper bd reactivate
              }
            }  # end test on free params
          }  # end reactivate constraints loop
        }  # end if bounds
 # ?? try step reset here
        if (oldstep < acctol) { oldstep <- acctol }  #   steplength
        if (oldstep > 1) { oldstep <- 1 }
      }  # end of inner loop (cycle)
      if (trace > 1) cat("End inner loop, cycle =", cycle, "\n")

    }  # end of outer loop
    msg <- "ncg seems to have converged"
    if (trace > 0) 
        cat(msg, "\n")
    #  par: The best set of parameters found.
    #  value: The value of 'fn' corresponding to 'par'.
    #  counts: number of calls to 'fn' and 'gr' (2 elements)
    # convergence: An integer code. '0' indicates successful
    #   convergence.
    #  message: A character string or 'NULL'.
    ans <- list(par, fmin, c(ifn, ig), 0, msg, bdmsk)
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message", "bdmsk")
    return(ans)
}  ## end of ncg
