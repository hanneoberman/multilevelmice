# HessQuality.R -- test initial and snewtonm final Hessians on funconstrain problems
sfname <- readline("Sink name=")
if (length(sfname)>0) sink(sfname, split=TRUE)
library(funconstrain)  # get the functions
# now translate for use with optimx
ffn <- function(n=NULL, fnum=NULL){
  # return list with tfn=function, tgr=gradient given fn number and n
  if (is.null(fnum)) stop("ffn needs a function number fnum")
  if ((fnum < 1) || (fnum > 35)) stop("fnum must be in [1, 35]")
  # select function
  funnam <- c("rosen", "freud_roth", "powell_bs", "brown_bs", "beale", 
              "jenn_samp", "helical", "bard", "gauss", "meyer", "gulf",
              "box_3d", "powell_s", "wood", "kow_osb", "brown_den", 
              "osborne_1", "biggs_exp6", "osborne_2", "watson", "ex_rosen", 
              "ex_powell", "penalty_1", "penalty_2", "var_dim", "trigon", 
              "brown_al", "disc_bv", "disc_ie", "broyden_tri", "broyden_band", 
              "linfun_fr", "linfun_r1", "linfun_r1z", "chebyquad")
  #  print(str(funnam))
  fname <- funnam[as.integer(fnum)]
  while (fnum %in% 1:35) {
    ameth <- optimx::ctrldefault(2)$allmeth
    mm <- 0 # in case m value needed
    if (fnum == 1) {
      n <- 2 # fixed
      mm <- 2
      tt <- rosen()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 2) {
      n <- 2 # fixed
      mm <- 2
      tt <- freud_roth()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 3) {
      n <- 2 # fixed
      mm <- 2
      tt <- powell_bs()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 4) {
      n <- 2 # fixed
      mm <- 3
      tt <- brown_bs()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 5) {
      n <- 2 # fixed
      mm <- 3
      tt <- beale()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 6) {
      n <- 2 # fixed
      mm <- 10
      tt <- jenn_samp()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 7) {
      n <- 3 # fixed
      tt <- helical()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 8) {
      n <- 3 # fixed
      mm <- 15
      tt <- bard()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 9) {
      n <- 3 # fixed
      mm <- 15
      tt <- gauss()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 10) {
      n <- 3 # fixed
      m <- 16 # ?? how to return
      tt <- meyer()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 11) {
      n <- 3
      mm <- 99
      tt <- gulf()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 12) {
      n <- 3
      mm <- 20
      tt <- box_3d()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 13) {
      n <- 4
      tt <- powell_s()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 14) {
      n <- 4
      tt <- wood()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 15) {
      mm <- 11
      n <- 4
      tt <- kow_osb()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 16) {
      mm <- 20
      n <- 4
      tt <- brown_den()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 17) {
      mm <- 33
      n <- 5
      tt <- osborne_1()
      ameth<-ameth[-which(ameth=="L-BFGS-B")] # remove L-BFGS-B from this case
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      lo[4] <- 0
      lo[5] <- 0
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 18) {
      mm <- 20
      n <- 6
      tt <- biggs_exp6()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 19) {
      mm <- 65
      n <- 11
      tt <- osborne_2()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 20) {
      n <-8
      mm <- 31
      tt <- watson()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 21) {
      n <- 10
      tt <- ex_rosen()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 22) {
      n <- 20
      tt <- ex_powell()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 23) {
      n <- 10
      mm <- n + 1
      tt <- penalty_1()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 24) {
      n <- 10
      mm <- n + 1
      tt <- penalty_2()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 25) {
      n <- 6
      mm <- n + 2
      tt <- var_dim()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 26) {
      n <- 8
      tt <- trigon()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 27) {
      n <- 8
      mm <- n
      tt <- brown_al()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 28) {
      n <- 6
      mm <- n
      tt <- disc_bv()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 29) {
      n <- 8
      mm <- n
      tt <- disc_ie()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 30) {
      n <- 8
      mm <- n
      tt <- broyden_tri()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 31) {
      n <- 8
      mm <- n
      tt <- broyden_band()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 32) {
      mm <- 10
      n <- 8
      tt <- linfun_fr()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 33) {
      mm <- 10
      n <- 8
      tt <- linfun_r1()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 34) {
      mm <- 10
      n <- 8
      tt <- linfun_r1z()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
    
    if (fnum == 35) {
      n <- 8
      m <- n
      tt <- chebyquad()
      if (is.function(tt$x0)) {
        xx0<-tt$x0(n)
      }
      else xx0 <- tt$x0
      lo <- rep((min(xx0)-0.1), n)
      up <- rep((max(xx0)+0.1), n)
      break }
  }    
  # ?? bounds are a trial only
  mask <- c(n)
  if (n > 3) mask <- c(1, mask)
  val <- list(npar = n, fffn=tt$fn, ffgr=tt$gr, ffhe=tt$he, x0=xx0, lo=lo, up=up, 
              mask=mask, fname=fname, ameth=ameth)
  #  cat("val:"); print(val); tmp<-readline('exit ffn')
  val
} # end ffn

library(optimx)
library(numDeriv) # probably comes in with optimx
# library(pracma) # uncomment to use jacobian() and hessian() from pracma
iprob <- as.numeric(readline("Prob #(0 for all):"))
if (iprob > 0) {
  miniprob<-iprob
  maxiprob<-iprob
} else {
  miniprob<-1
  maxiprob<-35
}
ctrl <- ctrldefault(2) # to get methods
allmeth <- ctrl$allmeth 

# bded <- readline("Bounds (L, U, B, blank for none):")
# bded <- strtrim(bded,1) # trim to 1 character
# if (bded == " ") bded<-"" # caution
# 
for (iprob in miniprob:maxiprob) {
    cat("iprob=",iprob,"\n")
    tfun <- ffn(fnum=iprob)
    # print(tfun)
    cat("Problem:", tfun$fname,"\n")
    x0 <- tfun$x0
    tfn <- tfun$fffn
    tgr <- tfun$ffgr
    the <- tfun$ffhe
    # ?? masking?
    lo<-NULL
    up<-NULL
    # if ((bded=="L") || (bded=="B")) lo <- tfun$lo
    # if ((bded=="U") || (bded=="B")) up <- tfun$up
    # # cat("about to call opm\n")

    soln <- snewtonm(x0, tfn, tgr, hess=the, lower=lo, upper=up, control=list(trace=0))
    xs <- soln$par

    f0<-tfn(x0)
    fs<-tfn(xs)
    g0<-tgr(x0)
    gs<-tgr(xs)
    mg0<-max(abs(g0))
    mgs<-max(abs(gs))
    cat("  f0\t\t  fs\t\t  mg0\t\t  mgs\n")
    cat(f0,"\t\t",fs,"\t\t",mg0,"\t\t",mgs,"\n\n")

    HH0 <- the(x0)
    Hn0 <- hessian(tfn, x0)
    HJ0 <- jacobian(tgr, x0)
    HHs <- the(xs)
    Hns <- hessian(tfn, xs)
    HJs <- jacobian(tgr, xs)

    cat("isSymmetric:\t","HH0\t\t Hn0\t\t HJ0\t\t HHs\t\t Hns\t\t HJs\n")
   options(digits=6)

   ttol<-100*.Machine$double.eps # default
   cat(ttol,":\t",isSymmetric(HH0),"\t\t",
       isSymmetric(Hn0, tol=ttol),"\t\t",
       isSymmetric(HJ0, tol=ttol),"\t\t",
       isSymmetric(HHs, tol=ttol),"\t\t",
       isSymmetric(Hns, tol=ttol),"\t\t",
       isSymmetric(HJs, tol=ttol),"\n")
   
       ttol<-(abs(f0)+10000)*.Machine$double.eps # 
       cat(ttol,":\t",isSymmetric(HH0),"\t\t",
           isSymmetric(Hn0, tol=ttol),"\t\t",
           isSymmetric(HJ0, tol=ttol),"\t\t",
           isSymmetric(HHs, tol=ttol),"\t\t",
           isSymmetric(Hns, tol=ttol),"\t\t",
           isSymmetric(HJs, tol=ttol),"\n")

       ttol<-(abs(f0)+1)*1e-8
       cat(ttol,":\t",isSymmetric(HH0),"\t\t",
               isSymmetric(Hn0, tol=ttol),"\t\t",
               isSymmetric(HJ0, tol=ttol),"\t\t",
               isSymmetric(HHs, tol=ttol),"\t\t",
               isSymmetric(Hns, tol=ttol),"\t\t",
               isSymmetric(HJs, tol=ttol),"\n")
       cat("       \t",
           max(abs(HH0-t(HH0))),"\t\t",
	         max(abs(Hn0-t(Hn0))),"\t\t",
           max(abs(HJ0-t(HJ0))),"\t\t",
           max(abs(HHs-t(HHs))),"\t\t",
           max(abs(Hns-t(Hns))),"\t\t",
           max(abs(HJs-t(HJs))),"\n")

    tmp <- readline("continue")

    
} # end problem loop

if (length(sfname) > 0) sink()
# End program
