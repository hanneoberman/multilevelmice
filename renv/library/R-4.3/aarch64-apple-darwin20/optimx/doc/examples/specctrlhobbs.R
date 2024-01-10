# rm(list=ls()) # comment out this line if you do not want the workspace cleared
##  author: John C. Nash
require(optimx)
sessionInfo()
## Optimization test function HOBBS
## ?? refs (put in .doc??)
## Nash and Walker-Smith (1987, 1989) ...
hobbs.f<- function(x){ # # Hobbs weeds problem -- function
    if (abs(120*x[3]) > 500) { # check computability
       fbad<-.Machine$double.xmax
       return(fbad)
    }
    res<-hobbs.res(x)
    f<-sum(res*res)
}
hobbs.res<-function(x){ # Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    t<-1:12
    if(abs(12*x[3])>50) {
       res<-rep(Inf,12)
    } else {
       res<-100*x[1]/(1+10*x[2]*exp(-0.1*x[3]*t)) - y
    }
}
x0 <- c(1,1,1)
cat("specctrlhobbs.R -- show special controls passed through optimr to bobyqa\n")
cat("default rhobeg and rhoend in bobyqa\n")
tdef <- optimr(x0, hobbs.f, method="bobyqa", control=list(trace=1))
proptimr(tdef)
tdef1 <- optimr(x0, hobbs.f, method="bobyqa", control=list(rhobeg=-1, trace=1))
proptimr(tdef1)
tbdef <- minqa::bobyqa(par=x0, fn=hobbs.f, control=list(iprint=1))
print(tbdef)

# optimr choices -- depend on elements of bounds in problem
tsp0 <- optimr(x0, hobbs.f, method="bobyqa", control=list(rhobeg=0, trace=1))
proptimr(tsp0)

tsp2 <- optimr(x0, hobbs.f, method="bobyqa", control=list(rhobeg=1, rhoend=1e-5, trace=1))
proptimr(tsp2)
tbsp2 <- minqa::bobyqa(par=x0, fn=hobbs.f, control=list(rhobeg=1, rhoend=1e-5, iprint=1))
print(tbsp2)

tsp3 <- optimr(x0, hobbs.f, method="bobyqa", control=list(rhobeg=10, rhoend=1e-4, trace=1))
proptimr(tsp3)
tbsp3 <- minqa::bobyqa(par=x0, fn=hobbs.f, control=list(rhobeg=10, rhoend=1e-4, iprint=1))
print(tbsp3)


