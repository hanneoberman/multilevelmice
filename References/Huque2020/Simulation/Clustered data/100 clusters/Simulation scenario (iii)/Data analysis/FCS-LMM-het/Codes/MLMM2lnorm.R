library(mitml)
library(reshape2)
library(mgcv)
library(lme4)
library(plyr)
#install.packages("mice")
library(mice)
library(foreign)
library(jomo)

strindex=Sys.getenv("PBS_ARRAY_INDEX")
myindex=strtoi(strindex)
lowerbound=(myindex-1) * 5  + 1
upperbound=(myindex) * 5

for (i in lowerbound:upperbound)
{
  dataname <- paste0("data", i,".dta")
  dataFCS2<-read.dta(dataname)
  ini2<-mice(dataFCS2,maxit=0)
  pred<-ini2$predictorMatrix
  pred["bmiz",]<-c(-2,2,2,0,2)
  #pred["QoL",]<-c(-2,2,2,2,0)
  meth2<-ini2$method
  meth2[c("bmiz")]<-"2l.norm"
  
  imp2<-mice(dataFCS2,m=30,meth=meth2,pred=pred)
  
  
  impData3L <- as.list(1:30)
  for(i in 1:30){
    impData3L[[i]] <- complete(imp2, action=i)
  }
  
  imp2.fitMM <- lapply(impData3L, FUN=function(x){
    lmer(QoL~bmiz+seifa+(1|ID)+(0+bmiz|ID),data=x)
  })
  
  fcs2.fitMM<-testEstimates(imp2.fitMM,var.comp = TRUE)
  
  
  beta.fcs2MM<-fcs2.fitMM$estimates[,1]
  se.fcs2MM<-fcs2.fitMM$estimates[,2]
  varcom.fcs2<-c(sqrt(fcs2.fitMM$var.comp[1,1]),sqrt(fcs2.fitMM$var.comp[2,1]),sqrt(fcs2.fitMM$var.comp[3,1]))
  pm3<-0
  pm3[(fcs2.fitMM$estimates[1,1]+qt(0.025, fcs2.fitMM$estimates[1,4])*fcs2.fitMM$estimates[1,2]) < 0.50 &  (fcs2.fitMM$estimates[1,1]+qt(0.975, fcs2.fitMM$estimates[1,4])*fcs2.fitMM$estimates[1,2]) > 0.50]<-1
  pm1<-0
  pm1[(fcs2.fitMM$estimates[2,1]+qt(0.025, fcs2.fitMM$estimates[2,4])*fcs2.fitMM$estimates[2,2]) < -0.20  &  (fcs2.fitMM$estimates[2,1]+qt(0.975, fcs2.fitMM$estimates[2,4])*fcs2.fitMM$estimates[2,2]) > -0.20]<-1
  pm2<-0
  pm2[(fcs2.fitMM$estimates[3,1]+qt(0.025, fcs2.fitMM$estimates[3,4])*fcs2.fitMM$estimates[3,2]) < 0.25 &  (fcs2.fitMM$estimates[3,1]+qt(0.975, fcs2.fitMM$estimates[3,4])*fcs2.fitMM$estimates[3,2]) > 0.25]<-1
  
  output<-as.data.frame(t(c(dataname,beta.fcs2MM,varcom.fcs2,se.fcs2MM,pm1=pm1,pm2=pm2,pm3=pm3)))
  write.table(output, "result2lnormRISM.csv",append = TRUE,sep=",",col.names=FALSE)
}
