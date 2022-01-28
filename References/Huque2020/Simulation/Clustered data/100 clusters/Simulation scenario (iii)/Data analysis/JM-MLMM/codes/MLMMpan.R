library(mitml)
library(reshape2)
library(mgcv)
library(lme4)
library(plyr)
library(mice)
library(foreign)
library(jomo)

strindex=Sys.getenv("PBS_ARRAY_INDEX")
myindex=strtoi(strindex)
lowerbound=(myindex-1) * 10  + 1
upperbound=(myindex) * 10

for (i in lowerbound:upperbound)
{
  dataname <- paste0("data", i,".dta")
  datapan<-read.dta(dataname)
  # Setup imputation model
  fm1<-bmiz~1+QoL+seifa+sep+(1+QoL|ID)
  imp<-panImpute(data=datapan,formula=fm1,m=30,seed=12345)
  #summary(imp)
  impList<-mitmlComplete(imp,print="all")
  fit.lmm<-with(impList,lmer(QoL~bmiz+seifa+(1|ID) +(0+bmiz|ID)))
  
  lmm.fit<-testEstimates(fit.lmm,var.comp = TRUE)
  
  beta.lmm<-lmm.fit$estimates[,1]
  se.lmm<-lmm.fit$estimates[,2]
  varcom.lmm<-c(sqrt(lmm.fit$var.comp[1,1]),sqrt(lmm.fit$var.comp[2,1]),sqrt(lmm.fit$var.comp[3,1]))
  
  pm3<-0
  pm3[(lmm.fit$estimates[1,1]+qt(0.025, lmm.fit$estimates[1,4])*lmm.fit$estimates[1,2]) < 0.50 &  (lmm.fit$estimates[1,1]+qt(0.975, lmm.fit$estimates[1,4])*lmm.fit$estimates[1,2]) > 0.50]<-1
  pm1<-0
  pm1[(lmm.fit$estimates[2,1]+qt(0.025, lmm.fit$estimates[2,4])*lmm.fit$estimates[2,2]) < -0.20  &  (lmm.fit$estimates[2,1]+qt(0.975, lmm.fit$estimates[2,4])*lmm.fit$estimates[2,2]) > -0.20]<-1
  pm2<-0
  pm2[(lmm.fit$estimates[3,1]+qt(0.025, lmm.fit$estimates[3,4])*lmm.fit$estimates[3,2]) < 0.25 &  (lmm.fit$estimates[3,1]+qt(0.975, lmm.fit$estimates[3,4])*lmm.fit$estimates[3,2]) > 0.25]<-1
  
  output<-as.data.frame(t(c(dataname,beta.lmm,varcom.lmm,se.lmm,pm1=pm1,pm2=pm2,pm3=pm3)))
  write.table(output, "resultpanRISM.csv",append = TRUE,sep=",",col.names=FALSE)
}
