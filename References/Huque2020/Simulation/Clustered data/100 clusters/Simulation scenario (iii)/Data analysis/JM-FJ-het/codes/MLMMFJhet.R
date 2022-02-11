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
lowerbound=(myindex-1) * 10  + 1
upperbound=(myindex) * 10

for (i in lowerbound:upperbound)
{
  dataname <- paste0("data", i,".dta")
  datajomo<-read.dta(dataname)
  # Setup imputation model
  fm2<-bmiz+QoL+seifa+sep~(1|ID)
  jomoimp<-jomoImpute(data=datajomo,formula=fm2,random.L1 ="full",m=30,seed=12345)
  jimpList<-mitmlComplete(jomoimp,print="all")
  ## Analysis model 2
  
  imp2.fitMM <- lapply(jimpList, FUN=function(x){
    lmer(QoL~bmiz+seifa+(1|ID)+(0+bmiz|ID),data=x)
  })
  
  Jomo2l.fit<-testEstimates(imp2.fitMM,var.comp = TRUE)
  beta.Jomo2l<-Jomo2l.fit$estimates[,1]
  se.Jomo2l<-Jomo2l.fit$estimates[,2]
  varcom.Jomo2l<-c(sqrt(Jomo2l.fit$var.comp[1,1]),sqrt(Jomo2l.fit$var.comp[2,1]),sqrt(Jomo2l.fit$var.comp[3,1]))
  
  pm3<-0
  pm3[(Jomo2l.fit$estimates[1,1]+qt(0.025, Jomo2l.fit$estimates[1,4])*Jomo2l.fit$estimates[1,2]) < 0.50 &  (Jomo2l.fit$estimates[1,1]+qt(0.975, Jomo2l.fit$estimates[1,4])*Jomo2l.fit$estimates[1,2]) > 0.50]<-1
  pm1<-0
  pm1[(Jomo2l.fit$estimates[2,1]+qt(0.025, Jomo2l.fit$estimates[2,4])*Jomo2l.fit$estimates[2,2]) < -0.20  &  (Jomo2l.fit$estimates[2,1]+qt(0.975, Jomo2l.fit$estimates[2,4])*Jomo2l.fit$estimates[2,2]) > -0.20]<-1
  pm2<-0
  pm2[(Jomo2l.fit$estimates[3,1]+qt(0.025, Jomo2l.fit$estimates[3,4])*Jomo2l.fit$estimates[3,2]) < 0.25 &  (Jomo2l.fit$estimates[3,1]+qt(0.975, Jomo2l.fit$estimates[3,4])*Jomo2l.fit$estimates[3,2]) > 0.25]<-1
  
  output<-as.data.frame(t(c(dataname,beta.Jomo2l,varcom.Jomo2l,se.Jomo2l,pm1=pm1,pm2=pm2,pm3=pm3)))
  write.table(output, "resultFJhet.csv",append = TRUE,sep=",",col.names=FALSE)
}
