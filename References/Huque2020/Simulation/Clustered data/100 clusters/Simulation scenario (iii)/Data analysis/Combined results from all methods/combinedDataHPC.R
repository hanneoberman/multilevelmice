setwd("~/Simulation/Clustered data/100 clusters/Simulation scenario (iii)/Data analysis")

### JM-MLMM
resultMLMM1<-as.data.frame(read.csv("JM-MLMM/Results/resultpanRISM.csv",sep=",",header = FALSE))
resultMLMM<-resultMLMM1[,3:ncol(resultMLMM1)]
names(resultMLMM)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultMLMM$Method<-rep(4,nrow(resultMLMM))
### JM-FJ
resultFJ1<-as.data.frame(read.csv("JM-FJ/Results/resultFJ.csv",sep=",",header = FALSE))
resultFJ<-resultFJ1[,3:ncol(resultFJ1)]
names(resultFJ)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultFJ$Method<-rep(5,nrow(resultFJ))
### JM-FJ-het
resultFJhet1<-as.data.frame(read.csv("JM-FJ-het/Results/resultFJhet.csv",sep=",",header = FALSE))
resultFJhet<-resultFJhet1[,3:ncol(resultFJhet1)]
names(resultFJhet)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultFJhet$Method<-rep(6,nrow(resultFJhet))

### JM-SMC
resultSMC1<-as.data.frame(read.csv("JM-SMC/Results/resultJomoSMC.csv",sep=",",header = FALSE))
resultSMC<-resultSMC1[,3:ncol(resultSMC1)]
names(resultSMC)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultSMC$Method<-rep(7,nrow(resultSMC))
### JM-SMC-het
resultSMChet1<-as.data.frame(read.csv("JM-SMC-het/Results/resultJomoSMChet.csv",sep=",",header = FALSE))
resultSMChet<-resultSMChet1[,3:ncol(resultSMChet1)]
names(resultSMChet)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultSMChet$Method<-rep(8,nrow(resultSMChet))
###FCS-LMM
resultFCS1<-as.data.frame(read.csv("FCS-LMM/Results/result2lpanRISM.csv",sep=",",header = FALSE))
resultFCS<-resultFCS1[,3:ncol(resultFCS1)]
names(resultFCS)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultFCS$Method<-rep(10,nrow(resultFCS))
###FCS-LMM-het
resultFCS2lnorm1<-as.data.frame(read.csv("FCS-LMM-het/Results/result2lnormRISM.csv",sep=",",header = FALSE))
resultFCS2lnorm<-resultFCS2lnorm1[,3:ncol(resultFCS2lnorm1)]
names(resultFCS2lnorm)<-c("cons","bmiz","seifa","varIntercept","varSlope","varResid","secons","sebmiz","seseifa","pbmiz","pSeifa","pCons")
resultFCS2lnorm$Method<-rep(11,nrow(resultFCS2lnorm))

result<-rbind(resultFCS,resultFCS2lnorm,resultFJ,resultFJhet, resultMLMM,resultSMC,resultSMChet)
write.table(result, "Combined results from all methods/Results.csv",sep=",",col.names=TRUE,row.names=FALSE)
