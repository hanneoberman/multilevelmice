library(data.table)
library(truncnorm)
library(nlme)
library(devtools)
library(dplyr)
library(ggplot2)
library(mice)

#0. Load imputation method and original dataset ----
# Source imputation function
source_url("https://github.com/johamunoz/Statsmed_Heckman/blob/main/4.Codes/mice.impute.2l.heckman.R?raw=TRUE")
# Dataset obtained from "https://www.kaggle.com/datasets/fabinmndez/obesitydata?select=ObesityDataSet_raw_and_data_sinthetic.csv"
# and it is described here https://www.sciencedirect.com/science/article/pii/S2352340919306985?via%3Dihub

# data<-as.data.table(read.csv2('./Data/obesity.csv', stringsAsFactors=TRUE, fileEncoding="latin1"))
data <- as.data.table(read.csv('./Data/obesity_data.csv', stringsAsFactors=FALSE))
#1. Modify dataset in order to get a variable weight with MNAR missing mechanism ----
# Generate countries where survey was taken
data<-data[,c("Weight","Gender","Age","Height","FAVC")] # Select some of the variables of the original dataset
data <- data %>% mutate(across(c("Weight", "Age", "Height"), ~as.numeric(.x)))
data$Cluster<-c(rep(1,500),rep(2,500),rep(3,500),rep(4,500),rep(5,111)) # in the imputation method is necessary cluster is included as numeric
summary(data)
#  Generation of the exclusion restriction (ERV): The time variable describes the time it takes for a person to answer the weight question.
data$Time<-rtruncnorm(n=nrow(data), a = 1, b = 10, mean = 5, sd = 2) 
data[,Gender:=ifelse(Gender=="Female",1,0)]
data[,FAVC:=ifelse(FAVC=="yes",1,0)] 
MAR_model<- lme(Weight~Gender+Age+Height+FAVC,random = ~1|Cluster, data = data) 
# Generate selection and outcome variable (weight)
data[,XOBO:=predict(MAR_model)]
data[,XSBS:=15-0.7*Gender-0.06*Age-2*Time] # Simulated selection model

# Simulate bivariate normal correlated errors
rho=-0.6 # We assume that non-responders are more likely to have a high weight.
d <- diag(2)
d[1,1]<-MAR_model$sigma^2
d[2,1] <- d[1,2] <- rho*MAR_model$sigma
e_norm <- mvtnorm::rmvnorm(n =nrow(data), mean = rep(0,2), sigma = d)

# Generate latent outcome and selection variables
data[,Weight.star:=XOBO+e_norm[,1]]
data[,ry.star:=XSBS+e_norm[,2]]

# Generate observed variables
data[,ry:=ifelse(ry.star>0&Cluster!=5,1,0)] #Include systematic missingness on cluster 5
data[,Weight :=ifelse(ry==1,Weight.star,NA)] 
summary(data)

# Datafinal 
data[ ,c("XOBO","XSBS","Weight.star","ry.star","ry") := NULL]


#2. Descriptive analysis simulated dataset ----

# Count missingness per group
dataNA<-data[, .(nNA = sum(is.na(Weight)),n=.N), by = Cluster]
dataNA[, propNA:=nNA/n]

#Plot weight
ggplot(data, aes(x = Weight, group=as.factor(Cluster), color = as.factor(Cluster),fill= as.factor(Cluster))) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.1)


#3. Imputation model
# Data passed into the mice function has to specify the cluster variable as a numeric variable, and in case that incomplete
# variable is a binary variable, it should specified as factor variable with 2 levels. 

ini <- mice(data, maxit = 0)
meth <- ini$method
#It is necessary to specify for the MNAR missing variable the method "2l.heckman" in the mice methods vector.     
meth["Weight"] <- "2l.heckman" 

#Furthermore, in the prediction matrix, the group or cluster variable should be specified as "-2", all predictor variables
#belonging to the selection and outcome as "1", the exclusion restrictions or predictor variables that are only included in the selection equation 
#as "-3" and those that are only included in the outcome equation as "-4".
pred <- ini$predictorMatrix
pred[,"Cluster"] <- -2 # Cluster variable
pred["Weight","Time"]  <- -3 # Exclusion restriction
pred["Weight",c("Height","FAVC")]  <- -4 # Exclusion restriction

# Imputation of weight variable without pmm 
imp <- mice(data = data, meth = meth, pred = pred, maxit = 1, m = 5, seed = 1)
densityplot(imp)


# Imputation of weight variable with (predictive mean matching) pmm: we can also use the ppm approach, by setting pmm= TRUE and providing a vector of donnors ypmm,
# for instance in this case the vector or donnors are weight from 25 to 200 kilograms and the imputed values are given in this range of values.

imp_pmm <- mice(data = data, meth = meth, pred = pred, maxit = 1, m = 5, seed = 1,pmm=TRUE,ypmm=seq(25,200,0.1))

densityplot(imp_pmm)

obesity <- data %>% janitor::clean_names() %>% mutate(weight = round(weight, 2), rt = round(time, 2), .keep = "unused")
obesity$bmi = round(b$weight / (b$height)^2, 2)
obesity = obesity[, c("cluster", "gender", "age", "height", "weight", "bmi", "rt")]
save(obesity, file = "./Data/obesity.Rdata")
