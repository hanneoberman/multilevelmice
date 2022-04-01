## code to prepare `countExample` dataset

library(tidyverse, warn.conflicts = FALSE)
library(MASS)
library(truncnorm)
library(magrittr)
library(reshape2)
library(fastDummies)
library(dplyr)
library(optmatch) #match
library(MatchIt)
library(WeightIt) #IPW
library(cobalt) # imbalance covariates
library(ggdag)
library(dagitty)
library(mice)
library(data.table)
library(PSweight)
library(sandwich) 
library(lmtest)

# Function to simulate MS data
simcountdata <- function(n, seed = 999,
                         beta = c(-0.5, -0.25, 0, 0.25, 0.5),
                         beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003),
                         percentiles = seq(0, 1, by = 0.2)){
  
  #' Generate simulated count data with settings based on real-world data
  #' Assume randomized treatment and independent covariates
  #'
  #' @param n sample size; integer
  #' @param seed randomization seed; integer
  #' @param beta coefficients characterizing treatment effect heterogeneity; vector of length 5
  #'             beta[1]*trt*I(high responder to DMF) +
  #'             beta[2]*trt*I(moderate responder to DMF) +
  #'             beta[3]*trt*I(neutral) +
  #'             beta[4]*trt*I(moderate responder to TERI) +
  #'             beta[5]*trt*I(high responder to TERI)
  #'             In the absence of treatment effect heterogeneity, set all beta[1:5] with the same values
  #' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
  #'               beta.x[1] (intercept)
  #'               beta.x[2]*ageatindex_centered
  #'               beta.x[3]*female
  #'               beta.x[4]*prerelapse_num
  #'               beta.x[5]*prevDMTefficacy== mediumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return data - a data frame of \code{n} rows and 8 columns of simulated data; data.frame
  
  
  set.seed(seed)
  
  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6){message("Wrong values of percentiles!")}
  
  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 10))
  colnames(ds) <- c("treatment", "ageatindex_centered", "female", "prerelapse_num",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms",
                    "postrelapse_num", "finalpostdayscount", "group")
  
  # Define X, A, and time
  ds %<>%
    mutate(trt =                 rbinom(n = n, size = 1, prob = 0.75),
           treatment =           as.factor(ifelse(trt == 1, "DMF", "TERI")),
           female =              rbinom(n = n, size = 1, prob = 0.75),
           ageatindex_centered = round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48, # rounded to integers
           prerelapse_num =      rpois(n = n, lambda = 0.44),
           prevDMTefficacy =     sample(x = c("None", "Low efficacy", "Medium and high efficacy"),
                                        size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11)),
           premedicalcost =      pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000), # rounded to 2 decimal points
           numSymptoms =         sample(x = c("0", "1", ">=2"),
                                        size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09)), # nuisance variable; do not include in the score
           finalpostdayscount =  ceiling(rgamma(n = n, shape = 0.9, scale = 500)), # rounded up to integers
           finalpostdayscount =  ifelse(finalpostdayscount > 2096, 2096, finalpostdayscount), # truncate at the max follow up day, 2096
           finalpostdayscount =  ifelse((finalpostdayscount > 2090) & (runif(1, 0, 1) < .5), 29, finalpostdayscount), # mimic the 1 month peak;  move roughly half of the large values to 29
           group =              "Simulated")
  
  # Define Y
  xmat.score <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost, ds) %>% as.matrix()
  gamma <- matrix(c(-0.33, # Intercept
                    -0.001, # Age
                    0.05, # female
                    -0.002, # prerelapse_num
                    0.33, 0.02, # Medium/high and none DMT efficacy
                    -0.0000005), nrow = 7) # premedicalcost
  ds <- ds %>% mutate(score = exp(xmat.score %*% gamma),
                      Iscore = cut(score, quantile(score, percentiles), include.lowest = T, labels = seq(1, 5)))
  
  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trt + trt*Iscore, ds) %>% as.matrix()
  
  betas <- matrix(c(beta.x[1], # Intercept
                    beta.x[2], # Age
                    beta.x[3], # female
                    beta.x[4], # prerelapse_num
                    beta.x[5], beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    # Previous beta.x: -1.54 Intercept, -0.01 Age, 0.06 female, 0.47 prerelapse_num, 0.68, 0.13 Medium/high and none DMT efficacy, 0.000003 premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1]))
  rate <-  exp(xmat.rate %*% betas)
  ds <- ds %>% mutate(postrelapse_num = rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3))
  
  return(list(data = ds, betas = betas, percentiles = percentiles, rate=rate))
}


# Generate complete data ----
generate_data<-function(beta,n,seed){
  data <- simcountdata(n = n, 
                       seed = seed,
                       beta = beta,
                       #beta=c(-0.2, -0.2, -0.2, -0.2, -0.2) 
                       #beta = c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)), 
                       # this beta is medium level of heterogeneity (see other levels in "Simulation design.pptx")
                       # other options:
                       #  none - beta = c(-0.2. -0.2, -0.2, -0.2, -0.2)
                       #  low - beta = c(log(0.7), log(0.75), log(1), log(1.05), log(1.1))
                       #  high - beta = c(log(0.3), log(0.5), log(1), log(1.1), log(1.5))
                       beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003)
  )$data %>% dplyr::rename(previous_treatment = prevDMTefficacy,
                           age = ageatindex_centered,
                           y = postrelapse_num,
                           years = finalpostdayscount,
                           previous_number_relapses = prerelapse_num,
                           previous_number_symptoms = numSymptoms,
                           previous_cost = premedicalcost) %>%
    mutate(previous_treatment = factor(previous_treatment, labels = c("drugA", "drugB", "drugC")),
           previous_number_symptoms = factor(previous_number_symptoms, labels = c("0", "1", ">=2")),
           years = years / 365.25,
           age = age + 48)%>%
    dplyr::select(age, female, previous_treatment, previous_cost, previous_number_symptoms, previous_number_relapses, treatment, y, years)
  return(data)}


# Get estimands ATE/ATT ----
get_estimand <- function(data,approach,estimand="ATE",wmethod="ps"){
  
  data<-data%>%mutate( DMF = ifelse(treatment == "DMF",1,0),
                       TERI = ifelse(treatment == "TERI",0,1))
  
  formula <-DMF ~ age + female + previous_treatment + previous_cost+  previous_number_symptoms + previous_number_relapses
  
  if(approach=="Match"){
    
    method <- ifelse(estimand=="ATE","full","nearest")
    # Step 1: apply matching
    mout <- matchit(formula, 
                    data = data,
                    family = binomial,
                    method = method,
                    caliper = 0.2,
                    std.caliper = TRUE,
                    estimand = estimand,
                    distance = "glm",
                    link = "logit",
                    replace = FALSE)
    
    # Step 2: retrieve matched sample
    if(estimand=="ATE"){
      mdata <- match.data(mout)
    }else{
      mdata<-get_matches(mout)
    }
  }else {
    mdata<-data
    W.out <- weightit(formula,data = mdata, estimand = estimand, method = wmethod)
    mdata$ipw<-W.out$weights
    
  }
  # Step 3: estimate ATT or ATE-----
  
  
  if(approach=="Match"){
    model_tx <- glm("y ~ DMF + offset(log(years))", 
                    family = poisson(link = "log"),
                    data = mdata)
    se.fit <- summary(model_tx)$coefficients[, 2]["DMF"]
  }else{
    model_tx <- glm("y ~ DMF + offset(log(years))", 
                    family = poisson(link = "log"),
                    data = mdata,
                    weights = ipw)
    se.fit<- coeftest(model_tx, vcov = vcovHC)["DMF","Std. Error"] #Robust estimate 
  }
  
  # Step 4: extract incidence rate ratio for DMF versus TERI
  za <- qnorm(1 - 0.05 / 2) # for 95% CIs
  fit <- summary(model_tx)$coefficients[, 1]["DMF"]
  lci <-(fit-za*se.fit)
  uci <-(fit+za*se.fit)
  return(c(approach=approach,estimand=estimand,fit=fit,se.fit=se.fit,lci=lci,uci=uci))
}



# Function to transform database with dummy to categorical variables----
databack<-function(data){
  data<-setDT(data)
  data[,previous_treatment:=as.factor(ifelse(previous_treatmentdrugA==1,"drugA",ifelse(previous_treatmentdrugB==1,"drugB",ifelse(is.na(previous_treatmentdrugA),NA,"drugC"))))]
  data[,previous_number_symptoms:=as.factor(ifelse(previous_number_symptoms1==1,"1",ifelse(previous_number_symptoms..2==1,">=2",ifelse(is.na(previous_number_symptoms1),NA,"0"))))]
  data[,treatment:=as.factor(ifelse(treatmentTERI,"TERI",ifelse(is.na(treatmentTERI),NA,"DMF")))]
  data[,previous_number_symptoms:=factor(previous_number_symptoms,levels=c("0","1",">=2"))]
  data[,treatmentTERI:=NULL]
  data[,previous_treatmentdrugA:=NULL]
  data[,previous_treatmentdrugB:=NULL]
  data[,previous_number_symptoms1:=NULL]
  data[,previous_number_symptoms..2:=NULL]
  return(as.data.frame(data))
}

#Get missing data ---
getmissdata <- function(data){
  
  data<-data.frame(model.matrix(~.-1,data))
  data$previous_treatmentdrugC<-NULL
  

  cost_mcar<-ampute(data, patterns=c(1,1,1,1,0,1,1,1,1,1,1),prop=0.1, mech="MCAR")$amp
  cost_mar<-ampute(data, patterns = c(1,1,1,1,0,1,1,1,1,1,1),weights=c(1,1,0,0,0,0,0,0,0,0,0), prop = 0.1, mech = "MAR")$amp #f(age,sex)
  presym_mcar<-ampute(data, patterns=c(1,1,1,1,1,0,0,1,1,1,1),prop=0.08, mech="MCAR")$amp
  prey_mcar<-ampute(data, patterns=c(1,1,1,1,1,1,1,0,1,1,1),prop=0.05, mech="MCAR")$amp
  prey_mar<-ampute(data, patterns = c(1,1,1,1,1,1,1,0,1,1,1),weights=c(0,0,0,0,0,1,1,0,0,0,0), prop = 0.05, mech = "MAR",type="LEFT")$amp #prev_y=f(previous_cost, previous symptoms)
  y_mcar<-ampute(data, patterns=c(1,1,1,1,1,1,1,1,1,0,1),prop=0.25, mech="MCAR")$amp
  y_mar<-ampute(data, patterns = c(1,1,1,1,1,1,1,1,1,0,1),weights=c(0,0,0,0,1,1,1,0,0,0,0), prop = 0.25, mech = "MAR",type="LEFT")$amp #prev_y=f(previous_cost, previous symptoms)
  y_mart<-ampute(data,patterns=c(1,1,1,1,1,1,1,1,1,0,1),weights=c(0,0,0,0,1,1,1,0,1,0,0),prop=0.25,mech="MAR")$amp #post_y=f(previous_cost, previous symptoms,treatment)
  y_mnar<-ampute(data,patterns=c(1,1,1,1,1,1,1,1,1,0,1),weights=c(0,0,0,0,1,1,1,0,0,1,0),prop=0.25,mech="MNAR")$amp #post_y=f(post_y values)
  
  #MCAR
  ampdata1<-cost_mcar
  ampdata1[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata1$previous_number_relapses<-prey_mcar$previous_number_relapses
  ampdata1$y<-y_mcar$y
  
  #MAR
  ampdata2<-cost_mar
  ampdata2[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata2$previous_number_relapses<-prey_mar$previous_number_relapses
  ampdata2$y<-y_mar$y
  
  #MAR-y(treatment)
  ampdata3<-cost_mar
  ampdata3[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata3$previous_number_relapses<-prey_mar$previous_number_relapses
  ampdata3$y<-y_mart$y
  
  #MAR-y(MNAR)
  ampdata4<-cost_mar
  ampdata4[,c("previous_number_symptoms1","previous_number_symptoms..2")]<-presym_mcar[,c("previous_number_symptoms1","previous_number_symptoms..2")]
  ampdata4$previous_number_relapses<-prey_mar$previous_number_relapses
  ampdata4$y<-y_mnar$y
  
  return(list(ampdata1=databack(ampdata1),ampdata2=databack(ampdata2),
              ampdata3=databack(ampdata3),ampdata4=databack(ampdata4)))
}

#Pool mice results-----
pool_est <- function(est,se_est){  
  est<-as.numeric(est)
  se_est<-as.numeric(se_est)
  mean<-mean(est)
  m <- length(est)
  w <-mean(se_est^2) # within variance
  b<-var(est) # between variance
  tv <- w + (1 + (1/m)) * b # total variance
  se_total <-sqrt(tv)
  r <- (1 + 1 / m) * (b / w)
  v <- (m - 1) * (1 + (1/r))^2
  t <- qt(0.975, v)
  lc<-mean-se_total*t
  uc<-mean+se_total*t
  return(c(mean,se_total,lc,uc))
}

#MICE+ balance+ estimand-----

mice_estimand<-function(miceout, approach,estimand){
  imputed.dfs <- complete(miceout, "all")
  limputed.dfs <- lapply(imputed.dfs, get_estimand,approach=approach,estimand=estimand)
  imputed.res<-do.call(rbind,limputed.dfs)
  pool_sest <- pool_est(est=imputed.res[,3],se_est=imputed.res[,4])
  return(c(approach,estimand,pool_sest))
}

# Get estimands for all scenarios----
getall<-function(mdata,miceout,estimand){
  
  
  #MICE + matchapproach="NONE","
  
  nmm1<-mice_estimand(miceout=miceout$ampdata1,approach="Match",estimand=estimand)
  nmm2<-mice_estimand(miceout=miceout$ampdata2,approach="Match",estimand=estimand)
  nmm3<-mice_estimand(miceout=miceout$ampdata3,approach="Match",estimand=estimand)
  nmm4<-mice_estimand(miceout=miceout$ampdata4,approach="Match",estimand=estimand)
  
  #MICE + ipw
  nmip1<-mice_estimand(miceout=miceout$ampdata1,approach="IPW",estimand=estimand)
  nmip2<-mice_estimand(miceout=miceout$ampdata2,approach="IPW",estimand=estimand)
  nmip3<-mice_estimand(miceout=miceout$ampdata3,approach="IPW",estimand=estimand)
  nmip4<-mice_estimand(miceout=miceout$ampdata4,approach="IPW",estimand=estimand)
  
  
  return(as.data.table(rbind(nmm1,nmm2,nmm3,nmm4,
                             nmip1,nmip2,nmip3,nmip4)))
}





# Procedure-----

#1. Generate full datasets according to heterogeneity in treatment specification ----
data_noHTE <- generate_data(n=10000,beta=c(-0.2, -0.2, -0.2, -0.2, -0.2),seed=1234) # no treatment effect
data_medHTE <- generate_data(n=10000,beta=c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)),seed=1234) #medium treatment effect


#2. Generate missingness in full data
#2.1 MCAR
#2.2 Cost MAR=f(age,female),(prev_num and post_num(Y) =f(previous_cost, previous symptoms) MAR at the same time) No...
#2.2 Cost MAR=f(age,female), (prev_num and post_num(Y) MAR at different observations same time) 
#2.3 As 2.2 but post_num(Y) is dependable on treament post_num(Y) =f(previous_cost, previous symptoms,treatment)
#2.4.As 2.3 but post_num(Y) is MNAR ampdata_y_mnar=f(post_Y)

mdata_noHTE<-getmissdata(data_noHTE)
mdata_medHTE<-getmissdata(data_medHTE)


#3. Estimate treatment effect estimand

#3.1. NoHTE
true.noHTE.ATE.Match<-get_estimand(data=data_noHTE,approach="Match",estimand="ATE") #Real data
true.noHTE.ATE.IPW<-get_estimand(data=data_noHTE,approach="IPW",estimand="ATE") #Real data
true.noHTE.ATT.Match<-get_estimand(data=data_noHTE,approach="Match",estimand="ATT") 
true.noHTE.ATT.IPW<-get_estimand(data=data_noHTE,approach="IPW",estimand="ATT") 
true.noHTE<-as.data.table(rbind(true.noHTE.ATE.Match,true.noHTE.ATE.IPW,true.noHTE.ATT.Match,true.noHTE.ATT.IPW))
true.noHTE[,HTE:="noHTE"]
miceout_noHTE<-lapply(mdata_noHTE,mice,m=5)
est.noHTE.ATE<-getall(mdata=mdata_noHTE,miceout=miceout_noHTE,estimand="ATE")
est.noHTE.ATT<-getall(mdata=mdata_noHTE,miceout=miceout_noHTE,estimand="ATT")

#3.1. medHTE
true.medHTE.ATE.Match<-get_estimand(data=data_medHTE,approach="Match",estimand="ATE") #Real data
true.medHTE.ATE.IPW<-get_estimand(data=data_medHTE,approach="IPW",estimand="ATE") #Real data
true.medHTE.ATT.Match<-get_estimand(data=data_medHTE,approach="Match",estimand="ATT") 
true.medHTE.ATT.IPW<-get_estimand(data=data_medHTE,approach="IPW",estimand="ATT") 
true.medHTE<-as.data.table(rbind(true.medHTE.ATE.Match,true.medHTE.ATE.IPW,true.medHTE.ATT.Match,true.medHTE.ATT.IPW))
true.medHTE[,HTE:="medHTE"]
miceout_medHTE<-lapply(mdata_medHTE,mice,m=5)
est.medHTE.ATE<-getall(mdata=mdata_medHTE,miceout=miceout_medHTE,estimand="ATE")
est.medHTE.ATT<-getall(mdata=mdata_medHTE,miceout=miceout_medHTE,estimand="ATT")
save.image(file='MS_missing10000.RData')

load(file='MS_missing.RData')


### Added by Hanne
saveRDS(miceout_noHTE, "MS_imps.RDS")
