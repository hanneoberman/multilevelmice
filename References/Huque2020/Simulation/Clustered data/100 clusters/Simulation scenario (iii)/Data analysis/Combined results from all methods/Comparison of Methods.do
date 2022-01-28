cd "~/Simulation/Clustered data/100 clusters/Simulation scenario (iii)/Data analysis"
use "available data/Results/Available data.dta",clear
renvars p1 p2 p3/pbmiz pSeifa pCons
gen Method=2 
saveold "Combined results from all methods/combinedStata.dta",replace 

import delimited "Combined results from all methods\Results.csv", clear case(preserve)  numericcols(_all)
duplicates drop
byso Method:sample 1000,count
saveold "Combined results from all methods/combinedR.dta",replace 

 
use "Combined results from all methods/combinedStata.dta",clear
append using "Combined results from all methods/combinedR.dta"
label define meth 1 "Full data" 2 "Available data" 3 "JM-MVN" 4 "JM-MLMM" 5 "JM-FJ" 6 "JM-FJ-het" 7 "JM-SMC" 8 "JM-SMC-het" 9 "FCS-standard"  10 "FCS-LMM" 11 "FCS-LMM-het"   
label value Method meth
gen bCons=(0.50-cons)/0.5 
gen bbmiz=(-0.20-bmiz)/0.2
gen bSeifa=(0.25-seifa)/0.25
gen bInt=(0.4-varIntercept)/0.4
gen bSlop=(0.2-varSlope)/0.2
gen bresid=(0.9-varResid)/0.9
sort Method
saveold "Combined results from all methods/Coefbias.dta",replace version(12)

***Finaly the final results for paper
byso Method: su 

