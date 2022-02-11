cd "C:\Users\h.huque\Documents\MI BiomJ\Codes for bimj.201900051\Simulation\Clustered data\100 clusters\Simulation scenario (iii)\Data simulation\Data"
use data1.dta,clear
 mixed QoL bmiz seifa ||ID:bmiz,stddev iterate(20) 
matrix define t=r(table)
matrix list t
matrix define coef=e(b)
matrix define var=vecdiag(e(V))
matrix list coef
matrix list var
drop ID seifa sep bmiz QoL
svmat coef 
svmat var
gen varIntercept=exp(coef5) 
gen varSlope=exp(coef4)  
gen varResid= exp(coef6)
gen secons=sqrt(var3) 
gen sebmiz=sqrt(var1) 
gen seseifa=sqrt(var2) 
gen p1=0
replace p1=1 if t[5,1]< -0.20 & t[6,1]> -0.20
gen p2=0
replace p2=1 if t[5,2]< 0.25 & t[6,2]> 0.25
gen p3=0
replace p3=1 if t[5,3]< 0.50 & t[6,3]> 0.50
gen data=1
renvars coef1 coef2 coef3\bmiz seifa cons 
keep cons bmiz seifa secons sebmiz seseifa varIntercept varSlope varResid p1 p2 p3 data
save "../../Data analysis/Available data/Results/Available data",replace


forvalues i=2/1000{
use data`i'.dta,clear
di `i'
 mixed QoL bmiz seifa ||ID:bmiz,stddev iterate(20) 
matrix define t=r(table)
matrix list t
matrix define coef=e(b)
matrix define var=vecdiag(e(V))
matrix list coef
matrix list var
drop ID seifa sep bmiz QoL
svmat coef 
svmat var
gen varIntercept=exp(coef5) 
gen varSlope=exp(coef4)  
gen varResid= exp(coef6)
gen secons=sqrt(var3) 
gen sebmiz=sqrt(var1) 
gen seseifa=sqrt(var2) 
gen p1=0
replace p1=1 if t[5,1]< -0.20 & t[6,1]> -0.20
gen p2=0
replace p2=1 if t[5,2]< 0.25 & t[6,2]> 0.25
gen p3=0
replace p3=1 if t[5,3]< 0.50 & t[6,3]> 0.50
gen data=1
renvars coef1 coef2 coef3\bmiz seifa cons 
keep cons bmiz seifa secons sebmiz seseifa varIntercept varSlope varResid p1 p2 p3 data
append using "../../Data analysis/Available data/Results/Available data.dta"
save "../../Data analysis/Available data/Results/Available data.dta",replace
}


 

 
