program define dataGen, rclass
version 13.0
syntax
drop _all
set obs 100
gen ID = _n
gen seifa=rnormal(0,1)
gen clusterSize=2+int((10-2+1)*runiform())
*** Generate cluster specific random intercept and slope
gen b_0=rnormal(0,0.4)
gen b_1=rnormal(0,0.2)
gen bm_0=rnormal(0,0.15)
expand clusterSize
drop clusterSize
gen lang=uniform() <= 0.90
qui gen moedu = uniform()<=0.60
*** socio-economic position
gen sep=-4.7+0.8*moedu+0.01*seifa+0.01*lang+rnormal(0,0.9)

*** Generate age
generate age=(168 +rnormal(11,1.5))/12
gen sex= uniform() <= 0.5
*** generate bmiz scores
gen bmiz=(-1.00+bm_0)+0.10*age+0.05*sex+rnormal(0,1)

* Generation of the QoL
generate QoL=(0.50+b_0)+(-0.2+b_1)*bmiz+(0.25)*seifa+rnormal(0,0.9)
quietly drop  age sex b_0 b_1 bm_0
end

program define dataGenM, rclass
version 13.0
syntax
quietly dataGen

/***********************************************
****	Missing data using MAR assumption	****
************************************************/
replace bmiz=. if runiform()<invlogit(-2.2+1.0*QoL+0.2*seifa-0.2*sep) 

end

cd "~\Simulation\Clustered data\100 clusters\Simulation scenario (iii)\sample 100\Data simulation\Data"
set seed 12345
forvalues i = 1(1)1000{
quietly dataGenM
drop moedu lang
saveold data`i'.dta,replace v(11)
}

 