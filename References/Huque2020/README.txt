Title of the manuscript: Multiple imputation methods for handling incomplete longitudinal and clustered data where the target analysis is a linear mixed effects model

Authors: Md Hamidul Huque, Margarita Moreno-Betancur, Matteo Quartagno, Julie A. Simpson, John B. Carlin, and Katherine J. Lee

Correponsing author for codes: Md Hamidul Huque, hamidul_b7@yahoo.com
Configuration: 

library(mitml)
library(reshape2)
library(mgcv)
library(lme4)
library(plyr)
#install.packages("mice")
library(mice)
library(foreign)
library(jomo)

> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17763)

Matrix products: default

locale:
[1] LC_COLLATE=English_Australia.1252  LC_CTYPE=English_Australia.1252    LC_MONETARY=English_Australia.1252
[4] LC_NUMERIC=C                       LC_TIME=English_Australia.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] jomo_2.6-7      foreign_0.8-71  mice_3.5.0      lattice_0.20-38 plyr_1.8.4      lme4_1.1-21     Matrix_1.2-17  
 [8] mgcv_1.8-28     nlme_3.1-139    reshape2_1.4.3  mitml_0.3-7    

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1        compiler_3.6.0    pillar_1.4.0      nloptr_1.2.1      tools_3.6.0       rpart_4.1-15     
 [7] boot_1.3-22       tibble_2.1.1      pkgconfig_2.0.2   rlang_0.3.4       rstudioapi_0.10   parallel_3.6.0   
[13] dplyr_0.8.1       stringr_1.4.0     generics_0.0.2    nnet_7.3-12       grid_3.6.0        tidyselect_0.2.5 
[19] glue_1.3.1        R6_2.4.0          survival_2.44-1.1 minqa_1.2.4       purrr_0.3.2       tidyr_0.8.3      
[25] magrittr_1.5      backports_1.1.4   MASS_7.3-51.4     splines_3.6.0     assertthat_0.2.1  stringi_1.4.3    
[31] broom_0.5.2       crayon_1.3.4      pan_1.6

Stata:

Stata/SE 16.0 for windows (64-bit x86-64)
Physical memory 16.00 GB

Data directory:
The main Simulation folder contains two subfolder: Clustered data and Longitudinal data

Cluster data folder again contains two subfolder based on sample sizes 1) 100 clusters: this folders provides analysis and results for Appendix table B3 and B4, corresponding to scenario (iii) and (iv), respectively (2)300 clusters: this folders provides analysis and results for Table 3 and 4 in the main paper. 

Longitudinal data folder contains two subfolder based on sample sizes 1) 1000 samples: this folders provides analysis and results for Appendix table B1 and B2, corresponding to scenario (i) and (ii), respectively (2)5000 samples: this folder provides analysis and results for Table 1 and 2 in the main paper.

JM-MVNI and FCS-Standard methods were evaluated using Stata (This may take while to complete 1000 replications).
All others methods were evaluated using high performance computing cluster.

For all the methods except JM-MVNI, FCS-standard and Available data analysis, generated data in Stata format (in "\Data simulation\Data" folder) were analyzed using R software. These datasets are located in the `\Data simulation\Data` folder under each of the simulation scenarios i to iv. 

As indicated in the main paper that the scenario (i) and (ii) corresponds to longitudinal data analysis for Table 1 and 2, respectively in the main paper for 5000 samples and appendix Table B1 and B2, respectively for 1000 samples. Similarly, scenarios (iii) and (iv) corresponds to the cluster data analysis results showed in Table 3 and 4 with 300 clusters and resutls for smaller sample sizes (100 clusters) as shown in the appedix table B3 and B4. Due to computation time, we recommend user should reproduce all the results with sample sample sizes (i.e., results presented in the appendix first) before moving to the large sample.

We have provided all the codes to reproduce the results in the `Data analysis` which contains subfolder as the name of the methods (heading of each Table) evaluated.In order to reproduce the results presented in Tables. The user needs to analysis datasets ("\Data simulation\Data") using all the codes in the `Data analysis` folder under corresponding scenarios, then combine those results using the combinedDataHPC.R followed by 'comparison of methods.do file` available in the `Combined results from all methods













All the analysis results were provided in the "Data analysis folder using subfolder in accordance with the name of the method evaluated.