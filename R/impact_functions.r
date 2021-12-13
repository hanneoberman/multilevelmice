library(foreign)
library(dplyr)

fit_missingness_models <- function(fname) {
  data <- read.spss(fname, to.data.frame = TRUE)
  
  # Names of incomplete variables
  vnames_incomplete <- c("d_pupil", "ct_class", "hypoxia", "hypotens", "tsah", "edh") #i_pupil 
  vnames_complete <- c("age", "d_motor")
  
  data$name <- ""
  data$name[data$trial == "74"] <- "TINT" #1118
  data$name[data$trial == "75"] <- "TIUS" #1041
  data$name[data$trial == "76"] <- "SLIN" #409
  data$name[data$trial == "77"] <- "SAP" #919
  data$name[data$trial == "79"] <- "PEG" #1510
  data$name[data$trial == "80"] <- "HIT I" # 350
  data$name[data$trial == "81"] <- "UK4" #791
  data$name[data$trial == "83"] <- "TCDB" #603
  data$name[data$trial == "84"] <- "SKB" #126
  data$name[data$trial == "85"] <- "EBIC" #822
  data$name[data$trial == "86"] <- "HIT II" #819
  data$name[data$trial == "88"] <- "NABIS" #385
  data$name[data$trial == "89"] <- "CSTAT" # 517
  data$name[data$trial == "90"] <- "PHARMOS" #856
  data$name[data$trial == "91"] <- "APOE" # 756
  
  
  results <- data.frame(trial = character(),
                        variable_missing = character(),
                        percent_missing = numeric(),
                        "beta_Intercept" = numeric(),
                        beta_d_pupil =  numeric(),
                        beta_ct_class = numeric(),
                        beta_hypoxia = numeric(),
                        beta_hypotens = numeric(),
                        beta_tsah = numeric(),
                        beta_edh = numeric(),
                        beta_age = numeric(),
                        beta_d_motor = numeric())
  



  # Iterate over the studies
  studyids <- unique(data$trial)
  
  # 
  for (sid in studyids) {
    dat <- subset(data, trial == sid)[,c("name", vnames_complete, vnames_incomplete)]
    
    nmissing <- colSums(apply(dat, 2, is.na))
    pct_missing <- 100*nmissing/nrow(dat)
    
    for (vname in vnames_incomplete) {
      dat$observed <- ifelse(is.na(dat[,vname]),0, 1)
      vnames_predictors <- setdiff(c(vnames_complete, vnames_incomplete), vname)
      
      # Remove variables that are systematically missing from the predictor space
      vnames_predictors <- setdiff(vnames_predictors, names(nmissing)[pct_missing == 100])
      
      fmla <- paste("observed ~", paste(vnames_predictors, collapse = "+"))
      
      fit <- glm(fmla, data = dat, family = binomial())
      
      cnames <- names(coef(fit))
      cnames[1] <- "Intercept"
      cnames <- paste("beta_", cnames, sep = "")
      
      results <- results %>% add_row(data.frame(trial = unique(dat$name),
                                                variable_missing = vname,
                                                percent_missing = pct_missing[vname]))
      
      results[nrow(results),cnames] <- coef(fit)
      
    }
  }
  
  return(results)
  
}