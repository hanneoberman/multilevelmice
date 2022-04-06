# env
# library(dplyr)
# logistic <- function(x){exp(x)/(1+exp(x))} #equal to plogis()

fit_missingness_models <- function(fname) {
  data <- read.spss(fname, to.data.frame = TRUE)
  
  # Names of incomplete variables
  vnames_incomplete <- c("pupil", "ct", "hypoxia", "hypotens", "tsah", "edh")  
  vnames_complete <- c("age", "motor_score", "d_mort")
  
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
  
  # Change coding of variables to match the coding in metamisc
  data$motor_score <- factor(data$d_motor, levels = c(1, 2, 3, 4, 5, 6),
                             labels = c("1/2", "1/2", "3", "4", "5/6", "5/6"))
  data$pupil <- factor(data$d_pupil, levels = c(1, 2, 3),
                       labels = c("Both", "One", "None"))
  data$ct <- factor(data$ct_class, levels = c(1, 2, 3, 4, 5),
                    labels = c("I/II", "I/II", "III", "IV/V", "IV/V"))
  
  
  results <- data.frame(trial = character(),
                        variable_missing = character(),
                        percent_missing = numeric(),
                        "beta_Intercept" = numeric(),
                        beta_pupilOne =  numeric(),
                        beta_pupilNone =  numeric(),
                        beta_ctIII = numeric(),
                        beta_ctIV.V = numeric(),
                        beta_hypoxia = numeric(),
                        beta_hypotens = numeric(),
                        beta_tsah = numeric(),
                        beta_edh = numeric(),
                        beta_age = numeric(),
                        beta_motor_score3 = numeric(),
                        beta_motor_score4 = numeric(),
                        "beta_motor_score5.6" = numeric(),
                        beta_d_mort = numeric())
  



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
      cnames[which(cnames == "motor_score5/6")] <-  "motor_score5.6"
      cnames[which(cnames == "ctIV/V")] <-  "ctIV.V"
      cnames <- paste("beta_", cnames, sep = "")
      
      results <- results %>% add_row(data.frame(trial = unique(dat$name),
                                                variable_missing = vname,
                                                percent_missing = pct_missing[vname]))
      
      if (pct_missing[vname] > 0 &  pct_missing[vname] < 100) {
        results[nrow(results),cnames] <- coef(fit)
      }
    }
  }
  
  return(results)

}

ampute_impact <- function(impact, p_obs){
  # just for function development 
  # data("impact", package = "metamisc")
  # load("Data/p_observed_impact.rda")
  # p_obs <- results
  
  # pre-process coefficients
  p_obs[is.na(p_obs)] <- 0
  betas <- split(p_obs, p_obs$trial)
  betas <- purrr::map(betas, ~{.[, 4:17] %>% as.data.frame(row.names = unique(p_obs$variable_missing)) %>% relocate(beta_pupilOne, .after = beta_pupilNone)})

  # create model matrices
  mm <- purrr::map(split(impact, impact$name), ~{model.matrix(name ~ pupil + ct + hypox + hypots + tsah + edh + age + motor_score + mort, .x) %>% as.data.frame()})
  
  # calculate linear predictor for missingness probability and convert to prob_obs
  pr <- purrr::map2_dfr(betas, mm, ~{as.matrix(.y) %*% t(.x[]) %>% 
      plogis(.) %>% 
      as.data.frame() %>% 
      mutate(across(everything(), ~rbinom(n = nrow(.y), size = 1, prob = .))) %>% 
      setNames(paste0("R_", names(.))) %>% 
      cbind(id = rownames(.), .)})
  
  # select lowest values to make incomplete
  impact_NA <- impact %>% 
    cbind(id = rownames(.), .) %>% 
    left_join(pr, by = "id") %>% 
    mutate(pupil = ifelse(R_pupil == 1, pupil, NA),
           pupil = factor(pupil, levels = 1:3, labels = c("Both", "None", "One")),
           ct = ifelse(R_ct == 1, ct, NA),
           ct = factor(ct, levels = 1:3, labels = c("I/II", "III", "IV/V")),
           hypox = ifelse(R_hypoxia == 1, hypox, NA),
           hypots = ifelse(R_hypotens == 1, hypots, NA),
           tsah = ifelse(R_tsah == 1, tsah, NA),
           edh = ifelse(R_edh == 1, edh, NA),
           .keep = "unused")
  return(impact_NA[, -1])
}

# saveRDS(impact_NA, "Data/impact_NA.RDS")
