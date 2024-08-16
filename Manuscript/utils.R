# some functions used in the manuscript

#' Visualize missing data (syst. vs. spor.)
#'
#' @return An object of class `ggplot`
plot_na <- function() {
  dat <- expand.grid(rows = 1:7, cols = 1:6) |>
    cbind(
      text = c("1", "1", "2", "2", "3", "", "N", rep("", 35)),
      miss = c(rep("", 16), "NA", "NA", "", "", "", "NA", "", "", "NA", rep("", 17))
    )
  ggplot(dat, aes(x = cols, y = rows)) +
    geom_tile(fill = "white",
              color = "black",
              linewidth = 0.5) +
    geom_text(aes(label = text), color = "black", size = 3) +
    geom_text(
      aes(label = miss),
      color = mice:::mdc(2),
      family = "mono",
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = 1:6,
      labels = c(
        "cluster",
        expression(X[1]),
        expression(X[2]),
        expression(X[3]),
        "...",
        expression(X[p])
      ),
      name = NULL,
      position = "top"
    ) +
    scale_y_continuous(
      breaks = 1:7,
      labels = c(1:5, "...", "n"),
      name = NULL,
      trans = "reverse"
    ) +
    # coord_cartesian(expand = c(0,0)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt"))
}

#' Amputation of the popmis datset
#' 
#' @return The popmis dataset from mice, with missing values.
ampute.popmis <- function(dat){
  set.seed(1234)
  dat$prob <- 1 - plogis(-5 + 0.3 * (dat$sex==1) + 0.7 * dat$texp)
  dat$teachpop <- ifelse(rbinom(nrow(dat), 1, dat$prob),
                         NA, dat$teachpop) # MAR teach popularity
  dat$texp[sample(1:nrow(dat), 0.05*nrow(dat))] <- NA #MCAR texp
  dat$prob <- NULL
  dat$const <- NULL
  return(dat)
}

plot_models <- function(list_models, mod_name){
  library(insight);library(tidyverse);library(broom.mixed);library(mitml)
  mod_table<-var_table<-NULL
  for(i in 1:length(list_models)){
    mod <- list_models[[i]]  
    if(class(mod)[[1]]=="mira"){
      # For model estimated on imputed dataset ----
      ## Get coefficient estimates -----
      mod_tab <-summary(pool(mod), conf.int = TRUE)[,c("term","estimate","2.5 %","97.5 %","p.value")] 
      ## Get variance estimates -----
      mitval <- data.frame(testEstimates(as.mitml.result(mod), 
                                         extra.pars = TRUE)[["extra.pars"]])
      names_var <- rownames(mitval)
      
      vartable <-data.frame(Term=character(0),Value=numeric(0))
      if (length(names_var[grepl("^(Residual~~Residual)", names_var)])!=0){
        vartable <- vartable|>
          add_row(Term="Residual variance",
                  Value=estim[grepl("^(Residual~~Residual)", names_var)])
      }
      if (length(names_var[grepl("^(Intercept~~Intercept)",
                                 names_var)])!=0){
        vartable <- vartable|>
          add_row(Term="Random intercept\nvariance",
                  Value=estim[grepl("^(Intercept~~Intercept)",
                                    names_var)]) }
      if (length(names_var[grepl("^(Intercept~~Slope)", 
                                 names_var)])!=0){
        vartable <- vartable|>
          add_row(Term=paste0("Random slope\nvariance ",
                              sub(".*\\|", "",
                                  names_var[grepl("^(Intercept~~Intercept)", names_var)])),
                  Value=estim[grepl("^(Intercept~~Slope)", names_var)])
      } 
      
      if (length(names_var[grepl("^(Slope~~Slope)", names_var)])!=0)
      {
        vartable <- vartable|>
          add_row(Term=paste0("Random slope-intercept correlation ",
                              
                              sub(".*\\|", "",
                                  names_var[grepl("^(Slope~~Slope)", names_var)])),
                  Value=estim[grepl("^(Slope~~Slope)",
                                    names_var)])
      }
      
      if ((length(names_var[grepl("^(ICC)", names_var)])!=0)){
        vartable <- vartable|>
          add_row(Term="ICC",
                  Value=estim[grepl("^(ICC)", names_var)])
      }
      
    }else{
      
      # For model on CC -----
      ## Get coefficient estimates -----
      mod_tab<-tidy(mod,effects="fixed",conf.int=TRUE)
      mod_tab <-mod_tab[,c("term","estimate","conf.low",
                           "conf.high","p.value")]
      rownames(mod_tab)<-NULL
      ## Get variance estimates -----
      varcomp <- insight::get_variance(mod)
      vartable <- data.frame(Term="Residual variance",
                             Value=as.numeric(varcomp$var.residual))
      if (!is.null(varcomp$var.intercept)){
        vartable <- vartable|>
          add_row(Term="Random intercept\nvariance",
                  Value=as.numeric(varcomp$var.intercept))
      }
      if (!is.null(varcomp$var.var.slope)){
        vartable <- vartable|>
          add_row(Term=paste0("Random slope var",
                              attr(varcomp$var.slope, "names")),
                  Value=as.numeric(varcomp$var.intercept))
      }
      
      
      if (!is.null(varcomp$cor.slope_intercept)){
        vartable <- vartable|>
          add_row(Term=paste0("Random slope-intercept corr",
                              attr(varcomp$cor.slope_intercept, "names")),
                  Value=as.numeric(varcomp$cor.slope_intercept))
      }
      if (!is.null(varcomp$var.intercept)){
        vartable <- vartable|>
          add_row(Term="ICC",
                  Value=as.numeric(varcomp$var.intercept)/
                    (as.numeric(varcomp$var.intercept)+
                       as.numeric(varcomp$var.residual)))
      }
    }
    colnames(mod_tab)<-c("Coefficient","Estimate","lower","upper",
                         "p.value")
    mod_tab$model<-mod_name[i]
    mod_table<-rbind(mod_table,mod_tab)  
    vartable$model <- mod_name[i]
    var_table <- rbind(var_table,vartable)
  }
  
  plot_coef <- ggplot(mod_table, aes(x = Coefficient, y = Estimate,
                                     ymin = lower, ymax = upper, group=model, color=model)) + 
    geom_hline(yintercept=0,linetype="dashed",color="black",size=0.1)+
    geom_linerange( position = position_dodge(width = 0.8)) + 
    geom_pointrange(size=0.1, position = position_dodge(width = 0.8))+
    coord_flip()+
    theme_classic()+theme(legend.position = "bottom")
  +scale_color_viridis_d()
  if(randisplay==T){
    var_table<-var_table|>filter(str_detect(Term, "^Random")) }
  plot_var <- ggplot(var_table, aes(x=Term, y=Value, fill=model))+
    geom_bar(stat="identity", position="dodge")+
    coord_flip()+scale_fill_viridis_d()+
    xlab("Term")+ylab("Estimate")+ 
    theme_light()+theme(legend.position = "bottom")
  
  
  return(list(plot_coef,plot_var))
}