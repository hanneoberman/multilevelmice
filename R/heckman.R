mice.impute.HeckmanIpdma<- function(y,ry,x, type, shrinkage_o=TRUE,shrinkage_s=TRUE,...){
  ##first 3 arguments are mandatory in mice
  ##type = how a var is used in the predictor matrix, new type is -3 for exclusion restriction var
  ##shrinkage in outcome/selection, may not be necessary 
  
  # 0. Define functions ----
  
  # 0.1. CopulaIPD: Apply Binomial or Gaussian model depending on y type
  ## this is fitting the copula equation in each cluster
  copulaIPD <- function(data=data,formula=formula,y=y){
    data$y<-as.numeric(as.character(data$y))
    
    if(class(y) == "numeric"){ #y continuos
      fit <- try(GJRM::copulaSampleSel(formula=formula, data=data), silent=TRUE)
      if(class(fit) == "try-error"){ # model is not estimable
        fit <- NA}
    }else if (class(y)=="factor" & nlevels(y)==2){ #y binomial
      fit <- try(GJRM::SemiParBIV(formula=formula, data=data, Model="BSS"), silent=TRUE)
      if(class(fit) == "try-error"){
        fit <- NA}
    }
    return(fit)
  }
  
  
  # 0.2. cov_mat_vector: Transform covariance matrix in a ordered vector
  ## additional function to transform the outcome for the copula, used later to calculate the marginals for the random effect param
  cov_mat_vector <- function(cov_mat,vnames){
    
    cov_mat[upper.tri(cov_mat)] <- "Up"
    cov_vec <- as.vector(cov_mat[vnames,vnames])
    cov_vec <- as.numeric(cov_vec[cov_vec!="Up"])
    return(cov_vec)
  }
  
  # 0.3. draw_theta_psi_mar: Estimate true effect size and draw a marginal theta and psi=var(theta) .
  ## psi is the standard error
  draw_theta_psi_mar<-function(coef_mat_s,Vb_list,vnames=NULL){
    
    if (is.null(vnames)){
      vnames <- colnames(coef_mat_s)} #use all set of parameters
    
    # Get covariance matrix
    coef_mat_s <- coef_mat_s[,vnames]
    cov_mat_s  <- do.call(rbind,lapply(Vb_list, cov_mat_vector, vnames=vnames))
    
    Marg_ind <- ifelse(length(Vb_list)!=1,0,1) #Indicator if Heckamn marginal model should be applied
    
    # Apply multivariate random-effects meta-analysis
    if(!Marg_ind){
      mvma <- try(mixmeta::mixmeta(coef_mat_s, cov_mat_s, method="reml", control=list(hessian=TRUE)),silent=TRUE)
      if(class(mvma)=="try_error"){
        mvma <- try(mixmeta::mixmeta(coef_mat_s, cov_mat_s, method="ml", control=list(hessian=TRUE)),silent=TRUE)
        if(class(mvma)=="try-error"){
          Marg_ind <- 1
        }}}
    
    if( Marg_ind ){ # marginal or only 1 study with copula estimation
      theta_star <- MASS::mvrnorm(n=1,mu=coef_mat_s, Sigma=Vb_list[[1]][vnames,vnames])
      psi_star<-Vb_list[[1]][vnames,vnames]
    }else{
      # Draw effects theta_star
      theta_star <- MASS::mvrnorm(n=1,mu=coef(mvma), Sigma=vcov(mvma))
      
      # Draw random effect, psi_star
      if(length(vnames)==1){
        qrsh <- 1/mvma$hessian
      }else{
        Hes <- as.matrix(Matrix::forceSymmetric(mvma$hessian))
        qrsh<-lqmm::make.positive.definite(solve(-Hes,tol = 1e-300))
      }
      
      rpar <- mvtnorm::rmvnorm(1, mean=mvma$par, sigma=qrsh, method="svd")
      
      if(length(vnames)==1){
        psi_star <- rpar^2
      }else{
        psi <- matrix(0, ncol( mvma$Psi), ncol(mvma$Psi))
        psi[lower.tri(psi,diag=TRUE)] <- rpar
        psi_star <- tcrossprod(psi)
      }
    }
    return(list(theta_star,psi_star))
  }
  
  
  #Conditional posterior distribution for shrinkage option
  ##estimate with or without shrinkage
  draw_cond_theta <- function(theta_mar,theta_k,var_theta_k,vnames,shrinkage){
    W_m <- ginv(theta_mar[[2]])
    W_k <- ginv(var_theta_k[vnames,vnames])
    S <- ginv(W_m+W_k)
    mu <- S%*%( W_k%*%as.vector(theta_k[vnames]) + W_m%*%as.vector(theta_mar[[1]]))
    if(shrinkage==TRUE){
      theta_star_i <- MASS::mvrnorm(n=1, mu=as.vector(mu), Sigma=S)
    }else{
      theta_star_i <- MASS::mvrnorm(n=1, mu=as.vector(theta_k[vnames]), Sigma=as.matrix(var_theta_k[vnames,vnames]))
    }
    return(theta_star_i)
  }
  
  ## there the actual function starts
  # 1. Define variables and dataset----
  #Rename covariates
  ## convert to design matrix (with dummies for factors)
  colnames(x) <- paste0("xv",1:length(colnames(x))) #change the covariates name for avoiding conflicts when y is covariate
  excl <- colnames(x)[type==-3] #names of exclusion variables
  nam<-colnames(x)
  
  # Define outcome and selection equation
  ## outcome equation = main model
  ## selection equation = exclusion restriction model
  out <- as.formula(paste0("y","~", paste(nam[type == 1],collapse="+")))
  sel <- as.formula(paste0("ry","~", paste(nam[type %in%c(1,-3)],collapse="+")))
  
  # Define data & prediction matrix
  data <- data.frame(ry, y, x)
  X <- data.frame(cbind(rep(1, nrow(x)),x))
  
  ## when the copula in one of the clusters fails, use this
  if (length(colnames(x)[type==-2])==0){
    warning("No clusterized incomplete variable, it will be imputed with no study specific model")
    names.clust<-"0"
    groupv<-rep(0,nrow(x))
    X$group_all<-groupv
    group_name<-"group_all"}
  else{
    group_name <- nam[type==-2]
    names.clust <- as.character(unique(x[,group_name]))
    groupv <- x[,group_name] # Define group vector
  }
  
  
  # Separate dataframe by group
  
  df_list <-  split.data.frame(data,groupv)
  
  #2. Step 1: Get the theta estimators for each study ----
  
  res_list <- suppressWarnings( lapply( df_list, copulaIPD, formula = list(sel,out), y = y))
  val<-sapply(res_list,function(x) typeof(x)) # identify systematic missing studies.
  res_list<-res_list[val=="list"] #keep only the clusters where copula function was calculated
  
  
  if(!any(val=="list")){
    warning("There is not enough information to impute at the study level, so the imputation is based on the total database without study specification.")
    res_list <- list("M"=copulaIPD(data=data,formula=list(sel,out),y=y))
    
    if(typeof(res_list[["M"]])!="list"){
      stop("There is not enough infomation to impute the Heckman model at marginal or study level")
    }
  }
  
  
  # Get effect vector for each study
  coef_list <- lapply(res_list, `[[`, c('coefficients'))
  
  
  # Get covariance matrix for each study
  Vb_list <- lapply(res_list, `[[`, c('Vb'))
  
  # Rename coeficient parameters in case that there are factor covariates
  nameslist <- names(coef_list) # get the names of a dply list.
  
  nameseq <- names(coef_list[[nameslist[1]]]) #equation names
  ostart <- which(nameseq =="(Intercept)")[2] #index where start outcome equation
  tstart <- which(nameseq =="theta.star") #index of theta parameter
  oend <- ifelse( class(y) == "numeric",tstart-2,tstart-1)
  
  selnaml <- paste(nameseq[1:(ostart-1)],"s",sep="_") #selection names with factor levels
  outnaml <- nameseq[ostart:oend] #outcome names with factor levels
  
  for(i in nameslist) { #rename columns and rows of selection variables
    names(coef_list[[i]])[1:(ostart-1)] <- selnaml
    colnames(Vb_list[[i]])[1:(ostart-1)] <- selnaml
    rownames(Vb_list[[i]])[1:(ostart-1)] <- selnaml
  }
  
  # Transform effects vector in a matrix with studies in rows
  coef_mat_s <- do.call(rbind, coef_list)
  
  
  #3. Step2: Get marginal theta and var(theta) ----
  
  #vnames=[selnam, outnam, sigma.star,theta.star]
  
  #3.1 Draw thetas from the universe of populations ----
  beta_s_mar <- draw_theta_psi_mar(coef_mat_s = coef_mat_s,Vb_list = Vb_list, vnames=selnaml)
  beta_o_mar <- draw_theta_psi_mar(coef_mat_s = coef_mat_s,Vb_list = Vb_list, vnames=outnaml)
  if (class(y) == "numeric"){
    sigma_t_mar <- draw_theta_psi_mar(coef_mat_s = coef_mat_s,Vb_list = Vb_list,vnames="sigma.star")} #copula package calls log(sigma) as sigma.star
  rho_t_mar <- draw_theta_psi_mar(coef_mat_s = coef_mat_s,Vb_list = Vb_list,vnames="theta.star")  #copula package calls atanh(rho) as theta.star
  
  #4. Get theta_k and var(theta_k) from full conditional distribution ----
  
  for (i in names.clust){ #Loop across studies
    
    if( !i%in%names(res_list)| i=="0"){ #systematically missing
      beta_s_star <- MASS::mvrnorm(n=1, mu=beta_s_mar[[1]], Sigma=beta_s_mar[[2]])
      beta_o_star <- MASS::mvrnorm(n=1, mu=beta_o_mar[[1]], Sigma=beta_o_mar[[2]])
      if (class(y) == "numeric"){
        sigma_star <- exp(MASS::mvrnorm(n=1, mu=sigma_t_mar[[1]], Sigma=sigma_t_mar[[2]]))}
      rho_star <- tanh(MASS::mvrnorm(n=1, mu=rho_t_mar[[1]], Sigma=rho_t_mar[[2]]))
      
    }else{ # sporadically missing
      beta_s_star <- draw_cond_theta(beta_s_mar,coef_list[[i]],Vb_list[[i]],selnaml,shrinkage_s) # no shrinkage.
      beta_o_star <- draw_cond_theta(beta_o_mar,coef_list[[i]],Vb_list[[i]],outnaml,shrinkage_o) #all the time shrinkage better.
      if (class(y) == "numeric"){
        sigma_star <- exp(draw_cond_theta(sigma_t_mar,coef_list[[i]],Vb_list[[i]],"sigma.star",shrinkage_o))}#copula package calls log(sigma) as sigma.star
      rho_star <- tanh(draw_cond_theta(rho_t_mar,coef_list[[i]],Vb_list[[i]],"theta.star",shrinkage_s&shrinkage_o)) #copula package calls atanh(rho) as theta.star
    }
    
    XSBS <- data.matrix(X[!ry & groupv==as.numeric(i), !colnames(X) %in%c(group_name)])%*%as.vector(beta_s_star)
    XOBO <- data.matrix(X[!ry & groupv==as.numeric(i), !colnames(X)%in%c(group_name,excl)])%*%as.vector(beta_o_star)
    
    
    if (class(y) == "numeric"){
      Ratio<-(-dnorm(XSBS)/(pnorm(-XSBS)))
      Ratio[is.na(Ratio)|is.infinite(Ratio)]<-0.0
      
      y.star <- XOBO + as.numeric(sigma_star)*as.numeric(rho_star)*Ratio + rnorm( nrow(XSBS), 0, sd=sigma_star)
      
    }else if (class(y) == "factor" & nlevels(y) == 2){
      
      p.star <- pbivnorm::pbivnorm(as.vector(XOBO), -as.vector(XSBS), -as.numeric(rho_star))/pnorm(-XSBS)
      #in case p.star return exactly 0 or 1 values
      p.star[p.star=="1"] <- 1.0
      p.star[p.star=="0"] <- 0.0
      
      ## the value that we'll impute
      y.star <- rbinom(nrow(XSBS), 1, as.numeric(p.star))
    }
    
    y[ !ry & groupv==as.numeric(i)] <- y.star
  }
  
  return(y[!ry])
}