### 2. FUNCTION LIBRARY ###

### GENERAL FUNCTIONS ### 

### In this script - all general function can be found that apply to both the unpenalized logistic regression and the penalized logistic regression scenario ### 

#logit function to be able to transform AUC and pool AUC 
logit <- function(x){log(x / (1-x))}

#estimate log c-statistic
#we estimate the log c-statistic in every imp set, to subsequently pool this log c-statistic and then transform back. 
compute_c_stat_log_imp <- function(obs_outcome, pred_outcome, method="delong", direction="<"){
  
  #estimate c_stat
  c_stat <- pROC::auc( response = obs_outcome,
                       predictor = pred_outcome,
                       direction = direction,
                       method = method)
  
  #estimate confidence intervals appropriately
  variance <- pROC::var(c_stat)
  logit_se <- sqrt(variance)/(c_stat*(1-c_stat))
  return(c (logitauc = logit(as.numeric(c_stat)),
            logit_se=logit_se))
}


# Rubin's rules for SE c-statistic---
rubins_rules_var_auc <- function( estimates,
                                  ses,
                                  n_imputed_sets){
  # within study variance
  within_var <- mean( ses^2)
  
  # between study variance
  between_var <- (1 + ( 1/ n_imputed_sets)) * var( estimates)
  
  # total variance
  total_var <- within_var + between_var
  
  return( total_var)
}

# Rubin's Rules to pool the standard errors of coefficients 

rubins_rules_var <- function(estimates, ses, n_imputed_sets){
  
  #within study variance
  within_var <- mean(ses^2)
  
  #between study variance
  theta_bar <- mean(estimates)
  
  between_var <- sum( (estimates - theta_bar)^2)/(n_imputed_sets-1)
  
  #total variance
  total_var <- within_var + between_var + between_var/n_imputed_sets 
  
  return(total_var)
}


#function for a majority vote 
majorityvote <- function(bindedset, impN, freq){
  nimp <- impN
  prop <- freq 
  calculatefreq <- function(nimp, prop) { V <- nimp * prop 
  return(V) } 
  vote <- calculatefreq(nimp, prop)
  
  
  finalmodel <- data.frame()
  
  count <- bindedset %>% group_by(variable) %>% mutate(freq=n()) %>% ungroup() %>% filter(freq > vote) %>% select(-freq)
  finalmodel <- count %>% group_by(variable) %>% summarise(Meancoef = mean(modeli.coefficients))
  return(finalmodel)
  
}

# Calculate model performance for the nullmodels
#calculate performance measures 
predcalcnull <- function(predimpsets, nullmodeldf){ 
  
  #input for the function: the imputationsets including the predicted values of the apparent model and the final nullmodel (pooled). 
  
  #store output: 
  totalpreds <- list()
  
  #specificy coefficients
  model_coefs <- with(nullmodeldf, setNames(Meancoef, variable))
  
  #loop over imputed datasets
  for(i in 1:length(predimpsets)){
    #specifiy each dataset
    predset <- predimpsets[[i]]
    
    #calculate null model predicted values 
    impsetprep <- predset %>% mutate_if(is.factor, as.numeric)
    impsetprep <- impsetprep %>% mutate(finalmodellp = as.numeric(model_coefs["(Intercept)"]))
    
    impsetpred <- predset
    impsetpred$lpnull <- impsetprep$finalmodellp
    impsetpred$prednull <- 1/(1+exp(-impsetpred$lpnull))
    totalpreds[[i]] <- impsetpred
  }
  return(totalpreds)
}


# Calibration plot: 
plot_smooth_calibration <- function( x_coordinates,
                                     y_coordinates){
  plot(x = 1,                 
       xlab = "Estimated probability", 
       ylab = "Observed risk",
       xlim = c(0, 1), 
       ylim = c(0, 1),
       main = "",
       type = "n")
  
  # add reference line
  abline( a = 0, b = 1, lty = 2)
  
  # add smooth calibration curve
  lines( x = x_coordinates,
         y = y_coordinates)
  
  # add histogram at bottom
  par(new = T)
  # from rms val.prob function
  lim <- c(0,1)
  bins <- seq(lim[1], lim[2], length=101)
  x <- x_coordinates[x_coordinates >= lim[1] & x_coordinates <= lim[2]]
  f <- table(cut(x, bins))
  j <- f > 0
  bins <- (bins[-101])[j]
  f <- f[j]
  f <- lim[1] + .15 * diff(lim) * f / max(f)
  segments(bins, 0, bins, f)
}

# Bootstrap function
#function for bootstrap: 
bootstrap <- function(dataimp, bootstrapn){
  
  #store output
  bootstraplist <- list()
  
  #set n for number of bootstrap samples 
  n <- bootstrapn
  
  #create loop
  for(i in 1:n){
    
    #make sure you sample with replacement 
    bootstrap <- sample_n(dataimp, nrow(dataimp), replace=TRUE)
    
    #store results in a list
    bootstraplist[[i]] <- bootstrap
  }
  return(bootstraplist)
}


# get bootstrap SE's for specific performance measures from dataframes and lists earlier constructed: 
getSEs <- function(performancelistbs, measure){
  
  #store output
  selist <- list()
  
  #loop over the imputed datasets with bootstraps
  for(i in 1:length(performancelistbs)){
    
    #take the list per imputation set and bind rows to create a stacked dataframe 
    performancelist <- performancelistbs[[i]]
    performancebinded <- dplyr::bind_rows(performancelist) 
    
    #chose the right performance measure 
    performancemeasure <- performancebinded %>% dplyr::filter(.data$performanceM == measure) 
    
    #filter on bootstrap values only (we use bootstrap model on bootstrap sample performance)
    performancemeasure2 <- performancemeasure %>% dplyr::select(performanceM, bootstrap)
    
    #calculate percentiles and subsequent the standard errors 
    upper <- quantile(performancemeasure2$bootstrap, probs=0.975)
    lower <- quantile(performancemeasure2$bootstrap, probs=0.025)
    se <- (upper-lower)/3.92
    selist[[i]] <- se
  }
  return(selist)
}





