### 2.3. FUNCTION LIBRARY FIRTH SPECIFIC FUNCTIONS ### 

### FIRTH SPECIFIC FUNCTIONS ###

backwardimpF <- function(impsets, formula){
  
  #store output
  models <- list()
  
  #loop over imputed datasets
  for (i in 1:length(impsets)){
    impset <<- impsets[[i]] #<<- is used because it seems that logistf cannot look into function environment so we had to save it in the general enviroment. NB: if you reuse this script, make sure you don't use impset in another function, because it is invisibly stored in general environment and R will look into that first.
    
    #create models using logistf
    full <- logistf(formula, data=impset)
    bwmodel <- backward(full, slstay=0.157, trace=F, printwork=F, data=impset) #specify AIC criterion, trace T will print everything, which can be usefull but it is faster to set to F. 
    bwmodelF <- flic(bwmodel)
    models[[i]] <- bwmodelF
  }
  return(models)
}

#extract models
extractmodels <- function(modelimplist){
  
  #store output 
  listmodeldfs <- list()
  
  #create loop over models of imputed datasets 
  for (i in 1:length(modelimplist)){
    modeli <- modelimplist[[i]]
    modeldf <- data.frame(modeli$coefficients)
    
    #get standard errors from logistf output
    modeldf2 <- as.data.frame(sqrt(diag(modeli$var)))
    
    #transformations to get the correct column names 
    modeldf3 <- modeldf %>% rownames_to_column("variable")
    modeldf4 <- modeldf2 %>% rownames_to_column("variable")
    modeldf4 <- modeldf4 %>% rename(se = "sqrt(diag(modeli$var))")
    
    #merge everything together and create output 
    modeldf5 <- merge(modeldf3, modeldf4, by="variable", all.x=TRUE)
    listmodeldfs[[i]] <- modeldf5
  }
  return(listmodeldfs)
}

#calculate predicted values 
predcalcF <- function(impsets, finalmodeldf){
  
  #store output
  predimps <- list()
  model_coefs <- with(finalmodeldf, setNames(Meancoef, variable))
  
  #create loop over imputed datasets, using the final model coefficients, multiplying them with the variable values 
  for(i in 1:length(impsets)){
    datasetimp <- impsets[[i]]
    
    #we need to set everything to numeric 
    impsetprep <- datasetimp %>% mutate_if(is.factor, as.numeric)
    
    #multiply all model coefficients with their respective variable value, not using the intercept variable (-1)
    impsetprep <- impsetprep %>% mutate(finalmodellp = model_coefs["(Intercept)"] + rowSums(across(all_of(names(model_coefs)[-1]),
                                                                                                   ~ get(cur_column()) * model_coefs[cur_column()]))) 
    
    #paste this value to the imputed dataset as predicted value 
    impsetpred <- datasetimp
    impsetpred$lp <- impsetprep$finalmodellp
    impsetpred$pred <- 1/(1+exp(-impsetpred$lp))
    predimps[[i]] <- impsetpred
  }
  return(predimps)
}

#nullmodels
nullmodels <- function(impsets){
  
  #store output
  nullmodels <- list()
  
  #loop over imputed datasets
  for (i in 1:length(impsets)){
    set <- impsets[[i]]
    #create models
    null <- logistf(outcome~1, data=set)
    nullmodels[[i]] <- null
  }
  return(nullmodels)
}

#calculate performance measures
impPerfF <- function(impsetspred){
  
  #store output
  performancetablesimp <- list()
  
  #create a loop over the imputed datasets, that now contains predicted values (pred) of the final model and predicted values of the nullmodel (prednull). 
  for(i in 1:length(impsetspred)){
    
    #calculate auc en logit_auc --> save the logit. We make sure we store these as a variable calles performanceM for all calculations, so we can later merge these sets together. 
    totalpred <- impsetspred[[i]] 
    performanceimp <- compute_c_stat_log_imp(obs_outcome=totalpred$outcome, pred_outcome=totalpred$pred)
    c_statlogdf <- as.data.frame(performanceimp)
    c_statlogdf <- c_statlogdf %>% rownames_to_column("performanceM")
    
    #calculate other performance statistics using val.prob, storing them as performance M variable in a seperate dataframe 
    performanceimp <- val.prob(totalpred$pred, as.numeric(totalpred$outcome), pl=FALSE)
    performancedf <- data.frame(performanceimp)
    performancedf <- performancedf %>% rownames_to_column("performanceM")
    performancedf <- performancedf %>% filter(!is.na(performanceM))
    
    #calculate performance statistics for the null model to calculate IPA further on
    performancenull <- val.prob(totalpred$prednull, as.numeric(totalpred$outcome), pl=FALSE)
    performancenulldf <- as.data.frame(performancenull)
    performancenulldf <- performancenulldf %>% rownames_to_column("performanceM")
    performancenulldf <- performancenulldf %>% filter(!is.na(performanceM)) ####dit moest ik erbij doen omdat hij om een of andere reden NA's kreeg, omdat er 1 variabele in kolom performanceM NA werd, en die nam hij mee met de Brier. 
    
    #now calculate the performance measures that are not returned by val.prop by hand (IPA, OE). We also transformed IPA when we were not yet sure how to pool IPA. 
    add.data.IPA <- data.frame(performanceM = "IPA", performanceimp = 1-((performancedf$performanceimp[performancedf$performanceM=="Brier"])/(performancenulldf$performancenull[performancenulldf$performanceM=="Brier"])))
    add.data.OE <- data.frame(performanceM = "OE", performanceimp = (mean(totalpred$outcome))/(mean(totalpred$pred)))
    add.data.BrierR <- data.frame(performanceM="BrierR", performanceimp = 1 - (mean((totalpred$outcome - totalpred$pred)^2) / (mean(totalpred$outcome) * (1 - mean(totalpred$outcome)))))
    dfperformance2 <- bind_rows(performancedf, c_statlogdf, add.data.IPA, add.data.OE, add.data.BrierR) #wrap it all toghether                          
    dfperformance3 <- dfperformance2 %>% filter(performanceM=="logitauc"|performanceM=="logit_se"|performanceM=="IPA"| performanceM=="Slope" | performanceM=="Intercept" | performanceM=="Brier" | performanceM=="BrierR" | performanceM=="OE")
    performancetablesimp[[i]] <- dfperformance3
    
  }
  return(performancetablesimp)
}


#calibration script Firth specific
# provide coordinates of the smooth calibration curve ----
draw_smooth_calibrationF <- function( obs_outcome, pred_outcome){
  # dit moet wss aangepast naar eigen data, moet alle geobserveerde obs en pred
  # waarden zijn uit stacked imputed set.
  y_input <- t( rbind( unlist( obs_outcome)))
  x_input <- t( rbind( unlist( pred_outcome)))
  
  y_coordinates <- loess( y_input ~ x_input)$fitted #-1 -1 was necessary for GLM but not Firth (already numeric). 
  
  # store all coordinates of smooth curve (ordered based on x)
  smooth_cal_slope_info <- list( x_coordinates = x_input[order( x_input)],
                                 y_coordinates = y_coordinates[order( x_input)])
  
  
  # output overall info calibration slope
  return( smooth_cal_slope_info)
}


# Calculate performance measures over bootstrap models: 
impPerfBSF <- function(imputatiesets, bootstrapsets, bsmodels, bsnullmodels){
  
  #store output as two seperate lists, one for performance in the bootstrap samples and one for performance in the original imputation sets. 
  performancetablesbs <- list()
  performancetablesimp <- list()
  
  
  #first create a loop length of imputation sets, creating the imputation sets, the list with their respective bootstrap samples (individual), list of models for each set and list of null models for each set.  
  for(i in 1:length(imputatiesets)){
    
    imputatieset <- imputatiesets[[i]]
    bootstraplist <- bootstrapsets[[i]]
    modelbslist <- bsmodels[[i]]
    nullmodelbslist <- bsnullmodels[[i]]
    
    #then loop within the bootstrap samples lists, using the models and the nullmodels. 
    for(j in 1:length(bootstraplist)){
      bootstrapdf <- bootstraplist[[j]]
      modelbs <- modelbslist[[j]]
      nullmodelbs <- nullmodelbslist[[j]]
      bootstrapdf$pred <- predict(modelbs, type="response", newdata=bootstrapdf)
      
      ### Bootstrap on bootstrap performance
      #auc bootstrapmodel in bootstrapsample
      bootstrap <- compute_c_stat_log_imp(obs_outcome=bootstrapdf$outcome, pred_outcome=bootstrapdf$pred)
      c_statlogdf <- as.data.frame(bootstrap)
      c_statlogdf <- c_statlogdf %>% rownames_to_column("performanceM")
      
      #use null models for calculating IPA (bootstrapmodel in bootstrap sample)
      bootstrapdf$prednull <- predict(nullmodelbs, type="response", newdata=bootstrapdf)
      null <- val.prob(bootstrapdf$prednull, as.numeric(bootstrapdf$outcome), pl=FALSE)
      performancenulldf <- as.data.frame(null)
      performancenulldf <- performancenulldf %>% rownames_to_column("performanceM")
      performancenulldf <- performancenulldf %>% filter(!is.na(performanceM))
      
      #andere performance measures bootstrapmodels in bootstrapsample 
      # first using val.prob
      bootstrap <- val.prob(bootstrapdf$pred, as.numeric(bootstrapdf$outcome), pl=FALSE)
      bootstrapdf2 <- data.frame(bootstrap)
      bootstrapdf2 <- bootstrapdf2 %>% rownames_to_column("performanceM")
      
      #then calculating IPA and O:E ratio. Storing all these measures as performanceM as well in order to later merge all the datasets. 
      add.data.IPAbs <- data.frame(performanceM = "IPA" , bootstrap = (1-((bootstrapdf2$bootstrap[bootstrapdf2$performanceM=="Brier"])/(performancenulldf$null[performancenulldf$performanceM=="Brier"])))) 
      add.data.IPAbslogit <- data.frame(performanceM="logitIPA", bootstrap = logit(as.numeric(add.data.IPAbs$bootstrap)))
      add.dataOE <- data.frame(performanceM = "OE", bootstrap = as.numeric((mean(as.numeric(bootstrapdf$outcome))))/(mean(bootstrapdf$pred)))
      add.dataBrierR <- data.frame(performanceM="BrierR", bootstrap = 1 - (mean((bootstrapdf$outcome - bootstrapdf$pred)^2) / (mean(bootstrapdf$outcome) * (1 - mean(bootstrapdf$outcome)))))
      bootstrapdf3 <- bind_rows(bootstrapdf2, c_statlogdf, add.data.IPAbs, add.data.IPAbslogit, add.dataOE, add.dataBrierR) 
      
      #merge everything together 
      bootstrapdf4 <- bootstrapdf3 %>% filter(performanceM=="C (ROC)"|performanceM=="Brier" | performanceM=="BrierR" | performanceM=="logitauc"|performanceM=="logit_se"|performanceM=="IPA"|performanceM=="logitIPA"|performanceM=="Slope" | performanceM=="Intercept" | performanceM=="OE") 
      
      
      ### IMPUTATION set
      #repeat this process for testing performance on the original imputation set
      imputatieset$pred <- predict(modelbs, type="response", newdata=imputatieset)
      
      #auc
      imp <- compute_c_stat_log_imp(obs_outcome=imputatieset$outcome, pred_outcome=imputatieset$pred)
      c_statlogdfimp <- as.data.frame(imp)
      c_statlogdfimp <- c_statlogdfimp %>% rownames_to_column("performanceM")
      
      #null model predicted values to calculate IPA later
      imputatieset$prednull <- predict(nullmodelbs, type="response", newdata=imputatieset)
      null <- val.prob(imputatieset$prednull, as.numeric(imputatieset$outcome), pl=FALSE)
      impnulldf <- as.data.frame(null)
      impnulldf <- impnulldf %>% rownames_to_column("performanceM")
      impnulldf <- impnulldf %>% filter(!is.na(performanceM))
      
      #performance measures using val.prob
      imp <- val.prob(imputatieset$pred, as.numeric(imputatieset$outcome), pl=FALSE)
      impdf <- data.frame(imp)
      impdf <- impdf %>% rownames_to_column("performanceM")
      
      #performance measures calculated by hand (IPA and O:E)
      add.data.IPAimp <- data.frame(performanceM = "IPA" , imp=(1-((impdf$imp[impdf$performanceM=="Brier"])/(impnulldf$null[impnulldf$performanceM=="Brier"])))) 
      add.data.IPAimplogit <- data.frame(performanceM="logitIPA", imp = logit(as.numeric(add.data.IPAimp$imp)))
      add.dataOEimp <- data.frame(performanceM = "OE", imp = as.numeric(mean(as.numeric(imputatieset$outcome)))/(mean(imputatieset$pred)))
      add.dataBrierR <- data.frame(performanceM="BrierR", imp = 1 - (mean((imputatieset$outcome - imputatieset$pred)^2) / (mean(imputatieset$outcome) * (1 - mean(imputatieset$outcome)))))
      impdf2 <- bind_rows(impdf, c_statlogdfimp, add.data.IPAimp, add.data.IPAimplogit, add.dataOEimp, add.dataBrierR) 
      
      #merge everything together              
      impdf3 <- impdf2 %>% filter(performanceM=="C (ROC)"|performanceM=="logitauc"|performanceM=="logit_se"| performanceM=="Brier" | performanceM=="BrierR" | performanceM=="IPA"|performanceM=="logitIPA"|performanceM=="Slope" | performanceM=="Intercept" | performanceM== "OE")
      
      performance <- merge(bootstrapdf4, impdf3, by="performanceM")
      performance$optimism <- performance$bootstrap - performance$imp
      
      #store the performance measures for individual bootstrap models
      performancetablesbs[[j]] <- performance
      
    }
    
    #store the bootstrap lists in the imputation list
    performancetablesimp[[i]] <- performancetablesbs
    
  }
  return(performancetablesimp)
}


