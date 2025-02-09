### 2.2. Function library GLM specific functions ### 

### GLM-SPECIFIC FUNCTIONS ###

#Function for GLM model development in every individual imputation set i using backward selection
backwardimpglm <- function(impsets, formula){
  
  #store output
  models <- list()
  
  #loop over imputed datasets
  for (i in 1:length(impsets)){
    impset <- impsets[[i]] 
    full <- glm(formula, family=binomial, data=impset)
    modeli <- step(full, direction="backward", trace=0)
    models[[i]] <- modeli
  }
  return(models)
}

#function to extract the models from the above generated output (i models for i imputation sets)
extractmodels <- function(modelimplist){
  
  #store output 
  listmodeldfs <- list()
  
  #create loop over models of imputed datasets 
  for (i in 1:length(modelimplist)){
    modeli <- modelimplist[[i]]
    modeldf <- data.frame(modeli$coefficients)
    
    #get standard errors from glm output
    modeldf2 <- as.data.frame(sqrt(diag(vcov(modeli))))
    
    #transformations to get the correct column names 
    modeldf3 <- modeldf %>% rownames_to_column("variable")
    modeldf4 <- modeldf2 %>% rownames_to_column("variable")
    modeldf4 <- modeldf4 %>% rename(se="sqrt(diag(vcov(modeli)))")
    
    #merge everything together and create output
    modeldf5 <- merge(modeldf3, modeldf4, by="variable", all.x=TRUE)
    listmodeldfs[[i]] <- modeldf5
  }
  return(listmodeldfs)
}

#refit model on all i imputation sets after majority vote for subsequent Wald test: 
refitmodelglm <- function(impsets, formula){
  
  #store output
  models <- list()
  
  #loop over imputed datasets
  for (i in 1:length(impsets)){
    impset <- impsets[[i]] 
    full <- glm(formula, family=binomial, data=impset)
    models[[i]] <- full
  }
  return(models)
}

# Calculation of predicted values in each imputed set i, to later be able to calculate performance measures: 
predcalcGLM <- function(impsets, finalmodeldf){
  
  #store output 
  predimps <- list()
  
  #GLM provides the model output of binary variables as "variablex1", since this is not the same name als the names of the variables in the dataset, we remove the extra "1" from this column names: 
  finalmodeldf$variable <- sub("(.*)1$", "\\1", finalmodeldf$variable) 
  model_coefs <- with(finalmodeldf, setNames(Meancoef, variable))
  
  #create loop over imputed datasets, using the final model coefficients, multiplying them with the variable values. 
  for(i in 1:length(impsets)){
    datasetimp <- impsets[[i]]
    
    #we need to set everything to numeric 
    names <- datasetimp %>% select_if(is.factor)
    names2 <- colnames(names)
    impsetprep <- datasetimp %>% mutate_if(is.factor, as.numeric)
    impsetprep <- impsetprep %>% mutate_at(names2, ~ . -1)
    
    #multiply all model coefficients with their respective variable value, not using the intercept variable (-1). --> calculate linear predictor 
    impsetprep <- impsetprep %>% mutate(finalmodellp = model_coefs["(Intercept)"] + rowSums(across(all_of(names(model_coefs)[-1]),
                                                                                                   ~ get(cur_column()) * model_coefs[cur_column()]))) 
    
    #paste this value to the imputed dataset as predicted value. 
    impsetpred <- datasetimp
    impsetpred$lp <- impsetprep$finalmodellp
    impsetpred$pred <- 1/(1+exp(-impsetpred$lp))
    predimps[[i]] <- impsetpred
  }
  return(predimps)
}

#nullmodels to calculate IPA and rescaled Brier scores (using GLM) 
nullmodels <- function(impsets){
  
  #store output
  nullmodels <- list()
  
  #loop over imputed datasets
  for (i in 1:length(impsets)){
    set <- impsets[[i]]
    #create models
    null <- glm(outcome~1, data=set, family=binomial)
    nullmodels[[i]] <- null
  }
  return(nullmodels)
}


#function to get all desired performance measures for apparent performance of the apparent model. 
### We choose AUC, IPA, Brier, OE ratio and calibration intercept and slope. 

impPerfGLM <- function(impsetspred){
  
  #store output
  performancetablesimp <- list()
  
  #create a loop over the imputed datasets, that now contain predicted values (pred) of the final model and predicted values of the nullmodel (prednull)
  for(i in 1:length(impsetspred)){
    
    #calculate auc en logit_auc --> save the logit. We make sure we store these as a variable called performanceM for all calculations, so we can later merge these sets together. 
    totalpred <- impsetspred[[i]] 
    performanceimp <- compute_c_stat_log_imp(obs_outcome=totalpred$outcome, pred_outcome=totalpred$pred)
    c_statlogdf <- as.data.frame(performanceimp)
    c_statlogdf <- c_statlogdf %>% rownames_to_column("performanceM")
    
    #calculate other performance statistics using val.prob, storing them as performanceM variable in a separate dataframe. 
    totalpred$outcome <- as.numeric(totalpred$outcome)-1
    performanceimp <- val.prob(totalpred$pred, totalpred$outcome, pl=FALSE)
    performancedf <- data.frame(performanceimp)
    performancedf <- performancedf %>% rownames_to_column("performanceM")
    performancedf <- performancedf %>% filter(!is.na(performanceM))
    
    #calculate performance statistics for the null model in order to calculate IPA further on
    performancenull <- val.prob(totalpred$prednull, as.numeric(totalpred$outcome), pl=FALSE)
    performancenulldf <- as.data.frame(performancenull)
    performancenulldf <- performancenulldf %>% rownames_to_column("performanceM")
    performancenulldf <- performancenulldf %>% filter(!is.na(performanceM)) #dit moest ik erbij doen omdat hij om een of andere reden NA's kreeg, omdat er 1 variabele in kolom performanceM NA werd, en die nam hij mee met de Brier. 
    
    #now calculate the performance measures that are not returned by val.prop by hand (IPA, OE). We also transformed IPA when we were not yet sure how to pool IPA. 
    add.data.IPA <- data.frame(performanceM = "IPA", performanceimp = 1-((performancedf$performanceimp[performancedf$performanceM=="Brier"])/(performancenulldf$performancenull[performancenulldf$performanceM=="Brier"])))
    add.data.IPA2 <- data.frame(performanceM = "logitIPA", performanceimp=logit(as.numeric(add.data.IPA$performanceimp)))
    add.data.OE <- data.frame(performanceM = "OE", performanceimp =  (mean((as.numeric(totalpred$outcome))))/(as.numeric(mean(as.numeric(totalpred$pred)))))
    add.data.BrierR <- data.frame(performanceM="BrierR", performanceimp = 1 - (mean((totalpred$outcome - totalpred$pred)^2) / (mean(totalpred$outcome) * (1 - mean(totalpred$outcome)))))
    dfperformance2 <- bind_rows(performancedf, c_statlogdf, add.data.IPA, add.data.IPA2, add.data.OE, add.data.BrierR) #wrap it all toghether                          
    dfperformance3 <- dfperformance2 %>% filter(performanceM=="logitauc"|performanceM=="logit_se"|performanceM=="IPA"|performanceM=="logitIPA"|performanceM=="Slope" | performanceM=="Intercept" | performanceM=="Brier" | performanceM=="BrierR" | performanceM=="OE")
    performancetablesimp[[i]] <- dfperformance3
    
  }
  return(performancetablesimp)
}


# For the calibration plot (GLM specific part): 
# provide coordinates of the smooth calibration curve ----
# This function is different from Firth because for GLM, 1 should be subtracted from the outcome because R transformed our outcome variable to 1 and 2 instead of 0 and 1
draw_smooth_calibration <- function( obs_outcome, pred_outcome){
  
  y_input <- t( rbind( unlist( obs_outcome)))
  x_input <- t( rbind( unlist( pred_outcome)))
  
  y_coordinates <- loess( (as.numeric(y_input)-1) ~ x_input)$fitted # this is where the 1 is subtracted for the GLM function whereas for the Firth function, this is not necessary. 
  
  # store all coordinates of smooth curve (ordered based on x)
  smooth_cal_slope_info <- list( x_coordinates = x_input[order( x_input)],
                                 y_coordinates = y_coordinates[order( x_input)])
  
  
  # output overall info calibration slope
  return( smooth_cal_slope_info)
}


# Function for model development in GLM bootstrapped samples: 
modelglm <- function(bootstrapsets, formula){
  
  #store output
  modelsbsglm <- list()
  
  #create loop over imputation sets 
  for(i in 1:length(bootstrapsets)){
    
    #first extract list with 500 bootstraps from list with imputation sets 
    impcat <- bootstrapsets[[i]]
    
    #store output of following loop
    modelsglm <- list()
    
    #then create loop within that imputation set number, over the 500 bootstrap samples it contains 
    for(j in 1:length(impcat)){
      setbsglm <- impcat[[j]]  
      
      #create models using GLM
      fullmodel <- glm(formula, data=setbsglm, family=binomial)
      bwmodelglm <- step(fullmodel, direction="backward", trace=0)
      modelsglm[[j]] <- bwmodelglm
    }
    modelsbsglm[[i]] <- modelsglm
  }
  return(modelsbsglm)
}

# Function for model development in GLM bootstrapped samples: 
modelglm_pl <- function(bootstrapsets, formula, nr_cores = 1, verbose = TRUE){
  nrI <- length(bootstrapsets)
  nrB <- length(bootstrapsets[[1]])
  if(verbose) print(paste("nr of imputations:", nrI))
  if(verbose) print(paste("nr of bootstraps:", nrB))
  
  cl_parms <- data.frame(nB = 1:nrB)
  # cl_parms <- merge(cl_parms, data.frame(nI = 1:nrI))
  cl_parms <- apply(cl_parms, 1, as.list)
  
  # bootstrapsets <- lapply(cl_parms, function(ee) bootstrapsets[[ee$nI]][[ee$nB]])
  # print(pryr::object_size(bootstrapsets))
  # Sys.info()[["sysname"]] == "Windows" && 
  
  
  if(nr_cores > 1){
    cl_modelglm <- makePSOCKcluster(nr_cores)
    on.exit({stopCluster(cl_modelglm); print("Cluster connections closed")}, add = TRUE)
    clusterEvalQ(cl_modelglm, library(stats))
    clusterExport(cl_modelglm, varlist = "formula", envir = environment())
  }
  .tic <- Sys.time()
  modelsbsglm <- vector(mode = "list", length = nrI)
  on.exit(return(modelsbsglm), add = TRUE)
  for (ii in 1:nrI){
    if(nr_cores > 1){
      if(verbose) cat(paste("------------ Now at imputation dataset nr:", ii, "---------------\n"))
      
      modelsbsglm[[ii]] <- parLapplyLB(cl_modelglm, bootstrapsets[[ii]], function(ee){
        step(glm(formula, data=ee, family=binomial), direction="backward", trace=0)
      })
      .toc <- Sys.time()
      if(verbose) cat("Elapsed time so far:", format(difftime(.toc, .tic)), "\n")
      if(verbose) cat("ETA:", format(.tic + difftime(.toc + (nrI-ii)*difftime(.toc, .tic)/ii, .tic)), "\n")
      if(verbose) print(clusterEvalQ(cl_modelglm, gc()))
    } else {
      modelsbsglm[[ii]] <- map(.progress = TRUE, bootstrapsets[[ii]], function(ee){
        step(glm(formula, data=ee, family=binomial), direction="backward", trace=0)
      })
    }
  }
  
  return(modelsbsglm)
}

# calculate performance over the bootstrap models (GLM specific) 
impPerfBSglm <- function(imputatiesets, bootstrapsets, bsmodels, bsnullmodels){
  
  #store output as two seperate lists, one for performance in the bootstrap samples and one for performance in the original imputation sets. 
  performancetablesbs <- list()
  performancetablesimp <- list()
  
  
  #first create a loop length of imputation sets, creating the imputation sets, the list with their respective bootstrap samples (individual), list of models for each set and list of null models for each set.  
  for(i in 1:length(imputatiesets)){
    
    imputatieset <- imputatiesets[[i]]
    bootstraplist <- bootstrapsets[[i]]
    modelbslist <- bsmodels[[i]]
    nullmodelbslist <- bsnullmodels[[i]]
    
    
    imputatieset$outcome <- as.numeric(imputatieset$outcome)-1
    
    #then loop within the bootstrap samples lists, using the models and the nullmodels. 
    for(j in 1:length(bootstraplist)){
      bootstrapdf <- bootstraplist[[j]]
      modelbs <- modelbslist[[j]]
      nullmodelbs <- nullmodelbslist[[j]]
      bootstrapdf$pred <- predict(modelbs, type="response", newdata=bootstrapdf)
      
      #auc bootstrapmodel in bootstrapsample
      bootstrap <- compute_c_stat_log_imp(obs_outcome=bootstrapdf$outcome, pred_outcome=bootstrapdf$pred)
      c_statlogdf <- as.data.frame(bootstrap)
      c_statlogdf <- c_statlogdf %>% rownames_to_column("performanceM")
      
      #use null models for calculating IPA (bootstrapmodel in bootstrap sample)
      bootstrapdf$outcome <- as.numeric(bootstrapdf$outcome)-1
      bootstrapdf$prednull <- predict(nullmodelbs, type="response", newdata=bootstrapdf)
      null <- val.prob(bootstrapdf$prednull, as.numeric(bootstrapdf$outcome), pl=FALSE)
      performancenulldf <- as.data.frame(null)
      performancenulldf <- performancenulldf %>% rownames_to_column("performanceM")
      performancenulldf <- performancenulldf %>% filter(!is.na(performanceM))
      
      #andere performance measures bootstrapmodels in bootstrapsample 
      # first using val.prob
      bootstrap <- val.prob(bootstrapdf$pred, bootstrapdf$outcome, pl=FALSE)
      bootstrapdf2 <- data.frame(bootstrap)
      bootstrapdf2 <- bootstrapdf2 %>% rownames_to_column("performanceM")
      
      #then calculating IPA and O:E ratio. Storing all these measures as performanceM as well in order to later merge all the datasets. 
      add.data.IPAbs <- data.frame(performanceM = "IPA" , bootstrap = (1-((bootstrapdf2$bootstrap[bootstrapdf2$performanceM=="Brier"])/(performancenulldf$null[performancenulldf$performanceM=="Brier"])))) 
      add.dataOE <- data.frame(performanceM = "OE", bootstrap = as.numeric((mean(as.numeric(bootstrapdf$outcome)))/(mean(bootstrapdf$pred))))
      add.data.BrierR <- data.frame(performanceM="BrierR", bootstrap = 1 - (mean((bootstrapdf$outcome - bootstrapdf$pred)^2) / (mean(bootstrapdf$outcome) * (1 - mean(bootstrapdf$outcome)))))
      bootstrapdf3 <- bind_rows(bootstrapdf2, c_statlogdf, add.data.IPAbs, add.dataOE, add.data.BrierR) 
      
      #merge everything together 
      bootstrapdf4 <- bootstrapdf3 %>% filter(performanceM=="C (ROC)"|performanceM=="Brier" | performanceM=="BrierR" | performanceM=="logitauc"|performanceM=="logit_se"|performanceM=="IPA"| performanceM=="Slope" | performanceM=="Intercept" | performanceM=="OE") 
      
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
      imp <- val.prob(imputatieset$pred, imputatieset$outcome, pl=FALSE)
      impdf <- data.frame(imp)
      impdf <- impdf %>% rownames_to_column("performanceM")
      
      #performance measures calculated by hand (IPA and O:E)
      add.data.IPAimp <- data.frame(performanceM = "IPA" , imp=(1-((impdf$imp[impdf$performanceM=="Brier"])/(impnulldf$null[impnulldf$performanceM=="Brier"])))) 
      add.dataOEimp <- data.frame(performanceM = "OE", imp = as.numeric(mean(imputatieset$outcome))/(mean(imputatieset$pred)))
      add.data.BrierR <- data.frame(performanceM="BrierR", imp = 1 - (mean((imputatieset$outcome - imputatieset$pred)^2) / (mean(imputatieset$outcome) * (1 - mean(imputatieset$outcome)))))
      impdf2 <- bind_rows(impdf, c_statlogdfimp, add.data.IPAimp, add.dataOEimp, add.data.BrierR) 
      
      #merge everything together              
      impdf3 <- impdf2 %>% filter(performanceM=="C (ROC)"|performanceM=="logitauc"|performanceM=="logit_se"| performanceM=="Brier" | performanceM=="BrierR" | performanceM=="IPA" |performanceM=="Slope" | performanceM=="Intercept" | performanceM== "OE")
      
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
