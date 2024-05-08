### 5.2. Apparent model - GLM (manually) ### 

# Of note: in this dataset, R has made all binary variables 1 and 2 instead of 0 and 1.
# It is important to keep this in mind when analyzing the data. 

#Create 10 models for each imputed dataset: 
modelsglm <- backwardimpglm(impsets=Final_imp, formula=outcome~age
                            +bmi
                            +gender
                            +open
                            +transhiatal
                            +dummy_Neotx1
                            +dummy_Neotx3
                            +T34
                            +comorb_dm
                            +comorb_cardiovasc
                            +comorb_hypertension
                            +smoking
                            +dummy_ASA34
                            +eGFR
                            +dummy_Hblow
                            +dummy_Hbhigh
                            +fev1_compl
                            +tiff_compl)

#function to extract the models from the above generated output
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

#get this modeldata (coefficients and standard errors) from the list with models 
modeldata <- extractmodels(modelsglm)

#first create long dataframe 
modeldatamerged <- bind_rows(modeldata)

#final model: 
finalmodelimpglm <- majorityvote(modeldatamerged, impN=10, freq=0.5) 

#make transformation of column names 
modeldatamerged2 <- as.data.frame(modeldatamerged %>% rename(coef = "modeli.coefficients"))

# do some transformations to be able to merge the standard errors with the final model dataframe 
standarderrors <- modeldatamerged2 %>% group_by(variable) %>% summarise(rubins_rules_var(estimates=coef, ses=se, n_imputed_sets=10))
standarderrors <- standarderrors %>% rename(ses = "rubins_rules_var(estimates = coef, ses = se, n_imputed_sets = 10)")

#merge with final model
finalmodelimpglm$OR <- exp(finalmodelimpglm$Meancoef) 
finalmodelimpglm2 <- merge(finalmodelimpglm, standarderrors, by="variable", all.x=TRUE) 

#calculate 95% CI 
finalmodelimpglm3 <- finalmodelimpglm2 
finalmodelimpglm3$lower <- finalmodelimpglm3$Meancoef - 1.96*sqrt(finalmodelimpglm3$ses)
finalmodelimpglm3$upper <- finalmodelimpglm3$Meancoef + 1.96*sqrt(finalmodelimpglm3$ses)
finalmodelimpglm3$lowerOR <- exp(finalmodelimpglm3$lower)
finalmodelimpglm3$upperOR <- exp(finalmodelimpglm3$upper)

#change order of variables to compare with psfmi 
finalmodelimpglm3 <- finalmodelimpglm3 %>% select(variable, Meancoef, ses, lower, upper, OR, lowerOR, upperOR)

saveRDS(finalmodelimpglm3, file=paste0("FinalmodelGLM", ".rds"))

#check model output: 
finalmodelimpglm3




