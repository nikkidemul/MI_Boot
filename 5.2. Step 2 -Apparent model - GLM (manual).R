### 5.2. Apparent model - Non-penalized logistic regression (GLM - manually) ### 

#In this script, we develop the apparent model over the i imputed datasets. 
### Check out 2. and 2.2. Function libraries for the code on the specific functions used. 

### Of note: in this dataset, R has will name all binary variables 1 and 2 instead of 0 and 1 after numerical transformation.
###### It is important to keep this in mind when analyzing the data. 

set.seed(21212)
#Create i models for each imputed dataset i (in our example: 10 models developed in 10 datasets): 
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

#get the modeldata (variables selected and their coefficients and standard errors) from the list with models 
modeldata <- extractmodels(modelsglm)

#first create long dataframe 
modeldatamerged <- bind_rows(modeldata)

#select the final model variables and coefficients using the majority vote: 
finalmodelglm <- majorityvote(modeldatamerged, impN=10, freq=0.5) 

#refit the models over all imputed datasets i with these selected variables: 
finalmodelglm$variable <- sub("(.*)1$", "\\1", finalmodelglm$variable) #glm puts 1 in the variable name for each binary variable, 
### in order to use these names for the refit, we need to remove this '1' from the variable name. 
variablesrefitglm <- setdiff(finalmodelglm$variable, "(Intercept)") #intercept is not a variable name - we remove these from the variable string. 
formula_stringglm <- paste("outcome ~", paste(variablesrefitglm, collapse= "+"))

#now refit the model using the variables selected by the majority vote: 
refitmodelsglm <- refitmodelglm(impsets=Final_imp, formula=formula_stringglm)

#extract the refitted models as above: 
modeldatarefitglm <- extractmodels(refitmodelsglm)
modeldatarefitmergedglm <- bind_rows(modeldatarefitglm)

#summarise the results: 
finalmodelrefitglm <- modeldatarefitmergedglm %>% group_by(variable) %>% summarise(Meancoef=mean(modeli.coefficients))

#calculate the variance of the coefficients using Rubin's Rules (see 2. Function library)
varsrefitglm <- modeldatarefitmergedglm %>% group_by(variable) %>% summarise(rubins_rules_var(estimates=modeli.coefficients, ses=se, n_imputed_sets=10)) %>% rename(var = "rubins_rules_var(estimates = modeli.coefficients, ses = se, n_imputed_sets = 10)")

#merge these vars with the final model
finalmodelrefitglm <- merge(finalmodelrefitglm, varsrefitglm, by="variable", all.x=TRUE)

#check if variables need to be dropped using Wald (using p > 0.05). 
#in case multiple variables are p > 0.5, consider dropping the least significant variable and then repeating these steps until no further variables can be removed
#calculate Wald statistic: 
finalmodelrefitglm$W <- (finalmodelrefitglm$Meancoef^2)/(finalmodelrefitglm$var)
#calculate p-value based on the Wald statistic: 
finalmodelrefitglm$p <- 1-pchisq(finalmodelrefitglm$W, df=1)

#check p-values: 
finalmodelrefitglm

#in our example, fev1_compl can be dropped from the model --> repeat refit. 
finalmodelrefitglm$variable <- sub("(.*)1$", "\\1", finalmodelrefitglm$variable) 
variablesrefit2glm <- setdiff(finalmodelrefitglm$variable, c("(Intercept)", "fev1_compl"))
formula_string2glm <- paste("outcome ~", paste(variablesrefit2glm, collapse= "+"))

refitmodels2glm <- refitmodelglm(impsets=Final_impF, formula=formula_string2glm)

modeldatarefit2glm <- extractmodels(refitmodels2glm)
modeldatarefitmerged2glm <- bind_rows(modeldatarefit2glm)

finalmodelrefit2glm <- modeldatarefitmerged2glm %>% group_by(variable) %>% summarise(Meancoef=mean(modeli.coefficients))

varsrefit2glm <- modeldatarefitmerged2glm %>% group_by(variable) %>% summarise(rubins_rules_var(estimates=modeli.coefficients, ses=se, n_imputed_sets=10)) %>% rename(var = "rubins_rules_var(estimates = modeli.coefficients, ses = se, n_imputed_sets = 10)")

#merge these vars with the final model
finalmodelrefit2glm <- merge(finalmodelrefit2glm, varsrefit2glm, by="variable", all.x=TRUE)

#check if variables need to be dropped using Wald (using p > 0.05). 
# in case multiple variables are p > 0.5, consider dropping the least significant variable and then repeating these steps until no further variables can be removed
finalmodelrefit2glm$W <- (finalmodelrefit2glm$Meancoef^2)/(finalmodelrefit2glm$var)
finalmodelrefit2glm$p <- 1-pchisq(finalmodelrefit2glm$W, df=1)

#check p-values: 
finalmodelrefit2glm
# all are p<0.05, this is the final model. 

finalmodelglm <- finalmodelrefit2glm

#calculate 95% CI 
finalmodelglm$lower <- finalmodelglm$Meancoef - 1.96*sqrt(finalmodelglm$var)
finalmodelglm$upper <- finalmodelglm$Meancoef + 1.96*sqrt(finalmodelglm$var)
finalmodelglm$OR <- exp(finalmodelglm$Meancoef) 
finalmodelglm$lowerOR <- exp(finalmodelglm$lower)
finalmodelglm$upperOR <- exp(finalmodelglm$upper)
finalmodelglm$se <- sqrt(finalmodelglm$var)

#change order of variables to be able to compare model output: 
finalmodelglm <- finalmodelglm %>% select(variable, Meancoef, se, lower, upper, OR, lowerOR, upperOR, p)
finalmodelglm


