### 5.3. APPARENT MODEL - FIRTH ###

set.seed(21212)

#Create i models for each imputed dataset i (in our example: 10 models developed in 10 datasets): 
modelsF <- backwardimpF(impsets=Final_impF, formula=outcome~age+
                         open+
                         T34+
                         dummy_Neotx1+
                         transhiatal+
                         fev1_compl+
                         dummy_Hblow+
                         comorb_dm+
                         gender+
                         eGFR+
                         dummy_ASA34+
                         dummy_Hbhigh+
                         comorb_cardiovasc+
                         comorb_hypertension+
                         smoking+
                         dummy_Neotx3+
                         bmi+
                         tiff_compl)

#get the modeldata (variables selected and their coefficients and standard errors) from the list with models 
modeldata <- extractmodelsF(modelsF)

#first create long dataframe 
modeldatamerged <- bind_rows(modeldata)

#select the final model variables and coefficients using the majority vote: 
finalmodelF <- majorityvote(modeldatamerged, impN=10, freq=0.5) 

#refit the models over all imputed datasets i with these selected variables: 
variablesrefit <- setdiff(finalmodelF$variable, "(Intercept)")
formula_string <- paste("outcome ~", paste(variablesrefit, collapse= "+"))

#now refit the model using the variables selected by the majority vote: 
refitmodels <- refitmodelF(impsets=Final_impF, formula=formula_string)

#extract the refitted models as above: 
modeldatarefit <- extractmodelsF(refitmodels)
modeldatarefitmerged <- bind_rows(modeldatarefit)

#summarise the results: 
finalmodelrefit <- modeldatarefitmerged %>% group_by(variable) %>% summarise(Meancoef=mean(modeli.coefficients))

#calculate the variance of the coefficients using Rubin's Rules (see 2. Function library)
varsrefit <- modeldatarefitmerged %>% group_by(variable) %>% summarise(rubins_rules_var(estimates=modeli.coefficients, ses=se, n_imputed_sets=10)) %>% rename(var = "rubins_rules_var(estimates = modeli.coefficients, ses = se, n_imputed_sets = 10)")

#merge these vars with the final model
finalmodelrefit <- merge(finalmodelrefit, varsrefit, by="variable", all.x=TRUE)

#check if variables need to be dropped using Wald (using p > 0.05). 
#in case multiple variables are p > 0.5, consider dropping the least significant variable and then repeating these steps until no further variables can be removed
#calculate Wald statistic: 
finalmodelrefit$W <- (finalmodelrefit$Meancoef^2)/(finalmodelrefit$var)
#calculate p-value based on the Wald statistic: 
finalmodelrefit$p <- 1-pchisq(finalmodelrefit$W, df=1)

#check p-values: 
finalmodelrefit

#in our example, fev1_compl can be dropped from the model --> repeat refit. 
variablesrefit2 <- setdiff(finalmodelrefit$variable, c("(Intercept)", "fev1_compl"))
formula_string2 <- paste("outcome ~", paste(variablesrefit2, collapse= "+"))

refitmodels2 <- refitmodelF(impsets=Final_impF, formula=formula_string2)

modeldatarefit2 <- extractmodelsF(refitmodels2)
modeldatarefitmerged2 <- bind_rows(modeldatarefit2)

finalmodelrefit2 <- modeldatarefitmerged2 %>% group_by(variable) %>% summarise(Meancoef=mean(modeli.coefficients))
finalmodelrefit2$OR <- exp(finalmodelrefit2$Meancoef) #add odds ratio

varsrefit2 <- modeldatarefitmerged2 %>% group_by(variable) %>% summarise(rubins_rules_var(estimates=modeli.coefficients, ses=se, n_imputed_sets=10)) %>% rename(var = "rubins_rules_var(estimates = modeli.coefficients, ses = se, n_imputed_sets = 10)")

#merge these vars with the final model
finalmodelrefit2 <- merge(finalmodelrefit2, varsrefit2, by="variable", all.x=TRUE)

#check if variables need to be dropped using Wald (using p > 0.05). 
# in case multiple variables are p > 0.5, consider dropping the least significant variable and then repeating these steps until no further variables can be removed
finalmodelrefit2$W <- (finalmodelrefit2$Meancoef^2)/(finalmodelrefit2$var)
finalmodelrefit2$p <- 1-pchisq(finalmodelrefit2$W, df=1)

#check p-values: 
finalmodelrefit2
# all are p<0.05, this is the final model. 

finalmodelF <- finalmodelrefit2

# For Firth, we calculate the DF by using the total sample size and subtract the number of predictors, -1 --> 619 - 6 - 1 = 612 
dfF <- 612

finalmodelF$se <- sqrt(finalmodelF$var)
finalmodelF$lower <- finalmodelF$Meancoef - (qt(0.975, dfF)*finalmodelF$se)
finalmodelF$upper <- finalmodelF$Meancoef + (qt(0.975, dfF)*finalmodelF$se)
finalmodelF$OR <- exp(finalmodelF$Meancoef)
finalmodelF$lowerOR <- exp(finalmodelF$Meancoef - (qt(0.975, dfF)*finalmodelF$se))
finalmodelF$upperOR <- exp(finalmodelF$Meancoef + (qt(0.975, dfF)*finalmodelF$se))


#change order of variables to be able to compare model output: 
finalmodelF <- finalmodelF %>% select(variable, Meancoef, se, lower, upper, OR, lowerOR, upperOR, p)
finalmodelF


saveRDS(finalmodelF, paste0("finalmodelF", ".rds"))


