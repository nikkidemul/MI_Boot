### 5.3. APPARENT MODEL - FIRTH ###

set.seed(21212)

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

#get the modeldata (coefficients and standard errors) from the list with models 
modeldata <- extractmodels(modelsF)

modeldatamerged <- bind_rows(modeldata)

#final model: 
finalmodelF <- majorityvote(modeldatamerged, impN=10, freq=0.5) 

#make transformation of column names 
modeldatamerged2 <- as.data.frame(modeldatamerged %>% rename(coef = "modeli.coefficients"))

# do some transformations to be able to merge the standard errors with the final model dataframe 
standarderrorsF <- modeldatamerged2 %>% group_by(variable) %>% summarise(rubins_rules_var(estimates=coef, ses=se, n_imputed_sets=10))
standarderrorsF <- standarderrorsF %>% rename(ses = "rubins_rules_var(estimates = coef, ses = se, n_imputed_sets = 10)")

#merge with final model
finalmodelF$OR <- exp(finalmodelF$Meancoef)
finalmodelF2 <- merge(finalmodelF, standarderrorsF, by="variable", all.x=TRUE)

#calculate 95% CI 
finalmodel3F <- finalmodelF2
finalmodel3F$lower <- finalmodel3F$Meancoef - 1.96*sqrt(finalmodel3F$ses)
finalmodel3F$upper <- finalmodel3F$Meancoef + 1.96*sqrt(finalmodel3F$ses)
finalmodel3F$lowerOR <- exp(finalmodel3F$lower)
finalmodel3F$upperOR <- exp(finalmodel3F$upper)

#change order of variables to be able to compare model output: 
finalmodel3F <- finalmodel3F %>% select(variable, Meancoef, ses, lower, upper, OR, lowerOR, upperOR)

saveRDS(finalmodel3F,file=paste0("FinalmodelFirth",".rds"))

#check model
finalmodel3F


