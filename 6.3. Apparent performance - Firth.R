### 6.3. Apparent performance - Firth ### 

#get the predicted values by final model for each imputed dataset 
impsetspredF <- predcalcF(Final_impF, finalmodelF)

#develop the nullmodels 
nullmodelsF <- nullmodels(Final_impF) #for logistic regression models, the output should be conform the prevalence of the outcome. 

#extract nullmodels
nulldataF <- extractmodels(nullmodelsF)

nulldatamergeF <- bind_rows(nulldataF)

#reuse this script to pool the estimates 
nullmodelimpF <- majorityvote(nulldatamergeF, impN=10, freq=0.5)  

#include the predicted value of the nullmodels in the datasets
completeF <- predcalcnull(impsetspredF, nullmodelimpF) 

#calculate the performance measures
performanceimpF <- impPerfF(completeF) 

#stack the sets, pool (take average) and transform back.
performancelongF <- bind_rows(performanceimpF)
performancesummaryF <- performancelongF %>% group_by(performanceM)
performancesummary2F <- performancesummaryF %>% summarise(Meanperformance = mean(performanceimp), .groups="drop")
performancesummary2F$nonlogmean <- plogis(performancesummary2F$Meanperformance)

#now we have a table both including the logit mean and the mean of the original values. The only performance measure that we wanted to pool as logit was AUC, therefore we create a new table where all performance measures are in their right format: 
finalperformanceF <- performancesummary2F %>% select(performanceM, Meanperformance)
finalperformanceF$Meanperformance[finalperformanceF$performanceM=="logitauc"] <- performancesummary2F$nonlogmean[performancesummary2F$performanceM=="logitauc"]

finalperformanceF <- finalperformanceF %>% filter(performanceM != "logitIPA")
finalperformanceF <- finalperformanceF %>% filter( performanceM != "logit_se")
finalperformanceF$performanceM <- ifelse(finalperformanceF$performanceM=="logitauc", "AUC", finalperformanceF$performanceM)

#calculate total variance for confidence intervals AUC (later stage): 
total_var_performance_impF <- rubins_rules_var_auc(estimates=performancelongF$performanceimp[performancelongF$performanceM=="logitauc"], ses=performancelongF$performanceimp[performancelongF$performanceM=="logit_se"], n_imputed_sets=10) #is opgeslagen als totale variantie. 

total_var_performance_impF #0.0106951

#Calibration plot
outcomestacked <- lapply(impsetspredF, function(x) x%>%select(outcome))
predstacked <- lapply(impsetspredF, function(x) x%>%select(pred))

smooth_info <- draw_smooth_calibrationF(obs_outcome=outcomestacked, pred_outcome=predstacked)
calibrationplot <- plot_smooth_calibration(x_coordinates=smooth_info$x_coordinates, y_coordinates=smooth_info$y_coordinates) 

#Save results: 
saveRDS(finalperformanceF, file=paste0("FinalperformanceF", ".rds"))

finalmodel3F

#ROC Curve: 
Myimp$lpF <- -1.81193367 + 0.04323434*Myimp$age + 0.42775800*as.numeric(Myimp$comorb_hypertension) + 
  -0.01057415*Myimp$fev1_compl + -0.45562568*as.numeric(Myimp$comorb_dm) + -0.52713502*as.numeric(Myimp$transhiatal) + 
  0.64277993*as.numeric(Myimp$open) + 0.39235827*as.numeric(Myimp$dummy_Neotx3)

Myimp$predF <- 1/(1+exp(-Myimp$lpF))
summary(Myimp$predF)

myROC <- roc(Myimp$outcome ~ Myimp$lpF, levels=c(0,1))
plot(myROC, xlim=c(1,0)) 
auc(myROC)



