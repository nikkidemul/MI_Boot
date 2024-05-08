### 6.2. Apparent performance GLM model ### 

#get the predicted values of the final model for each imputed dataset 
impsetspredGLM <- predcalcGLM(Final_imp, finalmodelimpglm3)

# Because we want to calculate the index of prediction accuracy (or rescaled Brier score for that matter) - we also have to develop nullmodels. 
#develop null models 
nullmodelsGLM <- nullmodels(Final_imp) #for logistic regression models, the output should be conform the prevalence of the outcome. 

#extract nullmodels 
nulldataGLM <- extractmodels(nullmodelsGLM)

# summarise the results: 
nulldatamergeGLM <- bind_rows(nulldataGLM)

#strictly we don't need to do a majority vote but we do want to summarise the results, for simplicity we use the same function: 
nullmodelimpGLM <- majorityvote(nulldatamergeGLM, impN=10, freq=0.5) 

#get new datasets, including the predicted value of the final model and the predicted values of the nullmodel. 
completeGLM <- predcalcnull(impsetspredGLM, nullmodelimpGLM) 

#Get the desired performance measures. 
#In our example we calculate AUC, IPA, Brier, BrierR, OE ratio and calibration intercept and slope. 
performanceimpGLM <- impPerfGLM(completeGLM) 

#stack de sets, take the average and then transform back 
performancelongGLM <- bind_rows(performanceimpGLM)
performancesummaryGLM <- performancelongGLM %>% group_by(performanceM)
performancesummary2GLM <- performancesummaryGLM %>% summarise(Meanperformance = mean(performanceimp), .groups="drop")
performancesummary2GLM$nonlogmean <- plogis(performancesummary2GLM$Meanperformance)

#now we have a table both including the logit mean and the mean of the original values. The only performance measure that we wanted to pool as logit was AUC, therefore we create a new table where all performance measures are in their right format: 

finalperformanceGLM <- performancesummary2GLM %>% select(performanceM, Meanperformance)
finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="logitauc"] <- performancesummary2GLM$nonlogmean[performancesummary2GLM$performanceM=="logitauc"]

finalperformanceGLM <- finalperformanceGLM %>% filter(performanceM != "logitIPA")
finalperformanceGLM <- finalperformanceGLM %>% filter( performanceM != "logit_se")
finalperformanceGLM$performanceM <- ifelse(finalperformanceGLM$performanceM=="logitauc", "AUC", finalperformanceGLM$performanceM)


# We have to pool the SE's of the AUC (logit). 
# We save these standard errors to later calculate the confidence intervals surrouding the AUC for both apparent performance and bootstrap performance. 

total_var_performance_impGLM <- rubins_rules_var_auc(estimates=performancelongGLM$performanceimp[performancelongGLM$performanceM=="logitauc"], ses=performancelongGLM$performanceimp[performancelongGLM$performanceM=="logit_se"], n_imputed_sets=10) 

total_var_performance_impGLM #0.01070087

# Calibration plot: 
outcomestacked <- lapply(impsetspredGLM, function(x) x%>%select(outcome))
predstacked <- lapply(impsetspredGLM, function(x) x%>%select(pred))

smooth_info <- draw_smooth_calibration(obs_outcome=outcomestacked, pred_outcome=predstacked)

calibrationplot <- plot_smooth_calibration(x_coordinates=smooth_info$x_coordinates, y_coordinates=smooth_info$y_coordinates) 

# For the ROC we can calculate the predicted values manually
Myimp$lpGLM <- -1.84368780 + 0.04389563*Myimp$age + 0.43335998*as.numeric(Myimp$comorb_hypertension) + 
  -0.01082717*Myimp$fev1_compl + -0.45949801*as.numeric(Myimp$comorb_dm) + 
  -0.53482468*as.numeric(Myimp$transhiatal) + 0.65965965*as.numeric(Myimp$open) + 
  0.39556042*as.numeric(Myimp$dummy_Neotx3)

Myimp$predGLM <- 1/(1+exp(-Myimp$lpGLM))
summary(Myimp$predGLM)
myROC <- roc(Myimp$outcome ~Myimp$lpGLM, levels=c(0,1))
plot(myROC, xlim=c(1,0))
auc(myROC)
