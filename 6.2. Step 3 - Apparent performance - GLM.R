### 6.2. Apparent performance non-penalized logistic regression (GLM) ### 

#In this script, we will use the model obtained in 5.2. (apparent model), and test its performance in the original data (the i imputed datasets). 
### Check out 2. and 2.2. Function libraries for the code on the specific functions used. 
#First, calculate the predicted values of the final model for each imputed dataset 
impsetspredGLM <- predcalcGLM(Final_imp, finalmodelglm)

# Because we want to calculate the index of prediction accuracy (or rescaled Brier score) - we also have to develop null models. 
# Null models are simply the event rate for binary outcomes (modelling: outcome ~ 1). 
# Develop null models:
nullmodelsGLM <- nullmodels(Final_imp) #for logistic regression models, the output should be conform the prevalence of the outcome. 

# Extract nullmodels from the output generated above: 
nulldataGLM <- extractmodels(nullmodelsGLM)

# Summarise the results into a long format dataframe: 
nulldatamergeGLM <- bind_rows(nulldataGLM)

# Summarise the results into a final null model:
### Strictly we don't need to do a majority vote but we do want to summarise the results into one final model, for simplicity we therefore used the same function: 
nullmodelimpGLM <- majorityvote(nulldatamergeGLM, impN=10, freq=0.5) 

# Create new datasets, including the predicted value of the final model and the predicted values of the nullmodel:  
completeGLM <- predcalcnull(impsetspredGLM, nullmodelimpGLM) 

# Calculate the desired performance measures. 
### In our example we calculate AUC, IPA, Brier, BrierR, OE ratio and calibration intercept and slope. 
performanceimpGLM <- impPerfGLM(completeGLM) 

# Then we have to pool the performance measures into one final model performance. For point estimates, taking the average suffices. 
### For the AUC: we have used a logit transformation to pool. After pooling, we can transform back. 
performancelongGLM <- bind_rows(performanceimpGLM)
performancesummaryGLM <- performancelongGLM %>% group_by(performanceM)
performancesummary2GLM <- performancesummaryGLM %>% summarise(Meanperformance = mean(performanceimp), .groups="drop")
performancesummary2GLM$nonlogmean <- plogis(performancesummary2GLM$Meanperformance)

# Now we have a table that includes both the logit mean and the mean of the original values. 
# The only performance measure for which we needed the logit was the AUC. Therefore we create a new table with only the final measures we need: 
finalperformanceGLM <- performancesummary2GLM %>% select(performanceM, Meanperformance)
finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="logitauc"] <- performancesummary2GLM$nonlogmean[performancesummary2GLM$performanceM=="logitauc"]

finalperformanceGLM <- finalperformanceGLM %>% filter(performanceM != "logitIPA")
finalperformanceGLM <- finalperformanceGLM %>% filter( performanceM != "logit_se")
finalperformanceGLM$performanceM <- ifelse(finalperformanceGLM$performanceM=="logitauc", "AUC", finalperformanceGLM$performanceM)

# Then, we also have to pool the standard error's of the logit AUC. 
# We save this calculated variance of the logit AUC to later calculate the confidence intervals surrounding the AUC (for both apparent and adjusted performance): 
total_var_performance_impGLM <- rubins_rules_var_auc(estimates=performancelongGLM$performanceimp[performancelongGLM$performanceM=="logitauc"], ses=performancelongGLM$performanceimp[performancelongGLM$performanceM=="logit_se"], n_imputed_sets=10) 

total_var_performance_impGLM #0.01037926

# Make a calibration plot:  
outcomestacked <- lapply(impsetspredGLM, function(x) x%>%select(outcome))
predstacked <- lapply(impsetspredGLM, function(x) x%>%select(pred))

smooth_info <- draw_smooth_calibration(obs_outcome=outcomestacked, pred_outcome=predstacked) 

calibrationplot <- plot_smooth_calibration(x_coordinates=smooth_info$x_coordinates, y_coordinates=smooth_info$y_coordinates) 

# Make a ROC: 
finalmodelglm 
# For the ROC we can calculate the predicted values manually as follows: 
Myimp$lpGLM <- -2.72780665 + 0.04244476*Myimp$age + 0.47924375*as.numeric(Myimp$comorb_hypertension) + 
 -0.41179816*as.numeric(Myimp$comorb_dm) + 
  -0.45425622*as.numeric(Myimp$transhiatal) + 0.65301430*as.numeric(Myimp$open) + 
  0.40846314*as.numeric(Myimp$dummy_Neotx3)

Myimp$predGLM <- 1/(1+exp(-Myimp$lpGLM))
summary(Myimp$predGLM)
myROC <- roc(Myimp$outcome ~Myimp$lpGLM, levels=c(0,1))
plot(myROC, xlim=c(1,0))
auc(myROC)
