### 6.3. Apparent performance - Firth's penalized logistic regression ### 

# In this script, we will use the model obtained in 5.3. (apparent model), and test its performance in the original data (the i imputed datasets). 
### Check out 2. and 2.3. Function libraries for the code on the specific functions used. 

# First, calculate the predicted values of the final model for each imputed dataset
impsetspredF <- predcalcF(Final_impF, finalmodelF)

# Because we want to calculate the index of prediction accuracy and the rescaled Brier score, we also have to develop null models.
# Null models are simply the prevalence of the outcome for logistic regression models. 
# Develop null models: 
nullmodelsF <- nullmodels(Final_impF) 

# Extract null models from the output generated above: 
nulldataF <- extractmodelsF(nullmodelsF)

# Summarise the results into a long format dataframe: 
nulldatamergeF <- bind_rows(nulldataF)

# Summarise the results into a final null model: 
### Strictly we don't need to do a majority vote here, but we do want to summarise the results into one final model, for simplicity we therefore used the same function: 
nullmodelimpF <- majorityvote(nulldatamergeF, impN=length(Final_impF), freq=0.5)  

# Create new datasets, including the predicted value of the final model, and the predicted values of the null models: 
completeF <- predcalcnull(impsetspredF, nullmodelimpF) 

# Calculate the desired performance measures: 
### In our exalmple we calculate AUC, IPA, Brier, rescaled Brier, O:E ratio and calibration intercept and slope. 
performanceimpF <- impPerfF(completeF) 

# Then we have to pool the performance measures into one final model performance. For point estimates, taking the average suffices. 
### For the AUC: we have used a logit transformation to pool (because it is a log-rank statistic). After pooling, we can transform back. 
performancelongF <- bind_rows(performanceimpF)
performancesummaryF <- performancelongF %>% group_by(performanceM)
performancesummary2F <- performancesummaryF %>% summarise(Meanperformance = mean(performanceimp), .groups="drop")
performancesummary2F$nonlogmean <- plogis(performancesummary2F$Meanperformance)

# Now we have a table that includes both the logit mean and the mean of the original values. 
# The only performance measure for which we needed the logit was the AUC. Therefore we create a new table with only the final measures we need: 
finalperformanceF <- performancesummary2F %>% select(performanceM, Meanperformance)
finalperformanceF$Meanperformance[finalperformanceF$performanceM=="logitauc"] <- performancesummary2F$nonlogmean[performancesummary2F$performanceM=="logitauc"]

finalperformanceF <- finalperformanceF %>% filter(performanceM != "logitIPA")
finalperformanceF <- finalperformanceF %>% filter( performanceM != "logit_se")
finalperformanceF$performanceM <- ifelse(finalperformanceF$performanceM=="logitauc", "AUC", finalperformanceF$performanceM)

# Then, we also have to pool the standard errors of the logit AUC.
# We save this calculated variance of the logit AUC to later calculate the confidence intervals surrounding the AUC (for both apparent and adjusted performance). 
total_var_performance_impF <- rubins_rules_var_auc(estimates=performancelongF$performanceimp[performancelongF$performanceM=="logitauc"], ses=performancelongF$performanceimp[performancelongF$performanceM=="logit_se"], n_imputed_sets=length(Final_impF))  

if(MIb.intrnl.verbose) total_var_performance_impF #0.01037926

# Make a calibration plot: 
outcomestacked <- lapply(impsetspredF, function(x) x%>%select(outcome))
predstacked <- lapply(impsetspredF, function(x) x%>%select(pred))

smooth_info <- draw_smooth_calibrationF(obs_outcome=outcomestacked, pred_outcome=predstacked)
calibrationplot <- plot_smooth_calibration(x_coordinates=smooth_info$x_coordinates, y_coordinates=smooth_info$y_coordinates) 

#Save results: 
#saveRDS(finalperformanceF, file=paste0("FinalperformanceF", ".rds"))

# Make a ROC: 
if(MIb.intrnl.verbose) print(finalmodelF)

# For the ROC we can calculate the predicted values manually as follows: 
model_coefs <- with(finalmodelF, setNames(Meancoef, variable))
# to numeric: 
Myimp <- Myimp %>% mutate_if(is.factor, as.numeric)

#multiply all model coefficients with their respective variable value, not using the intercept variable (-1)
Myimp <- Myimp %>% mutate(lpF = model_coefs["(Intercept)"] + rowSums(across(all_of(names(model_coefs)[-1]),
                                                                            ~ get(cur_column()) * model_coefs[cur_column()]))) 

#paste this value to the imputed dataset as predicted value 
Myimp$predF <- 1/(1+exp(-Myimp$lpF))

if(MIb.verbose) cat("\nA summary of ...:\n");print(summary(Myimp$predF))

if(MIb.verbose){
  myROC <- roc(Myimp$outcome ~ Myimp$lpF, levels=c(0,1))
  par(pty='s'); plot(myROC, xlim=c(1,0)); par(pty='m')
  auc(myROC)
}


