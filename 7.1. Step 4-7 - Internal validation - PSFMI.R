### 7.1. PSFMI VALIDATION ### 

# For PSFMI, we take the final model from script 5.1. 
modelD1_back$RR_model_final

# We select the variables from this model to do internal validation of that specific model. 
# The selected variables were: age, open, transhiatal, dummy_Neotx3, comorb_dm, comorb_hypertension

# Name this final model, the final model: 
finalmodelPSFMI <- psfmi_lr(data=Myimp, nimp=10, impvar=".imp", Outcome="outcome", predictors=c("age", 
                                                                                           "open", 
                                                                                           "transhiatal", 
                                                                                           "dummy_Neotx3", 
                                                                                           "comorb_dm",  
                                                                                           "comorb_hypertension"))
# Perform the bootstrap on this final model: 
set.seed(21212)
res_bootstrap <- psfmi_validate(finalmodelPSFMI, int_val=TRUE, p.crit=1, nboot=500, val_method="MI_boot", cal.plot=TRUE, plot.indiv=FALSE)

# Read the results: 
res_bootstrap$stats_val
