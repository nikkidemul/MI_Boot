### 7.1. PSFMI VALIDATION ### 

#The final model seleted by PSFMI contained the following variables
modelD1_back$RR_model_final
# age, open, transhiatal, dummy_Neotx3, comorb_dm, comorb_hypertension

# perform the internal validation using bootstrapping
finalmodelPSFMI <- psfmi_lr(data=Myimp, nimp=10, impvar=".imp", Outcome="outcome", predictors=c("age", 
                                                                                           "open", 
                                                                                           "transhiatal", 
                                                                                           "dummy_Neotx3", 
                                                                                           "comorb_dm",  
                                                                                           "comorb_hypertension"))

set.seed(21212)
res_bootstrap <- psfmi_validate(finalmodelPSFMI, int_val=TRUE, p.crit=1, nboot=500, val_method="MI_boot", cal.plot=TRUE, plot.indiv=FALSE)

res_bootstrap$stats_val
