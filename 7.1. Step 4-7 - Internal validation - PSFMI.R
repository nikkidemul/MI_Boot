### 7.1. PSFMI VALIDATION ### 

# For PSFMI, the script is as follows: 

# Name this final model, the final model: 
modelvalPSFMI <- psfmi_lr(data=Myimp, nimp=10, impvar=".imp", Outcome="outcome", predictors=c("age", 
                                                                                              "bmi",
                                                                                              "gender",
                                                                                              "open", 
                                                                                              "transhiatal", 
                                                                                              "dummy_Neotx1",
                                                                                              "dummy_Neotx3", 
                                                                                              "T34",
                                                                                              "comorb_dm",  
                                                                                              "comorb_cardiovasc",
                                                                                              "comorb_hypertension",
                                                                                              "smoking",
                                                                                              "dummy_ASA34",
                                                                                              "eGFR",
                                                                                              "dummy_Hblow",
                                                                                              "dummy_Hbhigh",
                                                                                              "fev1_compl",
                                                                                              "tiff_compl"))
# Perform the bootstrap on this final model: 
set.seed(21212)
res_bootstrap <- psfmi_validate(modelvalPSFMI, int_val=TRUE, p.crit=0.157, direction="BW", nboot=500, val_method="MI_boot", cal.plot=TRUE, plot.indiv=FALSE)

# Read the results: 
res_bootstrap$stats_val


