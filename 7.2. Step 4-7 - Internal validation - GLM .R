### 7.2. INTERNAL VALIDATION - Non-penalized logistic regression (GLM) ###

# In this script, we will perform the internal validation of the model obtained in script 5.2. using bootstrapping (500 bootstraps). 
### Check out 2. and 2.2. Function libraries for the code on the specific functions used. 


#First perform the bootstrap
set.seed(21212)
bootstrapsGLM <- lapply(Final_imp, bootstrap, bootstrapn=500)

#save results: 
#saveRDS(bootstrapsGLM, file=paste0("bootstrapGLM", ".rds"))

# Then, we repeat the model building process from step 2 and step 3 (script 5.2). 
### It is important that we use the same technique for the model development in the bootstrap, as in the original data (i imputed datasets). 
set.seed(21212)
tic <- Sys.time()
modelsbsimpGLM <- modelglm_pl(nr_cores = 2, 
                              bootstrapsets = bootstrapsGLM, 
                              formula=outcome~age+
                                              bmi+
                                              gender+
                                              open+
                                              transhiatal+
                                              dummy_Neotx1+
                                              dummy_Neotx3+
                                              T34+
                                              comorb_dm+
                                              comorb_cardiovasc+
                                              comorb_hypertension+
                                              smoking+
                                              dummy_ASA34+
                                              eGFR+
                                              dummy_Hblow+
                                              dummy_Hbhigh+
                                              fev1_compl+
                                              tiff_compl)
toc <- Sys.time()
print(toc - tic)
print(Sys.time())
#save results for later use 
#saveRDS(modelsbsimpGLM,file=paste0("modelsbsGLM",".rds"))

# Then we develop the nullmodels over the bootstrap samples as we did for the apparent model: 
nullmodelsbsglm <- lapply(bootstrapsGLM, nullmodels)

# Then we calculate the performance measures over the bootstrap samples: 
testperformanceglm <- impPerfBSglm(imputatiesets=Final_imp, bootstrapsets=bootstrapsGLM, bsmodels=modelsbsimpGLM, bsnullmodels=nullmodelsbsglm)

# This has created a list with 500 performance measures for 10 imputed datasets. We will stack this into one dataframe which we will need to calculate the 95% CI's surrounding the point estimates. 
performance_rowbindtotalglm <- bind_rows(testperformanceglm) 

# Rename some of the values to make it consistent: 
performance_rowbindtotalglm$performanceM <- ifelse(performance_rowbindtotalglm$performanceM == "C (ROC)", "AUC", performance_rowbindtotalglm$performanceM)

# Then we pool the performance measures. For point estimates, we take the average. 
mean_performanceglm <- performance_rowbindtotalglm %>% group_by(performanceM) %>% summarise(across(where(is.numeric), \(e) mean(e, na.rm = TRUE)), .groups="drop")

# Then, as for the apparent model performance, we have to transform back the logit AUC  (plogis function). 
mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="logitauc"] <- plogis(mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="logitauc"])
mean_performanceglm$imp[mean_performanceglm$performanceM=="logitauc"] <- plogis(mean_performanceglm$imp[mean_performanceglm$performanceM=="logitauc"]) 

#save final results
#saveRDS(mean_performanceglm,file=paste0("meanperformanceglm_val",".rds"))

# Then we will retrieve the standard errors using a function called 'getSEs' for which we refer to the function libraries. 
# This function provides us with the standard errors within each bootstrap procedure (for each imputationset i). 
### Calculation: we take the 500 standard errors for each imputation set, and take the 0.025 and 0.975 percentile. We divide this by 3.92 (2*1.96) to calculate the standard error within the imputation set. 
### This is pragmatically chosen because upper limit confidence interval = point estimate + 1.96*standard error, and lower limit = point estimate - 1.96*standard error. 

# After calculating the within imputation set standard errors through this function, we pool the within imputation set standard errors into the final standard errors. 
### Here, we have to take between imputation set variation into account as wel, and we use Rubin's Rules to do this. 
# We use the final standard error to calculate the 95% CI surrounding our point estimate of performance. 

# For the optimism corrected values: 
### The point estimate = apparent point estimate - optimism 
### The 95% CI can be calculated using the same formula's, but using the optimism corrected point estimate as a starting point. 

##### IPA #####
### Get within imputation set standard errors: 
IPAse <- getSEs(performancelistbs=testperformanceglm, measure="IPA")
IPAse2 <- bind_rows(IPAse)
colnames(IPAse2)[1] <- "se"

### Get between imputation set standard errors: 
IPApoint <- performancelongGLM %>% filter(performanceM=="IPA")
SEIPA <- sqrt(rubins_rules_var(estimates=IPApoint$performanceimp, ses=IPAse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
IPAupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"] + 1.96*SEIPA
IPAlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"] - 1.96*SEIPA

### Calculate the optimism correct point estimate + 95% CI. 
correctedIPAGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="IPA"]
IPAupperglmC <- correctedIPAGlm + 1.96*SEIPA
IPAlowerglmC <- correctedIPAGlm - 1.96*SEIPA

# This procedure is then repeated for the other performance measures

##### O:E ratio #####
### Get within imputation set standard errors: 
OEse <- getSEs(performancelistbs=testperformanceglm, measure="OE")
OEse2 <- bind_rows(OEse)
colnames(OEse2)[1] <- "se"

### Get between imputation set standard errors: 
OEpoint <- performancelongGLM %>% filter(performanceM=="OE")
SEOE <- sqrt(rubins_rules_var(estimates=OEpoint$performanceimp, ses=OEse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
OEupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"] + 1.96*SEOE
OElowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"] - 1.96*SEOE

### Calculate the optimism correct point estimate + 95% CI. 
correctedOEGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="OE"]
OEupperglmC <- correctedOEGlm + 1.96*SEOE
OElowerglmC <- correctedOEGlm - 1.96*SEOE

##### Brier score #####
### Get within imputation set standard errors: 
Brierse <- getSEs(performancelistbs=testperformanceglm, measure="Brier")
Brierse2 <- bind_rows(Brierse)
colnames(Brierse2)[1] <- "se"

### Get between imputation set standard errors: 
Brierpoint <- performancelongGLM %>% filter(performanceM=="Brier")
SEBrier <- sqrt(rubins_rules_var(estimates=Brierpoint$performanceimp, ses=Brierse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
BrierupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"] + 1.96*SEBrier
BrierlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"] - 1.96*SEBrier

### Calculate the optimism correct point estimate + 95% CI. 
correctedBrierGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="Brier"]
BrierupperglmC <- correctedBrierGlm + 1.96*SEBrier
BrierlowerglmC <- correctedBrierGlm - 1.96*SEBrier


##### Rescaled Brier Score #####
### Get within imputation set standard errors: 
BrierRse <- getSEs(performancelistbs=testperformanceglm, measure="BrierR")
BrierRse2 <- bind_rows(BrierRse)
colnames(BrierRse2)[1] <- "se"

### Get between imputation set standard errors: 
BrierRpoint <- performancelongGLM %>% filter(performanceM=="BrierR")
SEBrierR <- sqrt(rubins_rules_var(estimates=BrierRpoint$performanceimp, ses=BrierRse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
BrierRupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"] + 1.96*SEBrierR
BrierRlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"] - 1.96*SEBrierR

### Calculate the optimism correct point estimate + 95% CI. 
correctedBrierRGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="BrierR"]
BrierRupperglmC <- correctedBrierRGlm + 1.96*SEBrierR
BrierRlowerglmC <- correctedBrierRGlm - 1.96*SEBrierR


##### Calibration intercept #####
### Get within imputation set standard errors: 
Ise <- getSEs(performancelistbs=testperformanceglm, measure="Intercept")
Ise2 <- bind_rows(Ise)
colnames(Ise2)[1] <- "se"

### Get between imputation set standard errors: 
Ipoint <- performancelongGLM %>% filter(performanceM=="Intercept")
SEI <- sqrt(rubins_rules_var(estimates=Ipoint$performanceimp, ses=Ise2$se, n_imputed_sets=10))

### Calculate the confidence intervals
IupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"] + 1.96*SEI
IlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"] - 1.96*SEI

### Calculate the optimism correct point estimate + 95% CI. 
correctedIGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="Intercept"]
IupperglmC <- correctedIGlm + 1.96*SEI
IlowerglmC <- correctedIGlm - 1.96*SEI


##### Calibration slope #####
### Get within imputation set standard errors: 
Slopese <- getSEs(performancelistbs=testperformanceglm, measure="Slope")
Slopese2 <- bind_rows(Slopese)
colnames(Slopese2)[1] <- "se"

### Get between imputation set standard errors: 
Slopepoint <- performancelongGLM %>% filter(performanceM=="Slope")
SESlope <- sqrt(rubins_rules_var(estimates=Slopepoint$performanceimp, ses=Slopese2$se, n_imputed_sets=10))

### Calculate the confidence intervals
SlopeupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"] + 1.96*SESlope
SlopelowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"] - 1.96*SESlope

### Calculate the optimism correct point estimate + 95% CI. 
correctedSlopeGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="Slope"]
SlopeupperglmC <- correctedSlopeGlm + 1.96*SESlope
SlopelowerglmC <- correctedSlopeGlm - 1.96*SESlope



##### AUC #####
# For the AUC, we need the variance we save earlier in 6.2. 
# First, we rename and restructure our data. 
mean_performanceglm$performanceM <- ifelse(mean_performanceglm$performanceM=="auc", "AUC", mean_performanceglm$performanceM)

# We then calculate the corrected AUC using the optimism obtained from the bootstrap procedure.  
correctedAUCglm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"] - (mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"])

# We then calculate the confidence intervals using the variance calculated in script 6.2. 
AUClowerGLM <- plogis(logit(finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"])-1.96*sqrt(total_var_performance_impGLM))
AUCupperGLM <- plogis(logit(finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"])+1.96*sqrt(total_var_performance_impGLM))

# Then corrected confidence intervals: 
AUClowerGLMc <- AUClowerGLM - mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"]
AUCupperGLMc <- AUCupperGLM - mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"]


### FINAL PERFORMANCE TABLE ###
# Then we create the final performance table as output, combining all the above calculated measures. 

performanceM <- c("AUC", "OE", "Intercept", "Slope", "Brier", "BrierR", "IPA")
original <- c(finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"], 
              finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"], 
              finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"], 
              finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"], 
              finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"], 
              finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"], 
              finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"])
lowerO <- c(AUClowerGLM, 
            OElowerGLM, 
            IlowerGLM, 
            SlopelowerGLM, 
            BrierlowerGLM, 
            BrierRlowerGLM,
            IPAlowerGLM)
upperO <- c(AUCupperGLM, 
            OEupperGLM, 
            IupperGLM, 
            SlopeupperGLM, 
            BrierupperGLM, 
            BrierRupperGLM,
            IPAupperGLM)
bootstrap <- c(mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="logitauc"], 
               mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="OE"], 
               mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="Intercept"], 
               mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="Slope"], 
               mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="Brier"], 
               mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="BrierR"], 
               mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="IPA"])
imp <- c(mean_performanceglm$imp[mean_performanceglm=="logitauc"], 
         mean_performanceglm$imp[mean_performanceglm=="OE"], 
         mean_performanceglm$imp[mean_performanceglm=="Intercept"], 
         mean_performanceglm$imp[mean_performanceglm=="Slope"], 
         mean_performanceglm$imp[mean_performanceglm=="Brier"], 
         mean_performanceglm$imp[mean_performanceglm=="BrierR"], 
         mean_performanceglm$imp[mean_performanceglm=="IPA"])
optimism <- c(mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"], 
              mean_performanceglm$optimism[mean_performanceglm$performanceM=="OE"], 
              mean_performanceglm$optimism[mean_performanceglm$performanceM=="Intercept"], 
              mean_performanceglm$optimism[mean_performanceglm$performanceM=="Slope"], 
              mean_performanceglm$optimism[mean_performanceglm$performanceM=="Brier"], 
              mean_performanceglm$optimism[mean_performanceglm$performanceM=="BrierR"], 
              mean_performanceglm$optimism[mean_performanceglm$performanceM=="IPA"])
corrected <- c(correctedAUCglm, 
               correctedOEGlm, 
               correctedIGlm, 
               correctedSlopeGlm, 
               correctedBrierGlm, 
               correctedBrierRGlm,
               correctedIPAGlm)
lowerC <- c(AUClowerGLMc, 
            OElowerglmC, 
            IlowerglmC, 
            SlopelowerglmC, 
            BrierlowerglmC, 
            BrierRlowerglmC,
            IPAlowerglmC)
upperC <- c(AUCupperGLMc, 
            OEupperglmC, 
            IupperglmC, 
            SlopeupperglmC, 
            BrierupperglmC, 
            BrierRupperglmC,
            IPAupperglmC)


GLMperformancetable <- as.data.frame(cbind(performanceM, original, lowerO, upperO, bootstrap, imp, optimism, corrected, lowerC, upperC))
GLMperformancetable <- GLMperformancetable %>% mutate(original = as.numeric(original)) %>% mutate(lowerO = as.numeric(lowerO)) %>% mutate(upperO = as.numeric(upperO)) %>% mutate(bootstrap = as.numeric(bootstrap)) %>% mutate(imp = as.numeric(imp)) %>% mutate(optimism = as.numeric(optimism)) %>% mutate(corrected=as.numeric(corrected)) %>% mutate(lowerC = as.numeric(lowerC)) %>% mutate(upperC = as.numeric(upperC))

#saveRDS(GLMperformancetable, paste0("totalGLMperformance", ".rds"))

# Then we round to two decimals to enhance readibility: 
GLMperformancetable2 <- GLMperformancetable %>%mutate(across(is.numeric, round, digits=2))

#saveRDS(GLMperformancetable2, paste0("finalGLMperformance",".rds"))

# Check results: 
GLMperformancetable2

# As discussed in the manuscript: 
### In the left column we list the performance measures calculated. 
### Then 'original' shows the original performance measure (apparent performance) of the apparent model, on the original data. 
### The lowerO and upperO represent the confidence intervals (lower and upper bounds) of this apparent performance measure. 
### The bootstrap shows the performance  of the bootstrap model on the bootstrap sample the model was developed in. 
### The imp shows the performance of the bootstrap model on the imputation set their bootstrap sample was derived from (bootstrap on imputation set). 
### Optimism = imp - bootstrap. 
### Corrected performance measure is the original performance - optimism. 
### The lowerC and upperC represent the confidence intervals (lower and upper boundds) of this corrected performance measure. 





