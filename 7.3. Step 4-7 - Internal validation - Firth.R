### 7.3. INTERNAL VALIDATION - Penalized logistic regression (Firth) ### 

# In this script, we will perform the internal validation of the model obtained in script 5.3. using bootstrapping (500 bootstraps).

# First perform the bootstrap
set.seed(21212)
bootstrapsimpF <- lapply(Final_impF, bootstrap, bootstrapn=50)

# Save results: 
#saveRDS(bootstrapsimpF, file=paste0("bootstrapsimpF", ".rds"))

# Then, we repeat the model building proces from step 2 and 3 in the manuscript (script 5.2.).
### It is important that we use the same technique for the model development in the bootstrap, as in the original data (i imputed datasets). 

set.seed(21212)
modelsbsimpF <- map(.progress = TRUE, bootstrapsimpF, backwardimpF, formula=outcome~age+
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

# Then we develop the null models over the boostrap samples as we did for the apparent model: 
nullmodelsF <- lapply(bootstrapsimpF, nullmodels)

# Then we calculate the performance measures over the bootstrap samples: 
testperformanceF <- impPerfBSF(imputatiesets=Final_impF, bootstrapsets=bootstrapsimpF, bsmodels=modelsbsimpF, bsnullmodels=nullmodelsF)

#save results
#saveRDS(testperformanceF, file=paste0("listperformanceF",".rds"))

# This has created a list with 500 performance measures for 10 imputed datasets. We will stack this into one dataframe which we will need to calculate 95% CI's surrounding point estimates.  
performance_rowbindtotalF <- bind_rows(testperformanceF) 

# Rename some of the values to make it consistent: 
performance_rowbindtotalF$performanceM <- ifelse(performance_rowbindtotalF$performanceM == "C (ROC)", "AUC", performance_rowbindtotalF$performanceM)

# Then we pool the performance measures. For point estimates, we take the average. 
mean_performanceF <- performance_rowbindtotalF %>% group_by(performanceM) %>% summarise(across(where(is.numeric), mean, na.rm=TRUE), .groups="drop")

# Then, as for the apparent model performance, we have to transform back the logit AUC (plogis function). 
mean_performanceF$bootstrap[mean_performanceF$performanceM=="logitauc"] <-   plogis(mean_performanceF$bootstrap[mean_performanceF$performanceM=="logitauc"])
mean_performanceF$imp[mean_performanceF$performanceM=="logitauc"] <- plogis(mean_performanceF$imp[mean_performanceF$performanceM=="logitauc"]) 

# Save final results: 
#saveRDS(mean_performanceF,file=paste0("mean_performanceF",".rds"))



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
IPAse <- getSEs(performancelistbs=testperformanceF, measure = "IPA")
IPAse2 <- bind_rows(IPAse)
colnames(IPAse2)[1] <- "se"

### Get between imputation set standard errors: 
IPApoint <- performancelongF %>% filter(performanceM=="IPA")
SEIPA <- sqrt(rubins_rules_var(estimates=IPApoint$performanceimp, ses=IPAse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
IPAupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"] + 1.96*SEIPA
IPAlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"] - 1.96*SEIPA

### Calculate the optimism correct point estimate + 95% CI. 
correctedIPA <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"] - mean_performanceF$optimism[mean_performanceF$performanceM=="IPA"]
IPAupperFC <- correctedIPA + 1.96*SEIPA
IPAlowerFC <- correctedIPA - 1.96*SEIPA


##### O:E ratio #####
### Get within imputation set standard errors: 
OEse <- getSEs(performancelistbs=testperformanceF, measure="OE")
OEse2 <- bind_rows(OEse)
colnames(OEse2)[1] <- "se"

### Get between imputation set standard errors: 
OEpoint <- performancelongF %>% filter(performanceM=="OE")
SEOE <- sqrt(rubins_rules_var(estimates=OEpoint$performanceimp, ses=OEse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
OEupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"] + 1.96*SEOE
OElowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"] - 1.96*SEOE

### Calculate the optimism correct point estimate + 95% CI. 
correctedOE <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"] - mean_performanceF$optimism[mean_performanceF$performanceM=="OE"]
OEupperFC <- correctedOE + 1.96*SEOE
OElowerFC <- correctedOE - 1.96*SEOE


##### Brier score #####
### Get within imputation set standard errors: 
Brierse <- getSEs(performancelistbs=testperformanceF, measure="Brier")
Brierse2 <- bind_rows(Brierse)
colnames(Brierse2)[1] <- "se"

### Get between imputation set standard errors: 
Brierpoint <- performancelongF %>% filter(performanceM=="Brier")
SEBrier <- sqrt(rubins_rules_var(estimates=Brierpoint$performanceimp, ses=Brierse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
BrierupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"] + 1.96*SEBrier
BrierlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"] - 1.96*SEBrier

### Calculate the optimism correct point estimate + 95% CI. 
correctedBrier <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"] - mean_performanceF$optimism[mean_performanceF$performanceM=="Brier"]
BrierupperFC <- correctedBrier + 1.96*SEBrier
BrierlowerFC <- correctedBrier - 1.96*SEBrier


##### Rescaled Brier score #####
### Get within imputation set standard errors: 
BrierRse <- getSEs(performancelistbs=testperformanceF, measure="BrierR")
BrierRse2 <- bind_rows(BrierRse)
colnames(BrierRse2)[1] <- "se"

### Get between imputation set standard errors: 
BrierRpoint <- performancelongF %>% filter(performanceM=="BrierR")
SEBrierR <- sqrt(rubins_rules_var(estimates=BrierRpoint$performanceimp, ses=BrierRse2$se, n_imputed_sets=10))

### Calculate the confidence intervals
BrierRupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"] + 1.96*SEBrierR
BrierRlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"] - 1.96*SEBrierR

### Calculate the optimism correct point estimate + 95% CI. 
correctedBrierR <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"] - mean_performanceF$optimism[mean_performanceF$performanceM=="BrierR"]
BrierRupperFC <- correctedBrierR + 1.96*SEBrierR
BrierRlowerFC <- correctedBrierR - 1.96*SEBrierR

##### Calibration Intercept #####
### Get within imputation set standard errors: 
Ise <- getSEs(performancelistbs=testperformanceF, measure="Intercept")
Ise2 <- bind_rows(Ise)
colnames(Ise2)[1] <- "se"

### Get between imputation set standard errors: 
Ipoint <- performancelongF %>% filter(performanceM=="Intercept")
SEI <- sqrt(rubins_rules_var(estimates=Ipoint$performanceimp, ses=Ise2$se, n_imputed_sets=10))

### Calculate the confidence intervals
IupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"] + 1.96*SEI
IlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"] - 1.96*SEI

### Calculate the optimism correct point estimate + 95% CI. 
correctedI <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"] - mean_performanceF$optimism[mean_performanceF$performanceM=="Intercept"]
IupperFC <- correctedI + 1.96*SEI
IlowerFC <- correctedI - 1.96*SEI


##### Calibration Slope #####
### Get within imputation set standard errors: 
Slopese <- getSEs(performancelistbs=testperformanceF, measure="Slope")
Slopese2 <- bind_rows(Slopese)
colnames(Slopese2)[1] <- "se"

### Get between imputation set standard errors: 
Slopepoint <- performancelongF %>% filter(performanceM=="Slope")
SESlope <- sqrt(rubins_rules_var(estimates=Slopepoint$performanceimp, ses=Slopese2$se, n_imputed_sets=10))

### Calculate the confidence intervals
SlopeupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"] + 1.96*SESlope
SlopelowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"] - 1.96*SESlope

### Calculate the optimism correct point estimate + 95% CI. 
correctedSlope <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"] - mean_performanceF$optimism[mean_performanceF$performanceM=="Slope"]
SlopeupperFC <- correctedSlope + 1.96*SESlope
SlopelowerFC <- correctedSlope - 1.96*SESlope


##### AUC #####
# For the AUC, we need the variance we saved earlier in 6.3. 
# First, we rename and restructure our data: 
mean_performanceF$performanceM <- ifelse(mean_performanceF$performanceM=="auc", "AUC", mean_performanceF$performanceM)

# We then calculate the corrected AUC using the optimism obtained from the bootstrap procedure. 
correctedAUC <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"] - (mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"])

# We then calculate the confidence intervals using the variance calculated in script 6.3. 
AUClowerF <- plogis(logit(finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"])-1.96*sqrt(total_var_performance_impF))
AUCupperF <- plogis(logit(finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"])+1.96*sqrt(total_var_performance_impF))

# Then the confidence intervals surrounding the corrected performance measure (corrected for optimism)
AUClowerFC <- AUClowerF - mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"]
AUCupperFC <- AUCupperF - mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"]


### FINAL PERFORMANCE TABLE ###
# Then we create the final performance table as output, combining all the above calculated measures. 

performanceM <- c("AUC", "OE", "Intercept", "Slope", "Brier", "BrierR", "IPA")
original <- c(finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"], 
              finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"], 
              finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"], 
              finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"], 
              finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"], 
              finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"], 
              finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"])
lowerO <- c(AUClowerF, 
            OElowerF, 
            IlowerF, 
            SlopelowerF, 
            BrierlowerF, 
            BrierRlowerF, 
            IPAlowerF)
upperO <- c(AUCupperF, 
            OEupperF, 
            IupperF, 
            SlopeupperF, 
            BrierupperF, 
            BrierRupperF, 
            IPAupperF)
bootstrap <- c(mean_performanceF$bootstrap[mean_performanceF$performanceM=="logitauc"], 
               mean_performanceF$bootstrap[mean_performanceF$performanceM=="OE"], 
               mean_performanceF$bootstrap[mean_performanceF$performanceM=="Intercept"], 
               mean_performanceF$bootstrap[mean_performanceF$performanceM=="Slope"], 
               mean_performanceF$bootstrap[mean_performanceF$performanceM=="Brier"], 
               mean_performanceF$bootstrap[mean_performanceF$performanceM=="BrierR"],
               mean_performanceF$bootstrap[mean_performanceF$performanceM=="IPA"])
imp <- c(mean_performanceF$imp[mean_performanceF=="logitauc"], 
         mean_performanceF$imp[mean_performanceF=="OE"], 
         mean_performanceF$imp[mean_performanceF=="Intercept"], 
         mean_performanceF$imp[mean_performanceF=="Slope"], 
         mean_performanceF$imp[mean_performanceF=="Brier"], 
         mean_performanceF$imp[mean_performanceF=="BrierR"], 
         mean_performanceF$imp[mean_performanceF=="IPA"])
optimism <- c(mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"], 
              mean_performanceF$optimism[mean_performanceF$performanceM=="OE"], 
              mean_performanceF$optimism[mean_performanceF$performanceM=="Intercept"], 
              mean_performanceF$optimism[mean_performanceF$performanceM=="Slope"], 
              mean_performanceF$optimism[mean_performanceF$performanceM=="Brier"], 
              mean_performanceF$optimism[mean_performanceF$performanceM=="BrierR"], 
              mean_performanceF$optimism[mean_performanceF$performanceM=="IPA"])
corrected <- c(correctedAUC, 
               correctedOE, 
               correctedI, 
               correctedSlope, 
               correctedBrier, 
               correctedBrierR, 
               correctedIPA)
lowerC <- c(AUClowerFC, 
            OElowerFC, 
            IlowerFC, 
            SlopelowerFC, 
            BrierlowerFC, 
            BrierRlowerFC, 
            IPAlowerFC)
upperC <- c(AUCupperFC, 
            OEupperFC, 
            IupperFC, 
            SlopeupperFC, 
            BrierupperFC, 
            BrierRupperFC, 
            IPAupperFC)


Fperformancetable <- as.data.frame(cbind(performanceM, original, lowerO, upperO, bootstrap, imp, optimism, corrected, lowerC, upperC))
Fperformancetable <- Fperformancetable %>% mutate(original = as.numeric(original)) %>% mutate(lowerO = as.numeric(lowerO)) %>% mutate(upperO = as.numeric(upperO)) %>% mutate(bootstrap = as.numeric(bootstrap)) %>% mutate(imp = as.numeric(imp)) %>% mutate(optimism = as.numeric(optimism)) %>% mutate(corrected=as.numeric(corrected)) %>% mutate(lowerC = as.numeric(lowerC)) %>% mutate(upperC = as.numeric(upperC))

#saveRDS(Fperformancetable, paste0("totalperformanceFirth_bs", ".rds"))

# Then we round to two decimals to enhance readibility: 
Fperformancetable2 <- Fperformancetable %>%mutate(across(is.numeric, round, digits=2))

#saveRDS(Fperformancetable2, paste0("finalperformanceFirth_bs",".rds"))

# Check results: 
Fperformancetable2

# Interpretation: 
# As discussed in the manuscript: 
### In the left column we list the performance measures calculated. 
### Then 'original' shows the original performance measure (apparent performance) of the apparent model, on the original data. 
### The lowerO and upperO represent the confidence intervals (lower and upper bounds) of this apparent performance measure. 
### The bootstrap shows the performance  of the bootstrap model on the bootstrap sample the model was developed in. 
### The imp shows the performance of the bootstrap model on the imputation set their bootstrap sample was derived from (bootstrap on imputation set). 
### Optimism = imp - bootstrap. 
### Corrected performance measure is the original performance - optimism. 
### The lowerC and upperC represent the confidence intervals (lower and upper boundds) of this corrected performance measure.



