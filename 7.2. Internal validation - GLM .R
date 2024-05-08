### 7.2. INTERNAL VALIDATION - GLM ###


#First perform the bootstrap
set.seed(21212)
bootstrapsGLM <- lapply(Final_imp, bootstrap, bootstrapn=500)

#save results: 
saveRDS(bootstrapsGLM, file=paste0("bootstrapGLM", ".rds"))

set.seed(21212)
modelsbsimpGLM <- modelglm(bootstrapsGLM, formula=outcome~age
                           +bmi
                           +gender
                           +open
                           +transhiatal
                           +dummy_Neotx1
                           +dummy_Neotx3
                           +T34
                           +comorb_dm
                           +comorb_cardiovasc
                           +comorb_hypertension
                           +smoking
                           +dummy_ASA34
                           +eGFR
                           +dummy_Hblow
                           +dummy_Hbhigh
                           +fev1_compl
                           +tiff_compl) 

#save results for later use 
saveRDS(modelsbsimpGLM,file=paste0("modelsbsGLM",".rds"))

#develop the nullmodels over the bootstrap samples
nullmodelsbsglm <- lapply(bootstrapsGLM, nullmodels)

#calculate performance
testperformanceglm <- impPerfBSglm(imputatiesets=Final_imp, bootstrapsets=bootstrapsGLM, bsmodels=modelsbsimpGLM, bsnullmodels=nullmodelsbsglm)

#unlist to stacked dataframe, which we will need to calculate 95% confidence intervals.  
performance_rowbindtotalglm <- bind_rows(testperformanceglm) 
performance_rowbindtotalglm$performanceM <- ifelse(performance_rowbindtotalglm$performanceM == "C (ROC)", "AUC", performance_rowbindtotalglm$performanceM)

# Then pool performance measures
#take the mean of the performance measures 
mean_performanceglm <- performance_rowbindtotalglm %>% group_by(performanceM) %>% summarise(across(where(is.numeric), mean, na.rm=TRUE), .groups="drop")

#then convert back to plogis for AUC
mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="logitauc"] <- plogis(mean_performanceglm$bootstrap[mean_performanceglm$performanceM=="logitauc"])
mean_performanceglm$imp[mean_performanceglm$performanceM=="logitauc"] <- plogis(mean_performanceglm$imp[mean_performanceglm$performanceM=="logitauc"]) 

#save final results
saveRDS(mean_performanceglm,file=paste0("meanperformanceglm_val",".rds"))


#then use this function to get the standard errors, and bind the standard errors of the 10 imputation sets together. 
IPAse <- getSEs(performancelistbs=testperformanceglm, measure="IPA")
IPAse2 <- bind_rows(IPAse)
colnames(IPAse2)[1] <- "se"

#Then pool these ses using rubins rules
### IPA 
IPApoint <- performancelongGLM %>% filter(performanceM=="IPA")
SEIPA <- sqrt(rubins_rules_var(estimates=IPApoint$performanceimp, ses=IPAse2$se, n_imputed_sets=10))

#then use this value to calculate the confidence intervals 
IPAupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"] + 1.96*sqrt(SEIPA)
IPAlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"] - 1.96*sqrt(SEIPA)

#optimism corrected is original value IPA - corrected value IPA, using the same formula to calculate the confidence intervals. 
correctedIPAGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="IPA"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="IPA"]
IPAupperglmC <- correctedIPAGlm + 1.96*sqrt(SEIPA)
IPAlowerglmC <- correctedIPAGlm - 1.96*sqrt(SEIPA)

# The same for the other performance measures
### OE RATIO 
OEse <- getSEs(performancelistbs=testperformanceglm, measure="OE")
OEse2 <- bind_rows(OEse)
colnames(OEse2)[1] <- "se"

#pool ses using rubins rules
OEpoint <- performancelongGLM %>% filter(performanceM=="OE")
SEOE <- sqrt(rubins_rules_var(estimates=OEpoint$performanceimp, ses=OEse2$se, n_imputed_sets=10))

#calculate confidence intervals: 
OEupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"] + 1.96*sqrt(SEOE)
OElowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"] - 1.96*sqrt(SEOE)

#optimism corrected is original value OE - corrected value OE, using the same formula to calculate the confidence intervals. 
correctedOEGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="OE"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="OE"]
OEupperglmC <- correctedOEGlm + 1.96*sqrt(SEOE)
OElowerglmC <- correctedOEGlm - 1.96*sqrt(SEOE)

### BRIER SCORE
Brierse <- getSEs(performancelistbs=testperformanceglm, measure="Brier")
Brierse2 <- bind_rows(Brierse)
colnames(Brierse2)[1] <- "se"

Brierpoint <- performancelongGLM %>% filter(performanceM=="Brier")
SEBrier <- sqrt(rubins_rules_var(estimates=Brierpoint$performanceimp, ses=Brierse2$se, n_imputed_sets=10))

BrierupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"] + 1.96*sqrt(SEBrier)
BrierlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"] - 1.96*sqrt(SEBrier)

#optimism corrected is original value Brier - corrected value Brier, using the same formula to calculate the confidence intervals. 
correctedBrierGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Brier"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="Brier"]
BrierupperglmC <- correctedBrierGlm + 1.96*sqrt(SEBrier)
BrierlowerglmC <- correctedBrierGlm - 1.96*sqrt(SEBrier)


### Rescaled Brier score
BrierRse <- getSEs(performancelistbs=testperformanceglm, measure="BrierR")
BrierRse2 <- bind_rows(BrierRse)
colnames(BrierRse2)[1] <- "se"

BrierRpoint <- performancelongGLM %>% filter(performanceM=="BrierR")
SEBrierR <- sqrt(rubins_rules_var(estimates=BrierRpoint$performanceimp, ses=BrierRse2$se, n_imputed_sets=10))

BrierRupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"] + 1.96*sqrt(SEBrierR)
BrierRlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"] - 1.96*sqrt(SEBrierR)

#optimism corrected is original value Brier - corrected value Brier, using the same formula to calculate the confidence intervals. 
correctedBrierRGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="BrierR"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="BrierR"]
BrierRupperglmC <- correctedBrierRGlm + 1.96*sqrt(SEBrierR)
BrierRlowerglmC <- correctedBrierRGlm - 1.96*sqrt(SEBrierR)


### CALIBRATION INTERCEPT 
Ise <- getSEs(performancelistbs=testperformanceglm, measure="Intercept")
Ise2 <- bind_rows(Ise)
colnames(Ise2)[1] <- "se"

Ipoint <- performancelongGLM %>% filter(performanceM=="Intercept")
SEI <- sqrt(rubins_rules_var(estimates=Ipoint$performanceimp, ses=Ise2$se, n_imputed_sets=10))

IupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"] + 1.96*sqrt(SEI)
IlowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"] - 1.96*sqrt(SEI)

#optimism corrected is original value Intercept - corrected value Intercept, using the same formula to calculate the confidence intervals. 
correctedIGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Intercept"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="Intercept"]
IupperglmC <- correctedIGlm + 1.96*sqrt(SEI)
IlowerglmC <- correctedIGlm - 1.96*sqrt(SEI)


### CALIBRATION SLOPE 
Slopese <- getSEs(performancelistbs=testperformanceglm, measure="Slope")
Slopese2 <- bind_rows(Slopese)
colnames(Slopese2)[1] <- "se"

Slopepoint <- performancelongGLM %>% filter(performanceM=="Slope")
SESlope <- sqrt(rubins_rules_var(estimates=Slopepoint$performanceimp, ses=Slopese2$se, n_imputed_sets=10))

SlopeupperGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"] + 1.96*sqrt(SESlope)
SlopelowerGLM <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"] - 1.96*sqrt(SESlope)

#optimism corrected is original value Slope - corrected value Slope, using the same formula to calculate the confidence intervals. 
correctedSlopeGlm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="Slope"] - mean_performanceglm$optimism[mean_performanceglm$performanceM=="Slope"]
SlopeupperglmC <- correctedSlopeGlm + 1.96*sqrt(SESlope)
SlopelowerglmC <- correctedSlopeGlm - 1.96*sqrt(SESlope)



### AUC 
# restructure 
mean_performanceglm$performanceM <- ifelse(mean_performanceglm$performanceM=="auc", "AUC", mean_performanceglm$performanceM)

#corrected AUC 
correctedAUCglm <- finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"] - (mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"])

# confidence intervals 
AUClowerGLM <- plogis(logit(finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"])-1.96*sqrt(total_var_performance_impGLM))
AUCupperGLM <- plogis(logit(finalperformanceGLM$Meanperformance[finalperformanceGLM$performanceM=="AUC"])+1.96*sqrt(total_var_performance_impGLM))

#corrected
AUClowerGLMc <- AUClowerGLM - mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"]
AUCupperGLMc <- AUCupperGLM - mean_performanceglm$optimism[mean_performanceglm$performanceM=="AUC"]


### FINAL PERFORMANCE TABLE
#create final performance table 
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

saveRDS(GLMperformancetable, paste0("totalGLMperformance", ".rds"))
#nog even afronden op 2 decimalen
GLMperformancetable2 <- GLMperformancetable %>%mutate(across(is.numeric, round, digits=2))

saveRDS(GLMperformancetable2, paste0("finalGLMperformance",".rds"))

#check results: 
GLMperformancetable2





