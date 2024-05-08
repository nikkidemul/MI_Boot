### 7.3. INTERNAL VALIDATION - FIRTH ### 

#first perform the bootstrap
set.seed(21212)
bootstrapsimpF <- lapply(Final_impF, bootstrap, bootstrapn=500)

#save results: 
saveRDS(bootstrapsimpF, file=paste0("bootstrapsimpF", ".rds"))

#then redevelop the models

modelsbsimpF <- lapply(bootstrapsimpF, backwardimpF, formula=outcome~age+
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

# redevelop nullmodels: 
nullmodelsF <- lapply(bootstrapsimpF, nullmodels)

#calculate performance
testperformanceF <- impPerfBSF(imputatiesets=Final_impF, bootstrapsets=bootstrapsimpF, bsmodels=modelsbsimpF, bsnullmodels=nullmodelsF)

#save results
saveRDS(testperformanceF, file=paste0("listperformanceF",".rds"))

#unlist to stacked dataframe, which we will need to calculate 95% confidence intervals.  
performance_rowbindtotalF <- bind_rows(testperformanceF) 
performance_rowbindtotalF$performanceM <- ifelse(performance_rowbindtotalF$performanceM == "C (ROC)", "AUC", performance_rowbindtotalF$performanceM)

# pool performance measures

#take the mean of the performance measures 
mean_performanceF <- performance_rowbindtotalF %>% group_by(performanceM) %>% summarise(across(where(is.numeric), mean, na.rm=TRUE), .groups="drop")

#then convert back to plogis for AUC
mean_performanceF$bootstrap[mean_performanceF$performanceM=="logitauc"] <-   plogis(mean_performanceF$bootstrap[mean_performanceF$performanceM=="logitauc"])
mean_performanceF$imp[mean_performanceF$performanceM=="logitauc"] <- plogis(mean_performanceF$imp[mean_performanceF$performanceM=="logitauc"]) 

saveRDS(mean_performanceF,file=paste0("mean_performanceF",".rds"))



# Calculating confidence intervals 
#### IPA 
IPAse <- getSEs(performancelistbs=testperformanceF, measure = "IPA")
IPAse2 <- bind_rows(IPAse)
colnames(IPAse2)[1] <- "se"

#pool using Rubin's rules 
IPApoint <- performancelongF %>% filter(performanceM=="IPA")
SEIPA <- sqrt(rubins_rules_var(estimates=IPApoint$performanceimp, ses=IPAse2$se, n_imputed_sets=10))

#then use this value to calculate the confidence intervals 
IPAupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"] + 1.96*sqrt(SEIPA)
IPAlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"] - 1.96*sqrt(SEIPA)

#optimism corrected is original value IPA - corrected value IPA, using the same formula to calculate the confidence intervals. 
correctedIPA <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="IPA"] - mean_performanceF$optimism[mean_performanceF$performanceM=="IPA"]
IPAupperFC <- correctedIPA + 1.96*sqrt(SEIPA)
IPAlowerFC <- correctedIPA - 1.96*sqrt(SEIPA)


#### OE 
OEse <- getSEs(performancelistbs=testperformanceF, measure="OE")
OEse2 <- bind_rows(OEse)
colnames(OEse2)[1] <- "se"

#pool ses using rubins rules
OEpoint <- performancelongF %>% filter(performanceM=="OE")
SEOE <- sqrt(rubins_rules_var(estimates=OEpoint$performanceimp, ses=OEse2$se, n_imputed_sets=10))

#calculate confidence intervals: 
OEupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"] + 1.96*sqrt(SEOE)
OElowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"] - 1.96*sqrt(SEOE)

#corrected for optimism (same method as described for IPA): 
correctedOE <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="OE"] - mean_performanceF$optimism[mean_performanceF$performanceM=="OE"]
OEupperFC <- correctedOE + 1.96*sqrt(SEOE)
OElowerFC <- correctedOE - 1.96*sqrt(SEOE)


### Brier
Brierse <- getSEs(performancelistbs=testperformanceF, measure="Brier")
Brierse2 <- bind_rows(Brierse)
colnames(Brierse2)[1] <- "se"

Brierpoint <- performancelongF %>% filter(performanceM=="Brier")
SEBrier <- sqrt(rubins_rules_var(estimates=Brierpoint$performanceimp, ses=Brierse2$se, n_imputed_sets=10))

BrierupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"] + 1.96*sqrt(SEBrier)
BrierlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"] - 1.96*sqrt(SEBrier)

#corrected: 
correctedBrier <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Brier"] - mean_performanceF$optimism[mean_performanceF$performanceM=="Brier"]
BrierupperFC <- correctedBrier + 1.96*sqrt(SEBrier)
BrierlowerFC <- correctedBrier - 1.96*sqrt(SEBrier)


### Rescaled Brier
### Brier
BrierRse <- getSEs(performancelistbs=testperformanceF, measure="BrierR")
BrierRse2 <- bind_rows(BrierRse)
colnames(BrierRse2)[1] <- "se"

BrierRpoint <- performancelongF %>% filter(performanceM=="BrierR")
SEBrierR <- sqrt(rubins_rules_var(estimates=BrierRpoint$performanceimp, ses=BrierRse2$se, n_imputed_sets=10))

BrierRupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"] + 1.96*sqrt(SEBrierR)
BrierRlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"] - 1.96*sqrt(SEBrierR)

#corrected: 
correctedBrierR <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="BrierR"] - mean_performanceF$optimism[mean_performanceF$performanceM=="BrierR"]
BrierRupperFC <- correctedBrierR + 1.96*sqrt(SEBrierR)
BrierRlowerFC <- correctedBrierR - 1.96*sqrt(SEBrierR)

### Calibration intercept 
Ise <- getSEs(performancelistbs=testperformanceF, measure="Intercept")
Ise2 <- bind_rows(Ise)
colnames(Ise2)[1] <- "se"

Ipoint <- performancelongF %>% filter(performanceM=="Intercept")
SEI <- sqrt(rubins_rules_var(estimates=Ipoint$performanceimp, ses=Ise2$se, n_imputed_sets=10))
IupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"] + 1.96*sqrt(SEI)
IlowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"] - 1.96*sqrt(SEI)

#corrected: 
correctedI <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Intercept"] - mean_performanceF$optimism[mean_performanceF$performanceM=="Intercept"]
IupperFC <- correctedI + 1.96*sqrt(SEI)
IlowerFC <- correctedI - 1.96*sqrt(SEI)


### Calibration Slope
Slopese <- getSEs(performancelistbs=testperformanceF, measure="Slope")
Slopese2 <- bind_rows(Slopese)
colnames(Slopese2)[1] <- "se"

Slopepoint <- performancelongF %>% filter(performanceM=="Slope")
SESlope <- sqrt(rubins_rules_var(estimates=Slopepoint$performanceimp, ses=Slopese2$se, n_imputed_sets=10))
SlopeupperF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"] + 1.96*sqrt(SESlope)
SlopelowerF <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"] - 1.96*sqrt(SESlope)

#corrected: 
correctedSlope <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="Slope"] - mean_performanceF$optimism[mean_performanceF$performanceM=="Slope"]
SlopeupperFC <- correctedSlope + 1.96*sqrt(SESlope)
SlopelowerFC <- correctedSlope - 1.96*sqrt(SESlope)


### AUC
# restructure data
mean_performanceF$performanceM <- ifelse(mean_performanceF$performanceM=="auc", "AUC", mean_performanceF$performanceM)

# corrected AUC
correctedAUC <- finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"] - (mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"])

# confidence interval surrounding apparent AUC 
AUClowerF <- plogis(logit(finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"])-1.96*sqrt(total_var_performance_impF))
AUCupperF <- plogis(logit(finalperformanceF$Meanperformance[finalperformanceF$performanceM=="AUC"])+1.96*sqrt(total_var_performance_impF))

# confidence interval surrounding corrected AUC 
AUClowerFC <- AUClowerF - mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"]
AUCupperFC <- AUCupperF - mean_performanceF$optimism[mean_performanceF$performanceM=="AUC"]


#### FINAL PERFORMANCE TABLE ###
# Putting it all together: 
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

saveRDS(Fperformancetable, paste0("totalperformanceFirth_bs", ".rds"))
#round to 2 decimals: 
Fperformancetable2 <- Fperformancetable %>%mutate(across(is.numeric, round, digits=2))

saveRDS(Fperformancetable2, paste0("finalperformanceFirth_bs",".rds"))

Fperformancetable2
