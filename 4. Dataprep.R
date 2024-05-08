### 4. DATAPREP ###

# In this script we will prep the data (define dummy variables, put it in long format as well for the PSFMI analysis). 

### PSFMI 
Myimp <- as.data.frame(complete(imputedset, action="long"))

#dummy_Neotx
table(Myimp$Neotx)
Myimp <- Myimp %>% mutate(dummy_Neotx1 = ifelse(Neotx==1,1,0)) %>% mutate(dummy_Neotx3 = ifelse(Neotx==3,1,0)) %>% mutate(dummy_Neotx1 = as.factor(dummy_Neotx1)) %>% mutate(dummy_Neotx3 = as.factor(dummy_Neotx3))

#eGFRcat
table(Myimp$eGFRcat)
Myimp <- Myimp %>% mutate(eGFR = ifelse(eGFRcat==2 | eGFRcat==3,1,0))
Myimp <- Myimp %>% mutate(eGFR = as.factor(eGFR))
table(Myimp$eGFR)

#categorical ASAclass
table(Myimp$ASA)
Myimp <- Myimp %>% mutate(dummy_ASA12 = ifelse(ASA==1 | ASA==2,1,0)) %>% mutate(dummy_ASA34 = ifelse(ASA==3 | ASA==4,1,0)) %>% mutate(dummy_ASA12 = as.factor(dummy_ASA12)) %>% mutate(dummy_ASA34 = as.factor(dummy_ASA34))

#Cardiovascular variable
Myimp <- Myimp %>% mutate(comorb_per_vasc = comorb_per_vasc - 1) %>% mutate(comorb_per_vasc = as.factor(comorb_per_vasc))
Myimp <- Myimp %>% mutate(comorb_cardiovasc = ifelse((comorb_mi==1 | comorb_coronary==1 | comorb_heart_failure==1 | comorb_per_vasc==1 | comorb_cerebrovasc==1),1,0)) %>% mutate(comorb_cardiovasc = as.factor(comorb_cardiovasc)) 
table(Myimp$comorb_cardiovasc)

#dummy Hb
Myimp <- Myimp %>% mutate(dummy_Hblow = ifelse(Hbcat==1,1,0)) %>% mutate(dummy_Hbhigh = ifelse(Hbcat==2,1,0)) %>% mutate(dummy_Hblow = as.factor(dummy_Hblow)) %>% mutate(dummy_Hbhigh = as.factor(dummy_Hbhigh))



### MANUAL 
### Then repeat this to construct new imputed data in list format, which we need for the manual analysis 
Myimp2 <- as.data.frame(complete(imputedset, action="long", include=TRUE))

#dummy_Neotx
Myimp2 <- Myimp2 %>% mutate(dummy_Neotx1 = ifelse(Neotx==1,1,0)) %>% mutate(dummy_Neotx3 = ifelse(Neotx==3,1,0)) %>% mutate(dummy_Neotx1 = as.factor(dummy_Neotx1)) %>% mutate(dummy_Neotx3 = as.factor(dummy_Neotx3))
table(Myimp2$dummy_Neotx3)

#eGFRcat
table(Myimp2$eGFRcat)
Myimp2 <- Myimp2 %>% mutate(eGFR = eGFRcat) %>% mutate(eGFR = ifelse(eGFR==2 | eGFR==3,1,0))
Myimp2 <- Myimp2 %>% mutate(eGFR = as.factor(eGFR))
table(Myimp2$eGFR)

#categorical ASAclass
Myimp2 <- Myimp2 %>% mutate(dummy_ASA12 = ifelse(ASA==1 | ASA==2,1,0)) %>% mutate(dummy_ASA34 = ifelse(ASA==3 | ASA==4,1,0)) %>% mutate(dummy_ASA12 = as.factor(dummy_ASA12)) %>% mutate(dummy_ASA34 = as.factor(dummy_ASA34))

#Cardiovascular variable
Myimp2 <- Myimp2 %>% mutate(comorb_per_vasc = comorb_per_vasc - 1) %>% mutate(comorb_per_vasc = as.factor(comorb_per_vasc))
Myimp2 <- Myimp2 %>% mutate(comorb_cardiovasc = ifelse((comorb_mi==1 | comorb_coronary==1 | comorb_heart_failure==1 | comorb_per_vasc==1 | comorb_cerebrovasc==1),1,0)) %>% mutate(comorb_cardiovasc = as.factor(comorb_cardiovasc)) 

#dummy Hb
Myimp2 <- Myimp2 %>% mutate(dummy_Hblow = ifelse(Hbcat==1,1,0)) %>% mutate(dummy_Hbhigh = ifelse(Hbcat==2,1,0)) %>% mutate(dummy_Hblow = as.factor(dummy_Hblow)) %>% mutate(dummy_Hbhigh = as.factor(dummy_Hbhigh))

table(Myimp2$outcome)

Final_imp <- as.mids(Myimp2, where=NULL, .imp=".imp", .id=".id")

Final_imp <- complete(Final_imp, "all")

saveRDS(Final_imp, file=paste0("Final_imp", ".rds"))
saveRDS(Myimp, file=paste0("Myimp", ".rds"))



### MANUAL FIRTH #### 
### For Firth, we need all the data in numeric format (it does not handle factor variables). Therefore we repeat it one more time and rename these files as Firth.
Myimp2 <- as.data.frame(complete(imputedset, action="long", include=TRUE))

#dummy_Neotx
Myimp2 <- Myimp2 %>% mutate(dummy_Neotx1 = ifelse(Neotx==1,1,0)) %>% mutate(dummy_Neotx3 = ifelse(Neotx==3,1,0)) %>% mutate(dummy_Neotx1 = as.factor(dummy_Neotx1)) %>% mutate(dummy_Neotx3 = as.factor(dummy_Neotx3))

#eGFRcat
Myimp2 <- Myimp2 %>% mutate(eGFR = ifelse(eGFRcat==2 | eGFRcat==3,1,0))
Myimp2 <- Myimp2 %>% mutate(eGFR = as.factor(eGFR))
table(Myimp2$eGFR)
table(Myimp2$eGFRcat)

#categorical ASAclass
Myimp2 <- Myimp2 %>% mutate(dummy_ASA12 = ifelse(ASA==1 | ASA==2,1,0)) %>% mutate(dummy_ASA34 = ifelse(ASA==3 | ASA==4,1,0)) %>% mutate(dummy_ASA12 = as.factor(dummy_ASA12)) %>% mutate(dummy_ASA34 = as.factor(dummy_ASA34))

#Cardiovascular variable
Myimp2 <- Myimp2 %>% mutate(comorb_per_vasc = comorb_per_vasc - 1) %>% mutate(comorb_per_vasc = as.factor(comorb_per_vasc))
Myimp2 <- Myimp2 %>% mutate(comorb_cardiovasc = ifelse((comorb_mi==1 | comorb_coronary==1 | comorb_heart_failure==1 | comorb_per_vasc==1 | comorb_cerebrovasc==1),1,0)) %>% mutate(comorb_cardiovasc = as.factor(comorb_cardiovasc)) 
table(Myimp2$comorb_cardiovasc)

#dummy Hb
Myimp2 <- Myimp2 %>% mutate(dummy_Hblow = ifelse(Hbcat==1,1,0)) %>% mutate(dummy_Hbhigh = ifelse(Hbcat==2,1,0)) %>% mutate(dummy_Hblow = as.factor(dummy_Hblow)) %>% mutate(dummy_Hbhigh = as.factor(dummy_Hbhigh))

Myimp3 <- Myimp2
x <- sapply(Myimp3, is.factor)
Myimp3[ ,x] <- as.data.frame(apply(Myimp3[, x], 2, as.numeric))

Final_impF <- as.mids(Myimp3, where=NULL, .imp=".imp", .id=".id")

Final_impF <- complete(Final_impF, "all")

saveRDS(Final_impF, file=paste0("Final_impF", ".rds"))

