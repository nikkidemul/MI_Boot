### 3. MULTIPLE IMPUTATION SCRIPT ###

### We want to analyze Hb as a categorical variable, so we will make it a categorical variable before the imputation procedure. 
incompletedata$Hbcat <- as.factor(ifelse(incompletedata$Hemoglobine < 7.5, 1, ifelse(incompletedata$Hemoglobine>8.5,2,0)))
incompletedata <- incompletedata %>% select(!Hemoglobine)

# steps towards multiple imputation
set.seed(600)
system.time(imp.dry <- mice(incompletedata, maxit=0, printFlag=TRUE, visitsequence='monotone', seed=300))
testimp <- complete(imp.dry)

#check for strong correlations
cormytotaldata <- lapply(incompletedata, as.numeric)
cormytotaldata <- data.frame(cormytotaldata)
mycormatvar <- round(cor(cormytotaldata, use="pairwise.complete.obs"),3)
table(mycormatvar >=0.7 | mycormatvar <=-0.7) #so no correlations higher then 0.7. 
table(mycormatvar >=0.5 | mycormatvar <=-0.5) # so a few are correlated at the 0.5 level, which? 

###check these correlations in excel
write.table(mycormatvar, "correlationmatrix.csv", sep=";")

### the variables correlating > 0.5 are tiffeneau index and fev1 - which you would expect

mymethod <- imp.dry$method
mymethod #check whether all variables are imputed accordingly

#then continue by defining which variables are contributing to the prediction model and post processing if applicable. 

var2impute <- which(mymethod!="")
pred <- imp.dry$predictorMatrix
#don't let tiff and fev1 predict each other? 
pred["fev1_compl", "tiff_compl"] <- 0
pred["tiff_compl", "fev1_compl"] <- 0

system.time(imp.dry2 <- mice(incompletedata, meth=mymethod, pred=pred, maxit=0, printFlag =T, visitSequence = 'monotone', seed=300))
testimp2 <- complete(imp.dry2)
rm(imp.dry)

#postprocessing
postproc <- imp.dry2$post
postproc["fev1_compl"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(60,150))"
postproc["tiff_compl"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(60,150))"

m <- 10
maxit <- 10 

system.time(imp.dry3 <- mice(incompletedata, meth=mymethod, pred=pred, maxit=0, printFlag=T, visitSequence = 'monotone', seed=300))
test3 <- complete(imp.dry3)

#Run final imputation model
system.time(imputedset <- mice(data=incompletedata, 
                               m=m, 
                               maxit=maxit, 
                               method=mymethod, 
                               predictorMatrix = pred, 
                               printFlag=T, 
                               seed=300, 
                               post=postproc, 
                               visitSequence = 'monotone'))


#Then check your imputations
plot(imputedset)

densityplot(imputedset)























