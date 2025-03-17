### 3. MULTIPLE IMPUTATION SCRIPT ###

# steps towards multiple imputation
if(MIb.seed) set.seed(600)
system.time(imp.dry <- mice(incompletedata, maxit=0, printFlag=TRUE, visitsequence='monotone', seed=300))
testimp <- complete(imp.dry)

#check for strong correlations
cormytotaldata <- lapply(incompletedata, as.numeric)
cormytotaldata <- data.frame(cormytotaldata)
mycormatvar <- round(cor(cormytotaldata, use="pairwise.complete.obs"),3)
table(mycormatvar >=0.7 | mycormatvar <=-0.7) #so no correlations higher then 0.7. 
table(mycormatvar >=0.5 | mycormatvar <=-0.5) # so a few are correlated at the 0.5 level, which? 

###check these correlations in excel
# write.table(mycormatvar, "correlationmatrix.csv", sep=";")
if(MIb.verbose) View(mycormatvar)

### the variables correlating > 0.5 are tiffeneau index and fev1 - which you would expect

mymethod <- imp.dry$method
if(MIb.intrnl.verbose) print(mymethod)  #check whether all variables are imputed accordingly

#then continue by defining which variables are contributing to the prediction model and post processing if applicable. 

var2impute <- which(mymethod!="")
pred <- imp.dry$predictorMatrix
#don't allow for tiff and fev1 to predict each other: 
pred["fev1_compl", "tiff_compl"] <- 0
pred["tiff_compl", "fev1_compl"] <- 0

system.time(imp.dry2 <- mice(incompletedata, meth=mymethod, pred=pred, maxit=0, printFlag =T, visitSequence = 'monotone', seed=300))
testimp2 <- complete(imp.dry2)
rm(imp.dry)

#postprocessing
postproc <- imp.dry2$post
postproc["fev1_compl"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(60,150))"
postproc["tiff_compl"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(60,150))"

system.time(imp.dry3 <- mice(incompletedata, meth=mymethod, pred=pred, maxit=0, printFlag=T, visitSequence = 'monotone', seed=300))
test3 <- complete(imp.dry3)

#Run final imputation model
if(MIb.verbose) tic <- Sys.time()
imputedset <- mice(data=incompletedata, 
   m=MIb.MICE.m, 
   maxit=MIb.MICE.maxit, 
   method=mymethod, 
   predictorMatrix = pred, 
   print=MIb.intrnl.verbose, 
   seed=300, 
   post=postproc, 
   visitSequence = 'monotone')
if(MIb.verbose) cat("\nMICE imputation completed: ")
if(MIb.verbose) print(Sys.time() - tic)

#Then check your imputations
if(MIb.verbose) plot(imputedset)
if(MIb.verbose) densityplot(imputedset)























