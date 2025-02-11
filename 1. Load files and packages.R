### Load files and packages ###

rm(list=ls())

library(logistf)
library(mice)
library(tidyverse)
library(psfmi)
library(dplyr)


library(mgcv)
library(formula.tools)
library(Matrix)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
library(stats)
library(mice)
library(tibble)
library(pROC)
library(rms)

#The script is structured as follows: 
### Three scenario's are provided: PSFMI (always denoted as x.1), GLM (manual) unpenalized regression (always denoted as x.2), 
##### and Firth's penalized logistic regression (always denoted as x.3). 
### 2., 2.2., and 2.3. comprise as function library's. 2. contains all general functions used in all script from 3.x onwards, 
##### and 2.2 and 2.3. contain functions that are specific for the GLM and Firth scenarios. 
### In 3. MICE script, script can be found to perform multiple imputation.
### In 4. Dataprep, data is prepped to use for further analysis. For PSFMI, GLM (unpenalized logistic regression) and Firth, data should be in slightly different formats. 
### In 5.1 - 5.3, apparent model, including variable selection over the imputed datasets is scripted for respectively PSFMI, unpenalized logistic regression using GLM,
##### and Firth's penalized logistic regression. Functions used can be found in the corresponding function library's. 
### In 6.1 - 6.2, apparent peformance of the final models formulated in 5.x is estimated for respectively PSFMI, GLM and Firth. 
### In 7.1 - 7.3, the bootstrap procedure and calculation of optimism of performance measure can be found, for (again) respectively PSFMI, GLM and Firth. 
