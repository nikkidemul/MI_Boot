### APPARENT PERFORMANCE PSFMI ###

# The performance of the apparent model found in 5.1. can be test using: 
perf_model1 <- pool_performance(data=Myimp, nimp=10, impvar=".imp", formula= outcome ~ age + open + transhiatal + dummy_Neotx3 + comorb_dm +comorb_hypertension, model_type="binomial")

perf_model1
