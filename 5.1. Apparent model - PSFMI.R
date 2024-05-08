### 5.1. Apparent model - PSFMI ### 

# First we develop the model using the PSFMI package 

modelD1_back <- psfmi_lr(data=Myimp,
                         formula=outcome~age
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
                         +tiff_compl,
                         nimp=10,
                         impvar=".imp",
                         p.crit=0.157, method="D1", direction="BW")

modelD1_back$RR_model_final


