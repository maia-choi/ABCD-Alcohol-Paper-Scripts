
#Maia Choi 
#abcd alc analyses
#written for amarel server###
###script to run neuro sip analyses###
#####

#start script#

#read in packages, pray amarel cooperates
library(survival)
library(dplyr)
library(coxme)
# read in data
pheno <- read.csv("pheno_10_04.csv")
enviro <- read.csv("enviro.csv")
neuro <- read.csv("neuro.csv")
neuro <- subset(neuro, select = -c(smri_vol_scs_lesionrh, smri_vol_scs_lesionlh, smri_vol_scs_wmhintlh, 
smri_vol_scs_wmhintrh,smri_vol_scs_csf,smri_vol_scs_wmhint,smri_vol_scs_ccps,smri_vol_scs_ccmidps,
smri_vol_scs_ccct,smri_vol_scs_ccmidat,smri_vol_scs_ccat,smri_vol_scs_latventricles,
smri_vol_scs_allventricles,smri_vol_scs_intracranialv,smri_vol_scs_suprateialv,smri_vol_scs_subcorticalgv))

#merge your data
dat <- merge(pheno, enviro, by = c("src_subject_id"), all.x = TRUE)
dat <- merge(dat, neuro, by = c("src_subject_id"), all.x = TRUE)

#set refs for regressions
dat$pedu_c<- relevel(factor(dat$pedu_c), ref="low")
dat$income<- relevel(factor(dat$income), ref="lower")
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

#something with neuro stuff
neuro <- subset(neuro, select = -c(smri_vol_cdk_total, smri_vol_scs_wholeb))

######## first SIP###########################

########CORTICAL regions#################
#include total brain volume  as a covariate
# i had this pull the column names as my predictor list
cortical_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_cdk")]

reg_out <- list()
pvals<- list()
results_table<- NULL

for (i in 1:length(cortical_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- cortical_stx[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +scale(",y,")+ smri_vol_cdk_total+ factor(sex)+ pb_val +
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[8,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"sip_cortical_schoenfield.csv")

####SUBCORTICAL####################

subcort_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_scs")]

reg_out <- list()
pvals<- list()
results_table<- NULL

for (i in 1:length(subcort_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- subcort_stx[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +scale(",y,")+ smri_vol_scs_wholeb + factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[8,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"sip_subcortial_schoenfield.csv")

##############TASK BASED ######################################

task_var <- colnames(neuro)[startsWith(colnames(neuro), "tfmri")]

reg_out <- list()
pvals<- list()
results_table<- NULL
 for (i in 1:length(task_var)) {
   cat("Processing iteration:", i, "\n")
  y <- task_var[i]
   form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +scale(",y,")+ factor(sex)+ pb_val+
                          factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
   reg_out [[i]] <- coxme(form, data = dat)
   pvals[[i]]<-cox.zph(reg_out[[i]])
   ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[7,3]))
   results_table <- rbind(results_table, as.data.frame(ph_results))
 }

write.csv(results_table,"sip_task_schoenfield.csv")

######## first DRINK###########################

########CORTICAL regions#################
#include total brain volume  as a covariate
# i had this pull the column names as my predictor list
cortical_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_cdk")]

reg_out <- list()
pvals<- list()
results_table<- NULL

for (i in 1:length(cortical_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- cortical_stx[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_cdk_total+ factor(sex)+ pb_val +
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[8,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"drink_cortical_schoenfield.csv")

####SUBCORTICAL####################

subcort_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_scs")]

reg_out <- list()
pvals<- list()
results_table<- NULL

for (i in 1:length(subcort_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- subcort_stx[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_scs_wholeb + factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[8,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"drink_subcortial_schoenfield.csv")

##############TASK BASED ######################################

task_var <- colnames(neuro)[startsWith(colnames(neuro), "tfmri")]

reg_out <- list()
pvals<- list()
results_table<- NULL
for (i in 1:length(task_var)) {
  cat("Processing iteration:", i, "\n")
  y <- task_var[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +scale(",y,")+ factor(sex)+ pb_val+
                          factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[7,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"drink_task_schoenfield.csv")

######## PROGRESS from sip to drink###########################
dat <- dat %>%
  filter(sip_ever == 1)

# create age difference variable
dat$age_sip_drink_use<- dat$age_drink_use- dat$age_sip_use
dat$age_sip_drink_use[dat$age_sip_drink_use<0]<-0

########CORTICAL regions#################
#include total brain volume  as a covariate
# i had this pull the column names as my predictor list
cortical_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_cdk")]

reg_out <- list()
pvals<- list()
results_table<- NULL

for (i in 1:length(cortical_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- cortical_stx[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_cdk_total+ factor(sex)+ pb_val +
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[8,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"prog_cortical_schoenfield.csv")

####SUBCORTICAL####################

subcort_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_scs")]

reg_out <- list()
pvals<- list()
results_table<- NULL

for (i in 1:length(subcort_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- subcort_stx[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_scs_wholeb + factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[8,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"prog_subcortial_schoenfield.csv")

##############TASK BASED ######################################

task_var <- colnames(neuro)[startsWith(colnames(neuro), "tfmri")]

reg_out <- list()
pvals<- list()
results_table<- NULL
for (i in 1:length(task_var)) {
  cat("Processing iteration:", i, "\n")
  y <- task_var[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +scale(",y,")+ factor(sex)+ pb_val+
                          factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat)
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[7,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"prog_task_schoenfield.csv")

