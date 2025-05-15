
#Maia Choi 11/05####
#abcd alc analyses
#written for amarel server###
###script to run neuro sip analyses###
#####

#start script#

#read in packages
library(survival)
library(dplyr)
library(coxme)
# read in data
pheno <- read.csv("pheno_10_04.csv") 
enviro <- read.csv("enviro.csv") 
neuro <- read.csv("neuro.csv") 
neuro <- subset(neuro, select = -c(smri_vol_scs_lesionrh, smri_vol_scs_lesionlh, smri_vol_scs_wmhintlh, smri_vol_scs_wmhintrh))

#merge your data
dat <- merge(pheno, enviro, by = c("src_subject_id"), all.x = TRUE)
dat <- merge(dat, neuro, by = c("src_subject_id"), all.x = TRUE)

#set refs for regressions
dat$pedu_c<- relevel(factor(dat$pedu_c), ref="low")
dat$income<- relevel(factor(dat$income), ref="lower")
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

#something with neuro stuff
neuro <- subset(neuro, select = -c(smri_vol_cdk_total, smri_vol_scs_wholeb))

######## first sip###########################

# ########CORTICAL regions#################
#include total brain volume  as a covariate
# i had this pull the column names as my predictor list
cortical_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_cdk")]
cort_stxcov <- list()
for (i in 1:length(cortical_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- cortical_stx[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +scale(",y,")+ smri_vol_cdk_total+ factor(sex)+ pb_val +
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  cort_stxcov [[i]] <- coxme(form, data = dat)
}

# list the names of each predictor/ covar
var_names <- c("pheno", "smri_vol_cdk_total", "sex", "pb_val",
               "income higher", "income middle", "black",
               "hispanic", "other", "pedu high",
               "pedu mid", "pedu super high")
results_table <- NULL # initialize empty df

#loops through for as many predictors/ levels u have
num_vars <- length(var_names)

for (i in 1:length(cort_stxcov)) {
  cort_results <- data.frame("n" = cort_stxcov[[i]]$n[2])

  #define which stats to pull  # can modify to pull any stat from mode
  for (j in 1:num_vars) {
    cort_results[paste0("HR_",var_names[j])] <- summary(cort_stxcov[[i]])$coefficients[j,2]
    cort_results[paste0("HR lower_",var_names[j])] <- exp(confint(cort_stxcov[[i]]))[j, 1]
    cort_results[paste0("HR upper_",var_names[j])] <- exp(confint(cort_stxcov[[i]]))[j, 2]
    cort_results[paste0("p_",var_names[j])] <- summary(cort_stxcov[[i]])$coefficients[j, 5]
  }

  results_table <- rbind(results_table, cort_results)
}
results_table$pheno <- cortical_stx

write.csv(results_table,"sip_cortical_11_07.csv")

# ####SUBSCORTICAL####################
# 
subcort_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_scs")]
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

sc_stxcov <- list()
for (i in 1:length(subcort_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- subcort_stx[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +scale(",y,")+ smri_vol_scs_wholeb + factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  sc_stxcov [[i]] <- coxme(form, data = dat)
}

# results_table <- NULL # initialize empty df
#
# # Define variable names and indices
# var_names <- c("pheno", "smri_vol_scs_wholeb", "sex", "pb_val",
#                "income higher", "income middle", "black",
#                "hispanic", "other", "pedu high",
#                "pedu mid", "pedu super high")
#
# num_vars <- length(var_names)
#
# for (i in 1:length(sc_stxcov)) {
#   subcort_results <- data.frame("n" = sc_stxcov[[i]]$n[2])
#
#   for (j in 1:num_vars) {
#     subcort_results[paste0("HR_",var_names[j])] <- summary(sc_stxcov[[i]])$coefficients[j,2]
#     subcort_results[paste0("HR lower_",var_names[j])] <- exp(confint(sc_stxcov[[i]]))[j, 1]
#     subcort_results[paste0("HR upper_",var_names[j])] <- exp(confint(sc_stxcov[[i]]))[j, 2]
#     subcort_results[paste0("p_",var_names[j])] <- summary(sc_stxcov[[i]])$coefficients[j, 5]
#   }
#
#   results_table <- rbind(results_table, subcort_results)
# }
# results_table$pheno <- subcort_stx
#
# write.csv(results_table,"sip_subcortical_11_07.csv")
#
##############TASK BASED ######################################

task_var <- colnames(neuro)[startsWith(colnames(neuro), "tfmri")]
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

task_outcov <- list()
for (i in 1:length(task_var)) {
  cat("Processing iteration:", i, "\n")
  y <- task_var[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +scale(",y,")+ factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  task_outcov [[i]] <- coxme(form, data = dat)
}
results_table <- NULL # initialize empty df

# Define variable names and indices
var_names <- c("pheno", "sex", "pb_val",
               "income higher", "income middle", "black",
               "hispanic", "other", "pedu high",
               "pedu mid", "pedu super high")

num_vars <- length(var_names)

for (i in 1:length(task_outcov)) {
  task_results <- data.frame("n" = task_outcov[[i]]$n[2])

  for (j in 1:num_vars) {
    task_results[paste0("HR_",var_names[j])] <- summary(task_outcov[[i]])$coefficients[j,2]
    task_results[paste0("HR lower_",var_names[j])] <- exp(confint(task_outcov[[i]]))[j, 1]
    task_results[paste0("HR upper_",var_names[j])] <- exp(confint(task_outcov[[i]]))[j, 2]
    task_results[paste0("p_",var_names[j])] <- summary(task_outcov[[i]])$coefficients[j, 5]
  }

  results_table <- rbind(results_table, task_results)
}
results_table$pheno <- task_var

write.csv(results_table,"sip_taskbased_11_07.csv")


######## first DRINK###########################

########CORTICAL regions#################
#include total brain volume  as a covariate
# i had this pull the column names as my predictor list
cortical_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_cdk")]

cort_stxcov <- list()
for (i in 1:length(cortical_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- cortical_stx[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_cdk_total+ factor(sex)+ pb_val +
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  cort_stxcov [[i]] <- coxme(form, data = dat)
}

# list the names of each predictor/ covar
var_names <- c("pheno", "smri_vol_cdk_total", "sex", "pb_val",
               "income higher", "income middle", "black",
               "hispanic", "other", "pedu high",
               "pedu mid", "pedu super high")
results_table <- NULL # initialize empty df

#loops through for as many predictors/ levels u have
num_vars <- length(var_names)

for (i in 1:length(cort_stxcov)) {
  cort_results <- data.frame("n" = cort_stxcov[[i]]$n[2])

  #define which stats to pull  # can modify to pull any stat from mode
  for (j in 1:num_vars) {
    cort_results[paste0("HR_",var_names[j])] <- summary(cort_stxcov[[i]])$coefficients[j,2]
    cort_results[paste0("HR lower_",var_names[j])] <- exp(confint(cort_stxcov[[i]]))[j, 1]
    cort_results[paste0("HR upper_",var_names[j])] <- exp(confint(cort_stxcov[[i]]))[j, 2]
    cort_results[paste0("p_",var_names[j])] <- summary(cort_stxcov[[i]])$coefficients[j, 5]
  }

  results_table <- rbind(results_table, cort_results)
}
results_table$pheno <- cortical_stx

write.csv(results_table,"drink_cortical_11_07.csv")

####SUBSCORTICAL####################

subcort_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_scs")]
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

sc_stxcov <- list() 
for (i in 1:length(subcort_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- subcort_stx[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_scs_wholeb + factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  sc_stxcov [[i]] <- coxme(form, data = dat)  
}

results_table <- NULL # initialize empty df

# Define variable names and indices
var_names <- c("pheno", "smri_vol_scs_wholeb", "sex", "pb_val", 
               "income higher", "income middle", "black", 
               "hispanic", "other", "pedu high", 
               "pedu mid", "pedu super high")

num_vars <- length(var_names)

for (i in 1:length(sc_stxcov)) {
  subcort_results <- data.frame("n" = sc_stxcov[[i]]$n[2])
  
  for (j in 1:num_vars) {
    subcort_results[paste0("HR_",var_names[j])] <- summary(sc_stxcov[[i]])$coefficients[j,2]
    subcort_results[paste0("HR lower_",var_names[j])] <- exp(confint(sc_stxcov[[i]]))[j, 1]
    subcort_results[paste0("HR upper_",var_names[j])] <- exp(confint(sc_stxcov[[i]]))[j, 2]
    subcort_results[paste0("p_",var_names[j])] <- summary(sc_stxcov[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, subcort_results)
}
results_table$pheno <- subcort_stx

write.csv(results_table,"drink_subcortical_11_07.csv")


##############TASK BASED ######################################

task_var <- colnames(neuro)[startsWith(colnames(neuro), "tfmri")]
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

task_outcov <- list() 
for (i in 1:length(task_var)) {
  cat("Processing iteration:", i, "\n")
  y <- task_var[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +scale(",y,")+ factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  task_outcov [[i]] <- coxme(form, data = dat)  
}
results_table <- NULL # initialize empty df

# Define variable names and indices
var_names <- c("pheno", "sex", "pb_val", 
               "income higher", "income middle", "black", 
               "hispanic", "other", "pedu high", 
               "pedu mid", "pedu super high")

num_vars <- length(var_names)

for (i in 1:length(task_outcov)) {
  task_results <- data.frame("n" = task_outcov[[i]]$n[2])
  
  for (j in 1:num_vars) {
    task_results[paste0("HR_",var_names[j])] <- summary(task_outcov[[i]])$coefficients[j,2]
    task_results[paste0("HR lower_",var_names[j])] <- exp(confint(task_outcov[[i]]))[j, 1]
    task_results[paste0("HR upper_",var_names[j])] <- exp(confint(task_outcov[[i]]))[j, 2]
    task_results[paste0("p_",var_names[j])] <- summary(task_outcov[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, task_results)
}
results_table$pheno <- task_var

write.csv(results_table,"drink_taskbased_11_07.csv")


########PROGRESS from sip to drink###########################
dat <- dat %>%
  filter(sip_ever == 1)

# create age difference variable
dat$age_sip_drink_use<- dat$age_drink_use- dat$age_sip_use
dat$age_sip_drink_use[dat$age_sip_drink_use<0]<-0

########CORTICAL regions#################
#include total brain volume  as a covariate
# i had this pull the column names as my predictor list

cortical_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_cdk")]

cort_stxcov <- list() 
for (i in 1:length(cortical_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- cortical_stx[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_cdk_total+ factor(sex)+ 
                         pb_val + factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  cort_stxcov [[i]] <- coxme(form, data = dat)  
}

# list the names of each predictor/ covar 
var_names <- c("pheno", "smri_vol_cdk_total", "sex", "pb_val", 
               "income higher", "income middle", "black", 
               "hispanic", "other", "pedu high", 
               "pedu mid", "pedu super high")
results_table <- NULL # initialize empty df

#loops through for as many predictors/ levels u have
num_vars <- length(var_names)

for (i in 1:length(cort_stxcov)) {
  cort_results <- data.frame("n" = cort_stxcov[[i]]$n[2])
  
  #define which stats to pull  # can modify to pull any stat from mode
  for (j in 1:num_vars) {
    cort_results[paste0("HR_",var_names[j])] <- summary(cort_stxcov[[i]])$coefficients[j,2]
    cort_results[paste0("HR lower_",var_names[j])] <- exp(confint(cort_stxcov[[i]]))[j, 1]
    cort_results[paste0("HR upper_",var_names[j])] <- exp(confint(cort_stxcov[[i]]))[j, 2]
    cort_results[paste0("p_",var_names[j])] <- summary(cort_stxcov[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, cort_results)
}
results_table$pheno <- cortical_stx

write.csv(results_table,"prog_cortical_11_07.csv")

####SUBSCORTICAL####################

subcort_stx <- colnames(neuro)[startsWith(colnames(neuro), "smri_vol_scs")]
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

sc_stxcov <- list() 
for (i in 1:length(subcort_stx)) {
  cat("Processing iteration:", i, "\n")
  y <- subcort_stx[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +scale(",y,")+ smri_vol_scs_wholeb + factor(sex)+ 
                         pb_val+ factor(income) +factor(race_ethnicity) + factor(pedu_c)+(1|site_id_l/fid)"))
  sc_stxcov [[i]] <- coxme(form, data = dat)  
}

results_table <- NULL # initialize empty df

# Define variable names and indices
var_names <- c("pheno", "smri_vol_scs_wholeb", "sex", "pb_val", 
               "income higher", "income middle", "black", 
               "hispanic", "other", "pedu high", 
               "pedu mid", "pedu super high")

num_vars <- length(var_names)

for (i in 1:length(sc_stxcov)) {
  subcort_results <- data.frame("n" = sc_stxcov[[i]]$n[2])
  
  for (j in 1:num_vars) {
    subcort_results[paste0("HR_",var_names[j])] <- summary(sc_stxcov[[i]])$coefficients[j,2]
    subcort_results[paste0("HR lower_",var_names[j])] <- exp(confint(sc_stxcov[[i]]))[j, 1]
    subcort_results[paste0("HR upper_",var_names[j])] <- exp(confint(sc_stxcov[[i]]))[j, 2]
    subcort_results[paste0("p_",var_names[j])] <- summary(sc_stxcov[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, subcort_results)
}
results_table$pheno <- subcort_stx

write.csv(results_table,"prog_subcortical_11_07.csv")


##############TASK BASED ######################################

task_var <- colnames(neuro)[startsWith(colnames(neuro), "tfmri")]
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")

task_outcov <- list() 
for (i in 1:length(task_var)) {
  cat("Processing iteration:", i, "\n")
  y <- task_var[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +scale(",y,")+ factor(sex)+ pb_val+
                         factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  task_outcov [[i]] <- coxme(form, data = dat)  
}
results_table <- NULL # initialize empty df

# Define variable names and indices
var_names <- c("pheno", "sex", "pb_val", 
               "income higher", "income middle", "black", 
               "hispanic", "other", "pedu high", 
               "pedu mid", "pedu super high")

num_vars <- length(var_names)

for (i in 1:length(task_outcov)) {
  task_results <- data.frame("n" = task_outcov[[i]]$n[2])
  
  for (j in 1:num_vars) {
    task_results[paste0("HR_",var_names[j])] <- summary(task_outcov[[i]])$coefficients[j,2]
    task_results[paste0("HR lower_",var_names[j])] <- exp(confint(task_outcov[[i]]))[j, 1]
    task_results[paste0("HR upper_",var_names[j])] <- exp(confint(task_outcov[[i]]))[j, 2]
    task_results[paste0("p_",var_names[j])] <- summary(task_outcov[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, task_results)
}
results_table$pheno <- task_var

write.csv(results_table,"prog_taskbased_11_07.csv")





