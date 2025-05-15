
#Maia Choi
#abcd alc analyses
#written for amarel server###
###script to run univariate analyses###
##has sip,drink, progression###

#start script#

#read in packages, pray amarel cooperates
library(survival)
library(dplyr)
library(coxme)

#read in 3 data files
pheno <- read.csv("pheno_10_04.csv") 
enviro <- read.csv("enviro.csv") 
geno <- read.csv("geno_wt_resids.csv") 

#merge ur data
dat <- merge(pheno, enviro, by = c("src_subject_id"), all.x = TRUE)
dat <- merge(dat, geno, by = c("src_subject_id"), all.x = TRUE)

#set references for categorical variables
dat$pedu_c<- relevel(factor(dat$pedu_c), ref="low")
dat$income<- relevel(factor(dat$income), ref="lower")
dat$race_ethnicity<- relevel(factor(dat$race_ethnicity), ref="white")


##########First Sip###################
# Covars only model
covar_sip<- coxme(Surv(age_sip_use, sip_ever) ~  factor(sex)+ factor(income) + factor(race_ethnicity)+pedu_c +(1|site_id_l/fid), data = dat)

# list the names of each predictor with levels 4 output table
var_names <- c("sex", "income higher", "income middle", 
               "black","hispanic", "other", 
               "pedu high", "pedu mid", "pedu super high")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

univariate_mod <- data.frame("n" = covar_sip$n[2])

#define which stats to pull
for (j in 1:num_vars) {
  univariate_mod[paste0("HR_",var_names[j])] <- summary(covar_sip)$coefficients[j,2]
  univariate_mod[paste0("HR lower_",var_names[j])] <- confint(covar_sip)[j, 1]
  univariate_mod[paste0("HR upper_",var_names[j])] <- confint(covar_sip)[j, 2]
  univariate_mod[paste0("p_",var_names[j])] <- summary(covar_sip)$coefficients[j, 5]
}

results_table <- rbind(results_table, univariate_mod)

#write output table
write.csv(results_table,"sip_covars_only_11_05.csv")

#####Pheno models (imp, envrio, cog)####
#loop through univariate predictors
phenos <- c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
"upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance",
"coi_nat","coi_ed","coi_se",
"fes_y_ss_fc","pmq_y_ss_mean","fhx_su_bi", "macv_p_ss_r","prenat_su_exp",
"srpf_y_ss_ses","srpf_y_ss_iiss","srpf_y_ss_dfs",
"cog_score")

reg_out<- list()

for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +",y,"+ factor(sex)+ factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat)  
}

# list the names of each predictor with levels 4 output table
var_names <- c("pheno","sex", "income higher", "income middle", 
               "black","hispanic", "other", 
               "pedu high", "pedu mid", "pedu super high")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

for (i in 1:length(reg_out)) {
  univariate_mod <- data.frame("n" = reg_out[[i]]$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(reg_out[[i]])[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(reg_out[[i]])[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, univariate_mod)
}
#make sure to add pheno names to output table
results_table$pheno <- phenos

#write output table
write.csv(results_table,"sip_univariate_11_05.csv")

####PGS models#####
#loop through pgs
phenos <- c("alcC_resid","EXT_resid","DEP_resid","alcP_r_resid")

reg_out<- list()

for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +",y,"+ factor(sex)+ factor(income)+ pedu_c+factor(ancestry)+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat)  
}

# list the names of each predictor with levels 4 output table
var_names <- c("pheno","sex", "income higher", "income middle", 
               "pedu high", "pedu mid", "pedu super high",
               "ancestry")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

for (i in 1:length(reg_out)) {
  univariate_mod <- data.frame("n" = reg_out[[i]]$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(reg_out[[i]])[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(reg_out[[i]])[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, univariate_mod)
}
#make sure to add pheno names to output table
results_table$pheno <- phenos

#write output table
write.csv(results_table,"sip_PGS_11_05.csv")

##########First Drinky Drink ##################

# Covars only model
covar_drink<- coxme(Surv(age_drink_use, drink_ever) ~  factor(sex)+ factor(income) + pedu_c +factor(race_ethnicity)+(1|site_id_l/fid), data = dat)

# list the names of each predictor with levels 4 output table
var_names <- c("sex", "income higher", "income middle", 
               "black","hispanic", "other", 
               "pedu high", "pedu mid", "pedu super high")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

univariate_mod <- data.frame("n" = covar_drink$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(covar_drink)$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(covar_drink)[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(covar_drink)[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(covar_drink)$coefficients[j, 5]
  }
  
results_table <- rbind(results_table, univariate_mod)

#write output table
write.csv(results_table,"drink_covars_only_11_05.csv")

#####Pheno models (imp, envrio, cog)####
#loop through univariate predictors
phenos <- c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
            "upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance",
            "coi_nat","coi_ed","coi_se",
            "fes_y_ss_fc","pmq_y_ss_mean","fhx_su_bi", "macv_p_ss_r","prenat_su_exp",
            "srpf_y_ss_ses","srpf_y_ss_iiss","srpf_y_ss_dfs",
            "cog_score")

reg_out<- list()

for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat)  
}

# list the names of each predictor with levels 4 output table
var_names <- c("pheno","sex", "income higher", "income middle", 
               "black","hispanic", "other", 
               "pedu high", "pedu mid", "pedu super high")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

for (i in 1:length(reg_out)) {
  univariate_mod <- data.frame("n" = reg_out[[i]]$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(reg_out[[i]])[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(reg_out[[i]])[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, univariate_mod)
}
#make sure to add pheno names to output table
results_table$pheno <- phenos

#write output table
write.csv(results_table,"drink_univariate_11_05.csv")

####PGS models#####
#loop through pgs
phenos <- c("alcC_resid","EXT_resid","DEP_resid","alcP_r_resid")

reg_out<- list()

for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income)+ pedu_c+factor(ancestry)+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat)  
}

# list the names of each predictor with levels 4 output table
var_names <- c("pheno","sex", "income higher", "income middle", 
               "pedu high", "pedu mid", "pedu super high",
               "ancestry")

results_table <- NULL # initialize empty df

num_vars <- length(var_names) #loops thru predictors/levels

for (i in 1:length(reg_out)) {
  univariate_mod <- data.frame("n" = reg_out[[i]]$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(reg_out[[i]])[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(reg_out[[i]])[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, univariate_mod)
}
#make sure to add pheno names to output table
results_table$pheno <- phenos

#write output table
write.csv(results_table,"drink_PGS_11_05.csv")

########## Progression from Sippy Sip to Drinky Drink ###################

#filter to ever had a sip
dat_drink <- dat %>%
  filter(sip_ever == 1)

# create age difference variable
dat_drink$age_sip_drink_use<- dat_drink$age_drink_use- dat_drink$age_sip_use
dat_drink$age_sip_drink_use[dat_drink$age_sip_drink_use<0]<-0

# Covars only model
covar_prog<- coxme(Surv(age_drink_use, drink_ever) ~  factor(sex)+ factor(income) + pedu_c +factor(race_ethnicity)+(1|site_id_l/fid), data = dat_drink)

# list the names of each predictor with levels 4 output table
var_names <- c("sex", "income higher", "income middle", 
               "black","hispanic", "other", 
               "pedu high", "pedu mid", "pedu super high")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

univariate_mod <- data.frame("n" = covar_prog$n[2])

#define which stats to pull
for (j in 1:num_vars) {
  univariate_mod[paste0("HR_",var_names[j])] <- summary(covar_prog)$coefficients[j,2]
  univariate_mod[paste0("HR lower_",var_names[j])] <- confint(covar_prog)[j, 1]
  univariate_mod[paste0("HR upper_",var_names[j])] <- confint(covar_prog)[j, 2]
  univariate_mod[paste0("p_",var_names[j])] <- summary(covar_prog)$coefficients[j, 5]
}

results_table <- rbind(results_table, univariate_mod)

#write output table
write.csv(results_table,"prog_covars_only_11_05.csv")

#####Pheno models (imp, envrio, cog)####
#loop through univariate predictors
phenos <- c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
            "upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance",
            "coi_nat","coi_ed","coi_se",
            "fes_y_ss_fc","pmq_y_ss_mean","fhx_su_bi", "macv_p_ss_r","prenat_su_exp",
            "srpf_y_ss_ses","srpf_y_ss_iiss","srpf_y_ss_dfs",
            "cog_score")

reg_out<- list()

for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat_drink)  
}

# list the names of each predictor with levels 4 output table
var_names <- c("pheno","sex", "income higher", "income middle", 
               "black","hispanic", "other", 
               "pedu high", "pedu mid", "pedu super high")

results_table <- NULL # initialize empty df

num_vars <- length(var_names)#loops thru predictors/levels

for (i in 1:length(reg_out)) {
  univariate_mod <- data.frame("n" = reg_out[[i]]$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(reg_out[[i]])[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(reg_out[[i]])[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, univariate_mod)
}
#make sure to add pheno names to output table
results_table$pheno <- phenos

#write output table
write.csv(results_table,"prog_univariate_11_05.csv")

####PGS models#####
#loop through pgs
phenos <- c("alcC_resid","EXT_resid","DEP_resid","alcP_r_resid")

reg_out<- list()

for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income)+ pedu_c+factor(ancestry)+(1|site_id_l/fid)"))
  reg_out [[i]] <- coxme(form, data = dat_drink)  
}

# list the names of each predictor with levels 4 output table
var_names <- c("pheno","sex", "income higher", "income middle", 
               "pedu high", "pedu mid", "pedu super high",
               "ancestry")

results_table <- NULL # initialize empty df

num_vars <- length(var_names) #loops thru predictors/levels

for (i in 1:length(reg_out)) {
  univariate_mod <- data.frame("n" = reg_out[[i]]$n[2])
  
  #define which stats to pull
  for (j in 1:num_vars) {
    univariate_mod[paste0("HR_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j,2]
    univariate_mod[paste0("HR lower_",var_names[j])] <- confint(reg_out[[i]])[j, 1]
    univariate_mod[paste0("HR upper_",var_names[j])] <- confint(reg_out[[i]])[j, 2]
    univariate_mod[paste0("p_",var_names[j])] <- summary(reg_out[[i]])$coefficients[j, 5]
  }
  
  results_table <- rbind(results_table, univariate_mod)
}
#make sure to add pheno names to output table
results_table$pheno <- phenos

#write output table
write.csv(results_table,"prog_PGS_11_05.csv")



