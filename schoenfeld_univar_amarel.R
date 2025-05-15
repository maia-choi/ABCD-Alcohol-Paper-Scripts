
#Maia Choi
#abcd alc analyses
#written for amarel server###
###script to check for ph assumption/ schoenfield residuals###
##has sip,drink, progression for the univariate analyses###

#start script#

#read in packages
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

#####Pheno models (imp, envrio, cog)####
#loop through univariate predictors
phenos <- c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
            "upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance",
            "coi_nat","coi_ed","coi_se",
            "fes_y_ss_fc","pmq_y_ss_mean","fhx_su_bi", "macv_p_ss_r","prenat_su_exp",
            "srpf_y_ss_ses","srpf_y_ss_iiss","srpf_y_ss_dfs",
            "cog_score")

reg_out<- list()
pvals<- list()
results_table<- NA
for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +",y,"+ factor(sex)+ factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)  
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[6,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"sip_univariate_schoenfield.csv")

####PGS models#####
#loop through pgs
phenos <- c("alcC_resid","EXT_resid","DEP_resid","alcP_r_resid")
reg_out<- list()
pvals<- list()
results_table<- NA
for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_use, sip_ever) ~ +",y,"+ factor(sex)+ factor(income)+ pedu_c+factor(ancestry)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)  
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[6,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"sip_pgs_schoenfield.csv")

##########First Drink###################

#####Pheno models (imp, envrio, cog)####
#loop through univariate predictors
phenos <- c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
            "upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance",
            "coi_nat","coi_ed","coi_se",
            "fes_y_ss_fc","pmq_y_ss_mean","fhx_su_bi", "macv_p_ss_r","prenat_su_exp",
            "srpf_y_ss_ses","srpf_y_ss_iiss","srpf_y_ss_dfs",
            "cog_score")

reg_out<- list()
pvals<- list()
results_table<- NA
for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)  
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[6,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"drink_univariate_schoenfield.csv")

####PGS models#####
#loop through pgs
phenos <- c("alcC_resid","EXT_resid","DEP_resid","alcP_r_resid")
reg_out<- list()
pvals<- list()
results_table<- NA
for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income)+ pedu_c+factor(ancestry)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat)  
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[6,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"drink_pgs_schoenfield.csv")

########## Progression from Sippy Sip to Drinky Drink ###################

#filter to ever had a sip
dat_drink <- dat %>%
  filter(sip_ever == 1)

# create age difference variable
dat_drink$age_sip_drink_use<- dat_drink$age_drink_use- dat_drink$age_sip_use
dat_drink$age_sip_drink_use[dat_drink$age_sip_drink_use<0]<-0

#####Pheno models (imp, envrio, cog)####
#loop through univariate predictors
phenos <- c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
            "upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance",
            "coi_nat","coi_ed","coi_se",
            "fes_y_ss_fc","pmq_y_ss_mean","fhx_su_bi", "macv_p_ss_r","prenat_su_exp",
            "srpf_y_ss_ses","srpf_y_ss_iiss","srpf_y_ss_dfs",
            "cog_score")

reg_out<- list()
pvals<- list()
results_table<- NA
for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income) +factor(race_ethnicity) + pedu_c+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat_drink)  
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[6,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"prog_univariate_schoenfield.csv")

####PGS models#####
#loop through pgs
phenos <- c("alcC_resid","EXT_resid","DEP_resid","alcP_r_resid")
reg_out<- list()
pvals<- list()
results_table<- NA
for (i in 1:length(phenos)) {
  cat("Processing iteration:", i, "\n")
  y <- phenos[i]
  form <- formula(paste0("Surv(age_sip_drink_use, drink_ever) ~ +",y,"+ factor(sex)+ factor(income)+ pedu_c+factor(ancestry)+(1|site_id_l/fid)"))
  reg_out[[i]] <- coxme(form, data = dat_drink)  
  pvals[[i]]<-cox.zph(reg_out[[i]])
  ph_results <-("schoenfield_pvalue"= (pvals[[i]]$table[6,3]))
  results_table <- rbind(results_table, as.data.frame(ph_results))
}

write.csv(results_table,"prog_pgs_schoenfield.csv")

