##########################################
#Programmer: `Rachel Sendor
#Last update: Mar 2023
#Purpose:     DRC Non-falciparum Descriptive Epidemiology Analysis
#             CUMULATIVE INCIDENCE ESTIMATES 
#             
##########################################



#########

## Analyzing Cumulative Incidence of Pm, Po, and Pf 

#########


options(max.print=999999)
options(scipen = 999)
.libPaths()

#Loading R Packages:
library(tidyverse) 
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggbreak)
library(ggthemes)
library(survminer)
library(purrr)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)


##########################################################
# SET WORKING DIRECTORY

#### blinded file path
setwd(".../Dataset_AD/Results")



##########################################################
# Cumulative Incidence Analyses - Risk Curves - Overall and Stratified by different factors    
##########################################################

##Analysis Notes: 
#All datasets are comprised of subjects negative for the species-specific infection at baseline
#Active + Passive dataset includes the time to first active OR passive infection, only among those negative for infection at baseline 

#Analysis Datasets created in SAS: (program: AnalysisDatasetBuild_DRCNonPf_Mar22_RS.sas)


##########################################################
# IMPORTING DATASETS FOR ANALYSIS --- DATASET CREATED IN SAS (SEE PROGRAM: AnalysisDatasetBuild_DRCNonPf_Mar22_RS.sas;  Updated DRC Incidence Estimation_Jul23)

#Read in incidence-specific analytic datasets   (separate datasets created for each species, and each analysis population (ie 9 total datasets))


###PM INFECTIONS###

#Active only 

#Pm_Risk_Act<-read_excel(".../Pm_Risk_Cohort_Active_Sub_Jul22.xlsx")
Pm_Risk_Act<-read_excel(".../Pm_Risk_Cohort_Active_Sub_Jul22_Jul23.xlsx")

#Active & Passive  -VNP pop
#Pm_Risk_VNP<-read_excel(".../Pm_Risk_Cohort_VNP_Sub_Jul22.xlsx")
Pm_Risk_VNP<-read_excel(".../Pm_Risk_Cohort_VNP_Sub_Jul22_Jul23.xlsx")

#Passive only - changed Jul23 -VNP pop
Pm_Risk_VNP2<-read_excel(".../Pm_Risk_Cohort_VNP_Sub2_Jul22_Jul23.xlsx")

#Active & Passive  - TOtal pop
#Pm_Risk_Tot<-read_excel(".../Pm_Risk_Cohort_Total_Sub_Jul22.xlsx")
Pm_Risk_Tot<-read_excel(".../Pm_Risk_Cohort_Total_Sub_Jul22_Jul23.xlsx")



###PO INFECTIONS###

#Active only 
#Po_Risk_Act<-read_excel(".../Po_Risk_Cohort_Active_Sub_Jul22.xlsx")
Po_Risk_Act<-read_excel(".../Po_Risk_Cohort_Active_Sub_Jul22_Jul23.xlsx")

#Active & Passive  
#Po_Risk_VNP<-read_excel(".../Po_Risk_Cohort_VNP_Sub_Jul22.xlsx")
Po_Risk_VNP<-read_excel(".../Po_Risk_Cohort_VNP_Sub_Jul22_Jul23.xlsx")

#Passive only - changed Jul23 -VNP pop
Po_Risk_VNP2<-read_excel(".../Po_Risk_Cohort_VNP_Sub2_Jul22_Jul23.xlsx")

#Active & Passive  - TOtal pop
#Po_Risk_Tot<-read_excel(".../Po_Risk_Cohort_Total_Sub_Jul22.xlsx")
Po_Risk_Tot<-read_excel(".../Po_Risk_Cohort_Total_Sub_Jul22_Jul23.xlsx")



###PF INFECTIONS###

#Active only 
#Pf_Risk_Act<-read_excel(".../Pf_Risk_Cohort_Active_Sub_Jul22.xlsx")
Pf_Risk_Act<-read_excel(".../Pf_Risk_Cohort_Active_Sub_Jul22_Jul23.xlsx")

#Active & Passive  
#Pf_Risk_VNP<-read_excel(".../Pf_Risk_Cohort_VNP_Sub_Jul22.xlsx")
Pf_Risk_VNP<-read_excel(".../Pf_Risk_Cohort_VNP_Sub_Jul22_Jul23.xlsx")

#Passive only - changed Jul23 -VNP pop
Pf_Risk_VNP2<-read_excel(".../Pf_Risk_Cohort_VNP_Sub2_Jul22_Jul23.xlsx")

#Active & Passive  - TOtal pop
#Pf_Risk_Tot<-read_excel(".../Pf_Risk_Cohort_Total_Sub_Jul22.xlsx")
Pf_Risk_Tot<-read_excel(".../Pf_Risk_Cohort_Total_Sub_Jul22_Jul23.xlsx")





##########################################################

###DERIVING NEW VARIABLES FOR ANALYSIS###

##  Age_cat2 = baseline age: 1 = <5yrs;  2 = 5-14 yrs;  3= 15+ yrs
##  Time_active_mon = time to event (or censoring) in months 



Pm_Risk_Act<- Pm_Risk_Act%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_active_mon=Time_active/30.4)



Pm_Risk_VNP<- Pm_Risk_VNP%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)%>%
  mutate(Time_all_end_mon=Time_all_end/30.4)

Pm_Risk_VNP2<- Pm_Risk_VNP2%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                      ifelse(age_cat_FU0==2, 2, 
                                                             ifelse(age_cat_FU0==3, 2, 
                                                                    ifelse(age_cat_FU0==4, 3, 
                                                                           ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)


Pm_Risk_Tot<- Pm_Risk_Tot%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)%>%
  mutate(Time_all_end_mon=Time_all_end/30.4)%>%
  mutate(Time_all_end_mon_2=Time_all_end_2/30.4)%>%
  mutate(Pm_Event_All_CE= ifelse(Pm_Event_All==0, 0, 
                                 ifelse(Pm_Event_All==1, 1, 2)))

Po_Risk_Act<- Po_Risk_Act%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_active_mon=Time_active/30.4)

Po_Risk_VNP<- Po_Risk_VNP%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)%>%
  mutate(Time_all_end_mon=Time_all_end/30.4)

Po_Risk_VNP2<- Po_Risk_VNP2%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                      ifelse(age_cat_FU0==2, 2, 
                                                             ifelse(age_cat_FU0==3, 2, 
                                                                    ifelse(age_cat_FU0==4, 3, 
                                                                           ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)                


Po_Risk_Tot<- Po_Risk_Tot%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)%>%
  mutate(Time_all_end_mon=Time_all_end/30.4)%>%
  mutate(Time_all_end_mon_2=Time_all_end_2/30.4)%>%
  mutate(Po_Event_All_CE= ifelse(Po_Event_All==0, 0, 
                                 ifelse(Po_Event_All==1, 1, 2)))


Pf_Risk_Act<- Pf_Risk_Act%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_active_mon=Time_active/30.4)

Pf_Risk_VNP<- Pf_Risk_VNP%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)%>%
  mutate(Time_all_end_mon=Time_all_end/30.4)

Pf_Risk_VNP2<- Pf_Risk_VNP2%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                      ifelse(age_cat_FU0==2, 2, 
                                                             ifelse(age_cat_FU0==3, 2, 
                                                                    ifelse(age_cat_FU0==4, 3, 
                                                                           ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)                  

Pf_Risk_Tot<- Pf_Risk_Tot%>%mutate(age_cat2= ifelse(age_cat_FU0==1, 1, 
                                                    ifelse(age_cat_FU0==2, 2, 
                                                           ifelse(age_cat_FU0==3, 2, 
                                                                  ifelse(age_cat_FU0==4, 3, 
                                                                         ifelse(age_cat_FU0==5, 3, NA))))))%>%
  mutate(Time_VNP_end_mon=Time_VNP_end/30.4)%>%
  mutate(Time_all_end_mon=Time_all_end/30.4)%>%
  mutate(Time_all_end_mon_2=Time_all_end_2/30.4)%>%
  mutate(Pf_Event_All_CE= ifelse(Pf_Event_All==0, 0, 
                                 ifelse(Pf_Event_All==1, 1, 2)))





##########################################################
# COMBINED SURVIVAL CURVES OF PM, PO, AND PF INFECTION    
##########################################################

#Estimate the Kaplan Meier Crude Risk Function for All Datasets, where event = 1 and Time = time to event or censoring due to LTFU 
##Estimating cumulative incidence = 1-KM survival function (assuming no competing risks and non-informative censoring)

Pf_Act_Surv<-survfit(Surv(Time_active_mon, Pf_Event_Active) ~ 1, data = Pf_Risk_Act)
Pm_Act_Surv<-survfit(Surv(Time_active_mon, Pm_Event_Active) ~ 1, data = Pm_Risk_Act)
Po_Act_Surv<-survfit(Surv(Time_active_mon, Po_Event_Active) ~ 1, data = Po_Risk_Act)


fit <- list(pfs=Pf_Act_Surv, pos=Po_Act_Surv, pms=Pm_Act_Surv)
ggsurvplot(fit, combine = TRUE, 
           legend.title = "Survey-based Pop.",
           #  legend.labs = c(""),
           break.x.by = 6, 
           # Add p-value and tervals
           # pval = TRUE,
           conf.int = TRUE,
           #surv.median.line = c("hv"),
           censor = TRUE, 
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.17,
           cumevents = TRUE, 
           # cumcensor = TRUE, 
           #cumcensor.title = "Cumulative number censored",
           cumevent.title = "Cumulative infections",
           #risk.table.fontsize = 2.5,
           # font.x=2.5,
           tables.theme = theme_cleantable(),
           xlab = "Months from Baseline", 
           ylab = ("Prob. un-infected"), 
           xlim = c(0,24))



#Crude Pm Risk among Active Follow-ups






### FOR TABLE 2:  

##Look at Time in days and months to compare (months is likely easier to interpret in plots)
km <- survfit(Surv(Time_active_mon, Pm_Event_Active) ~ 1, data=Pm_Risk_Act, stype=2, id=Subject_ID)
delta_risk_Pm_act <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 
view(delta_risk_Pm_act)


####GRAPH CUMULATIVE INCIDENCE CURVES TO VISUALIZE 
ggplot(data = delta_risk_Pm_act) + 
  geom_step(aes(x = time, y = risk), color = "turquoise", size=1.4)+
  
  
  geom_ribbon(aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "turquoise", alpha=0.3)+
  ylab("Cumulative Incidence") +
  xlab("Months from Baseline")+
  labs(title="Risk of Pm Infxn during Active Visits")+
  theme_clean()



####CONFIRM NUMBER AT RISK AND CUMULATIVE EVENTS FOR TABLE 2 
ggsurvplot(
  fit = survfit(Surv(Time_active_mon, Pm_Event_Active) ~ 1, data = Pm_Risk_Act, stype=2, id=Subject_ID), 
  legend.title = "P.malariae",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,20))

#Summary with 95% CIs for survival estimates 
summary(survfit(Surv(Time_active_mon, Pm_Event_Active) ~ 1, data = Pm_Risk_Act, stype=2, id=Subject_ID))

#__________________________________________________________________________________________________#

#Crude Po Risk among Active Follow-ups

km <- survfit(Surv(Time_active_mon, Po_Event_Active) ~ 1, data=Po_Risk_Act, stype=2, id=Subject_ID)
delta_risk_Po_act <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 

view(delta_risk_Po_act)


ggplot(data = delta_risk_Po_act) + 
  geom_step(aes(x = time, y = risk), color = "red", size=1.4)+
  
  
  geom_ribbon(aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.3)+
  ylab("Cumulative Incidence") +
  xlab("Months from Baseline")+
  labs(title="Risk of Po Infxn during Active Visits")+
  theme_clean()

ggsurvplot(
  fit = survfit(Surv(Time_active_mon, Po_Event_Active) ~ 1, data = Po_Risk_Act, stype=2, id=Subject_ID), 
  legend.title = "P.ovale",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,20))


#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_active_mon, Po_Event_Active) ~ 1, data = Po_Risk_Act))
summary(delta_risk_Po_act)


#__________________________________________________________________________________________________#

#Crude Pf Risk among Active Follow-ups

#Summary with 95% CIs for survival estimates & cumulative incidence estimates
km <- survfit(Surv(Time_active_mon, Pf_Event_Active) ~ 1, data=Pf_Risk_Act, stype=2, id=Subject_ID)
summary(km)

delta_risk_Pf_act <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 
view(delta_risk_Pf_act)
summary(delta_risk_Pf_act)


ggplot(data = delta_risk_Pf_act) + 
  geom_step(aes(x = time, y = risk), color = "red", size=1.4)+
  
  
  geom_ribbon(aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "orange", alpha=0.3)+
  ylab("Cumulative Incidence") +
  xlab("Months from Baseline")+
  labs(title="Risk of Pf Infxn during Active Visits")+
  theme_clean()

ggsurvplot(
  fit = survfit(Surv(Time_active_mon, Pf_Event_Active) ~ 1, data = Pf_Risk_Act), 
  legend.title = "P.falciparum",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=14,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Probability un-infected"), 
  xlim = c(0,20))



#_________________________________________________________________________________________________#

# Pm risk in passive visits 

km <- survfit(Surv(Time_VNP_end_mon, Pm_event_VNP) ~ 1, data=Pm_Risk_VNP, stype=2, id=Subject_ID)
delta_risk_Pm_vnp <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 
view(delta_risk_Pm_vnp)


ggplot(data = delta_risk_Pm_vnp) + 
  geom_step(aes(x = time, y = risk), color = "blue", size=1.4)+
  ylab("Risk") +
  labs(title="Risk of Pm Infxn during Passive Visits")

ggsurvplot(
  fit = survfit(Surv(Time_VNP_end_mon, Pm_event_VNP) ~ 1, data = Pm_Risk_VNP, stype=2, id=Subject_ID), 
  legend.title = "P.malariae",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=14,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Prob. un-infected"), 
  xlim = c(0,36))

#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_VNP_end_mon, Pm_event_VNP) ~ 1, data = Pm_Risk_VNP, stype=2, id=Subject_ID))
summary(delta_risk_Pm_vnp)


###NOW with just Passive visits, and passive infections 

km <- survfit(Surv(Time_VNP_end_mon, Pm_event_VNP) ~ 1, data=Pm_Risk_VNP2, stype=2, id=Subject_ID)
delta_risk_Pm_vnp2 <- data.frame(time = summary(km)$time, 
                                 surv = summary(km)$surv, 
                                 risk = 1 - summary(km)$surv,
                                 LB_CI = 1 - summary(km)$upper,
                                 UB_CI = 1 - summary(km)$lower,
                                 risk_percent = (1 - summary(km)$surv)*100, 
                                 LB_CI_percent = (1- summary(km)$upper)*100,
                                 UB_CI_percent = (1- summary(km)$lower)*100, 
                                 nrisk = summary(km)$n.risk, 
                                 nevent = summary(km)$n.event, 
                                 n_censor = summary(km)$n.censor) 

view(delta_risk_Pm_vnp2)



ggsurvplot(
  fit = survfit(Surv(Time_VNP_end_mon, Pm_event_VNP) ~ 1, data = Pm_Risk_VNP2, stype=2, id=Subject_ID), 
  legend.title = "P.malariae",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=14,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Prob. un-infected"), 
  xlim = c(0,36))

#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_VNP_end_mon, Pm_event_VNP) ~ 1, data = Pm_Risk_VNP2, stype=2, id=Subject_ID))





#_________________________________________________________________________________________________#
# Po risk in passive visits 

km <- survfit(Surv(Time_VNP_end_mon, Po_event_VNP) ~ 1, data=Po_Risk_VNP, stype=2, id=Subject_ID)
delta_risk_Po_vnp <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor)  
view(delta_risk_Po_vnp)

ggplot(data = delta_risk_Po_vnp) + 
  geom_step(aes(x = time, y = risk), color = "red", size=1.4)+
  ylab("Risk") +
  labs(title="Risk of Po Infxn during Passive Visits")

ggsurvplot(
  fit = survfit(Surv(Time_VNP_end_mon, Po_event_VNP) ~ 1, data = Po_Risk_VNP, stype=2, id=Subject_ID), 
  legend.title = "P.ovale",
  legend.labs = c(""),
  break.x.by = 6, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,36))

#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_VNP_end_mon, Po_event_VNP) ~ 1, data = Po_Risk_VNP, stype=2, id=Subject_ID))
summary(delta_risk_Po_vnp)



#_________________________________________________________________________________________________#
# Pf risk in passive visits 

km <- survfit(Surv(Time_VNP_end_mon, Pf_event_VNP) ~ 1, data=Pf_Risk_VNP, stype=2, id=Subject_ID)
delta_risk_Pf_vnp <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 
view(delta_risk_Pf_vnp)

ggplot(data = delta_risk_Pf_vnp) + 
  geom_step(aes(x = time, y = risk), color = "orange", size=1.4)+
  ylab("Risk") +
  labs(title="Risk of Pf Infxn during Passive Visits")

ggsurvplot(
  fit = survfit(Surv(Time_VNP_end_mon, Pf_event_VNP) ~ 1, data = Pf_Risk_VNP, stype=2, id=Subject_ID), 
  legend.title = "P.falciparum",
  legend.labs = c(""),
  break.x.by = 6, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,36))

#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_VNP_end_mon, Pf_event_VNP) ~ 1, data = Pf_Risk_VNP, stype=2, id=Subject_ID))
summary(delta_risk_Pf_vnp)




#_________________________________________________________________________________________________#
# TOTAL STUDY POPULATION 


# Pm risk in total study visits 

km <- survfit(Surv(Time_all_end_mon, Pm_Event_All) ~ 1, data=Pm_Risk_Tot, stype=2, id=Subject_ID)
delta_risk_Pm_tot <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor)  
view(delta_risk_Pm_tot)


km <- survfit(Surv(Time_all_end_mon_2, Pm_Event_All) ~ 1, data=Pm_Risk_Tot, stype=2, id=Subject_ID)
delta_risk_Pm_tot_2 <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor)  
view(delta_risk_Pm_tot_2)



ggplot(data = delta_risk_Pm_tot) + 
  geom_step(aes(x = time, y = risk), color = "blue", size=1.4)+
  ylab("Risk") +
  labs(title="Risk of Pm Infxn during Active + Passive Visits")


ggsurvplot(
  fit = survfit(Surv(Time_all_end_mon, Pm_Event_All) ~ 1, data = Pm_Risk_Tot, stype=2, id=Subject_ID), 
  legend.title = "P.malariae",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,36))



ggsurvplot(
  fit = survfit(Surv(Time_all_end_mon_2, Pm_Event_All) ~ 1, data = Pm_Risk_Tot, stype=2, id=Subject_ID), 
  legend.title = "P.malariae",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,36))






#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_all_end_mon, Pm_Event_All) ~ 1, data = Pm_Risk_Tot))
summary(survfit(Surv(Time_all_end_mon_2, Pm_Event_All) ~ 1, data = Pm_Risk_Tot))



#_________________________________________________________________________________________________#
# Po risk in total study visits 

km <- survfit(Surv(Time_all_end_mon, Po_Event_All) ~ 1, data=Po_Risk_Tot, stype=2, id=Subject_ID)
delta_risk_Po_tot <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 
view(delta_risk_Po_tot)



ggplot(data = delta_risk_Po_tot) + 
  geom_step(aes(x = time, y = risk), color = "red", size=1.4)+
  ylab("Risk") +
  labs(title="Risk of Po Infxn during Active + Passive Visits")



ggsurvplot(
  fit = survfit(Surv(Time_all_end_mon, Po_Event_All) ~ 1, data = Po_Risk_Tot, stype=2, id=Subject_ID), 
  legend.title = "p. ovale",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,36))


#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_all_end_mon, Po_Event_All) ~ 1, data = Po_Risk_Tot))
summary(delta_risk_Po_act)





#_________________________________________________________________________________________________#
# Pf risk in total study visits 

km <- survfit(Surv(Time_all_end_mon, Pf_Event_All) ~ 1, data=Pf_Risk_Tot, stype=2, id=Subject_ID)
delta_risk_Pf_tot <- data.frame(time = summary(km)$time, 
                                surv = summary(km)$surv, 
                                risk = 1 - summary(km)$surv,
                                LB_CI = 1 - summary(km)$upper,
                                UB_CI = 1 - summary(km)$lower,
                                risk_percent = (1 - summary(km)$surv)*100, 
                                LB_CI_percent = (1- summary(km)$upper)*100,
                                UB_CI_percent = (1- summary(km)$lower)*100, 
                                nrisk = summary(km)$n.risk, 
                                nevent = summary(km)$n.event, 
                                n_censor = summary(km)$n.censor) 
view(delta_risk_Pf_tot)


ggsurvplot(
  fit = survfit(Surv(Time_all_end_mon, Pf_Event_All) ~ 1, data = Pf_Risk_Tot, stype=2, id=Subject_ID), 
  legend.title = "p. falciparum",
  legend.labs = c(""),
  break.x.by = 2, 
  # Add p-value and tervals
  pval = TRUE,
  conf.int = TRUE,
  #surv.median.line = c("hv"),
  censor = TRUE, 
  # Add risk table
  risk.table = TRUE,
  tables.height = 0.15,
  cumevents = TRUE, 
  cumcensor = TRUE, 
  risk.table.fontsize = 4,
  font.x=4,
  tables.theme = theme_cleantable(),
  xlab = "Months from Baseline", 
  ylab = ("Survival prob."), 
  xlim = c(0,36))


#Summary with 95% CIs for survival estimate
summary(survfit(Surv(Time_all_end_mon, Pf_Event_All) ~ 1, data = Pf_Risk_Tot))





#_________________________________________________________________________________________________#

#######  ACTIVE STUDY POPULATION ONLY 
colors <- c("Pm" = "blue", "Po" = "red", "Pf" = "orange")


p1_Act<- ggplot() + 
  geom_step(data = delta_risk_Pm_act, aes(x = time, y = risk, color = "Pm"), size=1.4)+
  geom_ribbon(data = delta_risk_Pm_act, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "blue", alpha=0.3)+
  geom_step(data = delta_risk_Po_act, aes(x = time, y = risk, color = "Po"), size=1.4)+
  geom_ribbon(data = delta_risk_Po_act, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.3)+
  geom_step(data = delta_risk_Pf_act, aes(x = time, y = risk, color = "Pf"), size=1.4)+
  geom_ribbon(data = delta_risk_Pf_act, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "orange", alpha=0.3)+
  ylab("Cumulative Incidence") +
  xlab("Months from Baseline")+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24))+
  labs(title="Risk of Malaria Infection by Species, Active Visits", 
       x = "Months from Baseline", 
       y = "Cumulative Incidence", 
       color = "Legend")+
  scale_color_manual(values = colors)+
  theme_clean()
#scale_color_manual(name = "Species", 
#          breaks= c("P.malariae", "P.ovale", "P.falciparum"), 
#         values= c("Pm" = "blue", "Po" = "red", "Pf" = "orange"))
plot(p1_Act)



#######  PASSIVE STUDY POPULATION ONLY 
p1_VNP<- ggplot() + 
  geom_step(data = delta_risk_Pm_vnp, aes(x = time, y = risk, color = "Pm"), size=1.4)+
  geom_ribbon(data = delta_risk_Pm_vnp, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "blue", alpha=0.3)+
  geom_step(data = delta_risk_Po_vnp, aes(x = time, y = risk, color = "Po"), size=1.4)+
  geom_ribbon(data = delta_risk_Po_vnp, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.3)+
  geom_step(data = delta_risk_Pf_vnp, aes(x = time, y = risk, color = "Pf"), size=1.4)+
  geom_ribbon(data = delta_risk_Pf_vnp, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "orange", alpha=0.3)+
  ylab("Cumulative Incidence") +
  xlab("Months from Baseline")+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(title="Risk of Malaria Infection by Species, Passive Visits", 
       x = "Months from Baseline", 
       y = "Cumulative Incidence", 
       color = "Species")+
  scale_color_manual(values = colors)+
  theme_clean()
#scale_color_manual(name = "Species", 
#          breaks= c("P.malariae", "P.ovale", "P.falciparum"), 
#         values= c("Pm" = "blue", "Po" = "red", "Pf" = "orange"))
plot(p1_VNP)



####FULL PLOT - TOTAL STUDY POPULATION - VNP and Active - THROUGH 31 DEC 17 - All Species
p1_Tot<-ggplot() + 
  geom_step(data=delta_risk_Po_tot, aes(x = time, y = risk, color = "Po"), size=0.4)+
  geom_ribbon(data=delta_risk_Po_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.2)+
  geom_step(data=delta_risk_Pm_tot, aes(x = time, y = risk, color = "Pm"), size=0.4)+
  geom_ribbon(data=delta_risk_Pm_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "blue", alpha=0.2)+
  geom_step(data=delta_risk_Pf_tot, aes(x = time, y = risk, color = "Pf"), size=0.4)+
  geom_ribbon(data=delta_risk_Pf_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "orange", alpha=0.2)+
  # geom_vline(xintercept = 365, color = "black", lwd = 1, alpha = 0.5)+
  # geom_vline(xintercept = 730, color = "black", lwd = 1, alpha = 0.5)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(title="", 
       x = "Months from Baseline", 
       y = "Cumulative Incidence", 
       color = "Species")+
  scale_color_manual(values = colors)+
  theme_set(
    #theme_clean(base_size = 15))
    theme_clean(base_size = 7))
    

p1_Tot <- p1_Tot + theme(legend.position = c(0.83, 0.40))
plot(p1_Tot)


#####TOTAL STUDY POPULATION: ZOOMED INSERT FOR PLOT WITH MULTIPLE CURVES
p2_Tot<-ggplot() + 
  geom_step(data=delta_risk_Po_tot, aes(x = time, y = risk, color = "Po"), size=0.4)+
  geom_ribbon(data=delta_risk_Po_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.2)+
  geom_step(data=delta_risk_Pm_tot, aes(x = time, y = risk, color = "Pm"), size=0.4)+
  geom_ribbon(data=delta_risk_Pm_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "blue", alpha=0.2)+
  #   geom_vline(xintercept = 365, color = "black", lwd = 1, alpha = 0.5)+
  #   geom_vline(xintercept = 730, color = "black", lwd = 1, alpha = 0.5)+
  ylim(0,0.20)+
  # labs(title="Pm and Po infection incidence, Total Population", fontsize = 4)+
  theme_clean(base_size = 15, panel.grid.major = element_line(size=0.2))+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months", 
       y = "Cum. Incid.")+
  scale_color_manual(values = c("Pm"= "blue", "Po"="red"), guide="none")

plot(p2_Tot)




#######NOW CREATING A ZOOMED INSERT

#p1_Tot_Insert<- 
  p1_Tot+
  annotation_custom(ggplotGrob(p2_Tot), xmin = -1.0, xmax=10.8, ymin=0.70, ymax=1.0)+
  geom_rect(aes(xmin = -1.0, xmax=10.8, ymin=0.70, ymax = 1.0), color = 'black', linetype = "solid", alpha = 0 )



####All species: CI Curves with 95% CIs 
jpeg("../Analysis/Figures/Draft Figures_Manuscript/KMRiskPlot_TotalPop_AllSpecies_Zoomed.jpg", quality = "100", 
     res = 600, bg = "transparent",  height=5, width=6, units = "in" )


### Exporting high res pdf for journal 
ggsave("AllSpecies_TotalPop_CumIncid.pdf", p1_Tot_Insert, width = 88, height = 90, units = "mm") 


### Exporting high res pdf for journal - Plot 1 no inset
ggsave("AllSpecies_TotalPop_CumIncid.pdf", p1_Tot, width = 90, height = 90, units = "mm") 
ggsave("ZoomPmPo_TotalPop_CumIncid.pdf", p2_Tot, width = 90, height = 90, units = "mm") 




##############################
## FIXING FIGURE SIZING FOR NATURE PUBLICATION 


####FULL PLOT - TOTAL STUDY POPULATION - VNP and Active - THROUGH 31 DEC 17 - All Species
p1_Tot<-ggplot() + 
  geom_step(data=delta_risk_Po_tot, aes(x = time, y = risk, color = "Po"), size=0.8)+
  geom_ribbon(data=delta_risk_Po_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.1)+
  geom_step(data=delta_risk_Pm_tot, aes(x = time, y = risk, color = "Pm"), size=0.8)+
  geom_ribbon(data=delta_risk_Pm_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "blue", alpha=0.1)+
  geom_step(data=delta_risk_Pf_tot, aes(x = time, y = risk, color = "Pf"), size=0.8)+
  geom_ribbon(data=delta_risk_Pf_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "orange", alpha=0.1)+
  # geom_vline(xintercept = 365, color = "black", lwd = 1, alpha = 0.5)+
  # geom_vline(xintercept = 730, color = "black", lwd = 1, alpha = 0.5)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(title="", 
       x = "Months from Baseline", 
       y = "Cumulative Incidence", 
       color = "Species")+
  scale_color_manual(values = colors)+
  theme_set(
    #theme_clean(base_size = 15))
    theme_clean(base_size = 7))


p1_Tot <- p1_Tot + theme(legend.position = c(0.84, 0.45))
plot(p1_Tot)


#####TOTAL STUDY POPULATION: ZOOMED INSERT FOR PLOT WITH MULTIPLE CURVES
p2_Tot<-ggplot() + 
  geom_step(data=delta_risk_Po_tot, aes(x = time, y = risk, color = "Po"), size=0.8)+
  geom_ribbon(data=delta_risk_Po_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.1)+
  geom_step(data=delta_risk_Pm_tot, aes(x = time, y = risk, color = "Pm"), size=0.8)+
  geom_ribbon(data=delta_risk_Pm_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "blue", alpha=0.1)+
  #   geom_vline(xintercept = 365, color = "black", lwd = 1, alpha = 0.5)+
  #   geom_vline(xintercept = 730, color = "black", lwd = 1, alpha = 0.5)+
  ylim(0,0.25)+
  # labs(title="Pm and Po infection incidence, Total Population", fontsize = 4)+
  theme_clean()+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months", 
       y = "Cum. Incid.")+
  scale_color_manual(values = c("Pm"= "blue", "Po"="red"), guide="none")

plot(p2_Tot)


#######NOW CREATING A ZOOMED INSERT

p1_Tot+
  annotation_custom(ggplotGrob(p2_Tot), xmin = -1.0, xmax=10.8, ymin=0.70, ymax=1.0)+
  geom_rect(aes(xmin = -1.0, xmax=10.8, ymin=0.70, ymax = 1.0), color = 'black', linetype = "solid", alpha = 0 )


##############################














####FULL PLOT - SURVEY POPULATION -Active ONLY - THROUGH FU3  - All Species
p1_Act<-ggplot() + 
  geom_step(data=delta_risk_Po_act, aes(x = time, y = risk, color = "Po"), size=1.3)+
  geom_ribbon(data=delta_risk_Po_act, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.1)+
  geom_step(data=delta_risk_Pm_act, aes(x = time, y = risk, color = "Pm"), size=1.3)+
  geom_ribbon(data=delta_risk_Pm_act, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "blue", alpha=0.1)+
  geom_step(data=delta_risk_Pf_act, aes(x = time, y = risk, color = "Pf"), size=1.3)+
  geom_ribbon(data=delta_risk_Pf_act, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "orange", alpha=0.1)+
  # geom_vline(xintercept = 365, color = "black", lwd = 1, alpha = 0.5)+
  # geom_vline(xintercept = 730, color = "black", lwd = 1, alpha = 0.5)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(title="", 
       x = "Months from Baseline", 
       y = "Cumulative Incidence", 
       color = "Species")+
  scale_color_manual(values = colors)+
  theme_set(
    theme_classic(base_size = 15))

plot(p1_Act)


#####TOTAL STUDY POPULATION: ZOOMED INSERT FOR PLOT WITH MULTIPLE CURVES
p2_Tot<-ggplot() + 
  geom_step(data=delta_risk_Po_tot, aes(x = time, y = risk, color = "Po"), size=0.8)+
  geom_ribbon(data=delta_risk_Po_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time),  fill = "red", alpha=0.1)+
  geom_step(data=delta_risk_Pm_tot, aes(x = time, y = risk, color = "Pm"), size=0.8)+
  geom_ribbon(data=delta_risk_Pm_tot, aes(ymin=LB_CI, ymax=UB_CI, x=time), fill = "blue", alpha=0.1)+
  #   geom_vline(xintercept = 365, color = "black", lwd = 1, alpha = 0.5)+
  #   geom_vline(xintercept = 730, color = "black", lwd = 1, alpha = 0.5)+
  ylim(0,0.25)+
  # labs(title="Pm and Po infection incidence, Total Population", fontsize = 4)+
  theme_clean()+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months", 
       y = "Cum. Incid.")+
  scale_color_manual(values = c("Pm"= "blue", "Po"="red"), guide="none")

plot(p2_Tot)


#######NOW CREATING A ZOOMED INSERT

p1_Tot+
  annotation_custom(ggplotGrob(p2_Tot), xmin = -1.0, xmax=10.8, ymin=0.70, ymax=1.0)+
  geom_rect(aes(xmin = -1.0, xmax=10.8, ymin=0.70, ymax = 1.0), color = 'black', linetype = "solid", alpha = 0 )



####All species: CI Curves with 95% CIs 
jpeg("../Analysis/Figures/Draft Figures_Manuscript/KMRiskPlot_TotalPop_AllSpecies_Zoomed.jpg", quality = "100", 
     res = 600, bg = "transparent",  height=5, width=6, units = "in" )






############################3
## Comparing the cumulative incidence plot against combined survival plots

Pf_Tot_Surv<-survfit(Surv(Time_all_end_mon, Pf_Event_All) ~ 1, data = Pf_Risk_Tot)
Pm_Tot_Surv<-survfit(Surv(Time_all_end_mon, Pm_Event_All) ~ 1, data = Pm_Risk_Tot)
Po_Tot_Surv<-survfit(Surv(Time_all_end_mon, Po_Event_All) ~ 1, data = Po_Risk_Tot)


fit <- list(pfs=Pf_Tot_Surv, pos=Po_Tot_Surv, pms=Pm_Tot_Surv)
ggsurvplot(fit, combine = TRUE, 
           legend.title = "Survey-based Pop.",
           #  legend.labs = c(""),
           break.x.by = 6, 
           # Add p-value and tervals
           # pval = TRUE,
           conf.int = TRUE,
           #surv.median.line = c("hv"),
           censor = TRUE, 
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.17,
           cumevents = TRUE, 
           # cumcensor = TRUE, 
           #cumcensor.title = "Cumulative number censored",
           cumevent.title = "Cumulative infections",
           #risk.table.fontsize = 2.5,
           # font.x=2.5,
           tables.theme = theme_cleantable(),
           xlab = "Months from Baseline", 
           ylab = ("Prob. un-infected"), 
           xlim = c(0,24))











##########################################################
# CRUDE RISK PLOTS ACROSS STRATA OF PREDICTOR VARS  
##########################################################

#Estimate the Kaplan Meier Crude Risk Function for All Datasets, where event = 1 and Time = time to event or censoring due to LTFU 
##Estimating cumulative incidence = 1-KM survival function (assuming no competing risks and non-informative censoring)
#____________________________________________________________________________________________________________________


##Coding a fake competing event in order to use ggcuminc package.  NA = 2. 
Pf_Risk_Tot <- Pf_Risk_Tot %>% mutate(Pf_Event_All_CE = as.factor(recode(Pf_Event_All_CE, `0` = 0, `1` = 1, `2` = 2)))
Pm_Risk_Tot <- Pm_Risk_Tot %>% mutate(Pm_Event_All_CE = as.factor(recode(Pm_Event_All_CE, `0` = 0, `1` = 1, `2` = 2)))
Po_Risk_Tot <- Po_Risk_Tot %>% mutate(Po_Event_All_CE = as.factor(recode(Po_Event_All_CE, `0` = 0, `1` = 1, `2` = 2)))




#####
#Age Categories (<5, 5-15, >15)  -- var: age_cat2
#####


#Labeling for plots 
Pf_Risk_Tot$age_cat2 <- ordered(Pf_Risk_Tot$age_cat2,
                                levels = c(1,2,3),
                                labels = c("Age <5", "Age 5-14", "Age 15+"))

Pm_Risk_Tot$age_cat2 <- ordered(Pm_Risk_Tot$age_cat2,
                                levels = c(1,2,3),
                                labels = c("Age <5", "Age 5-14", "Age 15+"))
Po_Risk_Tot$age_cat2 <- ordered(Po_Risk_Tot$age_cat2,
                                levels = c(1,2,3),
                                labels = c("Age <5", "Age 5-14", "Age 15+"))



#### Updated with new censoring information 

## P. malariae incidence  by age strata 
Pm_CumInc_AgeStrata <- cuminc(Surv(Time_all_end_mon_2, Pm_Event_All_CE) ~ age_cat2, data = Pm_Risk_Tot) %>% 
  ggcuminc(linetype_aes = TRUE, linewidth = 0.7) + 
  add_confidence_interval() +
  #add_risktable(risktable_group = "strata")+
  add_pvalue(location=c("annotation"), size=1.8)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months from baseline", 
       y = "Cumulative incidence")+
  scale_color_manual(values = c("Age <5"= "#3BBDD1", "Age 5-14"="#3252D1", "Age 15+"="#142052"), guide="none")+
  scale_fill_manual(values = c("Age <5"= "#3BBDD1", "Age 5-14"="#3252D1", "Age 15+"="#142052"))+
  theme(legend.position = c(0.72, 0.47), legend.text = element_text(size = 6), axis.text = element_text(size = 7), axis.title = element_text(size = 7) )

plot(Pm_CumInc_AgeStrata)


## P. ovale incidence  by age strata 
Po_CumInc_AgeStrata <- cuminc(Surv(Time_all_end_mon_2, Po_Event_All_CE) ~ age_cat2, data = Po_Risk_Tot) %>% 
  ggcuminc(linetype_aes = TRUE, linewidth = 0.7) + 
  add_confidence_interval() +
  add_risktable(risktable_group = "strata")+
  add_pvalue(location=c("annotation"), size =1.8)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months from baseline", 
       y = "Cumulative incidence")+
  scale_color_manual(values = c("Age <5"= "#DB5353", "Age 5-14"="#DA2420", "Age 15+"="#6E0D0C"), guide="none")+
  scale_fill_manual(values = c("Age <5"= "#DB5353", "Age 5-14"="#DA2420", "Age 15+"="#6E0D0C"))+
  theme(legend.position = c(0.72, 0.47), legend.text = element_text(size = 6), axis.text = element_text(size = 7), axis.title = element_text(size = 7) )

plot(Po_CumInc_AgeStrata)


## P. falciparum incidence  by age strata 
Pf_CumInc_AgeStrata <- cuminc(Surv(Time_all_end_mon_2, Pf_Event_All_CE) ~ age_cat2, data = Pf_Risk_Tot) %>% 
  ggcuminc(linetype_aes = TRUE, linewidth = 0.7) + 
  add_confidence_interval() +
  add_risktable(risktable_group = "strata")+
  add_pvalue(location=c("annotation"), size = 1.8)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months from baseline", 
       y = "Cumulative incidence")+
  scale_color_manual(values = c("Age <5"= "#f0d10c", "Age 5-14"="#f76b07", "Age 15+"="#874e04"), guide="none")+
  scale_fill_manual(values = c("Age <5"= "#f0d10c", "Age 5-14"="#f76b07", "Age 15+"="#874e04"))+
  theme(legend.position = c(0.72, 0.40), legend.text = element_text(size = 6), axis.text = element_text(size = 7), axis.title = element_text(size = 7) )

plot(Pf_CumInc_AgeStrata)







###Combining all species plots together for 3-panel figure 4A
Aggregate_AgeStrat_CI_Plot<-cowplot::plot_grid(Pm_CumInc_AgeStrata, Po_CumInc_AgeStrata, Pf_CumInc_AgeStrata, ncol=3)
plot(Aggregate_AgeStrat_CI_Plot)


#### Exporting High Res Figures for Publication 

ggsave("Aggregate_AgeStrat_CI_6Sep23.pdf", Aggregate_AgeStrat_CI_Plot, width = 180, height = 90, units = "mm") 




#####
#Sex Categories (1=Female; 0=Male)  -- var: sex
#####

#Labeling for plots 
Pf_Risk_Tot$gender <- ordered(Pf_Risk_Tot$gender,
                              levels = c(0,1),
                              labels = c("Male", "Female"))
Pm_Risk_Tot$gender <- ordered(Pm_Risk_Tot$gender,
                              levels = c(0,1),
                              labels = c("Male", "Female"))
Po_Risk_Tot$gender <- ordered(Po_Risk_Tot$gender,
                              levels = c(0,1),
                              labels = c("Male", "Female"))


## P. malariae incidence  by sex strata 
Pm_CumInc_SexStrata <- cuminc(Surv(Time_all_end_mon_2, Pm_Event_All_CE) ~ gender, data = Pm_Risk_Tot) %>% 
  ggcuminc(linetype_aes = TRUE, linewidth = 0.7) + 
  add_confidence_interval() +
  add_risktable(risktable_group = "strata")+
  add_pvalue(location=c("annotation"), size=1.8)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months from baseline", 
       y = "Cumulative incidence")+
  scale_color_manual(values = c("Female"= "#3BBDD1", "Male"="#3252D1"), guide="none")+
  scale_fill_manual(values = c("Female"= "#3BBDD1", "Male"="#3252D1"))+
  theme(legend.position = c(0.72, 0.47), legend.text = element_text(size = 6), axis.text = element_text(size = 7), axis.title = element_text(size = 7) )

plot(Pm_CumInc_SexStrata)


## P. ovale incidence  by age strata 
Po_CumInc_SexStrata <- cuminc(Surv(Time_all_end_mon_2, Po_Event_All_CE) ~ gender, data = Po_Risk_Tot) %>% 
  ggcuminc(linetype_aes = TRUE, linewidth = 0.7) + 
  add_confidence_interval() +
  add_risktable(risktable_group = "strata")+
  add_pvalue(location=c("annotation"), size=1.8)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months from baseline", 
       y = "Cumulative incidence")+
  scale_color_manual(values = c("Female"= "#DB5353", "Male"="#6E0D0C"), guide="none")+
  scale_fill_manual(values = c("Female"= "#DB5353", "Male"="#6E0D0C"))+
  theme(legend.position = c(0.72, 0.47), legend.text = element_text(size = 6), axis.text = element_text(size = 7), axis.title = element_text(size = 7) )

plot(Po_CumInc_SexStrata)


## P. falciparum incidence  by age strata 
Pf_CumInc_SexStrata <- cuminc(Surv(Time_all_end_mon_2, Pf_Event_All_CE) ~ gender, data = Pf_Risk_Tot) %>% 
  ggcuminc(linetype_aes = TRUE, linewidth = 0.7) + 
  add_confidence_interval() +
  add_risktable(risktable_group = "strata")+
  add_pvalue(location=c("annotation"), size=1.8)+
  ylim(0,1.0)+
  scale_x_continuous(limits = c(0, 36), breaks = c(0, 6, 12, 18, 24, 30, 36))+
  labs(x = "Months from baseline", 
       y = "Cumulative incidence")+
  scale_color_manual(values = c("Female"= "#ffc917", "Male"="#ff7817"), guide="none")+
  scale_fill_manual(values = c("Female"= "#ffc917", "Male"="#ff7817"))+
  theme(legend.position = c(0.72, 0.47), legend.text = element_text(size = 6), axis.text = element_text(size = 7), axis.title = element_text(size = 7) )

plot(Pf_CumInc_SexStrata)




###Combining all species plots together for 3-panel figure 4A
Aggregate_SexStrat_CI_Plot<-cowplot::plot_grid(Pm_CumInc_SexStrata, Po_CumInc_SexStrata, Pf_CumInc_SexStrata, ncol=3)
plot(Aggregate_SexStrat_CI_Plot)


#### Exporting High Res Figures for Publication 

ggsave("Aggregate_SexStrat_CI_6Sep23.pdf", Aggregate_SexStrat_CI_Plot, width = 180, height = 90, units = "mm") 

