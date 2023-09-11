##########################################
#Programmer: `Rachel Sendor
#Last update: Mar 2023
#Purpose:     DRC Non-falciparum Descriptive Epidemiology Analysis
#             
##########################################

##########################################################
#SETTING PARAMETERS / LOADING PACKAGES FOR USE

options(max.print=999999)
options(scipen = 999)
.libPaths()


 #Loading R Packages:
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(sf)
library(tableone)
library(devtools)
library(PropCIs)
library(ggplot2)
library(ggbreak)
library(ggthemes)
library(resample)
library(ggpubr)
library(haven)
library(lubridate)
library(plyr)
library(survminer)
library(purrr)
library(stringr)
library(forcats)
library(broom)
library(gee)
library(geepack)
library(ggExtra)
library(janitor)
library(linelist)
library(lme4)
library(lmerTest)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

#help(package="tidycmprsk")


##########################################################
# SET WORKING DIRECTORY

#### blinded file path
setwd(".../DRC Non-falciparum Project/Analysis/Dataset_AD")



# IMPORTING DATASETS FOR ANALYSIS --- DATASET CREATED IN SAS FOR DATA LINKAGE, CLEANING, & PRELIMINARY WRANGLING (SEE PROGRAM: AnalysisDatasetBuild_DRCNonPf_Mar22_RS.sas)


AD_Total<- read_excel('.../AD_Total_Jul22.xlsx', NULL)
AD_Active<- read_excel('.../AD_Active_Jul22.xlsx', NULL)
AD_Passive<- read_excel('.../AD_Passive_Jul22.xlsx',NULL)



#Active Population 
##Look by visit across same species, and between species
##exclude missing visits 
AD_Active_nonmissvisits<-AD_Active%>%filter(!is.na(Visit_date))

      
      
  
      
      
##########################################################
# ASSESSING FACTORS ASSOCIATED WITH PREVALENCE OF PM, PO, AND PF INFECTION    
##########################################################

## Method:  Modeling as Binomial GEE with Identity Link (for prev. diff), accounting for correlations between repeated testing across longitudinal follow-up visits within subjects. 
#           Using robust standard errors. 
#           Used Exch. correlation matrix

## 
## Relevant Sources: 1)  https://data.library.virginia.edu/getting-started-with-generalized-estimating-equations/


#### Total Population PDs ####

# N=1,595 subjects / 5,689 visits (BL, FU1, Fu2, Fu3) 
# Dataset: AD_Total

#Pm Infection (Any):  N=187 infections across 5688 non-missing Total visits   (572 FU visits missing as expected)
#Po Infection (Any):  N=78 infections across 5688 non-missing Total visits   (572 FU visits missing as expected)
#Pf Infection (Any):  N=1976 infections across 5688 non-missing Total visits   (572 FU visits missing as expected)


##Soring data by Subject_ID and then by Visit number to ensure it's formatted correctly for the clustering analysis. 
AD_Active <- AD_Active[order(AD_Active$Subject_ID, AD_Active$Visit, AD_Active$Visit_date),]


##Crude prevalence (accounting for repeated testing;  not accounting for household clustering)

##Pm

## Number of samples included in prevalence difference estimates (ie non-missing Pm)
AD_Active_Pmnomiss <- AD_Active%>% filter(Pm ==1 | Pm==0)

addmargins(table(AD_Active$Pm, useNA="always"))
addmargins(table(AD_Active_Pmnomiss$Pm, useNA="always"))

##PD = Prevalence Difference 
CrudePrev_GEE_Pm_Active_Pd<-geeglm(Pm ~ 1, 
                                   data = AD_Active, 
                                   id = Subject_ID, 
                                   family = binomial(link="identity"),
                                   corstr = "ex")

summary(CrudePrev_GEE_Pm_Active_Pd)
confint(fit) # 95% CI for the coefficients

##PR = Prevalence Ratio
CrudePrev_GEE_Pm_Active_PR<-geeglm(Pm ~ 1, 
                                   data = AD_Active, 
                                   id = Subject_ID, 
                                   family = binomial(link="log"),
                                   corstr = "ex")

summary(CrudePrev_GEE_Pm_Active_PR)


### Checking Mixed model accounting for within-subject and household clustering
CrudePrev_ME_Pm <- glmer(Pm ~ 1 + (1 | Village/Household/Subject_ID),
                        data = AD_Active,
                        family = binomial(link="log"))
summary(CrudePrev_ME_Pm)

CrudePrev_ME_Pm <- glmer(Pm ~ 1 + (Subject_ID | Village/Household/Subject_ID),
                        data = AD_Active,
                        family = binomial(link="log"))
summary(CrudePrev_ME_Pm)

CrudePrev_ME_Pm <- lmer(Pm ~ 1 + (Subject_ID | Household) + (Household | Village)  ,
                        data = AD_Active,
                        family = binomial(link="log"),
                        REML = FALSE)
summary(CrudePrev_ME_Pm)



##Po

AD_Active_Ponomiss <- AD_Active%>% filter(Po ==1 | Po==0)

addmargins(table(AD_Active$Po, useNA="always"))
addmargins(table(AD_Active_Ponomiss$Po, useNA="always"))



CrudePrev_GEE_Po<-geeglm(Po ~ 1, 
                         data = AD_Active, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(CrudePrev_GEE_Po)
confint(fit) # 95% CI for the coefficients



### Mixed model accounting for within-subject and household clustering
CrudePrev_ME_Po <- lmer(Po ~ 1 + (1 | Subject_ID) + (1 | Household),
                        data = AD_Active,
                        family = binomial(link="logit"),
                        REML = FALSE)
summary(CrudePrev_ME_Po)



##Pf

AD_Active_Pfnomiss <- AD_Active%>% filter(Pf ==1 | Pf==0)

addmargins(table(AD_Active$Pf, useNA="always"))
addmargins(table(AD_Active_Pfnomiss$Pf, useNA="always"))




CrudePrev_GEE_Pf<-geeglm(Pf ~ 1, 
                         data = AD_Active, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(CrudePrev_GEE_Pf)
confint(fit) # 95% CI for the coefficients


### Mixed model accounting forwithin-subject and household clutering
CrudePrev_ME_Pf <- lmer(Pf ~ 1 + (1 | Subject_ID) + (1 | Household),
                        data = AD_Active,
                        REML = FALSE)
summary(CrudePrev_ME_Pf)






##
##AGE CATEGORY (Age >=15 years = Ref)
##
AD_Active_AgeNoMiss<-AD_Active%>%filter(!is.na(AgeCat_Visit))
addmargins(table(AD_Active_AgeNoMiss$AgeCat_Visit, AD_Active_AgeNoMiss$Pm))
addmargins(table(AD_Active_AgeNoMiss$AgeCat_Visit, AD_Active_AgeNoMiss$Po))


#Pm#
Age_GEE_Pm_PD<-geeglm(Pm ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Active_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Pm_PD)

Age_GEE_Pm_PR<-geeglm(Pm ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Active_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Pm_PR)


#Age_GEE_Pm<-geeglm(Pm ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Active_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Pm)


#Po#

Age_GEE_Po_PD<-geeglm(Po ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Active_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Po_PD)

Age_GEE_Po_PR<-geeglm(Po ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Active_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Po_PR)


#Age_GEE_Po<-geeglm(Po ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Active_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Po)


#Pf#
Age_GEE_Pf_PD<-geeglm(Pf ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Active_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Pf_PD)

Age_GEE_Pf_PR<-geeglm(Pf ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Active_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Pf_PR)


#Age_GEE_Pf<-geeglm(Pf ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Active_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Pf)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Active_Age_Pm_PD<-tidy(Age_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Age_Po_PD<-tidy(Age_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Age_Pf_PD<-tidy(Age_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Active_Age_Pm_PR<-tidy(Age_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Age_Po_PR<-tidy(Age_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Age_Pf_PR<-tidy(Age_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##SEX CATEGORY (Male = Ref)
##
AD_Active_SexNoMiss<-AD_Active%>%filter(!is.na(Sex))

#Pm#
Sex_GEE_Pm_PD<-geeglm(Pm ~ Sex, 
                      data = AD_Active_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Pm_PD)

Sex_GEE_Pm_PR<-geeglm(Pm ~ Sex, 
                      data = AD_Active_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Pm_PR)

#Po#
Sex_GEE_Po_PD<-geeglm(Po ~ Sex, 
                      data = AD_Active_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Po_PD)

Sex_GEE_Po_PR<-geeglm(Po ~ Sex, 
                      data = AD_Active_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Po_PR)

#Pf#
Sex_GEE_Pf_PD<-geeglm(Pf ~ Sex, 
                      data = AD_Active_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Pf_PD)

Sex_GEE_Pf_PR<-geeglm(Pf ~ Sex, 
                      data = AD_Active_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Pf_PR)



##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator

#Prevalence Differences 
Coeff_Active_Sex_Pm_PD<-tidy(Sex_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Sex_Po_PD<-tidy(Sex_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Sex_Pf_PD<-tidy(Sex_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios
Coeff_Active_Sex_Pm_PR<-tidy(Sex_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Sex_Po_PR<-tidy(Sex_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Sex_Pf_PR<-tidy(Sex_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




##
##HEALTH AREA CATEGORY (1 = LINGWALA = REF CAT (urban);  Bu = 2 = Ref ; (rural);  3 = KIMPOKO  (peri-urban))
##
AD_Active_HealthAreaNoMiss<-AD_Active%>%filter(!is.na(HealthArea))

#Pm#
HealthArea_GEE_Pm_PD<-geeglm(Pm ~ as.factor(HealthArea), 
                             data = AD_Active_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Pm_PD)

HealthArea_GEE_Pm_PR<-geeglm(Pm ~ as.factor(HealthArea), 
                             data = AD_Active_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Pm_PR)

#Po#
HealthArea_GEE_Po_PD<-geeglm(Po ~ as.factor(HealthArea),  
                             data = AD_Active_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Po_PD)

HealthArea_GEE_Po_PR<-geeglm(Po ~ as.factor(HealthArea),  
                             data = AD_Active_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Po_PR)

#Pf#
HealthArea_GEE_Pf_PD<-geeglm(Pf ~ as.factor(HealthArea),  
                             data = AD_Active_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Pf_PD)

HealthArea_GEE_Pf_PR<-geeglm(Pf ~ as.factor(HealthArea),  
                             data = AD_Active_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Pf_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Active_HealthArea_Pm_PD<-tidy(HealthArea_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_HealthArea_Po_PD<-tidy(HealthArea_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_HealthArea_Pf_PD<-tidy(HealthArea_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios
Coeff_Active_HealthArea_Pm_PR<-tidy(HealthArea_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_HealthArea_Po_PR<-tidy(HealthArea_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_HealthArea_Pf_PR<-tidy(HealthArea_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")





##
##WEALTH CATEGORY (Average Wealth [3] = Ref)
##

AD_Active_WealthCatNoMiss<-AD_Active_Pfnomiss%>%filter(!is.na(WealthCat))
## comparing poor vs. average, and rich vs. average  (average = 1; poor = 2; rich = 3)

#Pm#
Wealth_GEE_Pm_PD<-geeglm(Pm ~ as.factor(WealthCat), 
                         data = AD_Active_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Pm_PD)

Wealth_GEE_Pm_PR<-geeglm(Pm ~ as.factor(WealthCat), 
                         data = AD_Active_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Pm_PR)


#Po#
Wealth_GEE_Po_PD<-geeglm(Po ~ as.factor(WealthCat), 
                         data = AD_Active_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Po_PD)

Wealth_GEE_Po_PR<-geeglm(Po ~ as.factor(WealthCat), 
                         data = AD_Active_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Po_PR)

#Pf#
Wealth_GEE_Pf_PD<-geeglm(Pf ~ as.factor(WealthCat), 
                         data = AD_Active_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Pf_PD)

Wealth_GEE_Pf_PR<-geeglm(Pf ~ as.factor(WealthCat), 
                         data = AD_Active_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Pf_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Difference 
Coeff_Active_Wealth_Pm_PD<-tidy(Wealth_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Wealth_Po_PD<-tidy(Wealth_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Wealth_Pf_PD<-tidy(Wealth_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Active_Wealth_Pm_PR<-tidy(Wealth_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Wealth_Po_PR<-tidy(Wealth_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Wealth_Pf_PR<-tidy(Wealth_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##FEVER CATEGORY (No Fever = Ref)
##
AD_Active_FeverNoMiss<-AD_Active%>%filter(!is.na(Feverpos))

addmargins(table(AD_Active_FeverNoMiss$Feverpos, AD_Active_FeverNoMiss$Pm, useNA = "always"))
addmargins(table(AD_Active_FeverNoMiss$Feverpos, AD_Active_FeverNoMiss$Po, useNA = "always"))

#Pm#
Fever_GEE_Pm_PD<-geeglm(Pm ~ Feverpos, 
                        data = AD_Active_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Pm_PD)

Fever_GEE_Pm_PR<-geeglm(Pm ~ Feverpos, 
                        data = AD_Active_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Pm_PR)


#Po#
Fever_GEE_Po_PD<-geeglm(Po ~ Feverpos, 
                        data = AD_Active_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Po_PD)

Fever_GEE_Po_PR<-geeglm(Po ~ Feverpos, 
                        data = AD_Active_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Po_PR)

#Pf#
Fever_GEE_Pf_PD<-geeglm(Pf ~ Feverpos, 
                        data = AD_Active_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Pf_PD)

Fever_GEE_Pf_PR<-geeglm(Pf ~ Feverpos, 
                        data = AD_Active_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Active_Fever_Pm_PD<-tidy(Fever_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Fever_Po_PD<-tidy(Fever_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Fever_Pf_PD<-tidy(Fever_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Active_Fever_Pm_PR<-tidy(Fever_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Fever_Po_PR<-tidy(Fever_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_Fever_Pf_PR<-tidy(Fever_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




write.table(Coeff_Active_Fever_Pm_PD, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Active_Fever_Po_PD, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Active_Fever_Pf_PD, "clipboard", sep="\t", col.names = F )

write.table(Coeff_Active_Fever_Pm_PR, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Active_Fever_Po_PR, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Active_Fever_Pf_PR, "clipboard", sep="\t", col.names = F )



##
##RDT-positive CATEGORY (Not RDT Positive = Ref)
##
AD_Active_RDTNoMiss<-AD_Active%>%filter(!is.na(RDTpos))

#Pm#
RDT_GEE_Pm_PD<-geeglm(Pm ~ RDTpos, 
                      data = AD_Active_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Pm_PD)

RDT_GEE_Pm_PR<-geeglm(Pm ~ RDTpos, 
                      data = AD_Active_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Pm_PR)


#Po#
RDT_GEE_Po_PD<-geeglm(Po ~ RDTpos, 
                      data = AD_Active_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Po_PD)

RDT_GEE_Po_PR<-geeglm(Po ~ RDTpos, 
                      data = AD_Active_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Po_PR)

#Pf#
RDT_GEE_Pf_PD<-geeglm(Pf ~ RDTpos, 
                      data = AD_Active_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Pf_PD)

RDT_GEE_Pf_PR<-geeglm(Pf ~ RDTpos, 
                      data = AD_Active_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Active_RDT_Pm_PD<-tidy(RDT_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_RDT_Po_PD<-tidy(RDT_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_RDT_Pf_PD<-tidy(RDT_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Active_RDT_Pm_PR<-tidy(RDT_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_RDT_Po_PR<-tidy(RDT_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_RDT_Pf_PR<-tidy(RDT_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##P.falciparum Co-infection (Neg Pf = Ref)
##
AD_Active_PfNoMiss<-AD_Active%>%filter(!is.na(Pf))
#Pm#
Pf_GEE_Pm_PD<-geeglm(Pm ~ Pf, 
                     data = AD_Active_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="identity"),
                     corstr = "ex")
summary(Pf_GEE_Pm_PD)

Pf_GEE_Pm_PR<-geeglm(Pm ~ Pf, 
                     data = AD_Active_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="log"),
                     corstr = "ex")
summary(Pf_GEE_Pm_PR)

#Po#
Pf_GEE_Po_PD<-geeglm(Po ~ Pf, 
                     data = AD_Active_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="identity"),
                     corstr = "ex")
summary(Pf_GEE_Po_PD)

Pf_GEE_Po_PR<-geeglm(Po ~ Pf, 
                     data = AD_Active_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="log"),
                     corstr = "ex")
summary(Pf_GEE_Po_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence DIfferences 
Coeff_Active_Pfcoinfec_Pm_PD<-tidy(Pf_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Pfcoinfec_Po_PD<-tidy(Pf_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Active_Pfcoinfec_Pm_PR<-tidy(Pf_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_Pfcoinfec_Po_PR<-tidy(Pf_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")



##
##SEASONALITY CATEGORY (Rainy Season = Ref)
##
AD_Active_SeasonNoMiss<-AD_Active%>%filter(!is.na(Season_dry))
#Pm#
DrySeason_GEE_Pm_PD<-geeglm(Pm ~ Season_dry, 
                            data = AD_Active_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Pm_PD)

DrySeason_GEE_Pm_PR<-geeglm(Pm ~ Season_dry, 
                            data = AD_Active_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Pm_PR)

#Po#
DrySeason_GEE_Po_PD<-geeglm(Po ~ Season_dry, 
                            data = AD_Active_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Po_PD)

DrySeason_GEE_Po_PR<-geeglm(Po ~ Season_dry, 
                            data = AD_Active_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Po_PR)

#Pf#
DrySeason_GEE_Pf_PD<-geeglm(Pf ~ Season_dry, 
                            data = AD_Active_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Pf_PD)

DrySeason_GEE_Pf_PR<-geeglm(Pf ~ Season_dry, 
                            data = AD_Active_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalecne Differneces 
Coeff_Active_DrySeason_Pm_PD<-tidy(DrySeason_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_DrySeason_Po_PD<-tidy(DrySeason_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_DrySeason_Pf_PD<-tidy(DrySeason_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalecne Ratios 
Coeff_Active_DrySeason_Pm_PR<-tidy(DrySeason_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_DrySeason_Po_PR<-tidy(DrySeason_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_DrySeason_Pf_PR<-tidy(DrySeason_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##BED NET USE PRIOR NIGHT (No Use = Ref)
##
AD_Active_BedNetNoMiss<-AD_Active%>%filter(!is.na(BedNetPrior))
#Pm#
BedNet_GEE_Pm_PD<-geeglm(Pm ~ BedNetPrior, 
                         data = AD_Active_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(BedNet_GEE_Pm_PD)

BedNet_GEE_Pm_PR<-geeglm(Pm ~ BedNetPrior, 
                         data = AD_Active_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(BedNet_GEE_Pm_PR)

#Po#
BedNet_GEE_Po_PD<-geeglm(Po ~ BedNetPrior, 
                         data = AD_Active_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(BedNet_GEE_Po_PD)

BedNet_GEE_Po_PR<-geeglm(Po ~ BedNetPrior, 
                         data = AD_Active_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(BedNet_GEE_Po_PR)
#Pf#
BedNet_GEE_Pf_PD<-geeglm(Pf ~ BedNetPrior, 
                         data = AD_Active_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(BedNet_GEE_Pf_PD)

BedNet_GEE_Pf_PR<-geeglm(Pf ~ BedNetPrior, 
                         data = AD_Active_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(BedNet_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Active_BedNet_Pm_PD<-tidy(BedNet_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_BedNet_Po_PD<-tidy(BedNet_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_BedNet_Pf_PD<-tidy(BedNet_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Active_BedNet_Pm_PR<-tidy(BedNet_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Active_BedNet_Po_PR<-tidy(BedNet_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Active_BedNet_Pf_PR<-tidy(BedNet_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")






######
#Active Visits:  Stacking all GEE Output for Estimate (95% CI) Plots

#Prevalence Differences

Active_GEE_Output_PD<-rbind(
  Coeff_Active_Age_Pm_PD, Coeff_Active_Age_Po_PD, Coeff_Active_Age_Pf_PD,
  Coeff_Active_Sex_Pm_PD, Coeff_Active_Sex_Po_PD, Coeff_Active_Sex_Pf_PD,
  Coeff_Active_HealthArea_Pm_PD, Coeff_Active_HealthArea_Po_PD, Coeff_Active_HealthArea_Pf_PD,
  #Coeff_Active_Village_Pm_PD, Coeff_Active_Village_Po_PD, Coeff_Active_Village_Pf_PD,
  Coeff_Active_Wealth_Pm_PD, Coeff_Active_Wealth_Po_PD, Coeff_Active_Wealth_Pf_PD,
  Coeff_Active_Fever_Pm_PD, Coeff_Active_Fever_Po_PD, Coeff_Active_Fever_Pf_PD,
  Coeff_Active_RDT_Pm_PD, Coeff_Active_RDT_Po_PD, Coeff_Active_RDT_Pf_PD,
  Coeff_Active_Pfcoinfec_Pm_PD, Coeff_Active_Pfcoinfec_Po_PD,
  Coeff_Active_DrySeason_Pm_PD, Coeff_Active_DrySeason_Po_PD, Coeff_Active_DrySeason_Pf_PD,
  Coeff_Active_BedNet_Pm_PD, Coeff_Active_BedNet_Po_PD, Coeff_Active_BedNet_Pf_PD
  # Coeff_Active_Visit_Pm_PD, Coeff_Active_Visit_Po_PD, Coeff_Active_Visit_Pf_PD
)

write.table(Active_GEE_Output_PD, "clipboard", sep="\t", col.names = T )



#Prevalence Ratios

Active_GEE_Output_PR<-rbind(
  Coeff_Active_Age_Pm_PR, Coeff_Active_Age_Po_PR, Coeff_Active_Age_Pf_PR,
  Coeff_Active_Sex_Pm_PR, Coeff_Active_Sex_Po_PR, Coeff_Active_Sex_Pf_PR,
  Coeff_Active_HealthArea_Pm_PR, Coeff_Active_HealthArea_Po_PR, Coeff_Active_HealthArea_Pf_PR,
  #Coeff_Active_Village_Pm_PD, Coeff_Active_Village_Po_PD, Coeff_Active_Village_Pf_PD,
  Coeff_Active_Wealth_Pm_PR, Coeff_Active_Wealth_Po_PR, Coeff_Active_Wealth_Pf_PR,
  Coeff_Active_Fever_Pm_PR, Coeff_Active_Fever_Po_PR, Coeff_Active_Fever_Pf_PR,
  Coeff_Active_RDT_Pm_PR, Coeff_Active_RDT_Po_PR, Coeff_Active_RDT_Pf_PR,
  Coeff_Active_Pfcoinfec_Pm_PR, Coeff_Active_Pfcoinfec_Po_PR,
  Coeff_Active_DrySeason_Pm_PR, Coeff_Active_DrySeason_Po_PR, Coeff_Active_DrySeason_Pf_PR,
  Coeff_Active_BedNet_Pm_PR, Coeff_Active_BedNet_Po_PR, Coeff_Active_BedNet_Pf_PR
  # Coeff_Active_Visit_Pm_PD, Coeff_Active_Visit_Po_PD, Coeff_Active_Visit_Pf_PD
)

write.table(Active_GEE_Output_PR, "clipboard", sep="\t", col.names = T )




#############
##  Plotting Prevalence Differences and 95% CIs by Population

Active_GEE_Output_PD$Ifxn <- factor(Active_GEE_Output_PD$Ifxn, levels=c("Pm", "Po", "Pf"))
Active_GEE_Output_PD$Ifxn <- factor(Active_GEE_Output_PD$Ifxn, levels=c("Pm", "Po", "Pf"))

Active_GEE_Output_PD$term <- recode_factor(Active_GEE_Output_PD$term, 
                                           "relevel(as.factor(AgeCat_Visit), ref = \"3\")1" = "Age <5 vs 15+", 
                                           "relevel(as.factor(AgeCat_Visit), ref = \"3\")2" = "Age 5-14 vs 15+",
                                           "Sex" = "Sex (F vs M)",
                                           "as.factor(HealthArea)2" = "Rural vs. Urban", 
                                           "as.factor(HealthArea)3" = "Peri-urban vs. Urban", 
                                           "as.factor(WealthCat)2" = "Poor vs Avg. Wealth",
                                           "as.factor(WealthCat)3" = "Wealthy vs Avg. Wealth",
                                           "Season_dry" =  "Dry vs Rainy Season",
                                           "RDTpos" = "RDT+ vs RDT-",
                                           "Pf" = "Pf coinfection (Y vs N)",
                                           "Feverpos" = "Fever (Y vs N)",
                                           "BedNetPrior" = "Bed Net Use (Y vs N)")

Active_GEE_Output_PR$Ifxn <- factor(Active_GEE_Output_PR$Ifxn, levels=c("Pm", "Po", "Pf"))
Active_GEE_Output_PR$Ifxn <- factor(Active_GEE_Output_PR$Ifxn, levels=c("Pm", "Po", "Pf"))

Active_GEE_Output_PR$term <- recode_factor(Active_GEE_Output_PR$term, 
                                           "relevel(as.factor(AgeCat_Visit), ref = \"3\")1" = "Age <5 vs 15+", 
                                           "relevel(as.factor(AgeCat_Visit), ref = \"3\")2" = "Age 5-14 vs 15+",
                                           "Sex" = "Sex (F vs M)",
                                           "as.factor(HealthArea)2" = "Rural vs. Urban", 
                                           "as.factor(HealthArea)3" = "Peri-urban vs. Urban", 
                                           "as.factor(WealthCat)2" = "Poor vs Avg. Wealth",
                                           "as.factor(WealthCat)3" = "Wealthy vs Avg. Wealth",
                                           "Season_dry" =  "Dry vs Rainy Season",
                                           "RDTpos" = "RDT+ vs RDT-",
                                           "Pf" = "Pf coinfection (Y vs N)",
                                           "Feverpos" = "Fever (Y vs N)",
                                           "BedNetPrior" = "Bed Net Use (Y vs N)")



Active_GEE_Output$term <- recode_factor(Active_GEE_Output$term, 
                                        "agecat_LT5" = "Age <5 vs 15+", 
                                        "agecat_5_14" = "Age 5-14 vs 15+",
                                        "Sex" = "Sex (F vs M)",
                                        "as.factor(HealthArea)2" = "Rural vs. Urban", 
                                        "as.factor(HealthArea)3" = "Peri-urban vs. Urban", 
                                        "Village_Impuru" = "Site-Impuru vs Bu",
                                        "Village_Pema" = "Site-Pema vs Bu",
                                        "Village_Kimpoko" = "Site-Kimpoko vs Bu",
                                        "Village_Ngamanzo" = "Site-Ngamanzo vs Bu",
                                        "Village_Iye" =  "Site-Iye vs Bu",
                                        "Village_Lingwala" = "Site-Voix de Peuple vs Bu",
                                        "Wealth_Poorest" = "Poorest vs Avg. Wealth",
                                        "Wealth_Poorer" = "Poorer vs Avg. Wealth",
                                        "Wealth_Wealthier" = "Wealthier vs Avg. Wealth",
                                        "Wealth_Wealthiest" = "Wealthiest vs Avg. Wealth",
                                        "Season_dry" =  "Dry vs Rainy Season",
                                        "RDTpos" = "RDT+ vs RDT-",
                                        "Pf" = "Pf coinfection (Y vs N)",
                                        "Feverpos" = "Fever (Y vs N)",
                                        "BedNetPrior" = "Bed Net Use (Y vs N)",
                                        "Visit" = "Calendar Time")


##checking the ordering of variables for gggplotting: 
Active_GEE_Output_PD$term %>%levels()
Active_GEE_Output_PR$term %>%levels()
#Passive_GEE_Output_PD$term %>%levels()
Active_GEE_Output_PD$term %>%levels()


#set order of levels to appear along y-axis
Active_GEE_Output_PD <- Active_GEE_Output_PD%>%
  mutate(term = term %>% 
           fct_relevel("Wealthy vs Avg. Wealth",
                       "Poor vs Avg. Wealth",
                       "Dry vs Rainy Season",
                       "Bed Net Use (Y vs N)",
                       "Fever (Y vs N)",
                       "Pf coinfection (Y vs N)",
                       "RDT+ vs RDT-",
                       "Peri-urban vs. Urban",
                       "Rural vs. Urban",
                       "Sex (F vs M)",
                       "Age 5-14 vs 15+",
                       "Age <5 vs 15+", 
           ))

Active_GEE_Output_PR <- Active_GEE_Output_PR%>%
  mutate(term = term %>% 
           fct_relevel("Wealthy vs Avg. Wealth",
                       "Poor vs Avg. Wealth",
                       "Dry vs Rainy Season",
                       "Bed Net Use (Y vs N)",
                       "Fever (Y vs N)",
                       "Pf coinfection (Y vs N)",
                       "RDT+ vs RDT-",
                       "Peri-urban vs. Urban",
                       "Rural vs. Urban",
                       "Sex (F vs M)",
                       "Age 5-14 vs 15+",
                       "Age <5 vs 15+", 
           ))



ActivePop_PD<-ggplot(Active_GEE_Output_PD, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.5,  width=0.4, position=position_dodge(width = 0.7))+
  geom_point(size=2, shape=21, colour="black", stroke = 0.6, position=position_dodge(width = 0.7)) +
  geom_vline(xintercept=0.0, lty=2, lwt=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Difference")+
  ylab("")+
  #labs(title = "Active Pop.",
  #     subtitle = "n=1,565 participants across 5,682 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=5, color="black"), plot.subtitle = element_text(size=5), axis.text.x=element_text(size=5, color="black"), axis.text.y.left=element_text(size=5,  color="black"), legend.text=element_text(size=5, color="black"), legend.title = element_text(size=5, color="black")
        #        , panel.background = element_rect(fill = "white", color = "black")
        #       , panel.grid.major.x = element_line(size = 0.3, linetype = 'solid', colour = "#D3D3D3"), panel.grid.minor.x = element_line(size = 0.2, linetype = 'solid', colour = "#D3D3D3"), panel.grid.major.y = element_blank()
        , legend.position = c(0.83, 0.88), legend.background = element_rect(size=0.3, linetype="solid", colour ="black")
        #          , legend.position = "none"
  )

#ActivePop_PD_noPf_2<-ActivePop_PD_noPf + expand_limits(x=c(-0.08, 0.08))
ActivePop_PD_2<-ActivePop_PD + xlim(-0.30, 0.75)
plot(ActivePop_PD_2) 



##Now without Pf also for a zoomed version on just Pm and Po differences
Active_GEE_Output_PD_NoPf<-Active_GEE_Output_PD%>%filter(Ifxn!= "Pf")
ActivePop_PD_noPf<-ggplot(Active_GEE_Output_PD_NoPf, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.4,  width=0.4, position=position_dodge(width = 0.3))+
  geom_point(size=2, shape=21, colour="black", stroke = 0.6, position=position_dodge(width = 0.3)) +
  geom_vline(xintercept=0.0, lty=2, lwt=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Difference")+
  ylab("")+
  #labs(title = "Active Pop.",
  #     subtitle = "n=1,565 participants across 5,682 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title = element_text(size=5, color="black"), plot.subtitle = element_text(size=5), axis.text.x=element_text(size=5, color="black"), axis.text.y.left=element_text(size=5,  color="black"),  legend.text=element_text(size=5, color="black"), legend.title = element_text(size=5, color="black")
        #                   , legend.position = c(0.82, 0.15), legend.background = element_rect(size=0.3, linetype="solid", colour ="black")
        , legend.position = "none"
  )
#ActivePop_PD_noPf_2<-ActivePop_PD_noPf + expand_limits(x=c(-0.08, 0.08))
ActivePop_PD_noPf_2<-ActivePop_PD_noPf + xlim(-0.05, 0.075)
plot(ActivePop_PD_noPf_2) 



ActivePop_PR <-ggplot(Active_GEE_Output_PR, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.7, width=0.6, position=position_dodge(width = 0.7))+
  geom_point(size=4, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.7)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  ggtitle("Active population (n=1,565 participants across 5,682 visits)")+
  xlab("Prevalence Ratio")+
  ylab("")+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=16, color="black"), axis.title = element_text(size =14), axis.text.x=element_text(size=11, color="black"), axis.text.y.left=element_text(size=11,  color="black"), legend.text=element_text(size=14, color="black"), legend.title = element_text(size=14, color="black"))
#ActivePop_PR_2<-ActivePop_PR + xlim(-0.9, 0.9)
plot(ActivePop_PR) 


##Now without Pf also for a zoomed version on just Pm and Po differences
Active_GEE_Output_PR_NoPf<-Active_GEE_Output_PR%>%filter(Ifxn!= "Pf")
ActivePop_PR_noPf<-ggplot(Active_GEE_Output_PR_NoPf, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.5,  width=0.4, position=position_dodge(width = 0.3))+
  geom_point(size=3, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.3)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Ratio")+
  ylab("")+
  labs(title = "Active Pop.",
       subtitle = "n=1,565 participants across 5,682 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=11, color="black"), plot.subtitle = element_text(size=9.5), axis.text.x=element_text(size=9, color="black"), axis.text.y.left=element_text(size=10,  color="black"),  legend.text=element_text(size=10, color="black"), legend.title = element_text(size=10, color="black"))
ActivePop_PR_noPf_2<-ActivePop_PR_noPf + expand_limits(x=c(-3, 6))
plot(ActivePop_PR_noPf_2) 






## Passive Visit Prevalence Association 


##Soring data by Subject_ID and then by Visit number to ensure it's formatted correctly for the clustering analysis. 
AD_Passive <- AD_Passive[order(AD_Passive$Subject_ID, AD_Passive$Visit, AD_Passive$Visit_date),]


##Crude prevalence (accounting for repeated testing;  not accounting for household clustering)

##Pm

AD_Passive_Pmnomiss <- AD_Passive%>% filter(Pm ==1 | Pm==0)

addmargins(table(AD_Passive$Pm, useNA="always"))
addmargins(table(AD_Passive_Pmnomiss$Pm, useNA="always"))



##PD
CrudePrev_GEE_Pm_Passive_Pd<-geeglm(Pm ~ 1, 
                                    data = AD_Passive, 
                                    id = Subject_ID, 
                                    family = binomial(link="identity"),
                                    corstr = "ex")

summary(CrudePrev_GEE_Pm_Passive_Pd)

##PR
CrudePrev_GEE_Pm_Passive_PR<-geeglm(Pm ~ 1, 
                                    data = AD_Passive, 
                                    id = Subject_ID, 
                                    family = binomial(link="log"),
                                    corstr = "ex")

summary(CrudePrev_GEE_Pm_Passive_PR)


### Mixed model accounting for within-subject and household clustering
CrudePrev_ME_Pm <- lmer(Pm ~ 1 + (1 | Subject_ID) + (1 | Household)  ,
                        data = AD_Passive,
                        family = binomial(link="identity"),
                        REML = FALSE)
summary(CrudePrev_ME_Pm)

CrudePrev_ME_Pm <- lmer(Pm ~ 1 + (1 | Subject_ID) + (1 | Household)  ,
                        data = AD_Passive,
                        family = binomial(link="log"),
                        REML = FALSE)
summary(CrudePrev_ME_Pm)



##Po

AD_Passive_Ponomiss <- AD_Passive%>% filter(Po ==1 | Po==0)

addmargins(table(AD_Passive$Po, useNA="always"))
addmargins(table(AD_Passive_Ponomiss$Po, useNA="always"))


CrudePrev_GEE_Po<-geeglm(Po ~ 1, 
                         data = AD_Passive, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(CrudePrev_GEE_Po)


### Mixed model accounting for within-subject and household clustering
CrudePrev_ME_Po <- lmer(Po ~ 1 + (1 | Subject_ID) + (1 | Household),
                        data = AD_Passive,
                        REML = FALSE)
summary(CrudePrev_ME_Po)



##Pf
AD_Passive_Pfnomiss <- AD_Passive%>% filter(Pf ==1 | Pf==0)

addmargins(table(AD_Passive$Pf, useNA="always"))
addmargins(table(AD_Passive_Pfnomiss$Pf, useNA="always"))



CrudePrev_GEE_Pf<-geeglm(Pf ~ 1, 
                         data = AD_Passive, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(CrudePrev_GEE_Pf)


### Mixed model accounting forwithin-subject and household clutering
CrudePrev_ME_Pf <- lmer(Pf ~ 1 + (1 | Subject_ID) + (1 | Household),
                        data = AD_Passive,
                        REML = FALSE)
summary(CrudePrev_ME_Pf)






##
##AGE CATEGORY (Age >=15 years = Ref)
##
AD_Passive_AgeNoMiss<-AD_Passive%>%filter(!is.na(AgeCat_Visit))
addmargins(table(AD_Passive_AgeNoMiss$AgeCat_Visit, AD_Passive_AgeNoMiss$Pm))
addmargins(table(AD_Passive$AgeCat_Visit, AD_Passive$Pm, useNA = "always"))


#Pm#
Age_GEE_Pm_PD<-geeglm(Pm ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Passive_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Pm_PD)

Age_GEE_Pm_PR<-geeglm(Pm ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Passive_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Pm_PR)



#Age_GEE_Pm<-geeglm(Pm ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Passive_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Pm)


#Po#

Age_GEE_Po_PD<-geeglm(Po ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Passive_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Po_PD)

Age_GEE_Po_PR<-geeglm(Po ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Passive_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Po_PR)


#Age_GEE_Po<-geeglm(Po ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Passive_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Po)


#Pf#
Age_GEE_Pf_PD<-geeglm(Pf ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Passive_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Pf_PD)

Age_GEE_Pf_PR<-geeglm(Pf ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Passive_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Pf_PR)


#Age_GEE_Pf<-geeglm(Pf ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Passive_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Pf)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Passive_Age_Pm_PD<-tidy(Age_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Age_Po_PD<-tidy(Age_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Age_Pf_PD<-tidy(Age_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_Age_Pm_PR<-tidy(Age_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Age_Po_PR<-tidy(Age_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Age_Pf_PR<-tidy(Age_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##SEX CATEGORY (Male = Ref)
##
AD_Passive_SexNoMiss<-AD_Passive%>%filter(!is.na(Sex))

#Pm#
Sex_GEE_Pm_PD<-geeglm(Pm ~ Sex, 
                      data = AD_Passive_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Pm_PD)

Sex_GEE_Pm_PR<-geeglm(Pm ~ Sex, 
                      data = AD_Passive_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Pm_PR)

#Po#
Sex_GEE_Po_PD<-geeglm(Po ~ Sex, 
                      data = AD_Passive_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Po_PD)

Sex_GEE_Po_PR<-geeglm(Po ~ Sex, 
                      data = AD_Passive_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Po_PR)

#Pf#
Sex_GEE_Pf_PD<-geeglm(Pf ~ Sex, 
                      data = AD_Passive_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Pf_PD)

Sex_GEE_Pf_PR<-geeglm(Pf ~ Sex, 
                      data = AD_Passive_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Pf_PR)



##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator

#Prevalence Differences 
Coeff_Passive_Sex_Pm_PD<-tidy(Sex_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Sex_Po_PD<-tidy(Sex_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Sex_Pf_PD<-tidy(Sex_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios
Coeff_Passive_Sex_Pm_PR<-tidy(Sex_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Sex_Po_PR<-tidy(Sex_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Sex_Pf_PR<-tidy(Sex_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




##
##HEALTH AREA CATEGORY (1 = LINGWALA = REF CAT (urban);  Bu = 2 = Ref ; (rural);  3 = KIMPOKO  (peri-urban))
##
AD_Passive_HealthAreaNoMiss<-AD_Passive%>%filter(!is.na(HealthArea))

#Pm#
HealthArea_GEE_Pm_PD<-geeglm(Pm ~ as.factor(HealthArea), 
                             data = AD_Passive_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Pm_PD)

HealthArea_GEE_Pm_PR<-geeglm(Pm ~ as.factor(HealthArea), 
                             data = AD_Passive_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Pm_PR)

#Po#
HealthArea_GEE_Po_PD<-geeglm(Po ~ as.factor(HealthArea),  
                             data = AD_Passive_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Po_PD)

HealthArea_GEE_Po_PR<-geeglm(Po ~ as.factor(HealthArea),  
                             data = AD_Passive_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Po_PR)

#Pf#
HealthArea_GEE_Pf_PD<-geeglm(Pf ~ as.factor(HealthArea),  
                             data = AD_Passive_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Pf_PD)

HealthArea_GEE_Pf_PR<-geeglm(Pf ~ as.factor(HealthArea),  
                             data = AD_Passive_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Pf_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Passive_HealthArea_Pm_PD<-tidy(HealthArea_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_HealthArea_Po_PD<-tidy(HealthArea_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_HealthArea_Pf_PD<-tidy(HealthArea_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios
Coeff_Passive_HealthArea_Pm_PR<-tidy(HealthArea_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_HealthArea_Po_PR<-tidy(HealthArea_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_HealthArea_Pf_PR<-tidy(HealthArea_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")





##
##WEALTH CATEGORY (Average Wealth [3] = Ref)
##

AD_Passive_WealthCatNoMiss<-AD_Passive%>%filter(!is.na(WealthCat))
## comparing poor vs. average, and rich vs. average  (average = 1; poor = 2; rich = 3)

#Pm#
Wealth_GEE_Pm_PD<-geeglm(Pm ~ as.factor(WealthCat), 
                         data = AD_Passive_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Pm_PD)

Wealth_GEE_Pm_PR<-geeglm(Pm ~ as.factor(WealthCat), 
                         data = AD_Passive_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Pm_PR)


#Po#
Wealth_GEE_Po_PD<-geeglm(Po ~ as.factor(WealthCat), 
                         data = AD_Passive_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Po_PD)

Wealth_GEE_Po_PR<-geeglm(Po ~ as.factor(WealthCat), 
                         data = AD_Passive_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Po_PR)

#Pf#
Wealth_GEE_Pf_PD<-geeglm(Pf ~ as.factor(WealthCat), 
                         data = AD_Passive_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Pf_PD)

Wealth_GEE_Pf_PR<-geeglm(Pf ~ as.factor(WealthCat), 
                         data = AD_Passive_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Pf_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Difference 
Coeff_Passive_Wealth_Pm_PD<-tidy(Wealth_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Wealth_Po_PD<-tidy(Wealth_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Wealth_Pf_PD<-tidy(Wealth_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_Wealth_Pm_PR<-tidy(Wealth_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Wealth_Po_PR<-tidy(Wealth_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Wealth_Pf_PR<-tidy(Wealth_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##FEVER CATEGORY (No Fever = Ref)
##
AD_Passive_FeverNoMiss<-AD_Passive%>%filter(!is.na(Feverpos))

addmargins(table(AD_Passive_FeverNoMiss$Feverpos, AD_Passive_FeverNoMiss$Pm, useNA = "always"))
addmargins(table(AD_Passive_FeverNoMiss$Feverpos, AD_Passive_FeverNoMiss$Po, useNA = "always"))

#Pm#
Fever_GEE_Pm_PD<-geeglm(Pm ~ Feverpos, 
                        data = AD_Passive_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Pm_PD)

Fever_GEE_Pm_PR<-geeglm(Pm ~ Feverpos, 
                        data = AD_Passive_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Pm_PR)


#Po#
Fever_GEE_Po_PD<-geeglm(Po ~ Feverpos, 
                        data = AD_Passive_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Po_PD)

Fever_GEE_Po_PR<-geeglm(Po ~ Feverpos, 
                        data = AD_Passive_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Po_PR)

#Pf#
Fever_GEE_Pf_PD<-geeglm(Pf ~ Feverpos, 
                        data = AD_Passive_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Pf_PD)

Fever_GEE_Pf_PR<-geeglm(Pf ~ Feverpos, 
                        data = AD_Passive_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Passive_Fever_Pm_PD<-tidy(Fever_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Fever_Po_PD<-tidy(Fever_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Fever_Pf_PD<-tidy(Fever_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_Fever_Pm_PR<-tidy(Fever_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Fever_Po_PR<-tidy(Fever_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Fever_Pf_PR<-tidy(Fever_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




write.table(Coeff_Passive_Fever_Pm_PD, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Passive_Fever_Po_PD, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Passive_Fever_Pf_PD, "clipboard", sep="\t", col.names = F )

write.table(Coeff_Passive_Fever_Pm_PR, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Passive_Fever_Po_PR, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Passive_Fever_Pf_PR, "clipboard", sep="\t", col.names = F )



##
##RDT-positive CATEGORY (Not RDT Positive = Ref)
##
AD_Passive_RDTNoMiss<-AD_Passive%>%filter(!is.na(RDTpos))

#Pm#
RDT_GEE_Pm_PD<-geeglm(Pm ~ RDTpos, 
                      data = AD_Passive_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Pm_PD)

RDT_GEE_Pm_PR<-geeglm(Pm ~ RDTpos, 
                      data = AD_Passive_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Pm_PR)


#Po#
RDT_GEE_Po_PD<-geeglm(Po ~ RDTpos, 
                      data = AD_Passive_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Po_PD)

RDT_GEE_Po_PR<-geeglm(Po ~ RDTpos, 
                      data = AD_Passive_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Po_PR)

#Pf#
RDT_GEE_Pf_PD<-geeglm(Pf ~ RDTpos, 
                      data = AD_Passive_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Pf_PD)

RDT_GEE_Pf_PR<-geeglm(Pf ~ RDTpos, 
                      data = AD_Passive_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Passive_RDT_Pm_PD<-tidy(RDT_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_RDT_Po_PD<-tidy(RDT_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_RDT_Pf_PD<-tidy(RDT_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_RDT_Pm_PR<-tidy(RDT_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_RDT_Po_PR<-tidy(RDT_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_RDT_Pf_PR<-tidy(RDT_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##P.falciparum Co-infection (Neg Pf = Ref)
##
AD_Passive_PfNoMiss<-AD_Passive%>%filter(!is.na(Pf))
#Pm#
Pf_GEE_Pm_PD<-geeglm(Pm ~ Pf, 
                     data = AD_Passive_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="identity"),
                     corstr = "ex")
summary(Pf_GEE_Pm_PD)

Pf_GEE_Pm_PR<-geeglm(Pm ~ Pf, 
                     data = AD_Passive_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="log"),
                     corstr = "ex")
summary(Pf_GEE_Pm_PR)

#Po#
Pf_GEE_Po_PD<-geeglm(Po ~ Pf, 
                     data = AD_Passive_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="identity"),
                     corstr = "ex")
summary(Pf_GEE_Po_PD)

Pf_GEE_Po_PR<-geeglm(Po ~ Pf, 
                     data = AD_Passive_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="log"),
                     corstr = "ex")
summary(Pf_GEE_Po_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence DIfferences 
Coeff_Passive_Pfcoinfec_Pm_PD<-tidy(Pf_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Pfcoinfec_Po_PD<-tidy(Pf_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_Pfcoinfec_Pm_PR<-tidy(Pf_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Pfcoinfec_Po_PR<-tidy(Pf_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")



##
##SEASONALITY CATEGORY (Rainy Season = Ref)
##
AD_Passive_SeasonNoMiss<-AD_Passive%>%filter(!is.na(Season_dry))
#Pm#
DrySeason_GEE_Pm_PD<-geeglm(Pm ~ Season_dry, 
                            data = AD_Passive_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Pm_PD)

DrySeason_GEE_Pm_PR<-geeglm(Pm ~ Season_dry, 
                            data = AD_Passive_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Pm_PR)

#Po#
DrySeason_GEE_Po_PD<-geeglm(Po ~ Season_dry, 
                            data = AD_Passive_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Po_PD)

DrySeason_GEE_Po_PR<-geeglm(Po ~ Season_dry, 
                            data = AD_Passive_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Po_PR)

#Pf#
DrySeason_GEE_Pf_PD<-geeglm(Pf ~ Season_dry, 
                            data = AD_Passive_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Pf_PD)

DrySeason_GEE_Pf_PR<-geeglm(Pf ~ Season_dry, 
                            data = AD_Passive_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalecne Differneces 
Coeff_Passive_DrySeason_Pm_PD<-tidy(DrySeason_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_DrySeason_Po_PD<-tidy(DrySeason_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_DrySeason_Pf_PD<-tidy(DrySeason_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalecne Ratios 
Coeff_Passive_DrySeason_Pm_PR<-tidy(DrySeason_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_DrySeason_Po_PR<-tidy(DrySeason_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_DrySeason_Pf_PR<-tidy(DrySeason_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



####
##Anemia(No Anemia= Ref)

##
AD_Passive_AnemiaNoMiss<-AD_Passive%>%filter(!is.na(anemia_WHO))


addmargins(table(AD_Passive$anemia_WHO, AD_Passive$Pm, useNA="always"))
addmargins(table(AD_Passive$anemia_WHO, AD_Passive$Po, useNA="always"))

addmargins(table(AD_Passive$anemia_WHO, AD_Passive$Pm_Species_Mono, useNA="always"))
addmargins(table(AD_Passive$anemia_WHO, AD_Passive$Po_Species_Mono, useNA="always"))

addmargins(table(AD_Passive_AnemiaNoMiss$anemia_WHO, AD_Passive_AnemiaNoMiss$Pm, useNA="always"))
addmargins(table(AD_Passive_AnemiaNoMiss$anemia_WHO, AD_Passive_AnemiaNoMiss$Po, useNA="always"))

AD_Passive_AnemiaNoMiss<-AD_Passive_AnemiaNoMiss%>%mutate(Anemia_ModSevere=ifelse(anemia_WHO==0|anemia_WHO==1, 0, 
                                                                                  ifelse(anemia_WHO==2|anemia_WHO==3, 1, NA)))%>%
                                                   mutate(Anemia_Any = ifelse(anemia_WHO==1 | anemia_WHO==2|anemia_WHO==3, 1,
                                                                                  ifelse(anemia_WHO==0, 0, NA)))

addmargins(table(AD_Passive_AnemiaNoMiss$Anemia_ModSevere, AD_Passive_AnemiaNoMiss$Pm, useNA="always"))
addmargins(table(AD_Passive_AnemiaNoMiss$Anemia_ModSevere, AD_Passive_AnemiaNoMiss$Po, useNA="always"))
addmargins(table(AD_Passive_AnemiaNoMiss$Anemia_Any, AD_Passive_AnemiaNoMiss$Pm, useNA="always"))
addmargins(table(AD_Passive_AnemiaNoMiss$Anemia_Any, AD_Passive_AnemiaNoMiss$Po, useNA="always"))


#Pm#
#ModSevere vs. Mild/None

Anemia_GEE_Passive_Pm_PD<-geeglm(Pm ~ Anemia_ModSevere, 
                                 data = AD_Passive_AnemiaNoMiss, 
                                 id = Subject_ID, 
                                 family = binomial(link="identity"),
                                 corstr = "ex")
summary(Anemia_GEE_Passive_Pm_PD)

Anemia_GEE_Passive_Pm_PR<-geeglm(Pm ~ Anemia_ModSevere, 
                                 data = AD_Passive_AnemiaNoMiss, 
                                 id = Subject_ID, 
                                 family = binomial(link="log"),
                                 corstr = "ex")
summary(Anemia_GEE_Passive_Pm_PR)

#Any Anemia vs. No Anemia
Anemia_GEE_Passive_Pm2_PD<-geeglm(Pm ~ Anemia_Any, 
                                  data = AD_Passive_AnemiaNoMiss, 
                                  id = Subject_ID, 
                                  family = binomial(link="identity"),
                                  corstr = "ex")
summary(Anemia_GEE_Passive_Pm2_PD)

Anemia_GEE_Passive_Pm2_PR<-geeglm(Pm ~ Anemia_Any, 
                                  data = AD_Passive_AnemiaNoMiss, 
                                  id = Subject_ID, 
                                  family = binomial(link="log"),
                                  corstr = "ex")
summary(Anemia_GEE_Passive_Pm2_PR)                                      


##Now with clustering accounted for
Anemia_GEE_Passive_Pm2_PD2 <- glmer(Pm ~ Anemia_Any + (1 | Subject_ID) + (1 | Household),
                                    data = AD_Passive_AnemiaNoMiss,
                                    family = binomial(link = "log"))
summary(Anemia_GEE_Passive_Pm2_PD2)


#Po#
#ModSevere vs. Mild/None

Anemia_GEE_Passive_Po_PD<-geeglm(Po ~ Anemia_ModSevere, 
                                 data = AD_Passive_AnemiaNoMiss, 
                                 id = Subject_ID, 
                                 family = binomial(link="identity"),
                                 corstr = "ex")
summary(Anemia_GEE_Passive_Po_PD)

Anemia_GEE_Passive_Po_PR<-geeglm(Po ~ Anemia_ModSevere, 
                                 data = AD_Passive_AnemiaNoMiss, 
                                 id = Subject_ID, 
                                 family = binomial(link="log"),
                                 corstr = "ex")
summary(Anemia_GEE_Passive_Po_PR)

#Any Anemia vs. No Anemia
Anemia_GEE_Passive_Po2_PD<-geeglm(Po ~ Anemia_Any, 
                                  data = AD_Passive_AnemiaNoMiss, 
                                  id = Subject_ID, 
                                  family = binomial(link="identity"),
                                  corstr = "ex")
summary(Anemia_GEE_Passive_Po2_PD)

Anemia_GEE_Passive_Po2_PR<-geeglm(Po ~ Anemia_Any, 
                                  data = AD_Passive_AnemiaNoMiss, 
                                  id = Subject_ID, 
                                  family = binomial(link="log"),
                                  corstr = "ex")
summary(Anemia_GEE_Passive_Po2_PR)    

#Pf#
#ModSevere vs. Mild/None

Anemia_GEE_Passive_Pf_PD<-geeglm(Pf ~ Anemia_ModSevere, 
                                 data = AD_Passive_AnemiaNoMiss, 
                                 id = Subject_ID, 
                                 family = binomial(link="identity"),
                                 corstr = "ex")
summary(Anemia_GEE_Passive_Pf_PD)

Anemia_GEE_Passive_Pf_PR<-geeglm(Pf ~ Anemia_ModSevere, 
                                 data = AD_Passive_AnemiaNoMiss, 
                                 id = Subject_ID, 
                                 family = binomial(link="log"),
                                 corstr = "ex")
summary(Anemia_GEE_Passive_Pf_PR)

#Any Anemia vs. No Anemia
Anemia_GEE_Passive_Pf2_PD<-geeglm(Pf ~ Anemia_Any, 
                                  data = AD_Passive_AnemiaNoMiss, 
                                  id = Subject_ID, 
                                  family = binomial(link="identity"),
                                  corstr = "ex")
summary(Anemia_GEE_Passive_Pf2_PD)

Anemia_GEE_Passive_Pf2_PR<-geeglm(Pf ~ Anemia_Any, 
                                  data = AD_Passive_AnemiaNoMiss, 
                                  id = Subject_ID, 
                                  family = binomial(link="log"),
                                  corstr = "ex")
summary(Anemia_GEE_Passive_Pf2_PR)      



##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Moderate/Severe vs. Mild / No Anemia
#Prevalence Differneces 
Coeff_Passive_Anemia_Pm_PD<-tidy(Anemia_GEE_Passive_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Po_PD<-tidy(Anemia_GEE_Passive_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Pf_PD<-tidy(Anemia_GEE_Passive_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_Anemia_Pm_PR<-tidy(Anemia_GEE_Passive_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Po_PR<-tidy(Anemia_GEE_Passive_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Pf_PR<-tidy(Anemia_GEE_Passive_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")


#Any Anemia vs. No Anemia 
#Prevalence Differneces 
Coeff_Passive_Anemia_Pm2_PD<-tidy(Anemia_GEE_Passive_Pm2_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Po2_PD<-tidy(Anemia_GEE_Passive_Po2_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Pf2_PD<-tidy(Anemia_GEE_Passive_Pf2_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Passive_Anemia_Pm2_PR<-tidy(Anemia_GEE_Passive_Pm2_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Po2_PR<-tidy(Anemia_GEE_Passive_Po2_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Passive_Anemia_Pf2_PR<-tidy(Anemia_GEE_Passive_Pf2_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



######
#Passive Visits:  Stacking all GEE Output for Estimate (95% CI) Plots

#Prevalence Differences

Passive_GEE_Output_PD<-rbind(
  Coeff_Passive_Age_Pm_PD, Coeff_Passive_Age_Po_PD, Coeff_Passive_Age_Pf_PD,
  Coeff_Passive_Sex_Pm_PD, Coeff_Passive_Sex_Po_PD, Coeff_Passive_Sex_Pf_PD,
  Coeff_Passive_HealthArea_Pm_PD, Coeff_Passive_HealthArea_Po_PD, Coeff_Passive_HealthArea_Pf_PD,
  Coeff_Passive_Wealth_Pm_PD, Coeff_Passive_Wealth_Po_PD, Coeff_Passive_Wealth_Pf_PD,
  Coeff_Passive_Fever_Pm_PD, Coeff_Passive_Fever_Po_PD, Coeff_Passive_Fever_Pf_PD,
  Coeff_Passive_RDT_Pm_PD, Coeff_Passive_RDT_Po_PD, Coeff_Passive_RDT_Pf_PD,
  Coeff_Passive_Pfcoinfec_Pm_PD, Coeff_Passive_Pfcoinfec_Po_PD,
  Coeff_Passive_DrySeason_Pm_PD, Coeff_Passive_DrySeason_Po_PD, Coeff_Passive_DrySeason_Pf_PD,
  Coeff_Passive_Anemia_Pm2_PD, Coeff_Passive_Anemia_Po2_PD, Coeff_Passive_Anemia_Pf2_PD,
  Coeff_Passive_Anemia_Pm_PD, Coeff_Passive_Anemia_Po_PD, Coeff_Passive_Anemia_Pf_PD
)

write.table(Passive_GEE_Output_PD, "clipboard", sep="\t", col.names = T )



#Prevalence Ratios

Passive_GEE_Output_PR<-rbind(
  Coeff_Passive_Age_Pm_PR, Coeff_Passive_Age_Po_PR, Coeff_Passive_Age_Pf_PR,
  Coeff_Passive_Sex_Pm_PR, Coeff_Passive_Sex_Po_PR, Coeff_Passive_Sex_Pf_PR,
  Coeff_Passive_HealthArea_Pm_PR, Coeff_Passive_HealthArea_Po_PR, Coeff_Passive_HealthArea_Pf_PR,
  #Coeff_Passive_Village_Pm_PD, Coeff_Passive_Village_Po_PD, Coeff_Passive_Village_Pf_PD,
  Coeff_Passive_Wealth_Pm_PR, Coeff_Passive_Wealth_Po_PR, Coeff_Passive_Wealth_Pf_PR,
  Coeff_Passive_Fever_Pm_PR, Coeff_Passive_Fever_Po_PR, Coeff_Passive_Fever_Pf_PR,
  Coeff_Passive_RDT_Pm_PR, Coeff_Passive_RDT_Po_PR, Coeff_Passive_RDT_Pf_PR,
  Coeff_Passive_Pfcoinfec_Pm_PR, Coeff_Passive_Pfcoinfec_Po_PR,
  Coeff_Passive_DrySeason_Pm_PR, Coeff_Passive_DrySeason_Po_PR, Coeff_Passive_DrySeason_Pf_PR,
  Coeff_Passive_Anemia_Pm2_PR, Coeff_Passive_Anemia_Po2_PR, Coeff_Passive_Anemia_Pf2_PR,
  Coeff_Passive_Anemia_Pm_PR, Coeff_Passive_Anemia_Po_PR, Coeff_Passive_Anemia_Pf_PR
)

write.table(Passive_GEE_Output_PR, "clipboard", sep="\t", col.names = T )




#############
##  Plotting Prevalence Differences and 95% CIs by Population

Passive_GEE_Output_PD$Ifxn <- factor(Passive_GEE_Output_PD$Ifxn, levels=c("Pm", "Po", "Pf"))
Passive_GEE_Output_PD$Ifxn <- factor(Passive_GEE_Output_PD$Ifxn, levels=c("Pm", "Po", "Pf"))

Passive_GEE_Output_PD$term <- recode_factor(Passive_GEE_Output_PD$term, 
                                            "relevel(as.factor(AgeCat_Visit), ref = \"3\")1" = "Age <5 vs 15+", 
                                            "relevel(as.factor(AgeCat_Visit), ref = \"3\")2" = "Age 5-14 vs 15+",
                                            "Sex" = "Sex (F vs M)",
                                            "as.factor(HealthArea)2" = "Rural vs Urban", 
                                            "as.factor(HealthArea)3" = "Peri-urban vs Urban", 
                                            "as.factor(WealthCat)2" = "Poor vs Avg. Wealth",
                                            "as.factor(WealthCat)3" = "Wealthy vs Avg. Wealth",
                                            "Season_dry" =  "Dry vs Rainy Season",
                                            "RDTpos" = "RDT+ vs RDT-",
                                            "Pf" = "Pf coinfection (Y vs N)",
                                            "Feverpos" = "Fever (Y vs N)",
                                            "Anemia_Any" = "Anemic vs Not Anemic",
                                            "Anemia_ModSevere" = "Mod./Severe vs Mild/No Anemia")

Passive_GEE_Output_PR$Ifxn <- factor(Passive_GEE_Output_PR$Ifxn, levels=c("Pm", "Po", "Pf"))
Passive_GEE_Output_PR$Ifxn <- factor(Passive_GEE_Output_PR$Ifxn, levels=c("Pm", "Po", "Pf"))

Passive_GEE_Output_PR$term <- recode_factor(Passive_GEE_Output_PR$term, 
                                            "relevel(as.factor(AgeCat_Visit), ref = \"3\")1" = "Age <5 vs 15+", 
                                            "relevel(as.factor(AgeCat_Visit), ref = \"3\")2" = "Age 5-14 vs 15+",
                                            "Sex" = "Sex (F vs M)",
                                            "as.factor(HealthArea)2" = "Rural vs Urban", 
                                            "as.factor(HealthArea)3" = "Peri-urban vs Urban", 
                                            "as.factor(WealthCat)2" = "Poor vs Avg. Wealth",
                                            "as.factor(WealthCat)3" = "Wealthy vs Avg. Wealth",
                                            "Season_dry" =  "Dry vs Rainy Season",
                                            "RDTpos" = "RDT+ vs RDT-",
                                            "Pf" = "Pf coinfection (Y vs N)",
                                            "Feverpos" = "Fever (Y vs N)",
                                            "Anemia_Any" = "Anemic vs Not Anemic", 
                                            "Anemia_ModSevere" = "Mod./Severe vs Mild/No Anemia")



##checking the ordering of variables for gggplotting: 
Passive_GEE_Output_PD$term %>%levels()
Passive_GEE_Output_PR$term %>%levels()
Passive_GEE_Output_PD$term %>%levels()
#  Total_GEE_Output_PD$term %>%levels()


#set order of levels to appear along y-axis
Passive_GEE_Output_PD <- Passive_GEE_Output_PD%>%
  mutate(term = term %>% 
           fct_relevel("Wealthy vs Avg. Wealth",
                       "Poor vs Avg. Wealth",
                       "Dry vs Rainy Season",
                       "Mod./Severe vs Mild/No Anemia",
                       "Anemic vs Not Anemic",
                       "Fever (Y vs N)",
                       "Pf coinfection (Y vs N)",
                       "RDT+ vs RDT-",
                       "Peri-urban vs Urban",
                       "Rural vs Urban",
                       "Sex (F vs M)",
                       "Age 5-14 vs 15+",
                       "Age <5 vs 15+", 
           ))

Passive_GEE_Output_PR <- Passive_GEE_Output_PR%>%
  mutate(term = term %>% 
           fct_relevel("Wealthy vs Avg. Wealth",
                       "Poor vs Avg. Wealth",
                       "Dry vs Rainy Season",
                       "Mod./Severe vs Mild/No Anemia",
                       "Anemic vs Not Anemic",
                       "Fever (Y vs N)",
                       "Pf coinfection (Y vs N)",
                       "RDT+ vs RDT-",
                       "Peri-urban vs Urban",
                       "Rural vs Urban",
                       "Sex (F vs M)",
                       "Age 5-14 vs 15+",
                       "Age <5 vs 15+", 
           ))



PassivePop_PD<-ggplot(Passive_GEE_Output_PD, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn)))+ 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.5,  width=0.4, position=position_dodge(width = 0.7))+
  geom_point(size=2, shape=21, colour="black", stroke = 0.6, position=position_dodge(width = 0.7)) +
  geom_vline(xintercept=0.0, lty=2, lwt=2 )+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Difference")+
  ylab("")+
  theme(plot.title =element_text(size=5, color="black"), plot.subtitle = element_text(size=5), axis.text.x=element_text(size=5, color="black"), axis.text.y.left=element_text(size=5,  color="black"), legend.text=element_text(size=5, color="black"), legend.title = element_text(size=5, color="black")
        #        , panel.background = element_rect(fill = "white", color = "black")
        #       , panel.grid.major.x = element_line(size = 0.3, linetype = 'solid', colour = "#D3D3D3"), panel.grid.minor.x = element_line(size = 0.2, linetype = 'solid', colour = "#D3D3D3"), panel.grid.major.y = element_blank()
        , legend.position = c(0.83, 0.88), legend.background = element_rect(size=0.3, linetype="solid", colour ="black")
        #          , legend.position = "none"
  )

#PassivePop_PD_noPf_2<-PassivePop_PD_noPf + expand_limits(x=c(-0.08, 0.08))
PassivePop_PD_2<-PassivePop_PD + xlim(-0.25, 0.75)
plot(PassivePop_PD_2) 


##Now without Pf also for a zoomed version on just Pm and Po differences
Passive_GEE_Output_PD_NoPf<-Passive_GEE_Output_PD%>%filter(Ifxn!= "Pf")
PassivePop_PD_noPf<-ggplot(Passive_GEE_Output_PD_NoPf, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn)))+ 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.4,  width=0.4, position=position_dodge(width = 0.3))+
  geom_point(size=2, shape=21, colour="black", stroke = 0.6, position=position_dodge(width = 0.3)) +
  geom_vline(xintercept=0.0, lty=2, lwt=2 )+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Difference")+
  ylab("")+
  #labs(title = "Clinic Sub-Pop.",
  #     subtitle = "n=1,050 participants across 3,407 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=5, color="black"), plot.subtitle = element_text(size=5), axis.text.x=element_text(size=5, color="black"), axis.text.y.left=element_text(size=5,  color="black"), legend.text=element_text(size=5, color="black"), legend.title = element_text(size=5, color="black")
        #        , panel.background = element_rect(fill = "white", color = "black")
        #       , panel.grid.major.x = element_line(size = 0.3, linetype = 'solid', colour = "#D3D3D3"), panel.grid.minor.x = element_line(size = 0.2, linetype = 'solid', colour = "#D3D3D3"), panel.grid.major.y = element_blank()
        , legend.position = c(0.89, 0.88), legend.background = element_rect(size=0.3, linetype="solid", colour ="black")
        #          , legend.position = "none"
  )

#PassivePop_PD_noPf_2<-PassivePop_PD_noPf + expand_limits(x=c(-0.08, 0.08))
PassivePop_PD_noPf_2<-PassivePop_PD_noPf + xlim(-0.05, 0.075)
plot(PassivePop_PD_noPf_2) 



PassivePop_PR <-ggplot(Passive_GEE_Output_PR, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.7, width=0.6, position=position_dodge(width = 0.7))+
  geom_point(size=4, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.7)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  #ggtitle("Survey-based population (n=1,565 participants across 5,682 visits)")+
  xlab("Prevalence Ratio")+
  ylab("")+
  labs(title = "Clinic Sub-Pop.",
       subtitle = "n=1,050 participants across 3,407 visits")+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=16, color="black"), axis.title = element_text(size =14), axis.text.x=element_text(size=11, color="black"), axis.text.y.left=element_text(size=11,  color="black"), legend.text=element_text(size=14, color="black"), legend.title = element_text(size=14, color="black"))
plot(PassivePop_PR) 


##Now without Pf also for a zoomed version on just Pm and Po differences
Passive_GEE_Output_PR_NoPf<-Passive_GEE_Output_PR%>%filter(Ifxn!= "Pf")
PassivePop_PR_noPf<-ggplot(Passive_GEE_Output_PR_NoPf, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.5,  width=0.4, position=position_dodge(width = 0.3))+
  geom_point(size=3, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.3)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Ratio")+
  ylab("")+
  labs(title = "Clinic Sub-Pop.",
       subtitle = "n=1,050 participants across 3,407 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=7, color="black"), plot.subtitle = element_text(size=9.5), axis.text.x=element_text(size=9, color="black"), axis.text.y.left=element_text(size=10,  color="black"),  legend.text=element_text(size=10, color="black"), legend.title = element_text(size=10, color="black"))
PassivePop_PR_noPf_2<-PassivePop_PR_noPf + expand_limits(x=c(-3, 6))
plot(PassivePop_PR_noPf_2) 




plot(ActivePop_PD_noPf_2)
plot(PassivePop_PD_noPf_2)













###Aggregating the non-pf and Pf Prevalence Differences Associated Factor Plots as 1 figure 
Aggregate_PrevDiffs_PmPo<-gridExtra::grid.arrange(ActivePop_PD_noPf_2, PassivePop_PD_noPf_2, ncol=2, nrow=1) 

Aggregate_PrevDiffs_PmPoPf<-gridExtra::grid.arrange(ActivePop_PD_2, PassivePop_PD, ncol=2, nrow=1) 




#### Exporting High Res Figures for Publication -- Active and Passive Populations, Non-Pf for manuscript,  with Pf for supplemental fig. 

## Non-Pf Figure - Active vs. Passive Strata 
ggsave("Active_PrevDiffs_PmPo_6Sep23.pdf", ActivePop_PD_noPf_2, width = 90, height = 120, units = "mm") 
ggsave("Passive_PrevDiffs_PmPo_6Sep23.pdf", PassivePop_PD_noPf_2, width = 90, height = 120, units = "mm")




## With Pf Figure - Active vs. Passive Strata 
ggsave("Active_PrevDiffs_PmPoPf_6Sep23.pdf", ActivePop_PD_2, width = 120, height = 120, units = "mm") 
ggsave("Passive_PrevDiffs_PmPoPf_6Sep23.pdf", PassivePop_PD_2, width = 120, height = 120, units = "mm")


























## TOTAL  Visit Prevalence Association 


##Soring data by Subject_ID and then by Visit number to ensure it's formatted correctly for the clustering analysis. 
AD_Total <- AD_Total[order(AD_Total$Subject_ID, AD_Total$Visit, AD_Total$Visit_date),]


##Crude prevalence (accounting for repeated testing;  not accounting for household clustering)

##Pm

##PD
CrudePrev_GEE_Pm_Total_Pd<-geeglm(Pm ~ 1, 
                                  data = AD_Total, 
                                  id = Subject_ID, 
                                  family = binomial(link="identity"),
                                  corstr = "ex")

summary(CrudePrev_GEE_Pm_Total_Pd)

##PR
CrudePrev_GEE_Pm_Total_PR<-geeglm(Pm ~ 1, 
                                  data = AD_Total, 
                                  id = Subject_ID, 
                                  family = binomial(link="log"),
                                  corstr = "ex")

summary(CrudePrev_GEE_Pm_Total_PR)


### Mixed model accounting for within-subject and household clustering
CrudePrev_ME_Pm <- lmer(Pm ~ 1 + (1 | Subject_ID) + (1 | Household)  ,
                        data = AD_Total,
                        family = binomial(link="identity"),
                        REML = FALSE)
summary(CrudePrev_ME_Pm)

CrudePrev_ME_Pm <- lmer(Pm ~ 1 + (1 | Subject_ID) + (1 | Household)  ,
                        data = AD_Total,
                        family = binomial(link="log"),
                        REML = FALSE)
summary(CrudePrev_ME_Pm)



##Po
CrudePrev_GEE_Po<-geeglm(Po ~ 1, 
                         data = AD_Total, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(CrudePrev_GEE_Po)


### Mixed model accounting for within-subject and household clustering
CrudePrev_ME_Po <- lmer(Po ~ 1 + (1 | Subject_ID) + (1 | Household),
                        data = AD_Total,
                        REML = FALSE)
summary(CrudePrev_ME_Po)



##Pf
CrudePrev_GEE_Pf<-geeglm(Pf ~ 1, 
                         data = AD_Total, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(CrudePrev_GEE_Pf)


### Mixed model accounting forwithin-subject and household clutering
CrudePrev_ME_Pf <- lmer(Pf ~ 1 + (1 | Subject_ID) + (1 | Household),
                        data = AD_Total,
                        REML = FALSE)
summary(CrudePrev_ME_Pf)






##
##AGE CATEGORY (Age >=15 years = Ref)
##
AD_Total_AgeNoMiss<-AD_Total%>%filter(!is.na(AgeCat_Visit))
addmargins(table(AD_Total_AgeNoMiss$AgeCat_Visit, AD_Total_AgeNoMiss$Pm))
addmargins(table(AD_Total$AgeCat_Visit, AD_Total$Pm))


#Pm#
Age_GEE_Pm_PD<-geeglm(Pm ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Total_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Pm_PD)

Age_GEE_Pm_PR<-geeglm(Pm ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Total_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Pm_PR)



#Age_GEE_Pm<-geeglm(Pm ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Total_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Pm)


#Po#

Age_GEE_Po_PD<-geeglm(Po ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Total_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Po_PD)

Age_GEE_Po_PR<-geeglm(Po ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Total_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Po_PR)


#Age_GEE_Po<-geeglm(Po ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Total_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Po)


#Pf#
Age_GEE_Pf_PD<-geeglm(Pf ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Total_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Age_GEE_Pf_PD)

Age_GEE_Pf_PR<-geeglm(Pf ~ relevel(as.factor(AgeCat_Visit), ref='3'), 
                      data = AD_Total_AgeNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Age_GEE_Pf_PR)


#Age_GEE_Pf<-geeglm(Pf ~ agecat_LT5 + agecat_5_14, 
#                   data = AD_Total_AgeNoMiss, 
#                   id = Subject_ID, 
#                   family = binomial(link="identity"),
#                   corstr = "ex")
#summary(Age_GEE_Pf)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Total_Age_Pm_PD<-tidy(Age_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Age_Po_PD<-tidy(Age_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Age_Pf_PD<-tidy(Age_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_Age_Pm_PR<-tidy(Age_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Age_Po_PR<-tidy(Age_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Age_Pf_PR<-tidy(Age_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##SEX CATEGORY (Male = Ref)
##
AD_Total_SexNoMiss<-AD_Total%>%filter(!is.na(Sex))

#Pm#
Sex_GEE_Pm_PD<-geeglm(Pm ~ Sex, 
                      data = AD_Total_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Pm_PD)

Sex_GEE_Pm_PR<-geeglm(Pm ~ Sex, 
                      data = AD_Total_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Pm_PR)

#Po#
Sex_GEE_Po_PD<-geeglm(Po ~ Sex, 
                      data = AD_Total_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Po_PD)

Sex_GEE_Po_PR<-geeglm(Po ~ Sex, 
                      data = AD_Total_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Po_PR)

#Pf#
Sex_GEE_Pf_PD<-geeglm(Pf ~ Sex, 
                      data = AD_Total_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(Sex_GEE_Pf_PD)

Sex_GEE_Pf_PR<-geeglm(Pf ~ Sex, 
                      data = AD_Total_SexNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(Sex_GEE_Pf_PR)



##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator

#Prevalence Differences 
Coeff_Total_Sex_Pm_PD<-tidy(Sex_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Sex_Po_PD<-tidy(Sex_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Sex_Pf_PD<-tidy(Sex_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios
Coeff_Total_Sex_Pm_PR<-tidy(Sex_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Sex_Po_PR<-tidy(Sex_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Sex_Pf_PR<-tidy(Sex_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




##
##HEALTH AREA CATEGORY (1 = LINGWALA = REF CAT (urban);  Bu = 2 = Ref ; (rural);  3 = KIMPOKO  (peri-urban))
##
AD_Total_HealthAreaNoMiss<-AD_Total%>%filter(!is.na(HealthArea))

#Pm#
HealthArea_GEE_Pm_PD<-geeglm(Pm ~ as.factor(HealthArea), 
                             data = AD_Total_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Pm_PD)

HealthArea_GEE_Pm_PR<-geeglm(Pm ~ as.factor(HealthArea), 
                             data = AD_Total_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Pm_PR)

#Po#
HealthArea_GEE_Po_PD<-geeglm(Po ~ as.factor(HealthArea),  
                             data = AD_Total_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Po_PD)

HealthArea_GEE_Po_PR<-geeglm(Po ~ as.factor(HealthArea),  
                             data = AD_Total_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Po_PR)

#Pf#
HealthArea_GEE_Pf_PD<-geeglm(Pf ~ as.factor(HealthArea),  
                             data = AD_Total_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="identity"),
                             corstr = "ex")
summary(HealthArea_GEE_Pf_PD)

HealthArea_GEE_Pf_PR<-geeglm(Pf ~ as.factor(HealthArea),  
                             data = AD_Total_HealthAreaNoMiss, 
                             id = Subject_ID, 
                             family = binomial(link="log"),
                             corstr = "ex")
summary(HealthArea_GEE_Pf_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Total_HealthArea_Pm_PD<-tidy(HealthArea_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_HealthArea_Po_PD<-tidy(HealthArea_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_HealthArea_Pf_PD<-tidy(HealthArea_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios
Coeff_Total_HealthArea_Pm_PR<-tidy(HealthArea_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_HealthArea_Po_PR<-tidy(HealthArea_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_HealthArea_Pf_PR<-tidy(HealthArea_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")





##
##WEALTH CATEGORY (Average Wealth [3] = Ref)
##

AD_Total_WealthCatNoMiss<-AD_Total%>%filter(!is.na(WealthCat))
## comparing poor vs. average, and rich vs. average  (average = 1; poor = 2; rich = 3)

#Pm#
Wealth_GEE_Pm_PD<-geeglm(Pm ~ as.factor(WealthCat), 
                         data = AD_Total_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Pm_PD)

Wealth_GEE_Pm_PR<-geeglm(Pm ~ as.factor(WealthCat), 
                         data = AD_Total_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Pm_PR)


#Po#
Wealth_GEE_Po_PD<-geeglm(Po ~ as.factor(WealthCat), 
                         data = AD_Total_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Po_PD)

Wealth_GEE_Po_PR<-geeglm(Po ~ as.factor(WealthCat), 
                         data = AD_Total_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Po_PR)

#Pf#
Wealth_GEE_Pf_PD<-geeglm(Pf ~ as.factor(WealthCat), 
                         data = AD_Total_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(Wealth_GEE_Pf_PD)

Wealth_GEE_Pf_PR<-geeglm(Pf ~ as.factor(WealthCat), 
                         data = AD_Total_WealthCatNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(Wealth_GEE_Pf_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Difference 
Coeff_Total_Wealth_Pm_PD<-tidy(Wealth_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Wealth_Po_PD<-tidy(Wealth_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Wealth_Pf_PD<-tidy(Wealth_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_Wealth_Pm_PR<-tidy(Wealth_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Wealth_Po_PR<-tidy(Wealth_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Wealth_Pf_PR<-tidy(Wealth_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##FEVER CATEGORY (No Fever = Ref)
##
AD_Total_FeverNoMiss<-AD_Total%>%filter(!is.na(Feverpos))

addmargins(table(AD_Total_FeverNoMiss$Feverpos, AD_Total_FeverNoMiss$Pm, useNA = "always"))
addmargins(table(AD_Total_FeverNoMiss$Feverpos, AD_Total_FeverNoMiss$Po, useNA = "always"))

#Pm#
Fever_GEE_Pm_PD<-geeglm(Pm ~ Feverpos, 
                        data = AD_Total_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Pm_PD)

Fever_GEE_Pm_PR<-geeglm(Pm ~ Feverpos, 
                        data = AD_Total_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Pm_PR)


#Po#
Fever_GEE_Po_PD<-geeglm(Po ~ Feverpos, 
                        data = AD_Total_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Po_PD)

Fever_GEE_Po_PR<-geeglm(Po ~ Feverpos, 
                        data = AD_Total_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Po_PR)

#Pf#
Fever_GEE_Pf_PD<-geeglm(Pf ~ Feverpos, 
                        data = AD_Total_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="identity"),
                        corstr = "ex")
summary(Fever_GEE_Pf_PD)

Fever_GEE_Pf_PR<-geeglm(Pf ~ Feverpos, 
                        data = AD_Total_FeverNoMiss, 
                        id = Subject_ID, 
                        family = binomial(link="log"),
                        corstr = "ex")
summary(Fever_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Total_Fever_Pm_PD<-tidy(Fever_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Fever_Po_PD<-tidy(Fever_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Fever_Pf_PD<-tidy(Fever_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_Fever_Pm_PR<-tidy(Fever_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Fever_Po_PR<-tidy(Fever_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Fever_Pf_PR<-tidy(Fever_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




write.table(Coeff_Total_Fever_Pm_PD, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Total_Fever_Po_PD, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Total_Fever_Pf_PD, "clipboard", sep="\t", col.names = F )

write.table(Coeff_Total_Fever_Pm_PR, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Total_Fever_Po_PR, "clipboard", sep="\t", col.names = F )
write.table(Coeff_Total_Fever_Pf_PR, "clipboard", sep="\t", col.names = F )



##
##RDT-positive CATEGORY (Not RDT Positive = Ref)
##
AD_Total_RDTNoMiss<-AD_Total%>%filter(!is.na(RDTpos))

#Pm#
RDT_GEE_Pm_PD<-geeglm(Pm ~ RDTpos, 
                      data = AD_Total_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Pm_PD)

RDT_GEE_Pm_PR<-geeglm(Pm ~ RDTpos, 
                      data = AD_Total_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Pm_PR)


#Po#
RDT_GEE_Po_PD<-geeglm(Po ~ RDTpos, 
                      data = AD_Total_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Po_PD)

RDT_GEE_Po_PR<-geeglm(Po ~ RDTpos, 
                      data = AD_Total_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Po_PR)

#Pf#
RDT_GEE_Pf_PD<-geeglm(Pf ~ RDTpos, 
                      data = AD_Total_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="identity"),
                      corstr = "ex")
summary(RDT_GEE_Pf_PD)

RDT_GEE_Pf_PR<-geeglm(Pf ~ RDTpos, 
                      data = AD_Total_RDTNoMiss, 
                      id = Subject_ID, 
                      family = binomial(link="log"),
                      corstr = "ex")
summary(RDT_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Total_RDT_Pm_PD<-tidy(RDT_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_RDT_Po_PD<-tidy(RDT_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_RDT_Pf_PD<-tidy(RDT_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_RDT_Pm_PR<-tidy(RDT_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_RDT_Po_PR<-tidy(RDT_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_RDT_Pf_PR<-tidy(RDT_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



##
##P.falciparum Co-infection (Neg Pf = Ref)
##
AD_Total_PfNoMiss<-AD_Total%>%filter(!is.na(Pf))
#Pm#
Pf_GEE_Pm_PD<-geeglm(Pm ~ Pf, 
                     data = AD_Total_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="identity"),
                     corstr = "ex")
summary(Pf_GEE_Pm_PD)

Pf_GEE_Pm_PR<-geeglm(Pm ~ Pf, 
                     data = AD_Total_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="log"),
                     corstr = "ex")
summary(Pf_GEE_Pm_PR)

#Po#
Pf_GEE_Po_PD<-geeglm(Po ~ Pf, 
                     data = AD_Total_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="identity"),
                     corstr = "ex")
summary(Pf_GEE_Po_PD)

Pf_GEE_Po_PR<-geeglm(Po ~ Pf, 
                     data = AD_Total_PfNoMiss, 
                     id = Subject_ID, 
                     family = binomial(link="log"),
                     corstr = "ex")
summary(Pf_GEE_Po_PR)


##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence DIfferences 
Coeff_Total_Pfcoinfec_Pm_PD<-tidy(Pf_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Pfcoinfec_Po_PD<-tidy(Pf_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_Pfcoinfec_Pm_PR<-tidy(Pf_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Pfcoinfec_Po_PR<-tidy(Pf_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")



##
##SEASONALITY CATEGORY (Rainy Season = Ref)
##
AD_Total_SeasonNoMiss<-AD_Total%>%filter(!is.na(Season_dry))
#Pm#
DrySeason_GEE_Pm_PD<-geeglm(Pm ~ Season_dry, 
                            data = AD_Total_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Pm_PD)

DrySeason_GEE_Pm_PR<-geeglm(Pm ~ Season_dry, 
                            data = AD_Total_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Pm_PR)

#Po#
DrySeason_GEE_Po_PD<-geeglm(Po ~ Season_dry, 
                            data = AD_Total_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Po_PD)

DrySeason_GEE_Po_PR<-geeglm(Po ~ Season_dry, 
                            data = AD_Total_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Po_PR)

#Pf#
DrySeason_GEE_Pf_PD<-geeglm(Pf ~ Season_dry, 
                            data = AD_Total_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="identity"),
                            corstr = "ex")
summary(DrySeason_GEE_Pf_PD)

DrySeason_GEE_Pf_PR<-geeglm(Pf ~ Season_dry, 
                            data = AD_Total_SeasonNoMiss, 
                            id = Subject_ID, 
                            family = binomial(link="log"),
                            corstr = "ex")
summary(DrySeason_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalecne Differneces 
Coeff_Total_DrySeason_Pm_PD<-tidy(DrySeason_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_DrySeason_Po_PD<-tidy(DrySeason_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_DrySeason_Pf_PD<-tidy(DrySeason_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalecne Ratios 
Coeff_Total_DrySeason_Pm_PR<-tidy(DrySeason_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_DrySeason_Po_PR<-tidy(DrySeason_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_DrySeason_Pf_PR<-tidy(DrySeason_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")




##
##BED NET USE PRIOR NIGHT (No Use = Ref)
##
AD_Total_BedNetNoMiss<-AD_Total%>%filter(!is.na(BedNetPrior))
#Pm#
BedNet_GEE_Pm_PD<-geeglm(Pm ~ BedNetPrior, 
                         data = AD_Total_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(BedNet_GEE_Pm_PD)

BedNet_GEE_Pm_PR<-geeglm(Pm ~ BedNetPrior, 
                         data = AD_Total_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(BedNet_GEE_Pm_PR)

#Po#
BedNet_GEE_Po_PD<-geeglm(Po ~ BedNetPrior, 
                         data = AD_Total_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(BedNet_GEE_Po_PD)

BedNet_GEE_Po_PR<-geeglm(Po ~ BedNetPrior, 
                         data = AD_Total_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(BedNet_GEE_Po_PR)
#Pf#
BedNet_GEE_Pf_PD<-geeglm(Pf ~ BedNetPrior, 
                         data = AD_Total_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="identity"),
                         corstr = "ex")
summary(BedNet_GEE_Pf_PD)

BedNet_GEE_Pf_PR<-geeglm(Pf ~ BedNetPrior, 
                         data = AD_Total_BedNetNoMiss, 
                         id = Subject_ID, 
                         family = binomial(link="log"),
                         corstr = "ex")
summary(BedNet_GEE_Pf_PR)

##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Prevalence Differences 
Coeff_Total_BedNet_Pm_PD<-tidy(BedNet_GEE_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_BedNet_Po_PD<-tidy(BedNet_GEE_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_BedNet_Pf_PD<-tidy(BedNet_GEE_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_BedNet_Pm_PR<-tidy(BedNet_GEE_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_BedNet_Po_PR<-tidy(BedNet_GEE_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_BedNet_Pf_PR<-tidy(BedNet_GEE_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")






####
##Anemia(No Anemia= Ref)

##
AD_Total_AnemiaNoMiss<-AD_Total%>%filter(!is.na(anemia_WHO))


addmargins(table(AD_Total$anemia_WHO, AD_Total$Pm, useNA="always"))
addmargins(table(AD_Total$anemia_WHO, AD_Total$Po, useNA="always"))

addmargins(table(AD_Total$anemia_WHO, AD_Total$Pm_Species_Mono, useNA="always"))
addmargins(table(AD_Total$anemia_WHO, AD_Total$Po_Species_Mono, useNA="always"))

addmargins(table(AD_Total_AnemiaNoMiss$anemia_WHO, AD_Total_AnemiaNoMiss$Pm, useNA="always"))
addmargins(table(AD_Total_AnemiaNoMiss$anemia_WHO, AD_Total_AnemiaNoMiss$Po, useNA="always"))

AD_Total_AnemiaNoMiss<-AD_Total_AnemiaNoMiss%>%mutate(Anemia_ModSevere=ifelse(anemia_WHO==0|anemia_WHO==1, 0, 
                                                                                  ifelse(anemia_WHO==2|anemia_WHO==3, 1, NA)))%>%
                                                mutate(Anemia_Any = ifelse(anemia_WHO==1 | anemia_WHO==2|anemia_WHO==3, 1,
                                                                           ifelse(anemia_WHO==0, 0, NA)))

#Pm#
#ModSevere vs. Mild/None

Anemia_GEE_Total_Pm_PD<-geeglm(Pm ~ Anemia_ModSevere, 
                               data = AD_Total_AnemiaNoMiss, 
                               id = Subject_ID, 
                               family = binomial(link="identity"),
                               corstr = "ex")
summary(Anemia_GEE_Total_Pm_PD)

Anemia_GEE_Total_Pm_PR<-geeglm(Pm ~ Anemia_ModSevere, 
                               data = AD_Total_AnemiaNoMiss, 
                               id = Subject_ID, 
                               family = binomial(link="log"),
                               corstr = "ex")
summary(Anemia_GEE_Total_Pm_PR)

#Any Anemia vs. No Anemia
Anemia_GEE_Total_Pm2_PD<-geeglm(Pm ~ Anemia_Any, 
                                data = AD_Total_AnemiaNoMiss, 
                                id = Subject_ID, 
                                family = binomial(link="identity"),
                                corstr = "ex")
summary(Anemia_GEE_Total_Pm2_PD)

Anemia_GEE_Total_Pm2_PR<-geeglm(Pm ~ Anemia_Any, 
                                data = AD_Total_AnemiaNoMiss, 
                                id = Subject_ID, 
                                family = binomial(link="log"),
                                corstr = "ex")
summary(Anemia_GEE_Total_Pm2_PR)                                      


#Po#
#ModSevere vs. Mild/None

Anemia_GEE_Total_Po_PD<-geeglm(Po ~ Anemia_ModSevere, 
                               data = AD_Total_AnemiaNoMiss, 
                               id = Subject_ID, 
                               family = binomial(link="identity"),
                               corstr = "ex")
summary(Anemia_GEE_Total_Po_PD)

Anemia_GEE_Total_Po_PR<-geeglm(Po ~ Anemia_ModSevere, 
                               data = AD_Total_AnemiaNoMiss, 
                               id = Subject_ID, 
                               family = binomial(link="log"),
                               corstr = "ex")
summary(Anemia_GEE_Total_Po_PR)

#Any Anemia vs. No Anemia
Anemia_GEE_Total_Po2_PD<-geeglm(Po ~ Anemia_Any, 
                                data = AD_Total_AnemiaNoMiss, 
                                id = Subject_ID, 
                                family = binomial(link="identity"),
                                corstr = "ex")
summary(Anemia_GEE_Total_Po2_PD)

Anemia_GEE_Total_Po2_PR<-geeglm(Po ~ Anemia_Any, 
                                data = AD_Total_AnemiaNoMiss, 
                                id = Subject_ID, 
                                family = binomial(link="log"),
                                corstr = "ex")
summary(Anemia_GEE_Total_Po2_PR)    

#Pf#
#ModSevere vs. Mild/None

Anemia_GEE_Total_Pf_PD<-geeglm(Pf ~ Anemia_ModSevere, 
                               data = AD_Total_AnemiaNoMiss, 
                               id = Subject_ID, 
                               family = binomial(link="identity"),
                               corstr = "ex")
summary(Anemia_GEE_Total_Pf_PD)

Anemia_GEE_Total_Pf_PR<-geeglm(Pf ~ Anemia_ModSevere, 
                               data = AD_Total_AnemiaNoMiss, 
                               id = Subject_ID, 
                               family = binomial(link="log"),
                               corstr = "ex")
summary(Anemia_GEE_Total_Pf_PR)

#Any Anemia vs. No Anemia
Anemia_GEE_Total_Pf2_PD<-geeglm(Pf ~ Anemia_Any, 
                                data = AD_Total_AnemiaNoMiss, 
                                id = Subject_ID, 
                                family = binomial(link="identity"),
                                corstr = "ex")
summary(Anemia_GEE_Total_Pf2_PD)

Anemia_GEE_Total_Pf2_PR<-geeglm(Pf ~ Anemia_Any, 
                                data = AD_Total_AnemiaNoMiss, 
                                id = Subject_ID, 
                                family = binomial(link="log"),
                                corstr = "ex")
summary(Anemia_GEE_Total_Pf2_PR)      



##Outputing GEE Results: Estimate + 95% CI from Robust Sandwich Estimator
#Moderate/Severe vs. Mild / No Anemia
#Prevalence Differneces 
Coeff_Total_Anemia_Pm_PD<-tidy(Anemia_GEE_Total_Pm_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Po_PD<-tidy(Anemia_GEE_Total_Po_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Pf_PD<-tidy(Anemia_GEE_Total_Pf_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_Anemia_Pm_PR<-tidy(Anemia_GEE_Total_Pm_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Po_PR<-tidy(Anemia_GEE_Total_Po_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Pf_PR<-tidy(Anemia_GEE_Total_Pf_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")


#Any Anemia vs. No Anemia 
#Prevalence Differneces 
Coeff_Total_Anemia_Pm2_PD<-tidy(Anemia_GEE_Total_Pm2_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Po2_PD<-tidy(Anemia_GEE_Total_Po2_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Pf2_PD<-tidy(Anemia_GEE_Total_Pf2_PD, exponentiate=FALSE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")

#Prevalence Ratios 
Coeff_Total_Anemia_Pm2_PR<-tidy(Anemia_GEE_Total_Pm2_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pm")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Po2_PR<-tidy(Anemia_GEE_Total_Po2_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Po")%>%filter(!term=="(Intercept)")
Coeff_Total_Anemia_Pf2_PR<-tidy(Anemia_GEE_Total_Pf2_PR, exponentiate=TRUE, conf.int=TRUE)%>%dplyr::mutate(Ifxn="Pf")%>%filter(!term=="(Intercept)")



######
#Total Visits:  Stacking all GEE Output for Estimate (95% CI) Plots

#Prevalence Differences

Total_GEE_Output_PD<-rbind(
  Coeff_Total_Age_Pm_PD, Coeff_Total_Age_Po_PD, Coeff_Total_Age_Pf_PD,
  Coeff_Total_Sex_Pm_PD, Coeff_Total_Sex_Po_PD, Coeff_Total_Sex_Pf_PD,
  Coeff_Total_HealthArea_Pm_PD, Coeff_Total_HealthArea_Po_PD, Coeff_Total_HealthArea_Pf_PD,
  Coeff_Total_Wealth_Pm_PD, Coeff_Total_Wealth_Po_PD, Coeff_Total_Wealth_Pf_PD,
  Coeff_Total_Fever_Pm_PD, Coeff_Total_Fever_Po_PD, Coeff_Total_Fever_Pf_PD,
  Coeff_Total_RDT_Pm_PD, Coeff_Total_RDT_Po_PD, Coeff_Total_RDT_Pf_PD,
  Coeff_Total_Pfcoinfec_Pm_PD, Coeff_Total_Pfcoinfec_Po_PD,
  Coeff_Total_DrySeason_Pm_PD, Coeff_Total_DrySeason_Po_PD, Coeff_Total_DrySeason_Pf_PD,
  Coeff_Total_Anemia_Pm2_PD, Coeff_Total_Anemia_Po2_PD, Coeff_Total_Anemia_Pf2_PD, 
  Coeff_Total_Anemia_Pm_PD, Coeff_Total_Anemia_Po_PD, Coeff_Total_Anemia_Pf_PD, 
  Coeff_Total_BedNet_Pm_PD, Coeff_Total_BedNet_Po_PD, Coeff_Total_BedNet_Pf_PD
)

write.table(Total_GEE_Output_PD, "clipboard", sep="\t", col.names = T )



#Prevalence Ratios

Total_GEE_Output_PR<-rbind(
  Coeff_Total_Age_Pm_PR, Coeff_Total_Age_Po_PR, Coeff_Total_Age_Pf_PR,
  Coeff_Total_Sex_Pm_PR, Coeff_Total_Sex_Po_PR, Coeff_Total_Sex_Pf_PR,
  Coeff_Total_HealthArea_Pm_PR, Coeff_Total_HealthArea_Po_PR, Coeff_Total_HealthArea_Pf_PR,
  #Coeff_Total_Village_Pm_PD, Coeff_Total_Village_Po_PD, Coeff_Total_Village_Pf_PD,
  Coeff_Total_Wealth_Pm_PR, Coeff_Total_Wealth_Po_PR, Coeff_Total_Wealth_Pf_PR,
  Coeff_Total_Fever_Pm_PR, Coeff_Total_Fever_Po_PR, Coeff_Total_Fever_Pf_PR,
  Coeff_Total_RDT_Pm_PR, Coeff_Total_RDT_Po_PR, Coeff_Total_RDT_Pf_PR,
  Coeff_Total_Pfcoinfec_Pm_PR, Coeff_Total_Pfcoinfec_Po_PR,
  Coeff_Total_DrySeason_Pm_PR, Coeff_Total_DrySeason_Po_PR, Coeff_Total_DrySeason_Pf_PR,
  Coeff_Total_Anemia_Pm2_PR, Coeff_Total_Anemia_Po2_PR, Coeff_Total_Anemia_Pf2_PR,
  Coeff_Total_Anemia_Pm_PR, Coeff_Total_Anemia_Po_PR, Coeff_Total_Anemia_Pf_PR,
  Coeff_Total_BedNet_Pm_PR, Coeff_Total_BedNet_Po_PR, Coeff_Total_BedNet_Pf_PR
)

write.table(Total_GEE_Output_PR, "clipboard", sep="\t", col.names = T )




#############
##  Plotting Prevalence Differences and 95% CIs by Population

Total_GEE_Output_PD$Ifxn <- factor(Total_GEE_Output_PD$Ifxn, levels=c("Pm", "Po", "Pf"))
Total_GEE_Output_PD$Ifxn <- factor(Total_GEE_Output_PD$Ifxn, levels=c("Pm", "Po", "Pf"))

Total_GEE_Output_PD$term <- recode_factor(Total_GEE_Output_PD$term, 
                                          "relevel(as.factor(AgeCat_Visit), ref = \"3\")1" = "Age <5 vs 15+", 
                                          "relevel(as.factor(AgeCat_Visit), ref = \"3\")2" = "Age 5-14 vs 15+",
                                          "Sex" = "Sex (F vs M)",
                                          "as.factor(HealthArea)2" = "Rural vs Urban", 
                                          "as.factor(HealthArea)3" = "Peri-urban vs Urban", 
                                          "as.factor(WealthCat)2" = "Poor vs Avg. Wealth",
                                          "as.factor(WealthCat)3" = "Wealthy vs Avg. Wealth",
                                          "Season_dry" =  "Dry vs Rainy Season",
                                          "BedNetPrior" = "Bed Net Use (Y vs N)",
                                          "RDTpos" = "RDT+ vs RDT-",
                                          "Pf" = "Pf coinfection (Y vs N)",
                                          "Feverpos" = "Fever (Y vs N)",
                                          "Anemia_Any" = "Anemic vs Not Anemic",
                                          "Anemia_ModSevere" = "Mod./Severe vs Mild/No Anemia")

Total_GEE_Output_PR$Ifxn <- factor(Total_GEE_Output_PR$Ifxn, levels=c("Pm", "Po", "Pf"))
Total_GEE_Output_PR$Ifxn <- factor(Total_GEE_Output_PR$Ifxn, levels=c("Pm", "Po", "Pf"))

Total_GEE_Output_PR$term <- recode_factor(Total_GEE_Output_PR$term, 
                                          "relevel(as.factor(AgeCat_Visit), ref = \"3\")1" = "Age <5 vs 15+", 
                                          "relevel(as.factor(AgeCat_Visit), ref = \"3\")2" = "Age 5-14 vs 15+",
                                          "Sex" = "Sex (F vs M)",
                                          "as.factor(HealthArea)2" = "Rural vs Urban", 
                                          "as.factor(HealthArea)3" = "Peri-urban vs Urban", 
                                          "as.factor(WealthCat)2" = "Poor vs Avg. Wealth",
                                          "as.factor(WealthCat)3" = "Wealthy vs Avg. Wealth",
                                          "Season_dry" =  "Dry vs Rainy Season",
                                          "BedNetPrior" = "Bed Net Use (Y vs N)",
                                          "RDTpos" = "RDT+ vs RDT-",
                                          "Pf" = "Pf coinfection (Y vs N)",
                                          "Feverpos" = "Fever (Y vs N)",
                                          "Anemia_Any" = "Anemic vs Not Anemic",
                                          "Anemia_ModSevere" = "Mod./Severe vs Mild/No Anemia")



##checking the ordering of variables for gggplotting: 
Total_GEE_Output_PD$term %>%levels()
Total_GEE_Output_PR$term %>%levels()
Total_GEE_Output_PD$term %>%levels()
Total_GEE_Output_PD$term %>%levels()


#set order of levels to appear along y-axis
Total_GEE_Output_PD <- Total_GEE_Output_PD%>%
  mutate(term = term %>% 
           fct_relevel("Wealthy vs Avg. Wealth",
                       "Poor vs Avg. Wealth",
                       "Dry vs Rainy Season",
                       "Bed Net Use (Y vs N)",
                       "Mod./Severe vs Mild/No Anemia",
                       "Anemic vs Not Anemic",
                       "Fever (Y vs N)",
                       "Pf coinfection (Y vs N)",
                       "RDT+ vs RDT-",
                       "Peri-urban vs Urban",
                       "Rural vs Urban",
                       "Sex (F vs M)",
                       "Age 5-14 vs 15+",
                       "Age <5 vs 15+", 
           ))

Total_GEE_Output_PR <- Total_GEE_Output_PR%>%
  mutate(term = term %>% 
           fct_relevel("Wealthy vs Avg. Wealth",
                       "Poor vs Avg. Wealth",
                       "Dry vs Rainy Season",
                       "Bed Net Use (Y vs N)",
                       "Mod./Severe vs Mild/No Anemia",
                       "Anemic vs Not Anemic",
                       "Fever (Y vs N)",
                       "Pf coinfection (Y vs N)",
                       "RDT+ vs RDT-",
                       "Peri-urban vs Urban",
                       "Rural vs Urban",
                       "Sex (F vs M)",
                       "Age 5-14 vs 15+",
                       "Age <5 vs 15+", 
           ))



TotalPop_PD <-ggplot(Total_GEE_Output_PD, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.7, width=0.6, position=position_dodge(width = 0.7))+
  geom_point(size=4, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.7)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  #          ggtitle("Survey-based population (n=1,565 participants across 5,682 visits)")+
  xlab("Prevalence Difference")+
  ylab("")+
  #labs(title = "Total Visits",
  #     subtitle = "n=1,565 participants across 9,089 visits")+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=16, color="black"), axis.title = element_text(size =14), axis.text.x=element_text(size=11, color="black"), axis.text.y.left=element_text(size=11,  color="black"), legend.text=element_text(size=14, color="black"), legend.title = element_text(size=14, color="black"))
plot(TotalPop_PD) 


##Now without Pf also for a zoomed version on just Pm and Po differences
Total_GEE_Output_PD_NoPf<-Total_GEE_Output_PD%>%filter(Ifxn!= "Pf")
TotalPop_PD_noPf<-ggplot(Total_GEE_Output_PD_NoPf, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.5,  width=0.4, position=position_dodge(width = 0.3))+
  geom_point(size=3, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.3)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Difference")+
  ylab("")+
  #labs(title = "Total Visits",
  #     subtitle = "n=1,565 participants across 9,089 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=11, color="black"), plot.subtitle = element_text(size=9.5), axis.text.x=element_text(size=9, color="black"), axis.text.y.left=element_text(size=10,  color="black"),  legend.text=element_text(size=10, color="black"), legend.title = element_text(size=10, color="black")
        #TotalPop_PD_noPf_2<-TotalPop_PD_noPf + expand_limits(x=c(-0.08, 0.08))
        , legend.position= "none"
  )

TotalPop_PD_noPf_2<-TotalPop_PD_noPf + xlim(-0.05, 0.075)
plot(TotalPop_PD_noPf_2) 



TotalPop_PR <-ggplot(Total_GEE_Output_PR, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.7, width=0.6, position=position_dodge(width = 0.7))+
  geom_point(size=4, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.7)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  #ggtitle("Survey-based population (n=1,565 participants across 5,682 visits)")+
  xlab("Prevalence Ratio")+
  ylab("")+
  labs(title = "Total Visits",
       subtitle = "n=1,565 participants across 9,089 visits")+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=16, color="black"), axis.title = element_text(size =14), axis.text.x=element_text(size=11, color="black"), axis.text.y.left=element_text(size=11,  color="black"), legend.text=element_text(size=14, color="black"), legend.title = element_text(size=14, color="black"))
plot(TotalPop_PR) 


##Now without Pf also for a zoomed version on just Pm and Po differences
Total_GEE_Output_PR_NoPf<-Total_GEE_Output_PR%>%filter(Ifxn!= "Pf")
TotalPop_PR_noPf<-ggplot(Total_GEE_Output_PR_NoPf, (aes(x=estimate, y=term, col=Ifxn, fill=Ifxn))) + 
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, col=Ifxn), size=0.5,  width=0.4, position=position_dodge(width = 0.3))+
  geom_point(size=3, shape=21, colour="black", stroke = 0.8, position=position_dodge(width = 0.3)) +
  geom_vline(xintercept=0.0, lty=2)+
  scale_fill_manual(values=c("blue", "red", "orange"))+
  scale_color_manual(values=c("black", "black", "black"))+
  guides(fill = guide_legend(title = "Species"), col=FALSE)+
  # ggtitle("Survey-based Pop. (N=1,565; 5,682 visits)")+
  xlab("Prevalence Ratio")+
  ylab("")+
  labs(title = "Total Visits",
       subtitle = "n=1,565 participants across 9,089 visits")+
  #guides(fill = guide_legend(title = "Species"), col=FALSE)+
  theme(plot.title =element_text(size=11, color="black"), plot.subtitle = element_text(size=9.5), axis.text.x=element_text(size=9, color="black"), axis.text.y.left=element_text(size=10,  color="black"),  legend.text=element_text(size=10, color="black"), legend.title = element_text(size=10, color="black"))
TotalPop_PR_noPf_2<-TotalPop_PR_noPf + expand_limits(x=c(-3, 6))
plot(TotalPop_PR_noPf_2) 



plot(TotalPop_PD_noPf_2)
plot(ActivePop_PD_noPf_2)
plot(PassivePop_PD_noPf_2)

###Aggregating the non-pf Risk Factor Plots as 1 figure 
gridExtra::grid.arrange(TotalPop_PD_noPf_2, ActivePop_PD_noPf_2, PassivePop_PD_noPf_2, ncol=3, nrow=1) 


###Aggregating the non-pf Risk Factor Plots as 1 figure 
gridExtra::grid.arrange(TotalPop_PD, ActivePop_PD, PassivePop_PD, ncol=3, nrow=1) 
















    
    
    