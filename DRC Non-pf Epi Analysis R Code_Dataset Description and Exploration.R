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
library(openxlsx)
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


##########################################################
# IMPORTING DATASETS FOR ANALYSIS 

## NOTE: DATASET CREATED IN SAS FOR DATA LINKAGE, CLEANING, & PRELIMINARY WRANGLING (SEE PROGRAM: AnalysisDatasetBuild_DRCNonPf_Mar22_RS.sas)

## Dataset contains 1 row per participant per study visit, including active surveillance visits and clinic visits. Exclusion criteria already applied (ie excluding extraneous clinic visits into 2018 and participants who didn't have a DBS at baseline)


AD_All_BL_FU_VNP_Long_Sub <- read_excel("../Dataset_AD/AD_All_BL_FU_VNP_Long_Sub_Jul22.xlsx", NULL)
View(AD_All_BL_FU_VNP_Long_Sub)




##########################################################
# DEFINING NEW VARIABLES FOR ANALYSES    
AD_All_BL_FU_VNP_Long_Sub$month_visit <- format(as.Date(AD_All_BL_FU_VNP_Long_Sub$Visit_date, format="%Y-%m-%d"),"%m")


AD_All_BL_FU_VNP_Long_Sub<-AD_All_BL_FU_VNP_Long_Sub %>%
  
  ## NEW VARS:  Pm_Species;  Po_Species;  Pf_Species = differentiating mono vs. mixed infections per species. 
  mutate(Pm_Species= ifelse(Pm==1&Po==0&Pf==0, "Pm Mono", 
                            ifelse(Pm==1&Po==1&Pf==0, "Pm Mixed", 
                                   ifelse(Pm==1&Po==0&Pf==1, "Pm Mixed",
                                          ifelse(Pm==1&Po==1&Pf==1, "Pm Mixed", NA)))))%>%
  mutate(Po_Species= ifelse(Pm==0&Po==1&Pf==0, "Po Mono", 
                            ifelse(Pm==1&Po==1&Pf==0, "Po Mixed", 
                                   ifelse(Pm==0&Po==1&Pf==1, "Po Mixed",
                                          ifelse(Pm==1&Po==1&Pf==1, "Po Mixed", NA)))))%>%
  mutate(Pf_Species= ifelse(Pm==0&Po==0&Pf==1, "Pf Mono", 
                            ifelse(Pm==1&Po==0&Pf==1, "Pf Mixed", 
                                   ifelse(Pm==0&Po==1&Pf==1, "Pf Mixed",
                                          ifelse(Pm==1&Po==1&Pf==1, "Pf Mixed", NA)))))%>%
  
  ## NEW VAR:  PSpecies_Detail = specifying each type of mono or mixed infection across all species 
  mutate(PSpecies_Detail= ifelse(Pm==1&Po==0&Pf==0, "Pm Mono", 
                                 ifelse(Pm==1&Po==1&Pf==0, "Pm Po", 
                                        ifelse(Pm==1&Po==0&Pf==1, "Pm Pf",
                                               ifelse(Pm==1&Po==1&Pf==1, "Pm Po Pf", 
                                                      ifelse(Pm==0&Po==1&Pf==0, "Po Mono", 
                                                             ifelse(Pm==0&Po==1&Pf==1, "Po Pf",
                                                                    ifelse(Pm==0&Po==0&Pf==1, "Pf Mono", NA))))))))%>%
  
  ## UPDATING VAR:  Updating "Rehydrated" variable to replace any missings with "0" to denote "not rehydrated".  1 = rehydrated DNA in sample; 0 = not rehydrated DNA in sample 
  mutate_at(vars("Rehydrated"), ~replace_na(.,0))%>%

  ## NEW VARS:  Pm_Sepcies_Mono; Po_Species_Mono; Pf_Species_Mono = specifying whether each species has a mono or mixed infection.
  mutate(Pm_Species_Mono=ifelse(Pm==1&Pm_Species=="Pm Mono", 1,
                                ifelse(Pm==1&Pm_Species=="Pm Mixed", 0, NA)))%>%
  mutate(Po_Species_Mono=ifelse(Po==1&Po_Species=="Po Mono", 1,
                                ifelse(Po==1&Po_Species=="Po Mixed", 0, NA)))%>%                
  
  mutate(Pf_Species_Mono=ifelse(Pf==1&Pf_Species=="Pf Mono", 1,
                                ifelse(Pf==1&Pf_Species=="Pf Mixed", 0, NA)))%>%
  
    ##NEW VAR: Season_dry = categorizing rainy season based on visit date: 
  mutate(Season_dry=ifelse(month_visit%in%c('01','02','03','04','10','11','12'), 0, 
                           ifelse(month_visit%in%c('05','06','07','08','09'), 1, NA)))%>%
  
  ## NEW VAR:  Sx_Prior_6mo_new = recategorizing variable whether subject had malaria symptoms in prior 6 months as dichotomous: 1= yes; 0 = no. 
  mutate(Sx_Prior6mo_new=ifelse(Sx_Prior6mo==1, 1, 
                                ifelse(Sx_Prior6mo==2, 1, 
                                       ifelse(Sx_Prior6mo==0, 0, NA))))%>%
  ## NEW VAR:  AntiMalTx_Prior6mo_new = recategorizing variable whether subject took antimalarial treatment in prior 6 months as dichotomous: 1= yes; 0 = no. 
  mutate(AntiMalTx_Prior6mo_new=ifelse(AntiMalTx_Prior6mo==1, 1, 
                                       ifelse(AntiMalTx_Prior6mo==2, 1, 
                                              ifelse(AntiMalTx_Prior6mo==0, 0, NA))))%>%
  ## NEW VAR:  WealthCat = recategorizing weatlh quintiles into three categories for poorer/poorest; average, or wealthier/wealthiest.  Setting referent cagetory (1 = average wealth); 2 = poor; 3 = wealthy. 
  mutate(WealthCat = ifelse(WealthQuintile==1, 2,
                            ifelse(WealthQuintile==2, 2,
                                   ifelse(WealthQuintile==3, 1, 
                                          ifelse(WealthQuintile==4, 3,
                                                 ifelse(WealthQuintile==5, 3,NA))))))%>%
  
  ## NEW VAR:  HealthArea = recategorizing the 7 different villages into their corresponding health areas (1 = urban village; 2,3,4 = rural villages; 5,6,7 = peri-urban villages. . 
  mutate(HealthArea=   ifelse(Village==1, 1,
                              ifelse(Village==2, 2, 
                                     ifelse(Village==3, 2, 
                                            ifelse(Village==4, 2,
                                                   ifelse(Village==5, 3, 
                                                          ifelse(Village==6, 3,
                                                                 ifelse(Village==7, 3, NA))))))))




##COnfirming New Variable Coding
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Pm_Species, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Po_Species, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Pf_Species, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$HealthArea, AD_All_BL_FU_VNP_Long_Sub$Village, useNA = "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$PSpecies_Detail, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))







##########################################################

# REVIEWING THE ANALYSIS DATASET  

##########################################################


##ASSESSING MISSING VISITS BY OUTCOME FREQUENCY AND EXCLUDING MISSING VISITS FROM ANALYSIS DATASET

# Missing Visits Summary: 
MissingVisits <-AD_All_BL_FU_VNP_Long_Sub%>%filter(is.na(Visit_date))

addmargins(table(MissingVisits$Visit_Type2, useNA = "always"))
addmargins(table(MissingVisits$Visit_Type2, MissingVisits$Village, useNA = "always"))
addmargins(table(MissingVisits$Visit_Type2, MissingVisits$HealthArea, useNA = "always"))
addmargins(table(MissingVisits$Visit, MissingVisits$HealthArea, useNA = "always"))


## RESULT:  578 visits are missing visit_date and related visit information --> all 100% missing are for subjects who were lost to follow-up and didn't return for their expected follow-up household visit. 
##          All villages had some participants who were lost to follow-up over time. 


MissingVisits_bysubject<-MissingVisits%>%distinct(Subject_ID, .keep_all = TRUE)

addmargins(table(MissingVisits_bysubject$Visit_Type2, MissingVisits_bysubject$HealthArea, useNA = "always"))
prop.table(MissingVisits_bysubject$HealthArea, margin=1)

## RESULT:  Lost to follow-up by health area: 
#     Urban   Rural   Peri-urban  <NA>  Sum
#FU   116   116     143           0     375
#<NA>   0     0       0           0       0
#Sum  116   116     143           0     375

#### IN-TEXT NUMBERS FOR LOST TO FOLLOW_UP BY REGION:  
#(n=116 lost eventually /385 Urban participants at baseline = 30.12%)
#(n=143 lost eventually /517 Peri-urban participants at bsaeline = 27.66%) 
#(n=116 lost eventually / 663 Rural participants at baseline = 17.496% )


# Malaria Outcomes among non-missing visits, stratified by visit type (BL = baseline; FU = follow-up; VNP = clinic visit )
AD_All_BL_FU_VNP_Long_Sub <-AD_All_BL_FU_VNP_Long_Sub%>%filter(!is.na(Visit_date))

addmargins(table(AD_All_BL_FU_VNP_Long_Sub$PSpecies_Detail, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))

## RESULTS:  
#           BL   FU  VNP <NA>  Sum
#Pf Mono   449 1339 1871    0 3659
#Pm Mono    16   33   46    0   95
#Pm Pf      31   98   85    0  214
#Pm Po       0    2    1    0    3
#Pm Po Pf    0    6    3    0    9
#Po Mono     2   15   43    0   60
#Po Pf       4   49   48    0  101
#<NA>     1063 2575 1310    0 4948
#Sum      1565 4117 3407    0 9089


addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Pf, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Pm, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Po, AD_All_BL_FU_VNP_Long_Sub$Visit_Type2, useNA = "always"))

## RESULTS: 

# Pf: 
#        BL   FU  VNP <NA>  Sum
# 0    1081 2603 1341    0 5025
# 1     484 1492 2009    0 3985
# <NA>    0   22   57    0   79
# Sum  1565 4117 3407    0 9089


# Pm: 
#        BL   FU  VNP <NA>  Sum
# 0    1518 3955 3214    0 8687
# 1      47  139  135    0  321
# <NA>    0   23   58    0   81
# Sum  1565 4117 3407    0 9089


# Po: 
#        BL   FU  VNP <NA>  Sum
# 0    1559 4022 3254    0 8835
# 1       6   72   95    0  173
# <NA>    0   23   58    0   81
# Sum  1565 4117 3407    0 9089



addmargins(table(AD_All_BL_FU_VNP_Long_Sub$HealthArea, AD_All_BL_FU_VNP_Long_Sub$Pf, useNA= "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$HealthArea, AD_All_BL_FU_VNP_Long_Sub$Po, useNA= "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$HealthArea, AD_All_BL_FU_VNP_Long_Sub$Pm, useNA= "always"))


### RESULTS:  there are is least 1 Pm, Po, and Pf infection in each health area. (fewest infections in urban village as expected)


addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Village, AD_All_BL_FU_VNP_Long_Sub$Pf, useNA= "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Village, AD_All_BL_FU_VNP_Long_Sub$Po, useNA= "always"))
addmargins(table(AD_All_BL_FU_VNP_Long_Sub$Village, AD_All_BL_FU_VNP_Long_Sub$Pm, useNA= "always"))

### RESULTS:  there are is least 1 Pm, Po, and Pf infection in each village site. (fewest infections in urban village as expected)






########################################################################

## Separate datasets for analyses into: 

#     1. Active Survey Visits Only = "AD_Active"
#     2. Passive Unscheduled Clinic (ie VNP) Visits Only, = "AD_Passive"
#     3. and Total (Active + Passive) Visits  for Analyses = "AD_Total"

########################################################################

#Each dataset constructed represents an analysis cohort - survey-based; clinic sub-population, and total (ie all visits)

#AD_Active   = By-visit dataset for BL, FU1, FU2, and FU3 (N=1,565 subs, 6260 visits possible-- but deleted out n=578 instances of missing FU visits, so 5,682 total visits (n=1,565 baseline, n=4,117 follow-ups))
#AD_Passive  = By-visit dataset for all VNP visits (N=1,050 subs, 3,407 total VNP visits through end of Dec 2017)
#AD_Total    = By-visit dataset for all Visits (Active and Passive) (N=1,565 subs, 9667 visits possible-- but deleted out n=578 instances of missing FU visits, so 9,089 total)

AD_Active <-AD_All_BL_FU_VNP_Long_Sub%>%filter(Visit_Type1=="Act")
AD_Passive <-AD_All_BL_FU_VNP_Long_Sub%>%filter(Visit_Type1=="Pas")
AD_Total <-AD_All_BL_FU_VNP_Long_Sub

##Quick look at malaria species infections by visit type in new datasets 
addmargins(table(AD_Active$PSpecies_Detail, AD_Active$Pm_Species, useNA = "always"))
addmargins(table(AD_Active$PSpecies_Detail, AD_Active$Po_Species, useNA = "always"))
addmargins(table(AD_Active$PSpecies_Detail, AD_Active$Pf_Species, useNA = "always"))


addmargins(table(AD_Passive$PSpecies_Detail, AD_Passive$Pm_Species, useNA = "always"))
addmargins(table(AD_Passive$PSpecies_Detail, AD_Passive$Po_Species, useNA = "always"))
addmargins(table(AD_Passive$PSpecies_Detail, AD_Passive$Pf_Species, useNA = "always"))






##########################################################

# ANALYSIS DATASET & STUDY POPULATION EXPLORATION

##########################################################



##############################################################################################################################################################################
# ASSESSING DIFFERENCES BETWEEN ACTIVE AND PASSIVE COHORTS TO DETERMINE VALIDITY OF INFERENCES FROM PASSIVE POP, & COMPARING TO ACTIVE POP.   
##############################################################################################################################################################################


#________________________________________________________________________________________________________________________________________________________________________________

# 1. Comparing demographics among the active vs. passive populations 

# Notes: 
# - Create comparison between Active and Passive, and those who did vs. did not ever have a Passive Visit
# - Differences in those who did vs. did not ever have a Passive Visit might necessitate IPSW to make more similar when comparing the active and passive populations
# - Limitation to note- we are unable to distinguish in our dataset whether participants did not ever visit a clinic because they did not have malaria symptoms, or alternatively, if they chose not to visit clinics for alternative reasons

#________________________________________________________________________________________________________________________________________________________________________________

##Comparing demographics at the visit level, stratified by visit type: 
addmargins(table(AD_Total$Visit_Type1, AD_Total$Visit_Type2, useNA = "always"))

Compare_VNP_Visits_ActvsPass<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                             strata = "Visit_Type1", 
                                             data = AD_Total, 
                                             factorVars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                             smd=TRUE, 
                                             addOverall=TRUE)

P_Compare_VNPpop1<-print(Compare_VNP_Visits_ActvsPass, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), smd=TRUE, showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)
#summary(Compare_VNP_Visits_ActvsPass)



##Comparing demographics at baseline between subjects who did vs. did not have at least 1 VNP visit during the study -- ie subject-level: 

#First flagging subjects at baseline who had at least 1 passive visit during the study
AD_Passive_Sub<-AD_Passive%>%mutate(PassPop_Flg=1)%>%dplyr::distinct(Subject_ID, PassPop_Flg)
AD_Active_BL<-AD_Active%>%filter(Visit_Type2=="BL")

PassivePop_YN<-left_join(AD_Active_BL, AD_Passive_Sub, by = c("Subject_ID"))
PassivePop_YN_Flg<-PassivePop_YN%>%mutate(PassPop_Flg = ifelse(is.na(PassPop_Flg), 0, PassPop_Flg)) 

addmargins(table(PassivePop_YN_Flg$PassPop_Flg, useNA = "always"))

## RESULTS: 
# 1 = subjects visited study clinic at least 1 time during study; 0 = subjects never visited study clinics during follow-up

#  0    1  <NA>  Sum 
#515 1050    0  1565 


# Second, compare baseline demographic characteristics between the 1,050 subjects who did visit clinics, vs. the 515 subjects who never visited a clinic
Compare_VNP_Pop_YesvsNo<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2", "elevation_m",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                        strata = "PassPop_Flg", 
                                        data = PassivePop_YN_Flg, 
                                        factorVars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                        smd=TRUE, 
                                        addOverall=TRUE)

P_Compare_VNPpop2<-print(Compare_VNP_Pop_YesvsNo, exact = c("AgeCat_Visit", "Sex", "Feverpos", "Pf", "Pm", "Po"),
                         nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), 
                         smd=TRUE, 
                         showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)
#summary(Compare_VNP_Pop_YesvsNo)


##########################################################
#OUTPUT TABLE:  SUPPLEMENTAL TABLE  - Comparison of Participant Characteristics between survey population and clinic subpopulation 
##########################################################




# Third, comparing demographics in VNP vists (ie at visit-level) between those that occurred during active follow-up period (last FU3 visit 10-20-16), vs. those that continued past active Follow-up period (Oct. 21 2016 - Dec 31 2017)

# Assessing range of visit dates per visit type (0 = BL; 1 = FU1; 2 = FU2; 3 = FU3)
table(AD_Active$Visit_date, AD_Active$Visit)
table(AD_Passive$Visit_date, AD_Passive$Visit)

# defining new variable to specify whether visit occurred after active FU visits ended ( 1= yes)
AD_Passive2<-AD_Passive%>%mutate(Flg_VNP_afteractive=ifelse(Visit_date<"2016-11-01", 0,
                                                            ifelse(Visit_date>"2016-10-30", 1, NA)))
addmargins(table(AD_Passive2$Visit_date, AD_Passive2$Flg_VNP_afteractive, useNA = "always"))


# Comparing charactieristics at each visit between those that occurred before vs. after active follow-up visits 
Compare_VNP_Visits_BeforevsAfterAct<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "Pf", "Pm", "Po", "anemia_WHO"), 
                                                    strata = "Flg_VNP_afteractive", data = AD_Passive2, factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2", "WealthQuintile", "Feverpos", "Pf", "Pm", "Po", "anemia_WHO"))

P_Compare_VNPpop3<-print(Compare_VNP_Visits_BeforevsAfterAct,  showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


# Fourth, summarizing demographics at baseline between active and passive surveillance populations: 
Compare_VNP_Pop_YesvsNo<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "Pf", "Pm", "Po"), 
                                        strata = "PassPop_Flg", data = PassivePop_YN_Flg, factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2", "WealthQuintile", "Feverpos", "Pf", "Pm", "Po"), includeNA = F, addOverall=TRUE)

PassivePop_YN_Flg_BL<-PassivePop_YN_Flg%>%filter(Visit == 0)
Compare_VNP_Pop_YesvsNo<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                        strata = "PassPop_Flg", 
                                        data = PassivePop_YN_Flg_BL, 
                                        #includeNA = T,
                                        factorVars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))


P_Compare_VNPpop2<-print(Compare_VNP_Pop_YesvsNo, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po", "elevation_m"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)







##############################################################################################################################################################################
# ASSESSING DEMOGRAPHIC AND CLINICAL CHARACTERISTICS OF ACTIVE AND PASSIVE STUDY POPULATIONS;  PRODUCING DESCRIPTIVE TABLES .   
##############################################################################################################################################################################


#________________________________________________________________________________________________________________________________________________________________________________

# 2. Summarizing demographic and clinical characteristics for study cohort populations 

# Notes: 
# - Create descriptive characteristic tables for active survey population, and separately for clinic sub-population 
# - Produce key summary metrics for in-text results in manuscript 
# - Describe participant loss-to-follow-up for Study Participant Flow Diagram Figure 

#________________________________________________________________________________________________________________________________________________________________________________



##Creating a new calendar time variable for passive visit to categorize according to 6-month intervals
#Merge in baseline visit date to passive visits, then calculate difference between baseline and passive visit (in months), and categorize into 6-month intervals

AD_Active_BL2<-AD_Active_BL%>%dplyr::select(Subject_ID, Visit_date)%>%dplyr::rename(Visit_date_BL=Visit_date)
AD_Passive3<-left_join(AD_Passive2, AD_Active_BL2, by = c("Subject_ID"))
AD_Passive4<-AD_Passive3%>%mutate(MonthDiff_VisDate_BL_VNP=((Visit_date-Visit_date_BL)/30.4))%>%
  mutate(Time_VNP_6moninterval_BL=ifelse(MonthDiff_VisDate_BL_VNP<6,  1, 
                                         ifelse(MonthDiff_VisDate_BL_VNP>=6&MonthDiff_VisDate_BL_VNP<12, 2, 
                                                ifelse(MonthDiff_VisDate_BL_VNP>=12&MonthDiff_VisDate_BL_VNP<18, 3,
                                                       ifelse(MonthDiff_VisDate_BL_VNP>=18&MonthDiff_VisDate_BL_VNP<24, 4,
                                                              ifelse(MonthDiff_VisDate_BL_VNP>=24&MonthDiff_VisDate_BL_VNP<30, 5,
                                                                     ifelse(MonthDiff_VisDate_BL_VNP>=30&MonthDiff_VisDate_BL_VNP<36, 6, NA)))))))%>%
  mutate(Time_VNP_12moninterval_BL=ifelse(Time_VNP_6moninterval_BL<3, 1, #visits within 6 and 12 months
                                          ifelse(Time_VNP_6moninterval_BL>2&Time_VNP_6moninterval_BL<5, 2, #visits between 12 and 24 months
                                                 ifelse(Time_VNP_6moninterval_BL>4, 3, NA)))) #visits after 25 months       





#### By Analysis Population - Overall and by Follow-up Visit & Species: 



## Active Population ##


##Look by visit across same species, and between species
##exclude missing visits 
AD_Active_nonmissvisits<-AD_Active%>%filter(!is.na(Visit_date))

No_FU3<-AD_Active%>%filter()



#overall - all species in active follow-up. 
Demog_Active_AllVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                       strata = "Visit", 
                                       data = AD_Active_nonmissvisits, 
                                       includeNA = F,
                                       factorVars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

P_Demog_Active_AllVisits<-print(Demog_Active_AllVisits, exact = c( "season_dry", "Sx_Prior6mo_new", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), showAllLevels = TRUE, is.na=F, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


##########################################################
#OUTPUT TABLE:  SUPPLEMENTAL TABLE  - Description of Total Participant Characteristics Across Active Study Visits DUring Follow-up
##########################################################



#overall - Pf species in active follow-up, by visit 
Demog_Active_PfVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                      strata = "Pf",
                                      # strata = "Pf_Species", 
                                      data = AD_Active_nonmissvisits, 
                                      includeNA = T,
                                      factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Pf_Demog_Active_AllVisits<-print(Demog_Active_PfVisits, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                 #exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo", "AntiMalTx_Prior6mo", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                 showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#overall - Pf species in active follow-up, baseline visit only 
AD_Active_nonmiss_BL<-AD_Active%>%filter(Visit_Type2=="BL")
Demog_Active_Pf_BL<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                   strata = "Pf",
                                   # strata = "Pf_Species", 
                                   data = AD_Active_nonmiss_BL, 
                                   #includeNA = T,
                                   factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Demog_Active_Pf_BL_print<-print(Demog_Active_Pf_BL, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                #exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo", "AntiMalTx_Prior6mo", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                showAllLevels = TRUE, 
                                #Turn on and off the below line to assess missingness by var, vs. calculate p-values excluding missingness.  
                                #is.na=T, 
                                noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)





#overall - Pm species in active follow-up, by visit 
Demog_Active_PmVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                      strata = "Pm", 
                                      data = AD_Active_nonmissvisits, 
                                      includeNA = T,
                                      factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Pm_Demog_Active_AllVisits<-print(Demog_Active_PmVisits, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#overall - Pm species in active follow-up, baseline visit only  
Demog_Active_Pm_BL<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                   strata = "Pm",
                                   # strata = "Pf_Species", 
                                   data = AD_Active_nonmiss_BL, 
                                   #includeNA = T,
                                   factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Demog_Active_Pm_BL_print<-print(Demog_Active_Pm_BL, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                #exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo", "AntiMalTx_Prior6mo", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                showAllLevels = TRUE, 
                                #Turn on and off the below line to assess missingness by var, vs. calculate p-values excluding missingness.  
                                #is.na=T, 
                                noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)




#overall - Po species in active follow-up, by visit 
Demog_Active_PoVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                      strata = "Po_Species", 
                                      data = AD_Active_nonmissvisits, 
                                      includeNA = T,
                                      factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Po_Demog_Active_AllVisits<-print(Demog_Active_PoVisits, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#overall - Po species in active follow-up, baseline visit only 
Demog_Active_Po_BL<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                   strata = "Po",
                                   # strata = "Po_Species", 
                                   data = AD_Active_nonmiss_BL, 
                                   includeNA = T,
                                   factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Demog_Active_Po_BL_print<-print(Demog_Active_Po_BL, nonnormal = c("Pd_Po", "Pd_Pm", "Pd_Po"), 
                                #exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo", "AntiMalTx_Prior6mo", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                showAllLevels = TRUE, 
                                #Turn on and off the below line to assess missingness by var, vs. calculate p-values excluding missingness.  
                                is.na=T, 
                                noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)






## Passive Population ##


# Looking at number of visits and number of subjects across 1-year, 2-years, and 2+ years from baseline 

# visits by timepoint from baseline 
addmargins(table(AD_Passive4$Time_VNP_6moninterval_BL, useNA = "always"))
addmargins(table(AD_Passive4$Time_VNP_12moninterval_BL, useNA = "always"))
#addmargins(table(AD_Passive4$Visit_date, useNA = "always"))
#addmargins(table(AD_Passive4$Visit_date_BL, useNA = "always"))


# subjects by timepoint from baseline 
AD_Passive4_subjectlvl_6mon<-AD_Passive4%>%
  group_by(Time_VNP_6moninterval_BL)%>%
  distinct(Subject_ID)
addmargins(table(AD_Passive4_subjectlvl_6mon$Time_VNP_6moninterval_BL))



AD_Passive4_subjectlvl_year<-AD_Passive4%>%
  group_by(Time_VNP_12moninterval_BL)%>%
  distinct(Subject_ID)
addmargins(table(AD_Passive4_subjectlvl_year$Time_VNP_12moninterval_BL))


##########################################################
#OUTPUT RESULTS:  SUPPLEMENTAL TABLE  - No. of subjects with clinic visits across time points from baseline 
##########################################################


#overall - all species in passive follow-up. 
Demog_Passive_AllVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "anemia_WHO", "season_dry",  "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                        strata = "Time_VNP_12moninterval_BL", 
                                        data = AD_Passive4, 
                                        includeNA = T,
                                        factorVars = c("AgeCat_Visit", "Sex", "Village", "HealthArea", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "anemia_WHO", "season_dry", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

P_Demog_Passive_AllVisits<-print(Demog_Passive_AllVisits, exact = c( "season_dry", "Feverpos", "RDTpos",  "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), showAllLevels = TRUE, is.na=F, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#summary(Demog_Passive_AllVisits)


##########################################################
#OUTPUT TABLE:  SUPPLEMENTAL TABLE  - Participant Characteristics Across Follow-up (Clinic Subpopulation)
##########################################################



# Median (IQR) number of passive visits by subject and by household -- in-text numbers for manuscript

AD_Passive_Subvisitcount<-AD_Passive4%>%
  group_by(Subject_ID) %>%
  dplyr::summarise(count=n())
summary(AD_Passive_Subvisitcount$count)

#RESULTS FOR MEDIAN # VNP VISITS PER SUBJECT IN PASSIVE POP (n=1050)
#Min.   1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   3.245   4.000  19.000 

AD_Passive_HHvisitcount<-AD_Passive4%>%
  group_by(Household) %>%
  dplyr::summarise(count=n())
summary(AD_Passive_HHvisitcount$count)

#RESULTS FOR MEDIAN # VNP VISITS PER HOUSEHOLD IN PASSIVE POP (n=1050)
#Min.  st Qu. Median  Mean    3rd Qu.    Max. 
#1.00  5.00   11.00   15.63   21.75     74.00 

##Looking at which household had 74 VNP visits (max)
print(AD_Passive_HHvisitcount$Household[AD_Passive_HHvisitcount$count==74])




######
#Checking distribution of visits over time 

#Active Pop
ggplot(AD_Active_nonmissvisits, aes(Visit_date))+
  geom_histogram()


ggplot(AD_Active_nonmissvisits, aes(Visit))+
  geom_histogram()

#Active Pop by Health Area
ggplot(AD_Active_nonmissvisits, aes(Visit))+
  geom_histogram()+
  facet_grid(rows = vars(Village))


AD_Active_missvisits<- AD_Active%>%filter(is.na(Visit_date))




#VNP Pop
AD_Passive4<-apply_labels(AD_Passive4, 
                          village=c("Bu" = 1, "Impuru" = 2, "Pema" = 3, "Kimpoko" = 4, "Ngamanzo" = 5, "Iye" = 6, "Voix du Peuple" = 7))

ggplot(AD_Passive4, aes(Time_VNP_6moninterval_BL))+
  geom_histogram()

ggplot(AD_Passive4, aes(Time_VNP_12moninterval_BL))+
  geom_histogram()

ggplot(AD_Passive4, aes(Visit_date))+
  geom_histogram()+
  facet_grid(rows = vars(Village))



##Look by visit across same species, and between species
##exclude missing visits 
PassivePop_Y_nonmissvisits<-PassivePop_YN_Flg%>%filter(PassPop_Flg==1)

#overall - all species in passive population (baseline visit only)
Demog_Passive_BL<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                 data = PassivePop_Y_nonmissvisits, 
                                 includeNA = T,
                                 factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

P_Demog_Passive_BL<-print(Demog_Passive_BL, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#overall - Pf species in active population 
Demog_Active_PfVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                      strata = "Pf",
                                      # strata = "Pf_Species", 
                                      data = PassivePop_Y_nonmissvisits, 
                                      includeNA = T,
                                      factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))


Pf_Demog_Active_AllVisits<-print(Demog_Active_PfVisits, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                 #exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo", "AntiMalTx_Prior6mo", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), 
                                 showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#overall - Pm species in active follow-up, by visit 
Demog_Active_PmVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                      strata = "Pm", 
                                      data = PassivePop_Y_nonmissvisits, 
                                      includeNA = T,
                                      factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Pm_Demog_Active_AllVisits<-print(Demog_Active_PmVisits, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)



#overall - Po species in active follow-up, by visit 
Demog_Active_PoVisits<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species", "Pd_Pf", "Pd_Pm", "Pd_Po"), 
                                      strata = "Po", 
                                      data = PassivePop_Y_nonmissvisits, 
                                      includeNA = T,
                                      factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"))

Po_Demog_Active_AllVisits<-print(Demog_Active_PoVisits, nonnormal = c("Pd_Pf", "Pd_Pm", "Pd_Po"), exact = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "RDTpos", "BedNetPrior", "season_dry", "Sx_Prior6mo_new", "AntiMalTx_Prior6mo_new", "Pf", "Pm", "Po", "Pf_Species", "Pm_Species", "Po_Species"), showAllLevels = TRUE, is.na=T, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)










##########################################################

# ARCHIVE STUDY ANALYSIS DATASETS FOR ANALYTIC USE IN SEPARATE R PROGRAMS 

##########################################################

write.xlsx(AD_Total, '.../AD_Total_Jul22.xlsx', colNames = TRUE)
write.xlsx(AD_Active, '.../AD_Active_Jul22.xlsx', colNames = TRUE)
write.xlsx(AD_Passive, '.../AD_Passive_Jul22.xlsx', colNames = TRUE)









