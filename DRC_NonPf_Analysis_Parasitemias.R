##########################################
#Programmer: `Rachel Sendor
#Last update: Mar 2023
#Purpose:     DRC Non-falciparum Descriptive Epidemiology Analysis
#             Parasite Density Analyses
#             
##########################################


# SET WORKING DIRECTORY
#### xxxxxxxxxxx

#Loading R Packages:
library(tidyverse) 
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggbreak)
library(ggthemes)



##########################################################
# SET WORKING DIRECTORY


##########################################################
# IMPORTING DATASETS FOR ANALYSIS 

## Dataset contains 1 row per participant per study visit, separated by total, active (survey), and passive (clinic) populations

AD_Total <- read_excel("../AD_Total_Jul22.xlsx")
View(AD_Total)

AD_Active <- read_excel("../AD_Active_Jul22.xlsx")
View(AD_Active)

AD_Passive <- read_excel("../AD_Passive_Jul22.xlsx")
View(AD_Passive)



################################################################################################

# STUDY ANALYSES: Assessing Parasitemia across species and comparing between species 

################################################################################################



##First, need to filter out the rehydrated samples and the missing visit dates to exclude those lost to FU  



##UPDATE: over-ride variable by multiplying all non-Pf PD's by 4 because of dilution factor during PCR assay - need to correct the Pd for non-pf only so it is more valid comparison with Pf Pd (since Pf was not diluted during PCR)

AD_Total<-AD_Total%>%mutate(Pd_Po_mult4 = Pd_Po*4)%>%
                     mutate(Pd_Pm_mult4 = Pd_Pm*4)


#Perform a quick comparison of rehydrated vs. non-rehydrated samples by species and demographics to see if there are important differences.  
#Looking in Total Pop.
AD_Total_nonmissvisits<-AD_Total%>%filter(!is.na(Visit_date))

addmargins(table(AD_Total_nonmissvisits$Rehydrated, AD_Total_nonmissvisits$Visit_Type2, useNA = "always"))

Compare_RehydratedSamples_YesvsNo<-CreateTableOne(vars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2",  "WealthQuintile", "Feverpos", "Pf", "Pm", "Po"), 
                                                  strata = "Rehydrated", data = AD_Total_nonmissvisits, factorVars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2", "WealthQuintile", "Feverpos", "Pf", "Pm", "Po"))

P_Compare_Rehydrate<-print(Compare_RehydratedSamples_YesvsNo, factorvars = c("AgeCat_Visit", "Sex", "Village", "Visit_Type2", "Wealthquintile", "Feverpos", "Pf", "Pm", "Po"), exact  = c("AgeCat_Visit", "Sex", "Visit_Type2", "Feverpos", "Pf", "Pm", "Po"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#Result: In Total Pop, the % of visits with samples that were rehydrated were higher for Bu, Baseline and VNP visit types, and Fever =N. 




##Now exclude rehydrated samples from parasitemia estimation since the concentration has been affected
#Also subset data for individual species' parasite density graphs
po_pd<-AD_Total_nonmissvisits%>%filter(Po==1&Rehydrated!=1)
pm_pd<-AD_Total_nonmissvisits%>%filter(Pm==1&Rehydrated!=1)
pm_pd<-pm_pd%>%filter(Pm_Species==as.factor('Pm Mono') | Pm_Species==as.factor('Pm Mixed'))
pf_pd<-AD_Total_nonmissvisits%>%filter(Pf==1&Rehydrated!=1)
pf_pd2<-pf_pd%>%filter(Pf_Species==as.factor('Pf Mono') | Pf_Species==as.factor('Pf Mixed'))
pall_pd<-AD_Total_nonmissvisits%>%filter(Rehydrated!=1)


##Looking at excluded due to rehydrated: 
po_pd_rehyd<-AD_Total_nonmissvisits%>%filter(Po==1&Rehydrated==1)
pm_pd_rehyd<-AD_Total_nonmissvisits%>%filter(Pm==1&Rehydrated==1)
pf_pd_rehyd<-AD_Total_nonmissvisits%>%filter(Pf==1&Rehydrated==1)




#####TOTAL POPULATION -- ALL VISITS NON-Rehydrated 

#####Basic Descriptive Stats (Median because non-normally distributed parasitemia)


#Po
po_pd %>% 
  dplyr::summarize(n(), median(Pd_Po_mult4), quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))



##### IN-TEXT RESULTS FOR MANUSCRIPT:   
#  `n()` `median(Pd_Po_mult4)` `quantile(Pd_Po_mult4, c(0.25))` `quantile(Pd_Po_mult4, c(0.75))` `min(Pd_Po_mult4)` `max(Pd_Po_mult4)`
#   164                  10.2                             2.73                             47.4              0.168            106476.



##rehydrated summary: 
po_pd_rehyd%>%group_by(Po_Species)%>%dplyr::summarise(n())


#mixed vs. mono summary: 
po_pd_mixed<-po_pd%>%filter(Po_Species=="Po Mixed")%>%dplyr::summarize(n(), median(Pd_Po_mult4),  quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))
po_pd_mono<-po_pd%>%filter(Po_Species=="Po Mono")%>%dplyr::summarize(n(), median(Pd_Po_mult4),  quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))
rbind(po_pd_mixed, po_pd_mono)

Compare_Pd_Po_mult4Species<-CreateTableOne(vars = c("Pd_Po_mult4"), strata = "Po_Species", data = po_pd)
print(Compare_Pd_Po_mult4Species, exact = c("Pd_Po_mult4"), nonnormal=c("Pd_Po_mult4"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


wilcox.test(Pd_Po_mult4 ~ Po_Species, data = po_pd)



#Pm
pm_pd %>% 
  dplyr::summarize(n(), median(Pd_Pm_mult4), quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))

##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#`n()`      `median(Pd_Pm_mult4)` `quantile(Pd_Pm_mult4, c(0.25))` `quantile(Pd_Pm_mult4, c(0.75))`    `min(Pd_Pm_mult4)` `max(Pd_Pm_mult4)`
# 307                  25.7                             7.69                             119.              0.618            358246.



##rehydrated summary: 
pm_pd_rehyd%>%group_by(Pm_Species)%>%dplyr::summarise(n())


pm_pd_mono<-pm_pd%>%filter(Pm_Species=="Pm Mono")%>%dplyr::summarize(n(), median(Pd_Pm_mult4), quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))
pm_pd_mixed<-pm_pd%>%filter(Pm_Species=="Pm Mixed")%>%dplyr::summarize(n(), median(Pd_Pm_mult4),quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))
rbind(pm_pd_mixed, pm_pd_mono)

Compare_Pd_Pm_mult4Species<-CreateTableOne(vars = c("Pd_Pm_mult4"), strata = "Pm_Species", data = pm_pd)
print(Compare_Pd_Pm_mult4Species, exact = c("Pd_Pm_mult4"), nonnormal=c("Pd_Pm_mult4"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

## COMPUTING MANUALLY TO CONFIRM AND OUTPUT TEST STATISTIC FOR NATURE REUQIREMENTS 
wilcox.test(Pd_Pm_mult4 ~ Pm_Species, data = pm_pd)

#Wilcoxon rank sum test with continuity correction

#data:  Pd_Pm_mult4 by Pm_Species
#W = 8486, p-value = 0.07
#alternative hypothesis: true location shift is not equal to 0



#Pf
pf_pd %>% 
  dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#`n()` `median(Pd_Pf)` `quantile(Pd_Pf, c(0.25))` `quantile(Pd_Pf, c(0.75))` `min(Pd_Pf)` `max(Pd_Pf)`

# 3730            267.                       18.8                      4526.        0.589      1165100



pf_pd_mono<-pf_pd%>%filter(Pf_Species=="Pf Mono")%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))
pf_pd_mixed<-pf_pd%>%filter(Pf_Species=="Pf Mixed")%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))
rbind( pf_pd_mixed, pf_pd_mono)

##rehydrated summary: 
pf_pd_rehyd%>%group_by(Pf_Species)%>%dplyr::summarise(n())


Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = pf_pd)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)







###comparing Pd's in Pf mono vs. Pf+Pm mixed, and Pf mono vs. Pf+Po mixed: 

#Pf mono vs. Pf+Po   (no Pm infections)
pf_pd_mono<-pf_pd%>%filter(Pf_Species=="Pf Mono")
pf_pd_mixed_Po<-pf_pd%>%filter(Pf_Species=="Pf Mixed"&Po==1&Pm!=1)
PfvsPfPo<- rbind( pf_pd_mixed_Po, pf_pd_mono)

PfvsPfPo%>%group_by(Pf_Species)%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = PfvsPfPo)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#RESULT - Pf mono vs. Pf+Po mixed -- TOTAL POP. 
#Stratified by Pf_Species
#                       level Pf Mixed                Pf Mono                 p     test   
#n                          96                      3418                                 
#Pd_Pf (median [IQR])       324.98 [47.52, 4103.50] 279.60 [17.38, 5014.62] 0.468 nonnorm




#Pf mono vs. Pf+Pm  (no Po infections)
pf_pd_mono<-pf_pd%>%filter(Pf_Species=="Pf Mono")
pf_pd_mixed_Pm<-pf_pd%>%filter(Pf_Species=="Pf Mixed"&Pm==1&Po!=1)
PfvsPfPm<- rbind( pf_pd_mixed_Pm, pf_pd_mono)

PfvsPfPm%>%group_by(Pf_Species)%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = PfvsPfPm)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#RESULT - Pf mono vs. Pf+Pm mixed -- TOTAL POP. 
#Stratified by Pf_Species
#                   level Pf Mixed                Pf Mono                 p     test   
#n                          206                     3418                                 
#Pd_Pf (median [IQR])       166.50 [33.39, 1041.74] 279.60 [17.38, 5014.62] 0.059 nonnorm







######### NOW REPEAT FOR PASSIVE AND ACTIVE POPULATIONS 

#PASSIVE POPULATION  VISITS NON-Rehydrated 

##Basic Descriptive Stats (Median because non-normally distributed parasitemia)

### filter to only passive pop
po_pd_passive<-po_pd%>%filter(Visit_Type2=="VNP")
po_pd_rehyd_passive<-po_pd_rehyd%>%filter(Visit_Type2=="VNP")

#Po
po_pd_passive %>% 
  dplyr::summarize(n(), median(Pd_Po_mult4), quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))


##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#`n()` `median(Pd_Po_mult4)` `quantile(Pd_Po_mult4, c(0.25))` `quantile(Pd_Po_mult4, c(0.75))` `min(Pd_Po_mult4)` `max(Pd_Po_mult4)`
# 93                  17.7                             4.59                             65.8              0.168              2875.


##rehydrated summary: 
po_pd_rehyd_passive%>%group_by(Po_Species)%>%dplyr::summarise(n())

#mixed vs. mono summary: 
po_pd_mixed<-po_pd_passive%>%filter(Po_Species=="Po Mixed")%>%dplyr::summarize(n(), median(Pd_Po_mult4),  quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))
po_pd_mono<-po_pd_passive%>%filter(Po_Species=="Po Mono")%>%dplyr::summarize(n(), median(Pd_Po_mult4),  quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))
rbind(po_pd_mixed, po_pd_mono)

Compare_Pd_Po_mult4Species<-CreateTableOne(vars = c("Pd_Po_mult4"), strata = "Po_Species", data = po_pd_passive)
print(Compare_Pd_Po_mult4Species, exact = c("Pd_Po_mult4"), nonnormal=c("Pd_Po_mult4"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#Pm



### filter to only passive pop
pm_pd_passive<-pm_pd%>%filter(Visit_Type2=="VNP")
pm_pd_rehyd_passive<-pm_pd_rehyd%>%filter(Visit_Type2=="VNP")    

pm_pd_passive %>% 
  dplyr::summarize(n(), median(Pd_Pm_mult4), quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))

##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#`n()` `median(Pd_Pm_mult4)` `quantile(Pd_Pm_mult4, c(0.25))` `quantile(Pd_Pm_mult4, c(0.75))` `min(Pd_Pm_mult4)` `max(Pd_Pm_mult4)`
# 132                  36.5                             4.68                             182.              0.618             23288.


##rehydrated summary: 
pm_pd_rehyd_passive%>%group_by(Pm_Species)%>%dplyr::summarise(n())


pm_pd_mono<-pm_pd_passive%>%filter(Pm_Species=="Pm Mono")%>%dplyr::summarize(n(), median(Pd_Pm_mult4), quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))
pm_pd_mixed<-pm_pd_passive%>%filter(Pm_Species=="Pm Mixed")%>%dplyr::summarize(n(), median(Pd_Pm_mult4),quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))
rbind(pm_pd_mixed, pm_pd_mono)

Compare_Pd_Pm_mult4Species<-CreateTableOne(vars = c("Pd_Pm_mult4"), strata = "Pm_Species", data = pm_pd_passive)
print(Compare_Pd_Pm_mult4Species, exact = c("Pd_Pm_mult4"), nonnormal=c("Pd_Pm_mult4"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#Pf

### filter to only passive pop
pf_pd_passive<-pf_pd%>%filter(Visit_Type2=="VNP")
pf_pd_rehyd_passive<-pf_pd_rehyd%>%filter(Visit_Type2=="VNP")    


pf_pd_passive %>% 
  dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#  `n()` `median(Pd_Pf)` `quantile(Pd_Pf, c(0.25))` `quantile(Pd_Pf, c(0.75))` `min(Pd_Pf)` `max(Pd_Pf)`
#  1970           2644.                       113.                     16931.        0.589      1165100


pf_pd_mono<-pf_pd_passive%>%filter(Pf_Species=="Pf Mono")%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))
pf_pd_mixed<-pf_pd_passive%>%filter(Pf_Species=="Pf Mixed")%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))
rbind( pf_pd_mixed, pf_pd_mono)


Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = pf_pd_passive)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


##rehydrated summary: 
pf_pd_rehyd_passive%>%group_by(Pf_Species)%>%dplyr::summarise(n())

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = pf_pd_rehyd_passive)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)





###comparing Pd's in Pf mono vs. Pf+Pm mixed, and Pf mono vs. Pf+Po mixed: 



#Pf mono vs. Pf+Po   (no Pm infections)
pf_pd_mono<-pf_pd_passive%>%filter(Pf_Species=="Pf Mono")
pf_pd_mixed_Po<-pf_pd_passive%>%filter(Pf_Species=="Pf Mixed"&Po==1&Pm!=1)
PfvsPfPo_passive<- rbind( pf_pd_mixed_Po, pf_pd_mono)

PfvsPfPo_passive%>%group_by(Pf_Species)%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = PfvsPfPo_passive)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#RESULT - Pf mono vs. Pf+Po mixed -- PASSIVE POP. 
#Stratified by Pf_Species
#                     level Pf Mixed                   Pf Mono                    p     test   
#n                          47                         1834                                    
#Pd_Pf (median [IQR])       1710.00 [165.18, 13082.50] 2897.00 [126.51, 18058.25] 0.596 nonnorm


#Pf mono vs. Pf+Pm   (no Po infections)
pf_pd_mono<-pf_pd_passive%>%filter(Pf_Species=="Pf Mono")
pf_pd_mixed_Pm<-pf_pd_passive%>%filter(Pf_Species=="Pf Mixed"&Pm==1&Po!=1)
PfvsPfPm_passive<- rbind( pf_pd_mixed_Pm, pf_pd_mono)

PfvsPfPm_passive%>%group_by(Pf_Species)%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = PfvsPfPm_passive)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#RESULT - Pf mono vs. Pf+Pm mixed -- PASSIVE POP. 
# Stratified by Pf_Species
#                     level Pf Mixed                Pf Mono                    p      test   
#n                          84                      1834                                     
#Pd_Pf (median [IQR])       287.78 [56.96, 4318.60] 2897.00 [126.51, 18058.25] <0.001 nonnorm






##### ACTIVE POPULATION VISITS NON-Rehydrated 

#####Basic Descriptive Stats (Median because non-normally distributed parasitemia)

### filter to only active pop
po_pd_active<-po_pd%>%filter(Visit_Type1=="Act")
po_pd_rehyd_active<-po_pd_rehyd%>%filter(Visit_Type1=="Act")

#Po
po_pd_active %>% 
  dplyr::summarize(n(), median(Pd_Po_mult4), quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))

##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#`n()` `median(Pd_Po_mult4)` `quantile(Pd_Po_mult4, c(0.25))` `quantile(Pd_Po_mult4, c(0.75))` `min(Pd_Po_mult4)` `max(Pd_Po_mult4)`
# 71                  5.77                             1.97                             28.0              0.647            106476.

##rehydrated summary: 
po_pd_rehyd_active%>%group_by(Po_Species)%>%dplyr::summarise(n())

#mixed vs. mono summary: 
po_pd_mixed<-po_pd_active%>%filter(Po_Species=="Po Mixed")%>%dplyr::summarize(n(), median(Pd_Po_mult4),  quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))
po_pd_mono<-po_pd_active%>%filter(Po_Species=="Po Mono")%>%dplyr::summarize(n(), median(Pd_Po_mult4),  quantile(Pd_Po_mult4, c(0.25)), quantile(Pd_Po_mult4, c(0.75)), min(Pd_Po_mult4), max(Pd_Po_mult4))
rbind(po_pd_mixed, po_pd_mono)

Compare_Pd_Po_mult4Species<-CreateTableOne(vars = c("Pd_Po_mult4"), strata = "Po_Species", data = po_pd_active)
print(Compare_Pd_Po_mult4Species, exact = c("Pd_Po_mult4"), nonnormal=c("Pd_Po_mult4"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#Pm

### filter to only active pop
pm_pd_active<-pm_pd%>%filter(Visit_Type1=="Act")
pm_pd_rehyd_active<-pm_pd_rehyd%>%filter(Visit_Type1=="Act")    

pm_pd_active %>% 
  dplyr::summarize(n(), median(Pd_Pm_mult4), quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))

##### IN-TEXT RESULTS FOR MANUSCRIPT: 
#`n()` `median(Pd_Pm_mult4)` `quantile(Pd_Pm_mult4, c(0.25))` `quantile(Pd_Pm_mult4, c(0.75))` `min(Pd_Pm_mult4)` `max(Pd_Pm_mult4)`
# 175                  22.4                             8.51                             72.2              0.778            358246.


##rehydrated summary: 
pm_pd_rehyd_active%>%group_by(Pm_Species)%>%dplyr::summarise(n())


pm_pd_mono<-pm_pd_active%>%filter(Pm_Species=="Pm Mono")%>%dplyr::summarize(n(), median(Pd_Pm_mult4), quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))
pm_pd_mixed<-pm_pd_active%>%filter(Pm_Species=="Pm Mixed")%>%dplyr::summarize(n(), median(Pd_Pm_mult4),quantile(Pd_Pm_mult4, c(0.25)), quantile(Pd_Pm_mult4, c(0.75)), min(Pd_Pm_mult4), max(Pd_Pm_mult4))
rbind(pm_pd_mixed, pm_pd_mono)

Compare_Pd_Pm_mult4Species<-CreateTableOne(vars = c("Pd_Pm_mult4"), strata = "Pm_Species", data = pm_pd_active)
print(Compare_Pd_Pm_mult4Species, exact = c("Pd_Pm_mult4"), nonnormal=c("Pd_Pm_mult4"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


#Pf

### filter to only active pop
pf_pd_active<-pf_pd%>%filter(Visit_Type1=="Act")
pf_pd_rehyd_active<-pf_pd_rehyd%>%filter(Visit_Type1=="Act")    


pf_pd_active %>% 
  dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

##### IN-TEXT RESULTS FOR MANUSCRIPT:
#`n()` `median(Pd_Pf)` `quantile(Pd_Pf, c(0.25))` `quantile(Pd_Pf, c(0.75))` `min(Pd_Pf)` `max(Pd_Pf)`
#1760            52.6                       8.18                       343.        0.612       268250


pf_pd_mono<-pf_pd_active%>%filter(Pf_Species=="Pf Mono")%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))
pf_pd_mixed<-pf_pd_active%>%filter(Pf_Species=="Pf Mixed")%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))
rbind( pf_pd_mixed, pf_pd_mono)

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = pf_pd_active)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)


##rehydrated summary: 
pf_pd_rehyd_active%>%group_by(Pf_Species)%>%dplyr::summarise(n())


Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = pf_pd_rehyd_active)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)





###comparing Pd's in Pf mono vs. Pf+Pm mixed, and Pf mono vs. Pf+Po mixed: 

#Pf mono vs. Pf+Po   (no Pm infections)
pf_pd_mono<-pf_pd_active%>%filter(Pf_Species=="Pf Mono")
pf_pd_mixed_Po<-pf_pd_active%>%filter(Pf_Species=="Pf Mixed"&Po==1&Pm!=1)
PfvsPfPo_active<- rbind( pf_pd_mixed_Po, pf_pd_mono)

PfvsPfPo_active%>%group_by(Pf_Species)%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = PfvsPfPo_active)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#RESULT - Pf mono vs. Pf+Po mixed -- ACTIVE POP. 
#Stratified by Pf_Species
#level                      Pf Mixed               Pf Mono              p     test   
#n                          49                     1584                              
#Pd_Pf (median [IQR])       122.11 [42.16, 485.90] 46.95 [7.40, 311.12] 0.010 nonnorm


#Pf mono vs. Pf+Pm   (no Po infections)
pf_pd_mono<-pf_pd_active%>%filter(Pf_Species=="Pf Mono")
pf_pd_mixed_Pm<-pf_pd_active%>%filter(Pf_Species=="Pf Mixed"&Pm==1&Po!=1)
PfvsPfPm_active<- rbind( pf_pd_mixed_Pm, pf_pd_mono)

PfvsPfPm_active%>%group_by(Pf_Species)%>%dplyr::summarize(n(), median(Pd_Pf), quantile(Pd_Pf, c(0.25)), quantile(Pd_Pf, c(0.75)), min(Pd_Pf), max(Pd_Pf))

Compare_Pd_PfSpecies<-CreateTableOne(vars = c("Pd_Pf"), strata = "Pf_Species", data = PfvsPfPm_active)
print(Compare_Pd_PfSpecies, exact = c("Pd_Pf"), nonnormal=c("Pd_Pf"), showAllLevels = TRUE, noSpaces = T)%>%write.table("clipboard", sep="\t", row.names = T)

#RESULT - Pf mono vs. Pf+Pm mixed -- ACTIVE POP. 
#Stratified by Pf_Species
#level                      Pf Mixed               Pf Mono              p     test   
#n                          122                    1584                              
#Pd_Pf (median [IQR])       126.88 [24.94, 496.40] 46.95 [7.40, 311.12] 0.001 nonnorm






##########################

#CREATING PARASITEMIA FIGURES FOR MANUSCRIPT 


#1. Density plots comparing Pf parasitemia in Pf mono vs. Pf+Po  and Pf mono vs. Pf+Pm  (with background histogram for raw counts)

## SEPARATELY BY VISIT TYPE (ACTIVE VS. PASSIVE POP.)   ALSO CREATE FOR TOTAL POP BUT DON'T INCLUDE IN MANUSCRIPT 

xlabels_manual<-c("0.01", "0.1", "1", "10", "100", "1,000", "10,000", "100,000", "1,000,000") 

##ACTIVE POP.

##Pf vs. Pm 

#Density plot using log scale 
dplot_PfbyPm<-ggplot(PfvsPfPm_active)+
  geom_density(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.5)+
  scale_x_log10(name="Estimated parasite density (p/uL)")
plot(dplot_PfbyPm)

barplot_PfPm<-ggplot(PfvsPfPm_active)+
  geom_bar(aes(x=Pd_Pf, y=(count)/sum(count), fill=Pf_Species))+
  scale_x_log10()
plot(barplot_PfPm)


#histogram plot using log scale
hplot_PfbyPm<-ggplot(PfvsPfPm_active)+
  geom_histogram(aes(x = Pd_Pf, fill = Pf_Species), position="dodge", alpha=0.3)+
  scale_x_log10(name="Estimated parasite density (p/uL)")
plot(hplot_PfbyPm)

hplot_PfbyPo<-ggplot(PfvsPfPo_active)+
  geom_histogram(aes(x = Pd_Pf, fill = Pf_Species), position="dodge", alpha=0.3)+
  scale_x_log10(name="Estimated parasite density (p/uL)")
plot(hplot_PfbyPo)

#overlaying density plot onto histogram to show the raw counts too
#un-stratified counts in background: 
PfbyPm_overlay<-ggplot(PfvsPfPm_active, aes(x=Pd_Pf, bins=10))+
  # geom_histogram(aes(y=..density..), alpha=0.5, col='gray', lwd=0.4)+
  # scale_x_log10(name="Estimated Pf Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  scale_x_log10(name="", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.4)
#labs(title="Pf mono vs. Pf+Pm mixed infections - Active visits", fontsize = 8) 
PfbyPm_overlay_fig <- PfbyPm_overlay + guides(fill=guide_legend(title=""))
plot(PfbyPm_overlay_fig)


#Stratified (but stacked) counts in background:
PfbyPm_overlay<-ggplot(PfvsPfPm_active, aes(x=Pd_Pf, fill=Pf_Species, bins=10))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.3, col='gray', lwd=0.2)+
  #scale_x_log10(name="Estimated Pf Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  scale_x_log10(name="", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.2)+
  scale_fill_manual(values=c("#5b9bf5", "#f5cc5b"))
#+labs(title="Pf mono vs. Pf+Pm mixed infections - Survey Pop.", fontsize = 8)

PfbyPm_overlay_fig <- PfbyPm_overlay + guides(fill=guide_legend(title="")) +
  # geom_vline(xintercept = 10, linetype="solid", color = "dark gray", size=1.3) +
  # annotate("text", label = "LoD", x = 7, y = 0.87, size = 5, colour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text=element_text(size=14), legend.position = c(0.79, 0.65))

PfbyPm_overlay_fig2 <- PfbyPm_overlay_fig + expand_limits(x=c(0.1,3000000), y=c(0.0, 0.80)) 
plot(PfbyPm_overlay_fig2)



#####Creating a scaled plot to depict differences in sample size
ggplot(PfvsPfPm_active) +
  geom_density(aes(x = Pd_Pf, y = after_stat(count), fill=Pf_Species), alpha = 0.8)+
  scale_x_log10(name="Estimated Pf parasite density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)



## Pf vs. Po

#Density plot on log scale
dplot_PfbyPo<-ggplot(PfvsPfPo_active)+
  geom_density(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.5)+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")
plot(dplot_PfbyPo)

#histogram plot using log scale
hplot_PfbyPo<-ggplot(PfvsPfPo_active)+
  geom_histogram(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.3)+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")
plot(hplot_PfbyPo)

#overlaying density plot onto histogram to show the raw counts too
#un-stratified counts in background: 
PfbyPo_overlay<-ggplot(PfvsPfPo_active, aes(x=Pd_Pf, bins=10))+
  geom_histogram(aes(y=..density..), alpha=0.5, col='gray', lwd=0.4)+
  scale_x_log10(name="Estimated Pf parasite density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.4)
#labs(title="Pf mono vs. Pf+Po mixed infections - Active visits", fontsize = 8)

PfbyPo_overlay_fig <- PfbyPo_overlay + guides(fill=guide_legend(title=""))
plot(PfbyPo_overlay_fig)

#Stratified (but stacked) counts in background:
PfbyPo_overlay<-ggplot(PfvsPfPo_active, aes(x=Pd_Pf, fill=Pf_Species, bins=10))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.3, col='gray', lwd=0.2)+
  #scale_x_log10(name="Estimated Pf Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  scale_x_log10(name="", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.2)+
  scale_fill_manual(values=c("#c44343", "#f5cc5b"))
#+labs(title="Pf mono vs. Pf+Po mixed infections - Survey Pop.", fontsize = 8)

PfbyPo_overlay_fig <- PfbyPo_overlay + guides(fill=guide_legend(title="")) +
  # geom_vline(xintercept = 10, linetype="solid", color = "dark gray", size=1.3) +
  # annotate("text", label = "LoD", x = 7, y = 0.87, size = 5, colour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text=element_text(size=14), legend.position = c(0.79, 0.65))

PfbyPo_overlay_fig2 <- PfbyPo_overlay_fig + expand_limits(x=c(0.1,3000000), y=c(0.0, 0.80))
plot(PfbyPo_overlay_fig2)




##PASSIVE POP.

##Pf vs. Pm 

#Density plot using log scale
dplot_PfbyPm<-ggplot(PfvsPfPm_passive)+
  geom_density(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.5)+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")
plot(dplot_PfbyPm)

#histogram plot using log scale
hplot_PfbyPm<-ggplot(PfvsPfPm_passive)+
  geom_histogram(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.3)+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")
plot(hplot_PfbyPm)

#overlaying density plot onto histogram to show the raw counts too
#un-stratified counts in background: 
PfbyPm_overlay_Pas<-ggplot(PfvsPfPm_passive, aes(x=Pd_Pf, bins=10))+
  geom_histogram(aes(y=..density..), alpha=0.5, col='gray', lwd=0.4)+
  #scale_x_log10(name="Estimated Pf Parasite Density (p/uL)")+
  scale_x_log10(name="Estimated Pf parasite density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.4)
#+labs(title="Pf mono vs. Pf+Pm mixed infections - Passive visits", fontsize = 8) 

PfbyPm_overlay_Pas_fig <- PfbyPm_overlay_Pas + guides(fill=guide_legend(title=""))
plot(PfbyPm_overlay_Pas_fig)

#Stratified (but stacked) counts in background:
PfbyPm_overlay_Pas<-ggplot(PfvsPfPm_passive, aes(x=Pd_Pf, fill=Pf_Species, bins=10))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.3, col='gray', lwd=0.2)+
  scale_x_log10(name="Estimated Pf parasite density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.2)+
  scale_fill_manual(values=c("#5b9bf5", "#f5cc5b"))
#+labs(title="Pf mono vs. Pf+Pm mixed infections - Clinic Sub-Pop.", fontsize = 8)


PfbyPm_overlay_Pas_fig <- PfbyPm_overlay_Pas + guides(fill=guide_legend(title="")) +
  # annotate("text", label = "LoD", x = 7, y = 0.87, size = 5, gend(title="")) +
  # geom_vline(xintercept = 10, linetype="solid", color = "dark grcolour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text=element_text(size=14), legend.position = c(0.79, 0.65))


PfbyPm_overlay_Pas_fig2 <- PfbyPm_overlay_Pas_fig + expand_limits(x=c(0.1 ,3000000), y=c(0.0, 0.80))
plot(PfbyPm_overlay_Pas_fig2)


## Pf vs. Po

#Density plot using log scale
dplot_PfbyPo<-ggplot(PfvsPfPo_passive)+
  geom_density(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.5)+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")
plot(dplot_PfbyPo)

#histogram plot using log scale
hplot_PfbyPo<-ggplot(PfvsPfPo_passive)+
  geom_histogram(aes(x = Pd_Pf, fill = Pf_Species), alpha=0.3)+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")
plot(hplot_PfbyPo)

#overlaying density plot onto histogram to show the raw counts too
#un-stratified counts in background: 
PfbyPo_overlay_Pas<-ggplot(PfvsPfPo_passive, aes(x=Pd_Pf, bins=10))+
  geom_histogram(aes(y=..density..), alpha=0.5, col='gray', lwd=0.4)+
  scale_x_log10(name="Estimated Pf parasite density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.4)
#+labs(title="Pf mono vs. Pf+Po mixed infections - Passive visits", fontsize = 8)

PfbyPo_overlay_Pas_fig <- PfbyPo_overlay_Pas + guides(fill=guide_legend(title=""))
plot(PfbyPo_overlay_Pas_fig)

#Stratified (but stacked) counts in background:
PfbyPo_overlay_Pas<-ggplot(PfvsPfPo_passive, aes(x=Pd_Pf, fill=Pf_Species, bins=10))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.3, col='gray', lwd=0.2)+
  scale_x_log10(name="Estimated Pf parasite density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pf_Species)), col='black', lwd=0.6, alpha=0.2)+
  scale_fill_manual(values=c("#c44343", "#f5cc5b"))
#+labs(title="Pf mono vs. Pf+Po mixed infections - Clinic Sub-Pop.", fontsize = 8)


PfbyPo_overlay_Pas_fig <- PfbyPo_overlay_Pas + guides(fill=guide_legend(title="")) +
  # annotate("text", label = "LoD", x = 7, y = 0.87, size = 5, gend(title="")) +
  # geom_vline(xintercept = 10, linetype="solid", color = "dark grcolour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text=element_text(size=14), legend.position = c(0.79, 0.65))


PfbyPo_overlay_Pas_fig2 <- PfbyPo_overlay_Pas_fig + expand_limits(x=c(0.1 ,3000000), y=c(0.0, 0.80))
plot(PfbyPo_overlay_Pas_fig2)



###Aggregating the Pf parasitemia figures
gridExtra::grid.arrange(PfbyPm_overlay_fig2, PfbyPo_overlay_fig2, 
                        PfbyPm_overlay_Pas_fig2, PfbyPo_overlay_Pas_fig2,   ncol=2, nrow=2)   #print Pf species density graphs in one frame

## Exporting aggregate image as a high-res file
png("Pf_mixedvsmono_ParasiteDensity_PublicationFig.png", width =14, height = 14, units = "cm", res = 600)
dev.off()



#######  Now Plotting Pm and Po density + Histogram plots comparing mixed vs. mono in Active and Passive visits 


##### ACTIVE VISITS

## Pm 
Pm_overlay_Act<-ggplot(pm_pd_active, aes(x=Pd_Pm_mult4, fill=Pm_Species, bins=10))+
  geom_histogram(aes(y=..density..), alpha=0.3, col='gray', lwd=0.2)+
  scale_x_log10(name="Estimated Pm Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pm_Species)), color=Pm_Species, lwd=1.0, alpha=0.4)+
  scale_fill_manual(values=c("#acbdf2", "#294ab3"))
#+labs(title="Pm mono vs. mixed infections - Active visits", fontsize = 6)

Pm_overlay_Act_fig <- Pm_overlay_Act + guides(fill=guide_legend(title="")) +
  geom_vline(xintercept = 10, linetype="solid", color = "dark gray", size=0.6) +
  annotate("text", label = "Pm LoD", x = 3.0, y = 1.8, size = 5, colour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text=element_text(size=20)) 
#legend.position = c(0.87, 0.25)))

Pm_overlay_Act_fig2 <- Pm_overlay_Act_fig + expand_limits(x=c(0.001,3000000))

plot(Pm_overlay_Act_fig2)



## Po 
Po_overlay_Act<-ggplot(po_pd_active, aes(x=Pd_Po_mult4, fill=Po_Species, bins=10))+
  # geom_histogram(aes(y=..density..), alpha=0.3, col='gray', lwd=0.2)+
  scale_x_log10(name="Estimated Po Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Po_Species)), col="black", lwd=1.2, alpha=0.6)+
  scale_fill_manual(values=c("#edbbbb", "#780101"))
#+labs(title="Po mono vs. mixed infections - Active visits", fontsize = 6)

Po_overlay_Act_fig <- Po_overlay_Act + guides(fill=guide_legend(title="")) +
  geom_vline(xintercept = 10, linetype="solid", color = "dark gray", size=0.6) +
  annotate("text", label = "Po LoD", x = 3.5, y = 1.8, size = 5, colour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text=element_text(size=20)) 
#legend.position = c(0.87, 0.25)))
Po_overlay_Act_fig2 <- Po_overlay_Act_fig + expand_limits(x=c(0.001,3000000))

plot(Po_overlay_Act_fig2)



##### PASSIVE VISITS

## Pm 
Pm_overlay_Pas<-ggplot(pm_pd_passive, aes(x=Pd_Pm_mult4, fill=Pm_Species, bins=10))+
  geom_histogram(aes(y=..density..), alpha=0.3, col='gray', lwd=0.2)+
  scale_x_log10(name="Estimated Pm Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Pm_Species)), col='black', lwd=1, alpha=0.4)
#+labs(title="Pm mono vs. mixed infections - Passive visits", fontsize = 8)

Pm_overlay_Pas_fig <- Pm_overlay_Pas + guides(fill=guide_legend(title="")) +
  geom_vline(xintercept = 10, linetype="solid", color = "dark gray", size=1.0) +
  annotate("text", label = "Pm LoD", x = 3.0, y = 1.8, size = 3, colour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  expand_limits(x=c(0,3000000))
#legend.position = c(0.87, 0.25)))

Pm_overlay_Pas_fig2 <- Pm_overlay_Pas_fig + expand_limits(x=c(0,3000000))

plot(Pm_overlay_Pas_fig2)


## Po 
Po_overlay_Pas<-ggplot(po_pd_passive, aes(x=Pd_Po_mult4, fill=Po_Species, bins=10))+
  geom_histogram(aes(y=..density..), alpha=0.3, col='gray', lwd=0.2)+
  scale_x_log10(name="Estimated Po Parasite Density (p/uL)", breaks=c(.01,.1,1,10, 100, 1000, 10000, 100000, 1000000),labels=xlabels_manual)+
  geom_density(aes(fill=factor(Po_Species)), col='black', lwd=1, alpha=0.4)
#+labs(title="Po mono vs. mixed infections - Passive visits", fontsize = 8)

Po_overlay_Pas_fig <- Po_overlay_Pas + guides(fill=guide_legend(title="")) +
  geom_vline(xintercept = 10, linetype="solid", color = "dark gray", size=1.0) +
  annotate("text", label = "Po LoD", x = 4.0, y = 1.8, size = 3, colour = "dark gray")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
#legend.position = c(0.87, 0.25)))
Po_overlay_Pas_fig2 <- Po_overlay_Pas_fig + expand_limits(x=c(0,3000000))

plot(Po_overlay_Pas_fig2)




###Aggregating ALL parasitemia figures
gridExtra::grid.arrange(PfbyPm_overlay_fig2, PfbyPo_overlay_fig2, 
                        PfbyPm_overlay_Pas_fig2, PfbyPo_overlay_Pas_fig2,   ncol=2, nrow=2)   #print Pf species density graphs in one frame


gridExtra::grid.arrange(Pm_overlay_Act_fig2, Po_overlay_Act_fig2, 
                        Pm_overlay_Pas_fig2, Po_overlay_Pas_fig2,   ncol=2, nrow=2)   #print Pm and Po species density graphs in one frame



##### Creating additional plots for visualization and descriptive review of data 
# parasitemia for total pop. 
#P.malariae - mixed vs. mono
pm<-ggplot(pm_pd)+
  geom_histogram(aes(Pd_Pm_mult4, fill=Pm_Species), color="gray80",bins=10)+
  scale_fill_brewer(palette = "Blues", direction=1, name="Species")+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")+
  guides(fill = guide_legend(reverse=T))+
  ggtitle(expression(italic("P. malariae (Total Pop.)")))+
  geom_vline(aes(xintercept = median(Pd_Pm_mult4[Pm_Species=="Pm Mono"]))) + 
  geom_vline(aes(xintercept = median(Pd_Pm_mult4[Pm_Species=="Pm Mixed"])))+
  theme_bw()+
  theme(legend.position = c(0.93, 0.8),
        plot.title = element_text(hjust = 0.5))

print(pm)

#P.ovale - mixed vs. mono
po<-ggplot(po_pd)+
  geom_histogram(aes(Pd_Po_mult4, fill=Po_Species), color="gray80",bins=10)+
  scale_fill_brewer(palette = "Reds", direction=1, name="Species")+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")+
  guides(fill = guide_legend(reverse=T))+
  geom_vline(aes(xintercept = median(Pd_Po_mult4[Po_Species=="Po Mono"]))) + 
  geom_vline(aes(xintercept = median(Pd_Po_mult4[Po_Species=="Po Mixed"])))+
  ggtitle(expression(italic("P. ovale spp. (Total Pop.)")))+
  theme_bw()+
  theme(legend.position = c(0.93, 0.8),
        plot.title = element_text(hjust = 0.5))

print(po)


#P.falciparum - mixed vs. mono
pf<-ggplot(pf_pd)+
  geom_histogram(aes(Pd_Pf, fill=Pf_Species), color="gray80",bins=10)+
  scale_fill_brewer(palette = "Greens", direction=1, name="Species")+
  scale_x_log10(name="Estimated Parasite Density (p/uL)")+
  guides(fill = guide_legend(reverse=T))+
  ggtitle(expression(italic("P. falciparum (Total Pop.)")))+
  geom_vline(aes(xintercept = median(Pd_Pf[Pf_Species=="Pf Mono"]))) + 
  geom_vline(aes(xintercept = median(Pd_Pf[Pf_Species=="Pf Mixed"])))+
  theme_bw()+
  theme(legend.position = c(0.93, 0.8),
        plot.title = element_text(hjust = 0.5))

print(pf)


gridExtra::grid.arrange(pm, po, pf, ncol=3)   #print all 3 species density graphs in one frame


ggsave(filename=".../DRC_nonfalcip_active_pd.png", device="png",
       height=4, width=8, units="in", dpi=600, bg="transparent")

