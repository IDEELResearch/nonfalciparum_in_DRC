##########################################
#Programmer: `Rachel Sendor
#Last update: Mar 2023
#Purpose:     DRC Non-falciparum Descriptive Epidemiology Analysis
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
library(openxlsx)
library(tableone)
library(purrr)



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

# STUDY ANALYSES: Assessing Repeat Pm, Po, and Pf Infections Over Study Period

################################################################################################


#________________________________________________________________________________________________________________________________________________________________________________

# SUMMARIZING NUMBER OF REPEATED PM, PO, AND PF INFECTIONS OVER TIME 

# Notes: 
#   - Producing in-text manuscript numbers for # of repeated non-pf infections per subject 
#   - Analyzing timing of repeat infections across study period (histograms by subject)

#________________________________________________________________________________________________________________________________________________________________________________



#Count total number of each species' infections by subject (non-missings)
AD_Total_Counts <- AD_Total %>% 
  group_by(Subject_ID) %>%
  dplyr::summarise(Pm_sum = sum((Pm), na.rm = TRUE), Po_sum = sum((Po), na.rm = TRUE), Pf_sum = sum((Pf), na.rm = TRUE))

AD_Active_Counts <- AD_Active %>% 
  group_by(Subject_ID) %>%
  dplyr::summarise(Pm_sum = sum((Pm), na.rm = TRUE), Po_sum = sum((Po), na.rm = TRUE), Pf_sum = sum((Pf), na.rm = TRUE))

AD_Passive_Counts <- AD_Passive %>% 
  group_by(Subject_ID) %>%
  dplyr::summarise(Pm_sum = sum((Pm), na.rm = TRUE), Po_sum = sum((Po), na.rm = TRUE), Pf_sum = sum((Pf), na.rm = TRUE))


AD_Total_Counts<-AD_Total_Counts%>%mutate(Pm_multiple=ifelse(Pm_sum>1, 1,0))%>%
  mutate(Pm_any=ifelse(Pm_sum>0, 1,0))%>%
  mutate(Po_multiple=ifelse(Po_sum>1, 1,0))%>%
  mutate(Po_any=ifelse(Po_sum>0, 1,0))%>%
  mutate(Pf_multiple=ifelse(Pf_sum>1, 1,0))%>%
  mutate(Pf_any=ifelse(Pf_sum>0, 1,0))%>%
  mutate(Pf_multiple_GT5=ifelse(Pf_sum>4, 1,0))%>%
  mutate(Pf_any=ifelse(Pf_sum>0, 1,0))

addmargins(table(AD_Total_Counts$Pm_sum, useNA = "always"))
addmargins(table(AD_Active_Counts$Pm_sum, useNA = "always"))
addmargins(table(AD_Passive_Counts$Pm_sum, useNA = "always"))


##Proportion of multiple infections out of total
addmargins(table(AD_Total_Counts$Pm_multiple[AD_Total_Counts$Pm_any==1]))
prop.table(table(AD_Total_Counts$Pm_multiple[AD_Total_Counts$Pm_any==1]))


addmargins(table(AD_Total_Counts$Po_sum, useNA = "always"))
addmargins(table(AD_Active_Counts$Po_sum, useNA = "always"))
addmargins(table(AD_Passive_Counts$Po_sum, useNA = "always"))


##Proportion of multiple infections out of total
addmargins(table(AD_Total_Counts$Po_multiple[AD_Total_Counts$Po_any==1]))
prop.table(table(AD_Total_Counts$Po_multiple[AD_Total_Counts$Po_any==1]))


addmargins(table(AD_Total_Counts$Pf_sum, useNA = "always"))
addmargins(table(AD_Active_Counts$Pf_sum, useNA = "always"))
addmargins(table(AD_Passive_Counts$Pf_sum, useNA = "always"))


##Proportion of multiple infections out of total
addmargins(table(AD_Total_Counts$Pf_multiple[AD_Total_Counts$Pf_any==1]))
prop.table(table(AD_Total_Counts$Pf_multiple[AD_Total_Counts$Pf_any==1]))

addmargins(table(AD_Total_Counts$Pf_multiple_GT5[AD_Total_Counts$Pf_any==1]))
prop.table(table(AD_Total_Counts$Pf_multiple_GT5[AD_Total_Counts$Pf_any==1]))



Pm_multiple<-ggplot(data=AD_Total_Counts)+
  geom_histogram(aes(Pm_sum))
plot(Pm_multiple)

Po_multiple<-ggplot(data=AD_Total_Counts)+
  geom_histogram(aes(Po_sum))
plot(Po_multiple)

Pf_multiple<-ggplot(data=AD_Total_Counts)+
  geom_histogram(aes(Pf_sum))
plot(Pf_multiple)




################################################################################################

#Plotting distribution of repeated infections over time 

################################################################################################


Infections_All <- AD_Total%>% 
   mutate(Visit_Type = ifelse(Visit == 0, "Baseline Visit",
                             ifelse(Visit>0&Visit<4, "Active FU Visit", 
                                    ifelse(Visit>3, "Passive FU Visit", NA))))%>%
  mutate(Village_Char = ifelse(Subject_ID>3070000&Subject_ID<3080000, "a_Lingwala",
                               ifelse(Subject_ID>1010000&Subject_ID<1020000, "g_Bu",
                                      ifelse(Subject_ID>1020000&Subject_ID<1030000, "f_Impuru",
                                             ifelse(Subject_ID>1030000&Subject_ID<1040000, "e_Pema",
                                                    ifelse(Subject_ID>2060000&Subject_ID<2070000, "b_Kimpoko",
                                                           ifelse(Subject_ID>2040000&Subject_ID<2050000, "d_Ngamanzo",
                                                                  ifelse(Subject_ID>2050000&Subject_ID<2060000, "c_Iye", NA))))))))

table(Infections_All$Village,Infections_All$Village, useNA="always")
##Confirmed that Village_New is correctly coded. 


#Sorting by SubjectID and Visit_date
Infections_All<- Infections_All%>%arrange(Subject_ID, Visit_date)

#Counting total number of subjects in the dataset 
counting_subs<-Infections_All%>%distinct(Subject_ID)


#Sum total number of infections per subject to have a by-subject dataset & create an artificial subject ID that is consecutively numbered as in Subject ID to see clustering by

total_counts <- Infections_All %>% 
  dplyr::group_by(Subject_ID) %>%
  dplyr::summarise(Pm_studytot = sum(Pm == 1, na.rm=T), Po_studytot=sum(Po == 1, na.rm = T), Pf_studytot=sum(Pf ==1, na.rm=T))%>%
  dplyr::mutate(Sub_ID_ConsecNum = row_number())


#Counting subjects by number of infections in the study (active and passive)
addmargins(table(total_counts$Pf_studytot))
addmargins(table(total_counts$Pm_studytot))
addmargins(table(total_counts$Po_studytot))



#Now link counts of Pm and Po totals to the full dataset to be able to descending sort by them 
#one to many merge
nonfalcip_persubvisit_totals <- merge(Infections_All, total_counts, by= 'Subject_ID')

#many to one merge too
Infections_All2<-Infections_All%>%distinct(Subject_ID, .keep_all = TRUE)
total_counts2 <-total_counts%>%left_join(Infections_All2, by='Subject_ID')


###Creating histograms of the number of infections per subject (by visit) to use in later plots for marginal plots
MarginalHist_InfectCount_Pm<-ggplot(total_counts2, aes(x=Pm_studytot)) +
  geom_bar(color="Black", fill="Blue")+
  geom_text(stat='count', aes(label=after_stat(count)), hjust= -0.3)+
  #xlim(-1, 4)+
  # ylim(0,400)+
  scale_x_continuous(name="# Pm Ifxns per subject", breaks=c(0, 1, 2, 3, 4), limits=c(-0.5,4.5))+
  # scale_x_discrete(name="# Pm Ifxns per subject", breaks=c("0","1","2","3","4"), limits=c("0","1","2","3","4"))+
  xlab("# Pm infections per subject")+
  ylab("# Subjects")+
  theme_bw()+
  coord_flip()+
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))

print(MarginalHist_InfectCount_Pm)


MarginalHist_InfectCount_Po<-ggplot(total_counts2, aes(x=Po_studytot)) +
  geom_bar(color="Black", fill="Red")+
  geom_text(stat='count', aes(label=after_stat(count)), hjust= -0.3)+
  scale_x_continuous(name="# Po Ifxns per subject", breaks=c(0, 1, 2, 3, 4), limits=c(-0.5,4.5))+
  # xlab("# Po infections per subject")+
  ylab("# Subjects")+
  theme_bw()+
  coord_flip()+
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))

print(MarginalHist_InfectCount_Po)



#NOW Sum total number of infections per Household to have a by-household dataset 
total_counts_HH <- Infections_All %>% 
  group_by(Household) %>%
  summarise(Pm_studytot_HH = sum(Pm == 1, na.rm=T), Po_studytot_HH=sum(Po == 1, na.rm = T), Pf_studytot_HH=sum(Pf ==1, na.rm=T))


#Now link counts of Pm and Po totals to the full dataset to be able to descending sort by them 
#one to many merge
nonfalcip_perHHvisit_totals <- merge(Infections_All, total_counts_HH, by= "Household")





##########################################################################

#CREATING HISTOGRAMS BY SUBJECT_ID, VISIT DATE, & SPECIES 

##########################################################################

##Basic scatterplot of calendar time on X-axis and new consecutive subject ID on Y axis.

####SCATTERPLOT OF ALL STUDY VISITS OVER TIME, BY VISIT TYPE
Plot_Visits_All <- ggplot(nonfalcip_persubvisit_totals, mapping= aes(x=Visit_date, y=Sub_ID_ConsecNum))+
  geom_point(aes(color=factor(Visit_Type))) +
  labs(x= 'Date of Visit', y='Subject ID')+ 
  facet_grid(Village_Char ~., scales="free", space="free")   

plot(Plot_Visits_All)


#Scatterplot of each visit by subject and species type 

#Pm

# First - sort descending by Pm_studytot so that those with multiple infections are at the top of the plot
#df_pm_any<-nonfalcip_persubvisit_totals %>% group_by(subject_id_num) %>% mutate(Pm_subtot=sum(Pm_studytot))
df_pm_any<-nonfalcip_persubvisit_totals
df_pm_any$Sub_ID_ConsecNum<-reorder(df_pm_any$Sub_ID_ConsecNum, df_pm_any$Pm_studytot)

#Sorted by total infections over study period - all subjects 
p_species_pm_any <- ggplot(df_pm_any, aes(x = Visit_date, y = Sub_ID_ConsecNum)) +
  geom_point(aes(color=factor(Pm)), alpha=0.6) +
  scale_color_manual(values = c("1" = 'blue', "0" = 'gray75'), labels = c("Yes", "No")) + 
  theme_clean(base_size =10) +
  theme(axis.line.x.bottom=element_line(color="black"),
        #axis.line.y.right=element_line(color="black"),
        axis.ticks.y =element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank() ) +
  labs(x= 'Date of Visit', y='Individual Study Participants') +
  #  geom_line(aes(group=Subject_ID))+
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))



p_species_pm_any <- p_species_pm_any + 
  guides(color=guide_legend(title = 'P.malariae PCR+'))+
  labs(title = "Pm Infections Throughout Study Follow-up  - Total Study Pop")
plot(p_species_pm_any)


##########TRYING TO FIRST PLOT ONLY THOSE WITH PM INFECTIONS AND COLOR BY HH ID TO IDENTIFY SAME HH.  THEN OVERLAY ONTO GRAPH OF ALL SUBJECTS TO SEE FULL DENOM OF THOSE WITHOUT AN INFECTION. 

##STEP 1
#subset to only those with at least 1 Pm infection during the study
#df_pm_only <-nonfalcip_persubvisit_totals%>%filter(Pm==1)

#sort descending by Pm_studytot so that those with multiple infections are at the top of the plot
##df_pm_only2<-df_pm_only %>% group_by(Subject_ID) %>% mutate(Pm_subtot=sum(Pm_studytot))
df_pm_only$Subject_ID<-reorder(df_pm_only$Subject_ID, df_pm_only$Pm_studytot)


#Sorted by total infections over study period - only  subjects with at least 1 Pm during study 
p_species_pm_only <- ggplot(df_pm_only, aes(x = Visit_date, y = Subject_ID)) +
  geom_point(aes(color=factor(Household)), alpha=0.5) +
  theme(text=element_text(size=5), axis.text.y=element_blank(), axis.ticks.y=element_blank() ) +
  labs(x= 'Date of Visit', y='Subject ID') +
  geom_line(aes(group=Subject_ID))+
  facet_grid(Village_Char ~., scales="free", space="free")


p_species_pm_only <- p_species_pm_only + guides(color=guide_legend(title = 'Households with Pm+'))
plot(p_species_pm_only)
p_species_pm_only <- p_species_pm_only + guides(color=guide_legend(title = 'P.malariae PCR+'))
plot(p_species_pm_only)



###TESTING SUBSET TO ONLY MULTIPLE PM INFECTIONS TO LOOK BY HOUSEHOLD

#subset to only those with at least 1 Pm infection during the study
df_pm_any2 <-nonfalcip_persubvisit_totals%>%filter(Pm_studytot>0)

#sort descending by Pm_studytot so that those with multiple infections are at the top of the plot
df_pm_any2$Subject_ID<-reorder(df_pm_any2$Subject_ID, df_pm_any2$Pm_studytot)


#Sorted by total infections over study period - only  subjects with at least 1 Pm during study 
# New facet label names for village var
Village_Char.labs <- c("Lingwala", "Bu", "Impuru", "Pema", "Kimpoko", "Ngamanzo", "Iye")
names(Village_Char.labs) <- c("a_Lingwala", "g_Bu", "f_Impuru", "e_Pema", "b_Kimpoko", "d_Ngamanzo", "c_Iye")


p_species_pm_all <- ggplot(df_pm_any2, aes(x = Visit_date, y = Subject_ID)) +
  geom_line(aes(group=Subject_ID), color="gray80")+
  geom_point(aes(color=factor(Pm_species)), size = 1.9, alpha=0.6) +
  # geom_point(shape=1, color="gray30", size=2)+
  theme_clean(base_size=10, ) +
  labs(x= 'Date of Visit', y='Subject ID') +
  scale_color_manual(values = c("Pm mono" = '#0080FF', "Pm mixed" = '#000080'), na.value="gray70") + 
  theme(axis.text.y=element_text(size=4), axis.ticks.y=element_blank() ) +
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))


p_species_pm_all <- p_species_pm_all + 
  guides(color=guide_legend(title = 'P.malariae PCR+'))+
  labs(title = "Subjects with >=1 Pm Infection During the Study - Total Study Pop")
plot(p_species_pm_all)






#subset to only those with at least 1 Pm infection during the study
df_pm_mult_only <-nonfalcip_persubvisit_totals%>%filter(Pm_studytot>1)

#sort descending by Pm_studytot so that those with multiple infections are at the top of the plot
##df_pm_only2<-df_pm_only %>% group_by(Subject_ID) %>% mutate(Pm_subtot=sum(Pm_studytot))
df_pm_mult_only$Subject_ID<-reorder(df_pm_mult_only$Subject_ID, df_pm_mult_only$Pm_studytot)


#Sorted by total infections over study period - only  subjects with at least 1 Pm during study 
# New facet label names for village var
Village_Char.labs <- c("Lingwala", "Bu", "Impuru", "Pema", "Kimpoko", "Ngamanzo", "Iye")
names(Village_Char.labs) <- c("a_Lingwala", "g_Bu", "f_Impuru", "e_Pema", "b_Kimpoko", "d_Ngamanzo", "c_Iye")


p_species_pm_only <- ggplot(df_pm_mult_only, aes(x = Visit_date, y = Subject_ID)) +
  geom_line(aes(group=Subject_ID), color="gray80")+
  geom_point(aes(color=factor(Pm)), size = 2, alpha=0.6) +
  # geom_point(shape=1, color="gray30", size=2)+
  theme_clean(base_size=10,)+
  theme(axis.ticks.y =element_blank(),
        axis.text.y=element_blank(),) +
  labs(x= 'Date of Visit', y='Individual Study Participants') +
  # scale_color_manual(values = c("Pm mono" = '#0080FF', "Pm mixed" = '#000080'), na.value="gray85") + 
  scale_color_manual(values = c("1" = 'blue', "0" = 'white'), labels = c("Yes", ""), na.value="white") + 
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))


p_species_pm_only <- p_species_pm_only + 
  guides(color=guide_legend(title = 'P.malariae PCR+'))+
  labs(title = "Subjects with Multiple Pm Infections - Total Study Pop")
plot(p_species_pm_only)





##############
###NOW CREATING AS ABOVE FOR PO INFECTIONS 

#Po

# First - sort descending by Po_studytot so that those with multiple infections are at the top of the plot
df_po_any<-nonfalcip_persubvisit_totals
df_po_any$Sub_ID_ConsecNum<-reorder(df_po_any$Sub_ID_ConsecNum, df_po_any$Po_studytot)

#Sorted by total infections over study period - all subjects 
p_species_po_any <- ggplot(df_po_any, aes(x = Visit_date, y = Sub_ID_ConsecNum)) +
  geom_point(aes(color=factor(Po)), alpha=0.70) +
  scale_color_manual(values = c("1" = '#FF3131', "0" = 'gray75'), labels = c("Yes", "No")) + 
  theme_clean(base_size =10) +
  theme(axis.line.x.bottom=element_line(color="black"),
        #axis.line.y.right=element_line(color="black"),
        axis.ticks.y =element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank() ) +
  labs(x= 'Date of Visit', y='Individual Study Participants') +
  #  geom_line(aes(group=Subject_ID))+
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))



p_species_po_any <- p_species_po_any + 
  guides(color=guide_legend(title = 'P.ovale PCR+'))+
  labs(title = "Po Infections Throughout Study Follow-up  - Total Study Pop")
plot(p_species_po_any)




##########TRYING TO FIRST PLOT ONLY THOSE WITH Po INFECTIONS AND COLOR BY HH ID TO IDENTIFY SAME HH.  THEN OVERLAY ONTO GRAPH OF ALL SUBJECTS TO SEE FULL DENOM OF THOSE WITHOUT AN INFECTION. 

##STEP 1
#subset to only those with at least 1 Po infection during the study
#df_po_only <-nonfalcip_persubvisit_totals%>%filter(Po==1)

#sort descending by Po_studytot so that those with multiple infections are at the top of the plot
##df_po_only2<-df_po_only %>% group_by(Subject_ID) %>% mutate(Po_subtot=sum(Po_studytot))
df_po_only$Subject_ID<-reorder(df_po_only$Subject_ID, df_po_only$Po_studytot)


#Sorted by total infections over study period - only  subjects with at least 1 Pm during study 
p_species_po_only <- ggplot(df_po_only, aes(x = Visit_date, y = Subject_ID)) +
  geom_point(aes(color=factor(Household)), alpha=0.5) +
  theme(text=element_text(size=5), axis.text.y=element_blank(), axis.ticks.y=element_blank() ) +
  labs(x= 'Date of Visit', y='Subject ID') +
  geom_line(aes(group=Subject_ID))+
  facet_grid(Village_Char ~., scales="free", space="free")


p_species_po_only <- p_species_po_only + guides(color=guide_legend(title = 'Households with Po+'))
plot(p_species_po_only)
p_species_po_only <- p_species_po_only + guides(color=guide_legend(title = 'P.ovale PCR+'))
plot(p_species_po_only)




#####LOOKING AT FOCUS OF THOSE WITH ONLY Po INFECTIONS 
#subset to only those with at least 1 Po infection during the study
df_po_any2 <-nonfalcip_persubvisit_totals%>%filter(Po_studytot>0)

#sort descending by Pm_studytot so that those with multiple infections are at the top of the plot
df_po_any2$Subject_ID<-reorder(df_po_any2$Subject_ID, df_po_any2$Po_studytot)


#Sorted by total infections over study period - only  subjects with at least 1 Pm during study 
# New facet label names for village var
Village_Char.labs <- c("Lingwala", "Bu", "Impuru", "Pema", "Kimpoko", "Ngamanzo", "Iye")
names(Village_Char.labs) <- c("a_Lingwala", "g_Bu", "f_Impuru", "e_Pema", "b_Kimpoko", "d_Ngamanzo", "c_Iye")


p_species_po_all <- ggplot(df_po_any2, aes(x = Visit_date, y = Subject_ID)) +
  geom_line(aes(group=Subject_ID), color="gray80")+
  geom_point(aes(color=factor(Po)), size = 2, alpha=0.9) +
  # geom_point(shape=1, color="gray30", size=2)+
  theme_clean(base_size=10) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank() ) +
  labs(x= 'Date of Visit', y='Individual Study Participants') +
  # scale_color_manual(values = c("Po mixed" = '#b30000'), na.value="gray85") + 
  scale_color_manual(values = c("1" = "#b30000", "0" = 'white'), labels = c("Yes", ""), na.value="white") + 
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))


p_species_po_all <- p_species_po_all + 
  guides(color=guide_legend(title = 'P.ovale PCR+'))+
  labs(title = "Subjects with >=1 Po Infection During the Study - Total Study Pop")
plot(p_species_po_all)






###TESTING SUBSET TO ONLY MULTIPLE Po INFECTIONS TO LOOK BY HOUSEHOLD
#subset to only those with at least 1 Po infection during the study
df_po_mult_only <-nonfalcip_persubvisit_totals%>%filter(Po_studytot>1)

#sort descending by Po_studytot so that those with multiple infections are at the top of the plot
df_po_mult_only$Subject_ID<-reorder(df_po_mult_only$Subject_ID, df_po_mult_only$Po_studytot)


#Sorted by total infections over study period - only  subjects with at least 1 Po during study 
# New facet label names for village var
Village_Char.labs <- c("Lingwala", "Bu", "Impuru", "Pema", "Kimpoko", "Ngamanzo", "Iye")
names(Village_Char.labs) <- c("a_Lingwala", "g_Bu", "f_Impuru", "e_Pema", "b_Kimpoko", "d_Ngamanzo", "c_Iye")


p_species_po_only <- ggplot(df_po_mult_only, aes(x = Visit_date, y = Subject_ID)) +
  geom_line(aes(group=Subject_ID), color="gray80")+
  geom_point(aes(color=factor(Po_species)), size = 2, alpha=0.9) +
  # geom_point(shape=1, color="gray30", size=2)+
  theme_clean(base_size=10, ) +
  labs(x= 'Date of Visit', y='Subject ID') +
  scale_color_manual(values = c("Po mono" = '#ff8080', "Po mixed" = '#b30000'), na.value="gray85") + 
  facet_grid(Village_Char ~., scales="free", space="free", labeller = labeller(Village_Char= Village_Char.labs))


p_species_po_only <- p_species_po_only + 
  guides(color=guide_legend(title = 'P.ovale PCR+'))+
  labs(title = "Subjects with Multiple Po Infections - Total Study Pop")
plot(p_species_po_only)






###################################3

#Create 2 plots next to each other for each non-pf species  

###Pm Figures 
Pm_Multiple_figure <- ggarrange(p_species_pm_any, p_species_pm_only,
                                labels = c("A", "B"),
                                ncol = 2)
plot(Pm_Multiple_figure)

Pm_Multiple_figure_2 <- ggarrange(p_species_pm_any, p_species_pm_all,
                                  labels = c("A", "B"),
                                  ncol = 2)
plot(Pm_Multiple_figure_2)



###Po Figures 
Po_Multiple_figure <- ggarrange(p_species_po_any, p_species_po_only,
                                labels = c("C", "D"),
                                ncol = 2)
plot(Po_Multiple_figure)


Po_Multiple_figure_2 <- ggarrange(p_species_po_any, p_species_po_all,
                                  labels = c("C", "D"),
                                  ncol = 2)
plot(Po_Multiple_figure_2)






###EXPLORING DATASETS FOR MULTIPLE INFECTIONS OVER TIME 

#Counting total number of Pf, Po, and Pm infections by visit type (infection totals, not subject totals)
addmargins(table(AD_Total$Pf, AD_Total$visit_type))
addmargins(table(AD_Total$Po, AD_Total$visit_type))
addmargins(table(AD_Total$Pm, AD_Total$visit_type))


#Counting number of repeat infections over time

addmargins(table(AD_Total$Pf, AD_Total$Visit))
addmargins(table(AD_Total$Pm, AD_Total$Visit))
addmargins(table(AD_Total$Po, AD_Total$Visit))




#Count total number of each species' infections by subject 
total_counts <- AD_Total %>% 
  group_by(Subject_ID) %>%
  summarise(Pm_studytot = sum(Pm, na.rm = T), Po_studytot = sum(Po, na.rm = T), Pf_studytot = sum(Pf, na.rm = T))


#Counting subjects by number of infections in the study (active and passive)
addmargins(table(total_counts$Pf_studytot))
addmargins(table(total_counts$Pm_studytot))
addmargins(table(total_counts$Po_studytot))

#Now link counts of Pm and Po totals to the full dataset to be able to descending sort by them 
#one to many merge
nonfalcip_persubvisit_totals <- merge(AD_Total, total_counts, by= 'Subject_ID')

#Reverse direction for a by_subject dataset 
HH_subset<-AD_Total%>%subset(select=c("Subject_ID", "Village", "Household"))
HH_subset2<-merge(x=total_counts, y=HH_subset, by="Subject_ID", all.x=TRUE)
HH_subset3=HH_subset2[which(!duplicated(HH_subset2$Subject_ID)),]


#Investigating if those with >1 Pm or >1 Po infection over time were from same village or same household
#Count total number of each species' infections by subject, across each village and HH

#Counting subjects by number of infections in the study (active and passive)
#Greater than 1 infection = multiple infections during study period 
Pf_multipleInfxn_Table<-addmargins(table(HH_subset3$Pf_studytot>1, HH_subset3$Village, useNA="always"))
print(Pf_multipleInfxn_Table)
prop.table(Pf_multipleInfxn_Table)

Po_multipleInfxn_Table<-addmargins(table(HH_subset3$Po_studytot>1, HH_subset3$Village, useNA="always"))
print(Po_multipleInfxn_Table)
prop.table(Po_multipleInfxn_Table)

Pm_multipleInfxn_Table<-addmargins(table(HH_subset3$Pm_studytot>1, HH_subset3$Village, useNA="always"))
print(Pm_multipleInfxn_Table)
prop.table(Pm_multipleInfxn_Table)

#Greater than 2 Infections 
Pf_multipleInfxn_Table<-addmargins(table(HH_subset3$Pf_studytot>2, HH_subset3$Village, useNA="always"))
print(Pf_multipleInfxn_Table)
prop.table(Pf_multipleInfxn_Table)

Po_multipleInfxn_Table<-addmargins(table(HH_subset3$Po_studytot>2, HH_subset3$Village, useNA="always"))
print(Po_multipleInfxn_Table)
prop.table(Po_multipleInfxn_Table)

Pm_multipleInfxn_Table<-addmargins(table(HH_subset3$Pm_studytot>2, HH_subset3$Village, useNA="always"))
print(Pm_multipleInfxn_Table)
prop.table(Pm_multipleInfxn_Table)

Pm_multipleInfxn_Table<-addmargins(table(HH_subset3$Pm_studytot>3, HH_subset3$Village, useNA="always"))
print(Pm_multipleInfxn_Table)
prop.table(Pm_multipleInfxn_Table)


##checking if same HH for subjects with 4 Pm infections 
Pm_multipleInfxn_Table2<-addmargins(table(HH_subset3$Pm_studytot>3, HH_subset3$Household, useNA="always"))
print(Pm_multipleInfxn_Table2)





#####################################################################################

###Scatterplots of Subject ID and Visit Date plotting infection occurrence over time 

#####################################################################################

#Pm

# First - sort descending by Pm_studytot so that those with multiple infections are at the top of the plot
nonfalcip_persubvisit_totals$Subject_ID<-reorder(nonfalcip_persubvisit_totals$Subject_ID, -nonfalcip_persubvisit_totals$Pm_studytot)

#Sorted by total infections over study period - all subjects 
p_species_pm_any <- ggplot(nonfalcip_persubvisit_totals, aes(x = Visit_date, y = Subject_ID)) +
  geom_point(aes(color=factor(Pm)), alpha=0.5) +
  scale_color_manual(values = c("0" = 'gray', "1" = 'red')) + 
  labs(x= 'Date of Visit', y='Subject ID')

p_species_pm_any <- p_species_pm_any + guides(color=guide_legend(title = 'P.malariae PCR+'))
plot(p_species_pm_any)

nonfalcip_persubvisit_totals <-nonfalcip_persubvisit_totals%>%
  mutate(HA=ifelse(Subject_ID>0&Subject_ID<2000000, "Bu", 
                   ifelse(Subject_ID>1999999&Subject_ID<3000000, "Kimpoko", 
                          ifelse(Subject_ID>2999999&Subject_ID<4000000, "Lingwala", NA))))%>%
  dplyr::select(Visit_date, Pm, Subject_ID, Pm_studytot, HA)





