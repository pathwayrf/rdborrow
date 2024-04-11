############################
# Creation of Part 2 Data Set
#############################

### Variables in the data set

### Variables we want:

# Treatment
# Time (Weeks)
# Age at enrollment (years)
# SMA Type (2 or 3)
# Ambulatory Status
# SMN2 Gene Copy Number
# Presence of scoliosis
# MFM at Baseline
# MFM at 12 months (primary endpoints+ all the followup visits)


library(dplyr)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(labelled)


load("Study/output/rct.bene.RData")
load("Study/output/rct.outcome.RData")
load("Study/output/ext.bene.RData")
load("Study/output/ext.outcome.RData")



# For the trial, we should use for those
# Age >=6: variable15, MFM 32 Total (recalculated)
# Age < 6: variable8, MFM 20 Total (recalculated)
# 
# For the EC, we should use
# variable18: MFM Total Score derived
# 
# DO NOT USE
# ‘Actual’, ‘forced MFM20’ or ‘forced MFM32’


######################### Create External control analytic dataset

## patient char file

# cols in adsub.sas7bdat
#cols=pool_metadata[[2]]%>%filter(Filename %in% "adsub.sas7bdat")%>%pull(`Name of\nVariable`)

# select patient char to keep for analyses
ext.bene<-(sap_2 %>% filter(PARAM == "Ambulatory Status")) %>%
  select(UNI_ID, STUDYID, ACTARM, Ambulatory_Status = AVALC) %>%
  left_join((sap_2 %>% filter(PARAM == "SMA Type"))) %>%
  select(UNI_ID, STUDYID,ACTARM,Ambulatory_Status, SMA_Type = AVALC) %>%
  left_join((sap_2 %>% filter(PARAM == "SMN2 Copy Number"))) %>%
  select(UNI_ID, STUDYID,ACTARM,Ambulatory_Status, SMA_Type, SMN2_Copy_Number = AVALC) %>%
  left_join((sap_2 %>% filter(PARAM == "Scoliosis at Baseline (Y/N)"))) %>%
  select(UNI_ID, STUDYID, ACTARM,Ambulatory_Status, SMA_Type, SMN2_Copy_Number, Scoliosis = AVALC, Age=AAGE)%>%
  filter(!STUDYID %in% "BP39055")%>%
  mutate(Study=ifelse(STUDYID=="BP29540","Natural History","Previous Placebo"))
 
#View(ext.bene)
 
## longitudinal outcome file

cols=pool_metadata[[2]]%>%filter(Filename %in% "adza.sas7bdat")%>%pull(`Name of\nVariable`)

print(pool_metadata[[2]]%>%filter(Filename %in% "adza.sas7bdat"),n=100)

table(sap_3$PARAM) # correct MFM used !
ext.outcome<-sap_3%>%
  filter(!STUDYID %in% "BP39055")%>%
  filter(PARAM == "MFM Total Score derived")%>%
  select(UNI_ID, STUDYID, MFM32=AVAL, AVISIT, AVISITN, change_from_baseline=CHG)%>%
  arrange(UNI_ID,AVISIT)%>%
  mutate(change_from_baseline=ifelse(is.na(change_from_baseline),0,change_from_baseline))%>%
  filter(AVISIT %in% c("BASELINE","Week 26",  "Week 52",  "Week 104", "Week 78"))%>%
  mutate(Study=ifelse(STUDYID=="BP29540","Natural History","Previous Placebo"))%>%
  mutate(time=tidyr::extract_numeric(AVISIT))

#View(ext.outcomes)



######################### Create RCT analytic dataset
 
# select patient char to keep for analyses
rct.bene<-(part2_data_adsub %>%
             filter(PARAM == "Ambulatory Status") %>%
             select(UNI_ID,STUDYID,ACTARM,
                    Ambulatory_Status = AVALC)) %>%
  left_join((part2_data_adsub %>%
               filter(PARAM == "SMA Type") %>%
               select(UNI_ID, SMA_Type = AVALC))) %>%
  left_join((part2_data_adsub %>%
               filter(PARAM == "SMN2 Copy Number (per medical records)") %>%
               select(UNI_ID, SMN2_Copy_Number = AVALC))) %>%
  left_join((part2_data_adsub %>%
               filter(PARAM == "Scoliosis") %>%
               select(UNI_ID, Scoliosis = AVALC, Age=AAGE)))%>%
  mutate(Study = "Sunfish Part 2")
 
table(part2_data_adza$PARAM[])
part2_data_adza$AAGE
head(part2_data_adza)


rct.outcome=part2_data_adza %>%
               filter(case_when(AGE<6 ~ PARAM %in% "MFM 20 Total (recalculated)", T ~ PARAM %in% "MFM 32 Total (recalculated)")) %>% 
               filter(AVISIT %in% c("BASELINE","Week 17","Week 35","Week 52","Week 78","Week 104")) %>%  
               select(UNI_ID, STUDYID, MFM32 = AVAL,AVISIT,AVISITN, change_from_baseline=CHG)%>% #, AGE, PARAM
  arrange(UNI_ID,AVISIT)%>%
  mutate(change_from_baseline=ifelse(is.na(change_from_baseline),0,change_from_baseline))%>%
  mutate(Study = "Sunfish Part 2")%>%
  mutate(time=tidyr::extract_numeric(AVISIT))
 
 
View(rct.bene)
View(rct.outcome)

bene=rbind(rct.bene, ext.bene)
outcome=rbind(rct.outcome, ext.outcome)

View(bene)
View(outcome)

save(bene, file = "Study/output/bene.RData")
save(outcome, file = "Study/output/outcome.RData")
 

