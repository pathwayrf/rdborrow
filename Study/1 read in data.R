# By running this script, analytic datasets (patient level and outcome level datasets) will be created for external controls and RCT



library(dplyr)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(labelled)
 
 
# Reading in meta data from Part 1
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part1/part1_pooled.rds", file = "part1_pooled.rds")

# List of size 2 with meta data
pool_metadata <- readRDS("part1_pooled.rds")
 

# First element is a description of the data sets
# Second element is a key of the variable names
pool_metadata[[1]]
View(pool_metadata[[2]])
print(pool_metadata[[2]]%>%filter(Filename %in% "adsub.sas7bdat"),n = 100)

#pooled data with NHS and RWE Trial
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part1/part1_sas_data_pooled.rds", file = "part1_sas_data_pooled.rds")
list.files("Study/code")
sas_data_pooled <- readRDS("part1_sas_data_pooled.rds")

 
# Extracting data frames
sap_1 <- sas_data_pooled[[1]] # Subject level Analysis Data Set
sap_2 <- sas_data_pooled[[2]] # Subcategory Analysis, contains the actual data of interest
sap_3 <- sas_data_pooled[[3]] # Specific clinical data
 

# Data from all 3 studies
table(sap_1$STUDYID) 
# BP29540: natural history study 81
# BP39055: RCT
# WN29836: 57
table(sap_2$STUDYID)
table(sap_3$STUDYID)

#View(sap_1)
#View(sap_2)
#View(sap_3)


#Natural History Study
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part1/part1_sas_data_nhs.rds", file = "part1_sas_data_nhs.rds")

sas_data_nhs <- readRDS("part1_sas_data_nhs.rds")
sas_data_nhs$aphy
 
#RWE Trial
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part1/part1_sas_data_rwetrial.rds", file = "part1_sas_data_rwetrial.rds")

sas_data_rwetrial <- readRDS("part1_sas_data_rwetrial.rds")

 

######################### Create External control analytic dataset

## patient char file

cols=pool_metadata[[2]]%>%filter(Filename %in% "adsub.sas7bdat")%>%pull(`Name of\nVariable`)
ext.bene<-sap_2%>%
  filter(!STUDYID %in% "BP39055")%>%
  select(cols[c(1,2,3,4,5,6,7,9,11,13,14,15,21,22,23,24,25,26)] )%>%
  filter(! PARAMCD %in% "MTCPFL")

table(ext.bene$PARAM,ext.bene$PARAMCD)
#View(ext.bene)

# from long to wide
dat1=ext.bene%>%filter(PARAMCD %in% "BM5D")%>%select(UNI_ID, PARAMCD, AVAL )
dat2=spread(dat1, key = PARAMCD, value =AVAL)
ext.bene<-ext.bene%>%left_join(dat2,by="UNI_ID")
 
dat1=ext.bene%>%filter(!PARAMCD %in% "BM5D")%>%select(UNI_ID, PARAMCD, AVALC )
dat2=spread(dat1, key = PARAMCD, value =AVALC)
ext.bene<-ext.bene%>%left_join(dat2,by="UNI_ID")

ext.bene<-ext.bene%>%distinct(UNI_ID, .keep_all = TRUE)%>%select(-PARAMCD,-PARAM,-PARCAT1,-AVAL,-AVALC,-AVALU)

sap_2%>%select(PARAM, PARAMCD)%>%filter(! PARAMCD %in% "MTCPFL")%>%distinct(PARAM, PARAMCD)

ext.bene<-ext.bene%>%set_variable_labels(AMSTAT="Ambulatory Status",
                                      BM5D="MFM Derived Total Score at Baseline",
                                      M5DS="MFM Scale for M5D",
                                      SCOLIO="Scoliosis at Baseline (Y/N)",
                                      SMAAMST="SMA Type (Ambulatory Status)",
                                      SMATYPE="SMA Type",
                                      TGENCONO="SMN2 Copy Number",
                                      SCOCURVE="Scoliosis Grade")


## longitudinal outcome file

cols=pool_metadata[[2]]%>%filter(Filename %in% "adza.sas7bdat")%>%pull(`Name of\nVariable`)

print(pool_metadata[[2]]%>%filter(Filename %in% "adza.sas7bdat"),n=100)

 
temp=sap_3%>%
  filter(PARAMCD %in% c("M5D"))%>%
  select(PARAM)%>%
  distinct()

# MFM Total Score derived
temp=sap_3%>%
  select(PARAM)%>%
  distinct()
View(temp)
  
ext.outcome<-sap_3%>%
  filter(!STUDYID %in% "BP39055")%>%
  filter(PARAMCD %in% c("M5D"))%>%
  select(cols[c(1,2,11,41,42,43,44,47,48,52,53,54,55,63)] )%>%arrange(UNI_ID,AVISIT)%>%mutate(CHG=ifelse(is.na(CHG),0,CHG))


######################### Create RCT analytic dataset
 

#Part2adsl
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part2/part2adsl_data.rds", file = "part2adsl_data.rds")

part2_data_adsl <- readRDS("part2adsl_data.rds")
View(part2_data_adsl)

table(part2_data_adsl$ARM)
# only RCT data, subject-level data describing baseline char and randomization

#Part2adsub
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part2/part2adsub_data.rds", file = "part2adsub_data.rds")

part2_data_adsub <- readRDS("part2adsub_data.rds")
View(part2_data_adsub)
 

#Part2adza
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part2/part2adza_data.rds", file = "part2adza_data.rds")

part2_data_adza <- readRDS("part2adza_data.rds")
View(part2_data_adza)

 
 

# make a visual of primary endpoint ovet time

# which variables used?
part2_data_adza%>%filter(PARAM %in% "MFM 32 Total (recalculated)")%>%
  select(PARAMCD)%>%distinct()
# what MFMs available?
temp=part2_data_adza%>%
  select(PARAM)%>%distinct()
View(temp)
 
agg=part2_data_adza%>%filter(PARAM %in% "MFM 32 Total (recalculated)")%>%
  group_by(ARM, AVISITN)%>%summarise(ave=mean(AVAL,na.rm = TRUE))%>%ungroup()%>%
  filter(AVISITN>=0 & AVISITN<29)

ind=part2_data_adza%>%filter(PARAM %in% "MFM 32 Total (recalculated)")%>%
  filter(AVISITN>=0 & AVISITN<29)

p1=agg%>%ggplot()+geom_point(aes(x=AVISITN,y=ave,group=ARM,col=ARM))+geom_line(aes(x=AVISITN,y=ave,group=ARM,col=ARM))+
  theme_bw()+scale_color_lancet()+
  scale_x_continuous(breaks=c(0,14,18,21,24,27),labels= c("Baseline","Week 17","Week 35","Week 52","Week 78","Week 104"))+
  theme(legend.position = "bottom")+
  geom_line(data=ind,aes(x=AVISITN,y=AVAL,group=UNI_ID,col=ARM),alpha=0.1)
p2=agg%>%ggplot()+geom_point(aes(x=AVISITN,y=ave,group=ARM,col=ARM))+geom_line(aes(x=AVISITN,y=ave,group=ARM,col=ARM))+
  theme_bw()+scale_color_lancet()+
  scale_x_continuous(breaks=c(0,14,18,21,24,27),labels= c("Baseline","Week 17","Week 35","Week 52","Week 78","Week 104"))+
  theme(legend.position = "bottom") 

grid.arrange(p2, p1, nrow = 1)
#ggsave("~/RCTRWD/xinercode/part2_data_adza.png", width = 30, height = 30, units = "cm")
 

#Part2index
aws.s3::save_object(object ="s3://pddata-euc1-ds-prd-01-fdanihu01rwe/part2/part2_index.rds", file = "part2_index.rds")

part2_index_xlsx <- readRDS("part2_index.rds")

 

## patient char file

 
cols=names(part2_data_adsub)
data.frame(row_number= 1:length(unlist(var_label(part2_data_adsub))),unlist(var_label(part2_data_adsub)))

# variable to label crosswalk
print(part2_data_adsub%>%filter(!PARAMCD %in% "INITSYMP")%>%select(PARAM, PARAMCD)%>%distinct(PARAM, PARAMCD),n=100)


rct.bene<-part2_data_adsub%>%
  select(cols[c(1,2,5,8,9,10,12,13,14,15,16,17,18,19,20,21,27,100,101,102,103,104)] )
View(rct.bene)

# from long to wide
contv=part2_data_adsub%>%filter(!is.na(AVAL) )%>%select(PARAMCD)%>%distinct(PARAMCD)%>%filter(!PARAMCD %in% "INITSYMP")
dat1=rct.bene%>%filter(PARAMCD %in% contv)%>%select(UNI_ID, PARAMCD, AVAL )
dat2=spread(dat1, key = PARAMCD, value =AVAL)
rct.bene<-rct.bene%>%left_join(dat2,by="UNI_ID")

dat1=rct.bene%>%filter(!PARAMCD %in% contv)%>%filter(!PARAMCD %in% "INITSYMP")%>%select(UNI_ID, PARAMCD, AVALC )
dat2=spread(dat1, key = PARAMCD, value =AVALC)
rct.bene<-rct.bene%>%left_join(dat2,by="UNI_ID")

rct.bene<-rct.bene%>%distinct(UNI_ID, .keep_all = TRUE)%>%select(-PARAMCD,-PARAM,-AVAL,-AVALC,-AVALU)

rct.bene<-rct.bene%>%set_variable_labels(AMSTAT="Ambulatory Status",
                                         BBMI="Baseline Body Mass Index (kg/m2)",
                                         BDIABP="Baseline Supine Diastolic Blood Pressure (mmHg)",
                                         BESTRESP="SMA Symptoms:Best",
                                         BHDCIRC="Baseline Head Circumference (cm)",
                                         BHFMSE20="Patient could walk at baseline (Yes/No)",
                                         BHT="Baseline Height (cm)",
                                         BMFM09="Patient could sit at baseline (Yes/No)",
                                         BMFM25="Patient could stand at baseline (Yes/No)",
                                         BRESP="Baseline Supine Respiratory Rate (breaths/min)",
                                         BSYSBP="Baseline Supine Systolic Blood Pressure (mmHg)",
                                         BWT="Baseline Weight (kg)",
                                         DURSRIS="Duration of disease prior to first dose of risdiplam (months)",
                                         DURSTRT="Duration of disease prior to first dose of study medication (months)",
                                         FRACT="Fractures",
                                         HIPSORD="Hip subluxation or dislocation",
                                         HIPSURG="Hip Surgery",
                                         MOTCURR="SMA Motor Function:Current",
                                         MOTHIGH="SMA Motor Function:Best",
                                         RESPDEV1="No Pulmonary Care (Non-Invasive or Invasive)",
                                         RULM="RULM Entry Item A score",
                                         SCOLIO="Scoliosis",
                                         SCOSURG="Surgery for Scoliosis",
                                         SMASAGE="Age at onset of Symptoms",
                                         SMATYPE="SMA Type",
                                         SMN2COGE="SMN2 Copy Number (Genotype)",
                                         SMN2COPY="SMN2 Copy Number (per medical records)",
                                         TRACH="Tracheotomy",
                                         BMFMDSEV="Disease severity",
                                         MFM32="MFM32 Item 9 Score",
                                         RESPDEV2="Cough Assist - Used Daily For Therapy, Not Illness Related",
                                         RESPDEV5="BiPAP Support For More Than 16 Hours Per Day",
                                         SCOCURVE="Scoliosis Grade",
                                         RESPDEV4="BiPAP Support For Less Than 16 Hours Per Day",
                                         RESPDEV3="Cough Assist - Used With An Illness",
                                         RESPDEV6="Airway Clearance Through Cough Assistance")


## longitudinal outcome file

# variable to label crosswalk
print(part2_data_adza%>%select(PARAM, PARAMCD)%>%distinct(PARAM, PARAMCD),n=100)

# cols
cols=names(part2_data_adza)
data.frame(row_number= 1:length(unlist(var_label(part2_data_adza))),unlist(var_label(part2_data_adza)))

 
rct.outcome<-part2_data_adza%>%
  filter(PARAMCD %in% c("MFMT20","MFMT32"))%>%
  filter(AVISITN>=0)%>%
  select(cols[c(1,2,38,124,125,126,127,128,129,130,131,132,133,145,146,154,157)] )%>%
  arrange(UNI_ID,AVISIT)%>%
  mutate(CHG=ifelse(is.na(CHG),0,CHG))


View(rct.outcome)

save(rct.bene, file = "rct.bene.RData")
save(rct.outcome, file = "rct.outcome.RData")
save(ext.bene, file = "ext.bene.RData")
save(ext.outcome, file = "ext.outcome.RData")
 
load("Study/output/rct.bene.RData")
load("Study/output/rct.outcome.RData")
load("Study/output/ext.bene.RData")
load("Study/output/ext.outcome.RData")

