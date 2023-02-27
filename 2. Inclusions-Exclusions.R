####2. INCLUSIONS/EXCLUSIONS####
pacman::p_load("haven","survey","jtools","remotes","svrepmisc","nhanesA",
               "knitr","tidyverse","plyr","dplyr","foreign","xml2","rvest")
load('raw_data.Rdata')  

#Respondents without sampling weights were excluded.n=15719
step1.data<-subset(data,!is.na(WT) & data$WT!=0)#5442 
ex.1<-nrow(data)-nrow(step1.data)

#age<20 n=2223
step2.data<-subset(step1.data, age >= 20  ) #3219
ex.2<-nrow(step1.data)-nrow(step2.data)

#Respondent with missing data(PAHs) were excluded.#n= 238
step3.data<-subset(step2.data, !is.na(OHNa_1) & !is.na(OHNa_2)& !is.na(OHFlu_2) &
                     !is.na(OHFlu_3) & !is.na(OHPh_1)& !is.na(OHPh_23) & !is.na(OHP_1)& !is.na(LBP))#2981
ex.3<-nrow(step2.data)-nrow(step3.data)

#Respondent with missing data(BMI) were excluded.n=90
step4.data<-subset(step3.data,!is.na(BMI))#2891
ex.4<-nrow(step3.data)-nrow(step4.data)

#Respondent with missing data(HbA1c) were excluded.n=103
step5.data<-subset(step4.data,!is.na(HbA1c))#2788
ex.5<-nrow(step4.data)-nrow(step5.data)

#Respondent with missing data(Cotinine) were excluded.n=45
step6.data<-subset(step5.data,!is.na(Cotinine))#2743
ex.6<-nrow(step5.data)-nrow(step6.data)

#Respondent with missing data(GGT) were excluded.n=9
step7.data<-subset(step6.data, !is.na(GGT))#2735
ex.7<-nrow(step6.data)-nrow(step7.data)
library(DataExplorer)
plot_missing(step7.data)

table(step6.data$OHNa_1_LOD)#0
table(step6.data$OHNa_2_LOD)#0
table(step6.data$OHFlu_2_LOD)#0
table(step6.data$OHFlu_3_LOD)#2
table(step6.data$OHPh_1_LOD)#1
table(step6.data$OHPh_2_LOD)#3
table(step6.data$OHPh_3_LOD)#12
table(step6.data$OHP_1_LOD)#12

COT.N<-which(step6.data$Cotinine == 0.011 | step6.data$Cotinine == 0.035)#600

analytic.data<-mutate(step1.data,inAnalysis= (step1.data$SEQN %in% step6.data$SEQN))
analytic.data <- select(analytic.data,-c("OHNa_1_LOD","OHNa_2_LOD","OHFlu_2_LOD","OHFlu_3_LOD",
                                         "OHPh_1_LOD","OHPh_2_LOD","OHPh_3_LOD","OHP_1_LOD"))
save(analytic.data,file='analytic.data.Rdata') 





