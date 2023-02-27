pacman::p_load("haven","survey","jtools","remotes","svrepmisc","nhanesA",
               "knitr","tidyverse","plyr","dplyr","foreign","xml2","rvest")
load('analytic.data.Rdata')
uw.data<-subset(analytic.data, inAnalysis == "TRUE" )#unweighted data
####Data type ####
str(analytic.data)
analytic.data[is.na(analytic.data)]<-0
n<-c('LBP','ALD100','RXD300','activity','Com')
for (i in n) {analytic.data[,i]<-as.factor(analytic.data[,i])
}
analytic.data$race<-relevel(analytic.data$race, ref = "Other Race")
analytic.data$gender<-relevel(analytic.data$gender, ref = "Female")
analytic.data$education<-relevel(analytic.data$education, ref = "!> High School")
analytic.data$marital<-relevel(analytic.data$marital, ref = "Other")
analytic.data$income.cut<-relevel(analytic.data$income.cut, ref = "< 130%")
analytic.data$LBP<-relevel(analytic.data$LBP, ref = "No")
analytic.data$activity<-relevel(analytic.data$activity, ref = "1")
analytic.data$Com<-relevel(analytic.data$Com, ref = "0")
analytic.data$ALD100<-relevel(analytic.data$ALD100, ref = "0")
analytic.data$RXD300<-relevel(analytic.data$RXD300, ref = "0")


####SURVEY DESIGN####
NHANES <- svydesign(id=~SDMVPSU, weights=~WT,strata=~SDMVSTRA, nest=TRUE, survey.lonely.psu = "adjust", data=analytic.data)
nhc <- subset(NHANES, inAnalysis == "TRUE")#weighted data
nhc.age <- update(nhc, age.cut=cut(age,c(-Inf,40,60,+Inf)))
nhc.BMI <- update(nhc, BMI.cut=cut(BMI,c(-Inf,20,25,30,+Inf)))
summary(uw.data[,'Cotinine'])#0.020    0.076   57.558   22.150
nhc.Cotinine<-update(nhc,Cotinine.cut=cut(Cotinine,c(-Inf,0.02,0.076,22.15,+Inf)))

#OHNa_1
logit1<-svyglm(LBP~OHNa_1*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHNa_1.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_1.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHNa_1.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_1.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_1.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHNa_1.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_1.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_1*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_1+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_1.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHNa_1.p<-c(OHNa_1.age,OHNa_1.gender,OHNa_1.BMI,OHNa_1.edu,OHNa_1.income,OHNa_1.Cot,OHNa_1.act,OHNa_1.drug)

#OHNa_2
logit1<-svyglm(LBP~OHNa_2*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHNa_2.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_2.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHNa_2.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_2.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_2.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHNa_2.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_2.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHNa_2*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHNa_2+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHNa_2.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHNa_2.p<-c(OHNa_2.age,OHNa_2.gender,OHNa_2.BMI,OHNa_2.edu,OHNa_2.income,OHNa_2.Cot,OHNa_2.act,OHNa_2.drug)

#OHFlu_2
logit1<-svyglm(LBP~OHFlu_2*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHFlu_2.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_2.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHFlu_2.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_2.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_2.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHFlu_2.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_2.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_2*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_2+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_2.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHFlu_2.p<-c(OHFlu_2.age,OHFlu_2.gender,OHFlu_2.BMI,OHFlu_2.edu,OHFlu_2.income,OHFlu_2.Cot,OHFlu_2.act,OHFlu_2.drug)

#OHFlu_3
logit1<-svyglm(LBP~OHFlu_3*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHFlu_3.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_3.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHFlu_3.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_3.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_3.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHFlu_3.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_3.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHFlu_3*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHFlu_3+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHFlu_3.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHFlu_3.p<-c(OHFlu_3.age,OHFlu_3.gender,OHFlu_3.BMI,OHFlu_3.edu,OHFlu_3.income,OHFlu_3.Cot,OHFlu_3.act,OHFlu_3.drug)

#OHPh_1
logit1<-svyglm(LBP~OHPh_1*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHPh_1.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_1.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHPh_1.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_1.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_1.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHPh_1.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_1.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_1*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_1+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_1.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHPh_1.p<-c(OHPh_1.age,OHPh_1.gender,OHPh_1.BMI,OHPh_1.edu,OHPh_1.income,OHPh_1.Cot,OHPh_1.act,OHPh_1.drug)

#OHPh_23
logit1<-svyglm(LBP~OHPh_23*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHPh_23.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_23.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHPh_23.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_23.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_23.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHPh_23.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_23.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHPh_23*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHPh_23+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHPh_23.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHPh_23.p<-c(OHPh_23.age,OHPh_23.gender,OHPh_23.BMI,OHPh_23.edu,OHPh_23.income,OHPh_23.Cot,OHPh_23.act,OHPh_23.drug)

#OHP_1
logit1<-svyglm(LBP~OHP_1*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OHP_1.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OHP_1.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OHP_1.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+education, family=quasibinomial,design =nhc, na.action = na.omit)
OHP_1.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OHP_1.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OHP_1.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OHP_1.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OHP_1*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OHP_1+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OHP_1.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OHP_1.p<-c(OHP_1.age,OHP_1.gender,OHP_1.BMI,OHP_1.edu,OHP_1.income,OHP_1.Cot,OHP_1.act,OHP_1.drug)

#OH_ALL
logit1<-svyglm(LBP~OH_ALL*age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+age.cut, family=quasibinomial,design =nhc.age, na.action = na.omit)
OH_ALL.age<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*gender, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+gender, family=quasibinomial,design =nhc, na.action = na.omit)
OH_ALL.gender<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+BMI.cut, family=quasibinomial,design =nhc.BMI, na.action = na.omit)
OH_ALL.BMI<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*education, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+education, family=quasibinomial,design =nhc, na.action = na.omit)
OH_ALL.edu<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+income.cut, family=quasibinomial,design =nhc, na.action = na.omit)
OH_ALL.income<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+Cotinine.cut, family=quasibinomial,design =nhc.Cotinine, na.action = na.omit)
OH_ALL.Cot<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*activity, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+activity, family=quasibinomial,design =nhc, na.action = na.omit)
OH_ALL.act<-sprintf("%0.3f", anova(logit1,logit2)$p)

logit1<-svyglm(LBP~OH_ALL*RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
logit2<-svyglm(LBP~OH_ALL+RXD300, family=quasibinomial,design =nhc, na.action = na.omit)
OH_ALL.drug<-sprintf("%0.3f", anova(logit1,logit2)$p)

OH_ALL.p<-c(OH_ALL.age,OH_ALL.gender,OH_ALL.BMI,OH_ALL.edu,OH_ALL.income,OH_ALL.Cot,OH_ALL.act,OH_ALL.drug)

inter.p<-data.frame(OHNa_1.p,OHNa_2.p,OHFlu_2.p,OHFlu_3.p,OHPh_1.p,OHPh_23.p,OHP_1.p,OH_ALL.p)
write.csv(inter.p,"inter.p.csv")
#"OHNa_1","OHNa_2","OHFlu_2","OHFlu_3","OHPh_1","OHPh_23","OHP_1","OH_ALL","LBP"


####
nhc.income0<-subset(nhc,income.cut=="0")
nhc.income1<-subset(nhc,income.cut=="< 130%")
nhc.income2<-subset(nhc,income.cut=="130%-349%")
nhc.income3<-subset(nhc,income.cut=="≥ 350%")

logit1<-svyglm(LBP~OHPh_23+age+BMI+gender+education+Cotinine+
                activity+RXD300, family=quasibinomial,design =nhc.income1, na.action = na.omit)
r1<-data.frame(OR = exp(coef(logit1)),CI =exp(confint(logit1)), wald = summary(logit1)$coefficients[,3]^2,p=summary(logit1)$coefficients[,4])
logit2<-svyglm(LBP~OHPh_23+age+BMI+gender+education+Cotinine+
                 activity+RXD300, family=quasibinomial,design =nhc.income2, na.action = na.omit)
r2<-data.frame(OR = exp(coef(logit2)),CI =exp(confint(logit2)), wald = summary(logit2)$coefficients[,3]^2,p=summary(logit2)$coefficients[,4])
logit3<-svyglm(LBP~OHPh_23+age+BMI+gender+education+Cotinine+
                 activity+RXD300, family=quasibinomial,design =nhc.income3, na.action = na.omit)
r3<-data.frame(OR = exp(coef(logit3)),CI =exp(confint(logit3)), wald = summary(logit3)$coefficients[,3]^2,p=summary(logit3)$coefficients[,4])

income1<-c("< 130%",sprintf("%0.2f", r1$OR[2]),sprintf("%0.2f",r1$CI.2.5..[2]),sprintf("%0.2f",r1$CI.97.5..[2]),sprintf("%0.3f",r1$p[2]))
income2<-c("130%-349%",sprintf("%0.2f", r2$OR[2]),sprintf("%0.2f",r2$CI.2.5..[2]),sprintf("%0.2f",r2$CI.97.5..[2]),sprintf("%0.3f",r2$p[2]))
income3<-c("≥ 350%",sprintf("%0.2f", r3$OR[2]),sprintf("%0.2f",r3$CI.2.5..[2]),sprintf("%0.2f",r3$CI.97.5..[2]),sprintf("%0.3f",r3$p[2]))
data.HR<-t(data.frame(income1,income2,income3))
colnames(data.HR)<-c("income","HR","LOW","HIGH","P")
write.csv(data.HR,"HR.csv")
