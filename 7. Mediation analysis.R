pacman::p_load("haven","survey","jtools","remotes","svrepmisc","nhanesA",
               "knitr","tidyverse","plyr","dplyr","foreign","xml2","rvest")
load('analytic.data.Rdata')
analytic.data$lgGGT<-log(analytic.data$GGT)
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

####Data hierarchy####
cut.model<-function(x){ q25<-summary(uw.data[,x])[2]
q50<-summary(uw.data[,x])[3]
q75<-summary(uw.data[,x])[5]
Q1<-round(median((analytic.data[,x])[which(analytic.data[,x]<= q25)]),3)
Q2<-round(median((analytic.data[,x])[which(analytic.data[,x] >q25 & analytic.data[,x]<=q50)]),3)
Q3<-round(median((analytic.data[,x])[which(analytic.data[,x] >q50 & analytic.data[,x]<=q75)]),3)
Q4<-round(median((analytic.data[,x])[which(analytic.data[,x] >q75)]),3)
cut(analytic.data[,x],c(-Inf,q25,q50,q75,+Inf),labels=c(Q1,Q2,Q3,Q4))
}
LM2<-c("OHNa_1","OHNa_2","OHFlu_2","OHFlu_3","OHPh_1","OHPh_23","OHP_1","OH_ALL")
for (x in LM2) {analytic.data[,paste0(x,".cut")]<-cut.model(x)
}

####SURVEY DESIGN####
NHANES <- svydesign(id=~SDMVPSU, weights=~WT,strata=~SDMVSTRA, nest=TRUE, survey.lonely.psu = "adjust", data=analytic.data)
nhc <- subset(NHANES, inAnalysis == "TRUE")#weighted data

####Exposed variables and mediation variables####
muvariable.names <- c("OHNa_1","OHNa_1.cut","OHNa_2","OHNa_2.cut","OHFlu_2","OHFlu_2.cut","OHFlu_3","OHFlu_3.cut",
                      "OHPh_1","OHPh_1.cut","OHPh_23","OHPh_23.cut","OHP_1","OHP_1.cut","OH_ALL","OH_ALL.cut")
lg2<-llply(muvariable.names, function(x) {
  svyglm(reformulate(c(x,"gender","age","education","income.cut","Cotinine","BMI","activity","RXD300"),
                     response = as.name("lgGGT")),design = nhc, na.action = na.omit)})
model<- function(x){
  B<-round(as.data.frame(lg2[x][[1]]$coefficients)[-1,],2)
  CI =round(as.data.frame(confint (lg2[x][[1]]))[-1,],2)
  id= rownames(as.data.frame(confint (lg2[x][[1]])))[-1]
  p=round(as.data.frame(summary(lg2[x][[1]])$coefficients[,4])[-1,],3)
  data<-data.frame(ID<-id,B<-B,CI<-CI,p<-p)
}
muuni.n<-c(1:(length(muvariable.names)))
lg2.1 <- lapply(muuni.n, model)

library(plyr)
lg2.2<- ldply(lg2.1,data.frame)
addkuo <- function(x,y){paste0("(",x,"-",y,")")}
lg2.2$B.CI<-paste0(sprintf("%0.2f", lg2.2$B....B),addkuo(lg2.2$X2.5..,lg2.2$X97.5..))
lg2.result<-select(lg2.2,ID="ID....id",B.CI="B.CI",p="p....p")
write.csv(lg2.result,"EM.result.csv")

####Trend test####
m<-c("OHNa_1.cut","OHNa_2.cut","OHFlu_2.cut","OHFlu_3.cut","OHPh_1.cut","OHPh_23.cut","OHP_1.cut", "OH_ALL.cut")
for (i in m) {analytic.data[,i]<-as.numeric(as.character(analytic.data[,i]))
}
NHANES <- svydesign(id=~SDMVPSU, weights=~WT,strata=~SDMVSTRA, nest=TRUE, survey.lonely.psu = "adjust", data=analytic.data)
nhc <- subset(NHANES, inAnalysis == "TRUE")#weighted data
muvariable.names <- c("OHNa_1.cut","OHNa_2.cut","OHFlu_2.cut","OHFlu_3.cut","OHPh_1.cut","OHPh_23.cut","OHP_1.cut", "OH_ALL.cut","LBP")
lg2p<-llply(muvariable.names, function(x) {
  svyglm(reformulate(c(x,"gender","age","education","income.cut","Cotinine","BMI","activity","RXD300"),
                     response = as.name("lgGGT")),design = nhc, na.action = na.omit)})
model<- function(x){
  B<-round(as.data.frame(lg2p[x][[1]]$coefficients)[-1,],2)
  CI =round(as.data.frame(confint (lg2p[x][[1]]))[-1,],2)
  id= rownames(as.data.frame(confint (lg2p[x][[1]])))[-1]
  p=round(as.data.frame(summary(lg2p[x][[1]])$coefficients[,4])[-1,],3)
  data<-data.frame(ID<-id,B<-B,CI<-CI,p<-p)
}
model(1)
muuni.n<-c(1:(length(muvariable.names)))
lg2.1p <- lapply(muuni.n, model)
lg2.2p<- ldply(lg2.1p,data.frame)
addkuo <- function(x,y){paste0("(",x,"-",y,")")}
lg2.2p$B.CI<-paste0(sprintf("%0.2f", lg2.2p$B....B),addkuo(lg2.2p$X2.5..,lg2.2p$X97.5..))
lg2.resultp<-select(lg2.2p,ID="ID....id",B.CI="B.CI",p="p....p")
write.csv(lg2.resultp,"EMP.result.csv")

####mediation variables and outcome variables####
logit1<-svyglm(LBP~lgGGT+gender+age+education+income.cut+Cotinine+BMI+activity+RXD300, family=quasibinomial,design = nhc, na.action = na.omit)
z<-data.frame(OR = exp(coef(logit1)),CI =exp(confint(logit1)), wald = summary(logit1)$coefficients[,3]^2,p=summary(logit1)$coefficients[,4])
z
summary(uw.data$lgGGT)
nhc.lgGGT <- update(nhc, lgGGT.cut=as.factor(cut(lgGGT,c(-Inf,2.639,2.996,3.466,+Inf))))

logit2<-svyglm(LBP~lgGGT.cut+gender+age+education+income.cut+Cotinine+BMI+activity+RXD300, family=quasibinomial,design = nhc.lgGGT, na.action = na.omit)
z<-data.frame(OR = exp(coef(logit2)),CI =exp(confint(logit2)), wald = summary(logit2)$coefficients[,3]^2,p=summary(logit2)$coefficients[,4])
z
nhc.lgGGT <- update(nhc, lgGGT.cut=as.numeric(cut(lgGGT,c(-Inf,2.639,2.996,3.466,+Inf))))
logit3<-svyglm(LBP~lgGGT.cut+gender+age+education+income.cut+Cotinine+BMI+activity+RXD300, family=quasibinomial,design = nhc.lgGGT, na.action = na.omit)
z<-data.frame(OR = exp(coef(logit3)),CI =exp(confint(logit3)), wald = summary(logit3)$coefficients[,3]^2,p=summary(logit3)$coefficients[,4])
z

####Mediation analysis####
library("mediation")
addkuo <- function(x,y){paste0("(",x,"-",y,")")}
##OHNa_1
fitM1 <- svyglm(lgGGT ~ OHNa_1 +gender+age+education+income.cut+Cotinine+BMI+activity+RXD300,
                design = nhc,na.action = na.omit)

fitO1 <- svyglm(LBP ~ OHNa_1 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+activity+RXD300,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHNa_1", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OHNa_1<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)

##OHNa_2
fitM1 <- svyglm(lgGGT ~ OHNa_2 +gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)

fitO1 <- svyglm(LBP ~ OHNa_2 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHNa_2", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OHNa_2<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)

##OHFlu_2
fitM1 <- svyglm(lgGGT ~ OHFlu_2 +gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)

fitO1 <- svyglm(LBP ~ OHFlu_2 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHFlu_2", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OHFlu_2<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)

##OHFlu_3
fitM1 <- svyglm(lgGGT ~ OHFlu_3 +gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)

fitO1 <- svyglm(LBP ~ OHFlu_3 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHFlu_3", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)

med.sum<-summary(med.out1)
OHFlu_3<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)

##OHPh_1
fitM1 <- svyglm(lgGGT ~ OHPh_1 +gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)
fitO1 <- svyglm(LBP ~ OHPh_1 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHPh_1", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OHPh_1<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)
##OHPh_23
fitM1 <- svyglm(lgGGT ~ OHPh_23 +gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)

fitO1 <- svyglm(LBP ~ OHPh_23 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHPh_23", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OHPh_23<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)
##OHP_1
fitM1 <- svyglm(lgGGT ~ OHP_1 +gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)
fitO1 <- svyglm(LBP ~ OHP_1 + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OHP_1", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OHP_1<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)
##OH_ALL
fitM1 <- svyglm(lgGGT ~ OH_ALL+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,na.action = na.omit)

fitO1 <- svyglm(LBP ~ OH_ALL + lgGGT+gender+age+education+income.cut+Cotinine+BMI+RXD300+activity,
                design = nhc,family = quasibinomial)
set.seed(567)
med.out1 <- mediate(fitM1, fitO1, treat = "OH_ALL", mediator = "lgGGT",boot = TRUE,
                    robustSE = TRUE, sims = 1000)
med.sum<-summary(med.out1)
OH_ALL<-c("DI"<-paste0(sprintf("%0.3f", med.sum$z.avg), addkuo(sprintf("%0.3f", med.sum$z.avg.ci[1]),sprintf("%0.3f", med.sum$z.avg.ci[2]))),
          "EI"<-paste0(sprintf("%0.3f", med.sum$d.avg), addkuo(sprintf("%0.3f", med.sum$d.avg.ci[1]),sprintf("%0.3f", med.sum$d.avg.ci[2]))),
          "pro"<-sprintf("%0.2f", 100*med.sum$n.avg),
          "p"<-med.sum$d.avg.p)


med<-t(data.frame("OHNa_1"=OHNa_1,"OHNa_2"=OHNa_2,"OHFlu_2"=OHFlu_2,
                  "OHFlu_3"=OHFlu_3,"OHPh_1"=OHPh_1,"OHPh_23"=OHPh_23,"OHP_1"=OHP_1,"OH_ALL"=OH_ALL))

write.csv(med,"med.csv")



