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
uw.data<-subset(analytic.data, inAnalysis == "TRUE" )#unweighted data

####Model 1####
variable.names <- c("OHNa_1","OHNa_1.cut","OHNa_2","OHNa_2.cut","OHFlu_2","OHFlu_2.cut","OHFlu_3","OHFlu_3.cut",
                    "OHPh_1","OHPh_1.cut","OHPh_23","OHPh_23.cut","OHP_1","OHP_1.cut","OH_ALL","OH_ALL.cut","LBP")
lg0<-llply(variable.names, function(x) {
  svyglm(reformulate(x,response = 'LBP'), family=quasibinomial,design = nhc, na.action = na.omit)})
unimodel<- function(x){
        OR<-round(as.data.frame(exp(coef(lg0[[x]])))[-1,],2)
        CI =round(as.data.frame(exp(confint(lg0[[x]])))[-1,],2)
        id= rownames(as.data.frame(exp(coef(lg0[[x]]))))[-1]
        wald = round(as.data.frame(summary(lg0[[x]])$coefficients[,3]^2)[-1,],3)
        p=round(as.data.frame(summary(lg0[[x]])$coefficients[,4])[-1,],3)
        data<-data.frame(ID<-id,OR<-OR,CI<-CI,wald<-wald,p<-p)
      }
uni.n<-c(1:(length(variable.names)-1))
lg0.1 <- lapply(uni.n, unimodel)
lg0.2<- ldply(lg0.1,data.frame)
addkuo <- function(x,y){paste0(" (",x,"-",y,")")}
lg0.2$OR.CI<-paste0(sprintf("%0.2f", lg0.2$OR....OR),addkuo(lg0.2$X2.5..,lg0.2$X97.5..))
lg0.result<-select(lg0.2,ID="ID....id",OR.CI="OR.CI",wald="wald....wald",p="p....p")
write.csv(lg0.result,"lg0.result.csv")

####Model 2#####
muvariable.names <- c("OHNa_1","OHNa_1.cut","OHNa_2","OHNa_2.cut","OHFlu_2","OHFlu_2.cut","OHFlu_3","OHFlu_3.cut",
                      "OHPh_1","OHPh_1.cut","OHPh_23","OHPh_23.cut","OHP_1","OHP_1.cut","OH_ALL","OH_ALL.cut","LBP")
lg1<-llply(muvariable.names, function(x) {
  svyglm(reformulate(c(x,"gender","age","BMI"),
                     response = as.name("LBP")), family=quasibinomial,design = nhc, na.action = na.omit)})
model<- function(x){
  OR<-round(as.data.frame(exp(coef(lg1[[x]])))[-1,],2)
  CI =round(as.data.frame(exp(confint(lg1[[x]])))[-1,],2)
  id= rownames(as.data.frame(exp(coef(lg1[[x]]))))[-1]
  wald = round(as.data.frame(summary(lg1[[x]])$coefficients[,3]^2)[-1,],3)
  p=round(as.data.frame(summary(lg1[[x]])$coefficients[,4])[-1,],3)
  data<-data.frame(ID<-id,OR<-OR,CI<-CI,wald<-wald,p<-p)
}

muuni.n<-c(1:(length(muvariable.names)-1))
lg1.1 <- lapply(muuni.n, model)

lg1.2<- ldply(lg1.1,data.frame)
addkuo <- function(x,y){paste0(" (",x,"-",y,")")}
lg1.2$OR.CI<-paste0(sprintf("%0.2f", lg1.2$OR....OR),addkuo(lg1.2$X2.5..,lg1.2$X97.5..))
lg1.result<-select(lg1.2,ID="ID....id",OR.CI="OR.CI",wald="wald....wald",p="p....p")
write.csv(lg1.result,"lg1.result.csv")

####Model 3####
muvariable.names <- c("OHNa_1","OHNa_1.cut","OHNa_2","OHNa_2.cut","OHFlu_2","OHFlu_2.cut","OHFlu_3","OHFlu_3.cut",
                      "OHPh_1","OHPh_1.cut","OHPh_23","OHPh_23.cut","OHP_1","OHP_1.cut","OH_ALL","OH_ALL.cut","LBP")
lg2<-llply(muvariable.names, function(x) {
  svyglm(reformulate(c(x,"gender","age","education","income.cut","Cotinine","BMI","activity","RXD300"),
                     response = as.name("LBP")), family=quasibinomial,design = nhc, na.action = na.omit)})
model<- function(x){
    OR<-round(as.data.frame(exp(coef(lg2[[x]])))[-1,],2)
    CI =round(as.data.frame(exp(confint(lg2[[x]])))[-1,],2)
    id= rownames(as.data.frame(exp(coef(lg2[[x]]))))[-1]
    wald = round(as.data.frame(summary(lg2[[x]])$coefficients[,3]^2)[-1,],3)
    p=round(as.data.frame(summary(lg2[[x]])$coefficients[,4])[-1,],3)
    data<-data.frame(ID<-id,OR<-OR,CI<-CI,wald<-wald,p<-p)
  }
muuni.n<-c(1:(length(muvariable.names)-1))
lg2.1 <- lapply(muuni.n, model)
lg2.2<- ldply(lg2.1,data.frame)
addkuo <- function(x,y){paste0(" (",x,"-",y,")")}
lg2.2$OR.CI<-paste0(sprintf("%0.2f", lg2.2$OR....OR),addkuo(lg2.2$X2.5..,lg2.2$X97.5..))
lg2.result<-select(lg2.2,ID="ID....id",OR.CI="OR.CI",wald="wald....wald",p="p....p")
write.csv(lg2.result,"lg2.result.csv")




