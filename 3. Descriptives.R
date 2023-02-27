pacman::p_load("haven","survey","jtools","remotes","svrepmisc","nhanesA",
               "knitr","tidyverse","plyr","dplyr","foreign","xml2","rvest")
load('analytic.data.Rdata')

####Data type ####
str(analytic.data)
analytic.data[is.na(analytic.data)]<-0
n<-c('LBP','ALD100','RXD300','activity','Com')
for (i in n) {analytic.data[,i]<-as.character(analytic.data[,i])
}
####SURVEY DESIGN####
NHANES <- svydesign(id=~SDMVPSU, weights=~WT,strata=~SDMVSTRA, nest=TRUE, survey.lonely.psu = "adjust", data=analytic.data)
nhc <- subset(NHANES, inAnalysis == "TRUE")#weighted data
uw.data<-subset(analytic.data, inAnalysis == "TRUE" )#unweighted data

plot(svysmooth(~age, design=nhc,bandwidth=1))
plot(svysmooth(~Cotinine, design=nhc,bandwidth=10))
plot(svysmooth(~HbA1c, design=nhc,bandwidth=1))
plot(svysmooth(~BMI, design=nhc,bandwidth=2))
plot(svysmooth(~GGT, design=nhc,bandwidth=2))
svyhist(~GGT, design=nhc)
svyhist(~Cotinine, design=nhc)
svyhist(~age, design=nhc)
####Table 1####
library(compareGroups)
require("tableone")
####weighted data####
colnames(uw.data)
wt.tab1<-svyCreateTableOne(vars = c("gender","age","race","education","marital","income.cut","activity",
                                    "Cotinine","BMI","HbA1c","Com","ALD100","RXD300",
                                    "OHNa_1","OHNa_2","OHFlu_2","OHFlu_3","OHPh_1","OHPh_23","OHP_1","OH_ALL"),
                           strata = "LBP",data = nhc,addOverall=TRUE)
wt.tab1.Con<-print(wt.tab1$ContTable,contDigits = 2, catDigits = 1,
                   pDigits = 3, showAllLevels = TRUE)%>%
  as.data.frame()
wt.tab1.Cat<-print(wt.tab1$CatTable,contDigits = 2, catDigits = 1,
                   pDigits = 3, showAllLevels = TRUE)%>%
  as.data.frame()
wt.tab1.Cat<-separate(data = wt.tab1.Cat, col = No, into = c("s", "No.2"), sep = "\\(")
wt.tab1.Cat<-separate(data = wt.tab1.Cat, col = Yes, into = c("s", "Yes.2"), sep = "\\(")
wt.tab1.Cat<-separate(data = wt.tab1.Cat, col = Overall, into = c("s", "Overall.2"), sep = "\\(")

####unweighted data####
ut.tab1 <- CreateTableOne(vars = c("gender","age","race","education","marital","income.cut","activity",
                                   "Cotinine","BMI","HbA1c","Com","ALD100","RXD300",
                                   "OHNa_1","OHNa_2","OHFlu_2","OHFlu_3","OHPh_1","OHPh_23","OHP_1","OH_ALL"),
                          strata = "LBP",data = uw.data,addOverall=TRUE)
ut.tab1.Cat<-print(ut.tab1$CatTable,contDigits = 2, catDigits = 1,
                   pDigits = 3, showAllLevels = TRUE)%>%
  as.data.frame()
ut.tab1.Cat<-separate(data = ut.tab1.Cat, col = No, into = c("No.1","s"), sep = "\\(")
ut.tab1.Cat<-separate(data = ut.tab1.Cat, col = Yes, into = c("Yes.1","s"), sep = "\\(")
ut.tab1.Cat<-separate(data = ut.tab1.Cat, col = Overall, into = c( "Overall.1","s"), sep = "\\(")
ut.tab1.Cat<-ut.tab1.Cat[,c("level","Overall.1","No.1","Yes.1")]
colnames(ut.tab1.Cat)

#####category data####
tabel1.1<-cbind(wt.tab1.Cat,ut.tab1.Cat)
tabel1.1<-unite(tabel1.1,"No",c("No.1","No.2"), sep="(", remove = F)
tabel1.1<-unite(tabel1.1,"Yes",c("Yes.1","Yes.2"), sep="(", remove = F)
tabel1.1<-unite(tabel1.1,"Overall",c("Overall.1","Overall.2"), sep="(", remove = F)
tabel1.1<-tabel1.1[,c("level","Overall","Yes","No","p")]
colnames(tabel1.1)

#####continuous data####
tabel1.2<-print(wt.tab1$ContTable,contDigits = 2, catDigits = 1,
                pDigits = 3, showAllLevels = TRUE, nonnormal=c("Cotinine","age"))%>%
  as.data.frame()
tabel1.2$level<-""
tabel1.2<-tabel1.2[,c("level","Overall","Yes","No","p")]

####Combine####
Tabel1<-rbind(tabel1.1,tabel1.2)
write.csv(Tabel1,"Tabel1.csv")

