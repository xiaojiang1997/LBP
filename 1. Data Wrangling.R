####1. Data Wrangling####
pacman::p_load("haven","survey","jtools","remotes","svrepmisc","nhanesA",
               "knitr","tidyverse","plyr","dplyr","foreign","xml2","rvest")
options(download.file.method = "auto")
## CYCLES 2001-2004  ##
# 1.1 DEMO ----------------------------------------------------------------
DEMO_b <- read.xport('./2001-2002/Demographics/demo_b.XPT')
DEMO_c <- read.xport('./2003-2004/Demographics/demo_c.XPT')
DEMO<-full_join(DEMO_c,DEMO_b)
DEMO<-DEMO[,c("SEQN","SDMVPSU" ,"SDMVSTRA", "RIAGENDR","RIDAGEYR",
              "RIDRETH1","DMDEDUC2","DMDMARTL","INDFMPIR")]
demo_vars<-names(DEMO)
demo<- nhanesA::nhanesTranslate("DEMO_b",demo_vars,data=DEMO)
DEMO <- plyr::rename(demo,c(RIAGENDR="gender",
                            RIDAGEYR="age",
                            RIDRETH1="race",
                            DMDEDUC2="education",
                            DMDMARTL="marital",
                            INDFMPIR="income"))

DEMO$race  <- car::recode(DEMO$race, "c('Other Race - Including Multi-Rac','Other Hispanic')='Other Race'")

DEMO$education <- car::recode(DEMO$education, 
                         "c('Less Than 9th Grade','9-11th Grade (Includes 12th grad','High School Grad/GED or Equivale')='!> High School';
                          c('Some College or AA degree','College Graduate or above')='> High School'; 
                          else = 'Unknown'")

DEMO$marital <- car::recode(DEMO$marital, "c('Married')='Married'; 
                                           else='Other'")
#Other includes widowed, divorced, separated, never married, living with partner, or unknown.

DEMO$income.cut <- cut(DEMO$income, breaks = c(-Inf,1.3,3.5,+Inf),labels = c(1,2,3), right = FALSE)
DEMO$income.cut <- car::recode(DEMO$income.cut,"c('1')='< 130%';
                                                c('2')='130%-349%'; 
                                                c('3')='≥ 350%';
                                                else = 'Unknown'")

# 1.2 Cotinine --------------------------------------------------------------
COT_b <- read.xport('./xpt/L06_b.XPT')
COT_c <- read.xport('./xpt/L06COT_c.XPT')
COT<-full_join(COT_c,COT_b)
COT$Cotinine<-COT$LBXCOT
COT<-COT[,c("SEQN","Cotinine")]

# 1.3 BMI -----------------------------------------------------------------
BMX_b <- read.xport('./xpt/bmx_b.XPT')
BMX_c <- read.xport('./xpt/bmx_c.XPT')
BMX<-full_join(BMX_c,BMX_b)
BMX<-BMX[,c("SEQN","BMXBMI" )]
BMX <- plyr::rename(BMX,c(BMXBMI="BMI"))

# 1.3 HbA1c ---------------------------------------------------------------
L10_b <- read.xport('./xpt/l10_b.XPT')
L10_c <- read.xport('./xpt/l10_c.XPT')
L10<-full_join(L10_b,L10_c)
L10<-L10[,c("SEQN","LBXGH" )]
L10 <- plyr::rename(L10,c(LBXGH="HbA1c"))

# 1.4 Physical demands of usual daily activities -------------------------------------------------------
PAQ_b <- read.xport('./xpt/paq_b.XPT')
PAQ_c <- read.xport('./xpt/paq_c.XPT')
PAQ<-full_join(PAQ_c,PAQ_b)
PAQ$PAQ180   <- car::recode(PAQ$PAQ180,"c('1')='1';
                                        c('2')='2';
                                        c('3')='3';
                                        c('4')='4';
                                        else = 'Unknown'")
PAQ<-PAQ[,c("SEQN","PAQ180" )]
PAQ <- plyr::rename(PAQ,c(PAQ180="activity"))

# 1.5 Comorbidities -------------------------------------------------------
mcq_b <- read.xport('./xpt/mcq_b.XPT')
mcq_c <- read.xport('./xpt/mcq_c.XPT')
mcq<-full_join(mcq_b,mcq_c)
mcq<-mcq[,c("SEQN","MCQ220","MCQ160K","MCQ160B","MCQ160C"
                 ,"MCQ160G","MCQ160E","MCQ160F","MCQ160L")]##cancer,chronic bronchitis,congestive heart failure,coronary heart disease,emphysema,heart attack,liver disease,stroke
mcq[is.na(mcq)]<-0
mcq[,c(2:9)]<-ifelse(mcq[,c(2:9)]==1,1,0)

bpq_b <- read.xport('./xpt/bpq_b.XPT')
bpq_c <- read.xport('./xpt/bpq_c.XPT')
bpq<-full_join(bpq_b,bpq_c)
bpq<-bpq[,c("SEQN","BPQ020","BPQ080")]#hypertension,high cholesterol level
bpq[is.na(bpq)]<-0
bpq[,c(2:3)]<-ifelse(bpq[,c(2:3)]==1,1,0)

Comorbidities<-full_join(mcq,bpq)
Comorbidities[is.na(Comorbidities)]<-0
Comorbidities$Com<-rowSums(Comorbidities[,c(2:ncol(Comorbidities))])
Comorbidities$Com[which(Comorbidities$Com>1)]='>1'
Comorbidities<-Comorbidities[,c("SEQN","Com")]

# 1.6 drinking -------------------------------------------------------------
alq_b <- read.xport('./xpt/alq_b.XPT')
alq_c <- read.xport('./xpt/alq_c.XPT')
alq<-full_join(alq_b,alq_c)
alq<-alq[,c("SEQN","ALD100")]
alq[is.na(alq)]<-0
alq[,2]<-ifelse(alq[,2]==1,1,0)

# 1.7 Analgesic Pain Relievers --------------------------------------------
drug_b <- read.xport('./xpt/rxqana_b.XPT')
drug_c <- read.xport('./xpt/rxqana_c.XPT')
drug<-full_join(drug_b,drug_c)
drug<-drug[,c("SEQN","RXD300")]
drug[is.na(drug)]<-0
drug[,2]<-ifelse(drug[,2]==1,1,0)
drug<-drug[-c(which(duplicated(drug))),]

# 1.8 PAHs ----------------------------------------------------------------
PAH_b <- read.xport('./2001-2002/Laboratory/phpypa_b.XPT')
PAH_b$WT<-1/2*(PAH_b$WTSPH2YR)
PAH_c <- read.xport('./2003-2004/Laboratory/l31pah_c.XPT')
PAH_c$WT<-1/2*(PAH_c$WTSB2YR)
PAH<-full_join(PAH_c,PAH_b)
PAH$OHNa_1<-log(PAH$URXP01/PAH$URXUCR+1)
PAH$OHNa_2<-log(PAH$URXP02/PAH$URXUCR+1)
PAH$OHFlu_2<-log(PAH$URXP04/PAH$URXUCR+1)
PAH$OHFlu_3<-log(PAH$URXP03/PAH$URXUCR+1)
PAH$OHPh_1<-log(PAH$URXP06/PAH$URXUCR+1)
PAH$OHPh_23<-log((PAH$URXP05+PAH$URXP07)/PAH$URXUCR+1)
PAH$OHP_1<-log(PAH$URXP10/PAH$URXUCR+1)
PAH$OH_ALL<-log((PAH$URXP01+PAH$URXP02+PAH$URXP04+PAH$URXP03+PAH$URXP06+PAH$URXP05+PAH$URXP07+
                   PAH$URXP10)/PAH$URXUCR+1)
PAH$OHNa_1_LOD<-PAH$URDP01LC
PAH$OHNa_2_LOD<-PAH$URDP02LC
PAH$OHFlu_2_LOD<-PAH$URDP04LC
PAH$OHFlu_3_LOD<-PAH$URDP03LC
PAH$OHPh_1_LOD<-PAH$URDP06LC
PAH$OHPh_2_LOD<-PAH$URDP05LC
PAH$OHPh_3_LOD<-PAH$URDP07LC
PAH$OHP_1_LOD<-PAH$URDP10LC
PAH<-PAH[,c("SEQN","OHNa_1","OHNa_2","OHFlu_2","OHFlu_3",
            "OHPh_1","OHPh_23","OHP_1","OH_ALL","OHNa_1_LOD",
            "OHNa_2_LOD","OHFlu_2_LOD","OHFlu_3_LOD",
            "OHPh_1_LOD","OHPh_2_LOD","OHPh_3_LOD","OHP_1_LOD","WT")]

# 1.9 LBP -----------------------------------------------------------------
MPQ_b <- read.xport('./2001-2002/Questionnaire/mpq_b.XPT')
MPQ_c <- read.xport('./2003-2004/Questionnaire/mpq_c.XPT')
MPQ<-full_join(MPQ_b,MPQ_c)
MPQ$LBP<-NA
MPQ$LBP[which(MPQ$MPQ070 == 1)]<-"Yes"
MPQ$LBP[which(MPQ$MPQ070 == 2)]<-"No"
MPQ<-MPQ[,c("SEQN","LBP")]

# 1.10 GGT ----------------------------------------------------------------
L40_b <- read.xport('./xpt/l40_b.XPT')
L40_c <- read.xport('./xpt/l40_c.XPT')
L40<-full_join(L40_b,L40_c)
L40<-L40[,c("SEQN","LBXSGTSI" )]
L40 <- plyr::rename(L40,c(LBXSGTSI ="GGT"))

# 1.11 Data Merge ---------------------------------------------------------------
data <- join_all(list(DEMO,PAH,MPQ,COT,BMX,L10,PAQ,Comorbidities,alq,drug,L40), by = "SEQN", type='full')

save(data,file='raw_data.Rdata') 


