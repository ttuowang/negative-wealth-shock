# Project: Causal Inference Final Project
#          Negative wealth shock and all-cause mortality
# Description: This R file is for data analysis.
# Author: Tuo Wang

# Preparation
library(tidyverse)
library(survival)
library(survminer)
library(grid)
library(ggplot2)
library(gridExtra)
library(scales)
library(ggthemes)

source("./workspace/plotting.R")
load("./data/data_box/data_11_27.RData")
load("./data/data_box/data_12_01.RData")

dim(rwshock)
dim(wealth.mortality.final)
dim(tracker.final)
dim(hrs.final)
dim(hrsimp.final)

names(wealth.mortality.final)

df <- mutate(wealth.mortality.final, 
             STATUS = ifelse(!is.na(tracker.final$KNOWNDECEASEDYR),1,0))

df$TIME <- df$YEAR - 1992

# A: Exploratory Data Analysis

# Chi-square test
table(df$SHOCK, df$MORTALITY)
chisq.test(df$SHOCK,df$MORTALITY)

# Total follow-up year:
sum(df$TIME)

# See the death rate
table(df$STATUS[df$SHOCK==1])
sum(df$STATUS[df$SHOCK==1] ==1)/sum(df$SHOCK==1)
table(df$STATUS[df$SHOCK==0])
sum(df$STATUS[df$SHOCK==0] ==1)/sum(df$SHOCK==0)

# Figure, number of death vs year
deathyear <- data.frame( table(df$DEATHYR))
names(deathyear) <- c("year", "death")
deathyear$year <- as.numeric(as.character(deathyear$year))
pdf("./results/age.pdf",height = 5, width = 6, family = "Times New Roman")
ggplot(deathyear[-24,], aes(x=year, y=death)) + 
  geom_point(size=2,alpha=0.7,color="#386cb0") +
  geom_line(alpha=0.9,color="#386cb0") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(breaks = seq(1992, 2016,2))
dev.off()
# Fit a kaplan meier curve. ( ? )
fit1 <- survfit(Surv(TIME,STATUS)~SHOCK, data=df)
ggsurvplot(fit1, data = df, palette = c("#386cb0","#fdb462"),
           risk.table = TRUE, ggtheme = theme_bw())

# B: covariates ( Time-independet covariates)

# (I). Baseline covariates
# Age at 1992 interview: 'AAGE' from tracker file [No NAs]
# self-reported sex: 'RAGENDER' from hrs file [No NAs]
# self-reported race/ethnicity: 'RARACEM' from hrs file [No NAs]
# Whether Hispanic: 'RAHISPAN' from hrs file [No NAs]
# educational attainment: RAEDYRS from hrs file [No NAs]
# smoking status: R1SMOKEN from hrs, 
#                 RwSMOKEV indicates whether the Respondent ever smoked cigarettes. 
#                 RwSMOKEN indicates whether the Respondent smokes now
# alcohol consumption: R1DRINKR from hrs
# physical activity: R1VIGACT from hrs
# BMI: hrs.final$R1BMI from hrs

# (II). Two indicators measuring financial dispostion at baseline.
# Financial risk aversion: R1RISK from hrs
# leave a bequest: R1BEQLRG from hrs

baseline_x <- select(hrs.final, HHIDPN,RAGENDER,RARACEM,
                     RAHISPAN,RAEDYRS,R1SMOKEN,R1DRINKR,
                     R1VIGACT,R1BMI,R1RISK,R1BEQLRG) 
baseline_x$AAGE <- tracker.final$AAGE

baseline_x <- drop_na(baseline_x)

df2 <- merge(df,baseline_x,by="HHIDPN")
df2 <- mutate(df2, RAGENDER = factor(RAGENDER),
              RARACEM  = factor(RARACEM),
              RAHISPAN = factor(RAHISPAN),
              R1SMOKEN = factor(R1SMOKEN),
              R1DRINKR = factor(R1DRINKR),
              R1VIGACT = factor(R1VIGACT),
              R1RISK   = factor(R1RISK),
              R1BEQLRG = factor(R1BEQLRG))
names(df2)
#rownames(df2) <- df2$HHIDPN
str(df2)
#table(df2$SHOCK, df2$MORTALITY)
#chisq.test(df2$SHOCK,df2$MORTALITY)

#table(df2$STATUS[df2$SHOCK==1])
#sum(df2$STATUS[df2$SHOCK==1] ==1)/sum(df2$SHOCK==1)
#table(df2$STATUS[df2$SHOCK==0])
#sum(df2$STATUS[df2$SHOCK==0] ==1)/sum(df2$SHOCK==0)

p1=ggplot(df2, aes(x=factor(SHOCK),y=R1BMI)) + 
  geom_boxplot(aes(fill = factor(SHOCK) ),alpha=00.7) +
  theme_bw()+
  theme(legend.position='bottom')+
  labs(x = "Negative weatlth shock", fill="Negative weatlth shock",y='BMI')+
  scale_fill_Publication()

p2=ggplot(df2, aes(x=factor(SHOCK),y=RAEDYRS)) + 
  geom_boxplot(aes(fill = factor(SHOCK) ),alpha=00.7) +
  theme_bw()+
  theme(legend.position='bottom')+
  labs(x = "Negative weatlth shock", fill="Negative weatlth shock",
       y='Education attainment')+
  scale_fill_Publication()

p3=ggplot(df2, aes(x=factor(SHOCK),y=AAGE)) + 
  geom_boxplot(aes(fill = factor(SHOCK) ),alpha=00.7) +
  theme_bw()+
  theme(legend.position='bottom')+
  labs(x = "Negative weatlth shock", fill="Negative weatlth shock",y='Age')+
  scale_fill_Publication()

pdf("./results/boxplot.pdf",height = 5, width = 8, family = "Times New Roman")
grid.arrange(p1,p2,p3, ncol=3)
dev.off()
#table(hrs.final$R1RISK, useNA = 'ifany')

# (III). Propensity Score Caliper Matching on the baseline covariates

source("./workspace/mahal_func.R")
library(optmatch)

# We need to cut the dataset into two parts in order to run the pairmatch function.
# The pairmatch function seems cant handle a very large search space.
# I omit the code for matching, please check matching.r for details.

head(matchedPairMat)
head(matchedPairMat2)

matchedpair <- rbind(matchedPairMat,matchedPairMat2)

propscore.model <- glm(SHOCK ~ RAGENDER+RARACEM+RAHISPAN+RAEDYRS+R1SMOKEN+
                         R1DRINKR+R1VIGACT+R1BMI+R1RISK+R1BEQLRG+AAGE,
                       x = TRUE, y = TRUE,
                       family = binomial, data = df2)

X = propscore.model$x[,-1] #remove intercept
A = propscore.model$y
logitps = predict(propscore.model)

X.treatment.before <- X[df2$SHOCK==1, ]
X.control.before <- X[df2$SHOCK==0, ]

X.treatment.after <- X[matchedpair[,1], ]
X.control.after <- X[matchedpair[,2], ]


treatmean.before <- apply(X.treatment.before, 2, mean)
treatvar.before <- apply(X.treatment.before, 2, var)
controlmean.before <- apply(X.control.before, 2, mean)
controlvar.before <- apply(X.control.before, 2, var)


treatmean.after <- apply(X.treatment.after, 2, mean)
treatvar.after <- apply(X.treatment.after, 2, var)
controlmean.after <- apply(X.control.after, 2, mean)
controlvar.after <- apply(X.control.after, 2, var)

stand.diff.before <- (treatmean.before - controlmean.before)/sqrt((treatvar.before+controlvar.before)/2)
stand.diff.after <- (treatmean.after - controlmean.after)/sqrt((treatvar.after+controlvar.after)/2)

standBeforeAfter = cbind(stand.diff.before,stand.diff.after)

colnames(standBeforeAfter ) = c("Before Match (Standardized Diff)",
                                "After Match (Standardized Diff)")
#knitr::kable(round(abs(standBeforeAfter),3), caption = "Differences in Covariates (Before and After)")

abs.stand.diff.before=abs(stand.diff.before)
abs.stand.diff.after=abs(stand.diff.after)
covariates=names(stand.diff.before)
plot.df=data.frame(abs.stand.diff=c(abs.stand.diff.before,abs.stand.diff.after),
                   covariates=rep(covariates,2),
                   type=c(rep("Before",length(covariates)),
                          rep("After",length(covariates))))
fig2 <- ggplot(plot.df,aes(x=abs.stand.diff,y=covariates))+
  geom_point(size=3,aes(color=factor(type)), alpha=0.8)+
  scale_shape_manual(values=c(4,1))+
  geom_vline(xintercept=c(.1,.2),lty=2) +
  scale_colour_Publication()+
  labs(color = "Before or After matching") +
  theme_bw() +
  theme(legend.position="top") 

require(extrafont)  
loadfonts()

png("./results/matchingplot.png",height = 500, width = 600, family = "Times New Roman")
fig2
dev.off()

# (IV). After matching analysis

df.treatment <- df2[matchedpair[,1],]
df.control <- df2[matchedpair[,2],]

mean(df.treatment$MORTALITY)
mean(df.control$MORTALITY)

t.test(df.treatment$MORTALITY, df.control$MORTALITY,paired=TRUE)

save(df.treatment, df.control, file = "./data/data_box/match_pair_12_01.RData")

load("./data/match_pair_12_01.rdata")
library(sensitivitymw)
Ymat.proj = data.frame(df.treatment$STATUS, df.control$STATUS)
# senmw does matching with multiple controls (same number of controls)
GammaSeq = seq(1,2,0.1)
upperBound = matrix(0,length(GammaSeq),3)
for(i in 1:length(GammaSeq)) {
  upperBound[i,1] =senmw(Ymat.proj, gamma = GammaSeq[i], method = "t")$pval
  upperBound[i,2] =senmw(Ymat.proj, gamma = GammaSeq[i], method = "p")$pval
  upperBound[i,3] = senmw(Ymat.proj, gamma = GammaSeq[i], method = "w")$pval
}
out = cbind(GammaSeq,upperBound)
colnames(out) = c("Gamma","t-test","trimmed mean test","weighted trimmed mean test")
round(out,3)
senmwCI(Ymat.proj, gamma = 1.25, method = "w", one.sided = TRUE)
library(sensitivitymv)
amplify(1.25, c(4 : 7)) #Here 4:7 is the \beta_{UA} amplitude
uniroot(function(x){amplify(1.25,x) - x},c(1.25+0.01,10))$root

#df.treatment$TREATED <- rep(1, nrow(df.treatment))
#df.control$TREATED <- rep(0, nrow(df.control))
#df.survival <- rbind.data.frame(df.treatment,df.control)
#fit2 <- survfit(Surv(TIME,STATUS)~TREATED, data=df.survival)
#ggsurvplot(fit2, data = df.survival, palette = c("#386cb0","#fdb462"),
#           risk.table = TRUE, ggtheme = theme_bw())

# (V.) Sensitivity analysis. 


# C: time-varying covariates ( Time-dependet covariates)

# (I). Time-varying covariates
# -marital status: 'RwMSTATH': Page 164, 
# -labor force status: 'RwLBRF': Page 1406
# -health insurance status: '' [...]
# -self report health: 'RwSHLT'
# -whether health limited the ability to work: 'RwHLTHLM'
# -hospitalization in the past 2 years: 'RwHOSP'
# -history of 8 chronic donditions:
# RwHIBPE, RwDIABE, RwCANCRE, RwLUNGE, RwHEARTE, RwSTROKE, RwPSYCHE, and RwARTHRE 
# -activities of daily living: [...]
# The ADLs include walking across a room (RwWALKR), 
# dressing (RwDRESS), bathing (RwBATH), eating (RwEAT), 
# getting in and out of bed (RwBED), and using the toilet (RwTOILT).

RwMSTATH <- paste0("R",c(1:13),"MSTATH") # 9 cats
RwLBRF <- paste0("R",c(1:13),"LBRF")     # 7 cats
RwSHLT <-paste0("R",c(1:13),"SHLT")      # 5 cats
RwHLTHLM <-paste0("R",c(1:13),"HLTHLM")  # 2 cats
RwHOSP <- paste0("R",c(1:13),"HOSP")     # 2 cats
RwHIBPE <- paste0("R",c(1:13),"HIBPE")   # 2 cats
RwDIABE <-paste0("R",c(1:13),"DIABE")    # 2 cats
RwCANCRE <- paste0("R",c(1:13),"CANCRE") # 2 cats
RwLUNGE <- paste0("R",c(1:13),"LUNGE")   # 2 cats
RwHEARTE <- paste0("R",c(1:13),"HEARTE") # 2 cats
RwSTROKE <- paste0("R",c(1:13),"STROKE") # 2 cats
RwPSYCHE <- paste0("R",c(1:13),"PSYCHE") # 2 cats
RwARTHRE <- paste0("R",c(1:13),"ARTHRE") # 2 cats

timevarying_x <- dplyr::select(
  hrs.final,HHIDPN, RwMSTATH,RwLBRF,RwSHLT,
  RwHLTHLM,RwHOSP,RwHIBPE,RwDIABE,RwCANCRE,
  RwLUNGE,RwHEARTE,RwSTROKE,RwPSYCHE,RwARTHRE)

# df2 is the dataframe containing all the baseline covariates without any NAs.

HwSHOCK <- c("H1SHOCK","H2SHOCK","H3SHOCK","H4SHOCK",
             "H5SHOCK","H6SHOCK","H7SHOCK","H8SHOCK",
             "H9SHOCK","H10SHOCK","H11SHOCK","H12SHOCK")

df5 <- merge(df2, timevarying_x,by="HHIDPN")
df5 <- merge(df5, dplyr::select(rwshock,HHIDPN,HwSHOCK), by="HHIDPN")
#df6 <- df5
# Factorize variables
for(i in c(20:32)){
  df5[, i] <- factor(df5[,i], levels = c(1:9))
}
for(i in c(33:45)){
  df5[, i] <- factor(df5[,i], levels = c(1:7))
}
for(i in c(46:58)){
  df5[, i] <- factor(df5[,i], levels = c(1:5))
}
for(i in c(59:188)){
  df5[, i] <- factor(df5[,i])
}
str(df5,list.len=ncol(df5))
#table(df$R1MSTATH, useNA = "ifany")


# (II). Risk set matching.

matchpair <- function(dataf){
  
  rownames(dataf) = c(1:nrow(dataf))
  
  formulaa <- as.formula( paste0(names(dataf)[ncol(dataf)],'~.-HHIDPN' ) )
  
  propscore.model <- glm(formulaa,
                         x = TRUE, y = TRUE,
                         family = binomial, data = dataf)
  X = propscore.model$x[,-1] #remove intercept
  A = propscore.model$y
  logitps = predict(propscore.model)
  distmat=smahal(A,X) 
  rownames(distmat)= rownames(dataf)[A == 1]
  colnames(distmat)= rownames(dataf)[A == 0]
  distmat_caliper=addcaliper(distmat,A,logitps,calipersd = 0.5)
  noControls = 1# Pair match 1 control to each treated
  matchvec=pairmatch(distmat_caliper,controls=noControls,data=df) #The last argument is used to re-order the output; it doesn't actually use the data itself.
  #summary(matchvec)
  #print(matchvec, grouped = TRUE)
  # each column represents who they are and the individual's logit PS value
  matchvec.num = as.numeric(substr(matchvec,start=3,stop=10))
  matchvec.num.notNA = matchvec.num[!is.na(matchvec.num)] #To remove individuals who didn't get matched.
  matchID = unique(matchvec.num.notNA)
  I = length(matchID)
  matchedPairMat = matrix(0,I,4)
  colnames(matchedPairMat) = c("SubjectID (Treated)","SubjectID (Control)","PS (Treated)","PS (Control)")
  treatedSubjID = rownames(dataf)[A==1]
  controlSubjID = rownames(dataf)[A==0]
  for(i in 1:I) {
    subjectIDs = which(matchvec.num == matchID[i])
    matchedPairMat[i,"SubjectID (Treated)"] = subjectIDs[subjectIDs %in% treatedSubjID]
    matchedPairMat[i,"SubjectID (Control)"] = subjectIDs[subjectIDs %in% controlSubjID]
    matchedPairMat[i,"PS (Treated)"] = round(logitps[matchedPairMat[i,"SubjectID (Treated)"]],3)
    matchedPairMat[i,"PS (Control)"] = round(logitps[matchedPairMat[i,"SubjectID (Control)"]],3)
  }
  
  matchedPairMat[,1] <- dataf$HHIDPN[matchedPairMat[,1]]
  matchedPairMat[,2] <- dataf$HHIDPN[matchedPairMat[,2]]
  return(matchedPairMat)
}

# Year 1994, wave 2:

table(df5$H1SHOCK, useNA = "ifany")
baseline_var <- names(df2)[9:19]
df_wave2 <- dplyr::select( df5, HHIDPN,baseline_var,R2MSTATH,R2LBRF,
                           R2SHLT,R2HLTHLM,R2HOSP,R2HIBPE,R2DIABE,R2CANCRE,
                           R2LUNGE,R2HEARTE,R2STROKE,R2PSYCHE,R2ARTHRE,H1SHOCK) %>%
  drop_na()

rownames(df_wave2) = c(1:nrow(df_wave2))
table(df_wave2$H1SHOCK, useNA = "ifany")
matchedPairMat_wave2 <- matchpair(df_wave2)

treatedID <- matchedPairMat_wave2[,1]

# Year 1996, wave 3:

df_wave3 <- dplyr::select( df5, HHIDPN,baseline_var,R3MSTATH,R3LBRF,
                           R3SHLT,R3HLTHLM,R3HOSP,R3HIBPE,R3DIABE,R3CANCRE,
                           R3LUNGE,R3HEARTE,R3STROKE,R3PSYCHE,R3ARTHRE,H2SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave3 <- matchpair(df_wave3)

treatedID <- c(treatedID, matchedPairMat_wave3[,1])

# Year 1998, wave 4:

df_wave4 <- dplyr::select( df5, HHIDPN,baseline_var,R4MSTATH,R4LBRF,
                           R4SHLT,R4HLTHLM,R4HOSP,R4HIBPE,R4DIABE,R4CANCRE,
                           R4LUNGE,R4HEARTE,R4STROKE,R4PSYCHE,R4ARTHRE,H3SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave4 <- matchpair(df_wave4)
treatedID <- c(treatedID, matchedPairMat_wave4[,1])

# Year 2000, wave 5:
df_wave5 <- dplyr::select( df5, HHIDPN,baseline_var,R5MSTATH,R5LBRF,
                           R5SHLT,R5HLTHLM,R5HOSP,R5HIBPE,R5DIABE,R5CANCRE,
                           R5LUNGE,R5HEARTE,R5STROKE,R5PSYCHE,R5ARTHRE,H4SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave5 <- matchpair(df_wave5)
treatedID <- c(treatedID, matchedPairMat_wave5[,1])

# Year 2002, wave 6:
df_wave6 <- dplyr::select( df5, HHIDPN,baseline_var,R6MSTATH,R6LBRF,
                           R6SHLT,R6HLTHLM,R6HOSP,R6HIBPE,R6DIABE,R6CANCRE,
                           R6LUNGE,R6HEARTE,R6STROKE,R6PSYCHE,R6ARTHRE,H5SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave6 <- matchpair(df_wave6)
treatedID <- c(treatedID, matchedPairMat_wave6[,1])

# Year 2004, wave 7:
df_wave7 <- dplyr::select( df5, HHIDPN,baseline_var,R7MSTATH,R7LBRF,
                           R7SHLT,R7HLTHLM,R7HOSP,R7HIBPE,R7DIABE,R7CANCRE,
                           R7LUNGE,R7HEARTE,R7STROKE,R7PSYCHE,R7ARTHRE,H6SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave7 <- matchpair(df_wave7)
treatedID <- c(treatedID, matchedPairMat_wave7[,1])

# Year 2006, wave 8:
df_wave8 <- dplyr::select( df5, HHIDPN,baseline_var,R8MSTATH,R8LBRF,
                           R8SHLT,R8HLTHLM,R8HOSP,R8HIBPE,R8DIABE,R8CANCRE,
                           R8LUNGE,R8HEARTE,R8STROKE,R8PSYCHE,R8ARTHRE,H7SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave8 <- matchpair(df_wave8)
treatedID <- c(treatedID, matchedPairMat_wave8[,1])

# Year 2008, wave 9:
df_wave9 <- dplyr::select( df5, HHIDPN,baseline_var,R9MSTATH,R9LBRF,
                           R9SHLT,R9HLTHLM,R9HOSP,R9HIBPE,R9DIABE,R9CANCRE,
                           R9LUNGE,R9HEARTE,R9STROKE,R9PSYCHE,R9ARTHRE,H8SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave9 <- matchpair(df_wave9)
treatedID <- c(treatedID, matchedPairMat_wave9[,1])

# Year 2010, wave 10:
df_wave10 <- dplyr::select( df5, HHIDPN,baseline_var,R10MSTATH,R10LBRF,
                           R10SHLT,R10HLTHLM,R10HOSP,R10HIBPE,R10DIABE,R10CANCRE,
                           R10LUNGE,R10HEARTE,R10STROKE,R10PSYCHE,R10ARTHRE,H9SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave10 <- matchpair(df_wave10)
treatedID <- c(treatedID, matchedPairMat_wave10[,1])

# Year 2012, wave 11:
df_wave11 <- dplyr::select( df5, HHIDPN,baseline_var,R11MSTATH,R11LBRF,
                            R11SHLT,R11HLTHLM,R11HOSP,R11HIBPE,R11DIABE,R11CANCRE,
                            R11LUNGE,R11HEARTE,R11STROKE,R11PSYCHE,R11ARTHRE,H10SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave11 <- matchpair(df_wave11)
treatedID <- c(treatedID, matchedPairMat_wave11[,1])

# Year 2014, wave 12:
df_wave12 <- dplyr::select( df5, HHIDPN,baseline_var,R12MSTATH,R12LBRF,
                            R12SHLT,R12HLTHLM,R12HOSP,R12HIBPE,R12DIABE,R12CANCRE,
                            R12LUNGE,R12HEARTE,R12STROKE,R12PSYCHE,R12ARTHRE,H11SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave12 <- matchpair(df_wave12)
treatedID <- c(treatedID, matchedPairMat_wave12[,1])

# Year 2016, wave 13:
df_wave13 <- dplyr::select( df5, HHIDPN,baseline_var,R13MSTATH,R13LBRF,
                            R13SHLT,R13HLTHLM,R13HOSP,R13HIBPE,R13DIABE,R13CANCRE,
                            R13LUNGE,R13HEARTE,R13STROKE,R13PSYCHE,R13ARTHRE,H12SHOCK) %>%
  filter( !(HHIDPN %in% treatedID)) %>%
  drop_na()

matchedPairMat_wave13 <- matchpair(df_wave13)
treatedID <- c(treatedID, matchedPairMat_wave13[,1])

matchedPairMat_all <- rbind(matchedPairMat_wave2,matchedPairMat_wave3,
                            matchedPairMat_wave4,matchedPairMat_wave5,
                            matchedPairMat_wave6,matchedPairMat_wave7,
                            matchedPairMat_wave8,matchedPairMat_wave9,
                            matchedPairMat_wave10,matchedPairMat_wave11,
                            matchedPairMat_wave12,matchedPairMat_wave13)

save(matchedPairMat_all, file = "./data/data_box/risksetmatch_mat_all.RData")
a1 <- matchedPairMat_all[,1]
a2 <- matchedPairMat_all[,2]

rsmdf.treatment <- filter(df5, HHIDPN %in% matchedPairMat_all[,1])


rsmdf.control <- filter(df5, HHIDPN %in% matchedPairMat_all[,2])

rsmdf.control <- NULL
for(i in c(1:length(a2))){
  rsmdf.control <- rbind.data.frame(rsmdf.control,df5[df5$HHIDPN==a2[i], ])
}

save(rsmdf.treatment,rsmdf.control, file = "./data/data_box/risksetmatch_df_all.RData")

mean(rsmdf.treatment$MORTALITY)
mean(rsmdf.control$MORTALITY)


t.test(rsmdf.treatment$MORTALITY, rsmdf.control$MORTALITY,paired=TRUE)
#length(rsmdf.treatment$MORTALITY)
#length(rsmdf.control$MORTALITY)

