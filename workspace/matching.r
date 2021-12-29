## Propensity score caliper matching
#-------------------------------------------------------------------------------#
df3 <- df2[1:4000, ]

names(df3)
rownames(df3)

propscore.model <- glm(SHOCK ~ RAGENDER+RARACEM+RAHISPAN+RAEDYRS+R1SMOKEN+
                         R1DRINKR+R1VIGACT+R1BMI+R1RISK+R1BEQLRG+AAGE,
                       x = TRUE, y = TRUE,
                       family = binomial, data = df3)

X = propscore.model$x[,-1] #remove intercept
A = propscore.model$y
logitps = predict(propscore.model)

distmat=smahal(A,X) 
# Label the rows and columns of the distance matrix by the subject id numbers in the data; you can set it something else as well.
rownames(distmat)= rownames(df3)[A == 1]
colnames(distmat)= rownames(df3)[A == 0]
distmat[1:4,1:7]

# Penalize this matrix
distmat_caliper=addcaliper(distmat,A,logitps,calipersd = 0.5)
distmat_caliper[1:4,1:7]

noControls = 1# Pair match 1 control to each treated
matchvec=pairmatch(distmat_caliper,controls=noControls,data=df3) #The last argument is used to re-order the output; it doesn't actually use the data itself.
summary(matchvec)
print(matchvec, grouped = TRUE)
# each column represents who they are and the individual's logit PS value
matchvec.num = as.numeric(substr(matchvec,start=3,stop=10))
matchvec.num.notNA = matchvec.num[!is.na(matchvec.num)] #To remove individuals who didn't get matched.
matchID = unique(matchvec.num.notNA)
I = length(matchID)

matchedPairMat = matrix(0,I,4)
colnames(matchedPairMat) = c("SubjectID (Treated)","SubjectID (Control)","PS (Treated)","PS (Control)")
treatedSubjID = rownames(df3)[A==1]
controlSubjID = rownames(df3)[A==0]


for(i in 1:I) {
  subjectIDs = which(matchvec.num == matchID[i])
  matchedPairMat[i,"SubjectID (Treated)"] = subjectIDs[subjectIDs %in% treatedSubjID]
  matchedPairMat[i,"SubjectID (Control)"] = subjectIDs[subjectIDs %in% controlSubjID]
  matchedPairMat[i,"PS (Treated)"] = round(logitps[matchedPairMat[i,"SubjectID (Treated)"]],3)
  matchedPairMat[i,"PS (Control)"] = round(logitps[matchedPairMat[i,"SubjectID (Control)"]],3)
}

#-------------------------------------------------------------------------------#
df4 <- df2[4001:8409, ]

names(df4)
rownames(df4)

propscore.model <- glm(SHOCK ~ RAGENDER+RARACEM+RAHISPAN+RAEDYRS+R1SMOKEN+
                         R1DRINKR+R1VIGACT+R1BMI+R1RISK+R1BEQLRG+AAGE,
                       x = TRUE, y = TRUE,
                       family = binomial, data = df4)

X = propscore.model$x[,-1] #remove intercept
A = propscore.model$y
logitps = predict(propscore.model)

distmat=smahal(A,X) 
# Label the rows and columns of the distance matrix by the subject id numbers in the data; you can set it something else as well.
rownames(distmat)= rownames(df4)[A == 1]
colnames(distmat)= rownames(df4)[A == 0]
distmat[1:4,1:7]

# Penalize this matrix
distmat_caliper=addcaliper(distmat,A,logitps,calipersd = 0.5)
distmat_caliper[1:4,1:7]

noControls = 1# Pair match 1 control to each treated
matchvec2=pairmatch(distmat_caliper,controls=noControls,data=df4) #The last argument is used to re-order the output; it doesn't actually use the data itself.
summary(matchvec2)
print(matchvec2, grouped = TRUE)
# each column represents who they are and the individual's logit PS value
matchvec.num2 = as.numeric(substr(matchvec2,start=3,stop=10))
matchvec.num.notNA2 = matchvec.num2[!is.na(matchvec.num2)] #To remove individuals who didn't get matched.
matchID2 = unique(matchvec.num.notNA2)
I = length(matchID2)

matchedPairMat2 = matrix(0,I,4)
colnames(matchedPairMat2) = c("SubjectID (Treated)","SubjectID (Control)","PS (Treated)","PS (Control)")

rownames(df4) = as.numeric(rownames(df4)) - 4000
treatedSubjID = rownames(df4)[A==1]
controlSubjID = rownames(df4)[A==0]


for(i in 1:I) {
  subjectIDs = which(matchvec.num2 == matchID2[i])
  matchedPairMat2[i,"SubjectID (Treated)"] = subjectIDs[subjectIDs %in% treatedSubjID]
  matchedPairMat2[i,"SubjectID (Control)"] = subjectIDs[subjectIDs %in% controlSubjID]
  matchedPairMat2[i,"PS (Treated)"] = round(logitps[matchedPairMat2[i,"SubjectID (Treated)"]],3)
  matchedPairMat2[i,"PS (Control)"] = round(logitps[matchedPairMat2[i,"SubjectID (Control)"]],3)
}

class(matchedPairMat2)
matchedPairMat2[,1] <- matchedPairMat2[,1]+4000
matchedPairMat2[,2] <- matchedPairMat2[,2]+4000


save(matchedPairMat,matchedPairMat2, file = "./data/data_box/match_data_12_01.RData")

## Risk set matching

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

