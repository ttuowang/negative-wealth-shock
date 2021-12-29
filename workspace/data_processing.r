# Project: Causal Inference Final Project
#          Negative wealth shock and all-cause mortality
# Description: This R file is for data processing.
# Author: Tuo Wang

# Preparation
setwd("Documents/UW-Madison Fall 2019/STAT 992/992finalproject/data")
library(tidyverse)
library(haven)
library(readr)

# Import these two dataset in Rstudio by using the import dataset button in 
# the environment section.
randhrsimp1992_2016v1 <- read_sas(
  "data/randhrsimp1992_2016v1_SAS/randhrsimp1992_2016v1.sas7bdat", 
  NULL)
randhrs1992_2016v1 <- read_sas(
  "data/randhrs1992_2016v1_SAS/randhrs1992_2016v1.sas7bdat", 
  NULL)
dim(randhrsimp1992_2016v1)
dim(randhrs1992_2016v1)

# Note that the order of randhrs1992_2016v1 and randhrsimp1992_2016v1 are the same
# The unique identifier for each case is HHIDPN.

# First, select all the people in the original HRS cohort, born in 1931 through
# 1941. Use the "cohort" variable to subset the HRS cases. 
# Variable name: HACOHORT, 3: Hrs
table(randhrs1992_2016v1$HACOHORT)

hrs.id <- randhrs1992_2016v1$HHIDPN[randhrs1992_2016v1$HACOHORT == 3]
hrs <- randhrs1992_2016v1[randhrs1992_2016v1$HACOHORT == 3, ]
hrsimp <- randhrsimp1992_2016v1[randhrs1992_2016v1$HACOHORT == 3,]

#------------------------------------------------------------------------------#
# Extract information in the tracker file 
# Set path to the data file "*.DA"
data.file <- "data/tracker file/trk2016/TRK2016TR_R.da"
# Set path to the dictionary file "*.DCT"
dict.file <- "data/tracker file/trk2016/TRK2016TR_R.dct"
# Read the dictionary file
df.dict <- read.table(dict.file, skip = 1, fill = TRUE, stringsAsFactors = FALSE)
# Set column names for dictionary dataframe
colnames(df.dict) <- c("col.num","col.type","col.name","col.width","col.lbl")
# Remove first and last row which only contains a closing }
df.dict <- df.dict[-nrow(df.dict),]
df.dict <- df.dict[-1, ]
# Extract numeric value from column width field
df.dict$col.width <- as.integer(
  sapply(
    df.dict$col.width, 
    gsub, 
    pattern = "[^0-9\\.]", 
    replacement = "")
  )
# Convert column types to format to be used with read_fwf function
df.dict$col.type <- sapply(
  df.dict$col.type, 
  function(x) ifelse(x %in% c("int","byte","long"), 
                     "i", ifelse(x == "float", "n", 
                                 ifelse(x == "double", "d", "c")))
  )
# Read the data file into a dataframe
tracker <- read_fwf(
  file = data.file, 
  fwf_widths(widths = df.dict$col.width, col_names = df.dict$col.name), 
  col_types = paste(df.dict$col.type, collapse = "")
  )
# Add column labels to headers
attributes(tracker)$variable.labels <- df.dict$col.lbl

# Add HHIDPN, HHIDPN = 1000* HHID + PN
tracker$HHIDPN <- as.numeric(tracker$HHID) * 1000 + as.numeric(tracker$PN)

# Save the tracker file
# write.csv(df, file="tracker.csv", row.names = FALSE, col.names = TRUE)
#------------------------------------------------------------------------------#

# Select the subset and arrange by HHIDPN
tracker2 <- filter(tracker, HHIDPN %in% hrs.id) %>% arrange(HHIDPN)

# There is one case that doesn't have tracker file, we can delete it.
# Case HHIDPN=56354031
setdiff(hrsimp$HHIDPN,tracker2$HHIDPN)
hrs <- filter(hrs, HHIDPN!=56354031)
hrsimp <- filter(hrsimp, HHIDPN!=56354031)

# **Important Note:** Check if the three datasets have the same order on HHIDPN
sum(tracker2$HHIDPN==hrs$HHIDPN)
sum(tracker2$HHIDPN==hrsimp$HHIDPN)

# Check the first interview year. We can see most of people (12543) have their first 
# interview in 1992, only 1072 cases have their first interview after 1992. Thus,
# we only include cases that have first interview in 1992 in the analysis.
# also we only include cases that were born between 1931 and 1941
table(tracker2$FIRSTIW)

cond <- (tracker2$FIRSTIW==1992 & tracker2$BIRTHYR >=1931 & tracker2$BIRTHYR<=1941)
hrs.id.1992 <- tracker2$HHIDPN[cond]

tracker.1992 <- tracker2[cond, ]
hrs.1992 <- hrs[cond,]
hrsimp.1992 <- hrsimp[cond,]

# Check ID.
#sum(tracker.1992$HHIDPN == hrs.id.1992)
#sum(hrs.1992$HHIDPN == hrs.id.1992)
#sum(hrsimp.1992$HHIDPN == hrs.id.1992)

#-----------------------#
## Negative Wealth Shock#
#-----------------------#
# The variable description starts from page 50 in randhrsimp1992_2016v1.pdf 

# Section B: Income
# 1. Individual Earnings Page 50
# 'RwIEARN' is the sum of Respondent's wage/salary income, bonuses/overtime pay/
# commissions/tips, 2nd job or military reserve earnings, and professional practice
# or trage income. 

# 2. Household Capital Income Page 79
# 'HwICAP' is the sum of household business or farm income, self-employment earnings, 
# business income, gross rent, dividend and interest income, trust funds or royalties, 
# and other asset income.

# 3. Pension and Annuity Page 132
# 'RwIPENA' is the sum of the Respondent’s income from all pensions and annuities.

# 4. Individual Income from Social Security DI or SSI Page 195
# 'RwISSDI' is the sum of the Respondent’s income from Social Security 
# disability (SDI) and Supplemental Security income (SSI).

# 5. Individual Income from Social Security Retirement Page 213 
# 'RwISRET' is the Respondent’s income from Social Security retirement, 
# spouse or widow benefits.

# 6. Individual Unemployment or Workers Compensation Page 238
# 'RwIUNWC' sums the Respondent’s income from unemployment and worker’s compensation.

# 7. Individual income from other government transfers Page 255
# 'RwIGXFR' sums the Respondent’s income from veterans’ benefits, welfare, and food stamps.

# 8. All other household income Page 288
# 'HwIOTHR' sums alimony, lump sums, and other income received.

# 9. Total household income (Respondent & spouse) Page 323
# 'HwITOT' reflects total income for the last calendar year. 
# HwITOT is set to the sum of Respondent and spouse earnings (RwIEARN, SwIEARN), 
# pensions and annuities (RwIPENA, SwIPENA), SSI and Social Security Disability 
# (RwISSDI, SwISSDI), Social Security retirement (RwISRET, SwISRET), unemployment 
# and workers compensation (RwIUNWC, SwIUNWC), other government transfers 
# (RwIGXFR, SwIGXFR), household capital income (HwICAP), and other income (HwIOTHR).

# Section C: Financial and Housing Wealth
# 1. Net value of real estate
# 'HwARLES'

# 2. Net value of vehicles
# 'HwATRAN'. The reported or imputed net value of vehicles is assigned to 'HwATRAN'

# 3. Net value of businesses
# 'HwABSNS'. The reported or imputed net value of businesses is assigned to HwABSNS

# 4. Net value of IRA, Keogh accounts
# 'HwAIRA.'. The reported or imputed net value of all IRA and Keogh accounts is 
# assigned to HwAIRA.

# 5. Net value of stocks, mutual funds, and investment trusts
# 'HwASTCK'. The reported or imputed net value of stocks and mutual funds is 
# assigned to HwASTCK.

# 6. Value of checking, savings, or money market accounts
# 'HwACHCK'. The reported or imputed value of checking, savings, and money market 
# accounts is assigned to HwACHCK.

# 7. Value of CD, government savings bonds, and T-bills
# 'HwACD'. The reported or imputed value of CDs, government savings bonds, and 
# treasury bills is assigned to HwACD.

# 8. Net value of all other savings, 'HwAOTHR'
# 9. Value of other debt, 'HwADEBT'
# 10. Value of truse, 'HwATRST'. 
# 11. Value of primary residence, 'HwAHOUS'.
# 12. Value of all mortgages/land contracts (primary residence), 'HwAMORT'.
# 13. Value of other home loans (primary residence)
# 14.* Net value of primary residence 'HwATOTH'
# HwAHOUS = value of primary residence
# HwAMORT = value of all mortgages/land contracts (primary residence) 
# HwAHMLN = value of other home loans (primary residence)

# 15. Value of secondary residence, 'HwAHOUB'
# 16. Value of all mortgages/land contracts (secondary residence), 'HwAMRTB'
# 17.* Net value of secondary residence, 'HwANETHB'
# 18.* Net value of non-housing financial wealth, 'HwATOTF'
# Sum (HwASTCK, HwACHCK, HwACD, HwABOND, HwAOTHR) - HwADEBT
# Note: This total does NOT include the value of IRAs and Keogh plans, 
# nor does it include the value of any real estate, vehicles, or businesses.
# 19. Total Wealth (Excluding Secondary Residence), 'HwATOTA'
# 20.* Total Wealth (Including Secondary Residence), 'HwATOTB'
# 21.* Total Wealth (Excluding IRAs), 'HwATOTW'
# 22.* Total Non-housing Wealth, 'H1ATOTN'
load("./data/data_1992.RData")
dim(tracker.1992)
dim(hrs.1992)
dim(hrsimp.1992)

table(is.na(hrsimp.1992$R2IEARN))

#save(tracker.1992,hrs.1992, hrsimp.1992, file = "data_1992.RData")


HwITOT <- c("H1ITOT","H2ITOT","H3ITOT","H4ITOT","H5ITOT",
            "H6ITOT","H7ITOT","H8ITOT","H9ITOT","H10ITOT",
            "H11ITOT","H12ITOT","H13ITOT")

HwITOTA <- c("H1ATOTA","H2ATOTA","H3ATOTA","H4ATOTA","H5ATOTA",
             "H6ATOTA","H7ATOTA","H8ATOTA","H9ATOTA","H10ATOTA",
             "H11ATOTA","H12ATOTA","H13ATOTA")

HwATOTW <- c("H1ATOTW","H2ATOTW","H3ATOTW","H4ATOTW","H5ATOTW",
             "H6ATOTW","H7ATOTW","H8ATOTW","H9ATOTW","H10ATOTW",
             "H11ATOTW","H12ATOTW","H13ATOTW")

HwNET <- c("H1NET","H2NET","H3NET","H4NET","H5NET",
           "H6NET","H7NET","H8NET","H9NET","H10NET",
           "H11NET","H12NET","H13NET")

HwSHOCK <- c("H1SHOCK","H2SHOCK","H3SHOCK","H4SHOCK",
             "H5SHOCK","H6SHOCK","H7SHOCK","H8SHOCK",
             "H9SHOCK","H10SHOCK","H11SHOCK","H12SHOCK")

inflation <- c(1.71, 1.62, 1.53, 1.47, 1.39,
               1.33, 1.27, 1.19, 1.11, 1.10,
               1.05, 1.01, 1.0)

wealth <- dplyr::select(hrsimp.1992,HHIDPN, HwITOTA) %>%
  mutate( H1NET  =  1.71*(H1ATOTA),
          H2NET  =  1.62*(H2ATOTA),
          H3NET  =  1.53*(H3ATOTA),
          H4NET  =  1.47*(H4ATOTA),
          H5NET  =  1.39*(H5ATOTA),
          H6NET  =  1.33*(H6ATOTA),
          H7NET  =  1.27*(H7ATOTA),
          H8NET  =  1.19*(H8ATOTA),
          H9NET  =  1.11*(H9ATOTA),
          H10NET =  1.10*(H10ATOTA),
          H11NET =  1.05*(H11ATOTA),
          H12NET =  1.01*(H12ATOTA),
          H13NET =  1.00* H13ATOTA) %>%
  dplyr::select(HHIDPN,HwNET) %>%
  mutate( H1DIFF = H1NET - H2NET,
          H2DIFF = H2NET - H3NET,
          H3DIFF = H3NET - H4NET,
          H4DIFF = H4NET - H5NET,
          H5DIFF = H5NET - H6NET,
          H6DIFF = H6NET - H7NET,
          H7DIFF = H7NET - H8NET,
          H8DIFF = H8NET - H9NET,
          H9DIFF = H9NET - H10NET,
          H10DIFF = H10NET - H11NET,
          H11DIFF = H11NET - H12NET,
          H12DIFF = H12NET - H13NET) %>%
  mutate( H1SHOCK = ifelse(H1DIFF >= 0.75*H1NET, 1, 0),
          H2SHOCK = ifelse(H2DIFF >= 0.75*H2NET, 1, 0),
          H3SHOCK = ifelse(H3DIFF >= 0.75*H3NET, 1, 0),
          H4SHOCK = ifelse(H4DIFF >= 0.75*H4NET, 1, 0),
          H5SHOCK = ifelse(H5DIFF >= 0.75*H5NET, 1, 0),
          H6SHOCK = ifelse(H6DIFF >= 0.75*H6NET, 1, 0),
          H7SHOCK = ifelse(H7DIFF >= 0.75*H7NET, 1, 0),
          H8SHOCK = ifelse(H8DIFF >= 0.75*H8NET, 1, 0),
          H9SHOCK = ifelse(H9DIFF >= 0.75*H9NET, 1, 0),
          H10SHOCK = ifelse(H10DIFF >= 0.75*H10NET, 1, 0),
          H11SHOCK = ifelse(H11DIFF >= 0.75*H11NET, 1, 0),
          H12SHOCK = ifelse(H12DIFF >= 0.75*H12NET, 1, 0))

wealth$SHOCKTOT <- NA
cond1 <- rowSums(is.na(wealth[,HwSHOCK])) != ncol(wealth[,HwSHOCK])
wealth$SHOCKTOT[cond1] <- rowSums(wealth[cond1,HwSHOCK], na.rm = TRUE)
table(wealth$SHOCKTOT,useNA = "ifany")

wealth$SHOCK <- ifelse(wealth$SHOCKTOT==0,0,1)
#-----------------------#
## All-Cause Mortality  #
#-----------------------#

# In tracker file, there is a variable EXDEATHYR

death <- dplyr::select(tracker.1992, HHIDPN, KNOWNDECEASEDYR) %>%
  mutate(MORTALITY = ifelse(is.na(KNOWNDECEASEDYR), 0, 1)) %>%
  rename( DEATHYR = KNOWNDECEASEDYR) 

death$YEAR <- death$DEATHYR
death$YEAR[is.na(death$DEATHYR)] <- tracker.1992$LASTALIVEYR[is.na(death$DEATHYR)]

wealth.mortality <- dplyr::select(
  cbind.data.frame(wealth,death[,-1]), HHIDPN, 
  SHOCKTOT, SHOCK,
  DEATHYR, MORTALITY, YEAR)

#-----------------------#
## Save final data      #
#-----------------------#
cond2 <- !is.na(wealth.mortality$SHOCK) 

rwshock <- wealth[cond2, ]
wealth.mortality.final <- wealth.mortality[cond2, ]
tracker.final <- tracker.1992[cond2,]
hrs.final <- hrs.1992[cond2, ]
hrsimp.final <- hrsimp.1992[cond2, ]

save(rwshock, wealth.mortality.final,tracker.final,
     hrs.final, hrsimp.final, file = "data_12_01.RData")


write_csv(wealth.mortality.final, path="wealth_mortality.csv")
write_csv(tracker.final, path="tracker_final.csv")



