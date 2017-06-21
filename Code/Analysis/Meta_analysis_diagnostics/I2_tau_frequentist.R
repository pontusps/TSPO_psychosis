
                ##############################################################
    # ##########   I^2: Diagnostic of study heterogenity using MLE estimates  ########## #
                ##############################################################
    
                         # Pontus P. Sigray, KI, Stockholm, May 2017
                
#Load packages   
library(tidyverse)
library(broom)
library(xlsx)

#Read in data
dfMaster<-read.xlsx("./DerivedData/All_studies_VT_no_overlap.xlsx",sheetIndex = 1)

#Rm Kenk's failed fit and Bloomfield's 2 MAB pats
dfMaster$HIP.VT[dfMaster$HIP.VT>70]<-NA
dfMaster<-dfMaster[!(dfMaster$Study=="Bloomfield" & dfMaster$ID=="MAB.pat"),]

#Z-score VT within genotype, within Studies
dfModelDat <-dfMaster %>%
  group_by(Study,Genotype) %>%
  mutate(FC.VT.z=as.numeric(scale(DLPFC.VT)),
         TC.VT.z=as.numeric(scale(TC.VT)),
         HIP.VT.z=as.numeric(scale(HIP.VT)))

#Make sure HC are 0 and pats are 1 (change order of factor-levels)
dfModelDat$HC_pat<-factor(dfModelDat$HC_pat,levels = c("HC","pat") )


#Estimate frequentist tau (MLE) using model M3 (actually weird to use these "mode" values, since the distribution over tau is so scewed - Does the frequentists know about this??)
#FC
m3.freq.FC<-lmerTest::lmer(data = dfModelDat,formula = (FC.VT.z ~ HC_pat + (1|Genotype) + (HC_pat|Study) ) )
tau.freq.FC<-m3.freq.FC %>%
  broom::tidy() %>%
  select(term,estimate) %>%
  filter(term=="sd_HC_patpat.Study") %>%
  .$estimate

tau.freq.FC<-tau.freq.FC

#TC
m3.freq.TC<-lmerTest::lmer(data = dfModelDat,formula = (TC.VT.z ~ HC_pat + (1|Genotype) + (HC_pat|Study) ) )
tau.freq.TC<-m3.freq.TC %>%
  broom::tidy() %>%
  select(term,estimate) %>%
  filter(term=="sd_HC_patpat.Study") %>%
  .$estimate

tau.freq.TC<-tau.freq.TC

#HIP
m3.freq.HIP<-lmerTest::lmer(data = dfModelDat,formula = (HIP.VT.z ~ HC_pat + (1|Genotype) + (HC_pat|Study) ) )
tau.freq.HIP<-m3.freq.HIP %>%
  broom::tidy() %>%
  select(term,estimate) %>%
  filter(term=="sd_HC_patpat.Study") %>%
  .$estimate

tau.freq.HIP<-tau.freq.HIP

#Center HAB and MAB VT values
dat<-dfMaster %>% 
  group_by(Study,Genotype) %>%
  mutate(HIP.VT.cent = as.numeric(scale((HIP.VT),scale = F)),
         TC.VT.cent = as.numeric(scale((TC.VT),scale = F)),
         FC.VT.cent = as.numeric(scale((DLPFC.VT),scale = F)))

#Z-scores centered VT values for HABs and MABs together
dat<-dat %>% 
  group_by(Study) %>%
  mutate(HIP.VT.cent.z = as.numeric(scale((HIP.VT.cent),scale = T)),
         TC.VT.cent.z = as.numeric(scale((TC.VT.cent),scale = T)),
         FC.VT.cent.z = as.numeric(scale((FC.VT.cent),scale = T)))

#Obtain lm estimates for each study by itself
coefs.FC<-dfModelDat %>% 
  group_by(Study) %>%
  do(tidy(lm(FC.VT.z ~ HC_pat ,data = .))) %>%
  filter( !grepl("Intercept",term) ) %>%
  select(Study,estimate,std.error)  %>%
  mutate(weights=1/std.error^2)

coefs.TC<-dfModelDat %>% 
  group_by(Study) %>%
  do(tidy(lm(TC.VT.z ~ HC_pat ,data = .))) %>%
  filter( !grepl("Intercept",term) ) %>%
  select(Study,estimate,std.error)  %>%
  mutate(weights=1/std.error^2)

coefs.HIP<-dfModelDat %>% 
  group_by(Study) %>%
  do(tidy(lm(HIP.VT.z ~ HC_pat ,data = .))) %>%
  filter( !grepl("Intercept",term) ) %>%
  select(Study,estimate,std.error)  %>%
  mutate(weights=1/std.error^2)

#Calculate I2, based on formula in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4410499/
i2.FC<- round( tau.freq.FC^2/(tau.freq.FC^2 + mean(coefs.FC$std.error^2) ),2)
i2.TC<- round(tau.freq.TC^2/(tau.freq.TC^2 + mean(coefs.TC$std.error^2) ),2)
i2.HIP<- round(tau.freq.HIP^2/(tau.freq.HIP^2 + mean(coefs.HIP$std.error^2) ),2)

df = length(unique(dfMaster$Study))-1

#Calculate Cochrane's Q
Q.FC=-df/(i2.FC-1)
Q.TC=-df/(i2.TC-1)
Q.HIP=-df/(i2.HIP-1)





