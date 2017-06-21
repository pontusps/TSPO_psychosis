

                        ############################################
            # ##########    LME hypothesis testing - frequentist    ########## #
                        ############################################

                        # Pontus P. Sigray, KI, Stockholm, May 2017

##############
# Preperation 
##############

#Load packages
library(tidyverse)
library(xlsx)
library(lmerTest)

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


##########################################################################
# Frequentist LME models with REML using Satterthwaite approximation of df
##########################################################################

m1.FC <- lmer(data = dfModelDat, formula = (FC.VT.z ~ HC_pat + (1|Genotype) + (1|Study) ))
sum.m1.FC<-summary(m1.FC)
deltaVT.FC<-sum.m1.FC$coefficients[2,]

m1.TC <- lmer(data = dfModelDat, formula = (TC.VT.z ~ HC_pat + (1|Genotype) + (1|Study) ))
sum.m1.TC<-summary(m1.TC)
deltaVT.TC<-sum.m1.TC$coefficients[2,]

m1.HIP <- lmer(data = dfModelDat, formula = (HIP.VT.z ~ HC_pat + (1|Genotype) + (1|Study) ))
sum.m1.HIP<-summary(m1.HIP)
deltaVT.HIP<-sum.m1.HIP$coefficients[2,]


######################
# Write output table 
######################

Regions<-c("FC","TC","HIP")

out.tab <- as.data.frame(rbind(deltaVT.FC,deltaVT.TC,deltaVT.HIP),row.names = Regions)

write.xlsx(out.tab,file = "./Results/Hypothesis_testing_frequentists/Frequentist_models.xlsx")

