

                    ############################################
        # ##########            LME medication status           ########## #
                    ############################################
        
                    # Pontus P. Sigray, KI, Stockholm, May 2017

##############
# Preperation 
##############

#Load packages
library(tidyverse)
library(xlsx)
library(rstan)
library(brms)

#Stan options                  
rstan_options(auto_write=T)
options(mc.cores=parallel::detectCores())

#Read in data
dfMaster<-read.xlsx("./DerivedData/All_studies_VT_no_overlap.xlsx",sheetIndex = 1)

#Rm Kenk's failed fit and Bloomfield's 2 MAB pats
dfMaster$HIP.VT[dfMaster$HIP.VT>70]<-NA
dfMaster<-dfMaster[!(dfMaster$Study=="Bloomfield" & dfMaster$ID=="MAB.pat"),]

#Z-score VT within genotype, within Studies
dfModelDat<-dfMaster %>%
group_by(Study,Genotype) %>%
mutate(FC.VT.z=as.numeric(scale(DLPFC.VT)),
   TC.VT.z=as.numeric(scale(TC.VT)),
   HIP.VT.z=as.numeric(scale(HIP.VT)))

#Change Hafizi and Collste drug-naive pats to frug-free
dfModelDat$Med.status<-as.character(dfModelDat$Med.status)
dfModelDat$Med.status[dfModelDat$Med.status=="naive"]<-'free'


#Add medication-status variable differentiating between on-meds pats, drug-free pats and controls
dfModelDat$Med.status.dumPrep<-dfModelDat$Med.status
dfModelDat$Med.status.dumPrep[dfModelDat$HC_pat=="HC"] <- "control"
dfModelDat$Med.status.dumPrep<-as.factor(dfModelDat$Med.status.dumPrep)

#n of each group
table(dfModelDat$Med.status)
table(dfModelDat$Med.status.dumPrep)

#Create dummies for medication effect
dfModelDat$Med.stat.control <- ifelse(dfModelDat$Med.status.dumPrep=="control",1,0)
dfModelDat$Med.stat.free <- ifelse(dfModelDat$Med.status.dumPrep=="free",1,0)
dfModelDat$Med.stat.medicated <- ifelse(dfModelDat$Med.status.dumPrep=="medicated",1,0)

##########
# Priors 
##########

# Uninformative priors over the fixed effects: patient-control difference and medication status (and dummy variables)
prior.FE <- set_prior("normal(0,10)", class = "b") 

# Prior for the "grand" sigma term in the likelihood (brms default now)
prior.RE.sigma <- set_prior("student_t(3,0,10)", class = "sigma") 

#  Prior on freely varying intercept for genotype group  
prior.RE.Intc.gene <-set_prior("cauchy(0,0.707)", class = "sd",group = "Genotype", coef = "Intercept") 

# Prior on freely varying intercept for study group
prior.RE.Intc.study <-set_prior("cauchy(0,0.707)", class = "sd",group = "Study", coef = "Intercept") 

## Prior on freely varying slope for genotype
#prior.RE.slope.gene <- set_prior("cauchy(0,0.707)", class = "sd",group = "Genotype")

## Prior on freely varying slope for study
#prior.RE.slope.study <- set_prior("cauchy(0,0.707)", class = "sd",group = "Study") 

######################################
# Base arguments common for all models
######################################

base_args <- list(
data = dfModelDat,
control = list(adapt_delta = 0.99),
iter = 21000*4, 
warmup = 1000,
thin=1
)

###################
# Run models 
###################

#Granville model: extra predictor for medicated and drug-free (DF pats and controls), controlling for pat-HC status
priors.M1<-c(prior.FE,
       prior.RE.sigma,
       prior.RE.Intc.gene,
       prior.RE.Intc.study)

#FC
args.M1.FC<-c(list(prior = priors.M1),base_args)
args.M1.FC$formula <- FC.VT.z ~ 1 + HC_pat + Med.status + (1|Genotype) + (1|Study)
capture.output(fit1.FC <- do.call(brm, args.M1.FC))
#plot(fit1.FC)
#stanplot(fit1.FC)

#TC
args.M1.TC<-c(list(prior = priors.M1),base_args)
args.M1.TC$formula <- TC.VT.z ~ 1 + HC_pat + Med.status + (1|Genotype) + (1|Study)
capture.output(fit1.TC <- do.call(brm, args.M1.TC))
#plot(fit1.TC)
#stanplot(fit1.TC)

#HIP
args.M1.HIP<-c(list(prior = priors.M1),base_args)
args.M1.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + Med.status + (1|Genotype) + (1|Study)
capture.output(fit1.HIP <- do.call(brm, args.M1.HIP))
#plot(fit1.HIP)


##############
# Save models
##############

save(fit1.FC,fit1.TC,fit1.HIP,file = "./Results/Medication_status/Med_status_models.RData")

########################################
# Plot posteriors of medication effect
########################################

#Mean and credInt for FC
samp.FC<-rstan::extract(fit1.FC$fit)
mean(samp.FC$b_Med.status)
quantile(as.numeric(rstan::extract(fit1.FC$fit)$b_Med.status),probs = c(0.05,0.95))

#Mean and credInt for TC
samp.TC<-rstan::extract(fit1.TC$fit)
mean(samp.TC$b_Med.status)
quantile(as.numeric(rstan::extract(fit1.TC$fit)$b_Med.status),probs = c(0.05,0.95))

#Mean and credInt for HIP
samp.HIP<-rstan::extract(fit1.HIP$fit)
mean(samp.HIP$b_Med.status)
quantile(as.numeric(rstan::extract(fit1.HIP$fit)$b_Med.status),probs = c(0.05,0.95))

