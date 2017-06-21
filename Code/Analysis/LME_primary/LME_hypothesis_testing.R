
        
                    ############################################
        # ##########            LME hypothesis testing          ########## #
                    ############################################
        
                    # Pontus P. Sigray, KI, Stockholm, May 2017

##############
# Preperation 
##############

#Load packages
library(tidyverse)
library(xlsx)
library(brms)
library(polspline)
    
#Stan options                  
rstan_options(auto_write=T)
options(mc.cores=parallel::detectCores())

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

##########
# Priors - Model 1 showed marginally better fit compared to the others
##########

# Prior over the fixed effect (deltaVT)
## 1. lower bound set to 0.  
prior.FE.lb <- set_prior("normal(0,0.5)", class = "b", lb = 0) 
prior.FE.lb.Robust.Large <- set_prior("normal(0,0.8)", class = "b", lb = 0) 
prior.FE.lb.Robust.Small <- set_prior("normal(0,0.2)", class = "b", lb = 0) 

## 2. upper bound set to 0. 
prior.FE.ub <- set_prior("normal(0,0.5)", class = "b", ub = 0) 
prior.FE.ub.Robust.Large <- set_prior("normal(0,0.8)", class = "b", ub = 0) 
prior.FE.ub.Robust.Small <- set_prior("normal(0,0.2)", class = "b", ub = 0) 

# Prior for the "grand" sigma term in the likelihood (brms default now)
prior.RE.sigma <- set_prior("student_t(3,0,10)", class = "sigma") 

# Prior on freely varying intercept for genotype group  
prior.RE.Intc.gene <-set_prior("cauchy(0,0.707)", class = "sd",group = "Genotype",coef = "Intercept") 

# Prior on freely varying intercept for study group
prior.RE.Intc.study <-set_prior("cauchy(0,0.707)", class = "sd",group = "Study", coef = "Intercept") 

##################
# Base arguments #
##################

base_args <- list(
  data = dfModelDat,
  control = list(adapt_delta = 0.99),
  iter = 21000*10, 
  warmup = 1000*4,
  thin=5)

###################################################
# Retrieve posteriors for lb (increase in patients)
###################################################

priors.lb<-c(prior.FE.lb, prior.RE.sigma, 
             prior.RE.Intc.gene, prior.RE.Intc.study)
priors.lb.small<-c(prior.FE.lb.Robust.Small, prior.RE.sigma, 
                   prior.RE.Intc.gene, prior.RE.Intc.study)
priors.lb.large<-c(prior.FE.lb.Robust.Large, prior.RE.sigma, 
                   prior.RE.Intc.gene, prior.RE.Intc.study)

### Frontal Cortex
args.lb.FC<-c(list(prior = priors.lb),base_args)
args.lb.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.FC <- do.call(brm, args.lb.FC))
#plot(fit.lb.FC)

#Robustness checks
args.lb.small.FC<-c(list(prior = priors.lb.small),base_args)
args.lb.small.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.small.FC <- do.call(brm, args.lb.small.FC))
#plot(fit.lb.small.FC)

args.lb.large.FC<-c(list(prior = priors.lb.large),base_args)
args.lb.large.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.large.FC <- do.call(brm, args.lb.large.FC))
#plot(fit.lb.large.FC)

### Temporal Cortex
args.lb.TC<-c(list(prior = priors.lb),base_args)
args.lb.TC$formula <- TC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.TC <- do.call(brm, args.lb.TC))
#plot(fit.lb.TC)

#Robustness checks
args.lb.small.TC<-c(list(prior = priors.lb.small),base_args)
args.lb.small.TC$formula <- TC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.small.TC <- do.call(brm, args.lb.small.TC))
#plot(fit.lb.small.TC)

args.lb.large.TC<-c(list(prior = priors.lb.large),base_args)
args.lb.large.TC$formula <- TC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.large.TC <- do.call(brm, args.lb.large.TC))
#plot(fit.lb.large.TC)

### Hippocampus
args.lb.HIP<-c(list(prior = priors.lb),base_args)
args.lb.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.HIP <- do.call(brm, args.lb.HIP))
#plot(fit.lb.HIP)

#Robustness checks
args.lb.small.HIP<-c(list(prior = priors.lb.small),base_args)
args.lb.small.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.small.HIP <- do.call(brm, args.lb.small.HIP))
#plot(fit.lb.small.HIP)

args.lb.large.HIP<-c(list(prior = priors.lb.large),base_args)
args.lb.large.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.lb.large.HIP <- do.call(brm, args.lb.large.HIP))
#plot(fit.lb.large.HIP)

####################################################
# Retrieve posteriors for ub (decrease in patients)
####################################################

priors.ub<-c(prior.FE.ub, prior.RE.sigma, 
             prior.RE.Intc.gene, prior.RE.Intc.study)
priors.ub.small<-c(prior.FE.ub.Robust.Small, prior.RE.sigma, 
                   prior.RE.Intc.gene, prior.RE.Intc.study)
priors.ub.large<-c(prior.FE.ub.Robust.Large, prior.RE.sigma, 
                   prior.RE.Intc.gene, prior.RE.Intc.study)

### Frontal Cortex
args.ub.FC<-c(list(prior = priors.ub),base_args)
args.ub.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.FC <- do.call(brm, args.ub.FC))
#plot(fit.ub.FC)

#Robustness checks
args.ub.small.FC<-c(list(prior = priors.ub.small),base_args)
args.ub.small.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.small.FC <- do.call(brm, args.ub.small.FC))
#plot(fit.ub.small.FC)

args.ub.large.FC<-c(list(prior = priors.ub.large),base_args)
args.ub.large.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.large.FC <- do.call(brm, args.ub.large.FC))
#plot(fit.ub.large.FC)

### Temporal Cortex
args.ub.TC<-c(list(prior = priors.ub),base_args)
args.ub.TC$formula <- TC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.TC <- do.call(brm, args.ub.TC))
#plot(fit.ub.TC)

#Robustness checks
args.ub.small.TC<-c(list(prior = priors.ub.small),base_args)
args.ub.small.TC$formula <- TC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.small.TC <- do.call(brm, args.ub.small.TC))
#plot(fit.ub.small.TC)

args.ub.large.TC<-c(list(prior = priors.ub.large),base_args)
args.ub.large.TC$formula <- TC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.large.TC <- do.call(brm, args.ub.large.TC))
#plot(fit.ub.large.TC)

### Hippocampus
args.ub.HIP<-c(list(prior = priors.ub),base_args)
args.ub.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.HIP <- do.call(brm, args.ub.HIP))
#plot(fit.ub.HIP)

#Robustness checks
args.ub.small.HIP<-c(list(prior = priors.ub.small),base_args)
args.ub.small.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.small.HIP <- do.call(brm, args.ub.small.HIP))
#plot(fit.ub.small.HIP)

args.ub.large.HIP<-c(list(prior = priors.ub.large),base_args)
args.ub.large.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit.ub.large.HIP <- do.call(brm, args.ub.large.HIP))
#plot(fit.ub.large.HIP)

##############
# Save models
##############

save(fit.lb.FC,fit.lb.TC,fit.lb.HIP,
     fit.lb.large.FC,fit.lb.large.TC,fit.lb.large.HIP,
     fit.lb.small.FC,fit.lb.small.TC,fit.lb.small.HIP,
     file = "./Results/Hypothesis_testing/HypothesisTesting_H1_models.RData")

save(fit.ub.FC,fit.ub.TC,fit.ub.HIP,
     fit.ub.large.FC,fit.ub.large.TC,fit.ub.large.HIP,
     fit.ub.small.FC,fit.ub.small.TC,fit.ub.small.HIP,
     file = "./Results/Hypothesis_testing/HypothesisTesting_H2_models.RData")

#####################################################
# Savage-Dickey approximation of BFs for lb and ub
####################################################


#Load truncated normal probability distribution package 
library(truncnorm)

### Lower bound hypothesis testing: Increase in patients v.s. null 
attach(what = "./Results/Hypothesis_testing/HypothesisTesting_H1_models.RData")

#Extract the samples from fits
samp.lb.FC<- rstan::extract(fit.lb.FC$fit)
samp.lb.TC<- rstan::extract(fit.lb.TC$fit)
samp.lb.HIP<- rstan::extract(fit.lb.HIP$fit)

#Extract the samples for the deltaVT parameter
deltaVT.lb.FC <- samp.lb.FC$b_HC_patpat
deltaVT.lb.TC <- samp.lb.TC$b_HC_patpat
deltaVT.lb.HIP <- samp.lb.HIP$b_HC_patpat

#Fit a density curve to the posterior
fit.posterior.lb.FC <- logspline(deltaVT.lb.FC,lbound = 0) #lb set to zero
fit.posterior.lb.TC <- logspline(deltaVT.lb.TC,lbound = 0) #lb set to zero
fit.posterior.lb.HIP <- logspline(deltaVT.lb.HIP,lbound = 0) #lb set to zero

#Height at posterior at point deltaVT=0
posterior.lb.FC <- dlogspline(0, fit.posterior.lb.FC)
posterior.lb.TC <- dlogspline(0, fit.posterior.lb.TC)
posterior.lb.HIP <- dlogspline(0, fit.posterior.lb.HIP)

#Height of order restricted prior at deltaVT = 0
prior.lb <- dnorm(0,0,0.5)*2 #truncnorm::dtruncnorm(x = 0,a = 0,mean = 0,sd = 0.5)          

#Calculate BF01: null over increase in pat
BF01.lb.FC <- posterior.lb.FC/prior.lb
BF01.lb.TC <- posterior.lb.TC/prior.lb
BF01.lb.HIP <- posterior.lb.HIP/prior.lb

#Calculate BF10: increase in pat over null
BF10.lb.FC <- 1/BF01.lb.FC
BF10.lb.TC <- 1/BF01.lb.TC
BF10.lb.HIP <- 1/BF01.lb.HIP

### Upper bound bound hypothesis testing: Decrease in patients v.s. null 
attach(what = "./Results/Hypothesis_testing/HypothesisTesting_H2_models.RData")


#Extract the samples from fits
samp.ub.FC<- rstan::extract(fit.ub.FC$fit)
samp.ub.TC<- rstan::extract(fit.ub.TC$fit)
samp.ub.HIP<- rstan::extract(fit.ub.HIP$fit)

#Extract the samples for the deltaVT parameter
deltaVT.ub.FC <- samp.ub.FC$b_HC_patpat
deltaVT.ub.TC <- samp.ub.TC$b_HC_patpat
deltaVT.ub.HIP <- samp.ub.HIP$b_HC_patpat

#Fit a density curve to the posterior
fit.posterior.ub.FC <- logspline(deltaVT.ub.FC,ubound = 0) #ub set to zero
fit.posterior.ub.TC <- logspline(deltaVT.ub.TC,ubound = 0) #ub set to zero
fit.posterior.ub.HIP <- logspline(deltaVT.ub.HIP,ubound = 0) #ub set to zero

#Height at posterior at point deltaVT=0
posterior.ub.FC <- dlogspline(0, fit.posterior.ub.FC)
posterior.ub.TC <- dlogspline(0, fit.posterior.ub.TC)
posterior.ub.HIP <- dlogspline(0, fit.posterior.ub.HIP)

#Height of order restricted prior at deltaVT = 0
prior.ub <- truncnorm::dtruncnorm(x = 0,b = 0,mean = 0,sd = 0.5)          

#Calculate BF01 
BF02.ub.FC <- posterior.ub.FC/prior.ub
BF02.ub.TC <- posterior.ub.TC/prior.ub
BF02.ub.HIP <- posterior.ub.HIP/prior.ub

#Calculate BF10
BF20.ub.FC <- 1/BF02.ub.FC
BF20.ub.TC <- 1/BF02.ub.TC
BF20.ub.HIP <- 1/BF02.ub.HIP

### Hypothesis testing: Increase in patients v.s. decrease in patients 
BF12.FC <- BF10.lb.FC/BF20.ub.FC
BF12.TC <- BF10.lb.TC/BF20.ub.TC
BF12.HIP <- BF10.lb.HIP/BF20.ub.HIP

BF21.FC <- 1/BF12.FC
BF21.TC <- 1/BF12.TC
BF21.HIP <- 1/BF12.HIP

###################################
# Robustness check of Bayes Factor of increase in patients v.s. null
###################################

#Extract the samples from fits
samp.lb.small.FC<- rstan::extract(fit.lb.small.FC$fit)
samp.lb.small.TC<- rstan::extract(fit.lb.small.TC$fit)
samp.lb.small.HIP<- rstan::extract(fit.lb.small.HIP$fit)

#Extract the samples for the deltaVT parameter
deltaVT.lb.small.FC <- samp.lb.small.FC$b_HC_patpat
deltaVT.lb.small.TC <- samp.lb.small.TC$b_HC_patpat
deltaVT.lb.small.HIP <- samp.lb.small.HIP$b_HC_patpat

#Fit a density curve to the posterior
fit.posterior.lb.small.FC <- logspline(deltaVT.lb.small.FC,lbound = 0) #lb set to zero
fit.posterior.lb.small.TC <- logspline(deltaVT.lb.small.TC,lbound = 0) #lb set to zero
fit.posterior.lb.small.HIP <- logspline(deltaVT.lb.small.HIP,lbound = 0) #lb set to zero

#Height at posterior at point deltaVT=0
posterior.lb.small.FC <- dlogspline(0, fit.posterior.lb.small.FC)
posterior.lb.small.TC <- dlogspline(0, fit.posterior.lb.small.TC)
posterior.lb.small.HIP <- dlogspline(0, fit.posterior.lb.small.HIP)

#Height of order restricted prior at deltaVT = 0
prior.lb.small <- dnorm(0,0,0.2)*2#truncnorm::dtruncnorm(x = 0,b = 0,mean = 0,sd = 0.2)          

#Calculate BF01 
BF01.lb.small.FC <- posterior.lb.small.FC/prior.lb.small
BF01.lb.small.TC <- posterior.lb.small.TC/prior.lb.small
BF01.lb.small.HIP <- posterior.lb.small.HIP/prior.lb.small

#Calculate BF10
BF10.lb.small.FC <- 1/BF01.lb.small.FC
BF10.lb.small.TC <- 1/BF01.lb.small.TC
BF10.lb.small.HIP <- 1/BF01.lb.small.HIP


###################################
# Robustness check of Bayes Factor of decrease in patients v.s. null
###################################

### Large effect (SD=0.8):  

#Extract the samples from fits
samp.ub.large.FC<- rstan::extract(fit.ub.large.FC$fit)
samp.ub.large.TC<- rstan::extract(fit.ub.large.TC$fit)
samp.ub.large.HIP<- rstan::extract(fit.ub.large.HIP$fit)

#Extract the samples for the deltaVT parameter
deltaVT.ub.large.FC <- samp.ub.large.FC$b_HC_patpat
deltaVT.ub.large.TC <- samp.ub.large.TC$b_HC_patpat
deltaVT.ub.large.HIP <- samp.ub.large.HIP$b_HC_patpat

#Fit a density curve to the posterior
fit.posterior.ub.large.FC <- logspline(deltaVT.ub.large.FC,ubound = 0) #ub set to zero
fit.posterior.ub.large.TC <- logspline(deltaVT.ub.large.TC,ubound = 0) #ub set to zero
fit.posterior.ub.large.HIP <- logspline(deltaVT.ub.large.HIP,ubound = 0) #ub set to zero

#Height at posterior at point deltaVT=0
posterior.ub.large.FC <- dlogspline(0, fit.posterior.ub.large.FC)
posterior.ub.large.TC <- dlogspline(0, fit.posterior.ub.large.TC)
posterior.ub.large.HIP <- dlogspline(0, fit.posterior.ub.large.HIP)

#Height of order restricted prior at deltaVT = 0
prior.ub.large <- dnorm(0,0,0.8)*2 #truncnorm::dtruncnorm(x = 0,b = 0,mean = 0,sd = 0.8)          

#Calculate BF01 
BF02.ub.large.FC <- posterior.ub.large.FC/prior.ub.large
BF02.ub.large.TC <- posterior.ub.large.TC/prior.ub.large
BF02.ub.large.HIP <- posterior.ub.large.HIP/prior.ub.large

#Calculate BF10
BF20.ub.large.FC <- 1/BF02.ub.large.FC
BF20.ub.large.TC <- 1/BF02.ub.large.TC
BF20.ub.large.HIP <- 1/BF02.ub.large.HIP

### Small effect (SD=0.2) decrease in patients

#Extract the samples from fits
samp.ub.small.FC<- rstan::extract(fit.ub.small.FC$fit)
samp.ub.small.TC<- rstan::extract(fit.ub.small.TC$fit)
samp.ub.small.HIP<- rstan::extract(fit.ub.small.HIP$fit)

#Extract the samples for the deltaVT parameter
deltaVT.ub.small.FC <- samp.ub.small.FC$b_HC_patpat
deltaVT.ub.small.TC <- samp.ub.small.TC$b_HC_patpat
deltaVT.ub.small.HIP <- samp.ub.small.HIP$b_HC_patpat

#Fit a density curve to the posterior
fit.posterior.ub.small.FC <- logspline(deltaVT.ub.small.FC,ubound = 0) #ub set to zero
fit.posterior.ub.small.TC <- logspline(deltaVT.ub.small.TC,ubound = 0) #ub set to zero
fit.posterior.ub.small.HIP <- logspline(deltaVT.ub.small.HIP,ubound = 0) #ub set to zero

#Height at posterior at point deltaVT=0
posterior.ub.small.FC <- dlogspline(0, fit.posterior.ub.small.FC)
posterior.ub.small.TC <- dlogspline(0, fit.posterior.ub.small.TC)
posterior.ub.small.HIP <- dlogspline(0, fit.posterior.ub.small.HIP)

#Height of order restricted prior at deltaVT = 0
prior.ub.small <- dnorm(0,0,0.2)*2#truncnorm::dtruncnorm(x = 0,b = 0,mean = 0,sd = 0.2)          

#Calculate BF01 
BF02.ub.small.FC <- posterior.ub.small.FC/prior.ub.small
BF02.ub.small.TC <- posterior.ub.small.TC/prior.ub.small
BF02.ub.small.HIP <- posterior.ub.small.HIP/prior.ub.small

#Calculate BF10
BF20.ub.small.FC <- 1/BF02.ub.small.FC
BF20.ub.small.TC <- 1/BF02.ub.small.TC
BF20.ub.small.HIP <- 1/BF02.ub.small.HIP



####################################
# Output table with Bayes Factors
####################################

Regions<-c("FC","TC","HIP")
BF01_increase<-c(BF01.lb.FC,BF01.lb.TC,BF01.lb.HIP)
BF10_increase<-c(BF10.lb.FC,BF10.lb.TC,BF10.lb.HIP)

BF02_decrease<-c(BF02.ub.FC,BF02.ub.TC,BF02.ub.HIP)
BF20_decrease<-c(BF20.ub.FC,BF20.ub.TC,BF20.ub.HIP)

BF12_increase_vs_decrease <- c(BF12.FC,BF12.TC,BF12.HIP)
BF21_decrease_vs_increase <- c(BF21.FC,BF21.TC,BF21.HIP)

BF.tab<-data.frame(Region=Regions,
                   BF01_increase=BF01_increase, BF10_increase=BF10_increase,
                   BF02_decrease=BF02_decrease,BF20_decrease,
                   BF12_increase_vs_decrease=BF12_increase_vs_decrease,
                   BF21_decrease_vs_increase=BF21_decrease_vs_increase)

write.xlsx(x = BF.tab,file = './Results/Hypothesis_testing/Bayes_factors.xlsx',row.names = F)

### Robustness check

Regions<-c("FC","TC","HIP")

BF02_large<-c(BF02.ub.large.FC,BF02.ub.large.TC,BF02.ub.large.HIP)
BF20_large<-c(BF20.ub.large.FC,BF20.ub.large.TC,BF20.ub.large.HIP)

BF02_small<-c(BF02.ub.small.FC,BF02.ub.small.TC,BF02.ub.small.HIP)
BF20_small<-c(BF20.ub.small.FC,BF20.ub.small.TC,BF20.ub.small.HIP)

BF01_small<-c(BF01.lb.small.FC,BF01.lb.small.TC,BF01.lb.small.HIP)
BF10_small<-c(BF10.lb.small.FC,BF10.lb.small.TC,BF10.lb.small.HIP)


BF.tab.Robust<-data.frame(Region=Regions,
                   BF02_large=BF02_large, BF20_large=BF20_large,
                   BF02_small=BF02_small, BF20_small=BF20_small,
                   BF01_small=BF01_small, BF10_small=BF10_small)

write.xlsx(x = BF.tab.Robust,file = './Results/Hypothesis_testing/Bayes_factors_robustCheck.xlsx',row.names = F)



