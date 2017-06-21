

            ###############################################
# ##########  LME hypothesis testing - gender as covariate ########## #
            ###############################################

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


###################
# Priors - Model 1 
###################

# Prior over the fixed effects (deltaVT and gender): no truncation - lb and ub will have to be extracted from the STAN fits 
prior.FE <- set_prior("normal(0,0.5)", class = "b") 

# Prior for the "grand" sigma term in the likelihood (brms default now)
prior.RE.sigma <- set_prior("student_t(3,0,10)", class = "sigma") 

# Prior on freely varying intercept for genotype group  
prior.RE.Intc.gene <-set_prior("cauchy(0,0.707)", class = "sd",group = "Genotype",coef = "Intercept") 

# Prior on freely varying intercept for study group
prior.RE.Intc.study <-set_prior("cauchy(0,0.707)", class = "sd",group = "Study", coef = "Intercept") 

##################
# Base arguments #
##################

#Run many iters since the posterior from will be truncated (ub and lb) after MCMC (brms does not allow truncation of only one of the predictors)

base_args <- list(
  data = dfModelDat,
  control = list(adapt_delta = 0.99),
  iter = 21000*10, 
  warmup = 1000*2,
  thin=2)

###################################################
# Retrieve posteriors for non-truncated hypothesis 
###################################################

priors<-c(prior.FE, prior.RE.sigma, 
             prior.RE.Intc.gene, prior.RE.Intc.study)

### Frontal Cortex
args.FC<-c(list(prior = priors),base_args)
args.FC$formula <- FC.VT.z ~ 1 + HC_pat + Sex + (1|Genotype) + (1|Study)
capture.output(fit.FC <- do.call(brm, args.FC))
#plot(fit.FC)

### Temporal Cortex
args.TC<-c(list(prior = priors),base_args)
args.TC$formula <- TC.VT.z ~ 1 + HC_pat + Sex + (1|Genotype) + (1|Study)
capture.output(fit.TC <- do.call(brm, args.TC))
#plot(fit.TC)

### Hippocampus
args.HIP<-c(list(prior = priors),base_args)
args.HIP$formula <- HIP.VT.z ~ 1 + HC_pat + Sex + (1|Genotype) + (1|Study)
capture.output(fit.HIP <- do.call(brm, args.HIP))
#plot(fit.HIP)

##############
# Save models
##############

save(fit.FC,fit.TC,fit.HIP,
     file = "./Results/Hypothesis_testing/HypothesisTesting_H1_H2_models_gender.RData")

#####################################################
# Savage-Dickey approximation of BFs for lb and ub (extract lb and ub from the same non-truncated model)
####################################################

### Lower bound hypothesis testing: Increase in patients v.s. null 
attach(what = "./Results/Hypothesis_testing/HypothesisTesting_H1_H2_models_gender.RData")

#Extract the samples from fits
samp.FC<- rstan::extract(fit.FC$fit)
samp.TC<- rstan::extract(fit.TC$fit)
samp.HIP<- rstan::extract(fit.HIP$fit)

#Extract the samples for the deltaVT parameter
deltaVT.lb.FC <- samp.FC$b_HC_patpat[samp.FC$b_HC_patpat>0]
deltaVT.lb.TC <- samp.TC$b_HC_patpat[samp.TC$b_HC_patpat>0]
deltaVT.lb.HIP <- samp.HIP$b_HC_patpat[samp.HIP$b_HC_patpat>0]

#Fit a density curve to the posterior
fit.posterior.lb.FC <- logspline(deltaVT.lb.FC,lbound = 0) #lb set to zero
fit.posterior.lb.TC <- logspline(deltaVT.lb.TC,lbound = 0) #lb set to zero
fit.posterior.lb.HIP <- logspline(deltaVT.lb.HIP,lbound = 0) #lb set to zero

#Height at posterior at point deltaVT=0
posterior.lb.FC <- dlogspline(0, fit.posterior.lb.FC)
posterior.lb.TC <- dlogspline(0, fit.posterior.lb.TC)
posterior.lb.HIP <- dlogspline(0, fit.posterior.lb.HIP)

#Height of order restricted prior at deltaVT = 0
prior.lb <- dnorm(0,0,0.5)*2          

#Calculate BF01: null over increase in pat
BF01.lb.FC <- posterior.lb.FC/prior.lb
BF01.lb.TC <- posterior.lb.TC/prior.lb
BF01.lb.HIP <- posterior.lb.HIP/prior.lb

#Calculate BF10: increase in pat over null
BF10.lb.FC <- 1/BF01.lb.FC
BF10.lb.TC <- 1/BF01.lb.TC
BF10.lb.HIP <- 1/BF01.lb.HIP

### Upper bound bound hypothesis testing: Decrease in patients v.s. null 

#Extract the samples from fits

#Extract the samples for the deltaVT parameter
deltaVT.ub.FC <- samp.FC$b_HC_patpat[samp.FC$b_HC_patpat<0]
deltaVT.ub.TC <- samp.TC$b_HC_patpat[samp.TC$b_HC_patpat<0]
deltaVT.ub.HIP <- samp.HIP$b_HC_patpat[samp.HIP$b_HC_patpat<0]

#Fit a density curve to the posterior
fit.posterior.ub.FC <- logspline(deltaVT.ub.FC,ubound = 0) #ub set to zero
fit.posterior.ub.TC <- logspline(deltaVT.ub.TC,ubound = 0) #ub set to zero
fit.posterior.ub.HIP <- logspline(deltaVT.ub.HIP,ubound = 0) #ub set to zero

#Height at posterior at point deltaVT=0
posterior.ub.FC <- dlogspline(0, fit.posterior.ub.FC)
posterior.ub.TC <- dlogspline(0, fit.posterior.ub.TC)
posterior.ub.HIP <- dlogspline(0, fit.posterior.ub.HIP)

#Height of order restricted prior at deltaVT = 0
prior.ub <- dnorm(x = 0,mean = 0,sd = 0.5)*2          

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

write.xlsx(x = BF.tab,file = './Results/Hypothesis_testing/Bayes_factors_covar_gender.xlsx',row.names = F)

