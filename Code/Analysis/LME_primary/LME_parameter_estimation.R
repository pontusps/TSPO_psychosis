

            ############################################
# ##########          LME parameter estimation          ########## #
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

##########
# Priors 
##########

# Uninformative Prior over the fixed effect: deltaVT. NO TRUNCATION for parameter estimation 
prior.FE <- set_prior("normal(0,10)", class = "b") 

# Prior for the "grand" sigma term in the likelihood (brms default now)
prior.RE.sigma <- set_prior("student_t(3,0,10)", class = "sigma") 

#  Prior on freely varying intercept for genotype group  
prior.RE.Intc.gene <-set_prior("cauchy(0,0.707)", class = "sd",group = "Genotype", coef = "Intercept") 

# Prior on freely varying intercept for study group
prior.RE.Intc.study <-set_prior("cauchy(0,0.707)", class = "sd",group = "Study", coef = "Intercept") 

# Prior on freely varying slope for genotype
prior.RE.slope.gene <- set_prior("cauchy(0,0.707)", class = "sd",group = "Genotype")

# Prior on freely varying slope for study
prior.RE.slope.study <- set_prior("cauchy(0,0.707)", class = "sd",group = "Study") 

######################################
# Base arguments common for all models
######################################

base_args <- list(
  data = dfModelDat,
  control = list(adapt_delta = 0.99),
  iter = 21000*4, 
  warmup = 1000*4,
  thin=4
)

############################
# Parameter estimation M1 - the pat-HC difference is assumed to be equally large between studies, and between genotype groups
############################

priors.M1<-c(prior.FE,
             prior.RE.sigma,
             prior.RE.Intc.gene,
             prior.RE.Intc.study)

args.M1.FC<-c(list(prior = priors.M1),base_args)
args.M1.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1|Study)
capture.output(fit1.FC <- do.call(brm, args.M1.FC))
#plot(fit1.FC)

args.M1.TC<-c(list(prior = priors.M1),base_args)
args.M1.TC$formula <- TC.VT.z ~ 1 + HC_pat  + (1|Genotype) + (1|Study)
capture.output(fit1.TC <- do.call(brm, args.M1.TC))
#plot(fit1.TC)

args.M1.HIP<-c(list(prior = priors.M1),base_args)
args.M1.HIP$formula <- HIP.VT.z ~ 1 + HC_pat  + (1|Genotype) + (1|Study)
capture.output(fit1.HIP <- do.call(brm, args.M1.HIP))
#plot(fit1.HIP)


############################################################
# Parameter estimation M3: get "Study-Tau" (as required by PRISMA Guidelines)
############################################################

#M3: random slope for study
priors.M3<-c(prior.FE,
             prior.RE.sigma,
             prior.RE.Intc.gene,
             prior.RE.Intc.study,
             prior.RE.slope.study)

args.M3.FC<-c(list(prior = priors.M3),base_args)
args.M3.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1|Genotype) + (1 + HC_pat|Study)
capture.output(fit3.FC <- do.call(brm, args.M3.FC))
#plot(fit3.FC)


args.M3.TC<-c(list(prior = priors.M3),base_args)
args.M3.TC$formula <- TC.VT.z ~ 1 + HC_pat  + (1|Genotype) + (1 + HC_pat|Study)
capture.output(fit3.TC <- do.call(brm, args.M3.TC))
#plot(fit3.TC)

args.M3.HIP<-c(list(prior = priors.M3),base_args)
args.M3.HIP$formula <- HIP.VT.z ~ 1 + HC_pat  + (1|Genotype) + (1 + HC_pat|Study)
capture.output(fit3.HIP <- do.call(brm, args.M3.HIP))
#plot(fit3.HIP)

##################
# Save models 
##################

save(fit3.FC,fit3.TC,fit3.HIP,file = "./Results/Parameter_estimation/par_estimation_M3.RData")
save(fit1.FC,fit1.TC,fit1.HIP,file = "./Results/Parameter_estimation/par_estimation_M1.RData")

#########################################
# Extract MAP and HPDI for grandDelta
#########################################

samp.M3.FC<-rstan::extract(fit3.FC$fit)
samp.M3.TC<-rstan::extract(fit3.TC$fit)
samp.M3.HIP<-rstan::extract(fit3.HIP$fit) 

#Extract mean of posteriors 
mean.M3.delta.FC <- mean(as.numeric(samp.M3.FC$b_HC_patpat))
mean.M3.delta.TC <- mean(as.numeric(samp.M3.TC$b_HC_patpat))
mean.M3.delta.HIP <- mean(as.numeric(samp.M3.HIP$b_HC_patpat))

mean.M3.tau.FC <- mean(as.numeric(samp.M3.FC$sd_Study__HC_patpat))
mean.M3.tau.TC <- mean(as.numeric(samp.M3.TC$sd_Study__HC_patpat))
mean.M3.tau.HIP <- mean(as.numeric(samp.M3.HIP$sd_Study__HC_patpat))

#Extract HPDI for grandDelta and tau (SD of random slopes of study)
HPDI.M3.deltaVT.FC <- rethinking::HPDI(as.numeric(samp.M3.FC$b_HC_patpat),prob = .95)
HPDI.M3.deltaVT.TC <- rethinking::HPDI(as.numeric(samp.M3.TC$b_HC_patpat),prob = .95)
HPDI.M3.deltaVT.HIP <- rethinking::HPDI(as.numeric(samp.M3.HIP$b_HC_patpat),prob = .95)

HPDI.M3.tau.FC <- rethinking::HPDI(as.numeric(samp.M3.FC$sd_Study__HC_patpat),prob = .95)
HPDI.M3.tau.TC <- rethinking::HPDI(as.numeric(samp.M3.TC$sd_Study__HC_patpat),prob = .95)
HPDI.M3.tau.HIP <- rethinking::HPDI(as.numeric(samp.M3.HIP$sd_Study__HC_patpat),prob = .95)

#Fit density kernels to find MAP of grand delta VT and tau

dens.M3.deltaVT.FC<-density(as.numeric(samp.M3.FC$b_HC_patpat))
dens.M3.deltaVT.TC<-density(as.numeric(samp.M3.TC$b_HC_patpat))
dens.M3.deltaVT.HIP<-density(as.numeric(samp.M3.HIP$b_HC_patpat))

dens.M3.tau.FC<-density(as.numeric(samp.M3.FC$sd_Study__HC_patpat))
dens.M3.tau.TC<-density(as.numeric(samp.M3.TC$sd_Study__HC_patpat))
dens.M3.tau.HIP<-density(as.numeric(samp.M3.HIP$sd_Study__HC_patpat))

#Get mode of posteriors
idx.deltaVT.FC<-which.max(dens.M3.deltaVT.FC$y)
mode.M1.deltaVT.FC<-dens.M3.deltaVT.FC$x[idx.deltaVT.FC]
idx.deltaVT.TC<-which.max(dens.M3.deltaVT.TC$y)
mode.M1.deltaVT.TC<-dens.M3.deltaVT.TC$x[idx.deltaVT.TC]
idx.deltaVT.HIP<-which.max(dens.M3.deltaVT.HIP$y)
mode.M1.deltaVT.HIP<-dens.M3.deltaVT.HIP$x[idx.deltaVT.HIP]

idx.tau.FC<-which.max(dens.M3.tau.FC$y)
mode.M1.tau.FC<-dens.M3.tau.FC$x[idx.tau.FC]
idx.tau.TC<-which.max(dens.M3.tau.TC$y)
mode.M1.tau.TC<-dens.M3.tau.TC$x[idx.tau.TC]
idx.tau.HIP<-which.max(dens.M3.tau.HIP$y)
mode.M1.tau.HIP<-dens.M3.tau.HIP$x[idx.tau.HIP]

#MAP and HPDI
sum.stat.M3.delta.FC <- c(mode.M1.deltaVT.FC,HPDI.M3.deltaVT.FC)
sum.stat.M3.delta.TC <- c(mode.M1.deltaVT.TC,HPDI.M3.deltaVT.TC)
sum.stat.M3.delta.HIP <- c(mode.M1.deltaVT.HIP,HPDI.M3.deltaVT.HIP)

sum.stat.M3.tau.FC <- c(mode.M1.tau.FC,HPDI.M3.tau.FC)
sum.stat.M3.tau.TC <- c(mode.M1.tau.TC,HPDI.M3.tau.TC)
sum.stat.M3.tau.HIP <- c(mode.M1.tau.HIP,HPDI.M3.tau.HIP)

#Sum stats into table
sum.stat.M3.delta.FC <- as.data.frame(cbind(mean.M3.delta.FC,t(as.data.frame(sum.stat.M3.delta.FC))))
colnames(sum.stat.M3.delta.FC)<-c('Mean','MAP','|0.95','0.95|')
sum.stat.M3.delta.TC <- as.data.frame(cbind(mean.M3.delta.TC,t(as.data.frame(sum.stat.M3.delta.TC))))
colnames(sum.stat.M3.delta.TC)<-c('Mean','MAP','|0.95','0.95|')
sum.stat.M3.delta.HIP <- as.data.frame(cbind(mean.M3.delta.HIP,t(as.data.frame(sum.stat.M3.delta.HIP))))
colnames(sum.stat.M3.delta.HIP)<-c('Mean','MAP','|0.95','0.95|')

sum.stat.M3.tau.FC <- as.data.frame(cbind(mean.M3.tau.FC,t(as.data.frame(sum.stat.M3.tau.FC))))
colnames(sum.stat.M3.tau.FC)<-c('Mean','MAP','|0.95','0.95|')
sum.stat.M3.tau.TC <- as.data.frame(cbind(mean.M3.tau.TC,t(as.data.frame(sum.stat.M3.tau.TC))))
colnames(sum.stat.M3.tau.TC)<-c('Mean','MAP','|0.95','0.95|')
sum.stat.M3.tau.HIP <- as.data.frame(cbind(mean.M3.tau.HIP,t(as.data.frame(sum.stat.M3.tau.HIP))))
colnames(sum.stat.M3.tau.HIP)<-c('Mean','MAP','|0.95','0.95|')

ROI_Model <- c('FC.M3','TC.M3','HIP.M3')

sum.stat.Table.M3.delta<-as.data.frame(cbind(ROI_Model,rbind(sum.stat.M3.delta.FC,sum.stat.M3.delta.TC,sum.stat.M3.delta.HIP)))
sum.stat.Table.M3.tau<-as.data.frame(cbind(ROI_Model,rbind(sum.stat.M3.tau.FC,sum.stat.M3.tau.TC,sum.stat.M3.tau.HIP)))

#Write table
write.xlsx(x = sum.stat.Table.M3.delta,file = './Results/Parameter_estimation/MAP_HPDI_M3_grandDelta.xlsx')



