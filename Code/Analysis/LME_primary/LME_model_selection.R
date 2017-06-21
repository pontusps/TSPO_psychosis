

            ################################################
# ##########               LME model selection              ########## #
            ################################################

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

# Uninformative Prior over the fixed effect: deltaVT. NO TRUNCATION for fit evaluation. 
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
  iter = 21000, 
  warmup = 1000
)

#######################
# Evaluate model fits
#######################

#M0: Model Intercept only
priors.M0<-c(prior.RE.sigma,
             prior.RE.Intc.gene,
             prior.RE.Intc.study)

args.M0.FC<-c(list(prior = priors.M0),base_args)
args.M0.FC$formula <- FC.VT.z ~ 1  + (1|Genotype) + (1|Study)
capture.output(fit0.FC <- do.call(brm, args.M0.FC))
#plot(fit0.FC)

args.M0.TC<-c(list(prior = priors.M0),base_args)
args.M0.TC$formula <- TC.VT.z ~ 1  + (1|Genotype) + (1|Study)
capture.output(fit0.TC <- do.call(brm, args.M0.TC))
#plot(fit0.TC)

args.M0.HIP<-c(list(prior = priors.M0),base_args)
args.M0.HIP$formula <- HIP.VT.z ~ 1  + (1|Genotype) + (1|Study)
capture.output(fit0.HIP <- do.call(brm, args.M0.HIP))
#plot(fit0.HIP)

fit0.FC.LOO<-brms::LOO(fit0.FC); fit0.FC.WAIC<-brms::WAIC(fit0.FC)

fit0.TC.LOO<-brms::LOO(fit0.TC); fit0.TC.WAIC<-brms::WAIC(fit0.TC)

fit0.HIP.LOO<-brms::LOO(fit0.HIP); fit0.HIP.WAIC<-brms::WAIC(fit0.HIP)

#M1: fixed effect deltaVT
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

fit1.FC.LOO<-brms::LOO(fit1.FC); fit1.FC.WAIC<-brms::WAIC(fit1.FC)

fit1.TC.LOO<-brms::LOO(fit1.TC); fit1.TC.WAIC<-brms::WAIC(fit1.TC)

fit1.HIP.LOO<-brms::LOO(fit1.HIP); fit1.HIP.WAIC<-brms::WAIC(fit1.HIP)

#M2: random slope for genotype
priors.M2<-c(prior.FE,
             prior.RE.sigma,
             prior.RE.Intc.gene,
             prior.RE.Intc.study,
             prior.RE.slope.gene)

args.M2.FC<-c(list(prior = priors.M2),base_args)
args.M2.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1 + HC_pat|Genotype) + (1|Study)
capture.output(fit2.FC <- do.call(brm, args.M2.FC))
#plot(fit2.FC)

args.M2.TC<-c(list(prior = priors.M2),base_args)
args.M2.TC$formula <- TC.VT.z ~ 1 + HC_pat  + (1 + HC_pat|Genotype) + (1|Study)
capture.output(fit2.TC <- do.call(brm, args.M2.TC))
#plot(fit2.TC)

args.M2.HIP<-c(list(prior = priors.M2),base_args)
args.M2.HIP$formula <- HIP.VT.z ~ 1 + HC_pat  + (1 + HC_pat|Genotype) + (1|Study)
capture.output(fit2.HIP <- do.call(brm, args.M2.HIP))
#plot(fit2.HIP)

fit2.FC.LOO<-brms::LOO(fit2.FC); fit2.FC.WAIC<-brms::WAIC(fit2.FC)

fit2.TC.LOO<-brms::LOO(fit2.TC); fit2.TC.WAIC<-brms::WAIC(fit2.TC)

fit2.HIP.LOO<-brms::LOO(fit2.HIP); fit2.HIP.WAIC<-brms::WAIC(fit2.HIP)


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

fit3.FC.LOO<-brms::LOO(fit3.FC); fit3.FC.WAIC<-brms::WAIC(fit3.FC)

fit3.TC.LOO<-brms::LOO(fit3.TC); fit3.TC.WAIC<-brms::WAIC(fit3.TC)

fit3.HIP.LOO<-brms::LOO(fit3.HIP); fit3.HIP.WAIC<-brms::WAIC(fit3.HIP)

#M4: random slope for genotype and study
priors.M4<-c(prior.FE,
             prior.RE.sigma,
             prior.RE.Intc.gene,
             prior.RE.Intc.study,
             prior.RE.slope.gene,
             prior.RE.slope.study)

args.M4.FC<-c(list(prior = priors.M4),base_args)
args.M4.FC$formula <- FC.VT.z ~ 1 + HC_pat + (1 + HC_pat|Genotype) + (1 + HC_pat|Study)
capture.output(fit4.FC <- do.call(brm, args.M4.FC))
#plot(fit4.FC)

args.M4.TC<-c(list(prior = priors.M4),base_args)
args.M4.TC$formula <- TC.VT.z ~ 1 + HC_pat  + (1 + HC_pat|Genotype) + (1 + HC_pat|Study)
capture.output(fit4.TC <- do.call(brm, args.M4.TC))
#plot(fit4.TC)


#Obsereved a small drift in one of the chains: Add thinning (thin=5) and ramp up samples (iter*5)
M4.HIP.args<-base_args
M4.HIP.args$iter<-base_args$iter*5
M4.HIP.args$thin<-5
M4.HIP.args$warmup<-base_args$warmup*5

args.M4.HIP<-c(list(prior = priors.M4),M4.HIP.args)
args.M4.HIP$formula <- HIP.VT.z ~ 1 + HC_pat  + (1 + HC_pat|Genotype) + (1 + HC_pat|Study)
capture.output(fit4.HIP <- do.call(brm, args.M4.HIP))
#plot(fit4.HIP)

fit4.FC.LOO<-brms::LOO(fit4.FC); fit4.FC.WAIC<-brms::WAIC(fit4.FC)

fit4.TC.LOO<-brms::LOO(fit4.TC); fit4.TC.WAIC<-brms::WAIC(fit4.TC)

fit4.HIP.LOO<-brms::LOO(fit4.HIP); fit4.HIP.WAIC<-brms::WAIC(fit4.HIP)

################################
# Output tables for model fit 
################################

Model<-c(0:4)

LOO_score.FC<-c(fit0.FC.LOO$looic,fit1.FC.LOO$looic,fit2.FC.LOO$looic,
                fit3.FC.LOO$looic,fit4.FC.LOO$looic)
WAIC_score.FC<-c(fit0.FC.WAIC$waic,fit1.FC.WAIC$waic,fit2.FC.WAIC$waic,
                 fit3.FC.WAIC$waic,fit4.FC.WAIC$waic)

LOO_score.TC<-c(fit0.TC.LOO$looic,fit1.TC.LOO$looic,fit2.TC.LOO$looic,
                fit3.TC.LOO$looic,fit4.TC.LOO$looic)
WAIC_score.TC<-c(fit0.TC.WAIC$waic,fit1.TC.WAIC$waic,fit2.TC.WAIC$waic,
                 fit3.TC.WAIC$waic,fit4.TC.WAIC$waic)

LOO_score.HIP<-c(fit0.HIP.LOO$looic,fit1.HIP.LOO$looic,fit2.HIP.LOO$looic,
                 fit3.HIP.LOO$looic,fit4.HIP.LOO$looic)
WAIC_score.HIP<-c(fit0.HIP.WAIC$waic,fit1.HIP.WAIC$waic,fit2.HIP.WAIC$waic,
                  fit3.HIP.WAIC$waic,fit4.HIP.WAIC$waic)


#Calculate difference to best fitting model and Akaike Weights
#FC
FC.WAIC.delta<-WAIC_score.FC-min(WAIC_score.FC)
FC.LOO.delta<-LOO_score.FC-min(LOO_score.FC)

FC.AkaikeWeights<- (exp(-.5*FC.LOO.delta))/(sum(exp(-.5*FC.LOO.delta)))
fit.tab.FC<-data.frame(Model = Model, LOO_score.FC = LOO_score.FC, WAIC_score.FC = WAIC_score.FC,
                        dLOO_score.FC = FC.LOO.delta, dWAIC_score.FC=FC.WAIC.delta, 
                        AkaikeWeights.FC = FC.AkaikeWeights)
knitr::kable(fit.tab.FC,digits = 2)


#TC
TC.WAIC.delta<-WAIC_score.TC-min(WAIC_score.TC)
TC.LOO.delta<-LOO_score.TC-min(LOO_score.TC)

TC.AkaikeWeights<- (exp(-.5*TC.LOO.delta))/(sum(exp(-.5*TC.LOO.delta)))
fit.tab.TC<-data.frame(Model = Model, LOO_score.TC = LOO_score.TC, WAIC_score.TC = WAIC_score.TC,
                        dLOO_score.TC = TC.LOO.delta, dWAIC_score.TC=TC.WAIC.delta, 
                        AkaikeWeights.TC = TC.AkaikeWeights)
knitr::kable(fit.tab.TC,digits = 2)

#HIP
HIP.WAIC.delta<-WAIC_score.HIP-min(WAIC_score.HIP)
HIP.LOO.delta<-LOO_score.HIP-min(LOO_score.HIP)
  
HIP.AkaikeWeights<- (exp(-.5*HIP.LOO.delta))/(sum(exp(-.5*HIP.LOO.delta)))
fit.tab.HIP<-data.frame(Model = Model, LOO_score.HIP = LOO_score.HIP, WAIC_score.HIP = WAIC_score.HIP,
                        dLOO_score.HIP = HIP.LOO.delta, dWAIC_score.HIP=HIP.WAIC.delta, 
                        AkaikeWeights.HIP = HIP.AkaikeWeights)
knitr::kable(fit.tab.HIP,digits = 2)

#Write fit tables
write.xlsx(x = fit.tab.FC,file = './Results/Model_fit/FC_model_fits.xlsx',row.names = F)
write.xlsx(x = fit.tab.TC,file = './Results/Model_fit/TC_model_fits.xlsx',row.names = F)
write.xlsx(x = fit.tab.HIP,file = './Results/Model_fit/HIP_model_fits.xlsx',row.names = F)

#Make omnibuss table
header.names<-c("Region",'Model','dLOO','dWAIC','Akaike_weights')

tmp.tab.FC<-fit.tab.FC[,c('Model','dLOO_score.FC','dWAIC_score.FC','AkaikeWeights.FC')]
tmp.tab.FC<-cbind(rep('FC',nrow(tmp.tab.FC)),tmp.tab.FC)
names(tmp.tab.FC)<-header.names

tmp.tab.TC<-fit.tab.TC[,c('Model','dLOO_score.TC','dWAIC_score.TC','AkaikeWeights.TC')]
tmp.tab.TC<-cbind(rep('TC',nrow(tmp.tab.TC)),tmp.tab.TC)
names(tmp.tab.TC)<-header.names

tmp.tab.HIP<-fit.tab.HIP[,c('Model','dLOO_score.HIP','dWAIC_score.HIP','AkaikeWeights.HIP')]
tmp.tab.HIP<-cbind(rep('HIP',nrow(tmp.tab.HIP)),tmp.tab.HIP)
names(tmp.tab.HIP)<-header.names

#Make fit-table for all ROIs
allRois<-rbind(tmp.tab.FC,tmp.tab.TC,tmp.tab.HIP)

#Write omnibus fit table
write.xlsx(x = allRois,file = './Results/Model_fit/All_rois_fit.xlsx',row.names = F)
knitr::kable(allRois,digits = 2)



