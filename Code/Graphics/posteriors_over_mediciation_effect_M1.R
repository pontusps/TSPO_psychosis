    
                ######################################################################
    # ##########    Graphical output: Posteriors over medication effect from model 1  ########## #
                ######################################################################
    
                           # Pontus P. Sigray, KI, Stockholm, May 2017


library(tidyverse)
library(bayesplot)
library(xlsx)

#Load brms fits from model 1 with medication effect as additional predictor 
attach('./Results/Medication_status/Med_status_models.RData')

######################
# Posteriors over medication effects
######################

#Colors taken from bayesplot package 
colScheme.blue<-bayesplot::color_scheme_set(scheme = "blue")
fill.col<-colScheme.blue$light
highlight.col<-colScheme.blue$light_highlight

samp.FC<- rstan::extract(fit1.FC$fit)
samp.TC<- rstan::extract(fit1.TC$fit)
samp.HIP<- rstan::extract(fit1.HIP$fit)

stat.med.FC<-samp.FC$b_Med.statusmedicated
stat.med.TC<-samp.TC$b_Med.statusmedicated
stat.med.HIP<-samp.HIP$b_Med.statusmedicated

deltaVT.FC<-samp.FC$b_HC_patpat
deltaVT.TC<-samp.TC$b_HC_patpat
deltaVT.HIP<-samp.HIP$b_HC_patpat

#Make posterior samples into a matrix for the bayesplot-package 
mat.samps.FC<-as.matrix(as.data.frame(cbind(samp.FC$b_HC_patpat,samp.FC$b_Med.statusmedicated)))
colnames(mat.samps.FC)<-c("Patient-control effect","Medication effect")

mat.samps.TC<-as.matrix(as.data.frame(cbind(samp.TC$b_HC_patpat,samp.TC$b_Med.statusmedicated)))
colnames(mat.samps.TC)<-c("Patient-control effect","Medication effect")

mat.samps.HIP<-as.matrix(as.data.frame(cbind(samp.HIP$b_HC_patpat,samp.HIP$b_Med.statusmedicated)))
colnames(mat.samps.HIP)<-c("Patient-control effect","Medication effect")

##################################################################
# Plot of deltaVT and the additional effect of medication status
##################################################################

#Basedrop the bayesplot of FC, as a tribute to prof. Peng
plot.med.stat.FC<-bayesplot::mcmc_areas(x = mat.samps.FC, 
                                        pars = c("Patient-control effect","Medication effect"),
                                        prob = 0.95,point_est = "mean") +
                              xlab(expression(paste("Frontal Cortex ",Delta,V[T]))) + 
                              theme(panel.border = element_rect(colour = "black", fill=NA, size=1))



#Basedrop the bayesplot of TC, as a tribute to prof. Peng
plot.med.stat.TC<-bayesplot::mcmc_areas(x = mat.samps.TC, 
                                        pars = c("Patient-control effect","Medication effect"),
                                        prob = 0.95,point_est = "mean")  +
                              xlab(expression(paste("Temporal Cortex ",Delta,V[T])))  + 
                              theme(panel.border = element_rect(colour = "black", fill=NA, size=1))



#Basedrop the bayesplot of HIP, as a tribute to prof. Peng
plot.med.stat.HIP<-bayesplot::mcmc_areas(x = mat.samps.HIP, 
                                         pars = c("Patient-control effect","Medication effect"),
                                         prob = 0.95,point_est = "mean") +
                              xlab(expression(paste("Hippocampus ",Delta,V[T])))  + 
                              theme(panel.border = element_rect(colour = "black", fill=NA, size=1))



gridOut <- gridExtra::grid.arrange(plot.med.stat.FC,plot.med.stat.TC,plot.med.stat.HIP)
ggsave(plot = gridOut,filename = paste0("./Results/Graphics/Effect_medication_deltaVT.png"),
       width = 20*0.75,height = 30*0.75,units = "cm",device = "png",dpi = 300)




