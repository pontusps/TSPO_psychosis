

            #####################################################################
# ##########    Graphical output: Posteriors over Random study slopes sd (tau)   ########## #
            #####################################################################

                        # Pontus P. Sigray, KI, Stockholm, May 2017

                          
library(tidyverse)
library(bayesplot)
library(xlsx)

#Load brms fits from model 3 
attach('./Results/Parameter_estimation/par_estimation_M3.RData')

######################
# Posteriors over tau
######################

#Colors taken from bayesplot package 
colScheme.blue<-bayesplot::color_scheme_set(scheme = "blue")
fill.col<-colScheme.blue$light
highlight.col<-colScheme.blue$light_highlight

samp.M3.FC<- rstan::extract(fit3.FC$fit)
samp.M3.TC<- rstan::extract(fit3.TC$fit)
samp.M3.HIP<- rstan::extract(fit3.HIP$fit)

tau.FC<-samp.M3.FC$sd_Study__HC_patpat
tau.TC<-samp.M3.TC$sd_Study__HC_patpat
tau.HIP<-samp.M3.HIP$sd_Study__HC_patpat

#Fit density kernels
dens.M3.tau.FC<-polspline::logspline(tau.FC,lbound = 0)
dens.M3.tau.TC<-polspline::logspline(tau.TC,lbound = 0)
dens.M3.tau.HIP<-polspline::logspline(tau.HIP,lbound = 0)

#Get mode of posteriors
#idx.tau.FC<-which.max(dens.M3.tau.FC$y)
#mode.M3.tau.FC<-dens.M3.tau.FC$x[idx.tau.FC]

#Save plot 
png(file = './Results/Graphics/Tau_posteriors.png',units = "cm",width = 7.5*1.2,height = 12.5*1.2, res=300 )

#Plot posterior
par(xaxs="i",yaxs="i") #Make axis connect
par(mfrow=c(3,1) ) #Number of facets in plot
par(mar=c(3.2,4,1,1) ) #Number of facets in plot

line1.thick<-6.5
line2.thick<-line1.thick-3

###Plot posterior for FC tau
hist(tau.FC, col=fill.col, freq = F, 
     ylim = c(0,3.5), xlim = c(0,1),breaks = 200, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Frontal Cortex ',tau)),line = 2.3,cex.lab=1.2)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1.2 )


#Add smoothed line for posterior
plot(dens.M3.tau.FC, lwd=line1.thick, add=T)
plot(dens.M3.tau.FC, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) dcauchy(x,0,.707),xlim=c(0,1),add=T,lwd=line1.thick)
plot(function(x) dcauchy(x,0,.707),xlim=c(0,1),add=T,lwd=line2.thick,col="white")

legend(.4,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(.4,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

box()

###Plot posterior for TC tau
hist(tau.TC, col=fill.col, freq = F, 
     ylim = c(0,3.5), xlim = c(0,1),breaks = 200, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Temporal Cortex ',tau)),line = 2.2,cex.lab=1.2)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1.2 )

#Add smoothed line for posterior
plot(dens.M3.tau.TC, lwd=line1.thick, add=T)
plot(dens.M3.tau.TC, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) dcauchy(x,0,.707),xlim=c(0,1),add=T,lwd=line1.thick)
plot(function(x) dcauchy(x,0,.707),xlim=c(0,1),add=T,lwd=line2.thick,col="white")

legend(.4,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(.4,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

box()

###Plot posterior for HIP tau
hist(tau.HIP, col=fill.col, freq = F, 
     ylim = c(0,4.5), xlim = c(0,1),breaks = 200, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Hippocampus ',tau)),line = 2.3,cex.lab=1.2)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1.2 )

#Add smoothed line for posterior
plot(dens.M3.tau.HIP, lwd=line1.thick, add=T)
plot(dens.M3.tau.HIP, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) dcauchy(x,0,.707),xlim=c(0,1),add=T,lwd=line1.thick)
plot(function(x) dcauchy(x,0,.707),xlim=c(0,1),add=T,lwd=line2.thick,col="white")

legend(.4,3.5, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(.4,3.5, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

box()

dev.off()
