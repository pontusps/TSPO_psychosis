

            #########################################################################
# ##########    Graphical output: Posteriors over deltaVT from Hypothesis testing    ########## #
            ########################################################################

                          # Pontus P. Sigray, KI, Stockholm, May 2017


library(tidyverse)
library(bayesplot)
library(xlsx)

################################
# Posteriors over deltaVT - H1 
###############################

#Load brms fits for Hypothesis testing 
attach('./Results/Hypothesis_testing/HypothesisTesting_H1_models.RData')

#Colors taken from bayesplot package 
colScheme.blue<-bayesplot::color_scheme_set(scheme = "blue")
fill.col<-colScheme.blue$light
highlight.col<-colScheme.blue$light_highlight

samp.lb.FC<- rstan::extract(fit.lb.FC$fit)
samp.lb.TC<- rstan::extract(fit.lb.TC$fit)
samp.lb.HIP<- rstan::extract(fit.lb.HIP$fit)

deltaVT.lb.FC<-samp.lb.FC$b_HC_patpat
deltaVT.lb.TC<-samp.lb.TC$b_HC_patpat
deltaVT.lb.HIP<-samp.lb.HIP$b_HC_patpat

#Fit density kernels
dens.lb.deltaVT.FC<-polspline::logspline(deltaVT.lb.FC,lbound = 0)
dens.lb.deltaVT.TC<-polspline::logspline(deltaVT.lb.TC,lbound = 0)
dens.lb.deltaVT.HIP<-polspline::logspline(deltaVT.lb.HIP,lbound = 0)

#Get mode of posteriors
#idx.deltaVT.FC<-which.max(dens.lb.deltaVT.FC$y)
#mode.lb.deltaVT.FC<-dens.lb.deltaVT.FC$x[idx.deltaVT.FC]

##############
# Save plot 
##############

#pdf(file = './Results/Graphics/deltaVT_posteriors_H1.pdf',width = (7-4),height = (7-2))
png(file = './Results/Graphics/deltaVT_posteriors_H1.png',units = "cm",width = 7.5*1.2,height = 12.5*1.2, res=300 )

#Plot posterior
par(xaxs="i",yaxs="i") #Make axis connect
par(mfrow=c(3,1) ) #Number of facets in plot
par(mar=c(3.2,4,1,1) ) #Number of facets in plot

line1.thick<-6
line2.thick<-line1.thick-3

x = 0 # x axis position for Savage-Dickey markers
y.prior = dnorm(0,0,.5)*2 # y axis position for Savage-Dickey prior marker

###Plot posterior for FC deltaVT
hist(deltaVT.lb.FC, col=fill.col, freq = F, 
     ylim = c(0,25), xlim = c(0,0.2),breaks = 100, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Frontal Cortex ',Delta,'V'[T])),line = 2.3,cex.lab=1)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1 )


#Add smoothed line for posterior
plot(dens.lb.deltaVT.FC, lwd=line1.thick, add=T)
plot(dens.lb.deltaVT.FC, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) truncnorm::dtruncnorm(x,a = 0,mean = 0,sd = 0.5),xlim = c(0,0.2),add=T,lwd=line1.thick)
plot(function(x) truncnorm::dtruncnorm(x,a = 0,mean = 0,sd = 0.5),xlim = c(0,0.2),add=T,lwd=line2.thick,col="white")

y.post = polspline::dlogspline(0,fit = dens.lb.deltaVT.FC)
points(x, y.prior, pch=21, bg='white', col='black',cex=2)
points(x, y.post, pch=21, bg=highlight.col, col='black',cex=2)

legend(.1,20, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(.1,20, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width


###Plot posterior for TC deltaVT
hist(deltaVT.lb.TC, col=fill.col, freq = F, 
     ylim = c(0,25), xlim = c(0,0.2),breaks = 100, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Temporal Cortex ',Delta,'V'[T])),line = 2.2,cex.lab=1)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1 )

#Add smoothed line for posterior
plot(dens.lb.deltaVT.TC, lwd=line1.thick, add=T)
plot(dens.lb.deltaVT.TC, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) truncnorm::dtruncnorm(x,a = 0,mean = 0,sd = 0.5),xlim = c(0,0.2),add=T,lwd=line1.thick)
plot(function(x) truncnorm::dtruncnorm(x,a = 0,mean = 0,sd = 0.5),xlim = c(0,0.2),add=T,lwd=line2.thick,col="white")

y.post = polspline::dlogspline(0,fit = dens.lb.deltaVT.TC)
points(x, y.prior, pch=21, bg='white', col='black',cex=2)
points(x, y.post, pch=21, bg=highlight.col, col='black',cex=2)

legend(.1,20, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(.1,20, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

###Plot posterior for HIP deltaVT
hist(deltaVT.lb.HIP, col=fill.col, freq = F, 
     ylim = c(0,30), xlim = c(0,0.2),breaks = 100, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Hippocampus ',Delta,'V'[T])),line = 2.3,cex.lab=1)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1 )

#Add smoothed line for posterior
plot(dens.lb.deltaVT.HIP, lwd=line1.thick, add=T)
plot(dens.lb.deltaVT.HIP, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) truncnorm::dtruncnorm(x,a = 0,mean = 0,sd = 0.5),xlim = c(0,0.2),add=T,lwd=line1.thick)
plot(function(x) truncnorm::dtruncnorm(x,a = 0,mean = 0,sd = 0.5),xlim = c(0,0.2),add=T,lwd=line2.thick,col="white")

y.post = polspline::dlogspline(0,fit = dens.lb.deltaVT.HIP)
points(x, y.prior, pch=21, bg='white', col='black',cex=2)
points(x, y.post, pch=21, bg=highlight.col, col='black',cex=2)

legend(.1,25, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(.1,25, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width


dev.off()


################################
# Posteriors over deltaVT - H2
###############################

#Load brms fits for Hypothesis testing 
attach('./Results/Hypothesis_testing/HypothesisTesting_H2_models.RData')

#Colors taken from bayesplot package 
colScheme.blue<-bayesplot::color_scheme_set(scheme = "blue")
fill.col<-colScheme.blue$light
highlight.col<-colScheme.blue$light_highlight

samp.ub.FC<- rstan::extract(fit.ub.FC$fit)
samp.ub.TC<- rstan::extract(fit.ub.TC$fit)
samp.ub.HIP<- rstan::extract(fit.ub.HIP$fit)

deltaVT.ub.FC<-samp.ub.FC$b_HC_patpat
deltaVT.ub.TC<-samp.ub.TC$b_HC_patpat
deltaVT.ub.HIP<-samp.ub.HIP$b_HC_patpat

#Fit density kernels
dens.ub.deltaVT.FC<-polspline::logspline(deltaVT.ub.FC,ubound = 0)
dens.ub.deltaVT.TC<-polspline::logspline(deltaVT.ub.TC,ubound = 0)
dens.ub.deltaVT.HIP<-polspline::logspline(deltaVT.ub.HIP,ubound = 0)

#Get mode of posteriors
#idx.deltaVT.FC<-which.max(dens.ub.deltaVT.FC$y)
#mode.ub.deltaVT.FC<-dens.ub.deltaVT.FC$x[idx.deltaVT.FC]

#Save plot 
#pdf(file = './Results/Graphics/deltaVT_posteriors_H2.pdf',width = 7.5*1.2,height = 12.5*1.2)
png(file = './Results/Graphics/deltaVT_posteriors_H2.png',units = "cm",width = 7.5*1.2,height = 12.5*1.2, res=200 )

#Plot posterior
par(xaxs="i",yaxs="i") #Make axis connect
par(mfrow=c(3,1) ) #Number of facets in plot
par(mar=c(3.2,4,1,1) ) #Number of facets in plot

line1.thick<-6
line2.thick<-line1.thick-3

###Plot posterior for FC deltaVT
hist(deltaVT.ub.FC, col=fill.col, freq = F, 
     ylim = c(0,4), xlim = c(-1,0),breaks = 75, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Frontal Cortex ',Delta,'V'[T])),line = 2.3,cex.lab=1)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1 )


#Add smoothed line for posterior
plot(dens.ub.deltaVT.FC, lwd=line1.thick, add=T)
plot(dens.ub.deltaVT.FC, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) truncnorm::dtruncnorm(x,b = 0,mean = 0,sd = 0.5),xlim = c(-1,0),add=T,lwd=line1.thick)
plot(function(x) truncnorm::dtruncnorm(x,b = 0,mean = 0,sd = 0.5),xlim = c(-1,0),add=T,lwd=line2.thick,col="white")

y.post = polspline::dlogspline(0,fit = dens.ub.deltaVT.FC)
points(x, y.prior, pch=21, bg='white', col='black',cex=2)
points(x, y.post, pch=21, bg=highlight.col, col='black',cex=2)

legend(-0.95,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(-0.95,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

###Plot posterior for TC deltaVT
hist(deltaVT.ub.TC, col=fill.col, freq = F, 
     ylim = c(0,3.5), xlim = c(-1,0),breaks = 75, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Temporal Cortex ',Delta,'V'[T])),line = 2.2,cex.lab=1)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1 )

#Add smoothed line for posterior
plot(dens.ub.deltaVT.TC, lwd=line1.thick, add=T)
plot(dens.ub.deltaVT.TC, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) truncnorm::dtruncnorm(x,b = 0,mean = 0,sd = 0.5),xlim = c(-1,0),add=T,lwd=line1.thick)
plot(function(x) truncnorm::dtruncnorm(x,b = 0,mean = 0,sd = 0.5),xlim = c(-1,0),add=T,lwd=line2.thick,col="white")

y.post = polspline::dlogspline(0,fit = dens.ub.deltaVT.TC)
points(x, y.prior, pch=21, bg='white', col='black',cex=2)
points(x, y.post, pch=21, bg=highlight.col, col='black',cex=2)

legend(-0.95,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(-0.95,3, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

###Plot posterior for HIP deltaVT
hist(deltaVT.ub.HIP, col=fill.col, freq = F, 
     ylim = c(0,4), xlim = c(-1,0),breaks = 75, 
     xlab="" ,
     ylab="",
     main="",
     cex.lab = 1,las=1)
title(xlab = expression(paste('Hippocampus ',Delta,'V'[T])),line = 2.3,cex.lab=1)
title(ylab = expression(paste('Probability')),line = 2.3,cex.lab=1 )

#Add smoothed line for posterior
plot(dens.ub.deltaVT.HIP, lwd=line1.thick, add=T)
plot(dens.ub.deltaVT.HIP, lwd=line2.thick, col=highlight.col, add=T)

plot(function(x) truncnorm::dtruncnorm(x,b = 0,mean = 0,sd = 0.5),xlim = c(-1,0),add=T,lwd=line1.thick)
plot(function(x) truncnorm::dtruncnorm(x,b = 0,mean = 0,sd = 0.5),xlim = c(-1,0),add=T,lwd=line2.thick,col="white")

y.post = polspline::dlogspline(0,fit = dens.ub.deltaVT.HIP)
points(x, y.prior, pch=21, bg='white', col='black',cex=2)
points(x, y.post, pch=21, bg=highlight.col, col='black',cex=2)


legend(-0.95,3.75, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line1.thick,line1.thick),col=c("black","black"))# gives the legend lines the correct color and width

legend(-0.95,3.75, bty = "n",# places a legend at the appropriate place 
       c("Posterior","Prior"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(line2.thick,line2.thick),col=c(highlight.col,"white"))# gives the legend lines the correct color and width

dev.off()
