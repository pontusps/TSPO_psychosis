
            #############################################################
# ##########    Graphical output: Forest plots and posteriors over Tau   ########## #
            #############################################################

                  # Pontus P. Sigray, KI, Stockholm, May 2017

#Load packages
library(tidyverse)
library(bayesplot)
library(xlsx)

#Source plot script
source("./R/Graphics/brms_forest/brms_forest_HB.R")

#Load brms fits for model 1 and model 3 
attach('./Results/Parameter_estimation/par_estimation_M1.RData')
attach('./Results/Parameter_estimation/par_estimation_M3.RData')
            

##############################################
# Graphics: Brms-forest-plots of random study slopes and grand delta
##############################################
#For brms forest plot see: https://mvuorre.github.io/post/2017/better-brms-forest-plots/

#Plot random slopes and grand delta using Model 3 estimates
plot.FC<-bayes_forest_HB(data = dfModelDat,model = fit3.FC,ROI = "FC",
         level = 0.95,show_data = T, dens_fill = "#86aec2" )

#Customize labels 
plot.FC<-plot.FC + 
scale_y_continuous(expression(paste('Frontal Cortex ',Delta,V[T])),limits = c(-2,3)) + 
theme(axis.title.x = element_text(size=14,face="bold")
,axis.text.y = element_text(size=11) 
)

##Save plot
ggsave(filename = './Results/Graphics/Bayes_forest_FC.png',plot = plot.FC,dpi = 300)

#Plot random slopes and grand delta using Model 3 estimates
plot.TC<-bayes_forest_HB(data = dfModelDat,model = fit3.TC,ROI = "TC", 
         level = 0.95,show_data = T,dens_fill = "#6497b1" )
#Customize labels 
plot.TC<-plot.TC + 
scale_y_continuous(expression(paste('Temporal Cortex ',Delta,V[T])),limits = c(-2,3)) + 
theme(axis.title.x = element_text(size=14,face="bold")
,axis.text.y=element_blank()
) 

##Save plot
ggsave(filename = './Results/Graphics/Bayes_forest_TC.png',plot = plot.TC,dpi = 300)

#Plot random slopes and grand delta using Model 3 estimates
plot.HIP<-bayes_forest_HB(data = dfModelDat,model = fit3.HIP,ROI = "HIP",
          level = 0.95,show_data = T,dens_fill = "#4c7d96" )
#Customize labels 
plot.HIP<-plot.HIP + 
scale_y_continuous(expression(paste('Hippocampus ',Delta,V[T])),limits = c(-2,3)) +
theme(axis.title.x = element_text(size=14,face="bold")
,axis.text.y=element_blank()
) 

##Save plot
ggsave(filename = './Results/Graphics/Bayes_forest_HIP.png',plot = plot.HIP,dpi = 300)

gridGrand<-gridExtra::grid.arrange(plot.FC, plot.TC, plot.HIP, ncol=3, widths=c(4.89,4,4))
ggsave(filename = './Results/Graphics/Bayes_forest_All_Rois.png',
plot = gridGrand,dpi = 300,width = 12,height = 5)


                