
                          ##############################################
              # ##########  Graphical output: Beeswarm plot of raw data ########## #
                          ##############################################
              
                          # Pontus P. Sigray, KI, Stockholm, March 2017
#Load packages
library(tidyverse)
library(xlsx)
library(gridExtra)
library(grid)
library(ggbeeswarm)
library(RColorBrewer)

#Set color brewer
myColors<-c("#6497b1","#b17e64")

#Read in data
dfMaster<-read.xlsx("./DerivedData/All_studies_VT_no_overlap.xlsx",sheetIndex = 1)

#Remove Kenk hippocampus model fit failure  
dfMaster$HIP.VT[(dfMaster$HIP.VT>70)]<-NA

#Remove Bloomfield MAB pats
dfMaster<-dfMaster[!(dfMaster$Study=="Bloomfield" & dfMaster$ID=="MAB.pat"),]

#Add mean VT per study and per cell (ID) to each row in dataframe  
dfMaster<-dfMaster %>%
  group_by(Study,ID) %>%
  mutate(meanHIP=mean(HIP.VT,na.rm = T),
         meanTC=mean(TC.VT),
         meanDLPFC=mean(DLPFC.VT))

############
# Plot HIP #
############
dfMaster.hip<-dfMaster

#Perform "do()" on grouped data
plotCode.HIP<-dfMaster.hip %>%
  group_by(Study) %>%
  do(plot.code = ggplot(data=.) + aes(x = as.factor(ID), y =  HIP.VT ) 
     + geom_quasirandom(width=0.3,size=2.2,aes(shape=as.factor(HC_pat),colour = as.factor(Genotype))) 
     + geom_errorbar(width= 0.3,size=1.1,aes(y=meanHIP,ymin=meanHIP,ymax=meanHIP))
     + theme_bw()
     + scale_colour_manual(name='',labels=c('HABs  ','MABs'), values = myColors)  #fix legend
     + scale_shape(name='',labels=c('HC    ','SZ') )  #fix legend
     + theme(axis.title.x =element_blank(),
             axis.title.y = element_blank(),
             axis.ticks = element_blank(), 
             axis.text.x = element_blank(),
             axis.text.y = element_text(size=12),
             legend.position='none',
             legend.key = element_blank()
             #,plot.margin = unit(c(0.5,0.2,0.2,0.1),'cm') 
             )
     + expand_limits(y=0)
  ) 

############
# Plot TC #
############
dfMaster.TC<-dfMaster

#Perform "do()" on grouped data
plotCode.TC<-dfMaster.TC %>%
  group_by(Study) %>%
  do(plot.code = ggplot(data=.) + aes(x = as.factor(ID), y =  TC.VT ) 
     + geom_quasirandom(width=0.3,size=2.2,aes(shape=as.factor(HC_pat),colour = as.factor(Genotype))) 
     + geom_errorbar(width= 0.3,size=1.1,aes(y=meanTC,ymin=meanTC,ymax=meanTC))
     + theme_bw()
     + scale_colour_manual(name='',labels=c('HABs  ','MABs'), values = myColors)  #fix legend
     + scale_shape(name='',labels=c('HC    ','SZ') )  #fix legend
     + theme(axis.title.x =element_blank(),
             axis.title.y = element_blank(),
             axis.ticks = element_blank(), 
             axis.text.x = element_blank(),
             axis.text.y = element_text(size=12),
             legend.position='none',
             legend.key = element_blank()
             #,plot.margin = unit(c(0.5,0.2,0.2,0.1),'cm') 
             )
     + expand_limits(y=0)
  ) 

############
# Plot DLPFC #
############
dfMaster.DLPFC<-dfMaster

#Perform "do()" on grouped data
plotCode.DLPFC<-dfMaster.DLPFC %>%
  group_by(Study) %>%
  do(plot.code = ggplot(data=.) + aes(x = as.factor(ID), y =  DLPFC.VT ) 
     + geom_quasirandom(width=0.3,size=2.2,aes(shape=as.factor(HC_pat),colour = as.factor(Genotype))) 
     + geom_errorbar(width= 0.3,size=1.1,aes(y=meanDLPFC,ymin=meanDLPFC,ymax=meanDLPFC))
     + theme_bw()
     + scale_colour_manual(name='',labels=c('HABs  ','MABs'), values = myColors)  #fix legend
     + scale_shape(name='',labels=c('HC    ','SZ') )  #fix legend
     + theme(axis.title.x =element_blank(),
             axis.title.y = element_blank(),
             axis.ticks = element_blank(), 
             axis.text.x = element_blank(),
             axis.text.y = element_text(size=12),
             legend.position='none',
             legend.key = element_blank()
             #,plot.margin = unit(c(0.5,0.2,0.2,0.1),'cm') 
             )
     + expand_limits(y=0)
  ) 

########################
# Arrange into one plot
########################

# Rename to make grid table easier
plotCode.HIP$ROI<-'C.HIP'
plotCode.TC$ROI<-'B.TC'
plotCode.DLPFC$ROI<-'A.DLPFC'

#Merge all plot code into a tibble 
allPlotCode<-rbind(plotCode.DLPFC,plotCode.TC,plotCode.HIP)

#Rearrange order of Study
allPlotCode$Study <- factor(allPlotCode$Study,levels=c("Kenk","Bloomfield","Coughlin","Hafizi","Collste"))

allPlotCode<-allPlotCode %>%
  group_by(ROI) %>% 
  arrange(Study) %>%
  arrange(ROI)

### Set the widths of each plot-area to be the same (aligning the y-axes): 
grobs <- list()
widths <- list()

# Collect the widths for each grob of each plot
for (i in 1:length(allPlotCode$plot.code)){
  grobs[[i]] <- ggplotGrob(allPlotCode$plot.code[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

#Asign the max width to each grob
for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

#Arrange all grobs for all ROI into a 3x5 grid (ROI x Study). Add labels, legends and such in InkScape. 
#Layout for grid plot
layout<-rbind(c(1,2,3,4,5), 
              c(6,7,8,9,10),
              c(11,12,13,14,15))

#Save plot
grid.each.study<-gridExtra::grid.arrange(grobs=grobs, layout_matrix=layout)
ggsave(filename =paste0("./Results/Graphics/VT_scatter_all_studies.svg"),plot = grid.each.study,device = svg,width =25, height = 30, units = "cm")

##########################################
# Pooled plot of z-transformed VT values #
##########################################

#Z-transform values within studies, within genotype (set all HAB and MAB groups on "equal scaling").

dfMaster.z<-dfMaster %>%
  group_by(Study,Genotype) %>%
  mutate(HIP.z_wS=as.numeric(scale(HIP.VT)),
         TC.z_wS=as.numeric(scale(TC.VT)),
         DLPFC.z_wS=as.numeric(scale(DLPFC.VT)))

#Calculate Grand mean across all studies for each ID-group
dfMaster.z.plot <- dfMaster.z %>%
  group_by(ID) %>%
  mutate(meansHIP.z = mean(HIP.z_wS,na.rm = T),
         meansTC.z = mean(TC.z_wS,na.rm = T),
         meansDLPFC.z = mean(DLPFC.z_wS,na.rm = T)) 

dfMaster.z.plot<-as.data.frame(dfMaster.z.plot)

#Add column with color specs for Genotype groups
dfMaster.z.plot$MyCol[dfMaster.z.plot$Genotype=="HAB"]<-myColors[1]
dfMaster.z.plot$MyCol[dfMaster.z.plot$Genotype=="MAB"]<-myColors[2]

###########################
# Pooled plot of z-values #
###########################

#HIP
plot.code.HIP.All<- dfMaster.z.plot %>% 
  na.omit() %>%
  group_by(Genotype) %>%
  do( plot.code.All = ggplot(data = ., aes(x = as.factor(ID),y = HIP.z_wS)) +
        geom_quasirandom(color=unique(.$MyCol),width=0.2,size=2.2,aes(shape=as.factor(HC_pat) ) ) + 
        geom_errorbar(width= 0.3,size=1.1,
                      aes(y=meansHIP.z,ymin=meansHIP.z,ymax=meansHIP.z)) + 
        theme_bw() + 
        theme(legend.position='none',legend.key = element_blank()) + 
        scale_shape(name='',labels=c('HC    ','SZ')) +   
        theme(axis.title.y = element_blank())  + 
        theme(axis.title.x =element_blank()) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size=12) ) + 
        ylim(-2,3.5)
  )

#TC
plot.code.TC.All<- dfMaster.z.plot %>% 
  na.omit() %>%
  group_by(Genotype) %>%
  do( plot.code.All = ggplot(data = ., aes(x = as.factor(ID),y = TC.z_wS)) +
        geom_quasirandom(color=unique(.$MyCol),width=0.2,size=2.2,aes(shape=as.factor(HC_pat) ) ) + 
        geom_errorbar(width= 0.3,size=1.1,
                      aes(y=meansTC.z,ymin=meansTC.z,ymax=meansTC.z)) + 
        theme_bw() + 
        theme(legend.position='none',legend.key = element_blank()) + 
        scale_shape(name='',labels=c('HC    ','SZ')) +   
        theme(axis.title.y = element_blank())  + 
        theme(axis.title.x =element_blank()) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size=12) ) + 
        ylim(-2,3.5)
  )

#DLPFC
plot.code.DLPFC.All<- dfMaster.z.plot %>% 
  na.omit() %>%
  group_by(Genotype) %>%
  do( plot.code.All = ggplot(data = ., aes(x = as.factor(ID),y = DLPFC.z_wS)) +
        geom_quasirandom(color=unique(.$MyCol),width=0.2,size=2.2,aes(shape=as.factor(HC_pat) ) ) + 
        geom_errorbar(width= 0.3,size=1.1,
                      aes(y=meansDLPFC.z,ymin=meansDLPFC.z,ymax=meansDLPFC.z)) + 
        theme_bw() + 
        theme(legend.position='none',legend.key = element_blank()) + 
        scale_shape(name='',labels=c('HC    ','SZ')) +   
        theme(axis.title.y = element_blank())  + 
        theme(axis.title.x =element_blank()) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size=12) ) + 
        ylim(-2,3.5)
  )

########################
# Arrange into one plot
########################

plot.code.DLPFC.All$ROI<-"A.DLPFC"
plot.code.TC.All$ROI<-"B.TC"
plot.code.HIP.All$ROI<-"C.HIP"

#Merge all plot code into a tibble 
z.allPlotCode<-rbind(plot.code.DLPFC.All,plot.code.TC.All,plot.code.HIP.All)

#Rearrange order
z.allPlotCode<-z.allPlotCode %>%
  arrange(Genotype) %>%
  arrange(ROI)

#Layout for grid plot
layout.All<-rbind(c(1,2), 
                  c(3,4),
                  c(5,6))
#Save plot
gridPool<-gridExtra::grid.arrange(grobs=z.allPlotCode$plot.code.All,layout_matrix=layout.All)
ggsave(filename =paste0("./Results/Graphics/Pooled_VT_all_studies_row.svg"),plot = gridPool,device = svg,width =10, height = 30, units = "cm")

#######################################
# Merge scatterplot of studies and pooled plot
#######################################

#Save plot
gridGrand<-gridExtra::grid.arrange(grid.each.study,gridPool,ncol=2,widths=c(7,3))
ggsave(filename =paste0("./Results/Graphics/Grand_VT_scatter_plot.svg"),plot = gridGrand,device = svg,width =35, height = 25, units = "cm")


