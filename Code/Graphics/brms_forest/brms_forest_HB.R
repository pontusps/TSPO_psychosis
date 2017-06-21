#Forest plot using posteriors of random slopes fitted with the brms-package. 
#Written by M. Vuorre  https://mvuorre.github.io/post/2017/better-brms-forest-plots/
#Modified by Pontus P. Sigray, April 2017, Karolinska Institutet
#Requires David Robinson's and Ben Marwick's geom_flat_violin function to work: https://gist.github.com/benmarwick/2a1bb0133ff568cbe28d

bayes_forest_HB <- function (data, model, ROI, level = 0.95, xlim = NULL, ylim = c(-3,3), show_data = FALSE, 
          sort_estimates = FALSE, dens_fill = "#6497b1", dens_col = 'NA'){
  
  library(tidyverse)
  library(broom)
  source(paste0("./R/Graphics/","geom_flat_violin.R"))
  #bayesplot::color_scheme_set(scheme = "blue")
  
  if (!exists("geom_flat_violin", mode = "function")) {
    stop("geom_flat_violin needed for this function to work. Please source it.", 
         call. = FALSE)
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  if (!requireNamespace("broom", quietly = TRUE)) {
    stop("broom needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  
  data$study<-data$Study
  #Extract posteriors from brms-fit
  samps <- brms::posterior_samples(model, pars = c("b_", "r_"))
  
  #Keep only fixed effect and Study-slopes
  samps<-samps %>%
    select(-contains("Intercept") )
  
  names(samps) <- c("All", levels(data$study)) #
  samps.m <- as.matrix(samps)
  samps.m[, 2:dim(samps.m)[2]] <- samps.m[, 2:dim(samps.m)[2]] + 
    samps.m[, 1]
  sampsdf <- as.data.frame(samps.m)
  names(sampsdf) <- names(samps)
  
  #Make df into long format for plotting
  samps.l <- tidyr::gather_(sampsdf, key_col = quote(study), 
                            value_col = "value", gather_cols = names(sampsdf))
  
  #Summary stats of each posterior (mean, 95% CredInt): Consider using ´rethinking´ to get MAP and HPDI here instead
  samples_s <- group_by_(samps.l, "study")
  samples_s <- summarise_(samples_s, mean = ~mean(value), lwr = ~quantile(value,
                                                                          probs = 0.5 - level/2), upr = ~quantile(value, probs = 0.5 + level/2))
  
  samples_s$s <- paste0(round(samples_s$mean, 2), " [", round(samples_s$lwr,2), ", ", round(samples_s$upr, 2), "]")
  
  ##Get Study labels from original df
  samples_s$label<- samples_s$study
  
  if (sort_estimates) {
    samples_s <- arrange_(samples_s, "mean")
    samples_s$order <- 1:nrow(samples_s)
    samples_s$order <- ifelse(samples_s$study == "All", -Inf, 
                              samples_s$order)
    samples_s$study <- reorder(samples_s$study, samples_s$order)
  }else {
    order.indx<-data.frame(study=c("Kenk","Bloomfield","Coughlin","Hafizi","Collste"),order=c(5,4,3,2,1))
    samples_s <- merge(samples_s,order.indx,by="study",all.x = T)
    samples_s$order <- ifelse(samples_s$study == "All", -Inf, 
                              samples_s$order)
    samples_s$study <- reorder(samples_s$study, samples_s$order)
  }
  
  #Bold the "All" label
  xlabels<-c()
  xbreaks <- samples_s$study
  xlabels[1] <- expression(bold(All)) 
  xlabels[2:length(samples_s$study)] <- samples_s$label[2:length(samples_s$study)]
  
  #Plot forest. ylim must be set manually as input argument 
  p1 <- ggplot(samples_s, aes_(quote(study), quote(mean))) + 
    geom_linerange(aes_(ymin = quote(lwr), ymax = quote(upr)),
                   position = position_nudge(x = -0.1)) + 
    geom_flat_violin(data = samps.l, fill = dens_fill, col = dens_col, 
                     width = 0.95, alpha = 0.5, aes_(y = quote(value)),
                     position = position_nudge(x = -0.1)) + 
    geom_point(aes_(y = quote(mean)),position = position_nudge(x = -0.1)) + 
    geom_vline(xintercept = 1.5) + 
    geom_vline(xintercept = max(as.integer(samples_s$study)) + 0.5) + 
    geom_hline(yintercept = 0, lty = 2, size = 0.25,alpha = 0.6) + 
    geom_text(data = filter_(samples_s, quote(study !="All")), 
              hjust = "inward", vjust = "middle", aes_(label = quote(s), y = quote(ylim[2]))) + 
    geom_text(data = filter_(samples_s, quote(study == "All")), 
              hjust = "inward", vjust = "middle", aes_(label = quote(s), y = quote(ylim[2])), fontface = "bold") + 
    scale_x_discrete(breaks = samples_s$study, labels = xlabels, expand = c(0, 0.5)) + 
    scale_y_continuous(limits = ylim,expand = c(0, 0)) + 
    coord_flip() + 
    theme(axis.title.y = element_blank(), panel.border = element_rect(fill = NA, colour = NA, size = 0.6), 
          axis.line.x = element_line(size = 0.7), plot.title = element_text(face = "bold", hjust = 0.5), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA))
  
  #Show "raw-difference" between pats and controls within each study and the betas' corresponding 95% CI: 
  if (show_data) {
    
    #Center HAB and MAB VT values
    dat<-data %>% 
      group_by(Study,Genotype) %>%
      mutate(HIP.VT.cent = as.numeric(scale((HIP.VT),scale = F)),
             TC.VT.cent = as.numeric(scale((TC.VT),scale = F)),
             FC.VT.cent = as.numeric(scale((DLPFC.VT),scale = F)))
    
    #Z-scores centered VT values for HABs and MABs together
    dat<-dat %>% 
      group_by(Study) %>%
      mutate(HIP.VT.cent.z = as.numeric(scale((HIP.VT),scale = T)),
             TC.VT.cent.z = as.numeric(scale((TC.VT),scale = T)),
             FC.VT.cent.z = as.numeric(scale((DLPFC.VT),scale = T)))
    
    #Check which ROI is bering plotted
    if (ROI=="FC") { 
      
      #Obtain lm estimates
      coefs.FC<-dat %>% 
        group_by(Study) %>%
        do(tidy(lm(FC.VT.cent.z ~ HC_pat ,data = .))) %>%
        filter( !grepl("Intercept",term) ) %>%
        select(Study,estimate,std.error)
      
      #Merge with dataframe in ggplot2 object
      p1$data<-merge(p1$data,coefs.FC,by.x="study",by.y = "Study",all.x = T)
      
    }else if(ROI=="TC"){
      
      coefs.TC<-dat %>% 
        group_by(Study) %>%
        do(tidy(lm(TC.VT.cent.z ~ HC_pat,data = .))) %>%
        filter( !grepl("Intercept",term) ) %>%
        select(Study,estimate,std.error)
      
      #Merge with dataframe in ggplot2 object
      p1$data<-merge(p1$data,coefs.TC,by.x="study",by.y = "Study",all.x = T)
      
    }else if(ROI=="HIP"){
      
      coefs.HIP<-dat %>% 
        group_by(Study) %>%
        do(tidy(lm(HIP.VT.cent.z ~ HC_pat,data = .))) %>%
        filter( !grepl("Intercept",term) ) %>%
        select(Study,estimate,std.error)
      
      #Merge with dataframe in ggplot2 object<
      p1$data<-merge(p1$data,coefs.HIP,by.x="study",by.y = "Study",all.x = T)
      
    }else{
      stop("Specify ROI argument: FC, TC or HIP", 
           call. = FALSE)
    }
    
    #Add estimates and their  stand.error to plot
    p1 <- p1 + geom_point(size = 1.2, aes_(y = quote(estimate)), 
                          shape = 4, position = position_nudge(x = -0.2), na.rm = T) + 
      geom_linerange(size = 0.2, aes_(ymin = quote(estimate - std.error * 2 ), 
                                      ymax = quote(estimate + std.error * 2 )), 
                                      position = position_nudge(x = -0.2), na.rm = T)
  }
  return(p1)
}

