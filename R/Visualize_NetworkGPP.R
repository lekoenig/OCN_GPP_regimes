#################################################################################
#######   FUNCTIONS FOR USE IN VISUALIZATION/FIGURES: River NETWORK GPP   ####### 
#################################################################################


##----------------------------------------------##
##     Plot regimes from Savoy et al. (2019)     
##                Figure A.1
##----------------------------------------------##

  # Smoothing function for plotting reach-scale GPP regimes (from Phil Savoy and Jordan Read)
    smooth_5 <- function(data){
      smth <- data
      for (i in 2:(length(data)-3)){
        smth[i] <- mean(data[(i-2):(i+2)])
      }
      return(smth)
    }
    
  # Create data frame of smoothed reach-scale GPP regimes
    create.smoothed.dat <- function(data){
                  data <- data %>%
                          mutate(doy = c(1:365),
                                 medianGPP = apply(data,1,median,na.rm=TRUE),
                                 Quant25 = apply(data,1,quantile,probs=c(0.25),na.rm=TRUE),
                                 Quant75 = apply(data,1,quantile,probs=c(0.75),na.rm=TRUE))
                  data$medianGPP_smooth <- smooth_5(data$medianGPP)
                  data$Quant25_smooth <- smooth_5(data$Quant25)
                  data$Quant75_smooth <- smooth_5(data$Quant75)
                  return(data)
    }

  # Define plot for Figure A.1
    plot.regime <- function(dat,regime.title.chr){
      dat %>% 
        ggplot() + 
        geom_ribbon(aes(x=doy,ymin=Quant25_smooth,ymax=Quant75_smooth),fill='gray35',alpha=.4)+
        geom_line(aes(x=doy,y=medianGPP_smooth),color='black',size=0.9)+
        ggtitle(regime.title.chr)+coord_cartesian(ylim=c(0,10))+
        scale_x_continuous(breaks=c(0,100,200,300),labels=c(0,100,200,300))+
        scale_y_continuous(breaks=c(0,2,4,6,8,10),labels=c(0,2,4,6,8,10))+
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              plot.title = element_text(size = 13, face = "plain",family="serif"))
  }



##----------------------------------------------------##
##   Grab customized text/titles from indiv. figures        
##          for placement in a larger grid
##               used for Figures 1, 2
##----------------------------------------------------##

  get_yaxis <- function(plot){
    grobs <- plot_to_gtable(plot)$grobs
    yaxis <- grobs[[13]]
  }
  
  get_xaxis <- function(plot){
    grobs <- plot_to_gtable(plot)$grobs
    xaxis <- grobs[[12]]
  }
  
  get_ytext <- function(plot){
    grobs <- plot_to_gtable(plot)$grobs
    xaxis <- grobs[[3]]
  }
  
  get_xtext <- function(plot){
    grobs <- plot_to_gtable(plot)$grobs
    xaxis <- grobs[[7]]
  }
  
  get_arrow <- function(plot){
    grobs <- plot_to_gtable(plot)$grobs
    xaxis <- grobs[[6]]
  }
  
  get_ramp <- function(plot){
    grobs <- plot_to_gtable(plot)$grobs
    xaxis <- grobs[[1]]
  }



##----------------------------------------------##
##     Define panels for Figure 2 grid plot
##----------------------------------------------##

  # Define plotting function for 40 km^2 panels:
    plot.grids40 <- function(Sims){
      ribbon.colors <- rep("gray78",length(b.40))
      Sims %>% ggplot() + 
        geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=dfname))+
        geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname),size=0.7)+
        labs(x=expression(""),y=expression("")) +coord_cartesian(ylim=c(0,150))+
        scale_y_continuous(breaks=c(0,50,100,150))+
        scale_color_viridis(discrete=TRUE,direction = -1,option = "D",begin=0,end=1)+
        scale_fill_manual(values=ribbon.colors)+
        theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
              axis.text.x = element_blank(),axis.text.y=element_blank(),
              legend.position='none')
    }

  # Define plotting function for 160 km^2 panels:
    plot.grids160 <- function(Sims){
      ribbon.colors <- rep("gray78",length(b.160))
      Sims %>% ggplot() + 
        geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=dfname))+
        geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname),size=0.7)+
        labs(x=expression(""),y=expression("")) +coord_cartesian(ylim=c(0,600))+
        scale_y_continuous(breaks=c(0,200,400,600))+
        scale_color_viridis(discrete=TRUE,direction = -1,option = "D",begin=0,end=0.6)+
        scale_fill_manual(values=ribbon.colors)+
        theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
              axis.text.x = element_blank(),axis.text.y=element_blank(),
              legend.position='none')
    }

  # Define plotting function for 450 km^2 panels:
    plot.grids450 <- function(Sims){
      ribbon.colors <- rep("gray78",length(b.450))
      Sims %>% ggplot() + 
        geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=dfname))+
        geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname),size=0.7)+
        labs(x=expression(""),y=expression("")) +coord_cartesian(ylim=c(0,1800))+
        scale_y_continuous(breaks=c(0,600,1200,1800))+
        scale_color_viridis(discrete=TRUE,direction = -1,option = "D",begin=0,end=0.4)+
        scale_fill_manual(values=ribbon.colors)+
        theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
              axis.text.x = element_blank(),axis.text.y=element_blank(),
              legend.position='none')
    }

  # Define plotting function for 2,600 km^2 panels:
    plot.grids2600 <- function(Sims){
      ribbon.colors <- rep("gray78",length(b.2600))
      Sims %>% ggplot() + 
        geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=dfname))+
        geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname),size=0.7)+
        labs(x=expression(""),y=expression("")) +coord_cartesian(ylim=c(0,12300))+
        scale_y_continuous(breaks=c(0,4000,8000,12000))+
        scale_color_viridis(discrete=TRUE,direction = 1,option = "D",begin = 0,end=.8)+
        scale_fill_manual(values=ribbon.colors)+
        theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
              axis.text.x = element_blank(),axis.text.y=element_blank(),
              legend.position='none')
    }



##----------------------------------------------##
##      Define table A.1 in manuscript text 
##    (physical attributes of sub-catchments)
##----------------------------------------------##

  Fill.in.TableA.1 <- function(size,b.size){
    
    GPP.ReachArea <- data.frame(matrix(NA,nrow=(length(b.size)),ncol=7))
    colnames(GPP.ReachArea) <- c('WatershedSize','Rep','ReachArea_km2','Sum_length_m','Mean_length_m','DrainageArea_km2','Mean_width_m')
    
    for(i in 1:length(b.size)){
      c.ws <- get_ancestors(net2500,b.size[i])
      subws <- induced.subgraph(net2500,c.ws)
      GPP.ReachArea[i,1] <- size
      GPP.ReachArea[i,2] <- paste("Rep",i,sep="")
      GPP.ReachArea[i,3] <- (sum(V(subws)$ReachArea))/(1000*1000)
      GPP.ReachArea[i,4] <- sum((V(subws)$ReachArea/V(subws)$width)) # total length
      GPP.ReachArea[i,5] <- mean((V(subws)$ReachArea/V(subws)$width)) # mean length
      GPP.ReachArea[i,6] <- vertex.attributes(net2500,b.size[i])$DrainArea_km2
      GPP.ReachArea[i,7] <- mean(V(subws)$width) # mean width
    }
    return(GPP.ReachArea)
  }



##---------------------------------------------------------##
##           Define table 1 in manuscript text 
##---------------------------------------------------------##

  Fill.in.Table1.allsims <- function(df,chr.scenario,size,b.size){
    #  For example, Productive model scenario for a 450 km2 watershed:
    #  df = Sim450.P
    #  chr.scenario = "Productive"
    #  size = 450
    #  b.size = b.450
    
    # Summarize annual network GPP across simulations and sub-catchment reps:
    Sims = list()
    
    for(i in 1:length(b.size)){
      dat   <- as.data.frame(df[[i]])
      dat$Scenario <- chr.scenario
      dat$doy <- c(1:365)
      dat$dfname <- rep(paste("Rep",as.character(i),sep=""),365)
      dat_long <- tidyr::gather(dat,simulation.number,network.GPP,V1:V1000,factor_key=FALSE)
      
      Sims[[i]] <- dat_long
    }
    
    GPP.Sims = do.call(rbind, Sims)
    
    GPP.Sims.df <- GPP.Sims %>%
      group_by(Scenario,dfname,simulation.number) %>%
      summarize(AnnualGPP_kgO2.yr = sum(network.GPP),
                MeanGPP_kgO2.d = mean(network.GPP)) %>%
      mutate(AnnualGPP_kgC.yr = (AnnualGPP_kgO2.yr*12/32),
             MeanGPP_kgC.d = (MeanGPP_kgO2.d*12/32))
    
    # Gather network structure information (e.g. benthic area, length):
    GPP.ReachArea <- data.frame(matrix(NA,nrow=(length(b.size)),ncol=7))
    colnames(GPP.ReachArea) <- c('WatershedSize','Rep','ReachArea_km2','Sum_length_m','Mean_length_m','DrainageArea_km2','Mean_width_m')
    
    for(i in 1:length(b.size)){
      c.ws <- get_ancestors(net2500,b.size[i])
      subws <- induced.subgraph(net2500,c.ws)
      GPP.ReachArea[i,1] <- size
      GPP.ReachArea[i,2] <- paste("Rep",i,sep="")
      GPP.ReachArea[i,3] <- (sum(V(subws)$ReachArea))/(1000*1000)
      GPP.ReachArea[i,4] <- sum(V(subws)$ReachArea/V(subws)$width) # total length
      GPP.ReachArea[i,5] <- mean(V(subws)$ReachArea/V(subws)$width) # mean length
      GPP.ReachArea[i,6] <- vertex.attributes(net2500,b.size[i])$DrainArea_km2
      GPP.ReachArea[i,7] <- mean(V(subws)$width) # mean width
    }
    
    Results.Table1 <- left_join(GPP.ReachArea,GPP.Sims.df,by=c("Rep"="dfname")) %>%
      mutate(MeanDailyGPP_gCm2d = (MeanGPP_kgC.d*1000)/(ReachArea_km2*1000*1000))
    
    GPP.cum.sums <- GPP.Sims %>% group_by(Scenario,dfname,simulation.number) %>%
      mutate(cumGPP_kgO2 = cumsum(network.GPP),
             cumPropGPP = cumGPP_kgO2/max(cumGPP_kgO2)) %>%
      summarize(CHK_AnnualGPP_kgO2yr = max(cumGPP_kgO2),
                doy_50pct = min(which(cumPropGPP>0.5)),
                doy_max = which.max(network.GPP))
    
    Results.Table <- left_join(Results.Table1, GPP.cum.sums, by = c("Rep" = "dfname", "simulation.number" = "simulation.number","Scenario" = "Scenario"))
    
    if(all.equal(Results.Table$AnnualGPP_kgO2.yr,Results.Table$CHK_AnnualGPP_kgO2yr)=="TRUE") {Results.Table <- Results.Table} else{print("warning: annual sums not equal")}
    
    return(Results.Table)
  }



##---------------------------------------------------------##
##     Create data subset for plotting network regimes  
##          (medians and 95% confidence intervals)     
##---------------------------------------------------------##

  Subset.replicate.regimes <- function(df,chr.scenario,which.reps,b.size){
    #  For example, Productive model scenario for a 450 km2 watershed:
    #  df = Sim450.P
    #  chr.scenario = "Productive"
    #  which.reps = 1:10
    #  b.size = b.450
    
    Reps <- paste("Rep",which.reps,sep="")
    
    ## Summarize annual network GPP across simulations and sub-catchment reps:
    Sims = list()
    
    for(i in 1:length(b.size)){
      dat   <- as.data.frame(df[[i]])
      dat$Scenario <- chr.scenario
      dat$doy <- c(1:365)
      dat$dfname <- rep(paste("Rep",as.character(i),sep=""),365)
      dat_long <- tidyr::gather(dat,simulation.number,network.GPP,V1:V1000,factor_key=FALSE)
      
      Sims[[i]] <- dat_long
    }
    
    GPP.Sims = do.call(rbind, Sims)
    
    # Calculate median and 95% confidence intervals for 1,000 simulations x each rep:
    GPPSims.plot <- GPP.Sims %>% group_by(Scenario,dfname,doy) %>%
      summarize(meanvalue = mean(network.GPP),
                Q2.5 = quantile(network.GPP,probs=0.025),
                Q97.5 = quantile(network.GPP,probs=0.975))
    GPPSims.plot <- filter(GPPSims.plot,dfname %in% Reps)
    
    return(GPPSims.plot)
  }



##---------------------------------------------------------##
##     Create data subset for plotting network regimes
##                Riparian clearing scenario
##          (medians and 95% confidence intervals)     
##---------------------------------------------------------##

  Subset.RC.regimes <- function(df,chr.scenario,b.size){
    #  For example, Productive model scenario for a 2600 km2 watershed:
    #  df = Sim2600.P.RC
    #  chr.scenario = "Productive"
    #  size = 2600
    #  b.size = b.2600[10] # only using full network
    
    # Summarize annual network GPP across simulations and sub-catchment reps:
    Sims = list()
    
    for(i in 1:length(rc.sims)){
      dat   <- as.data.frame(df[[i]])
      dat$Scenario <- chr.scenario
      dat$dfname <- rep(paste("subsample",rc.sims[i],sep="."),365)
      dat$doy <- c(1:365)
      dat_long <- tidyr::gather(dat,simulation.number,network.GPP,V1:V1000,factor_key=FALSE)
      
      Sims[[i]] <- dat_long
    }
    
    GPP.Sims = do.call(rbind, Sims)
    
    # Calculate median and 95% confidence intervals for 1,000 simulations x each rep:
    GPPSims.plot <- GPP.Sims %>% 
      group_by(Scenario,dfname,doy) %>%
      summarize(meanvalue = mean(network.GPP),
                Q2.5 = quantile(network.GPP,probs=0.025),
                Q97.5 = quantile(network.GPP,probs=0.975)) 
    
    return(GPPSims.plot)
  }
  
  

##--------------------------------------------##
##     Define table A.3 in manuscript text       
##--------------------------------------------##

  Fill.in.TableA3.allsims <- function(df,chr.scenario,size,b.size){
    #  For example, Productive model scenario for a 2600 km2 watershed:
    #  df = Sim2600.P.RC
    #  chr.scenario = "Productive"
    #  size = 2600
    #  b.size = b.2600[10] # only using full network
    
    # Summarize annual, network GPP across simulations and sub-catchment reps:
    Sims = list()
    
    for(i in 1:length(rc.sims)){
      dat   <- as.data.frame(df[[i]])
      dat$Scenario <- chr.scenario
      dat$dfname <- rep(paste("subsample",rc.sims[i],sep="."),365)
      dat$doy <- c(1:365)
      dat_long <- tidyr::gather(dat,simulation.number,network.GPP,V1:V1000,factor_key=FALSE)
      
      Sims[[i]] <- dat_long
    }
    
    GPP.Sims = do.call(rbind, Sims)
    
    GPP.Sims.df <- GPP.Sims %>%
      group_by(Scenario,dfname,simulation.number) %>%
      summarize(AnnualGPP_kgO2.yr = sum(network.GPP)) %>%
      mutate(AnnualGPP_kgC.yr = (AnnualGPP_kgO2.yr*12/32))
    
    return(GPP.Sims.df)
  }




