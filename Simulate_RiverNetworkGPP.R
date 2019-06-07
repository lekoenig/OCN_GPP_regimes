
## Analysis steps for simulating river-network GPP using optimal channel networks
## Last updated 26 February 2019


library(igraph)
library(tidyr)
library(forcats)
library(broom)
library(dplyr)
library(ggplot2)
library(cowplot)
library(knitr)
library(scales)
library(purrr)
library(randomizr)
library(viridis)
library(raster)
library(GGally)
library(gridGraphics)
library(grid)
library(animation)
#install.packages("colourschemes", repos="http://R-Forge.R-project.org")   # note that R package "colourschemes" now has to be downloaded from R-Forge
library(colourschemes)
library(intergraph)
library(zoo)


# Analysis script calls on functions in R folder:
source("./R/NetworkGPP_Functions.R")
source("./R/Visualize_NetworkGPP.R")


##-----------------------------------------------------##
##     INPUT DATA - Optimal Channel Networks (OCNs)     
##-----------------------------------------------------##

  # For this study, we used a 512x512 pixel optimal channel network (OCN) to simulate a 2,621 km^2 river network.
  # This OCN ("thres 50") corresponds to a drainage density of approximately 1 km/km^2

  # Load OCN thres file (used for plotting networks):
  ocn512<- read.table("data/OCN_512_512_thres_50.txt",header=FALSE)
  colnames(ocn512) <- c("ID","X","Y","Downstream_ID","Length","Drainage_Area","IDreach","Direct_DrainageArea_Pixel","Direct_DrainageArea_Reach")

  # Load OCN reach_ntw file:
  net_2500 <- read.csv("data/OCN_512_512_thres_50_reach_ntw.csv",header=T,stringsAsFactors=FALSE)
  #net_2500$Width.Empirical <- 0.001772556*net_2500$DRAINAREA_m2^0.46326  ## NOT USED: width-area relationship all 47 sites
  net_2500$Width.Empirical <- 0.0013*net_2500$DRAINAREA_m2^0.479          ## width-area relationship excluding nwis_11273400 from the relationship
  net_2500$ReachArea <- net_2500$Width.Empirical * net_2500$ReachLength_m
  
  # define nodes:
  net_2500_nodes <- net_2500[,c("ID","XCOORD","YCOORD","Width.Empirical","ReachArea","DRAINAREA_m2")] 
  colnames(net_2500_nodes) <- c("ID","x","y","width","ReachArea","DrainArea_m2")
  net_2500_nodes$DrainArea_km2 <- net_2500_nodes$DrainArea_m2/(1000*1000)

  # define links:
  net_2500_links <- net_2500[-nrow(net_2500),c("ID","TO_ID")]
  colnames(net_2500_links) <- c("FROM_ID","TO_ID")

  # define graph and create an igraph object that represents the stream channel network:
  net2500 <- graph.data.frame (d=net_2500_links, vertices=net_2500_nodes, directed=T)  
  V(net2500)$dist.from.outlet <- distances(net2500, v=V(net2500)[length(V(net2500)$name)], to=V(net2500), weights=NA)
  V(net2500)$ReachArea[length(V(net2500)$name)] <- 1  # define area of most-downstream reach (1 m2)


##-----------------------------------------------------##
##     INPUT DATA - Reach-scale productivity regimes 
##              from Savoy et al. (2019)     
##-----------------------------------------------------##

  # Load site information:
  SiteInfo <- read.csv("data/Savoy_hydroshare/site_basic.csv",header=T,stringsAsFactors = FALSE)
  SiteInfo.Table <- SiteInfo %>% dplyr::select(Site_ID,Site_name,Lat,Lon,WS_area,Width,four_clus)
  
  # Import gap-filled productivity regimes for each site following a given regime (downloaded from CUAHSI Hydroshare: https://www.hydroshare.org/resource/eba152073b4046178d1a2ffe9a897ebe/).
  site.regimes <- read.csv("data/Savoy_hydroshare/avg_gpp_filled.csv",header=T,stringsAsFactors = FALSE)    # rows are sites and columns are day of year
  
  # Which regime does each site (identified by NWIS ID) belong to?
  sites_summerdecline <- which(colnames(site.regimes) %in% (filter(SiteInfo.Table,four_clus=="summer decline")$Site_ID))
  sites_springpeak <- which(colnames(site.regimes) %in% (filter(SiteInfo.Table,four_clus=="spring peak")$Site_ID))
  sites_summerpeak <- which(colnames(site.regimes) %in% (filter(SiteInfo.Table,four_clus=="summer peak")$Site_ID))
  sites_aseasonal <- which(colnames(site.regimes) %in% (filter(SiteInfo.Table,four_clus=="aseasonal")$Site_ID))
  
  # Create subsets of the data that correspond with each regime:
  summerdecline.dat <- site.regimes[,sites_summerdecline]
  springpeak.dat <- site.regimes[,sites_springpeak]
  summerpeak.dat <- site.regimes[,sites_summerpeak]
  aseasonal.dat <- site.regimes[,sites_aseasonal]
  
  # Find the probability of a given reach-scale productivity regime with each width quantile of the empirical dataset (informs Stochastic model scenario):
    #Find quantiles of width distribution:
      quantile(SiteInfo$Width)
    
    #Assign each site to a width quantile:
      SiteInfo$Width_quantile <- assign.width.class(SiteInfo$Width)
    
    #Find probabilities of each regime (from column "four_clus") within each width quantile (printed in Table A.2):
      Stoch.prob <- SiteInfo %>% 
        group_by(Width_quantile,four_clus) %>%
        summarise (n = n()) %>%
        mutate(freq = (n / sum(n))*100)
      print(Stoch.prob)

    
##-----------------------------------------------------##
##     Calculate river-network productivity (GPP)    
##-----------------------------------------------------##

  # Calculate network-scale GPP for a range of watershed size given three modeled scenarios: Productive rivers, Unproductive rivers, and Stochastic.

  # Define the number of simulations to use throughout (used across scenarios/watershed size simulations):
  N=1000    

  
#--------- Analyze GPP for 40 km^2 sub-catchments ---------   

  # Define sub-catchment replicates:
  b.40 <- which(net_2500_nodes$DrainArea_km2<42.1 & net_2500_nodes$DrainArea_km2>37.9)
  b.40 <- b.40[c(1,2,5,6,8,9,11,12,14,15,16,17,18)]
  
  # Define empty lists to house annual, network GPP estimates for each replicate 40 km^2 sub-catchment:
  Sim40.P = list()
  Sim40.A = list()
  Sim40.S = list()    
  
  # Perform network calculations for all replicate 40 km^2 sub-catchments:   
  for(j in 1:length(b.40)){
    
    c.40 <- get_ancestors(net2500,b.40[j])
      subws <- induced.subgraph(net2500,c.40)
    
    set.seed(j)
    
    stoch.40 <- data.frame(V(subws)$name,V(subws)$width,V(subws)$ReachArea)
    colnames(stoch.40) <- c('node','width','ReachArea') 
    stoch.40$width.class <- assign.width.class(stoch.40$width)
    
    ifelse(length(unique(stoch.40$width.class))==1,
           block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
           ifelse(length(unique(stoch.40$width.class))==2,
                  block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                        c(0.0,0.18,0.46,0.36)),
                  ifelse(length(unique(stoch.40$width.class))==3,
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33)),
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33),
                                               c(0.67,0.25,0.08,0.0)))))
    
    stoch.40$Regime <- block_ra(blocks = stoch.40$width.class, 
                                block_prob_each = block_m_each,
                                conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
    
    Network.Simulations.Prod <- replicate(N,Network.Prod(network = subws))
    Network.Simulations.Unprod <- replicate(N,Network.Unprod(network = subws))
    Network.Simulations.Stoch <- replicate(N, Network.Stoch(network = stoch.40))
    
    Sim40.P[[j]] = Network.Simulations.Prod
    Sim40.A[[j]] = Network.Simulations.Unprod
    Sim40.S[[j]] = Network.Simulations.Stoch
    
  } 

#--------- Analyze GPP for 160 km^2 sub-catchments ---------   

  # Define sub-catchment replicates:
  b.160 <- which(net_2500_nodes$DrainArea_km2<168.1 & net_2500_nodes$DrainArea_km2>151.9)
  b.160 <- b.160[c(2,3,6,11,14)]

  # Define empty lists to house annual, network GPP estimates for each replicate 160 km^2 sub-catchment:
  Sim160.P = list()
  Sim160.A = list()
  Sim160.S = list()    
  
  # Perform network calculations for all replicate 160 km^2 sub-catchments:   
  for(j in 1:length(b.160)){
    
    c.160 <- get_ancestors(net2500,b.160[j])
    subws <- induced.subgraph(net2500,c.160)
    
    set.seed(j)
    
    stoch.160 <- data.frame(V(subws)$name,V(subws)$width,V(subws)$ReachArea)
    colnames(stoch.160) <- c('node','width','ReachArea') 
    stoch.160$width.class <- assign.width.class(stoch.160$width)
  
    ifelse(length(unique(stoch.160$width.class))==1,
           block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
           ifelse(length(unique(stoch.160$width.class))==2,
                  block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                        c(0.0,0.18,0.46,0.36)),
                  ifelse(length(unique(stoch.160$width.class))==3,
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33)),
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33),
                                               c(0.67,0.25,0.08,0.0)))))
    
    stoch.160$Regime <- block_ra(blocks = stoch.160$width.class, 
                                 block_prob_each = block_m_each,
                                 conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
    
    Network.Simulations.Prod <- replicate(N,Network.Prod(network = subws))
    Network.Simulations.Unprod <- replicate(N,Network.Unprod(network = subws))
    Network.Simulations.Stoch <- replicate(N,Network.Stoch(network = stoch.160))
    
    Sim160.P[[j]] = Network.Simulations.Prod
    Sim160.A[[j]] = Network.Simulations.Unprod
    Sim160.S[[j]] = Network.Simulations.Stoch
    
  }

#--------- Analyze GPP for 450 km^2 sub-catchments ---------   

  # Define sub-catchment replicates:
  b.450 <- which(net_2500_nodes$DrainArea_km2<472.6 & net_2500_nodes$DrainArea_km2>427.4)
  b.450 <- b.450[c(1,8)]

  # Define empty lists to house annual, network GPP estimates for each replicate 450 km^2 sub-catchment:
  Sim450.P = list()
  Sim450.A = list()
  Sim450.S = list()    
  
  # Perform network calculations for all replicate 450 km^2 sub-catchments:   
  for(j in 1:length(b.450)){
    
    c.450 <- get_ancestors(net2500,b.450[j])
    subws <- induced.subgraph(net2500,c.450)
    
    set.seed(j)
    
    stoch.450 <- data.frame(V(subws)$name,V(subws)$width,V(subws)$ReachArea)
    colnames(stoch.450) <- c('node','width','ReachArea') 
    stoch.450$width.class <- assign.width.class(stoch.450$width)
    
    ifelse(length(unique(stoch.450$width.class))==1,
           block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
           ifelse(length(unique(stoch.450$width.class))==2,
                  block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                        c(0.0,0.18,0.46,0.36)),
                  ifelse(length(unique(stoch.450$width.class))==3,
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33)),
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33),
                                               c(0.67,0.25,0.08,0.0)))))
    
    stoch.450$Regime <- block_ra(blocks = stoch.450$width.class, 
                                 block_prob_each = block_m_each,
                                 conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
    
    Network.Simulations.Prod <- replicate(N,Network.Prod(network = subws))
    Network.Simulations.Unprod <- replicate(N,Network.Unprod(network = subws))
    Network.Simulations.Stoch <- replicate(N,Network.Stoch(network = stoch.450))
    
    Sim450.P[[j]] = Network.Simulations.Prod
    Sim450.A[[j]] = Network.Simulations.Unprod
    Sim450.S[[j]] = Network.Simulations.Stoch
    
  }
  
#--------- Analyze GPP for 2600 km^2 sub-catchments ---------   

  # Define sub-catchment replicates:
  b.2600 <- which(net_2500_nodes$DrainArea_km2<2752 & net_2500_nodes$DrainArea_km2>2490)
  b.2600 <- b.2600[10]
  
  # Define empty lists to house annual network GPP estimates for each replicate 2600 km^2 sub-catchment:
  Sim2600.P = list()
  Sim2600.A = list()
  Sim2600.S = list()    
  
  # Perform network calculations for all replicate 2600 km^2 sub-catchments:   
  for(j in 1:length(b.2600)){
    
    c.2600 <- get_ancestors(net2500,b.2600[j])
    subws <- induced.subgraph(net2500,c.2600)
    
    set.seed(j)
    
    stoch.2600 <- data.frame(V(subws)$name,V(subws)$width,V(subws)$ReachArea)
    colnames(stoch.2600) <- c('node','width','ReachArea') 
    stoch.2600$width.class <- assign.width.class(stoch.2600$width)
    
    ifelse(length(unique(stoch.2600$width.class))==1,
           block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
           ifelse(length(unique(stoch.2600$width.class))==2,
                  block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                        c(0.0,0.18,0.46,0.36)),
                  ifelse(length(unique(stoch.2600$width.class))==3,
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33)),
                         block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                               c(0.0,0.18,0.46,0.36),
                                               c(0.25,0.09,0.33,0.33),
                                               c(0.67,0.25,0.08,0.0)))))
    
    stoch.2600$Regime <- block_ra(blocks = stoch.2600$width.class, 
                                  block_prob_each = block_m_each,
                                  conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
    
    Network.Simulations.Prod <- replicate(N,Network.Prod(network = subws))
    Network.Simulations.Unprod <- replicate(N,Network.Unprod(network = subws))
    Network.Simulations.Stoch <- replicate(N,Network.Stoch(network = stoch.2600))
    
    Sim2600.P[[j]] = Network.Simulations.Prod
    Sim2600.A[[j]] = Network.Simulations.Unprod
    Sim2600.S[[j]] = Network.Simulations.Stoch
  } 
    
    
##---------------------------------------------------##
##   Create Table 1: mean and sd across simulations       
##---------------------------------------------------##  
  
  # Call function to populate table for the three scenarios (Productive rivers, Stochastic, Unproductive rivers) and each watershed size (40, 160, 450, 2600):
  P.40.table <- Fill.in.Table1.allsims(Sim40.P,"Productive",40,b.40)
  S.40.table <- Fill.in.Table1.allsims(Sim40.S,"Stochastic",40,b.40)
  A.40.table <- Fill.in.Table1.allsims(Sim40.A,"Unproductive",40,b.40)
  
  P.160.table <- Fill.in.Table1.allsims(Sim160.P,"Productive",160,b.160)
  S.160.table <- Fill.in.Table1.allsims(Sim160.S,"Stochastic",160,b.160)
  A.160.table <- Fill.in.Table1.allsims(Sim160.A,"Unproductive",160,b.160)
  
  P.450.table <- Fill.in.Table1.allsims(Sim450.P,"Productive",450,b.450)
  S.450.table <- Fill.in.Table1.allsims(Sim450.S,"Stochastic",450,b.450)
  A.450.table <- Fill.in.Table1.allsims(Sim450.A,"Unproductive",450,b.450)
  
  P.2600.table <- Fill.in.Table1.allsims(Sim2600.P,"Productive",2600,b.2600)
  S.2600.table <- Fill.in.Table1.allsims(Sim2600.S,"Stochastic",2600,b.2600)
  A.2600.table <- Fill.in.Table1.allsims(Sim2600.A,"Unproductive",2600,b.2600)
  
  Table1.allsims <- bind_rows(P.40.table,S.40.table,A.40.table,
                              P.160.table,S.160.table,A.160.table,
                              P.450.table,S.450.table,A.450.table,
                              P.2600.table,S.2600.table,A.2600.table)
  
  # Manuscript Table 1:
  Table1 <- Table1.allsims %>%
    group_by(WatershedSize,Scenario) %>%
    summarize(MeanAnnualGPP.kgCyr = mean(AnnualGPP_kgC.yr),
              sd.AnnualGPP.kgCyr = sd(AnnualGPP_kgC.yr),
              MeanDailyGPP.gCm2d = mean(MeanDailyGPP_gCm2d),
              sd.MeanDailyGPP.gCm2d = sd(MeanDailyGPP_gCm2d),
              Mean.doy50pct = mean(doy_50pct),
              sd.Mean.doy50pct = sd(doy_50pct),
              doy.maxGPP = mean(doy_max),
              sd.doymaxGPP = sd(doy_max)) %>%
    mutate(CV.AnnualGPP.kgCyr = (sd.AnnualGPP.kgCyr/MeanAnnualGPP.kgCyr)*100)
  print(Table1)
  
  # Report mean annual GPP (kg C yr^-1) in Table 1 as 3 significant digits:
  print(formatC(signif(Table1$MeanAnnualGPP.kgCyr,digits=3), digits=3,format="fg", flag="#"))
  
  # Report the standard deviation of mean annual GPP (kg C yr^-1) in Table 1 as 3 significant digits:
  print(formatC(signif(Table1$sd.AnnualGPP.kgCyr,digits=3), digits=3,format="fg", flag="#"))
  

##------------------------------------------------------------------------##
##   Figure 1: envelope of productivity regimes for 2,621 km^2 watershed         
##------------------------------------------------------------------------##  

  # Plot network-scale GPP regimes and the spatial distribution of GPP on day of peak productivity for all three modeled scenarios.  
  
  # 1) Generate spatial patterns in reach-scale GPP (g O2 m^-2 d^-1): 
  V(net2500)$width.class <- assign.width.class(V(net2500)$width)
  
  ifelse(length(unique(V(net2500)$width.class))==1,
         block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
         ifelse(length(unique(V(net2500)$width.class))==2,
                block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                      c(0.0,0.18,0.46,0.36)),
                ifelse(length(unique(V(net2500)$width.class))==3,
                       block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                             c(0.0,0.18,0.46,0.36),
                                             c(0.25,0.09,0.33,0.33)),
                       block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                             c(0.0,0.18,0.46,0.36),
                                             c(0.25,0.09,0.33,0.33),
                                             c(0.67,0.25,0.08,0.0)))))
  
  V(net2500)$Regime <- block_ra(blocks = V(net2500)$width.class, 
                                block_prob_each = block_m_each,
                                conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
  
  for(j in 1:length(V(net2500))){
    V(net2500)$GPP.doy.P[j] <- ifelse(V(net2500)$width[j]>9,sample(summerpeak.dat,size=1,replace=TRUE)[207,1],
                                      sample(springpeak.dat,size=1,replace=TRUE)[207,1])
    
    V(net2500)$GPP.doy.A[j] <- ifelse(V(net2500)$width[j]>9,sample(aseasonal.dat,size=1,replace=TRUE)[95,1],
                                      sample(springpeak.dat,size=1,replace=TRUE)[95,1])
    
    V(net2500)$GPP.doy.S[j] <- ifelse(V(net2500)$Regime[j] == "SummerPeak",
                                      sample(summerpeak.dat,size=1,replace=TRUE)[109,1],
                                      ifelse(V(net2500)$Regime[j] == "Aseasonal",
                                             sample(aseasonal.dat,size=1,replace=TRUE)[109,1],
                                             ifelse(V(net2500)$Regime[j] == "SummerDecline",
                                                    sample(summerdecline.dat,size=1,replace=TRUE)[109,1],
                                                    sample(springpeak.dat,size=1,replace=TRUE)[109,1])))
  }
  
  # 2) Plot spatial patterns in reach-scale GPP:
  tramp = multiRamp(rbind(c(-3.2,2.3)),
                    list(c("darkgoldenrod2","#014636")))
  trampA = multiRamp(rbind(c(-3.2,1.1)),
                     list(c("darkgoldenrod2","#014636")))
  trampS = multiRamp(rbind(c(-3.2,2.1)),
                     list(c("darkgoldenrod2","#014636")))
  
  my_colors.P = tramp(log(V(net2500)$GPP.doy.P*12/32))
  my_colors.A = trampA(log(V(net2500)$GPP.doy.A*12/32))
  my_colors.S = trampS(log(V(net2500)$GPP.doy.S*12/32))
  
  my_area = 0.9*log(V(net2500)$width) 
  
  V(net2500)$my_colors.P <- my_colors.P
  V(net2500)$my_colors.A <- my_colors.A
  V(net2500)$my_colors.S <- my_colors.S
  
  # color streamlines by GPP:
  list.2500.GPP <- data.frame(V(net2500)$name,V(net2500)$width,V(net2500)$ReachArea,V(net2500)$GPP.doy.P,V(net2500)$GPP.doy.A,V(net2500)$GPP.doy.S,V(net2500)$my_colors.P,V(net2500)$my_colors.A,V(net2500)$my_colors.S)
  colnames(list.2500.GPP) <- c('node','width','ReachArea','GPP.doy.P','GPP.doy.A','GPP.doy.S',
                               'my_colors.P','my_colors.A','my_colors.S') 
  
  list.2500.GPP$IDreach     <- as.numeric(as.character(list.2500.GPP$node))
  list.2500.GPP$my_colors.P <- as.character(list.2500.GPP$my_colors.P)
  list.2500.GPP$my_colors.A <- as.character(list.2500.GPP$my_colors.A)
  list.2500.GPP$my_colors.S <- as.character(list.2500.GPP$my_colors.S)
  
  ocn512.GPP <- left_join(ocn512,list.2500.GPP,by="IDreach")
  
  # Plot spatial arrangement of GPP throughout the network:
  par(mfrow=c(3,1),mar=c(0.75,0.3,0,0.3),bg="transparent")  #mar: bottom,left,top,right
  plot(ocn512.GPP$X, ocn512.GPP$Y, pch=16, col=ocn512.GPP$my_colors.P, cex=0.09*log(ocn512.GPP$Drainage_Area), xlab="", ylab='', axes=F)
  
  plot(ocn512.GPP$X, ocn512.GPP$Y, pch=16, col=ocn512.GPP$my_colors.S, cex=0.09*log(ocn512.GPP$Drainage_Area), xlab="", ylab='', axes=F)
  
  plot(ocn512.GPP$X, ocn512.GPP$Y, pch=16, col=ocn512.GPP$my_colors.A, cex=0.09*log(ocn512.GPP$Drainage_Area), xlab="", ylab='', axes=F)
  grid.echo()
  net.grids <- grid.grab()     # from grid/gridGraphics packages

  # 3) Plot network-scale regimes: 
  GPP.Fig1b.P <- Subset.replicate.regimes(df=Sim2600.P,chr.scenario="Productive",which.reps=1,b.size=b.2600)
  GPP.Fig1b.S <- Subset.replicate.regimes(df=Sim2600.S,chr.scenario="Stochastic",which.reps=1,b.size=b.2600)
  GPP.Fig1b.A <- Subset.replicate.regimes(df=Sim2600.A,chr.scenario="Unproductive",which.reps=1,b.size=b.2600)
  
  Fig1b.P <- GPP.Fig1b.P %>% 
    ggplot() + 
    geom_ribbon(aes(x=doy, ymin=Q2.5*12/32, ymax=Q97.5*12/32),fill="gray78")+
    geom_vline(xintercept = 207,color="gray25",linetype="dashed")+
    geom_line(aes(x=doy,y=meanvalue*12/32),size=.9,color="#014636")+ 
    coord_cartesian(ylim=c(0,12100))+labs(x="Day of year",y=expression(River~network~GPP~(kg~C~d^-1)))+
    scale_y_continuous(breaks=c(0,4000,8000,12000))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    annotate("text", label = "a)", x = 15, y = 11948.8, size = 5.25)+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  Fig1b.S <- GPP.Fig1b.S %>% 
    ggplot() + 
    geom_ribbon(aes(x=doy, ymin=Q2.5*12/32, ymax=Q97.5*12/32),fill="gray78")+
    geom_vline(xintercept = 109,color="gray25",linetype="dashed")+
    geom_line(aes(x=doy,y=meanvalue*12/32),size=.9,color="#014636")+
    coord_cartesian(ylim=c(0,12100))+labs(x="Day of year",y=expression(River~network~GPP~(kg~C~d^-1)))+
    scale_y_continuous(breaks=c(0,4000,8000,12000))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    annotate("text", label = "b)", x = 15, y = 11948.8, size = 5.25)+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  Fig1b.A <- GPP.Fig1b.A %>% 
    ggplot() + 
    geom_ribbon(aes(x=doy, ymin=Q2.5*12/32, ymax=Q97.5*12/32),fill="gray78")+
    geom_vline(xintercept = 95,color="gray25",linetype="dashed")+
    geom_line(aes(x=doy,y=meanvalue*12/32),size=.9,color="#014636")+
    coord_cartesian(ylim=c(0,12100))+labs(x="Day of year",y=expression(River~network~GPP~(kg~C~d^-1)))+
    scale_y_continuous(breaks=c(0,4000,8000,12000))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    annotate("text", label = "c)", x = 15, y = 11948.8, size = 5.25)+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())
    
  # Define legends and lab axis text/titles:
  Fig1b.lab <- GPP.Fig1b.A %>% ggplot() + 
    geom_ribbon(aes(x=doy, ymin=Q2.5*12/32, ymax=Q97.5*12/32),fill="gray78",inherit.aes = FALSE)+
    geom_point(aes(x=doy,y=meanvalue*12/32),size=0.8,color="#014636")+
    coord_cartesian(ylim=c(0,12100))+
    labs(x="Day of year",y=expression(River~network~GPP~(kg~C~d^-1)))+
    scale_x_continuous(breaks=c(0,100,200,300))
  
  legend.lab <- GPP.Fig1b.A %>% ggplot() + 
    geom_ribbon(aes(x=doy, ymin=Q2.5*12/32, ymax=Q97.5*12/32),fill="gray78",inherit.aes = FALSE)+
    geom_point(aes(x=doy,y=meanvalue*12/32),size=0.8,color="#014636")+
    coord_cartesian(ylim=c(0,12100))+
    labs(y="Day of year",x=expression(GPP~(g~C~m^-2~d^-1)))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    theme(axis.title.x=element_text(size=11))
  
  Fig1b.axistext <- GPP.Fig1b.A %>% ggplot() + 
    geom_ribbon(aes(x=doy, ymin=Q2.5*12/32, ymax=Q97.5*12/32),fill="gray78",inherit.aes = FALSE)+
    geom_point(aes(x=doy,y=meanvalue*12/32),size=0.8,color="#014636")+
    coord_cartesian(ylim=c(0,12100))+
    labs(x="Day of year",y=expression(River~network~GPP~(kg~C~d^-1)))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),             
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  # Extract titles/text needed for grid figure:
  doy_labels   <- get_xtext(Fig1b.axistext)
  title.x      <- get_xaxis(Fig1b.lab)
  title.y      <- get_yaxis(Fig1b.lab)
  y.text       <- get_ytext(Fig1b.axistext)
  legend.title <- get_xaxis(legend.lab)
  
  labD = grid::grobTree(textGrob("d)", gp=gpar(col="black", fontsize=14, fontface="plain")))
  labE = grid::grobTree(textGrob("e)", gp=gpar(col="black", fontsize=14, fontface="plain")))
  labF = grid::grobTree(textGrob("f)", gp=gpar(col="black", fontsize=14, fontface="plain")))
  
  # 4) Build multi-panel plot:
  title <- grid::grobTree(textGrob("Productive rivers", gp=gpar(col="black", fontsize=14, fontface="bold")))
  title.S <- grid::grobTree(textGrob("Stochastic", gp=gpar(col="black", fontsize=14, fontface="bold")))
  title.A <- grid::grobTree(textGrob("Unproductive rivers", gp=gpar(col="black", fontsize=14, fontface="bold")))
  
  Fig1.bare <- plot_grid(Fig1b.P,Fig1b.S,Fig1b.A,ncol=1,align="hv")
  
  Fig1.combo <- plot_grid(Fig1.bare,NULL,net.grids,NULL,ncol=4,rel_widths=c(.88,.05,.7,.05))
  print(Fig1.combo)
  
  Fig1.xlab <- plot_grid(Fig1.combo,NULL,ncol=1,rel_heights = c(1,0.05))
  Fig1.labs <- plot_grid(NULL,Fig1.xlab,ncol=2,rel_widths = c(0.06,1))
  
  Fig1.print <- (Fig1.labs + 
                   draw_grob(title.x, .095, -.47, .5, 1)+
                   draw_grob(title.y,-0.217,0.0395,.5,1)+
                   draw_grob(labD,0.048,.470,1,1)+
                   draw_grob(labE,0.048,.153,1,1)+
                   draw_grob(labF,0.048,-0.167,1,1))
  
  Fig1b.print <- plot_grid(NULL,Fig1.print,nrow=2,rel_heights=c(0.015,1))
  
  Fig1c.print <- (Fig1b.print +
                    draw_grob(title,-0.15,0.4915,1,1)+       
                    draw_grob(title.S,-0.15,0.177,1,1)+     
                    draw_grob(title.A,-0.15,-0.1345,1,1)+
                    NULL)    
  
  print(Fig1c.print)
  dev.off()
  
  # Create raster scale gradient to use in legend:
  r <- raster(matrix(runif(100), ncol=10)) 
  
  # using raw time series:
  cuts=seq(log(0.0407622),log(9.974182),by=0.2)  # bounds of tramp function for P scenario (g C m-2 d-1 in log units)
  cutsA=seq(log(0.0407622),log(3.004166),by=0.2) # bounds of tramp function for A scenario (g C m-2 d-1 in log units)
  cutsS=seq(log(0.0407622),log(8.16617),by=0.2)  # bounds of tramp function for S scenario (g C m-2 d-1 in log units)
  
  par(bg = "transparent")  # switch off background to avoid obscuring adjacent plots 
  
  #Productive rivers:
  plot(r,fun=log, legend.only=TRUE, legend.width=1, legend.shrink=.65, 
       smallplot=c(0.1,0.20, 0.72,0.94),
       col=(tramp(seq(-3.2,2.3,by=0.2))),
       breaks=cuts, axis.args=list(at=c((-3.2),(-0.6931472),(0.6931472),(2.014903)), 
                                   labels=c("0.0",round(exp(-0.6931472),1),"2.0",round(exp(2.014903),1)), cex.axis=.9),
       legend.args=list(text=expression(""), side=3, font=2, 
                        line=0.5, cex=.9))  
  
  #Stochastic:
  plot(r,fun=log, legend.only=TRUE, legend.width=1, legend.shrink=.65, 
       smallplot=c(0.1,0.20, 0.4,0.64),
       col=(trampS(seq(-3.2,2.1,by=0.2))),
       breaks=cutsS, axis.args=list(at=c((-3.2),(-0.6931472),(0.6931472),(1.791759)), 
                                    labels=c("0.0",round(exp(-0.6931472),1),"2.0","6.0"), cex.axis=.9),
       legend.args=list(text=expression(""), side=3, font=2, 
                        line=0.5, cex=.9)) 
  
  #Unproductive rivers:
  plot(r,fun=log, legend.only=TRUE, legend.width=1, legend.shrink=.65, 
       smallplot=c(0.1,0.20, 0.09,0.32),
       col=(trampA(seq(-3.2,1.1,by=0.2))),
       breaks=cutsA, axis.args=list(at=c((-3.2),(-0.6931472),0,(0.6931472)), 
                                    labels=c("0.0",round(exp(-0.6931472),1),"1.0","2.0"), cex.axis=.9),
       legend.args=list(text=expression(""), side=3, font=2, 
                        line=0.5, cex=.9))  
  
  recordlegend <- recordPlot() # record the previous plot
  
  # Print Figure 1:
  Final.Fig1 <- plot_grid(Fig1c.print,recordlegend,rel_widths=c(1.65,.32))+
    draw_grob(legend.title, 0.656, 0.566, 0.5, 0.8)
  print(Final.Fig1)  #680width x 735height
  
  # Export plot:
  png("figures/Fig1_EnvelopeRegimes.png", width = 6.8, height = 7.32,units = "in",res = 300)
  print(Final.Fig1)
  dev.off()
  
  pdf("figures/Fig1_EnvelopeRegimes.pdf", width = 6.8, height = 7.32)
  print(Final.Fig1)
  dev.off()
  
  
##------------------------------------------------------------------------##
##   Figure 2: Network-scale productivity regimes across watershed size            
##------------------------------------------------------------------------##   
    
  # Plot network-scale GPP regimes for 3 modeled scenarios across a range in watershed size
  
  GPP40.Fig2.P   <- Subset.replicate.regimes(df=Sim40.P,chr.scenario="Productive",which.reps=1:length(b.40),b.size=b.40)
  GPP40.Fig2.S   <- Subset.replicate.regimes(df=Sim40.S,chr.scenario="Stochastic",which.reps=1:length(b.40),b.size=b.40)
  GPP40.Fig2.A   <- Subset.replicate.regimes(df=Sim40.A,chr.scenario="Unproductive",which.reps=1:length(b.40),b.size=b.40)
  
  GPP160.Fig2.P  <- Subset.replicate.regimes(df=Sim160.P,chr.scenario="Productive",which.reps=1:length(b.160),b.size=b.160)
  GPP160.Fig2.S  <- Subset.replicate.regimes(df=Sim160.S,chr.scenario="Stochastic",which.reps=1:length(b.160),b.size=b.160)
  GPP160.Fig2.A  <- Subset.replicate.regimes(df=Sim160.A,chr.scenario="Unproductive",which.reps=1:length(b.160),b.size=b.160)
  
  GPP450.Fig2.P  <- Subset.replicate.regimes(df=Sim450.P,chr.scenario="Productive",which.reps=1:length(b.450),b.size=b.450)
  GPP450.Fig2.S  <- Subset.replicate.regimes(df=Sim450.S,chr.scenario="Stochastic",which.reps=1:length(b.450),b.size=b.450)
  GPP450.Fig2.A  <- Subset.replicate.regimes(df=Sim450.A,chr.scenario="Unproductive",which.reps=1:length(b.450),b.size=b.450)
  
  GPP2600.Fig2.P <- Subset.replicate.regimes(df=Sim2600.P,chr.scenario="Productive",which.reps=1:length(b.2600),b.size=b.2600)
  GPP2600.Fig2.S <- Subset.replicate.regimes(df=Sim2600.S,chr.scenario="Stochastic",which.reps=1:length(b.2600),b.size=b.2600)
  GPP2600.Fig2.A <- Subset.replicate.regimes(df=Sim2600.A,chr.scenario="Unproductive",which.reps=1:length(b.2600),b.size=b.2600)
  
  # Define axis labels for 40 km^2 panels:
  LABELS.40 <- GPP40.Fig2.P %>% ggplot() + 
    geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname))+
    labs(x=expression("Day of year"),y=expression(River~network~GPP~(kg~C~d^-1)))+
    coord_cartesian(ylim=c(0,150))+
    scale_y_continuous(breaks=c(0,50,100,150))+
    scale_color_viridis(discrete=TRUE,direction = -1,option = "D")+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)), 
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  # Define axis labels for 160 km^2 panels:
  LABELS.160 <- GPP160.Fig2.P %>% ggplot() + 
    geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname))+
    labs(x=expression("Day of year"),y=expression(River~network~GPP~(kg~C~d^-1)))+
    coord_cartesian(ylim=c(0,600))+
    scale_y_continuous(breaks=c(0,200,400,600))+
    scale_color_viridis(discrete=TRUE,direction = -1,option = "D",begin=0,end=0.6)+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),       
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  # Define axis labels for 450 km^2 panels:
  LABELS.450 <- GPP450.Fig2.P %>% ggplot() + 
    geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname))+
    labs(x=expression("Day of year"),y=expression(River~network~GPP~(kg~C~d^-1)))+
    coord_cartesian(ylim=c(0,1800))+
    scale_y_continuous(breaks=c(0,600,1200,1800))+
    scale_color_viridis(discrete=TRUE,direction = -1,option = "D",begin=0,end=0.4)+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  # Define axis labels for 2,600 km^2 panels:
  LABELS.2600 <- GPP2600.Fig2.P %>% ggplot() + 
    geom_line(aes(x=doy,y=meanvalue*12/32,color=dfname))+
    labs(x=expression("Day of year"),y=expression(River~network~GPP~(kg~C~d^-1)))+
    coord_cartesian(ylim=c(0,12300))+
    scale_y_continuous(breaks=c(0,4000,8000,12000))+
    scale_color_viridis(discrete=TRUE,direction = 1,option = "D",begin = 0,end=.8)+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),       
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  # Create customized arrow to place in right margin:
  arrow.plot <- ggplot() + 
    geom_segment(aes(x=4, xend=4, y=11, yend=10),
                 arrow = arrow(length = unit(0.5, "cm")),size=1,color="gray30")+
    labs(x=expression('40'~km^2),y=expression(Increasing~watershed~size))
  arrow.plot2 <- ggplot() + 
    geom_segment(aes(x=4, xend=4, y=11, yend=10),
                 arrow = arrow(length = unit(0.5, "cm")),size=1,color="gray30")+
    labs(x=expression('2,600'~km^2),y=expression(Increasing~watershed~size))
  
  arrow <- get_arrow(arrow.plot)
  arrow.ws.small <- get_xaxis(arrow.plot)
  arrow.ws.large <- get_xaxis(arrow.plot2)
  arrow.incr <- get_yaxis(arrow.plot)
    
  # Generate individual grid plots for each watershed size:
  plot40.P   <- plot.grids40(GPP40.Fig2.P)
  plot40.S   <- plot.grids40(GPP40.Fig2.S)
  plot40.A   <- plot.grids40(GPP40.Fig2.A)
  
  plot160.P   <- plot.grids160(GPP160.Fig2.P)
  plot160.S   <- plot.grids160(GPP160.Fig2.S)
  plot160.A   <- plot.grids160(GPP160.Fig2.A)
  
  plot450.P   <- plot.grids450(GPP450.Fig2.P)
  plot450.S   <- plot.grids450(GPP450.Fig2.S)
  plot450.A   <- plot.grids450(GPP450.Fig2.A)
  
  plot2600.P   <- plot.grids2600(GPP2600.Fig2.P)
  plot2600.S   <- plot.grids2600(GPP2600.Fig2.S)
  plot2600.A   <- plot.grids2600(GPP2600.Fig2.A)
  
  # Extract titles/text needed for grid figure:
  doy_labels      <- get_xtext(LABELS.2600)
  GPP_labels.40   <- get_ytext(LABELS.40)
  GPP_labels.160  <- get_ytext(LABELS.160)
  GPP_labels.450  <- get_ytext(LABELS.450)
  GPP_labels.2500 <- get_ytext(LABELS.2600)
  title.y <- get_yaxis(LABELS.2600)
  title.x <- get_xaxis(LABELS.2600)
  
  title   <- ggdraw() + draw_label("Productive rivers", fontface='bold')
  title.S <- ggdraw() + draw_label("Stochastic", fontface='bold')
  title.A <- ggdraw() + draw_label("Unproductive rivers", fontface='bold')
  
  # Build multi-panel grid plot:
  Multi.plot.skeleton <- plot_grid(
    plot40.P,plot40.S,plot40.A,
    plot160.P,plot160.S,plot160.A,
    plot450.P,plot450.S,plot450.A,
    plot2600.P,plot2600.S,plot2600.A,
    ncol=3,align="hv")
  
  titles <- plot_grid(title,title.S,title.A,ncol=3)
  Multi.plot.titles <- plot_grid(titles,Multi.plot.skeleton, ncol=1, rel_heights=c(0.05, 1)) # rel_heights values control title margins
  
  Multi.plot.x <- plot_grid(Multi.plot.titles,NULL,NULL,ncol=1,rel_heights = c(1,0.05,0.04))
  Multi.plot.y <- plot_grid(NULL,Multi.plot.x,NULL,ncol=3,rel_widths = c(0.15,1,0.2))
  
  Multi.plot.print <- (Multi.plot.y + 
                         draw_grob(title.x,x=-0.013,y=-0.464,scale=1)+
                         draw_grob(title.y,x=-0.46,y=0.045,scale=1) + #y=0.0
                         draw_grob(GPP_labels.40,x=-0.468,y=0.351,scale=.185)+
                         draw_grob(GPP_labels.160,x=-0.468,y=0.1325,scale=.185)+
                         draw_grob(GPP_labels.450,x=-0.468,y=-0.085,scale=.185)+
                         draw_grob(GPP_labels.2500,x=-0.468,y=-0.304,scale=.185)+
                         draw_grob(doy_labels,x=.232,y=-0.513,scale=.218)+
                         draw_grob(doy_labels,x=-0.014,y=-0.513,scale=.218)+
                         draw_grob(doy_labels,x=-0.262,y=-0.513,scale=.218)+
                         draw_grob(arrow,x=.435,y=0.03,scale=.75)+
                         draw_grob(arrow.ws.small,x=.435,y=.40)+
                         draw_grob(arrow.ws.large,x=.435,y=-0.335)+
                         draw_grob(arrow.incr,x=.420,y=0.04))
  print(Multi.plot.print)
  
  # Export network grid (manuscript figure 2, 700x600 pixels):
  png("figures/Fig2_NetworkGPP_Grid.png", width = 8, height = 6.85,units = "in",res = 300)
  print(Multi.plot.print)
  dev.off()
  
  pdf("figures/Fig2_NetworkGPP_Grid.pdf", width=8, height=6.85)
  print(Multi.plot.print)
  dev.off()
  
  
##------------------------------------------------------------------##
##   Calculate role of small streams in network-scale productivity       
##------------------------------------------------------------------##  

  # Define 2,621 km^2 network:
  c.2600  <- get_ancestors(net2500,b.2600[1])
  subws <- induced.subgraph(net2500,c.2600)
  
  set.seed(1)
  
  stoch.2600 <- data.frame(V(subws)$name,V(subws)$width,V(subws)$ReachArea)
  colnames(stoch.2600) <- c('node','width','ReachArea') 
  stoch.2600$width.class <- assign.width.class(stoch.2600$width)
  
  ifelse(length(unique(stoch.2600$width.class))==1,
         block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
         ifelse(length(unique(stoch.2600$width.class))==2,
                block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                      c(0.0,0.18,0.46,0.36)),
                ifelse(length(unique(stoch.2600$width.class))==3,
                       block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                             c(0.0,0.18,0.46,0.36),
                                             c(0.25,0.09,0.33,0.33)),
                       block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                             c(0.0,0.18,0.46,0.36),
                                             c(0.25,0.09,0.33,0.33),
                                             c(0.67,0.25,0.08,0.0)))))
  
  stoch.2600$Regime <- block_ra(blocks = stoch.2600$width.class, 
                                block_prob_each = block_m_each,
                                conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
  
  # Calculate annual, network-scale GPP (kg O2 yr^-1):
  # creates matrix with dimensions 43 nodes (rows) x 1000 simulations (cols)
  Node.GPP.Prod   <- replicate(N,SmallStream.ntwFunctionP(network = subws))   
  Node.GPP.Unprod <- replicate(N,SmallStream.ntwFunctionA(network = subws))
  Node.GPP.Stoch  <- replicate(N,SmallStream.ntwFunctionS(network = stoch.2600))
  
  # Define levels of drainage area for which to calculate the proportion of total, network-scale GPP:
  areas <- as.numeric(c(0,20,50,100,150,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2621))
  
  # Calculate proportion of each area (defined above) to total network-scale GPP:
    # 1) Productive rivers:
    frac.GPP.P.matrix <- data.frame(matrix(NA, nrow = length(areas), ncol = 10))
    
    for(j in 1:length(areas)){
      which.nodes <- which(V(net2500)$DrainArea_km2<areas[j])
      for(i in 1:10){
        frac.GPP.P.size <- sum(Node.GPP.Prod[which.nodes,i])/sum(Node.GPP.Prod[,i])
        frac.GPP.P.matrix[j,i] <- frac.GPP.P.size 
      }
    }
    
    # 2) Unproductive rivers:
    frac.GPP.A.matrix <- data.frame(matrix(NA, nrow = length(areas), ncol = 10))
    
    for(j in 1:length(areas)){
      which.nodes <- which(V(net2500)$DrainArea_km2<areas[j])
      for(i in 1:10){
        frac.GPP.A.size <- sum(Node.GPP.Unprod[which.nodes,i])/sum(Node.GPP.Unprod[,i])
        frac.GPP.A.matrix[j,i] <- frac.GPP.A.size 
      }
    }
    
    # 3) Stochastic:
    frac.GPP.S.matrix <- data.frame(matrix(NA, nrow = length(areas), ncol = 10))
    
    for(j in 1:length(areas)){
      which.nodes <- which(V(net2500)$DrainArea_km2<areas[j])
      for(i in 1:10){
        frac.GPP.S.size <- sum(Node.GPP.Stoch[which.nodes,i])/sum(Node.GPP.Stoch[,i])
        frac.GPP.S.matrix[j,i] <- frac.GPP.S.size 
      }
    }
    
    Mean.P <- apply(frac.GPP.P.matrix,1,mean)
    Q025.P <- apply(frac.GPP.P.matrix,1,quantile,probs=0.025)
    Q975.P <- apply(frac.GPP.P.matrix,1,quantile,probs=0.975)
    Scenario.P <- "Productive"
    SmallStreams.P <- data.frame(Mean.P,Q025.P,Q975.P,Scenario.P)
    SmallStreams.P$Area <- areas
    colnames(SmallStreams.P) <- c("Mean","Q2.5","Q97.5","Scenario","WatershedSize")
    
    Mean.A <- apply(frac.GPP.A.matrix,1,mean)
    Q025.A <- apply(frac.GPP.A.matrix,1,quantile,probs=0.025)
    Q975.A <- apply(frac.GPP.A.matrix,1,quantile,probs=0.975)
    Scenario.A <- "Unproductive"
    SmallStreams.A <- data.frame(Mean.A,Q025.A,Q975.A,Scenario.A)
    SmallStreams.A$Area <- areas
    colnames(SmallStreams.A) <- c("Mean","Q2.5","Q97.5","Scenario","WatershedSize")
    
    Mean.S <- apply(frac.GPP.S.matrix,1,mean)
    Q025.S <- apply(frac.GPP.S.matrix,1,quantile,probs=0.025)
    Q975.S <- apply(frac.GPP.S.matrix,1,quantile,probs=0.975)
    Scenario.S <- "Stochastic"
    SmallStreams.S <- data.frame(Mean.S,Q025.S,Q975.S,Scenario.S)
    SmallStreams.S$Area <- areas
    colnames(SmallStreams.S) <- c("Mean","Q2.5","Q97.5","Scenario","WatershedSize")
    
  # Calculate proportion that each drainage area contributes to total streambed surface area in the river network:
  area.cumulative2600.0 <-    0
  area.cumulative2600.20 <-   sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<20)])/(1000*1000)
  area.cumulative2600.50 <-   sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<50)])/(1000*1000)
  area.cumulative2600.100 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<100)])/(1000*1000)
  area.cumulative2600.150 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<150)])/(1000*1000)
  area.cumulative2600.200 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<200)])/(1000*1000)
  area.cumulative2600.300 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<300)])/(1000*1000)
  area.cumulative2600.400 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<400)])/(1000*1000)
  area.cumulative2600.500 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<500)])/(1000*1000)
  area.cumulative2600.600 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<600)])/(1000*1000)
  area.cumulative2600.700 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<700)])/(1000*1000)
  area.cumulative2600.800 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<800)])/(1000*1000)
  area.cumulative2600.900 <-  sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<900)])/(1000*1000)
  area.cumulative2600.1000 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<1000)])/(1000*1000)
  area.cumulative2600.1200 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<1200)])/(1000*1000)
  area.cumulative2600.1400 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<1400)])/(1000*1000)
  area.cumulative2600.1600 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<1600)])/(1000*1000)
  area.cumulative2600.1800 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<1800)])/(1000*1000)
  area.cumulative2600.2000 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<2000)])/(1000*1000)
  area.cumulative2600.2200 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<2200)])/(1000*1000)
  area.cumulative2600.2400 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<2400)])/(1000*1000)
  area.cumulative2600.2600 <- sum(V(net2500)$ReachArea[which(V(net2500)$DrainArea_km2<2621.5)])/(1000*1000)
  
  SmallStreams.P$BenthicArea <- c(area.cumulative2600.0,area.cumulative2600.20,area.cumulative2600.50,area.cumulative2600.100,
                                  area.cumulative2600.150,area.cumulative2600.200,area.cumulative2600.300,area.cumulative2600.400,
                                  area.cumulative2600.500,area.cumulative2600.600,area.cumulative2600.700,area.cumulative2600.800,
                                  area.cumulative2600.900,area.cumulative2600.1000,area.cumulative2600.1200,area.cumulative2600.1400,
                                  area.cumulative2600.1600,area.cumulative2600.1800,area.cumulative2600.2000,area.cumulative2600.2200,
                                  area.cumulative2600.2400,area.cumulative2600.2600)
  SmallStreams.A$BenthicArea <- c(area.cumulative2600.0,area.cumulative2600.20,area.cumulative2600.50,area.cumulative2600.100,
                                  area.cumulative2600.150,area.cumulative2600.200,area.cumulative2600.300,area.cumulative2600.400,
                                  area.cumulative2600.500,area.cumulative2600.600,area.cumulative2600.700,area.cumulative2600.800,
                                  area.cumulative2600.900,area.cumulative2600.1000,area.cumulative2600.1200,area.cumulative2600.1400,
                                  area.cumulative2600.1600,area.cumulative2600.1800,area.cumulative2600.2000,area.cumulative2600.2200,
                                  area.cumulative2600.2400,area.cumulative2600.2600)
  SmallStreams.S$BenthicArea <- c(area.cumulative2600.0,area.cumulative2600.20,area.cumulative2600.50,area.cumulative2600.100,
                                  area.cumulative2600.150,area.cumulative2600.200,area.cumulative2600.300,area.cumulative2600.400,
                                  area.cumulative2600.500,area.cumulative2600.600,area.cumulative2600.700,area.cumulative2600.800,
                                  area.cumulative2600.900,area.cumulative2600.1000,area.cumulative2600.1200,area.cumulative2600.1400,
                                  area.cumulative2600.1600,area.cumulative2600.1800,area.cumulative2600.2000,area.cumulative2600.2200,
                                  area.cumulative2600.2400,area.cumulative2600.2600)
  
  # Create table:
  frac.datatable1 <- rbind(SmallStreams.P,SmallStreams.A,SmallStreams.S)
  frac.datatable1$Scenario <- factor(frac.datatable1$Scenario, levels = c("Productive", "Stochastic", "Unproductive"))
  
  # Re-format data table for plotting:
  frac.datatable <- data.frame(SmallStreams.P$WatershedSize,SmallStreams.P$BenthicArea,
                               SmallStreams.P$Mean,SmallStreams.P$Q2.5,SmallStreams.P$Q97.5,
                               SmallStreams.S$Mean,SmallStreams.S$Q2.5,SmallStreams.S$Q97.5,
                               SmallStreams.A$Mean,SmallStreams.A$Q2.5,SmallStreams.A$Q97.5)
  colnames(frac.datatable) <- c("WatershedSize","BenthicArea","Mean.P","Q_2.5_P","Q_97.5_P",
                                "Mean.S","Q_2.5_S","Q_97.5_S","Mean.A","Q_2.5_A","Q_97.5_A")
  
  StreamSize.Stats <- frac.datatable1 %>%
    filter(WatershedSize==100) %>%
    mutate(PercentGPP.contribution = (Mean*100),
           PercentBSA.contribution = (BenthicArea/9.010476))
  print(StreamSize.Stats)
  
  
##-----------------------------------------------------------------##
##   Figure 3: Contribution of small streams to network-scale GPP             
##-----------------------------------------------------------------## 
  
  frac.GPP <- frac.datatable %>% ggplot() + 
    geom_ribbon(aes(x=WatershedSize,ymin=Q_2.5_P,ymax=Q_97.5_P,fill="aProductive"),alpha=.35)+
    geom_line(aes(x=WatershedSize,y=Mean.P,color="aProductive"),size=1) +
    geom_ribbon(aes(x=WatershedSize,ymin=Q_2.5_A,ymax=Q_97.5_A,fill="cUnproductive"),alpha=.35)+
    geom_line(aes(x=WatershedSize,y=Mean.A,color="cUnproductive"),size=1) +
    geom_ribbon(aes(x=WatershedSize,ymin=Q_2.5_S,ymax=Q_97.5_S,fill="bStochastic"),alpha=.35)+
    geom_line(aes(x=WatershedSize,y=Mean.S,color="bStochastic"),size=1) +
    labs(x=expression(""),y=("Proportion of cumulative\nriver network GPP"))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+ 
    scale_x_continuous(breaks=c(0,500,1000,1500,2000,2500),labels=c(0,500,1000,1500,2000,2500))+
    scale_fill_manual(name="Model scenario",
                      breaks=c("aProductive","bStochastic","cUnproductive"),
                      labels=c("Productive rivers","Stochastic","Unproductive rivers"),
                      values=c("darkblue","#21908CFF","goldenrod1"))+
    scale_color_manual(name="Model scenario",
                       breaks=c("aProductive","bStochastic","cUnproductive"),
                       labels=c("Productive rivers","Stochastic","Unproductive rivers"),
                       values=c("darkblue","#21908CFF","goldenrod1"))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_blank(),
          legend.position="none")
  
  frac.GPP.lab <- frac.datatable1 %>% ggplot() + 
    geom_line(aes(x=WatershedSize,y=Mean,color=Scenario),size=1) +
    geom_ribbon(aes(x=WatershedSize,ymin=Q2.5,ymax=Q97.5,fill=Scenario),alpha=.35)+
    labs(x=expression(""),y=expression(Fraction~of~cumulative~river~network~GPP))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+ 
    scale_x_continuous(breaks=c(0,500,1000,1500,2000,2500),labels=c(0,500,1000,1500,2000,2500))+
    scale_color_manual(name="Model scenario",
                       breaks=c("Productive","Stochastic","Unproductive"),
                       labels=c("Productive rivers","Stochastic","Unproductive rivers"),
                       values=c("darkblue","#21908CFF","goldenrod1"))+
    scale_fill_manual(name="Model scenario",
                      breaks=c("Productive","Stochastic","Unproductive"),
                      labels=c("Productive rivers","Stochastic","Unproductive rivers"),
                      values=c("darkblue","#21908CFF","goldenrod1"))+  #purple:#440154FF
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_blank())
  
  cum.Area <- frac.datatable %>% ggplot() + geom_line(aes(x=WatershedSize,y=BenthicArea),color='gray30',size=1) + labs(x=expression(Contributing~watershed~area~(km^2)),y=expression(Streambed~surface~area~(km^2)))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 10))+ 
    scale_x_continuous(breaks=c(0,500,1000,1500,2000,2500),labels=c(0,500,1000,1500,2000,2500))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  cum.Area.lab <- frac.datatable %>% ggplot() + geom_line(aes(x=WatershedSize,y=BenthicArea),color='gray30',size=1) + labs(x=expression(Contributing~watershed~area~(km^2)),y=expression(Streambed~surface~area~(km^2)))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 10))+ 
    scale_x_continuous(breaks=c(0,500,1000,1500,2000,2500),labels=c(0,500,1000,1500,2000,2500))+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),     
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  plot1 <- plot_grid(frac.GPP,cum.Area,ncol=1,align="hv",labels=c('a)','b)'),label_fontface = "plain",label_x = 0.22,label_y = 0.98)
  plot1 <- plot_grid(frac.GPP,cum.Area,ncol=1,align="hv")
  plot1.xtitle <- get_xaxis(cum.Area.lab)
  plot1.x <- get_xtext(cum.Area.lab)
  legendA <- get_legend(frac.GPP.lab)
  plot2 <- plot_grid(plot1,NULL,ncol=1,rel_heights=c(1,0.055))
  plot3 <- plot_grid(plot2,NULL,ncol=2,rel_widths=c(1,0.06))
  plot4 <- (plot3 + draw_grob(plot1.xtitle,0.33,-0.465,0.5,1)+
              draw_grob(plot1.x,0.179,-0.915,0.74,1)+
              draw_grob(legendA,0.64,0.18,0.5,1))
  plot5 <- plot_grid(NULL,plot4,NULL,ncol=1,rel_heights=c(0.01,1,0.01))
  print(plot5)
  
  pdf("figures/Fig3_CumulativeGPP.pdf", width=5.2, height=7.8) ## png: 499 x 737
  print(plot5)
  dev.off()  
  
  png("figures/Fig3_CumulativeGPP.png", width = 5.2, height = 7.8,units = "in",res = 300)
  print((plot3 + draw_grob(plot1.xtitle,0.33,-0.465,0.5,1)+
           draw_grob(plot1.x,0.179,-0.915,0.74,1)+
           draw_grob(legendA,0.64,0.18,0.5,1)))
  dev.off()
  

##------------------------------------------------------------------##
##   Analyze light-availability model scenarios (objective 2)       
##------------------------------------------------------------------##
  
#--------- Vernal window scenarios ---------   
  
  # Modify empirical datasets to reflect extended vernal window (+7, +14 days):  
  springpeak.V7.dat <- springpeak.dat
  springpeak.V7.dat[,]<-NA
  for(i in 1:15){
    springpeak.V7.dat[c(0:39),i] <- springpeak.dat[c(0:39),i]
    springpeak.V7.dat[c(40:88),i] <- springpeak.dat[c(47:95),i]
    springpeak.V7.dat[c(96:365),i] <- springpeak.dat[c(96:365),i]
    springpeak.V7.dat[,i] <- na.approx(springpeak.V7.dat[,i])
  }
  
  springpeak.V14.dat <- springpeak.V7.dat
  springpeak.V14.dat[,]<-NA
  for(i in 1:15){
    springpeak.V14.dat[c(0:19),i] <- springpeak.V7.dat[c(0:19),i]
    springpeak.V14.dat[c(20:83),i] <- springpeak.V7.dat[c(26:89),i]
    springpeak.V14.dat[c(96:365),i] <- springpeak.V7.dat[c(96:365),i]
    springpeak.V14.dat[,i] <- na.approx(springpeak.V14.dat[,i])
  }
  
  # Define full network:
  c.2600 <- get_ancestors(net2500,b.2600[1])
  subws <- induced.subgraph(net2500,c.2600)
  
  stoch.RC <- data.frame(V(subws)$name,V(subws)$width,V(subws)$ReachArea)
  colnames(stoch.RC) <- c('node','width','ReachArea') 
  stoch.RC$width.class <- assign.width.class(stoch.RC$width)
  
  ifelse(length(unique(stoch.RC$width.class))==1,
         block_m_each <- rbind(c(.08,0.08,0.25,0.59)),
         ifelse(length(unique(stoch.RC$width.class))==2,
                block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                      c(0.0,0.18,0.46,0.36)),
                ifelse(length(unique(stoch.RC$width.class))==3,
                       block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                             c(0.0,0.18,0.46,0.36),
                                             c(0.25,0.09,0.33,0.33)),
                       block_m_each <- rbind(c(.08,0.08,0.25,0.59),
                                             c(0.0,0.18,0.46,0.36),
                                             c(0.25,0.09,0.33,0.33),
                                             c(0.67,0.25,0.08,0.0)))))
  
  stoch.RC$Regime <- block_ra(blocks = stoch.RC$width.class, 
                              block_prob_each = block_m_each,
                              conditions=c("SummerPeak", "Aseasonal", "SummerDecline","SpringPeak"))
  
  set.seed(1) 
  
  Network.Simulations.Prod.V7 <- replicate(N,Network.Prod.V7(network = subws))
  Network.Simulations.Unprod.V7 <- replicate(N,Network.Unprod.V7(network = subws))
  Network.Simulations.Stoch.V7 <- replicate(N,Network.Stoch.V7(network = stoch.RC))
  
  Network.Simulations.Prod.V14 <- replicate(N,Network.Prod.V14(network = subws))
  Network.Simulations.Unprod.V14 <- replicate(N,Network.Unprod.V14(network = subws))
  Network.Simulations.Stoch.V14 <- replicate(N,Network.Stoch.V14(network = stoch.RC))
  
  testV7.P <- apply(Network.Simulations.Prod.V7,2,sum)
  testV7.A <- apply(Network.Simulations.Unprod.V7,2,sum)
  testV7.S <- apply(Network.Simulations.Stoch.V7,2,sum)
  
  testV14.P <- apply(Network.Simulations.Prod.V14,2,sum)
  testV14.A <- apply(Network.Simulations.Unprod.V14,2,sum)
  testV14.S <- apply(Network.Simulations.Stoch.V14,2,sum)
  
  GPP.Simulations.all.vernal <- data.frame(testV7.P,testV7.S,testV7.A,testV14.P,testV14.S,testV14.A)
  vernal.mean <- apply(GPP.Simulations.all.vernal,2,mean)
  vernal.sd <- apply(GPP.Simulations.all.vernal,2,sd)
  
  GPP.Simulations.all.vernal <- data.frame(c(rep("vernal.7",3),rep("vernal.14",3)),rep(c("Productive","Stochastic","Unproductive"),2))
  colnames(GPP.Simulations.all.vernal) <- c("dfname","Scenario")
  GPP.Simulations.all.vernal$meanvalue <- vernal.mean*12/32
  GPP.Simulations.all.vernal$sd <- vernal.sd*12/32
  colnames(GPP.Simulations.all.vernal) <- c("dfname","Scenario","MeanAnnualGPP.kgCyr","sd.AnnualGPP.kgCyr")
  
  TableA.3.Vernal <- GPP.Simulations.all.vernal 
  print(TableA.3.Vernal)
  
  # Report mean annual GPP (kg C yr^-1) in Table A.3 as 3 significant digits:
  print(formatC(signif(TableA.3.Vernal$MeanAnnualGPP.kgCyr,digits=3), digits=3,format="fg", flag="#"))
  
  
  
#--------- Riparian clearing scenarios ---------   

  set.seed(1)
  
  subsample.20  <- sample(x=which(V(subws)$width<9),size=(0.2*(length(which(V(subws)$width<9)))))
  subsample.40  <- sample(x=which(V(subws)$width<9),size=(0.4*(length(which(V(subws)$width<9)))))
  subsample.60  <- sample(x=which(V(subws)$width<9),size=(0.6*(length(which(V(subws)$width<9)))))
  subsample.80  <- sample(x=which(V(subws)$width<9),size=(0.8*(length(which(V(subws)$width<9)))))
  subsample.100 <- sample(x=which(V(subws)$width<9),size=(1*(length(which(V(subws)$width<9)))))
  
  ## Define empty lists to house annual network GPP estimates for each RC scenario:
  Sim2600.P.RC = list()
  Sim2600.A.RC = list()
  Sim2600.S.RC = list()
  
  rc.sims <- c(0,20,40,60,80,100)
  
  for(i in 1:(length(rc.sims))){
    
    rc <- rc.sims[i]
    subsample.name <- paste("subsample",rc,sep=".")
    subsample.ntw <- V(subws)$name[sample(which(V(subws)$width<9),((rc/100)*(length(which(V(subws)$width<9)))))]
    
    set.seed(1)
    
    Network.Simulations.Prod.RC   <- replicate(N,Network.Prod.RC(network = subws))
    Network.Simulations.Unprod.RC <- replicate(N,Network.Unprod.RC(network = subws))
    Network.Simulations.Stoch.RC  <- replicate(N,Network.Stoch.RC(network = stoch.RC))
    
    Sim2600.P.RC[[i]] = Network.Simulations.Prod.RC
    Sim2600.A.RC[[i]] = Network.Simulations.Unprod.RC
    Sim2600.S.RC[[i]] = Network.Simulations.Stoch.RC
    
  }
  
  # Populate table that displays annual, network-scale GPP under the riparian clearing scenario
  P.2600.RC.table <- Fill.in.TableA3.allsims(df=Sim2600.P.RC,chr.scenario="Productive",size=2600,b.size=b.2600[1])
  S.2600.RC.table <- Fill.in.TableA3.allsims(df=Sim2600.S.RC,chr.scenario="Stochastic",size=2600,b.size=b.2600[1])
  A.2600.RC.table <- Fill.in.TableA3.allsims(df=Sim2600.A.RC,chr.scenario="Unproductive",size=2600,b.size=b.2600[1])
  
  TableA3.allsims <- bind_rows(P.2600.RC.table,S.2600.RC.table,A.2600.RC.table)
  
  TableA.3.RC <- TableA3.allsims %>%
    group_by(Scenario,dfname) %>%
    filter(dfname == "subsample.0"|dfname=="subsample.20"|dfname=="subsample.60"|dfname=="subsample.100") %>%
    summarize(MeanAnnualGPP.kgCyr = mean(AnnualGPP_kgC.yr),
              sd.AnnualGPP.kgCyr = sd(AnnualGPP_kgC.yr))
  print(TableA.3.RC)
  
  # Report mean annual GPP (kg C yr^-1) in Table A.3 as 3 significant digits:
  TableA.3.RC$MeanAnnualGPP_round <- formatC(signif(TableA.3.RC$MeanAnnualGPP.kgCyr,digits=3), digits=3,format="fg", flag="#")
  
  
##-------------------------------------------------##
##     Figure 4: Plot riparian clearing results                     
##-------------------------------------------------## 
  
  GPP2600.RC.Fig4.P   <- Subset.RC.regimes(df=Sim2600.P.RC,chr.scenario="Productive",b.size=b.2600[1])
  GPP2600.RC.Fig4.S   <- Subset.RC.regimes(df=Sim2600.S.RC,chr.scenario="Stochastic",b.size=b.2600[1])
  GPP2600.RC.Fig4.A   <- Subset.RC.regimes(df=Sim2600.A.RC,chr.scenario="Unproductive",b.size=b.2600[1])
  
  RC.levels <- c("subsample.0","subsample.20","subsample.40","subsample.60","subsample.80","subsample.100")
  
  mainpiece <- GPP2600.RC.Fig4.P %>% ggplot() +
    geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=factor(dfname,levels=RC.levels)),alpha=0.4)+
    geom_line(aes(x=doy,y=meanvalue*12/32,color=factor(dfname,levels=RC.levels)),size=0.7)+
    scale_color_manual(values=c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF"),guide=guide_legend(title= "% Riparian clearing"))+
    scale_fill_manual(values=c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF"),guide=guide_legend(title= "% Riparian clearing"))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,8000,16000,24000))+
    coord_cartesian(ylim=c(0,24000))+
    labs(x=expression(Day~of~year),y=expression(Network~GPP~(kg~C~d^-1)))+
    theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
          axis.text.y = element_blank(),axis.text.x = element_blank(),
          legend.position='none')
  
  mainpiece.A <- GPP2600.RC.Fig4.A %>% ggplot() +
    geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=factor(dfname,levels=RC.levels)),alpha=.4)+
    geom_line(aes(x=doy,y=meanvalue*12/32,color=factor(dfname,levels=RC.levels)),size=0.7)+
    scale_color_manual(values=c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF"),guide=guide_legend(title= "% Riparian clearing"))+
    scale_fill_manual(values=c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF"),guide=guide_legend(title= "% Riparian clearing"))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,8000,16000,24000))+
    coord_cartesian(ylim=c(0,24000))+
    labs(x=expression(Day~of~year),y=expression(Network~GPP~(kg~C~d^-1)))+
    theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
          axis.text.y = element_blank(),axis.text.x = element_blank(),
          legend.position='none')
  
  mainpiece.S <- GPP2600.RC.Fig4.S %>% ggplot() +
    geom_ribbon(aes(x=doy,ymin=Q2.5*12/32,ymax=Q97.5*12/32,fill=factor(dfname,levels=RC.levels)),alpha=.4)+
    geom_line(aes(x=doy,y=meanvalue*12/32,color=factor(dfname,levels=RC.levels)),size=0.7)+
    scale_color_manual(values=c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF"),guide=guide_legend(title= "% Riparian clearing"))+
    scale_fill_manual(values=c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF"),guide=guide_legend(title= "% Riparian clearing"))+
    scale_y_continuous(breaks=c(0,8000,16000,24000))+
    scale_x_continuous(breaks=c(0,100,200,300))+
    coord_cartesian(ylim=c(0,24000))+
    labs(x=expression(Day~of~year),y=expression(Network~GPP~(kg~C~d^-1)))+
    theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
          axis.text.y = element_blank(),axis.text.x = element_blank(),
          legend.position='none')
  
  GPP2600.RC.Fig4.S$Level <- as.integer(case_when(
    GPP2600.RC.Fig4.S$dfname=="subsample.0"    ~ 0,
    GPP2600.RC.Fig4.S$dfname=="subsample.20"   ~ 20,
    GPP2600.RC.Fig4.S$dfname=="subsample.40"   ~ 40,
    GPP2600.RC.Fig4.S$dfname=="subsample.60"   ~ 60,
    GPP2600.RC.Fig4.S$dfname=="subsample.80"   ~ 80,
    GPP2600.RC.Fig4.S$dfname=="subsample.100"  ~ 100))
  
  my_palette <- c("goldenrod1","#7AD151FF","#22A884FF","#2A788EFF","#414487FF","#440154FF")
  
  piece1 <- GPP2600.RC.Fig4.S %>% ggplot() + geom_line(aes(x=doy,y=meanvalue*12/32,color=Level),size=1)+
    scale_colour_gradientn(name = "", 
                           #low = "darkblue", high = "goldenrod1",
                           colors = (my_palette),                        
                           breaks=c(0,20,40,60,80,100),
                           labels=c("0", "20","40", "60","80","100"))+
    guides(color = guide_colorbar(barwidth = 1, barheight = 6),title.vjust=35000)+
    coord_cartesian(ylim=c(0,24000))+
    scale_y_continuous(breaks=c(0,8000,16000,24000))+
    labs(x="% small streams,\nriparian cleared  ")+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=14))
  
  piece2 <- GPP2600.RC.Fig4.S %>% ggplot() +  
    geom_line(aes(x=doy,y=meanvalue*12/32,color=Level),size=1)+
    scale_colour_gradient(name = "", 
                          low = "darkblue", high = "goldenrod",
                          breaks=c(0,20,40,60,80,100),
                          labels=c("0", "20","40", "60","80","100"))+
    guides(color = guide_colorbar(barwidth = 1, barheight = 6),title.vjust=35000)+
    coord_cartesian(ylim=c(0,24000))+
    scale_y_continuous(breaks=c(0,8000,16000,24000))+
    labs(x=expression(Day~of~year),y=expression(River~network~GPP~(kg~C~d^-1)))+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=14),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=14),           
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=13),
          axis.text.y=element_text(size=13),
          axis.ticks=element_blank(),
          legend.position="none")
  
  legend.landuse <- get_legend(piece1)
  legendtitle.landuse <- get_xaxis(piece1)
  axis.text.x <- get_xtext(piece2)
  axis.text.y <- get_ytext(piece2)
  axis.title.y <- get_yaxis(piece2)
  axis.title.x <- get_xaxis(piece2)
  
  title <- ggdraw() + draw_label("Productive rivers", fontface='bold')
  title.S <- ggdraw() + draw_label("Stochastic", fontface='bold')
  title.A <- ggdraw() + draw_label("Unproductive rivers", fontface='bold')
  
  titles <- plot_grid(title,title.S,title.A,ncol=3)
  
  plot.bare <- plot_grid(mainpiece,mainpiece.S,mainpiece.A,ncol=3,rel_widths=c(1,1,1))
  plot.titles <- plot_grid(titles,plot.bare,ncol=1,rel_heights=c(0.05,0.70))
  plot.spacey <- plot_grid(NULL,plot.titles,NULL,ncol=3,rel_widths=c(0.09,1,.13))
  plot.spacex <- plot_grid(NULL,plot.spacey,NULL,nrow=3,rel_heights=c(0.02,1,.15))
  
  Figtogether <- (plot.spacex + draw_grob(legend.landuse, x=0.79, y=0, scale=.8)+
                    draw_grob(legendtitle.landuse,x=0.43,y=0.25,scale=.8)+
                    draw_grob(axis.title.y,x=-0.48,y=0.048,scale=.3)+
                    draw_grob(axis.text.y,x=-0.780,y=0.033,scale=.736)+
                    draw_grob(axis.title.x,x=-0.008,y=-0.435,scale=.3)+
                    draw_grob(axis.text.x,x=-0.0142,y=-0.46,scale=.250)+
                    draw_grob(axis.text.x,x=-0.287,y=-0.46,scale=.250)+
                    draw_grob(axis.text.x,x=0.260,y=-0.46,scale=.250))
  print(Figtogether)
  
  png("figures/Figure4_LandUseScenarios.png",width = 10.5, height = 3.5,units="in", res=350)
  print(Figtogether)
  dev.off()
  
  pdf("figures/Figure4_LandUseScenarios.pdf",width = 10.5, height = 3.5)
  print(Figtogether)
  dev.off()
  

##--------------------------------------------------------------##
##  FIGURE S1: Plot optimal channel network and sub-catchments           
##--------------------------------------------------------------## 

  # Define sub-catchments in network file relative to igraph networks:
  list.40 <- get_ancestors(net2500,1554)    
  ex.40 <- (ocn512[ocn512$IDreach %in% list.40,])
  ex.40.ID <- ex.40$ID
  list.160 <- get_ancestors(net2500,1675)   
  ex.160 <- (ocn512[ocn512$IDreach %in% list.160,])
  ex.160.ID <- ex.160$ID
  list.450 <- get_ancestors(net2500,2377)   
  ex.450 <- (ocn512[ocn512$IDreach %in% list.450,])
  ex.450.ID <- ex.450$ID
  list.2500 <- get_ancestors(net2500,2834)  
  ex.2500 <- (ocn512[ocn512$IDreach %in% list.2500,])
  ex.2500.ID <- ex.2500$ID
    
  # Define colors for each sub-catchment:
  ocn512$ecol <- "gray10"
  ocn512$ecol[ocn512$ID %in% ex.40.ID] <- "#1CCC9B"       
  ocn512$ecol[ocn512$ID %in% ex.160.ID] <- "dodgerblue"   
  ocn512$ecol[ocn512$ID %in% ex.450.ID] <- "orange2"      
  #ocn512$ecol[ocn512$ID %in% ex.2500.ID] <- "gray25"
      
  # Create Figure A.1 in Network GPP MS:
  pdf(file="figures/FigureS1_subcatchments.pdf",width=5,height=5)
  par(mfrow=c(1,1), mar=c(0.1,0.45,0.1,0.1))
  plot(ocn512$X, ocn512$Y, pch=16, col=ocn512$ecol, cex=0.06*log(ocn512$Drainage_Area), xlab="", ylab='', axes=F) #cex originally set at 0.09
  dev.off()
    
        
##--------------------------------------------------##
##  FIGURE S2: Reproduce regimes from Savoy et al.        
##--------------------------------------------------##    

  # For reproducing median regimes from Savoy et al. (2019), use non-gap-filled data.(downloaded from CUAHSI Hydroshare: https://www.hydroshare.org/resource/eba152073b4046178d1a2ffe9a897ebe/).
  site.regimes.avg <- read.csv("data/Savoy_hydroshare/avg_gpp.csv",header=T,stringsAsFactors = FALSE)      # rows are sites and columns are day of year
    
  sites_summerdecline.avg <- which(colnames(site.regimes.avg) %in% (filter(SiteInfo.Table,four_clus=="summer decline")$Site_ID))
  sites_springpeak.avg <- which(colnames(site.regimes.avg) %in% (filter(SiteInfo.Table,four_clus=="spring peak")$Site_ID))
  sites_summerpeak.avg <- which(colnames(site.regimes.avg) %in% (filter(SiteInfo.Table,four_clus=="summer peak")$Site_ID))
  sites_aseasonal.avg <- which(colnames(site.regimes.avg) %in% (filter(SiteInfo.Table,four_clus=="aseasonal")$Site_ID))
    
  summerdecline.dat.avg <- site.regimes.avg[,sites_summerdecline]
  springpeak.dat.avg <- site.regimes.avg[,sites_springpeak]
  summerpeak.dat.avg <- site.regimes.avg[,sites_summerpeak]
  aseasonal.dat.avg <- site.regimes.avg[,sites_aseasonal]
    
  aseasonal.dat.plot <- create.smoothed.dat(data = aseasonal.dat.avg) 
  summerdecline.dat.plot <- create.smoothed.dat(data = summerdecline.dat.avg) 
  springpeak.dat.plot <- create.smoothed.dat(data = springpeak.dat.avg) 
  summerpeak.dat.plot <- create.smoothed.dat(data = summerpeak.dat.avg) 
    
  y.labels.aseasonal <- springpeak.dat.plot %>% 
    ggplot() + geom_line(aes(x=doy,y=medianGPP),color='darkgreen',size=1)+
    xlab("Day of year")+ylab(expression(GPP~(g~O[2]~m^-2~d^-1)))+ggtitle('Aseasonal')+
    coord_cartesian(ylim=c(0,10)) + 
    scale_x_continuous(breaks=c(0,100,200,300),labels=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,2,4,6,8,10),labels=c(0,2,4,6,8,10))+
    geom_ribbon(aes(x=doy,ymin=Quant25,ymax=Quant75),fill='green4',alpha=.4)
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin=margin(t=20,r=0,b=0,l=0)))
    
  labels <- springpeak.dat.plot %>% 
    ggplot() + geom_line(aes(x=doy,y=medianGPP),color='darkgreen',size=1)+
    xlab("Day of year")+ylab(expression(GPP~(g~O[2]~m^-2~d^-1)))+ggtitle('Aseasonal')+
    coord_cartesian(ylim=c(0,10)) +
    scale_x_continuous(breaks=c(0,100,200,300),labels=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,2,4,6,8,10),labels=c(0,2,4,6,8,10))+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),        
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
    
  main.plot <- plot_grid(
    plot.regime(dat = aseasonal.dat.plot,
                regime.title.chr = "Aseasonal"),
    plot.regime(dat = summerpeak.dat.plot,
                regime.title.chr = "Summer peak"),
    plot.regime(dat = summerdecline.dat.plot,
                regime.title.chr = "Summer decline"),
    plot.regime(dat = springpeak.dat.plot,
                regime.title.chr = "Spring peak"),
    ncol=2,align='hv',labels=c("a)","b)","c)","d)"),label_fontface = "plain", 
    label_x = 0.039,label_y = 0.98)
    
  main.xtitle <- get_xaxis(y.labels.aseasonal)
  main.ytitle <- get_yaxis(y.labels.aseasonal)
  main.xlabels <- get_xtext(labels)
  main.ylabels <- get_ytext(labels)
    
  main.plot <- plot_grid(NULL,main.plot,ncol=2,rel_widths=c(0.13,1))
  main.plot <- plot_grid(main.plot,NULL,ncol=1,rel_heights = c(1,0.13))
    
  main.1b <- (main.plot + 
                draw_grob(main.xtitle,x=0.070,y=-0.43,scale=0.2)+
                draw_grob(main.ytitle,x=-0.45,y=0.1,scale=0.2)+
                draw_grob(main.xlabels,x=-0.16,y=-0.558,scale=0.393)+
                draw_grob(main.xlabels,x=0.283,y=-0.558,scale=0.393)+
                draw_grob(main.ylabels,x=-0.53,y=0.256,scale=0.335)+
                draw_grob(main.ylabels,x=-0.53,y=-0.185,scale=0.335))
  #print(main.1b) #512x443
    
  png("figures/FigureS2_SavoyRegimes.png",width = 5.5, height = 4.761,units="in", res=350)
  print(main.1b)
  dev.off()
    
  pdf("figures/FigureS2_SavoyRegimes.pdf",width = 5.5, height = 4.761)
  print(main.1b)
  dev.off()

  
##--------------------------------------------------##
##  FIGURE S3: Plot reach-scale regimes (median)
##       light-availability model scenarios
##--------------------------------------------------##

  # Recproduce vernal window simulations from median productivity regime lines 
  # the median smoothed data were used only for the purposes of plotting the vernal window regimes as in Figure S3 of the manuscript text
  # +7 days:
  springpeak.dat.plot$medianGPP_smooth_V7 <- NA
  springpeak.dat.plot$medianGPP_smooth_V7[c(0:39)] <- springpeak.dat.plot$medianGPP_smooth[c(0:39)]
  springpeak.dat.plot$medianGPP_smooth_V7[c(40:88)] <- springpeak.dat.plot$medianGPP_smooth[c(47:95)]
  springpeak.dat.plot$medianGPP_smooth_V7[c(96:365)] <- springpeak.dat.plot$medianGPP_smooth[c(96:365)]
  springpeak.dat.plot$medianGPP_smooth_V7 <- na.approx(springpeak.dat.plot$medianGPP_smooth_V7)
  
  springpeak.dat.plot$Quant25_smooth_V7 <- NA
  springpeak.dat.plot$Quant25_smooth_V7[c(0:39)] <- springpeak.dat.plot$Quant25_smooth[c(0:39)]
  springpeak.dat.plot$Quant25_smooth_V7[c(40:88)] <- springpeak.dat.plot$Quant25_smooth[c(47:95)]
  springpeak.dat.plot$Quant25_smooth_V7[c(96:365)] <- springpeak.dat.plot$Quant25_smooth[c(96:365)]
  springpeak.dat.plot$Quant25_smooth_V7 <- na.approx(springpeak.dat.plot$Quant25_smooth_V7)
  
  springpeak.dat.plot$Quant75_smooth_V7 <- NA
  springpeak.dat.plot$Quant75_smooth_V7[c(0:39)] <- springpeak.dat.plot$Quant75_smooth[c(0:39)]
  springpeak.dat.plot$Quant75_smooth_V7[c(40:88)] <- springpeak.dat.plot$Quant75_smooth[c(47:95)]
  springpeak.dat.plot$Quant75_smooth_V7[c(96:365)] <- springpeak.dat.plot$Quant75_smooth[c(96:365)]
  springpeak.dat.plot$Quant75_smooth_V7 <- na.approx(springpeak.dat.plot$Quant75_smooth_V7)
  
  # +14 days
  springpeak.dat.plot$medianGPP_smooth_V14 <- NA
  springpeak.dat.plot$medianGPP_smooth_V14[c(0:19)] <- springpeak.dat.plot$medianGPP_smooth_V7[c(0:19)]
  springpeak.dat.plot$medianGPP_smooth_V14[c(20:83)] <- springpeak.dat.plot$medianGPP_smooth_V7[c(26:89)]
  springpeak.dat.plot$medianGPP_smooth_V14[c(96:365)] <- springpeak.dat.plot$medianGPP_smooth_V7[c(96:365)]
  springpeak.dat.plot$medianGPP_smooth_V14 <- na.approx(springpeak.dat.plot$medianGPP_smooth_V14)
  
  springpeak.dat.plot$Quant25_smooth_V14 <- NA
  springpeak.dat.plot$Quant25_smooth_V14[c(0:19)] <- springpeak.dat.plot$Quant25_smooth_V7[c(0:19)]
  springpeak.dat.plot$Quant25_smooth_V14[c(20:83)] <- springpeak.dat.plot$Quant25_smooth_V7[c(26:89)]
  springpeak.dat.plot$Quant25_smooth_V14[c(96:365)] <- springpeak.dat.plot$Quant25_smooth_V7[c(96:365)]
  springpeak.dat.plot$Quant25_smooth_V14 <- na.approx(springpeak.dat.plot$Quant25_smooth_V14)
  
  springpeak.dat.plot$Quant75_smooth_V14 <- NA
  springpeak.dat.plot$Quant75_smooth_V14[c(0:19)] <- springpeak.dat.plot$Quant75_smooth_V7[c(0:19)]
  springpeak.dat.plot$Quant75_smooth_V14[c(20:83)] <- springpeak.dat.plot$Quant75_smooth_V7[c(26:89)]
  springpeak.dat.plot$Quant75_smooth_V14[c(96:365)] <- springpeak.dat.plot$Quant75_smooth_V7[c(96:365)]
  springpeak.dat.plot$Quant75_smooth_V14 <- na.approx(springpeak.dat.plot$Quant75_smooth_V14)
  
  # Create plot:
  pdf("figures/FigureS3_ClimateLandUseRegimes.pdf",width = 5.5, height = 7.8)
  FigA <- ggplot() + xlab("")+ylab("")+
    coord_cartesian(ylim=c(0,2)) + 
    scale_x_continuous(breaks=c(0,100,200,300),labels=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,0.5,1,1.5,2))+
    geom_ribbon(data = springpeak.dat.plot, aes(x=doy,ymin=Quant25_smooth_V14*12/32,ymax=Quant75_smooth_V14*12/32,fill="e14 days advanced"),alpha=1)+  
    geom_ribbon(data = springpeak.dat.plot, aes(x=doy,ymin=Quant25_smooth_V7*12/32,ymax=Quant25_smooth_V7*12/32,fill="d7 days advanced"),alpha=1)+  
    geom_ribbon(data = springpeak.dat.plot, aes(x=doy,ymin=Quant25_smooth*12/32,ymax=Quant75_smooth*12/32,fill="baseline"),alpha=1)+
    geom_blank(data = summerpeak.dat.plot, aes(x=doy,ymin=Quant25_smooth*12/32,ymax=Quant75_smooth*12/32,fill="friparian"))+
    geom_blank(data = summerpeak.dat.plot, aes(x=doy,y=medianGPP_smooth*12/32,color="friparian"))+
    geom_line(data = springpeak.dat.plot, aes(x=doy,y=medianGPP_smooth_V14*12/32,color="e14 days advanced"),size=1)+
    geom_line(data = springpeak.dat.plot, aes(x=doy,y=medianGPP_smooth_V7*12/32,color="d7 days advanced"),size=1)+
    geom_line(data = springpeak.dat.plot, aes(x=doy,y=medianGPP_smooth*12/32,color="baseline"),size=1)+
    scale_color_manual(name= "Model scenario",
                       breaks=c("baseline","d7 days advanced","e14 days advanced","friparian"),
                       labels=c(" Base scenario"," 7 days advanced"," 14 days advanced"," Riparian clearing"),                                                         
                       values=c("darkblue","#ece7f2","#2b8cbe","forestgreen"))+  
    scale_fill_manual(name="Model scenario",
                      breaks=c("baseline","d7 days advanced","e14 days advanced","friparian"),
                      labels=c(" Base scenario"," 7 days advanced"," 14 days advanced"," Riparian clearing"),
                      values=c("#8E9BC6","gray80","#A0B6C9","#D5E7D5"))+  
    theme(plot.margin = margin(t=10,r=15,b=2,l=10),
          legend.position=c(.6,.8),
          legend.title = element_blank(),
          legend.key.size = unit(1.1,"line"),
          axis.text.x=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin=margin(t=0,r=0,b=0,l=0)))
  
  FigLABEL <- ggplot()+
    geom_blank(data = springpeak.dat.plot, aes(x=doy,y=medianGPP_smooth*12/32,color="friparian"))+
    xlab("")+
    ylab(expression(Reach-scale~GPP~(g~C~m^-2~d^-1)))+
    coord_cartesian(ylim=c(0,2))+
    scale_x_continuous(breaks=c(0,100,200,300),labels=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,0.5,1,1.5,2))+
    theme(panel.background=element_blank(),
          panel.border=element_rect(color="transparent"),
          plot.margin = unit(c(.2, 0, .2, .2), "cm"),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),        
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.ticks=element_blank(),
          legend.position="none")
  
  FigB <- ggplot() + 
    geom_ribbon(data = springpeak.dat.plot, aes(x=doy,ymin=Quant25_smooth*12/32,ymax=Quant75_smooth*12/32),fill="#8E9BC6",alpha=1)+
    geom_ribbon(data = summerpeak.dat.plot, aes(x=doy,ymin=Quant25_smooth*12/32,ymax=Quant75_smooth*12/32),fill="forestgreen",alpha=.23)+
    guides(fill=FALSE)+
    geom_line(data = summerpeak.dat.plot, aes(x=doy,y=medianGPP_smooth*12/32,color='Riparian clearing'),size=1)+
    geom_line(data = springpeak.dat.plot, aes(x=doy,y=medianGPP_smooth*12/32,color="Base scenario"),size=1)+
    scale_color_manual(values=c("darkblue", "forestgreen"))+
    xlab("Day of year")+ylab("")+
    coord_cartesian(ylim=c(0,3.5))+
    scale_x_continuous(breaks=c(0,100,200,300),labels=c(0,100,200,300))+
    scale_y_continuous(breaks=c(0,1,2,3),labels=c("0.0","1.0","2.0","3.0"))+
    theme(plot.margin = margin(t=2,r=15,b=15,l=10),
          legend.position="none",
          legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin=margin(t=15,r=0,b=0,l=0)))
  
  # Join multiple panels in one plot:
  together <- plot_grid(FigA,FigB,labels=c("a)","b)"),label_fontface = "plain",label_x=0.17,label_y=.98,ncol=1,align="v",rel_heights=c(.855,1))
  together2 <- plot_grid(NULL,together,ncol=2,rel_widths=c(0.03,1))
  
  title.y <- get_yaxis(FigLABEL)
  print(together2 + draw_grob(title.y, -.21, .05, .5, 1))
  
  dev.off()
  
  
##-----------------------------------------------------------##
##     Figure S4: Allometric scaling of river network GPP                
##-----------------------------------------------------------##
  
  Allometric <- Table1 %>% ggplot() + 
    geom_smooth(aes(x=WatershedSize,y=MeanAnnualGPP.kgCyr,color=Scenario),size=1.0,method="lm",se = FALSE)+
    geom_point(aes(x=WatershedSize,y=MeanAnnualGPP.kgCyr,color=Scenario),size=3)+
    scale_x_log10()+
    coord_cartesian(xlim=c(34,3000))+
    scale_y_log10(breaks=c(10000,100000,1000000),
                  labels=c("10,000","100,000","1,000,000"))+
    xlab(expression(Watershed~size~(km^2)))+ylab(expression(Mean~Annual~GPP~(kg~C~yr^-1)))+
    scale_color_manual(name="Model scenario",
                       breaks=c("Productive","Stochastic","Unproductive"),
                       labels=c("Productive rivers","Stochastic","Unproductive rivers"),
                       values=c("darkblue","#21908CFF","goldenrod1"))
  
  pdf("figures/FigureS4_AllometricScaling.pdf",width = 6.048, height = 3.75)
  print(Allometric)
  dev.off() ## size: 605 x 375
  
  # Print allometric scaling slopes
  Allometric.scaling.results <- Table1 %>% group_by(Scenario) %>%
    do(model = lm(log10(MeanAnnualGPP.kgCyr) ~ log10(WatershedSize), data = .))
  Allometric.scaling.Coef = broom::tidy(Allometric.scaling.results, model)
  print(Allometric.scaling.Coef)
  
  
##--------------------------------------------------------##
##     Table S1: Physical properties of sub-catchments                
##--------------------------------------------------------##  

  TableA1.40_allreps <- Fill.in.TableA.1(40,b.40)
  TableA1.40 <- TableA1.40_allreps %>% 
    mutate(DrainageDensity = ((Sum_length_m/1000)/DrainageArea_km2)) %>%
    group_by(WatershedSize) %>%
    summarize(MeanDrainDensity = mean(DrainageDensity),
              SD.DrainDensity = sd(DrainageDensity),
              Mean.BSA.km2 = mean(ReachArea_km2),
              SD.BSA.km2 = sd(ReachArea_km2),
              Mean.reach.width.m = mean(Mean_width_m),
              SD.reach.width.m = sd(Mean_width_m))
  
  TableA1.160_allreps <- Fill.in.TableA.1(160,b.160)
  TableA1.160 <- TableA1.160_allreps %>% 
    mutate(DrainageDensity = ((Sum_length_m/1000)/DrainageArea_km2)) %>%
    group_by(WatershedSize) %>%
    summarize(MeanDrainDensity = mean(DrainageDensity),
              SD.DrainDensity = sd(DrainageDensity),
              Mean.BSA.km2 = mean(ReachArea_km2),
              SD.BSA.km2 = sd(ReachArea_km2),
              Mean.reach.width.m = mean(Mean_width_m),
              SD.reach.width.m = sd(Mean_width_m))
  
  TableA1.450_allreps <- Fill.in.TableA.1(450,b.450)
  TableA1.450 <- TableA1.450_allreps %>% 
    mutate(DrainageDensity = ((Sum_length_m/1000)/DrainageArea_km2)) %>%
    group_by(WatershedSize) %>%
    summarize(MeanDrainDensity = mean(DrainageDensity),
              SD.DrainDensity = sd(DrainageDensity),
              Mean.BSA.km2 = mean(ReachArea_km2),
              SD.BSA.km2 = sd(ReachArea_km2),
              Mean.reach.width.m = mean(Mean_width_m),
              SD.reach.width.m = sd(Mean_width_m))
  
  TableA1.2600_allreps <- Fill.in.TableA.1(2600,b.2600)
  TableA1.2600 <- TableA1.2600_allreps %>% 
    mutate(DrainageDensity = ((Sum_length_m/1000)/DrainageArea_km2)) %>%
    group_by(WatershedSize) %>%
    summarize(MeanDrainDensity = mean(DrainageDensity),
              SD.DrainDensity = sd(DrainageDensity),
              Mean.BSA.km2 = mean(ReachArea_km2),
              SD.BSA.km2 = sd(ReachArea_km2),
              Mean.reach.width.m = mean(Mean_width_m),
              SD.reach.width.m = sd(Mean_width_m))
  
  TableA1 <- bind_rows(TableA1.40,TableA1.160,TableA1.450,TableA1.2600)
  print(TableA1)
  
  