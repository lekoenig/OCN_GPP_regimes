
################################################################
#####   FUNCTIONS FOR USE IN ANALYSIS: River Network GPP   ##### 
################################################################


##-----------------------------------------------------##
##     Create sub-catchments within an igraph object      
##-----------------------------------------------------##

  get_ancestors <- function(g,v) {
    return(igraph::subcomponent(g,v,"in"))
  }


  
##----------------------------------------------------------##
##   Assign GPP to individual reaches within river network    
##    under three different modeled scenarios (g O2 d^-1)
##----------------------------------------------------------##

  # Productive rivers:
    calcGPP.Prod <- function(width,area){
     GPP.doy <- ifelse(width>9,sample(summerpeak.dat,size=1,replace=TRUE)*area,sample(springpeak.dat,size=1,replace=TRUE)*area)
    }

  # Unproductive rivers:
    calcGPP.Unprod <- function(width,area){
     GPP.doy <- ifelse(width>9,sample(aseasonal.dat,size=1,replace=TRUE)*area,sample(springpeak.dat,size=1,replace=TRUE)*area)
    }  

  # Stochastic:
    calcGPP.Stoch <- function(Regime,area){  
     GPP.doy <- ifelse(Regime == "SummerPeak",sample(summerpeak.dat,size=1,replace=TRUE)*area,
                        ifelse(Regime=="Aseasonal",sample(aseasonal.dat,size=1,replace=TRUE)*area,
                               ifelse(Regime=="SummerDecline",sample(summerdecline.dat,size=1,replace=TRUE)*area,
                                      sample(springpeak.dat,size=1,replace=TRUE)*area)))
    }

  # Assign discrete width bin to individual reaches within the river network: 
    assign.width.class <- function(width){
     width.class <- as.factor(case_when(
        width < 5.31    ~ "width.1",
        width < 12.01   ~ "width.2",
        width < 32.41   ~ "width.3",
        width > 32.41   ~ "width.4"))
    }

    
    
##----------------------------------------------------------##
##   Assign GPP to individual reaches within river network    
##          riparian clearing scenario (g O2 d^-1)
##----------------------------------------------------------##

  # Productive rivers:
    calcGPP.Prod.RC <- function(width,area,name){
      GPP.doy <-  ifelse(name %in% subsample.ntw,sample(summerpeak.dat,size=1,replace=TRUE)*area,
                         ifelse(width>9,sample(summerpeak.dat,size=1,replace=TRUE)*area,
                                sample(springpeak.dat,size=1,replace=TRUE)*area))
    }

  # Unproductive rivers:
    calcGPP.Unprod.RC <- function(width,area,name){
      GPP.doy <-  ifelse(name %in% subsample.ntw,sample(summerpeak.dat,size=1,replace=TRUE)*area,
                         ifelse(width>9,sample(aseasonal.dat,size=1,replace=TRUE)*area,
                                sample(springpeak.dat,size=1,replace=TRUE)*area))
    }

  # Stochastic:
    calcGPP.Stoch.RC <- function(Regime,area,name){  
      GPP.doy <- ifelse(name %in% subsample.ntw,sample(summerpeak.dat,size=1,replace=TRUE)*area,
                        ifelse(Regime == "SummerPeak",sample(summerpeak.dat,size=1,replace=TRUE)*area,
                               ifelse(Regime=="Aseasonal",sample(aseasonal.dat,size=1,replace=TRUE)*area,
                                      ifelse(Regime=="SummerDecline",sample(summerdecline.dat,size=1,replace=TRUE)*area,sample(springpeak.dat,size=1,replace=TRUE)*area))))
    }

    
    
##----------------------------------------------------------##
##   Assign GPP to individual reaches within river network    
##     vernal window scenarios: +7, +14 days (g O2 d^-1)
##----------------------------------------------------------##    

  # Productive
    calcGPP.Prod.V7 <- function(width,area,name){
      GPP.doy <-  ifelse(width>9,sample(summerpeak.dat,size=1,replace=TRUE)*area,
                         sample(springpeak.V7.dat,size=1,replace=TRUE)*area)
    }
    calcGPP.Prod.V14 <- function(width,area,name){
      GPP.doy <-  ifelse(width>9,sample(summerpeak.dat,size=1,replace=TRUE)*area,
                         sample(springpeak.V14.dat,size=1,replace=TRUE)*area)
    }

  # Unproductive:
    calcGPP.Unprod.V7 <- function(width,area,name){
      GPP.doy <-  ifelse(width>9,sample(aseasonal.dat,size=1,replace=TRUE)*area,
                         sample(springpeak.V7.dat,size=1,replace=TRUE)*area)
    }
    calcGPP.Unprod.V14 <- function(width,area,name){
      GPP.doy <-  ifelse(width>9,sample(aseasonal.dat,size=1,replace=TRUE)*area,
                         sample(springpeak.V14.dat,size=1,replace=TRUE)*area)
    }

  # Stochastic:
    calcGPP.Stoch.V7 <- function(Regime,area){  
      GPP.doy <- ifelse(Regime == "SummerPeak",sample(summerpeak.dat,size=1,replace=TRUE)*area,
                        ifelse(Regime=="Aseasonal",sample(aseasonal.dat,size=1,replace=TRUE)*area,
                               ifelse(Regime=="SummerDecline",sample(summerdecline.dat,size=1,replace=TRUE)*area,sample(springpeak.V7.dat,size=1,replace=TRUE)*area)))
    }
    calcGPP.Stoch.V14 <- function(Regime,area){  
      GPP.doy <- ifelse(Regime == "SummerPeak",sample(summerpeak.dat,size=1,replace=TRUE)*area,
                        ifelse(Regime=="Aseasonal",sample(aseasonal.dat,size=1,replace=TRUE)*area,
                               ifelse(Regime=="SummerDecline",sample(summerdecline.dat,size=1,replace=TRUE)*area,sample(springpeak.V14.dat,size=1,replace=TRUE)*area)))
    }



##-------------------------------------------------##
##     Calculate network-scale GPP (kg O2 d^-1)     
##-------------------------------------------------##  

  # Productive:
    Network.Prod <- function(network){
      Network.GPP <- mapply(calcGPP.Prod,V(network)$width,V(network)$ReachArea) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

  # Unproductive:
    Network.Unprod <- function(network){
      Network.GPP <- mapply(calcGPP.Unprod,V(network)$width,V(network)$ReachArea) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

  # Stochastic:
    Network.Stoch <- function(network){
                Network.GPP <- mapply(calcGPP.Stoch,network$Regime,network$ReachArea) %>%
                               transpose() %>%
                               lapply(flatten_dbl) %>%
                               map(sum) %>%
                               flatten_dbl()/1000 
    }



##-------------------------------------------------##
##     Calculate network-scale GPP (kg O2 d^-1) 
##            riparian clearing scenario
##-------------------------------------------------##  

  # Productive:
    Network.Prod.RC <- function(network){
      Network.GPP <- mapply(calcGPP.Prod.RC,V(network)$width,V(network)$ReachArea,V(network)$name) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

  # Unproductive:
    Network.Unprod.RC <- function(network){
      Network.GPP <- mapply(calcGPP.Unprod.RC,V(network)$width,V(network)$ReachArea,V(network)$name) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

  # Stochastic:
    Network.Stoch.RC <- function(network){
      Network.GPP <- mapply(calcGPP.Stoch.RC,network$Regime,network$ReachArea,network$node) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

    
##-------------------------------------------------##
##     Calculate network-scale GPP (kg O2 d^-1) 
##       vernal window scenarios: +7, +14 days
##-------------------------------------------------## 
    
  # Productive:
    Network.Prod.V7 <- function(network){
      Network.GPP <- mapply(calcGPP.Prod.V7,V(network)$width,V(network)$ReachArea,V(network)$name) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }
    Network.Prod.V14 <- function(network){
      Network.GPP <- mapply(calcGPP.Prod.V14,V(network)$width,V(network)$ReachArea,V(network)$name) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

  # Unproductive:
    Network.Unprod.V7 <- function(network){
      Network.GPP <- mapply(calcGPP.Unprod.V7,V(network)$width,V(network)$ReachArea,V(network)$name) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }
    Network.Unprod.V14 <- function(network){
      Network.GPP <- mapply(calcGPP.Unprod.V14,V(network)$width,V(network)$ReachArea,V(network)$name) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }

  # Stochastic:
    Network.Stoch.V7 <- function(network){
      Network.GPP <- mapply(calcGPP.Stoch.V7,network$Regime,network$ReachArea) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }
    Network.Stoch.V14 <- function(network){
      Network.GPP <- mapply(calcGPP.Stoch.V14,network$Regime,network$ReachArea) %>%
        transpose() %>%
        lapply(flatten_dbl) %>%
        map(sum) %>%
        flatten_dbl()/1000 
    }



##---------------------------------------------------------##
##   Calculate role of small streams to network-scale GPP    
##       (proportion of total 2,621 km^2 network)
##                  used in figure 3
##---------------------------------------------------------##

  # Productive:
    SmallStream.ntwFunctionP  <- function(network){
      Network.GPP <- mapply(calcGPP.Prod,V(network)$width,V(network)$ReachArea) %>%   # daily rates of GPP for each node (g O2 d^-1)
        lapply(sum) %>%            # daily rates of GPP for each node (g O2 yr^-1)
        flatten_dbl()/1000         # annual rates of GPP for each node (kg O2 yr^-1)
    }

  # Unproductive:
    SmallStream.ntwFunctionA  <- function(network){
      Network.GPP <- mapply(calcGPP.Unprod,V(network)$width,V(network)$ReachArea) %>%   # daily rates of GPP for each node (g O2 d^-1)
        lapply(sum) %>%            # daily rates of GPP for each node (g O2 yr^-1)
        flatten_dbl()/1000         # annual rates of GPP for each node (kg O2 yr^-1)
    }

  # Stochastic:
    SmallStream.ntwFunctionS  <- function(network){
      Network.GPP <- mapply(calcGPP.Stoch,network$Regime,network$ReachArea) %>%   # daily rates of GPP for each node (g O2 d^-1)
        lapply(sum) %>%            # daily rates of GPP for each node (g O2 yr^-1)
        flatten_dbl()/1000         # annual rates of GPP for each node (kg O2 yr^-1)
    }


