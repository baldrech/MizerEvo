# Plotting methods for the MizerSim class
library(gridBase) # to get grid.newpage
# Biomass through time
#' Plot the biomass of each species through time
#'
#' After running a projection, the biomass of each species can be plotted against time. The biomass is calculated within user defined size limits (see \code{\link{getBiomass}}).
#' This plot is pretty easy to do by hand. It just gets the biomass using the \code{\link{getBiomass}} method and plots using the ggplot2 package. You can then fiddle about with colours and linetypes etc. Just look at the source code for details.

plotBiomass <- function(object, print_it=TRUE, start_time=as.numeric(dimnames(object@n)[[1]][1]), end_time = as.numeric(dimnames(object@n)[[1]][dim(object@n)[1]]), ...){
  b <- getBiomass(object, ...)
  names(dimnames(b))[names(dimnames(b))=="sp"] <- "Species"
  if(start_time >= end_time){
    stop("start_time must be less than end_time")
  }
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <= end_time),,drop=FALSE]
  bm <- melt(b)
  # Force Species column to be a character (if numbers used - may be interpreted as integer and hence continuous)
  bm$species <- as.character(bm$species)
  # Due to log10, need to set a minimum value, seems like a feature in ggplot
  min_value <- 1e-300
  bm <- bm[bm$value >= min_value,]
  p <- ggplot(bm) + 
    geom_line(aes(x=time,y=value, colour=species, linetype=species)) + 
    scale_y_continuous(trans="log10", name="Biomass") + 
    scale_x_continuous(name="Time") +
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
    ggtitle(NULL)
  if (nrow(object@params@species_params)>12){
    p <- ggplot(bm) + geom_line(aes(x=time,y=value, group=species)) + scale_y_continuous(trans="log10", name="Biomass") + scale_x_continuous(name="Time") 
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Plot the total yield of each species through time
#'
#' After running a projection, the total yield of each species across all
#' fishing gears can be plotted against time. 
#' This plot is pretty easy to do by hand. It just gets the biomass using the \code{\link{getYield}} method and plots using the ggplot2 package. You can then fiddle about with colours and linetypes etc. Just look at the source code for details.
plotYield <- function(object, print_it = TRUE, ...){
  y <- getYield(object, ...)
  names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
  ym <- melt(y)
  p <- ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=Species)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
  if (nrow(object@params@species_params)>12){
    p <- ggplot(ym) + geom_line(aes(x=time,y=value, group=Species)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Plot the total yield of each species by gear through time
#'
#' After running a projection, the total yield of each species by 
#' fishing gear can be plotted against time. 
#' This plot is pretty easy to do by hand. It just gets the biomass using the \code{\link{getYieldGear}} method and plots using the ggplot2 package. You can then fiddle about with colours and linetypes etc. Just look at the source code for details.
plotYieldGear <- function(object, print_it=TRUE, ...){
  y <- getYieldGear(object, ...)
  names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
  ym <- melt(y)
  p <- ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=gear)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
  if (nrow(object@params@species_params)>12){
    p <- ggplot(ym) + geom_line(aes(x=time,y=value, group=Species)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Plot the abundance spectra of each species and the background population
#'
#' After running a projection, the spectra of the abundance of each species and the background population can be plotted.
#' The abundance is averaged over the specified time range (a single value for the time range can be used to plot a single time step).
#' The abundance can be in terms of numbers or biomass, depending on the \code{biomass} argument.
plotSpectra <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/100, biomass = TRUE, print_it = TRUE, ...){
  time_elements <- get_time_elements(object,time_range)
  spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
  background_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
  y_axis_name = "Abundance"
  if (biomass){
    spec_n <- sweep(spec_n,2,object@params@w,"*")
    background_n <- background_n * object@params@w_full
    y_axis_name = "Biomass"
  }
  # Make data.frame for plot
  plot_dat <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
  plot_dat <- rbind(plot_dat, data.frame(value = c(background_n), Species = "Background", w = object@params@w_full))
  # lop off 0s in background and apply min_w
  plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w),]
  p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = y_axis_name, trans="log10")
  if (nrow(object@params@species_params)>12){
    p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = y_axis_name, trans="log10")
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Plot the feeding level of each species by size 
#'
#' After running a projection, plot the feeding level of each species by size.
#' The feeding level is averaged over the specified time range (a single value for the time range can be used).
plotFeedingLevel <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
  feed_time <- getFeedingLevel(object=object, time_range=time_range, drop=FALSE, ...)
  feed <- apply(feed_time, c(2,3), mean)
  plot_dat <- data.frame(value = c(feed), Species = dimnames(feed)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
  p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Feeding Level", lim=c(0,1))
  if (nrow(object@params@species_params)>12){
    p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Feeding Level", lim=c(0,1))
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Plot M2 of each species by size 
#'
#' After running a projection, plot M2 of each species by size.
#' M2 is averaged over the specified time range (a single value for the time range can be used to plot a single time step).
plotM2 <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
  m2_time <- getM2(object, time_range=time_range, drop=FALSE, ...)
  m2 <- apply(m2_time, c(2,3), mean)
  plot_dat <- data.frame(value = c(m2), Species = dimnames(m2)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
  p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "M2", lim=c(0,max(plot_dat$value)))
  if (nrow(object@params@species_params)>12){
    p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "M2", lim=c(0,max(plot_dat$value)))
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Plot total fishing mortality of each species by size 
#'
#' After running a projection, plot the total fishing mortality of each species by size.
#' The total fishing mortality is averaged over the specified time range (a single value for the time range can be used to plot a single time step).
plotFMort <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
  f_time <- getFMort(object, time_range=time_range, drop=FALSE, ...)
  f <- apply(f_time, c(2,3), mean)
  plot_dat <- data.frame(value = c(f), Species = dimnames(f)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
  p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Total fishing mortality", lim=c(0,max(plot_dat$value)))
  if (nrow(object@params@species_params)>12){
    p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Total fishing mortality", lim=c(0,max(plot_dat$value)))
  }
  if (print_it)
    print(p)
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

#' Summary plot for \code{MizerSim} objects
#'
#' After running a projection, produces 5 plots in the same window:
#' feeding level, abundance spectra, predation mortality and fishing
#' mortality of each species by size; and biomass of each species through
#' time.
#' This method just uses the other plotting methods and puts them
#' all in one window.
plotSummary <- function(x, ...){
  p1 <- plotFeedingLevel(x,print_it = FALSE,...)
  p2 <- plotSpectra(x,print_it = FALSE,...)
  p3 <- plotBiomass(x,print_it = FALSE,...)
  p4 <- plotM2(x,print_it = FALSE,...)
  p5 <- plotFMort(x,print_it = FALSE,...)
  grid.newpage()
  glayout <- grid.layout(3,2) # widths and heights arguments
  vp <- viewport(layout = glayout)
  pushViewport(vp)
  vplayout <- function(x,y)
    viewport(layout.pos.row=x, layout.pos.col = y)
  print(p1+ theme(legend.position="none"), vp = vplayout(1,1))
  print(p3+ theme(legend.position="none"), vp = vplayout(1,2))
  print(p4+ theme(legend.position="none"), vp = vplayout(2,1))
  print(p5+ theme(legend.position="none"), vp = vplayout(2,2))
  print(p2+ theme(legend.position="right", legend.key.size=unit(0.1,"cm")), vp = vplayout(3,1:2))
}

# My plots --------------------------------------------------------------------------------------------------------

# plot biomass

plotDynamics <- function(object, phenotype = TRUE, bloodline = NULL, light = FALSE, print_it = T){
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  biomass <- getBiomass(object) # n * w * dw and sum by species
  
  SpIdx = NULL # getting rid of the species that went extinct at the beginning 
  for (i in unique(object@params@species_params$species))
    if (sum(biomass[,i]) != 0 & dim(object@params@species_params[object@params@species_params$species == i,])[1] != 1) 
      SpIdx = c(SpIdx,i)
  
  biomassTot = NULL
  biomassTemp = biomass
  colnames(biomassTemp) = object@params@species_params$species
  for (i in SpIdx)
  {
    biomassSp = biomassTemp[,which(colnames(biomassTemp) == i)]
    biomassSp = apply(biomassSp,1,sum)
    biomassTot = cbind(biomassTot,biomassSp)
  }
  colnames(biomassTot) = SpIdx
  
  # biomassTot is the biomass of species
  # biomass is the biomass of phenotypes
  plotBiom <- function(x,light)
  {
    Biom <- melt(x) # melt for ggplot
    names(Biom) = list("time","sp","value")
    # Due to log10, need to set a minimum value, seems like a feature in ggplot
    min_value <- 1e-300
    Biom <- Biom[Biom$value >= min_value,]
    # take the first digit of the species column and put it in a new column
    Biom$bloodline = sapply(Biom[,2], function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
    if(light) Biom= Biom[seq(1,dim(Biom)[1],2),] # if the data is too heavy and it takes to much time
    return(Biom)
  }
  
  BiomPhen = NULL
  if (phenotype) BiomPhen <- plotBiom(x = biomass,light =light)
  BiomSp <- plotBiom(biomassTot,light)
  
  
  # colourCount = length(unique(Biom$sp))
  # getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used
  
  if (is.null(bloodline))
  {
    p <-
      ggplot(BiomSp) +
      geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1.2) +
      geom_line(data = BiomPhen, aes(x = time, y = value, colour = as.factor(bloodline), group = sp), alpha = 0.2) +
      scale_y_log10(name = "Biomass in g.m^-3", limits = c(1e-11, NA), breaks = c(1 %o% 10^(-11:-1))) +
      scale_x_continuous(name = "Time in years") +
      #geom_hline(yintercept = 1e-11) +
      #geom_text(aes(0,1e-11,label = "Extinction threshold", vjust = 1, hjust = -0.15), size = 2.5) +
      labs(color='Species') +
      theme(legend.key = element_rect(fill = "white"))+
      theme_bw()+
      scale_colour_manual(values=cbPalette)+ # colorblind
      #scale_color_grey(name = "Species")+ # grey
      ggtitle("Community biomass")
  }
  
  else 
  {
    Sbm= Sbm[Sbm$bloodline == bloodline,]
    
    
    p <-
      ggplot(Sbm) +
      geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp)) +#, linetype = sp)) +
      scale_y_log10(name = "Biomass in g.m^-3", limits = c(1e-11, NA), breaks = c(1 %o% 10^(-11:-1))) +
      scale_x_continuous(name = "Time in years") +
      labs(color='Species') +
      theme(legend.key = element_rect(fill = "white"))+
      theme_bw()+
      scale_color_grey()+
      ggtitle("Community biomass")
  }
  
  if(print_it) print(p)
  return(p)
}

plotDynamicsMulti <- function(folder)
{
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind

  #NO FISHERIES PART
  path_to_sim = paste(folder,"/normal",sep="")
  bigSim <- bunchLoad(folder = path_to_sim)
  sim <- superStich(bigSim)
  
  #top left corner biomass of 1 stochastic run without fisheries
  plot_dat <- biom(bigSim[[1]])
  p1 <- ggplot(plot_dat[[1]]) +
    geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1.2) +
    geom_line(data = plot_dat[[2]], aes(x = time, y = value, colour = as.factor(bloodline), group = sp), alpha = 0.2) +
    scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(1e-7, NA), breaks = c(1 %o% 10^(-11:-1))) +
    scale_x_continuous(name = NULL,labels = NULL) +
    #scale_x_continuous(name = "Time in years") +
    labs(color='Species') +
    theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    scale_colour_manual(values=cbPalette)+ # colorblind
    ggtitle("a)")
  
  # bottom left corner, biomass of averaged run without fisheries
  plot_dat <- biom(sim,phenotype = F)
  p2 <- ggplot(plot_dat[[1]]) +
    geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1.2) +
    scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(1e-05, NA), breaks = c(1 %o% 10^(-5:-1))) +
    scale_x_continuous(name = "Time in years") +
    labs(color='Species') +
    theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    scale_colour_manual(values=cbPalette)+ # colorblind
    ggtitle("c)")
  
  #FISHERIES PART
  path_to_sim = paste(folder,"/fisheries",sep="")
  bigSim <- bunchLoad(folder = path_to_sim)
  # dirContent <- dir(path_to_sim)
  # bigSim = list()
  # for(i in 1:length(dirContent))
  # {
  #   if (file.exists(paste(path_to_sim,"/",dirContent[i],"/run.Rdata",sep = "")))
  #   {
  #     sim <- get(load(paste(path_to_sim,"/",dirContent[i],"/run.Rdata",sep="")))
  #   } else {
  #     sim <- get(load(paste(path_to_sim,"/",dirContent[i],"/defaultRun.Rdata",sep="")))
  #     sim = superOpt(sim)
  #   } 
  #   bigSim[[i]] = sim
  #   rm(sim)
  # }
  sim <- superStich(bigSim)
  
  #top right corner biomass of 1 stochastic run with fisheries
  plot_dat <- biom(bigSim[[1]])
  p3 <- ggplot(plot_dat[[1]]) +
    geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1.2) +
    geom_line(data = plot_dat[[2]], aes(x = time, y = value, colour = as.factor(bloodline), group = sp), alpha = 0.2) +
    scale_y_log10(name = NULL, limits = c(1e-7, NA), labels = NULL, breaks = c(1 %o% 10^(-11:-1))) +
    scale_x_continuous(name = NULL,labels = NULL) +
    labs(color='Species') +
    theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    scale_colour_manual(values=cbPalette)+ # colorblind
    ggtitle("b)")
  
  # bottom left corner, biomass of averaged run without fisheries
  plot_dat <- biom(sim,phenotype = F)
  p4 <- ggplot(plot_dat[[1]]) +
    geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1.2) +
    scale_y_log10(name = NULL, limits = c(1e-5, NA), labels = NULL, breaks = c(1 %o% 10^(-11:-1))) +
    scale_x_continuous(name = "Time in years") +
    labs(color='Species') +
    theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    scale_colour_manual(values=cbPalette)+ # colorblind
    ggtitle("d)")
  
  #multiplot
  path_to_png = paste(folder,"/popDynamics.png",sep="")
  png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 600)
  
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(nrow=2, ncol=2, 
                                             widths = unit(c(10,10), "cm"), 
                                             heights = unit(c(15,5), "cm")))) 
  print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) 
  print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) 
  print(p3, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) 
  print(p4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) 
  
  dev.off() 
}

# plot size spectrum
plotSS <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/100, biomass = TRUE, print_it = TRUE, species = TRUE, ...){
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  min_w = 0.001
  
  time_elements <- get_time_elements(object,time_range)
  spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
  background_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
  
  y_axis_name = "Abundance"
  if (biomass){
    spec_n <- sweep(spec_n,2,object@params@w,"*")
    background_n <- background_n * object@params@w_full
    y_axis_name = "Biomass"
  }
  
  # Make data.frame for plot
  plot_datSP <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)), bloodline = object@params@species_params$species)
  plot_datPkt <- data.frame(value = c(background_n), Species = "Background", w = object@params@w_full)
  
  if (species)
  {
    dimnames(spec_n)$species = object@params@species_params$species
    SpIdx = unique(object@params@species_params$species)
    spec_sp = matrix(data = NA, ncol = dim(spec_n)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(spec_n)$size))
    names(dimnames(spec_sp))=list("species","size")
    
    for (i in 1:dim(spec_sp)[1])
    {
      temp = spec_n # save to manip
      temp[which(rownames(spec_n) != i), ] = 0 # make everything but the targeted species to go 0 to have correct normalisation
      temp = apply(temp, 2, sum)
      spec_sp[i, ] = temp
    }
    plot_datSP <- data.frame(value = c(spec_sp), Species = dimnames(spec_sp)[[1]], w = rep(object@params@w, each=length(SpIdx)))
  }
  
  # lop off 0s in background and apply min_w
  plot_datSP <- plot_datSP[(plot_datSP$value > 0) & (plot_datSP$w >= min_w),]
  plot_datPkt <- plot_datPkt[(plot_datPkt$value > 0) & (plot_datPkt$w >= min_w),]
  #getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used
  
  if (species)
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value, colour = as.factor(Species), group = Species)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, linetype = Species), size = 1.5) +
      scale_x_log10(name = "Size in g", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_log10(name = "Abundance density in individuals.m^-3", limits = c(1e-8,1e4)) +
      theme(panel.background = element_blank(),legend.key = element_rect(fill = "white"))+
      theme_bw()+
      scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
      #scale_color_grey(name = "Species")+ # grey
      ggtitle("Size spectrum")
  }
  
  else
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value, colour = as.factor(bloodline), group = Species)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, colour = Species), size = 1.5) +
      scale_x_log10(name = "Size in g", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_log10(name = "Abundance density in individuals.m^-3", limits = c(1e-35,1e4)) +
      theme(panel.background = element_blank(),legend.key = element_rect(fill = "white"))+
      theme_bw()+
      scale_colour_manual(values=cbPalette)+ # colorblind
      #scale_color_grey(name = "Species")+ # grey
      ggtitle("Size spectrum")
    
  }
  
  
  if (print_it)
    print(p)
  return(p)
}

# fisheries mortality
plotFM <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  f_time <- getFMort(object, time_range=time_range, drop=FALSE, ...)
  f <- apply(f_time, c(2,3), mean)
  plot_dat <- data.frame(value = c(f), Species = dimnames(f)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
  
  # take the first digit of the species column and put it in a new column
  plot_dat$bloodline = sapply(plot_dat$Species, function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, group = bloodline,color = as.factor(bloodline))) +
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-1:5)), limits = c(1e-2,1e5)) + 
    scale_y_continuous(name = "Total fishing mortality", lim=c(0,max(plot_dat$value))) +
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    scale_colour_manual(values=cbPalette)+ # colorblind
    ggtitle(NULL)
  
  if (print_it)
    print(p)
  return(p)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# plot feeding / satiation level
plotFood <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, throughTime = F, start = 1000, every = 1000, print_it = T, returnData = F, ...){
  
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  if (throughTime)
  {
    time_range = seq(start,max(as.numeric(dimnames(object@n)$time)),every)
    time_range = c(time_range,max(as.numeric(dimnames(object@n)$time))) # so it counts the last time step which is probably not even
    time_range = unique(time_range)
    feeding = array(data = NA, dim = c(length(unique(object@params@species_params$species)),100,length(time_range)),  
                    dimnames = list(as.character(unique(object@params@species_params$species)),object@params@w,time_range)) 
    Critfeeding = matrix(data=NA, nrow = length(time_range), ncol= 100, dimnames = list(time_range,object@params@w))
    for (i in time_range)
    {
      
      feed_time <- getFeedingLevel(object=object, time_range=i, drop=FALSE)#, ...) # get the feeding time
      feed <- apply(feed_time, c(2,3), mean) # average on the time frame
      
      Cfeed_time <- getCFeedingLevel(object=object, time_range=i, drop=FALSE)#, ...) # get the critical feeding level
      Critfeed <- apply(Cfeed_time, c(2,3), mean) # average on the time frame
      Critfeed <- Critfeed[1,] # all rows the same
      
      dimnames(feed)$sp = object@params@species_params$species
      SpIdx = unique(object@params@species_params$species) # get the species names
      feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(feed)$w)) # prepare the new object
      names(dimnames(feed_sp))=list("species","size")
      
      for (j in 1:dim(feed_sp)[1])
      {
        temp = feed # save to manip
        temp[which(rownames(feed) != j), ] = 0 # keep the ecotypes from the species only
        temp = apply(temp, 2, sum)
        temp = temp / length(which(rownames(feed)==j)) # do the mean (in 2 steps)
        feed_sp[j, ] = temp
      }
      feeding[,,which(dimnames(feeding)[[3]] == i)] = feed_sp
      Critfeeding[which(dimnames(Critfeeding)[[1]] == i),] = Critfeed
    }
    
    a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
    vlines <- data.frame(xint = a,grp = c(1:9))
    
    plot_dat = melt(feeding)
    colnames(plot_dat) = c("species","size","time","value")
    plot_crit = melt(Critfeeding)
    colnames(plot_crit) = c("time","size","value")
    p <- ggplot(plot_dat) + 
      geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
      geom_line(data = plot_crit, aes(x = size, y = value), linetype = "dashed") +
      scale_x_log10(name = "Size") + 
      scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
      geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
      facet_grid(time ~ .)+
      scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"))+
      
      ggtitle("Feeding level through time")
    
    if(print_it) print(p)
    return(p)
    
  }
  
  
  feed_time <- getFeedingLevel(object=object, time_range=time_range, drop=FALSE, ...) # get the feeding time
  feed <- apply(feed_time, c(2,3), mean) # average on the time frame
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(feed)$sp = object@params@species_params$species
    SpIdx = unique(object@params@species_params$species) # get the species names
    feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(feed)$w)) # prepare the new object
    names(dimnames(feed_sp))=list("species","size")
    
    for (i in 1:dim(feed_sp)[1])
    {
      temp = feed # save to manip
      temp[which(rownames(feed) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(feed)==i)) # do the mean (in 2 steps)
      feed_sp[i, ] = temp
    }
    feed = feed_sp
  }
  
  plot_dat <- data.frame(value = c(feed), Species = dimnames(feed)[[1]], w = rep(object@params@w, each=length(dimnames(feed)[[1]])))
  
  name = paste("Feeding level at time",time_range,sep=" ")
  
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10") + 
    scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
    ggtitle(name)
  
  if(print_it) print(p)
  
  if (returnData) return(plot_dat) else return(p)
}

# plot of the number of ecotype per time step

PlotNoSp <- function(object, print_it = T, returnData = F, init = T, dt = 0.1){
  
  if (is.list(object)) # this means that it gets a list of sim and needs to do the average
  {
    avgList = list()
    for (x in 1:length(object))
    {
      sim = object[[x]]
      # I really need to fix the time max in the sims ;(;(;(
      if(!init) {
        TimeMax = sim@params@species_params$timeMax[1] + 40000 
        timeStart = 40000} else {
          TimeMax = sim@params@species_params$timeMax[1]
          timeStart = 1
        }
      
      SumPar = cbind(sim@params@species_params$pop,sim@params@species_params$extinct)
      colnames(SumPar) = c("Apparition","Extinction")
      rownames(SumPar) = sim@params@species_params$ecotype
      SumPar = SumPar[order(SumPar[,1],decreasing=FALSE),]
      SumPar = as.data.frame(SumPar) # little dataframe with species apparition and extinction
      for (i in 1:dim(SumPar)[1]) if (SumPar$Extinction[i] == 0) SumPar$Extinction[i] = TimeMax # the not extinct ones get the end sim as extinction value
      
      avg = matrix(data = 0, nrow = TimeMax, ncol =4)
      avg[,1] = seq(1:TimeMax)*0.1
      ApCount = length(which(SumPar$Apparition == 0)) # number of starting species
      ExCount = 0
      SumParEx = SumPar[order(SumPar$Extinction,decreasing=F),] #need this order for extinction in loop
      
      for (i in timeStart:TimeMax) # for each time step
      {
        for (j in 1:dim(SumPar)[1])
          if (SumPar$Apparition[j] <= i & SumPar$Extinction[j] >= i) # how many phenotypes are alive right now?
            avg[i,2] = avg[i,2] + 1
          
          if (i %in% SumPar$Apparition) ApCount = ApCount + 1 # how many in total
          avg[i,3] = ApCount
          
          if (i %in% SumParEx$Extinction) ExCount = ExCount + 1 #how many extinct?
          avg[i,4] = ExCount
      }
      
      avg = as.data.frame(avg)
      dimnames(avg)[[2]] = list("time", "number","apparition","extinction") # dataframe of the number of species at each time step
      avg = avg[seq(1,dim(avg)[1],10),]
      avgList[[x]] = avg
    }
    
    # here I have a list of the demographic data of each sim 
    # now I need to do stats on them
    
    #DemoMean = do.call(mean, )
    
    #ans1 = aaply(laply(avgList, as.matrix), c(2, 3), mean)
    
    avgMatrix <- abind(avgList, along=3) # convert the list in an array
    avgMean = apply(avgMatrix, c(1,2), mean) # do the stat manip
    avgSd = apply(avgMatrix, c(1,2), sd)
    
    # I have to check if the sd function works properly, but not right now
    #    f = apply((sweep(yo1,2,avgMean,"-"))^2,2,sum)/length(avgList) # this is the variance at each time step
    #   g = sqrt(f) # this is the standard population deviation at each time step
    # yo1 = avgList[[1]]
    # yo2 = avgList[[2]]
    
    stat = data.frame(avgMean,avgSd[,-1])# prepare my data
    colnames(stat) = c("time","Mno","Mpop","Mext","Sno","Spop","Sext")
    
    # first build the standard deviation curve only 
    g1 <- ggplot(stat)+
      geom_smooth(aes(x = time, y = (Mno-Sno))) +
      geom_smooth(aes(x=time, y = (Mno+Sno))) +
      geom_smooth(aes(x=time, y = (Mext-Sext)))+
      geom_smooth(aes(x=time, y = (Mext+Sext)))
    
    gg1 <- ggplot_build(g1)
    
    dfRibbon <- data.frame(xNo = gg1$data[[1]]$x, yminNo = gg1$data[[1]]$y,  ymaxNo = gg1$data[[2]]$y, #and extract the smooth data for the ribbon
                           xExt = gg1$data[[3]]$x, yminExt = gg1$data[[3]]$y,  ymaxExt = gg1$data[[4]]$y)
    
    
    p <- ggplot(stat) +
      stat_smooth(aes(x = time, y = Mno, color = "green")) +
      #geom_smooth(aes(x = time, y = (Mno-Sno))) +
      #geom_smooth(aes(x=time, y = (Mno+Sno)))+
      geom_ribbon(data = dfRibbon, aes(x = xNo, ymin = yminNo, ymax = ymaxNo),
                  fill = "grey", alpha = 0.4)+
      geom_smooth(aes(x=time, y = Mext, color = "blue"))+
      #geom_smooth(aes(x=time, y = (Mext-Sext)))+
      #geom_smooth(aes(x=time, y = (Mext+Sext)))+
      geom_ribbon(data = dfRibbon, aes(x = xExt, ymin = yminExt, ymax = ymaxExt),
                  fill = "grey", alpha = 0.4)+
      # geom_smooth(data = yo1 ,aes(x=time, y = number))+
      # geom_smooth(data = yo2 ,aes(x=time, y = number))+
      scale_x_continuous(name = "Time (yr)") +
      scale_y_continuous(name = "Number of phenotypes")+
      scale_colour_discrete(labels = c("Extinct","Alive"))+
      theme(legend.title=element_blank(),
            legend.position=c(0.19,0.95),
            legend.justification=c(1,1),
            legend.key = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      ggtitle("Variation of phenotype's number throughout the simulation")
    
    if (print_it)  print(p)
    
    if (returnData) return(list(stat,dfRibbon)) else return(p)
  }
  
  
  # I really need to fix the time max in the sims ;(;(;(
  if(!init) {
    TimeMax = object@params@species_params$timeMax[1] + 40000 
    timeStart = 40000} else {
      TimeMax = object@params@species_params$timeMax[1]
      timeStart = 1
    }
  SumPar = cbind(object@params@species_params$pop,object@params@species_params$extinct) # weird things happen without the as.numeric
  colnames(SumPar) = c("Apparition","Extinction")
  rownames(SumPar) = object@params@species_params$ecotype
  SumPar = SumPar[order(SumPar[,1],decreasing=FALSE),]
  SumPar = as.data.frame(SumPar) # little dataframe with species apparition and extinction
  for (i in 1:dim(SumPar)[1]) if (SumPar$Extinction[i] == 0) SumPar$Extinction[i] = TimeMax # the not extinct ones get the end sim as extinction value
  
  avg = matrix(data = 0, nrow = TimeMax, ncol =4)
  avg[,1] = seq(1:(TimeMax))
  ApCount = length(which(SumPar$Apparition < timeStart)) # number of starting species
  ExCount = length(which(SumPar$Extinction < timeStart))
  SumParEx = SumPar[order(SumPar$Extinction,decreasing=F),] #need this order for extinction in loop
  
  for (i in (timeStart):(TimeMax)) # for each time step
  {
    for (j in 1:dim(SumPar)[1])
      if (SumPar$Apparition[j] <= i & SumPar$Extinction[j] >= i) # how many phenotypes are alive right now?
        avg[i,2] = avg[i,2] + 1
      
      if (i %in% SumPar$Apparition) ApCount = ApCount + 1 # how many in total
      avg[i,3] = ApCount
      
      if (i %in% SumParEx$Extinction) ExCount = ExCount + 1 #how many extinct?
      avg[i,4] = ExCount
  }
  
  avg = as.data.frame(avg)
  dimnames(avg)[[2]] = list("time", "number","apparition","extinction") # dataframe of the number of species at each time step
  
  p <- ggplot(avg) +
    geom_line(aes(x = time, y = number, color = "green")) +
    #geom_line(aes(x=time, y = apparition, color = "blue"))+
    geom_line(aes(x=time, y = extinction, color = "red"))+
    scale_x_continuous(name = "Time (yr)") +
    scale_y_continuous(name = "Number of phenotypes")+
    #scale_colour_discrete(labels = c("Total","Alive","Extinct"))+
    scale_colour_discrete(labels = c("Alive","Extinct"))+
    theme(legend.title=element_blank(),
          legend.position=c(0.19,0.95),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    ggtitle("Variation of phenotype's number throughout the simulation")
  
  if (print_it)  print(p)
  
  if (returnData) return(list(stat,dfRibbon)) else return(p)
}

# my plot of the dead (x 3 because that's why)
#' Plot M2 of each species by size 

plotUdead <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),species = T, throughTime = F, print_it = TRUE, returnData = F, ...){
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = c(1:9))
  
  if (throughTime)
  {
    time_range = seq(500,max(as.numeric(dimnames(object@n)$time)),500)
    time_range = c(time_range,max(as.numeric(dimnames(object@n)$time))) # so it counts the last time step which is probably not even
    time_range = unique(time_range)
    mortality = matrix(data = NA, nrow = length(time_range), ncol = 100, dimnames = list(time_range,object@params@w))
    
    for (i in time_range)
    {
      m2_time <- getM2(object, time_range=i, drop=FALSE, ...)
      m2 <- apply(m2_time, c(2,3), mean)
      mortality[which(rownames(mortality)==i),] <- m2[1,] #all rows are the same
    }
    
    
    plot_dat = melt(mortality)
    colnames(plot_dat) = list("time","size","value")
    
    p <- ggplot(plot_dat)+
      geom_line(aes(x=size, y=value, color = as.factor(time)))+
      scale_x_log10(name = "Prey size in g") +
      scale_y_continuous(name = "Predation mortality in g.m^-3") +
      #geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
      geom_vline(xintercept = c(object@params@species_params$w_inf[1:9]), linetype = "dashed")+
      scale_colour_manual(values=cbPalette, name = "Time in years")+ # colorblind
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"))+
      ggtitle("Predation mortality through time")
    
    if (print_it)
      print(p)
    return(p)
  }
  
  m2_time <- getM2(object, time_range=time_range, drop=FALSE, ...)
  m2 <- apply(m2_time, c(2,3), mean)
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(m2)$prey = object@params@species_params$species
    SpIdx = unique(object@params@species_params$species) # get the species names
    m2_sp = matrix(data = NA, ncol = dim(m2)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(m2)$w_prey)) # prepare the new object
    names(dimnames(m2_sp))=list("prey","w_prey")
    
    for (i in 1:dim(m2_sp)[1])
    {
      temp = m2 # save to manip
      temp[which(rownames(m2) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(m2)==i)) # do the mean (in 2 steps)
      m2_sp[i, ] = temp
    }
    m2 = m2_sp
  }
  name = paste("Predation Mortality at time",time_range,sep=" ")
  
  plot_dat <- data.frame(value = c(m2), Species = dimnames(m2)[[1]], w = rep(object@params@w, each=length(dimnames(m2)[[1]])))
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10") + 
    scale_y_continuous(name = "M2", lim=c(0,max(plot_dat$value))) +
    ggtitle(name)
  
  if (print_it) print(p)
  if(returnData) return(plot_dat) else return(p)
}

plotTotMort <- function(directory, predMort = NULL, time_range = NULL, print_it =F, returnData = F)
{
  # let's assume that object is a list of different sims
  object <- bunchLoad(folder = paste(directory,"/normal",sep=""))
  
  plotList <- list()
  for (i in 1:length(object))
  {
    sim <- object[[i]]
    if (is.null(time_range)) time_range =  max(as.numeric(dimnames(sim@n)$time))
    
    # predation mortality
    if (is.null(predMort)) 
    {
      m2_time <- getM2(sim, time_range=time_range, drop=FALSE)#, ...)
      m2 <- apply(m2_time, c(2,3), mean) 
    } else {
      m2 <- matrix(data = predMort[i,], nrow = length(sim@params@species_params$ecotype) , ncol = dim(predMort)[2],byrow = T, dimnames = list(sim@params@species_params$ecotype,sim@params@w) )
      names(dimnames(m2)) = c("prey","w_prey")
    }
    
    # fisheries mortality
    f_time <- getFMort(sim, time_range=time_range, drop=FALSE) #, ...)
    f <- apply(f_time, c(2,3), mean)
    # background mortality
    
    bm <- sim@params@species_params$z0
    
    # and let's assume I'm doing only single species runs for now
    m2S <- m2[1,]
    fS <- f[1,]
    bmS <- rep(bm[1],length(m2S))
    
    dead_dat <- data.frame(m2S,fS,bmS,sim@params@w,row.names = NULL)
    colnames(dead_dat) <- c("predation","fisheries","background","size")
    
    p <-  ggplot(dead_dat) +
      geom_line(aes(x = size, y = predation, color = "predation")) +
      geom_line(aes(x = size, y = fisheries, color = "fisheries")) +
      geom_line(aes(x = size, y = background, color = "background")) +
      scale_y_continuous(name = "mortality") +
      scale_x_continuous(name = "size", trans = "log10")
    
    temp <- ggplot_build(p) 
    temp$data[[1]]$group = i
    
    plotList[[i]] <- temp
  }
  
  dat = do.call(rbind,lapply(plotList,function(x) x$data[[1]]))
  
  p <- ggplot(dat) +
    geom_line(aes(x=x,y=y,color=as.factor(group)))+
    scale_x_continuous(name="size in log") +
    scale_y_continuous(name = "mortality in g", limits = c(0,1)) +
    scale_colour_grey(name = "Eta")+ 
    ggtitle("Predation mortality")
  
  if (print_it) print(p)
  
  if(returnData) return(dat) else return(p)
  
}

plotScythe <- function(object, whatTime = max(as.numeric(dimnames(object@n)$time)),print_it = TRUE, returnData = F, comments = T)
{

  z <- getZ(object = object@params, n = object@n[whatTime,,], n_pp = object@n_pp[whatTime,], effort = object@effort[whatTime])
  dimnames(z)$prey = object@params@species_params$species
  SpIdx = unique(object@params@species_params$species) # get the species names
  
  z_sp = matrix(data = NA, ncol = dim(z)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(z)$w_prey)) # prepare the new object
  names(dimnames(z_sp))=list("prey","w_prey")
    
    for (i in 1:dim(z_sp)[1])
    {
      temp = z # save to manip
      temp[which(rownames(z) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(z)==i)) # do the mean (in 2 steps)
      z_sp[i, ] = temp
    }
    z = z_sp
  
  name = paste("Predation Mortality at time",whatTime,sep=" ")
  
  plot_dat <- data.frame(value = c(z), Species = dimnames(z)[[1]], w = rep(object@params@w, each=length(dimnames(z)[[1]])))
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10") + 
    scale_y_continuous(name = "M2", lim=c(0,max(plot_dat$value))) +
    ggtitle(name)
  
  if (print_it) print(p)
  if(returnData) return(plot_dat) else return(p)

}

# plot growth at time t

plotGrowth <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F,...){
  
  time_elements <- get_time_elements(object,time_range)
  growth_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    growth <- getEGrowth(object@params, n=n, n_pp = object@n_pp[x,])
    return(growth)})
  
  #growth <- apply(growth_time, c(2,3), mean) # use this when I will have time_range on more than one time
  growth = growth_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(growth)$sp = object@params@species_params$species
    SpIdx = unique(object@params@species_params$species) # get the species names
    growth_sp = matrix(data = NA, ncol = dim(growth)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(growth)$w)) # prepare the new object
    names(dimnames(growth_sp))=list("species","size")
    
    for (i in 1:dim(growth_sp)[1])
    {
      temp = growth # save to manip
      temp[which(rownames(growth) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(growth)==i)) # do the mean (in 2 steps)
      growth_sp[i, ] = temp
    }
    growth = growth_sp
  }
  
  name = paste("Growth level at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(growth), Species = dimnames(growth)[[1]], w = rep(object@params@w, each=length(dimnames(growth)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10") + 
    scale_y_continuous(name = "growth Level", trans ="log10")+
    ggtitle(name)
  
  if(print_it) print(p)
  
  if (returnData) return(plot_dat) else return(p)
}

# growth curve
plotGrowthSpeed <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = T, ...){
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  # getEgrowth is my g at time t, but with varying feeding level (from the sim)
  # could do an option for constant feeding level
  time_elements <- get_time_elements(object,time_range)
  growth_time <- aaply(which(time_elements), 1, function(x){
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    growth <- getEGrowth(object@params, n=n, n_pp = object@n_pp[x,])
    return(growth)})
  
  rownames(growth_time) = object@params@species_params$species # phenotypes from same species have the same name
  
  SpIdx = NULL
  for (i in unique(object@params@species_params$species))
    if (sum(growth_time[,i]) != 0 & dim(object@params@species_params[object@params@species_params$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
      SpIdx = c(SpIdx,i)
  
  growth_timeSp=NULL
  for (i in SpIdx)
  {
    gtSp = growth_time # save to manipulate
    gtSp = gtSp[which(rownames(growth_time) == i),] # keep only the phenotypes from species i
    gtSpMean = apply(gtSp,2,mean) # I get the mean species at each size
    growth_timeSp = rbind(growth_timeSp,gtSpMean)
  }
  
  growth_time = growth_timeSp
  dimnames(growth_time)[[1]] = SpIdx
  
  arrivalM = matrix(nrow = dim(growth_time)[1], ncol = dim(growth_time)[2], dimnames = list(rownames(growth_time),colnames(growth_time)))
  for (i in 1:dim(arrivalM)[1]) # for every species
    for (j in 1:dim(arrivalM)[2]) # for every size
      if (growth_time[i,j] > 0) 
        arrivalM[i,j] = growth_time[i,j] + object@params@w[j] # matrix with the biomass put into growth (from birth I guess)
  
  arrivalM = arrivalM/object@params@species_params$w_inf[SpIdx] # normalise the growth value (but not the size bin)
  
  MSp = NULL
  for (i in SpIdx) # for every species
    MSp = rbind(MSp,cbind(rep(i,dim(arrivalM)[2]),arrivalM[which(i == rownames(arrivalM)),],as.numeric(colnames(arrivalM))/object@params@species_params$w_inf[i])) #normalise their size bins
  
  colnames(MSp) = list("species","value","wnorm") # this gets a 3 col matrix with all the info I need for ggplot
  rownames(MSp) = NULL
  MSP = as.data.frame(MSp)
  
  p <- ggplot(MSP)+
    geom_smooth(aes(x=wnorm,y=value,group=species, color = as.factor(species)), se=F, method = "loess", na.rm = T)+
    scale_x_continuous(limits = c(0,1), name = "Body size in g (scaled with M)")+
    scale_y_continuous(limits = c(0,1), name = "Body size in g at next time step (scaled with M)") +
    geom_abline(slope = 1, linetype = "dashed" )+
    theme(panel.background = element_blank(),legend.key = element_rect(fill = "white"))+
    theme_bw()+
    scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
    ggtitle("Growth speed")
  
  if (print_it) print(p)
  
  return(p)
  
}


plotGrowthCurve <- function(object,time_range = max(as.numeric(dimnames(object@n)$time)), print_it = T, generation = 30, SpIdx = NULL, ...){
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  #in case of one sp only
  #if (length(dim(object@n) <3)) object@n <- array(data = object@n,dim = c(dim(object@n)[1],1,dim(object@n)[2]),dimnames = list(dimnames(object@n)[[1]],"1",dimnames(object@n)[[2]]))
  
  # getEgrowth is my g at time t, but with varying feeding level (from the sim)
  # could do an option for constant feeding level
  time_elements <- get_time_elements(object,time_range)
  growth_time <- aaply(which(time_elements), 1, function(x){
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    growth <- getEGrowth(object@params, n=n, n_pp = object@n_pp[x,])
    return(growth)})
  
  rownames(growth_time) = object@params@species_params$species # phenotypes from same species have the same name
  
  # cat("growth\n")
  # print(growth_time)
  
  if(is.null(SpIdx))
    for (i in unique(object@params@species_params$species))
      if (sum(growth_time[,i]) != 0 & dim(object@params@species_params[object@params@species_params$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
        SpIdx = c(SpIdx,i)
  
  growth_timeSp=NULL
  for (i in SpIdx)
  {
    gtSp = growth_time # save to manipulate
    gtSp = gtSp[which(rownames(growth_time) == i),] # keep only the phenotypes from species i
    if (!is.null(dim(gtSp))) gtSpMean = apply(gtSp,2,mean) else gtSpMean = gtSp # I get the mean species at each size
    growth_timeSp = rbind(growth_timeSp,gtSpMean)
  }
  
  growth_time = growth_timeSp
  dimnames(growth_time)[[1]] = SpIdx
  
  # cat("growthSp\n")
  # print(growth_timeSp)
  
  arrivalM = matrix(nrow = dim(growth_time)[1], ncol = dim(growth_time)[2], dimnames = list(rownames(growth_time),colnames(growth_time)))
  for (i in 1:dim(arrivalM)[1]) # for every species
    for (j in 1:dim(arrivalM)[2]) # for every size
      if (growth_time[i,j] > 0) 
        arrivalM[i,j] = growth_time[i,j] + object@params@w[j] # matrix with the biomass put into growth (from birth I guess)
  
  
  #replace the NA in arrivalM by respective w_inf
  for (i in as.numeric(rownames(arrivalM)))
    for(j in 1:dim(arrivalM)[2])
      if (is.na(arrivalM[which(i == rownames(arrivalM)),j]))
        arrivalM[which(i == rownames(arrivalM)),j] = object@params@species_params$w_inf[i]
  
  time = seq(1:generation) #generation time to plot
  sizeAge = matrix(data=NA,nrow = length(SpIdx),ncol = length(time), dimnames = list(as.character(SpIdx),as.character(time)))
  
  sizeAge[,1] = arrivalM[,1]
  
  for(i in SpIdx)
  {
    for (j in time[-1])
    {
      if (sizeAge[which(i == rownames(sizeAge)), j - 1] != object@params@species_params$w_inf[which(object@params@species_params$ecotype == i)])
      {
        look = object@params@w[which.min(abs(object@params@w - sizeAge[which(i == rownames(sizeAge)), j - 1]))] #whats the size class the closest
        if ((sizeAge[which(i == rownames(sizeAge)), j - 1] - look) > 0)  {
          encadrant = c(which.min(abs(object@params@w - sizeAge[which(i == rownames(sizeAge)), j - 1])), which.min(abs(object@params@w - sizeAge[which(i == rownames(sizeAge)), j - 1])) + 1) #determine between which size class the size is
        } else if ((sizeAge[which(i == rownames(sizeAge)), j - 1] - look) < 0) {
          encadrant = c(which.min(abs(object@params@w - sizeAge[which(i == rownames(sizeAge)), j - 1])) - 1 , which.min(abs(object@params@w - sizeAge[which(i == rownames(sizeAge)), j - 1])))
        } else
        {
          encadrant = which.min(abs(object@params@w - sizeAge[which(i == rownames(sizeAge)), j - 1]))
        } # need to do something when this happens
        ratioSize = abs(object@params@w[encadrant] - sizeAge[which(i == rownames(sizeAge)), j - 1]) / (object@params@w[encadrant][2] - object@params@w[encadrant][1])# tells me the normalised distance to the bins
        
        nextsize = arrivalM[which(i == rownames(arrivalM)), encadrant] # bins growth at t+2
        sizeAge[which(i == rownames(sizeAge)), j] = (nextsize[2] - nextsize[1]) * ratioSize[1] + nextsize[1] #real size
        if (sizeAge[which(i == rownames(sizeAge)), j] > object@params@species_params$w_inf[which(object@params@species_params$ecotype == i)])
          sizeAge[which(i == rownames(sizeAge)), j] = object@params@species_params$w_inf[which(object@params@species_params$ecotype == i)] #cannot go over w_inf
        
      } else {
        sizeAge[which(i == rownames(sizeAge)), j] = object@params@species_params$w_inf[which(object@params@species_params$ecotype == i)]
      }
    }
  }
  # sizeAge summarise the size each species reach after j time steps
  
  sizeAge = cbind(rep(object@params@w[1],dim(sizeAge)[1]),sizeAge) # add first column (starting at min_w)
  sizeAgeN = sizeAge/object@params@species_params$w_inf[SpIdx] # normalise for comparison
  colnames(sizeAgeN) = c(0,time)
  
  # cat("size at age\n")
  # print(sizeAge)
  
  dataAge = melt(sizeAgeN)
  colnames(dataAge) = c("species","generation","size")
  
  p <- ggplot(dataAge) +
    geom_line(aes(x=generation,y=size,color = as.factor(species)), na.rm = T) +
    scale_x_continuous(name = "generation time") +
    scale_y_continuous(name = "normalised size") +
    geom_hline(yintercept = 0.25, linetype = "dashed") +
    #scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle("Size at age")
  
  if(print_it) print(p)
  return(p)
}

# traits plots -----------

plotTraits <- function(object, maturationSize = T, ppmr = T, dietBreadth = T, print_it = F)
{
  
  # prepare the data
  # little initialisation 
  SumPar = object@params@species_params #shortcut
  TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),0.1*SumPar$pop,0.1*SumPar$extinct,SumPar$w_mat,SumPar$beta,SumPar$sigma) # weird things happen without the as.numeric / *0.1 because dt
  colnames(TT) = c("Lineage","Ecotype","Apparition","Extinction","Maturation_size","PPMR","Diet_breadth")
  rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  TT = as.data.frame(TT) # I use TT later in the graphs (its like an artefact)
  for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = SumPar$timeMax[i]*0.1
  
  # Weighted mean
  # 1) matrix of summed abundance of mature ind at each time step
  truc = object@n
  # put 0 in object@n when w < w_mat
  for (i in 1:dim(truc)[1]) # for each time step
  {
    for (j in 1:dim(truc)[2]) # for each ecotypes
    {
      w_lim = SumPar$w_mat[j] # get the maturation size of the ecotype
      S <- numeric(length(object@params@w))
      S[sapply(w_lim, function(i) which.min(abs(i - object@params@w)))] <- 1 # find what w bin is the closest of the maturation size
      NoW_mat = which(S == 1) # what is the col number of that size bin
      truc[i,j,1:NoW_mat-1] <-0 # everything under this value become 0
    }
  }
  abundanceM = apply(truc, c(1,2),sum) # sum the abundance left 
  
  # 2) normalisation per species 
  colnames(abundanceM) = SumPar$species # phenotypes from same species have the same name
  abundanceNormal = matrix(0,nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])
  
  # I am getting rid of the species which went instinct at the begining and that went extinct without having mutants (no trait variation)
  SpIdx = NULL
  for (i in unique(SumPar$species))
    if (sum(abundanceM[,i]) != 0 & dim(SumPar[SumPar$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
      SpIdx = c(SpIdx,i)
  
  for (i in SpIdx)
  {
    abundanceSp = abundanceM # save to manipulate
    abundanceSp[,which(colnames(abundanceM) != i)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
    abundanceSp = sweep(abundanceSp,1,apply(abundanceSp,1,sum),"/") # normalise
    abundanceSp[is.nan(abundanceSp)] <-0 # when I divide by 0 I get nan
    abundanceNormal = abundanceNormal + abundanceSp # I just need to add them up to get the final matrix
  }
  
  # Now I have normalised abundance, I need to apply the trait I want to plot on them
  
  # Maturation size
  if (maturationSize)
  {
    abundanceT = sweep(abundanceNormal,2,SumPar$w_mat,"*") # I use the normalised abundance and multiply by the trait value
    
    # Calculate mean at each time step
    TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
    names(dimnames(TotMean)) = list("time","species")
    
    for (i in SpIdx)
    {
      AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
      if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
      TotMean[,i] = AMean
    }
    
    # Calculate variance and standard deviation
    # it is the sum of the difference between value and mean squared and multiplied by the weight
    
    statMS = list() # list with the stats of all species as I need to do this for each species separatly
    
    for (i in SpIdx)
    {
      meanSp = TotMean[,i] # take the mean of the species
      traitSp = SumPar$w_mat[SumPar$species == i] # take the traits of the ecotypes in the species
      weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
      stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","variance","sd"))) # initialise the matrix
      for (j in 1:length(meanSp)) # for each time step
      {
        variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
        stat[j,3] = variance
        stat[j,4] = sqrt(variance) # calculate the standard deviation
      }
      statMS[[i]] = as.data.frame(stat) # put in the list
    }
    
    # I have the stats for every species, just need to plot now
    plotMSstore = list()
    for (i in SpIdx)
    {
      stat = statMS[[i]] # take the stats of the right species
      stat$up = stat$mean + stat$sd
      stat$low = stat$mean - stat$sd
      
      phenotype = TT[TT$Lineage == i,3:5] # recup the traits time values
      Phen = melt(phenotype,"Maturation_size") # make the dataframe
      name = paste("Maturation size of species ",i, sep = "")
      
      #short cut the data frame when species does not reach end of simulation
      if (sum(which(Phen$value == SumPar$timeMax[[1]]*0.1)) == 0) # if there is at least one value equal to the end of the sime it means that the species survived until then
        stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
      
      p = ggplot() +
        geom_point(data = Phen, aes(x = value, y = Maturation_size, group = Maturation_size)) +
        geom_line(data = Phen, aes(x = value, y = Maturation_size, group = Maturation_size)) + 
        geom_smooth(data = stat, aes(x = time, y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = time, y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = time, y = up, color ="blue")) +
        scale_x_continuous(name = "Time", limits  = c(0,SumPar$timeMax[i]*0.1)) +
        scale_y_continuous(name = "Trait value") +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position=c(0.9,1), 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle(name)
      
      
      if (print_it) print(p)
      
      plotMSstore[[i]] = p
    }
  }
  # PPMR 
  if (ppmr)
  {
    abundanceT = sweep(abundanceNormal,2,SumPar$beta,"*") # I use the normalised abundance and multiply by the trait value
    
    # Calculate mean at each time step
    TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
    names(dimnames(TotMean)) = list("time","species")
    
    for (i in SpIdx)
    {
      AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
      if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
      TotMean[,i] = AMean
    }
    
    # Calculate variance and standard deviation
    # it is the sum of the difference between value and mean squared and multiplied by the weight
    
    statPPMR = list() # list with the stats of all species
    for (i in SpIdx)
    {
      meanSp = TotMean[,i] # take the mean of the species
      traitSp = SumPar$beta[SumPar$species == i] # take the traits of the ecotypes in the species
      weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
      stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","variance","sd"))) # initialise the matrix
      for (j in 1:length(meanSp)) # for each time step
      {
        variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
        stat[j,3] = variance
        stat[j,4] = sqrt(variance) # calculate the standard deviation
      }
      statPPMR[[i]] = as.data.frame(stat) # put in the list
    }
    
    # I have the stats for every species, just need to plot now
    plotPPMRstore = list()
    for (i in SpIdx)
    {
      stat = statPPMR[[i]] # take the stats of the right species
      stat$up = stat$mean + stat$sd
      stat$low = stat$mean - stat$sd
      
      phenotype = TT[TT$Lineage == i,c(3,4,6)] # recup the traits time values
      Phen = melt(phenotype,"PPMR") # make the dataframe
      name = paste("PPMR of species ",i, sep = "")
      
      #short cut the data frame when species does not reach end of simulation
      if (sum(which(Phen$value == SumPar$timeMax[[1]]*0.1)) == 0) # if there is at least one value equal to the end of the sime it means that the species survived until then
        stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
      
      p = ggplot() +
        geom_point(data = Phen, aes(x = value, y = PPMR, group = PPMR)) +
        geom_line(data = Phen, aes(x = value, y = PPMR, group = PPMR)) + 
        geom_smooth(data = stat, aes(x = time, y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = time, y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = time, y = up, color ="blue")) +
        scale_x_continuous(name = "Time", limits  = c(0,SumPar$timeMax[i]*0.1)) +
        scale_y_continuous(name = "Trait value") +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position=c(0.9,1), 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle(name)
      
      if (print_it) print(p)
      
      plotPPMRstore[[i]] = p
    }
  }
  # Diet breath
  if (dietBreadth)
  {
    abundanceT = sweep(abundanceNormal,2,SumPar$sigma,"*") # I use the normalised abundance and multiply by the trait value
    
    # Calculate mean at each time step
    TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
    names(dimnames(TotMean)) = list("time","species")
    
    for (i in SpIdx)
    {
      AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
      if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
      TotMean[,i] = AMean
    }
    
    # Calculate variance and standard deviation
    # it is the sum of the difference between value and mean squared and multiplied by the weight
    
    statDB = list() # list with the stats of all species
    for (i in SpIdx)
    {
      meanSp = TotMean[,i] # take the mean of the species
      traitSp = SumPar$sigma[SumPar$species == i] # take the traits of the ecotypes in the species
      weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
      stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","variance","sd"))) # initialise the matrix
      for (j in 1:length(meanSp)) # for each time step
      {
        variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
        stat[j,3] = variance
        stat[j,4] = sqrt(variance) # calculate the standard deviation
      }
      statDB[[i]] = as.data.frame(stat) # put in the list
      #statDB[[i]] = stat# put in the list
    }
    
    # I have the stats for every species, just need to plot now
    plotDBstore = list()
    for (i in SpIdx)
    {
      stat = statDB[[i]] # take the stats of the right species
      stat$up = stat$mean + stat$sd
      stat$low = stat$mean - stat$sd
      
      phenotype =  TT[TT$Lineage == i,c(3,4,7)] # recup the traits time values
      Phen = melt(phenotype,"Diet_breadth") # make the dataframe
      name = paste("Diet breadth of species ",i, sep = "")
      
      #short cut the data frame when species does not reach end of simulation
      if (sum(which(Phen$value == SumPar$timeMax[[1]]*0.1)) == 0) # if there is at least one value equal to the end of the sime it means that the species survived until then
        stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
      
      p = ggplot() +
        geom_line(data = Phen, aes(x = value, y = Diet_breadth, group = Diet_breadth)) +
        geom_point(data = Phen, aes(x = value, y = Diet_breadth, group = Diet_breadth)) + 
        geom_smooth(data = stat, aes(x = time, y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = time, y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = time, y = up, color ="blue")) +
        scale_x_continuous(name = "Time", limits  = c(0,SumPar$timeMax[i]*0.1)) +
        scale_y_continuous(name = "Trait value") +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position=c(0.9,1), legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle(name)
      
      if (print_it) print(p)
      
      plotDBstore[[i]] = p
    }
  }
  
  return(list(plotMSstore,plotPPMRstore,plotDBstore))
}


# multi plots, will be included in the previous function
plotTraitsMulti <- function(object, print_it = F, save_it = F,dir = NULL, res = 600, SpIdx = NULL, Mat = T, PPMR = T, Sig = T, returnData =F, comments = T, window = NULL, Wmat = NULL, TimeMax = NULL, Normalisation = T)
{
  dt = 0.1 # because it is always here somewhere
  
  if (comments) cat("plotTraitsMulti begins\n") 
  
  if (is.null(SpIdx)) SpIdx <- seq(1,9)
  if (comments) 
   {cat("SpIdx set to:")
  print(SpIdx)}
  
  # determine the initial maturation size for normalisation purpose
  # if (Mat && is.null(Wmat))
  # {
    #Wmat = c(25, 79.056942, 250, 790.569415, 2500, 7905.694150, 25000)
    #Wmat = object@params@species_params$w_mat[1:length(object@params@species_params$species)]
    # if (sum(SpIdx == seq(1,9)) == length(SpIdx)) Wmat = c(2.5, 7.905694, 25, 79.056942, 250, 790.569415, 2500, 7905.694150, 25000) #eta = 0.25
    if (sum(SpIdx == seq(1,9)) == length(SpIdx)) Wmat = c(5.00000, 15.81139, 50.00000, 158.11388, 500.00000, 1581.13883, 5000.00000, 15811.38830, 50000.00000) #eta = 0.5
    if (sum(SpIdx == seq(1,3)) == length(SpIdx)) Wmat = c(5,50,500)
    
    if (comments) {
      cat("Wmat is")
      print(Wmat)
    }
  #}

  
  # beforehand
  if (save_it & is.null(dir))  # if want to save but not specified where
    dir = paste(getwd(),"/temporary",sep = "")
  
  ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir)), FALSE) #create the file if it does not exists
  
  if (comments) cat(sprintf("Starting trait plot\n"))
  
  # prepare the data
  # little initialisation 
  SumPar = object@params@species_params #shortcut
  TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),0.1*SumPar$pop,0.1*SumPar$extinct,SumPar$w_mat,SumPar$beta,SumPar$sigma) # weird things happen without the as.numeric / *0.1 because dt
  colnames(TT) = c("Species","Phenotype","Apparition","Extinction","Maturation_size","PPMR","Width_Feed")
  rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  TT = as.data.frame(TT) # I use TT later in the graphs (its like an artefact)
  
  # because I am starting simulations with previous abundance but not updating the time max (because I forgot and I do not want to run all my sims again) I need to do it here (which will disapear once I do the update in model.r)
  
  if (is.null(TimeMax)) TimeMax = 4e4 #default initialisation period
  TimeMax <- (SumPar$timeMax[1] + TimeMax) * dt
  for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = TimeMax
  
  if (is.null(window)) window = c(-1,1)
  if (comments) cat(sprintf("windows set from %g to %g\n",window[1],window[2]))
  # Weighted mean
  # 1) matrix of summed abundance of mature ind at each time step
  truc = object@n
  # put 0 in object@n when w < w_mat
  for (i in 1:dim(truc)[1]) # for each time step
  {
    for (j in 1:dim(truc)[2]) # for each ecotypes
    {
      w_lim = SumPar$w_mat[j] # get the maturation size of the ecotype
      S <- numeric(length(object@params@w))
      S[sapply(w_lim, function(i) which.min(abs(i - object@params@w)))] <- 1 # find what w bin is the closest of the maturation size
      NoW_mat = which(S == 1) # what is the col number of that size bin
      truc[i,j,1:NoW_mat-1] <-0 # everything under this value become 0
    }
  }
  abundanceM = apply(truc, c(1,2),sum) # sum the abundance left 
  
  # 2) normalisation per species 
  colnames(abundanceM) = SumPar$species # phenotypes from same species have the same name
  abundanceNormal = matrix(0,nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])
  
  # I am getting rid of the species which went instinct at the begining and that went extinct without having mutants (no trait variation)
  # SpIdx = NULL
  # for (i in unique(SumPar$species))
  #   if (sum(abundanceM[,i]) != 0 & dim(SumPar[SumPar$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
  #     SpIdx = c(SpIdx,i)
  
  # SpIdx is annoying:
  if (length(SpIdx) > length(unique(SumPar$species))) SpIdx = unique(SumPar$species) # ok so I might have species not even reaching this point so I'm short cuting spidx automatcaly
  
  
  if(comments) cat(sprintf("Spidx set to\n"))
  if (comments) print(SpIdx)
  
  for (i in SpIdx)
  {
    abundanceSp = abundanceM # save to manipulate
    abundanceSp[,which(colnames(abundanceM) != i)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
    abundanceSp = sweep(abundanceSp,1,apply(abundanceSp,1,sum),"/") # normalise
    abundanceSp[is.nan(abundanceSp)] <-0 # when I divide by 0 I get nan
    abundanceNormal = abundanceNormal + abundanceSp # I just need to add them up to get the final matrix
  }
  if (comments) cat(sprintf("Abundance normalised\n"))
  # Now I have normalised abundance, I need to apply the trait I want to plot on them
  plotMatStore <- list()
  plotPPMRStore <- list()
  plotDBStore <- list()
  
  # MATURATION SIZE
  if (Mat)
  {
    abundanceT = sweep(abundanceNormal,2,SumPar$w_mat,"*") # I use the normalised abundance and multiply by the trait value
    
    # Calculate mean at each time step
    #ncol = length(unique(SumPar$species))
    TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = max(unique(SumPar$species)), dimnames = list(rownames(abundanceT),as.character(seq(1,max(unique(SumPar$species)))))) #sort(unique(SumPar$species))[length(unique(SumPar$species))]
    names(dimnames(TotMean)) = list("time","species")
    
    for (i in SpIdx)
    {
      AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
      if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
      TotMean[,i] = AMean
    }
    
    # Calculate variance and standard deviation
    # it is the sum of the difference between value and mean squared and multiplied by the weight
    
    statMS = list() # list with the stats of all species as I need to do this for each species separatly
    
    for (i in SpIdx)
    {
      meanSp = TotMean[,i] # take the mean of the species
      traitSp = SumPar$w_mat[SumPar$species == i] # take the traits of the ecotypes in the species
      weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
      stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 5,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","percentMean","percentSd"))) # initialise the matrix
      for (j in 1:length(meanSp)) # for each time step
      {
        if (is.null(dim(weightSp))) {variance = sum(((traitSp-meanSp[j])^2)*weightSp[j])} else {variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,])} # calculate the variance, condition if only one phen
        stat[j,3] = sqrt(variance) # calculate the standard deviation
        # stat[j,4] = (meanSp[j] - SumPar$w_mat[i])/SumPar$w_mat[i] # normalisation of mean
        # stat[j,5] = stat[j,3]/SumPar$w_mat[i] # normalised sd
        stat[j,4] = (meanSp[j] - Wmat[i])/Wmat[i] # normalisation of mean
        stat[j,5] = stat[j,3]/Wmat[i] # normalised sd
      }
      statMS[[i]] = as.data.frame(stat) # put in the list
    }
    
    #determine the y axis limit
    # do something here
    
    
    # I have the stats for every species, just need to plot now
    
    
    for (i in SpIdx)
    {
      stat = statMS[[i]] # take the stats of the right species
      
      phenotype = TT[TT$Species == i,3:5] # recup the traits time values
      Phen = melt(phenotype,"Maturation_size") # make the dataframe
      #name = paste("Maturation size of species ",i, sep = "")
      
      #short cut the data frame when species does not reach end of simulation
      if (sum(which(Phen$value == TimeMax)) == 0) # if there is at least one value equal to the end of the sim it means that the species survived until then
        stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
      
      name = paste("Species",i, sep = " ")
      
      # prepare the data for the ribbon
      g1 <- ggplot(stat)+
        geom_smooth(aes(x = time, y = percentMean-percentSd)) +
        geom_smooth(aes(x= time, y = percentMean+percentSd)) 
      
      gg1 <- ggplot_build(g1)
      dfRibbon <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y) #and extract the smooth data for the ribbon
      
      if (!Normalisation) stat$percentMean = stat$mean
      
      if(length(SpIdx)>1)
      {
        if (i == SpIdx[1]) # top of the multi plot
        {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean)) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(labels = NULL, name = NULL, limits  = c(0,TimeMax)) +
            #scale_y_continuous(name = name, limits = c(-0.3,0.3), breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)) + #how can seq() be so bad at his only job
            scale_y_continuous(name = name, limits = window) + #, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            
            ggtitle(NULL)
        }
        else if (i == SpIdx[length(SpIdx)]) # bottom of multiplots
          
        {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean)) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(name = "Time (yr)", limits  = c(0,TimeMax)) +
            scale_y_continuous(name = name, limits = window) + #, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle(NULL)
        }
        
        else {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean)) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(labels = NULL, name = NULL, limits  = c(0,TimeMax)) +
            scale_y_continuous(name = name, limits = window) + #, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle(NULL)
        }
      }
      
      if (length(SpIdx) == 1)
      {
        p = ggplot() +
          geom_smooth(data = stat, aes(x = time, y = percentMean)) +
          geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.4)+
          geom_hline(yintercept = 0, linetype = "dashed") +
          theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
          scale_x_continuous(name = "Time (yr)", limits  = c(0,TimeMax)) +
          scale_y_continuous(name = name, limits = window) + #, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
          ggtitle("Maturation size (g)")
      }
      
      plotMatStore[[i]] = p
      
      
    }
  }
  # PPMR 
  if (PPMR)
  {
    abundanceT = sweep(abundanceNormal,2,SumPar$beta,"*") # I use the normalised abundance and multiply by the trait value
    
    # Calculate mean at each time step
    TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
    names(dimnames(TotMean)) = list("time","species")
    
    for (i in SpIdx)
    {
      AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
      if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
      TotMean[,i] = AMean
    }
    
    # Calculate variance and standard deviation
    # it is the sum of the difference between value and mean squared and multiplied by the weight
    
    statPPMR = list() # list with the stats of all species
    for (i in SpIdx)
    {
      meanSp = TotMean[,i] # take the mean of the species
      traitSp = SumPar$beta[SumPar$species == i] # take the traits of the ecotypes in the species
      weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
      stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 5,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","percentMean","percentSd"))) # initialise the matrix
      for (j in 1:length(meanSp)) # for each time step
      {
        variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
        stat[j,3] = sqrt(variance) # calculate the standard deviation
        stat[j,4] = (meanSp[j] - SumPar$beta[i])/SumPar$beta[i] # normalisation of mean
        stat[j,5] = stat[j,3]/SumPar$beta[i] # normalised sd
      }
      statPPMR[[i]] = as.data.frame(stat) # put in the list
    }
    
    # I have the stats for every species, just need to plot now
    
    
    for (i in SpIdx)
    {
      stat = statPPMR[[i]] # take the stats of the right species
      
      phenotype = TT[TT$Species == i,c(3,4,6)] # recup the traits time values
      Phen = melt(phenotype,"PPMR") # make the dataframe
      #name = paste("PPMR of species ",i, sep = "")
      
      #short cut the data frame when species does not reach end of simulation
      if (sum(which(Phen$value == TimeMax)) == 0) # if there is at least one value equal to the end of the sime it means that the species survived until then
        stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
      
      # prepare the data for the ribbon
      g1 <- ggplot(stat)+
        geom_smooth(aes(x = time, y = percentMean-percentSd)) +
        geom_smooth(aes(x= time, y = percentMean+percentSd)) 
      
      gg1 <- ggplot_build(g1)
      dfRibbon <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y) #and extract the smooth data for the ribbon
      
      if(length(SpIdx)>1)
      {
        if (i == SpIdx[1]) # top of the multi plot
        {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean, colour ="#000000")) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(labels = NULL, name = NULL, limits  = c(0,TimeMax)) +
            scale_y_continuous(name = NULL, limits = window, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle("b)")
          
        }
        else if (i == SpIdx[length(SpIdx)]) # bottom of multiplots
          
        {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean, colour ="#000000")) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(name = "Time (yr)", limits  = c(0,TimeMax)) +
            scale_y_continuous(name = NULL, limits = window, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle(NULL)
        }
        
        else {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean, colour ="#000000")) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(labels = NULL, name = NULL, limits  = c(0,TimeMax)) +
            scale_y_continuous(name = NULL, limits = window, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle(NULL)
        }
      }
      if (length(SpIdx) == 1)
      {
        p = ggplot() +
          # geom_point(data = Phen, aes(x = value, y = PPMR, group = PPMR)) +
          # geom_line(data = Phen, aes(x = value, y = PPMR, group = PPMR)) + 
          geom_smooth(data = stat, aes(x = time, y = mean, color =" red")) +
          geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                      fill = "grey", alpha = 0.4)+
          geom_hline(yintercept = stat$mean[1], linetype = "dashed") +
          scale_x_continuous(name = "Time (yr)", limits  = c(0,TimeMax)) +
          scale_y_continuous(name = NULL) +
          scale_colour_discrete(labels = c("mean","standard deviation"))+
          theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          guides(color=guide_legend(override.aes=list(fill=NA)))+
          ggtitle("Preferred PPMR")
      }
      
      plotPPMRStore[[i]] <- p
      
      
      
    }
  }
  # WIDTH FEEDING KERNEL
  if (Sig)
  {
    abundanceT = sweep(abundanceNormal,2,SumPar$sigma,"*") # I use the normalised abundance and multiply by the trait value
    
    # Calculate mean at each time step
    TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
    names(dimnames(TotMean)) = list("time","species")
    
    for (i in SpIdx)
    {
      AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
      if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
      TotMean[,i] = AMean
    }
    
    # Calculate variance and standard deviation
    # it is the sum of the difference between value and mean squared and multiplied by the weight
    
    statDB = list() # list with the stats of all species
    for (i in SpIdx)
    {
      meanSp = TotMean[,i] # take the mean of the species
      traitSp = SumPar$sigma[SumPar$species == i] # take the traits of the ecotypes in the species
      weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
      stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 5,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","percentMean","percentSd"))) # initialise the matrix
      for (j in 1:length(meanSp)) # for each time step
      {
        variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
        stat[j,3] = sqrt(variance) # calculate the standard deviation
        stat[j,4] = (meanSp[j] - SumPar$sigma[i])/SumPar$sigma[i] # normalisation of mean
        stat[j,5] = stat[j,3]/SumPar$sigma[i] # normalised sd
      }
      statDB[[i]] = as.data.frame(stat) # put in the list
      #statDB[[i]] = stat# put in the list
    }
    
    
    
    # I have the stats for every species, just need to plot now
    
    for (i in SpIdx)
    {
      stat = statDB[[i]] # take the stats of the right species
      stat$up = stat$mean + stat$sd
      stat$low = stat$mean - stat$sd
      
      phenotype =  TT[TT$Species == i,c(3,4,7)] # recup the traits time values
      Phen = melt(phenotype,"Diet_breadth") # make the dataframe
      #name = paste("Diet breadth of species ",i, sep = "")
      
      #short cut the data frame when species does not reach end of simulation
      if (sum(which(Phen$value == TimeMax)) == 0) # if there is at least one value equal to the end of the sime it means that the species survived until then
        stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
      
      # prepare the data for the ribbon
      g1 <- ggplot(stat)+
        geom_smooth(aes(x = time, y = percentMean-percentSd)) +
        geom_smooth(aes(x= time, y = percentMean+percentSd)) 
      
      gg1 <- ggplot_build(g1)
      dfRibbon <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y) #and extract the smooth data for the ribbon
      
      if(length(SpIdx)>1)
      {
        if (i == SpIdx[1]) # top of the multi plot
        {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean, colour ="#000000")) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(labels = NULL, name = NULL, limits  = c(0,TimeMax)) +
            scale_y_continuous(name = NULL, limits = window, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle("c)")
          
        }
        else if (i == SpIdx[length(SpIdx)]) # bottom of multiplots
          
        {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean, colour ="#000000")) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(name = "Time (yr)", limits  = c(0,TimeMax)) +
            scale_y_continuous(name = NULL, limits = window, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle(NULL)
        }
        
        else {
          p = ggplot() +
            geom_smooth(data = stat, aes(x = time, y = percentMean, colour ="#000000")) +
            geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey", alpha = 0.4)+
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                  legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
            scale_x_continuous(labels = NULL, name = NULL, limits  = c(0,TimeMax)) +
            scale_y_continuous(name = NULL, limits = window, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
            ggtitle(NULL)
        }
      }
      if (length(SpIdx) == 1)
      {
        p = ggplot() +
          # geom_point(data = Phen, aes(x = value, y = Diet_breadth, group = Diet_breadth)) +
          # geom_line(data = Phen, aes(x = value, y = Diet_breadth, group = Diet_breadth)) + 
          geom_smooth(data = stat, aes(x = time, y = mean, color =" red")) +
          geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax),
                      fill = "grey", alpha = 0.4)+
          geom_hline(yintercept = stat$mean[1], linetype = "dashed") +
          scale_x_continuous(name = "Time (yr)", limits  = c(0,TimeMax)) +
          scale_y_continuous(name = NULL) +
          scale_colour_discrete(labels = c("mean","standard deviation"))+
          theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          guides(color=guide_legend(override.aes=list(fill=NA)))+
          ggtitle("Width feeding kernel")
      }
      
      plotDBStore[[i]] <- p
      
    }
  }
  
  if (comments) cat("plotTraitsMulti ends\n") 
  
  if (returnData) return(list(plotMatStore,plotPPMRStore,plotDBStore))
  
  
  # construction of the multiplot
  if (PPMR == F & Sig == F)
  {
    rowPos = rep(c(1,2,3),3)
    colPos = c(rep(1,3),rep(2,3),rep(3,3))
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=3, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(10,3), "cm"))))
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = rowPos[i], layout.pos.col = colPos[i])) 
    }
    
  } else {
    
    p <- multiplot(plotMatStore[[1]],plotMatStore[[2]], plotMatStore[[3]],
                   plotMatStore[[4]],plotMatStore[[5]],
                   plotMatStore[[6]],plotMatStore[[7]],plotMatStore[[8]],plotMatStore[[9]],
                   plotPPMRStore[[1]],plotPPMRStore[[2]], plotPPMRStore[[3]],
                   plotPPMRStore[[4]], plotPPMRStore[[5]],
                   plotPPMRStore[[6]],plotPPMRStore[[7]],plotPPMRStore[[8]],plotPPMRStore[[9]],
                   plotDBStore[[1]],plotDBStore[[2]],plotDBStore[[3]],
                   plotDBStore[[4]], plotDBStore[[5]],
                   plotDBStore[[6]],plotDBStore[[7]],plotDBStore[[8]],plotDBStore[[9]],
                   cols=3)
  }
  
  return(p)
  
}

#plot the mean of every runs and differentiate normal and fished
plotTraitOverlap <- function(directory, SpIdx = NULL, comments = T, Mat = T, PPMR = T, Sig = T,window = NULL, init =F)
{
  tic()
  if (is.null(SpIdx)) SpIdx = seq(1,9)
  
  normalData = F # to determine what data I have
  fisheriesData = F
  
 # get the initial wmat and Time max values
  
    if(dir.exists(file.path(paste(directory,"/init",sep=""))))
    {
      WmatSim <- get(load(paste(directory,"/init/run1/run.Rdata",sep="")))
      Wmat = WmatSim@params@species_params$w_mat[1:SpIdx[length(SpIdx)]]
      TimeMax = WmatSim@params@species_params$timeMax[1]
      if (comments) {
        cat("Wmat is")
        print(Wmat)}
    } else {cat("Could not find Wmat not TimeMax values, taking the default ones\n")}
  

  if(init) #if I want to plot the initialisation part
  {
    if(dir.exists(file.path(paste(directory,"/init",sep=""))))
    {
      if (comments) cat("Using initialisation simulations\n")
      normalData = T
      normalList <- bunchLoad(folder = paste(directory,"/init",sep=""))
      # get the plots
      normalTraitsList <- list()
      for (i in 1:length(normalList))
      {
        if (comments) cat(sprintf("Using run %g\n",i))
        normalTraitsList[[i]] <- plotTraitsMulti(object = normalList[[i]],PPMR = F,Sig = F,returnData = T, SpIdx = SpIdx, Wmat = Wmat)
      }
      if (comments) cat("Plots loaded")
      speciesListN <- list()
      for (j in SpIdx) # for every species
      {
        speciesData <- NULL
        for (i in 1:length(normalTraitsList)) # at each time
        {
          if (!is.null(normalTraitsList[[i]][[1]][[j]]))
          {
            a <- ggplot_build(normalTraitsList[[i]][[1]][[j]]) # take the plot data
            if (!is.null(a$data[[1]]$group))
            {
              a$data[[1]]$group <- i # change the group to the run number
              speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
            }
          }
        }
        speciesListN[[j]] <- speciesData #this is a list of all the species at each time
      }
    } 
  } else {
  # NO FISHERIES PART
  # need to load the data first
  if(dir.exists(file.path(paste(directory,"/normal",sep=""))))
  {
    if (comments) cat("Using simulations without fisheries\n")
    normalData = T
    normalList <- bunchLoad(folder = paste(directory,"/normal",sep=""))
    # get the plots
    normalTraitsList <- list()
    for (i in 1:length(normalList))
    {
      if (comments) cat(sprintf("Using run %g\n",i))
      normalTraitsList[[i]] <- plotTraitsMulti(object = normalList[[i]],PPMR = F,Sig = F,returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = TimeMax)
    }
    if (comments) cat("Plots loaded")
    speciesListN <- list()
    for (j in SpIdx) # for every species
    {
      speciesData <- NULL
      for (i in 1:length(normalTraitsList)) # at each time
      {
        if (!is.null(normalTraitsList[[i]][[1]][[j]]))
        {
          a <- ggplot_build(normalTraitsList[[i]][[1]][[j]]) # take the plot data
          if (!is.null(a$data[[1]]$group))
          {
            a$data[[1]]$group <- i # change the group to the run number
            speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
          }
        }
      }
      speciesListN[[j]] <- speciesData #this is a list of all the species at each time
    }
  }
  # FISHERIES PART
  # need to load the data first
  if(dir.exists(file.path(paste(directory,"/fisheries",sep=""))))
  {
    if (comments) cat("Using simulation with fisheries\n")
    fisheriesData = T
    fishList <- bunchLoad(folder = paste(directory,"/fisheries",sep=""))
    # get the plots
    fishTraitsList <- list()
    for (i in 1:length(fishList)){
      if (comments) cat(sprintf("Using run %g\n",i))
      fishTraitsList[[i]] <- plotTraitsMulti(object = fishList[[i]],PPMR = F,Sig = F,returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = TimeMax)
    }
    if (comments) cat("Plots loaded")
    speciesListF <- list()
    for (j in SpIdx) # for every species
    {
      speciesData <- NULL
      for (i in 1:length(fishTraitsList)) # at each time
      {
        if (!is.null(fishTraitsList[[i]][[1]][[j]]))
        {
          a <- ggplot_build(fishTraitsList[[i]][[1]][[j]]) # take the plot data
          if (!is.null(a$data[[1]]$group))
          {
            a$data[[1]]$group <- i # change the group to the run number
            speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
          }
        }
      }
      speciesListF[[j]] <- speciesData #this is a list of all the species at each time
    }
  }
  }
  # plot time
  plotMatStore = list()
  plotSigStore = list()
  plotPPMStore = list()
  
  if (is.null(window)) window = c(-1,1)
  if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
  
  # what data do I have ? Only normal, only fisheries or both?
  
  if (normalData && fisheriesData) # if I have both
  {
    temp = NULL
    for (i in SpIdx) if (!is.null(speciesListN[[i]])) temp = c(temp,i) #update SpIdx if species are extinct at all times
    SpIdx = temp
    
    for (i in SpIdx) # do the plots for each species and traits
    {
      if (comments) cat(sprintf("plot %i\n",i))
      # if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
      title = NULL
      
      if (Mat){
        plotMatStore[[i]] <-ggplot(speciesListN[[i]])+
          geom_line(aes(x=x,y=y, group = group)) +#colour=as.factor(group))) +
          geom_line(data = speciesListF[[i]], aes(x=x,y=y, group = group), linetype = "dashed", color = "red") +
          scale_x_continuous(name = NULL)+
          scale_y_continuous(name = sprintf("Species %i",i), limits = window)+
          scale_color_grey(name = "Species")+ # grey
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[1])
      }
      if (Sig)
      {
        plotSigStore[[i]] <-ggplot(speciesList[[i]])+
          geom_point(aes(x=x,y=y,colour=as.factor(group))) +
          geom_line(data =trendList[[i]], aes(x=x,y=y,group = run))+
          scale_x_continuous(name = NULL, limits = window)+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[2])
      }
      if(PPMR)
      {
        plotPPMStore[[i]] <-ggplot(speciesList[[i]])+
          geom_point(aes(x=x,y=y,colour=as.factor(group))) +
          scale_x_continuous(name = NULL, limits = window)+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[3])
      }
    }
  } else { # if I have only one of them
    if (fisheriesData) speciesListN = speciesListF # if I have only fisheries, I convert everything to normal as it is the same
    temp = NULL
    for (i in SpIdx) if (!is.null(speciesListN[[i]])) temp = c(temp,i) #update SpIdx if species are extinct at all times
    SpIdx = temp
    
    for (i in SpIdx) # do the plots for each species and traits
    {
      if (comments) cat(sprintf("plot %i\n",i))
      #if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
      title = NULL
      
      if (Mat){
        plotMatStore[[i]] <-ggplot(speciesListN[[i]])+
          geom_line(aes(x=x,y=y, group = group)) +#colour=as.factor(group))) +
          scale_x_continuous(name = NULL)+
          scale_y_continuous(name = sprintf("Species %i",i), limits = window)+
          scale_color_grey(name = "Species")+ # grey
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[1])
      }
      if (Sig)
      {
        plotSigStore[[i]] <-ggplot(speciesList[[i]])+
          geom_point(aes(x=x,y=y,colour=as.factor(group))) +
          geom_line(data =trendList[[i]], aes(x=x,y=y,group = run))+
          scale_x_continuous(name = NULL, limits = window)+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[2])
      }
      if(PPMR)
      {
        plotPPMStore[[i]] <-ggplot(speciesList[[i]])+
          geom_point(aes(x=x,y=y,colour=as.factor(group))) +
          scale_x_continuous(name = NULL, limits = window)+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[3])
      }
    }
  }
  
  
  
  if (comments) cat("plot done, starting the multiplot\n")
  # do the multiplot
   if (init) path_to_png = paste(directory,"/TraitInit.png",sep="") else path_to_png = paste(directory,"/Trait.png",sep="")
  
  if (PPMR == F & Sig == F)
  {
    gridID = matrix(c(1,rep(2,3),rep(3,5),rep(1,2),rep(2,4),rep(3,3)), nrow = 2, ncol = 9,byrow = T, dimnames = list(c("nRow","nCol"),seq(1,9))) # grid dimension depending on sp number
    rowPos = rep(1:gridID[1,length(SpIdx)],gridID[2,length(SpIdx)])
    colPos = NULL
    for (i in 1:gridID[2,length(SpIdx)]) colPos = c(colPos,rep(i,gridID[1,length(SpIdx)]))

    noRow = rowPos[length(rowPos)] # last value to determine number of row/col
    noCol = colPos[length(colPos)]
    
    png(filename=path_to_png, width = 10*noCol, height = 10*noRow, units = "cm",res = 1000)
    
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=noRow, ncol=noCol, 
                                               widths = unit(rep(10,noCol), "cm"), 
                                               heights = unit(rep(10,noRow), "cm"))))

    
    
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = rowPos[i], layout.pos.col = colPos[i])) 
    }
    
  } else {
    
    png(filename=path_to_png,width = 30, height = 45, units = "cm",res = 1000)
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=9, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(5,9), "cm"))))
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)) 
      print(plotSigStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 2)) 
      print(plotPPMStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 3)) 
    }
  }
  dev.off()
  if(comments) toc()
  
}

# biomass avg + trait below, species overlapping, this one to compare with and without predation
plotBiomAndTraits <- function(folder,Mat = T, PPMR = F,Sig = F, SpIdx = NULL, file = NULL, comments = T, window = NULL, what = c("no predation","predation"))
{
  tic()
  
  if(is.null(file)) file = "/normal"
  path_to_png = paste(folder,"/Biom&Traits.png",sep="")
  
  if (file == "/fisheries") path_to_png = paste(folder,"/Biom&TraitsFisheries.png",sep="")
  
  # get the initial stuff
  dt = 0.1
  if(dir.exists(file.path(paste(folder,"/noPred/init",sep=""))))
  {
    WmatSim <- get(load(paste(folder,"/noPred/init/run1/run.Rdata",sep=""))) # get one of the initial sims
    if (is.null(SpIdx)) SpIdx = unique(WmatSim@params@species_params$species) # determine spidx if not already given by user
    no_sp = length(SpIdx) # get the species number from this
    Wmat = WmatSim@params@species_params$w_mat[1:SpIdx[length(SpIdx)]] # get the mat size values
    TimeMax = WmatSim@params@species_params$timeMax[1] # and the length of the initialisation
    if (comments) cat("Initialisation done\n")
  } else {cat("Could not find initial values, plot is going to crash in 3 ... 2 ... 1...\n")}
  
  plotList <- list() # to store the plots
  for (column in what) # are we plotting pred or no pred sims ?
  {
    if (column == "no predation") 
      {
      directory <- paste(folder,"/noPred",sep = "")
      listPosition = 1
      title = "No interactions"
    } else if (column == "predation") 
      {
      directory <- paste(folder,"/Pred",sep = "")
      listPosition = 2
      title = "Predation interactions"
    } else stop("I don't know what to plot")

  multiSim <- bunchLoad(folder = paste(directory,file,sep="")) # load the folder (either normal or fisheries)
  
  if(comments) cat("sim loaded\n")
  
  # color gradient
  colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
  colGrad <- colfunc(no_sp)
  names(colGrad) <- seq(1,no_sp)

  
  # BIOMASS
  sim <- superStich(multiSim) # stich the abundance together

  plot_datBiom <- biom(multiSim[[1]],phenotype = T) # get the species and phenotypes abundance for one sim
  plot_datBiom$color
  
  
  p1 <- ggplot(plot_datBiom[[1]]) +
    geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1) +
    geom_line(data = plot_datBiom[[2]], aes(x = time, y = value, colour = as.factor(bloodline), group = sp), alpha = 0.2) +
    scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(1e-15, NA)) +
    scale_x_continuous(name = NULL, labels = NULL) +
    scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    ggtitle(title)

  plot_datBiom <- biom(sim,phenotype = F) # get the abundance data for all the sitched sim
  p2 <- ggplot(plot_datBiom[[1]]) +
    geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1) +
    scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(NA, NA)) + 
    scale_x_continuous(name = NULL) +
    scale_color_manual(name = "Species", values =colGrad)+ # color gradient
    theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    ggtitle(NULL)

  if(comments) cat("Biomass done\n")
  
  #TRAITS
   # get the plots
  TraitsList <- list()
  for (i in 1:length(multiSim))
  {
    if (comments) cat(sprintf("Using run %g\n",i))
    TraitsList[[i]] <- plotTraitsMulti(object = multiSim[[i]],Mat = Mat, PPMR = PPMR,Sig = Sig, returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = TimeMax, window = window)
  }
  if (comments) cat("Plots loaded")
  speciesList <- list()
  for (j in SpIdx) # for every species
  {
    speciesData <- NULL
    for (i in 1:length(TraitsList)) # at each time
    {
      if (!is.null(TraitsList[[i]][[1]][[j]]))
      {
        a <- ggplot_build(TraitsList[[i]][[1]][[j]]) # take the plot data
        if (!is.null(a$data[[1]]$group))
        {
          a$data[[1]]$species <- j # add species identity in the data
          a$data[[1]]$group <- 100*j + i # change the group to the run number (up to 99 different runs to keep species identity as hundreds)
          
          speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
        }
      }
    }
    speciesList[[j]] <- speciesData #this is a list of all the species at each time
  }
  
  # let's add the initial trait run
  multiSimInit <- bunchLoad(folder = paste(directory,"/init",sep="")) # load the init sims
  
  TraitsListInit <- list()
  for (i in 1:length(multiSimInit))
  {
    if (comments) cat(sprintf("Using run %g\n",i))
    TraitsListInit[[i]] <- plotTraitsMulti(object = multiSimInit[[i]],Mat = Mat, PPMR = PPMR,Sig = Sig, returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = 0, window = window)
  }
  if (comments) cat("Init plots loaded")
  speciesListInit <- list()
  for (j in SpIdx) # for every species
  {
    speciesDataInit <- NULL
    for (i in 1:length(TraitsListInit)) # at each time
    {
      if (!is.null(TraitsListInit[[i]][[1]][[j]]))
      {
        a <- ggplot_build(TraitsListInit[[i]][[1]][[j]]) # take the plot data
        if (!is.null(a$data[[1]]$group))
        {
          a$data[[1]]$species <- j # add species identity in the data
          a$data[[1]]$group <- 100*j + i # change the group to the run number (up to 99 different runs to keep species identity as hundreds)
          
          speciesDataInit <- rbind(speciesDataInit,a$data[[1]]) #bind the same species at different time
        }
      }
    }
    speciesListInit[[j]] <- speciesDataInit #this is a list of all the species at each time
  }
  
  if (is.null(window)) window = c(-1,1)
  if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
  plot_datTrait <- do.call("rbind",speciesList)
  plot_datTraitInit <- do.call("rbind",speciesListInit)
  
  # this is not even coding anymore, adding the time shift to the sim after init
  plot_datTrait$x <- plot_datTrait$x + TimeMax * dt
  
  plot_datMiracle <- rbind(plot_datTrait,plot_datTraitInit)
  
  p3 <- ggplot(plot_datMiracle)+
    geom_line(aes(x=x,y=y, group = group, color = as.factor(species))) +
    scale_x_continuous(name = "Time in years")+
    scale_y_continuous(name = "Trait relative proportion to initial value", limits = window)+
    scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom", 
          #legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"))+ 
    guides(color = guide_legend(nrow=1)) +
    ggtitle(NULL)
  
  plotList[[listPosition]] <- list(p1,p2,p3)
  }
  plotList <- plotList[lapply(plotList,length)>0] # if a slot is empty
  
  # MULTIPLOT

  # what do we get to plot?
  rowPos = rep(seq(1,3),2)  
  colPos = c(rep(1,3),rep(2,3))
  if (length(plotList) == 2) # if I have both plot list it means it is a 2 col plot
  {
  png(filename=path_to_png, width = 20, height = 30, units = "cm",res = 600)

  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(nrow=3, ncol=2, 
                                             widths = unit(10, "cm"), 
                                             heights = unit(10, "cm"))))
  for (i in rowPos)
    for (j in colPos)
    print(plotList[[j]][[i]], vp = viewport(layout.pos.row = i, layout.pos.col = j))
  
  # print(plotList[[1]][[1]], vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) 
  # print(plotList[[1]][[2]], vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) 
  # print(plotList[[1]][[3]], vp = viewport(layout.pos.row = 3, layout.pos.col = 1)) 
  # print(plotList[[2]][[1]], vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) 
  # print(plotList[[2]][[2]], vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) 
  # print(plotList[[2]][[3]], vp = viewport(layout.pos.row = 3, layout.pos.col = 2)) 
  
  } else {
    
    png(filename=path_to_png, width = 10, height = 30, units = "cm",res = 600)

    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=3, ncol=1, 
                                               widths = unit(10, "cm"), 
                                               heights = unit(10, "cm"))))
    
    for (i in rowPos)
      print(plotList[[1]][[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  
  dev.off()
  toc()
}

#this one to compare different sp number
plotDiffSp <- function(folder,Mat = T, PPMR = F,Sig = F, SpIdx = NULL, file = NULL, comments = T, window = NULL, what = c("3sp","9sp"), returnData = F)
{
  tic()
  
  if(is.null(file)) file = "/init"
  path_to_png = paste(folder,"/Biom&TraitsCompareSpecies.png",sep="")
  
  #if (file == "/fisheries") path_to_png = paste(folder,"/Biom&TraitsFisheries.png",sep="")
  
  plotList <- list() # to store the plots
  buildList <- list()
  dataList <- list() # to return the data
  for (column in what) # are we plotting pred or no pred sims ?
  {
    columnList <- list() # to store the data per column
    if (column == "3sp") 
    {
      directory <- paste(folder,"/3sp",sep = "")
      listPosition = 1
      title = c("(a)","(b)","(c)")
    } else if (column == "9sp") 
    {
      directory <- paste(folder,"/9sp",sep = "")
      listPosition = 2
      title = c("(d)","(e)","(f)")
    } else stop("I don't know what to plot")
    
    multiSim <- bunchLoad(folder = paste(directory,file,sep="")) # load the folder (either normal or fisheries)
    if(comments) cat("sim loaded\n")
    # get the initial stuff
    dt = 0.1

      SpIdx = sort(unique(multiSim[[1]]@params@species_params$species)) # determine spidx if not already given by user
      no_sp = length(SpIdx) # get the species number from this
      Wmat = multiSim[[1]]@params@species_params$w_mat[1:SpIdx[length(SpIdx)]] # get the mat size values
      TimeMax = multiSim[[1]]@params@species_params$timeMax[1] # and the length of the initialisation
      if (comments) cat("Initialisation done\n")

    # color gradient
    colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
    colGrad <- colfunc(no_sp)
    names(colGrad) <- seq(1,no_sp)

    # BIOMASS
    sim <- superStich(multiSim) # stich the abundance together

    plot_datBiom1 <- biom(multiSim[[1]],phenotype = T) # get the species and phenotypes abundance for one sim
    if (returnData) columnList[[1]] <- plot_datBiom1
      
    p1 <- ggplot(plot_datBiom1[[1]]) +
      geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1) +
      geom_line(data = plot_datBiom1[[2]], aes(x = time, y = value, colour = as.factor(bloodline), group = sp), alpha = 0.2) +
      scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(1e-15, NA)) +
      scale_x_continuous(name = NULL, labels = NULL) +
      scale_color_manual(name = "Species", values = colGrad)+ # color gradient
      theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
            legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
      ggtitle(title[1])

    plot_datBiom2 <- biom(sim,phenotype = F) # get the abundance data for all the sitched sim
    if (returnData) columnList[[2]] <- plot_datBiom2
    
    p2 <- ggplot(plot_datBiom2[[1]]) +
      geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1) +
      scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(NA, NA)) +
      scale_x_continuous(name = NULL, labels = NULL) +
      scale_color_manual(name = "Species", values =colGrad)+ # color gradient
      theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
            legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
      ggtitle(title[2])

    if(comments) cat("Biomass done\n")
    
    #TRAITS
    # get the plots
    TraitsList <- list()
    for (i in 1:length(multiSim))
    {
      if (comments) cat(sprintf("Using run %g\n",i))
      TraitsList[[i]] <- plotTraitsMulti(object = multiSim[[i]],Mat = Mat, PPMR = PPMR,Sig = Sig, returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = TimeMax, window = window)
    }
    if (comments) cat("Plots loaded")
    speciesList <- list()
    for (j in SpIdx) # for every species
    {
      speciesData <- NULL
      for (i in 1:length(TraitsList)) # at each time
      {
        if (!is.null(TraitsList[[i]][[1]][[j]]))
        {
          a <- ggplot_build(TraitsList[[i]][[1]][[j]]) # take the plot data
          if (!is.null(a$data[[1]]$group))
          {
            a$data[[1]]$species <- j # add species identity in the data
            a$data[[1]]$group <- 100*j + i # change the group to the run number (up to 99 different runs to keep species identity as hundreds)
            
            speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
          }
        }
      }
      speciesList[[j]] <- speciesData #this is a list of all the species at each time
    }
    
    if (is.null(window)) window = c(-1,1)
    if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
    
    plot_datTrait <- do.call("rbind",speciesList)
    if (returnData) columnList[[3]] <- plot_datTrait

    p3 <- ggplot(plot_datTrait)+
      geom_line(aes(x=x,y=y, group = group, color = as.factor(species))) +
      scale_x_continuous(name = "Time in years")+
      scale_y_continuous(name = "Trait relative proportion to initial value", limits = window)+
      scale_color_manual(name = "Species", values = colGrad)+ # color gradient
      theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom", 
            #legend.justification=c(1,1),
            legend.key = element_rect(fill = "white"))+ 
      guides(color = guide_legend(nrow=1)) +
      ggtitle(title[3])
    
    plotList[[listPosition]] <- list(p1,p2,p3)
    buildList[[listPosition]] <- list(ggplot_build(p1),ggplot_build(p2),ggplot_build(p3)) # the color gradient changes for the first loop if I don't do this. Why? I have no fucking idea
    dataList[[listPosition]] <- columnList
  }
  if (returnData) return(dataList)
  plotList <- plotList[lapply(plotList,length)>0] # if a slot is empty

  # MULTIPLOT
  
  # what do we get to plot?
  rowPos = rep(seq(1,3),2)  
  colPos = c(rep(1,3),rep(2,3))
  if (length(plotList) == 2) # if I have both plot list it means it is a 2 col plot
  {
    png(filename=path_to_png, width = 24, height = 30, units = "cm",res = 600)
    
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=3, ncol=2, 
                                               widths = unit(12, "cm"), 
                                               heights = unit(10, "cm"))))
    for (i in rowPos)
      for (j in colPos)
        print(plotList[[j]][[i]], vp = viewport(layout.pos.row = i, layout.pos.col = j))

  } else {
    
    png(filename=path_to_png, width = 12, height = 30, units = "cm",res = 600)
    
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=3, ncol=1, 
                                               widths = unit(12, "cm"), 
                                               heights = unit(10, "cm"))))
    
    for (i in rowPos)
      print(plotList[[1]][[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  
  dev.off()
  toc()
}

#this one when fished (so I have to load init and fisheries and maybe normal)
plotDiffSpFished <- function(folder,Mat = T, PPMR = F,Sig = F, SpIdx = NULL, comments = T, window = NULL, what = c("normal","fisheries"), option = "A")
{
  tic()
  dt = 0.1
  path_to_png = paste(folder,"/Biom&TraitsCompareFished.png",sep="")
  
  plotList <- list() # to store the plots
   multiSimInit <- bunchLoad(folder = paste(folder,"/init",sep=""))
  
  if(comments) cat("Init sims loaded\n")
  
  for (column in what) # are we plotting normal or fisheries sims ?
  {
    if (column == "normal") 
    {
      listPosition = 1
      multiSim <- bunchLoad(folder = paste(folder,"/normal",sep=""))
      title = c("(a) Without fisheries","(b)","(c)")
    } else if (column == "fisheries") 
    {
      listPosition = 2
      multiSim <- bunchLoad(folder = paste(folder,"/fisheries",sep=""))
      title = c("(d) With fisheries","(e)","(f)")
    } else stop("I don't know what to plot")
    
    if (option == "A")
    {
    #multiSimNormal <- bunchLoad(folder = paste(folder,"/normal",sep=""))
    template = multiSimInit[[1]]
    longSimList <- list()
    #prepare the data 
    for (x in 1:length(multiSimInit))
    {
      if(comments)  cat(sprintf("Using run %i\n",x))
      # get the right species params
      SummaryParams = template@params # structure and basics params
      SummaryParams@species_params = rbind(multiSimInit[[x]]@params@species_params,multiSim[[x]]@params@species_params) # get the sp ID from both sim
      SummaryParams@species_params$timeMax = sum(unique(SummaryParams@species_params$timeMax)) # update the timemax 
      a <- SummaryParams@species_params
      a <- a[order(a$ecotype, a$extinct, decreasing=TRUE),] # weird 3 lines to get rid of duplicates and keep the ones with the extinction value
      a <- a[!duplicated(a$ecotype),]
      SummaryParams@species_params = a[order(a$pop,a$ecotype),]
      
      if (comments) cat("Data handling\n")
      
      sim = list(multiSimInit[[x]],multiSim[[x]]) # get the data
      
      # stitiching the sims together
      Dtime = SummaryParams@species_params$timeMax[1] * dt
      Dsp = length(SummaryParams@species_params$ecotype)
      Dw = dim(sim[[1]]@n)[3]
      # put all the sim at the same dimension
      biomList <- list()
      for (i in 1:length(sim)) # for each sim
      {
        biom <- array(data = 0, dim = c(dim(sim[[i]]@n)[1],Dsp,Dw), dimnames = list(dimnames(sim[[i]]@n)$time, SummaryParams@species_params$ecotype, SummaryParams@w)) # create an array of the right dimension
        names(dimnames(biom)) = c("time","species","size")
        for (j in dimnames(sim[[i]]@n)$species) # fill it when necessary
          biom[,which(dimnames(biom)$species == j),] = sim[[i]]@n[,which(dimnames(sim[[i]]@n)$species == j),]
        
        biomList[[i]] <- biom # store it
      }
      biom <- do.call("abind",list(biomList, along = 1)) # abind the list
      names(dimnames(biom)) = list("time","species","size")
      dimnames(biom)$time = seq(1, SummaryParams@species_params$timeMax[1]* dt)  
      
      # I have to do the phyto aussi
      phyto <- do.call(abind, c(lapply(sim, function(isim) isim@n_pp),along = 1))
      
      # taking care of the effort
      effort <- do.call(rbind, lapply(sim, function(isim) isim@effort))
      names(dimnames(effort)) = list("time","effort")
      dimnames(effort)$time = seq(1, SummaryParams@species_params$timeMax[1]* dt) 
      
      # reconstruct the mizer object
      sim = template
      sim@params=SummaryParams
      sim@n = biom
      sim@effort = effort
      sim@n_pp = phyto
      
      longSimList[[x]] <- sim
      rm(list = "biom","phyto","effort","sim")
      gc()
      
    }
    }
    
    if (comments) cat("Data ready")  
    # get the initial stuff
    SpIdx = sort(unique(longSimList[[1]]@params@species_params$species)) # determine spidx if not already given by user
    no_sp = length(SpIdx) # get the species number from this
    Wmat = longSimList[[1]]@params@species_params$w_mat[1:SpIdx[length(SpIdx)]] # get the mat size values
    TimeMax = longSimList[[1]]@params@species_params$timeMax[1] # and the length of the initialisation
    if (comments) cat("Initialisation done\n")
    
    
    # color gradient
    colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
    colGrad <- colfunc(no_sp)
    names(colGrad) <- seq(1,no_sp)
    # colGrad[[column]] <- colfunc(no_sp)
    # names(colGrad[[column]]) <- seq(1,no_sp)
    
    
    #BIOMASS
    # sim <- superStich(longSimList) # stich the abundance together
    # 
    # plot_datBiom <- biom(longSimList[[1]],phenotype = T) # get the species and phenotypes abundance for one sim
    # p1 <- ggplot(plot_datBiom[[1]]) +
    #   geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1) +
    #   geom_line(data = plot_datBiom[[2]], aes(x = time, y = value, colour = as.factor(bloodline), group = sp), alpha = 0.2) +
    #   scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(1e-15, NA)) +
    #   geom_vline(xintercept = 4000, linetype = "dashed") +
    #   scale_x_continuous(name = NULL, labels = NULL) +
    #   scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    #   theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
    #         panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
    #         legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    #   ggtitle(title[1])
    # 
    # plot_datBiom <- biom(sim,phenotype = F) # get the abundance data for all the sitched sim
    # p2 <- ggplot(plot_datBiom[[1]]) +
    #   geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp), size = 1) +
    #   scale_y_log10(name = expression(paste("Biomass in g.m"^"-3")), limits = c(NA, NA)) +
    #   geom_vline(xintercept = 4000, linetype = "dashed") +
    #   scale_x_continuous(name = NULL) +
    #   scale_color_manual(name = "Species", values =colGrad)+ # color gradient
    #   theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
    #         panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
    #         legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
    #   ggtitle(title[2])
    # 
    # if(comments) cat("Biomass done\n")
    
    #TRAITS
    # get the plots
    TraitsList <- list()
    for (i in 1:length(longSimList))
    {
      if (comments) cat(sprintf("Using run %g\n",i))
      TraitsList[[i]] <- plotTraitsMulti(object = longSimList[[i]],Mat = Mat, PPMR = PPMR,Sig = Sig, returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = TimeMax, window = window)
    }
    if (comments) cat("Plots loaded")
    speciesList <- list()
    for (j in SpIdx) # for every species
    {
      speciesData <- NULL
      for (i in 1:length(TraitsList)) # at each time
      {
        if (!is.null(TraitsList[[i]][[1]][[j]]))
        {
          a <- ggplot_build(TraitsList[[i]][[1]][[j]]) # take the plot data
          if (!is.null(a$data[[1]]$group))
          {
            a$data[[1]]$species <- j # add species identity in the data
            a$data[[1]]$group <- 100*j + i # change the group to the run number (up to 99 different runs to keep species identity as hundreds)
            
            speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
          }
        }
      }
      speciesList[[j]] <- speciesData #this is a list of all the species at each time
    }
    
    #let's add the initial trait run
    #multiSimInit <- bunchLoad(folder = paste(directory,"/init",sep="")) # load the init sims
    
    # TraitsListInit <- list()
    # for (i in 1:length(multiSimInit))
    # {
    #   if (comments) cat(sprintf("Using run %g\n",i))
    #   TraitsListInit[[i]] <- plotTraitsMulti(object = multiSimInit[[i]],Mat = Mat, PPMR = PPMR,Sig = Sig, returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = 0, window = window)
    # }
    # if (comments) cat("Init plots loaded")
    # speciesListInit <- list()
    # for (j in SpIdx) # for every species
    # {
    #   speciesDataInit <- NULL
    #   for (i in 1:length(TraitsListInit)) # at each time
    #   {
    #     if (!is.null(TraitsListInit[[i]][[1]][[j]]))
    #     {
    #       a <- ggplot_build(TraitsListInit[[i]][[1]][[j]]) # take the plot data
    #       if (!is.null(a$data[[1]]$group))
    #       {
    #         a$data[[1]]$species <- j # add species identity in the data
    #         a$data[[1]]$group <- 100*j + i # change the group to the run number (up to 99 different runs to keep species identity as hundreds)
    # 
    #         speciesDataInit <- rbind(speciesDataInit,a$data[[1]]) #bind the same species at different time
    #       }
    #     }
    #   }
    #   speciesListInit[[j]] <- speciesDataInit #this is a list of all the species at each time
    # }
    
    if (is.null(window)) window = c(-1,1)
    if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
    
    plot_datTrait <- do.call("rbind",speciesList)
    
    # this is not even coding anymore, adding the time shift to the sim after init
    #plot_datTrait$x <- plot_datTrait$x + TimeMax * dt
    
    #plot_datMiracle <- rbind(plot_datTrait,plot_datTraitInit)
    #print(colGrad)
    p3 <- ggplot(plot_datTrait)+
      geom_line(aes(x=x,y=y, group = group, color = as.factor(species))) +
      scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
      scale_y_continuous(name = "Trait relative proportion to initial value", limits = window)+
      geom_vline(xintercept = 4000, linetype = "dashed") +
      scale_color_manual(name = "Species", values = colGrad)+ # color gradient
      theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom", 
            #legend.justification=c(1,1),
            legend.key = element_rect(fill = "white"))+ 
      guides(color = guide_legend(nrow=1)) +
      ggtitle(title[1])
    #print(p3)
    #plotList[[listPosition]] <-list(p1,p2,p3)
    # column = what[2]
    #  a <- list(p3)
    
    plotList[[listPosition]] <- p3
  }
  plotList <- plotList[lapply(plotList,length)>0] # if a slot is empty
  
  # MULTIPLOT
  
  png(filename=path_to_png, width = 24, height = 10, units = "cm",res = 600)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow=1, ncol=2,
                                             widths = unit(12, "cm"),
                                             heights = unit(10, "cm"))))
  
  print(plotList[[1]], vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plotList[[2]], vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  
  # what do we get to plot?
  # rowPos = rep(seq(1,3),2)
  # colPos = c(rep(1,3),rep(2,3))
  # if (length(plotList) == 2) # if I have both plot list it means it is a 2 col plot
  # {
  #   png(filename=path_to_png, width = 24, height = 30, units = "cm",res = 600)
  #   
  #   grid.newpage()
  #   pushViewport(viewport(layout = grid.layout(nrow=3, ncol=2,
  #                                              widths = unit(12, "cm"),
  #                                              heights = unit(10, "cm"))))
  #   for (i in rowPos)
  #     for (j in colPos)
  #       print(plotList[[j]][[i]], vp = viewport(layout.pos.row = i, layout.pos.col = j))
  #   
  # } else {
  #   
  #   png(filename=path_to_png, width = 12, height = 30, units = "cm",res = 600)
  #   
  #   grid.newpage()
  #   pushViewport(viewport(layout = grid.layout(nrow=3, ncol=1,
  #                                              widths = unit(12, "cm"),
  #                                              heights = unit(10, "cm"))))
  #   
  #   for (i in rowPos)
  #     print(plotList[[1]][[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  # }
  # 
  dev.off()
  toc()
}

 
plotTraitRelative <- function(folder,Mat = T, PPMR = F,Sig = F, SpIdx = NULL, comments = T, window = NULL, what = c("3sp","9sp"), returnData = F)
{
  tic()
  dt = 0.1
  path_to_png = paste(folder,"/Biom&TraitsCompareFishedRelative.png",sep="")
  
  plotList <- list() # to store the plots
  buildList <- list()
  
  for (column in what)
  {
    if (column == "3sp") 
    {
      listPosition = 1
      folder2 <- paste(folder,"/3sp",sep="")
      title = "3 species"
    } else if (column == "9sp") 
    {
      listPosition = 2
      folder2 <- paste(folder,"/9sp",sep="")
      title = "9 species"
    } else stop("I don't know what to plot")
    
  dataList <- list()
  init <- get(load(paste(folder2,"/init/run1/run.Rdata",sep="")))
  
  for (z in 1:2)
  {
if (z == 1) multiSim <- bunchLoad(folder = paste(folder2,"/normal",sep=""))
if (z == 2) multiSim <- bunchLoad(folder = paste(folder2,"/fisheries",sep=""))
      

   
    if (comments) cat("Data ready\n")  
    # get the initial stuff
    SpIdx = sort(unique(init@params@species_params$species)) # determine spidx if not already given by user
    no_sp = length(SpIdx) # get the species number from this
    Wmat = init@params@species_params$w_mat[1:SpIdx[length(SpIdx)]] # get the mat size values
    TimeMax = init@params@species_params$timeMax[1] # and the length of the initialisation
    if (comments) cat("Initialisation done\n")
    
    
    # color gradient
    colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
    colGrad <- colfunc(no_sp)
    names(colGrad) <- seq(1,no_sp)
    if (no_sp == 9) colLine <- c("solid","dashed","solid","dashed","solid","longdash","solid","longdash","solid") else colLine = rep("solid",3) # 3,6 and 8 are dashed
    
    names(colLine) <- seq(1,no_sp)
    
    #TRAITS
    # get the plots
    TraitsList <- list()
    for (i in 1:length(multiSim))
    {
      if (comments) cat(sprintf("Using run %g\n",i))
      TraitsList[[i]] <- plotTraitsMulti(object = multiSim[[i]],Mat = Mat, PPMR = PPMR,Sig = Sig, returnData = T, SpIdx = SpIdx, Wmat = Wmat, TimeMax = TimeMax, window = window, Normalisation = F)
    }
    if (comments) cat("Plots loaded")
    speciesList <- list()
    for (j in SpIdx) # for every species
    {
      speciesData <- NULL
      for (i in 1:length(TraitsList)) # at each time
      {
        if (!is.null(TraitsList[[i]][[1]][[j]]))
        {
          a <- ggplot_build(TraitsList[[i]][[1]][[j]]) # take the plot data
          if (!is.null(a$data[[1]]$group))
          {
            a$data[[1]]$species <- j # add species identity in the data
            a$data[[1]]$group <- 100*j + i # change the group to the run number (up to 99 different runs to keep species identity as hundreds)
            
            speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
          }
        }
      }
      speciesList[[j]] <- speciesData #this is a list of all the species at each time
    }
    
    if (is.null(window)) window = c(-1,1)
    if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
    
    plot_datTrait <- do.call("rbind",speciesList)
    
    dataList[[z]] <- plot_datTrait
    
  }

  dataList[[1]]$y <- dataList[[2]]$y / dataList[[1]]$y 


    p3 <- ggplot(dataList[[1]])+
      geom_line(aes(x=x,y=y, group = group, color = as.factor(species), linetype = as.factor(species))) +
      scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
      scale_y_continuous(name = "Fisheries eta value / normal eta value", limits = c(0.8,1.3))+
      scale_color_manual(name = "Species", values = colGrad)+ # color gradient
      scale_linetype_manual(name = "Species",values = colLine) +
      theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom", 
            #legend.justification=c(1,1),
            legend.key = element_rect(fill = "white"))+ 
      guides(color = guide_legend(nrow=1)) +
      ggtitle(title)

    plotList[[listPosition]] <- p3
    buildList[[listPosition]] <- list(ggplot_build(p3),dataList[[1]]) # the color gradient changes for the first loop if I don't do this. Why? I have no fucking idea
    
  }
  plotList <- plotList[lapply(plotList,length)>0] # if a slot is empty
  
  # MULTIPLOT
  
  png(filename=path_to_png, width = 24, height = 10, units = "cm",res = 600)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow=1, ncol=2,
                                             widths = unit(12, "cm"),
                                             heights = unit(10, "cm"))))
  
  print(plotList[[1]], vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plotList[[2]], vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  dev.off()
  toc()
  
  if (returnData)
  {
    myData <- do.call(list ,lapply(buildList, function(x) x[[2]])) # take the dataframe if want to take a look at the data or play with it
    return(myData)
  }
}
# bar plot of alive and extinct to compare
plotBarSp <- function (object, print_it = F, title = T)
{
  if (title) title = "Proportion of alive/extinct phenotypes" else title = NULL
  
  if (is.list(object)) # if I want to do average stuff
  {
    print("this is a list of sim")
    avgList = list()
    for (x in 1:length(object))
    {
      sim = object[[x]]
      ID = sim@params@species_params # handy shortcut
      Phen_dat = data.frame(ID$species, ID$extinct) # data to wrok on
      plot_dat = matrix(data = NA, nrow = length(unique(ID$species)), ncol = 3, dimnames = list(unique(ID$species), c("w_inf", "alive", "ext")))
      for (i in unique(ID$species)) # do species by species
      {
        Phen_ext = Phen_dat[which(Phen_dat$ID.species == i), 2] 
        no_alive = length(Phen_ext[Phen_ext == 0])
        no_ext = length(Phen_ext) - no_alive
        plot_dat[i, ] = c(ID$w_inf[i], no_alive, no_ext) # get the data, number of alive and extinct pehnotype at last time step by w_inf
      }
      avgList[[x]] = plot_dat # store everything
    }
    
    # time to do stats
    
    avgMatrix <- abind(avgList, along=3) # convert the list in an array
    avgMean = apply(avgMatrix, c(1,2), mean) # do the stat manip
    avgSd = apply(avgMatrix, c(1,2), sd)
    plot_dat = data.frame(avgMean,avgSd[,-1])# prepare my data
    colnames(plot_dat) = c("w_inf","Malive","Mext","Salive","Sext")
    
    #prepare the data for ggplot
    ggplot_dat = data.frame(c(plot_dat$w_inf,plot_dat$w_inf),c(rep("alive",length(unique(ID$species))),rep("ext",length(unique(ID$species)))),c(plot_dat$Malive,plot_dat$Mext),c(plot_dat$Salive,plot_dat$Sext))
    colnames(ggplot_dat) = list("w_inf","state","mean","sd")
    
    # trying to do something for standars deviation but weird
    ggplot_dat2 = data.frame(rep(plot_dat$w_inf,3),c(rep("sdalive",length(unique(ID$species))),rep("sdext",length(unique(ID$species))),rep("ext",length(unique(ID$species)))),
                             c(plot_dat$Salive,plot_dat$Sext,plot_dat$Mext-plot_dat$Sext))
    colnames(ggplot_dat2) = list("w_inf","state","value")
    ggplot_dat2$value[ggplot_dat2$value <0] = 0
    # standard deviation: first I have sd for extinct AND alive, so cannot just do a sd at the junction between the 2
    # I'm trying to plot mean ext - sd ext + sd ext + sd alive so it should do a sd not symetric at the junction
    # the result is colorfully horrible
    
    
    ggplot_dat2 = ggplot_dat
    ggplot_dat2 = ggplot_dat2[,-3]
    ggplot_dat2$mean = 2*ggplot_dat$sd 
    
    p <- ggplot(ggplot_dat) +
      geom_bar(aes(x = w_inf, y = mean, fill = state), stat = "identity") +
      #geom_point(aes(x=w_inf, y = mean))+
      #geom_bar(data = ggplot_dat2, aes(x = w_inf, y = value, fill = state), stat = "identity", alpha = 0.5)+
      scale_x_log10(breaks = c(round(plot_dat$w_inf)), name = "Size in g" ) +
      scale_y_continuous(name = "Number of phenotypes") +
      #scale_colour_discrete(labels = c("alive","extinct"))+
      theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
            legend.justification=c(1,1),legend.key = element_rect(fill = "white")) +
      ggtitle(title)
    
    if (print_it) print(p)
    
  } else {
    
    ID = object@params@species_params # handy shortcut
    Phen_dat = data.frame(ID$species, ID$extinct) # data to wrok on
    plot_dat = matrix(data = NA, nrow = length(unique(ID$species)), ncol = 3, dimnames = list(unique(ID$species), c("w_inf", "alive", "ext")))
    for (i in unique(ID$species)) # do species by species
    {
      Phen_ext = Phen_dat[which(Phen_dat$ID.species == i), 2] 
      no_alive = length(Phen_ext[Phen_ext == 0])
      no_ext = length(Phen_ext) - no_alive
      plot_dat[i, ] = c(ID$w_inf[i], no_alive, no_ext) # get the data, number of alive and extinct pehnotype at last time step by w_inf
    }
    
    ggplot_dat = data.frame(c(plot_dat[,1],plot_dat[,1]),c(rep("alive",length(unique(ID$species))),rep("ext",length(unique(ID$species)))),c(plot_dat[,2],plot_dat[,3]))
    colnames(ggplot_dat) = list("w_inf","state","value")
    p <- ggplot(ggplot_dat) +
      geom_bar(aes(x = w_inf, y = value, fill = state), stat = "identity") +
      scale_x_log10(breaks = c(round(ggplot_dat$w_inf)), name = "Size in g" ) +
      scale_y_continuous(name = "Number of phenotypes") +
      theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
            legend.justification=c(1,1),legend.key = element_rect(fill = "white")) +
      ggtitle(title)
    
    if (print_it) print(p)
  }
  return(p)
}

#plot of trait variation by biomass variation through time
plotPerformance <- function(object)
{
  # get the biomass variation
  biomass <- getBiomass(object) # n * w * dw and sum by species
  
  SpIdx = NULL # getting rid of the species that went extinct at the beginning 
  for (i in unique(object@params@species_params$species))
    if (sum(biomass[,i]) != 0 & dim(object@params@species_params[object@params@species_params$species == i,])[1] != 1) 
      SpIdx = c(SpIdx,i)
  
  biomassTot = NULL
  biomassTemp = biomass
  colnames(biomassTemp) = object@params@species_params$species
  for (i in SpIdx)
  {
    biomassSp = biomassTemp[,which(colnames(biomassTemp) == i)]
    biomassSp = apply(biomassSp,1,sum)
    biomassTot = cbind(biomassTot,biomassSp)
  }
  colnames(biomassTot) = SpIdx
  # biomassTot is the biomass of species
  
  # Biom <- melt(biomassTot) # melt for ggplot
  # names(Biom) = list("time","sp","value")
  # # Due to log10, need to set a minimum value, seems like a feature in ggplot
  # min_value <- 1e-300
  # Biom <- Biom[Biom$value >= min_value,]
  # # take the first digit of the species column and put it in a new column
  # Biom$bloodline = sapply(Biom[,2], function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  
  #get the trait variation
  # little initialisation 
  SumPar = object@params@species_params #shortcut
  TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),0.1*SumPar$pop,0.1*SumPar$extinct,SumPar$w_mat,SumPar$beta,SumPar$sigma) # weird things happen without the as.numeric / *0.1 because dt
  colnames(TT) = c("Lineage","Ecotype","Apparition","Extinction","Maturation_size","PPMR","Diet_breadth")
  rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  TT = as.data.frame(TT) # I use TT later in the graphs (its like an artefact)
  for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = SumPar$timeMax[i]*0.1
  
  # Weighted mean
  # 1) matrix of summed abundance of mature ind at each time step
  truc = object@n
  # put 0 in object@n when w < w_mat
  for (i in 1:dim(truc)[1]) # for each time step
  {
    for (j in 1:dim(truc)[2]) # for each ecotypes
    {
      w_lim = SumPar$w_mat[j] # get the maturation size of the ecotype
      S <- numeric(length(object@params@w))
      S[sapply(w_lim, function(i) which.min(abs(i - object@params@w)))] <- 1 # find what w bin is the closest of the maturation size
      NoW_mat = which(S == 1) # what is the col number of that size bin
      truc[i,j,1:NoW_mat-1] <-0 # everything under this value become 0
    }
  }
  abundanceM = apply(truc, c(1,2),sum) # sum the abundance left 
  
  # 2) normalisation per species 
  colnames(abundanceM) = SumPar$species # phenotypes from same species have the same name
  abundanceNormal = matrix(0,nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])
  
  # I am getting rid of the species which went instinct at the begining and that went extinct without having mutants (no trait variation)
  SpIdx = NULL
  for (i in unique(SumPar$species))
    if (sum(abundanceM[,i]) != 0 & dim(SumPar[SumPar$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
      SpIdx = c(SpIdx,i)
  
  for (i in SpIdx)
  {
    abundanceSp = abundanceM # save to manipulate
    abundanceSp[,which(colnames(abundanceM) != i)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
    abundanceSp = sweep(abundanceSp,1,apply(abundanceSp,1,sum),"/") # normalise
    abundanceSp[is.nan(abundanceSp)] <-0 # when I divide by 0 I get nan
    abundanceNormal = abundanceNormal + abundanceSp # I just need to add them up to get the final matrix
  }
  
  # Now I have normalised abundance, I need to apply the trait I want to plot on them
  
  # Maturation size
  abundanceT = sweep(abundanceNormal,2,SumPar$w_mat,"*") # I use the normalised abundance and multiply by the trait value
  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
  names(dimnames(TotMean)) = list("time","species")
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
    TotMean[,i] = AMean
  }
  # Calculate variance and standard deviation
  statMS = list() # list with the stats of all species as I need to do this for each species separatly
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = SumPar$w_mat[SumPar$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","delta"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
      stat[j,3] = sqrt(variance)# calculate the standard deviation
      stat[j,4] = abs(object@params@species_params$w_mat[i] - meanSp[j])/object@params@species_params$w_mat[i] # variation between starting value and actual one, normalised
    }
    statMS[[i]] = as.data.frame(stat) # put in the list
  }
  # now I need to plot delta trait against biomass
  plotMatStore <- list()
  for (i in SpIdx)
  {
    plot_dat = data.frame(statMS[[i]][,c(1, 4)], biomassTot[,i])
    colnames(plot_dat) = list("time", "delta", "biomass")
    
    if (i == SpIdx[length(SpIdx)])
    {
      p <- ggplot(plot_dat) +
        geom_point(aes(x = delta, y = biomass, colour = time )) +
        scale_y_log10(name = NULL, limits = c(1e-2,2e-1), breaks = c(1 %o% 10^(-4:1)))+
        scale_x_continuous(name = "Normalised maturation size (percent of change)", limits = c(0,0.12))+
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white")) #none remove legend
    } else {
      p <- ggplot(plot_dat) +
        geom_point(aes(x = delta, y = biomass, colour = time )) +
        scale_y_log10(name = NULL, limits = c(1e-2,2e-1), breaks = c(1 %o% 10^(-4:1)))+
        scale_x_continuous(name = NULL, limits = c(0,0.12))+
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white")) #none remove legend
      
      
    }
    
    plotMatStore[[i]] <- p
  }
  
  
  # PPMR
  abundanceT = sweep(abundanceNormal,2,SumPar$beta,"*") # I use the normalised abundance and multiply by the trait value
  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
  names(dimnames(TotMean)) = list("time","species")
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
    TotMean[,i] = AMean
  }
  # Calculate variance and standard deviation
  statPPMR = list() # list with the stats of all species
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = SumPar$beta[SumPar$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","delta"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
      stat[j,3] = sqrt(variance)# calculate the standard deviation
      stat[j,4] = abs(object@params@species_params$beta[i] - meanSp[j]) # variation between starting value and actual one
    }
    statPPMR[[i]] = as.data.frame(stat) # put in the list
  }
  
  # now I need to plot delta trait against biomass
  plotPPMRStore <- list()
  for (i in SpIdx)
  {
    plot_dat = data.frame(statPPMR[[i]][,c(1, 4)], biomassTot[,i])
    colnames(plot_dat) = list("time", "delta", "biomass")
    
    if (i == SpIdx[length(SpIdx)])
    {
      p <- ggplot(plot_dat) +
        geom_point(aes(x = delta, y = biomass, colour = time )) +
        scale_y_log10(name = NULL, limits = c(1e-2,2e-1), labels = NULL)+
        scale_x_continuous(name = "PPMR variation",limits = c(0,15))+
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))
    } else {
      
      p <- ggplot(plot_dat) +
        geom_point(aes(x = delta, y = biomass, colour = time )) +
        scale_y_log10(name = NULL, limits = c(1e-2,2e-1), labels = NULL)+
        scale_x_continuous(name = NULL,limits = c(0,15))+
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none",
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))
    }
    
    plotPPMRStore[[i]] <- p
    
  }
  
  # Diet breath
  abundanceT = sweep(abundanceNormal,2,SumPar$sigma,"*") # I use the normalised abundance and multiply by the trait value
  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(SumPar$species)), dimnames = list(rownames(abundanceT),unique(SumPar$species)))
  names(dimnames(TotMean)) = list("time","species")
  
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
    TotMean[,i] = AMean
  }
  
  # Calculate variance and standard deviation
  # it is the sum of the difference between value and mean squared and multiplied by the weight
  
  statDB = list() # list with the stats of all species
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = SumPar$sigma[SumPar$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","delta"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
      stat[j,3] = sqrt(variance)# calculate the standard deviation
      stat[j,4] = abs(object@params@species_params$sigma[i] - meanSp[j]) # variation between starting value and actual one
    }
    statDB[[i]] = as.data.frame(stat) # put in the list
    #statDB[[i]] = stat# put in the list
  }
  
  # now I need to plot delta trait against biomass
  plotDBStore <- list()
  for (i in SpIdx)
  {
    plot_dat = data.frame(statDB[[i]][,c(1, 4)], biomassTot[,i])
    colnames(plot_dat) = list("time", "delta", "biomass")
    
    if (i == SpIdx[length(SpIdx)])
    {
      p <- ggplot(plot_dat) +
        geom_point(aes(x = delta, y = biomass, colour = time )) +
        scale_y_log10(name = NULL, limits = c(1e-2,2e-1), labels = NULL)+
        scale_x_continuous(name = "Feeding kernel variation", limits = c(0,0.2))+
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white")) 
    } else {
      
      p <- ggplot(plot_dat) +
        geom_point(aes(x = delta, y = biomass, colour = time )) +
        scale_y_log10(name = NULL, limits = c(1e-2,2e-1), labels = NULL)+
        scale_x_continuous(name = NULL, limits = c(0,0.2))+
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white")) 
    }
    plotDBStore[[i]] <- p
  }
  
  
  p <- multiplot(plotMatStore[[1]],plotMatStore[[2]], plotMatStore[[3]],
                 plotMatStore[[4]],plotMatStore[[5]],
                 plotMatStore[[6]],plotMatStore[[7]],plotMatStore[[8]],plotMatStore[[9]],
                 plotPPMRStore[[1]],plotPPMRStore[[2]], plotPPMRStore[[3]],
                 plotPPMRStore[[4]], plotPPMRStore[[5]],
                 plotPPMRStore[[6]],plotPPMRStore[[7]],plotPPMRStore[[8]],plotPPMRStore[[9]],
                 plotDBStore[[1]],plotDBStore[[2]],plotDBStore[[3]],
                 plotDBStore[[4]], plotDBStore[[5]],
                 plotDBStore[[6]],plotDBStore[[7]],plotDBStore[[8]],plotDBStore[[9]],
                 cols=3)
  
} 


# plot survival en reproduction success and show the lifetime reprod sucess depending on w_mat and sigma
plotFitness <- function(listObject,whatTime = NULL , where = getwd(), returnData = F, save_it = F, SpIdx = NULL,Mat = T, PPMR = T, Sigma = T, comments = T, window = c(-1,1), Wmat = NULL)
{
  if (comments) cat("plotFitness begins\n")
  if (comments) cat("starting fitness plots\n")
  
  if(is.list(listObject) == F) # if I get a sim alone, I put it in a list of length one to standardise everything
    listObject = list(listObject)
  
  if(is.null(whatTime)) # get the time max from the first sim (as they are the same in any case)
    whatTime = max(as.numeric(dimnames(listObject[[1]]@n)$time))
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  # determine the initial maturation size for normalisation purpose
  if (Mat && is.null(Wmat))
  {
    #Wmat = c(25, 79.056942, 250, 790.569415, 2500, 7905.694150, 25000)
    #Wmat = object@params@species_params$w_mat[1:length(object@params@species_params$species)]
  # if (sum(SpIdx == seq(1,9)) == length(SpIdx)) Wmat = c(2.5, 7.905694, 25, 79.056942, 250, 790.569415, 2500, 7905.694150, 25000) #eta = 0.25
  if (sum(SpIdx == seq(1,9)) == length(SpIdx)) Wmat = c(5.00000, 15.81139, 50.00000, 158.11388, 500.00000, 1581.13883, 5000.00000, 15811.38830, 50000.00000) #eta = 0.5
  if (sum(SpIdx == seq(1,3)) == length(SpIdx)) Wmat = c(5,50,500)
  
  if (comments) {
    cat("Wmat is")
    print(Wmat)
  }
  }
  

  
  listOutput = list()
  for (x in 1:length(listObject))
  {
    if (comments) cat(sprintf("using simulation %g\n",x))
    object = listObject[[x]]
    SumPar = object@params@species_params #shortcut
    abundance = object@n #shortcut  
    dimnames(abundance)$species = SumPar$species #phenotypes get their species name
    time_range = seq((whatTime-50),whatTime,1)
    

    window = window
    
    if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
    
    abundance=abundance[which(dimnames(abundance)$time == whatTime-50):which(dimnames(abundance)$time == whatTime),,]  #shorten the abundance matrix already
    
    # get the growth at each time step
    growth_mat = array(data = NA,dim = c(length(time_range),length(SumPar$ecotype),length(object@params@w)),dimnames = list(time_range,SumPar$ecotype,object@params@w))
    names(dimnames(growth_mat)) = list("time","species","size")
    for(j in time_range) growth_mat[which(dimnames(growth_mat)$time == j),,] = getEGrowth(object@params,object@n[j,,],object@n_pp[j,])
    meanGrowth = apply(growth_mat,c(2,3),mean) # do the mean already
    if (comments) cat(sprintf("growth calculated\n"))
    
    #get the spawn at each time step
    spawn_mat = array(data = NA,dim = c(length(time_range),length(SumPar$ecotype),length(object@params@w)),dimnames = list(time_range,SumPar$ecotype,object@params@w))
    names(dimnames(spawn_mat)) = list("time","species","size")
    for(j in time_range) spawn_mat[which(dimnames(spawn_mat)$time == j),,] = getESpawning(object@params,object@n[j,,],object@n_pp[j,])
    meanSpawn = apply(spawn_mat,c(2,3),mean) # do the mean already
    if (comments) cat(sprintf("spawn calculated\n"))
    
    # which species are not extinct?
    # print("Spidx before")
    # print(SpIdx)
    SpIdx = NULL
    for (i in unique(SumPar$species))
      if (sum(abundance[,which(dimnames(abundance)$species == i),]) != 0) # if species is extinct, do not plot anything
      { a = abundance[,which(dimnames(abundance)$species == i),]   # if only one phen do not plot as apply will bug
      if (length(dim(a[,apply(a,2,sum) != 0,])) > 2)
        SpIdx = c(SpIdx,i)
      }
    SpIdx <- sort(SpIdx)
    
    if (comments) cat(sprintf("SpIdx is\n"))
    if (comments) print(SpIdx)
    # plots
    plotList = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) #need 9 slots in there, even if they're not filled
    dataList = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
    for (i in SpIdx)
    {
      abundanceSp = abundance[,which(dimnames(abundance)$species == i),] # get the abundance matrix only of the species concerned
      dimnames(abundanceSp)$species = SumPar[SumPar$species == i,]$ecotype # get their phen names back
      abundanceSp = abundanceSp[,apply(abundanceSp,2,sum) != 0,] # get rid of extinct species
      meanPhen = apply(abundanceSp,c(2,3),mean) #get mean over time
      sdPhen = apply(abundanceSp,c(2,3),sd) # get standard deviation to exclude species if vary too much (not implemented yet)
      growthSp = meanGrowth[c(dimnames(abundanceSp)$species),] # get the mean growth of the selected phenotypes
      
      # survival probability at size per phenotypes
      survp <- (meanPhen*growthSp)/(meanPhen[,1]*growthSp[,1])
      
      plot_dat1 = melt(survp) # prep data for ggplot
      
      p1 <- ggplot(plot_dat1) +
        geom_line(aes(x = size, y = value, colour = as.factor(species)))+
        scale_x_log10(breaks = c(1 %o% 10^(-3:5)))+
        scale_y_log10(name = "Survival probability")+        
        #scale_colour_manual(values=cbPalette, name = "Phenotypes")+ # colorblind
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle(NULL)
      
      # reproductive success
      spawnSp = meanSpawn[c(dimnames(abundanceSp)$species),] # get the mean spawn of the selected phenotypes
      
      repros <- survp*spawnSp/growthSp
      
      plot_dat2 = melt(repros)
      
      p2 <- ggplot(plot_dat2) +
        geom_line(aes(x = size, y = value, colour = as.factor(species)))+
        scale_x_log10(breaks = c(1 %o% 10^(-3:5)))+
        scale_y_log10(name = "Reproductive success")+        
        #scale_colour_manual(values=cbPalette, name = "Phenotypes")+ # colorblind
        theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle(NULL)
      
      # lifetime repro success
      
      r0<-rowSums(sweep(repros,2,object@params@dw,"*"),na.rm=T)
      
      plot_dat3 <- data.frame(r0, SumPar$w_mat[SumPar$ecotype %in% as.numeric(c(dimnames(abundanceSp)$species))],SumPar$sigma[SumPar$ecotype %in% as.numeric(c(dimnames(abundanceSp)$species))],SumPar$beta[SumPar$ecotype %in% as.numeric(c(dimnames(abundanceSp)$species))]) # %in% is awesome
      colnames(plot_dat3) = list("r0","w_mat","sigma","beta")
      
      #normalise the data here
      #plot_dat3$w_mat = (plot_dat3$w_mat- SumPar$w_mat[i])/plot_dat3$w_mat
      plot_dat3$w_mat = (plot_dat3$w_mat- Wmat[i])/plot_dat3$w_mat
      
      
      p3 <-ggplot(plot_dat3)+
        geom_point(aes(x=sigma,y=w_mat,colour=r0, size = r0))+
        scale_x_continuous(name = "Width of feeding kernel")+
        scale_y_continuous(name = "Maturation size")+
        scale_colour_gradientn(colours = terrain.colors(10), name = "r0")+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(NULL)
      
      
      pW_mat <-ggplot(plot_dat3)+
        geom_point(aes(x=w_mat,y=r0))+
        scale_x_continuous(name = "Maturation size", limits = window)+
        scale_y_continuous(name = "r0")+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(NULL)
      
      pSigma <-ggplot(plot_dat3)+
        geom_point(aes(x=((sigma-SumPar$sigma[i])/sigma),y=r0))+
        geom_smooth(aes(x=((sigma-SumPar$sigma[i])/sigma),y=r0), se = F) +
        scale_x_continuous(name = "Width of feeding kernel", limits = window)+
        scale_y_continuous(name = "r0")+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(NULL)
      
      
      pBeta <-ggplot(plot_dat3)+
        geom_point(aes(x=((beta-SumPar$beta[i])/beta),y=r0))+
        scale_x_continuous(name = "Preferred PPMR", limits = window)+
        scale_y_continuous(name = "r0")+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(NULL)
      
      #create plot_data4 with the trend line of feeding kernel
      
      plot_dat4 = NULL
      if(Sigma){
        if (sum(abundanceSp) != 0) # in the case the species is extinct and nothing has been plotted
        {
          plot_TL <- ggplot_build(pSigma)
          plot_dat4 <- plot_TL$data[[2]][,1:2] #is the trendline
          plot_dat4$run = rep(x,dim(plot_dat4)[1]) # add an identification col
        }
      }
      plotList[[i]] = list(p1,p2,p3,pW_mat,pSigma,pBeta)
      dataList[[i]] = list(plot_dat1,plot_dat2,plot_dat3,plot_dat4)
      # print("plotdat")
      # print(plot_dat3)
    }
    
    listOutput[[x]] = list(plotList,dataList)
  }
  if (returnData) {
    if (comments) cat("plotFitness ends\n")
    return(listOutput)
  }
  
  #save the plots
  if(save_it)
  {
    for (i in SpIdx) #does not really work at the moment, unless I want to save the figures on the disk, you can choose which 3 plots to use
    {
      path_to_png = paste(where,"/fitnessSp",i,".png",sep="")
      png(filename=path_to_png,width = 20, height = 30, units = "cm",res = 600)
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow=3, ncol=1,
                                                 widths = unit(20, "cm"),
                                                 heights = unit(c(10,10,10), "cm"))))
      print(plotList[[i]][[3]], vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(plotList[[i]][[4]], vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
      print(plotList[[i]][[5]], vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      
      dev.off()
    }
  }
  # plots w_mat vs r0 for all the phenotypes of 1 run
  # plot_dat <- do.call(rbind,lapply(dataList, function (x) x[[3]])) # did it first try, Im da best
  # plot_dat$species = sapply(rownames(plot_dat), function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1]) # create new column with species name (first digit of phen)
  # 
  # p<- ggplot(plot_dat)+
  #   geom_point(aes(x = w_mat, y = r0,colour = beta,size = sigma))+
  #   scale_x_log10() +
  #   scale_y_log10()
  
  
  if (comments) cat("plotFitness ends\n")
  return(p)
}

# do the multiplot of plotFitness across several run of the same type in the same folder
plotFitnessMulti <- function(folder, returnData = F, SpIdx = NULL, Mat = T, Sig = T, PPMR = T, comments = T, whatTime = NULL, window = c(-1,1), Wmat = NULL)
{
  if (comments) cat("plotFitnessMulti begins\n")
  
  if (comments) cat("starting fitness plot\n")
  listOfSim <- bunchLoad(folder) #load the sims from a folder
  if (comments) cat("sim loaded\n")
  
  window = window
  
  if(is.null(whatTime)) # get the time max from the first sim (as they are the same in any case)
    whatTime = max(as.numeric(dimnames(listOfSim[[1]]@n)$time))
  if (comments) cat(sprintf("time selected is %g\n",whatTime))
  
  if (is.null(SpIdx)) SpIdx = unique(listOfSim[[1]]@params@species_params$species)
  fitnessOutput = plotFitness(listOfSim, returnData = T, SpIdx = SpIdx, Sigma = Sig, PPMR = PPMR, whatTime = whatTime, window = window, Wmat = Wmat) #get the fitness data for plots
  if (comments) cat("fitness data obtained\n")
  SumPar = listOfSim[[1]]@params@species_params #shortcut
  
  #windows scale do not work for the moment
  # a = sort(SumPar$sigma) #get sigma values
  # b = c(head(a,1),tail(a,1)) # get both end of the spectrum
  # c = abs((b-SumPar$sigma[1])/b) # normalise
  # d = sort(c,decreasing = T)[1] # get the biggest one
  # window = c( - round(d, digits = 2), round(d, digits = 2)) # this is the maximum limit 
  window = window
  if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
  
  speciesList = list()
  trendList = list()
  for (i in SpIdx)
  {
    #print(i)
    speciesList[[i]] <- do.call(rbind,lapply(fitnessOutput, function (x) x[[2]][[i]][[3]])) #badaboum: create a list with data (traits + r0) for each species across runs
    trendList[[i]] <- do.call(rbind,lapply(fitnessOutput, function (x) x[[2]][[i]][[4]])) #get the trend line values
  }
  # if I fail somewere with SpIdx before, I correct it here
  temp = NULL
  for (i in SpIdx)
    if (!is.null(dim(speciesList[[i]])))
      if (dim(speciesList[[i]])[1] != 0)
        temp = c(temp, i)
  SpIdx = temp
  
  
  plotMatStore = list()
  plotSigStore = list()
  plotPPMStore = list()
  for (i in SpIdx) # do the plots for each species and traits
  {
    if (comments) cat(sprintf("plot %i\n",i))
    if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
    title = NULL
    
    if (Mat){
      plotMatStore[[i]] <-ggplot(speciesList[[i]])+
        geom_point(aes(x=w_mat,y=r0))+ 
        scale_x_continuous(name = NULL, limits = window)+
        scale_y_continuous(name = sprintf("Species %i",i))+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(title[1])
    }
    if (Sig)
    {
      plotSigStore[[i]] <-ggplot(speciesList[[i]])+
        geom_point(aes(x=((sigma-SumPar$sigma[i])/sigma),y=r0))+
        geom_line(data =trendList[[i]], aes(x=x,y=y,group = run))+
        scale_x_continuous(name = NULL, limits = window)+
        scale_y_continuous(name = NULL, labels = NULL)+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(title[2])
    }
    if(PPMR)
    {
      plotPPMStore[[i]] <-ggplot(speciesList[[i]])+
        geom_point(aes(x=((beta-SumPar$beta[i])/beta),y=r0))+
        scale_x_continuous(name = NULL, limits = window)+
        scale_y_continuous(name = NULL, labels = NULL)+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(title[3])
    }
  }
  
  if (returnData) return(list(plotMatStore,plotSigStore,plotPPMStore))
  
  if (comments) print("plot done, starting the multiplot")
  # do the multiplot
  path_to_png = paste(folder,"/fitnessMulti",whatTime,".png",sep="")
  
  if (PPMR == F & Sig == F)
  {
    png(filename=path_to_png,width = 30, height = 30, units = "cm",res = 1000)
    rowPos = rep(c(1,2,3),3)
    colPos = c(rep(1,3),rep(2,3),rep(3,3))
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=3, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(10,3), "cm"))))
    for (i in SpIdx)
    {print(i)
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = rowPos[i], layout.pos.col = colPos[i])) 
    }
    
  } else {
    
    png(filename=path_to_png,width = 30, height = 45, units = "cm",res = 1000)
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=9, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(5,9), "cm"))))
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)) 
      print(plotSigStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 2)) 
      print(plotPPMStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 3)) 
    }
  }
  dev.off()
}

plotFitnessMultiOverlap <- function(directory, SpIdx = NULL, Mat = T, PPMR = T, Sig = T, comments = T, window = c(-1,1), oneSpMode = F, whatTime = NULL, save_it = T)
{
  
  normalData = F # to determine what data I have
  fisheriesData = F
  
  if (Mat) # get the initial wmat values
  {
    if(dir.exists(file.path(paste(directory,"/init",sep=""))))
    {
      WmatSim <- get(load(paste(directory,"/init/run1/run.Rdata",sep="")))
    Wmat = WmatSim@params@species_params$w_mat[1:SpIdx[length(SpIdx)]]
    if (comments) {
      cat("Wmat is")
      print(Wmat)}
    } else {cat("Could not find Wmat values, taking the default ones\n")}
  }
  
  
  #NO FISHERIES PART
  if(dir.exists(file.path(paste(directory,"/normal",sep=""))))
  {
    normalData = T
    if (comments) cat("Simulations without fisheries\n")
    plotListNormal <- plotFitnessMulti(folder = paste(directory,"/normal",sep=""),returnData = T, SpIdx = SpIdx, Mat = Mat, PPMR = PPMR, Sig = Sig, window = window, Wmat = Wmat, whatTime = whatTime)
    # it should give me a list of 3: plotwmat,plotsigma,plotppmr
  }
  #FISHERIES PART
  if(dir.exists(file.path(paste(directory,"/fisheries",sep=""))))
  {
    fisheriesData = T
    if (comments) cat("Simulations with fisheries\n")
    plotListFish <- plotFitnessMulti(folder = paste(directory,"/fisheries",sep=""),returnData = T, SpIdx = SpIdx, Mat = Mat, PPMR = PPMR, Sig = Sig, window = window, Wmat = Wmat, whatTime = whatTime)
  }
  
  plotMatStore = list()
  plotSigStore = list()
  plotPPMStore = list()
  
  if (is.null(SpIdx)) SpIdx = seq(1,9,1) #no time now
  
  if (normalData && fisheriesData){ #if I have both
    for (i in SpIdx) # do the plots for each species and traits
    {
      # get data
      
      if (comments) cat(sprintf("Plot for species %g\n",i))
      
      #if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
      title = NULL
      
      if (Mat){
        if(!is.null(plotListNormal[[1]][[i]])) plot_datN_mat <- ggplot_build(plotListNormal[[1]][[i]]) else plot_datN_mat = NULL
        if(!is.null(plotListFish[[1]][[i]]))   plot_datF_mat <- ggplot_build(plotListFish[[1]][[i]]) else plot_datF_mat = NULL
        
        
        plotMatStore[[i]] <-ggplot(plot_datN_mat$data[[1]])+
          geom_point(aes(x=x,y=y))+
          geom_point(data = plot_datF_mat$data[[1]], aes(x=x,y=y), shape = 3, size = 3, color = "red") +
          scale_x_continuous(name = NULL, limits = window)+
          scale_y_continuous(name = sprintf("Species %i",i))+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[1])
      }
      if (Sig)
      {
        plot_datN_sig <- ggplot_build(plotListNormal[[2]][[i]])
        plot_datF_sig <- ggplot_build(plotListFish[[2]][[i]])
        
        plotSigStore[[i]] <-ggplot(plot_datN_sig$data[[1]])+
          geom_point(aes(x=x,y=y))+
          geom_point(data = plot_datF_sig$data[[1]], aes(x=x,y=y), shape = 3, size = 3, color = "red") +
          geom_line(data = plot_datN_sig$data[[2]],aes(x=x,y=y,group=group))+
          geom_line(data = plot_datF_sig$data[[2]], aes(x=x,y=y,group=group), linetype = "dashed", color = "red" ) +
          scale_x_continuous(name = NULL, limits = c(-0.25,0.25))+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[2])
      }
      if (PPMR)
      {
        plot_datN_ppm <- ggplot_build(plotListNormal[[3]][[i]])
        plot_datF_ppm <- ggplot_build(plotListFish[[3]][[i]])
        
        plotPPMStore[[i]] <-ggplot(plot_datN_ppm$data[[1]])+
          geom_point(aes(x=x,y=y))+
          geom_point(data = plot_datF_ppm$data[[1]], aes(x=x,y=y), shape = 3, size = 3, color = "red") +
          scale_x_continuous(name = NULL, limits = c(-0.25,0.25))+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[3])
      }
    }
  } else { # if I have only one
    
    if (fisheriesData) plotListNormal = plotListFish # if I have only fisheries, I convert everything to normal as it is the same
    
    for (i in SpIdx) # do the plots for each species and traits
    {
      # get data
      
      if (comments) cat(sprintf("Plot for species %g\n",i))
      
      if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
      
      # print("lastdat")
      # print(ggplot_build(plotListNormal[[1]][[i]]))
      
      
      if (Mat){
        plot_datN_mat <- ggplot_build(plotListNormal[[1]][[i]])
        if (oneSpMode) plot_datN_mat$data[[1]] =  plot_datN_mat$data[[1]][!duplicated(plot_datN_mat$data[[1]]$x),] # little cheat to get rid of duplicates that SHOULD have the same fucking fitness
        #plot_datF_mat <- ggplot_build(plotListFish[[1]][[i]])
        
        plotMatStore[[i]] <-ggplot(plot_datN_mat$data[[1]])+
          geom_point(aes(x=x,y=y))+
          #geom_point(data = plot_datF_mat$data[[1]], aes(x=x,y=y), shape = 3, size = 3, color = "red") +
          scale_x_continuous(name = NULL, limits = window)+
          scale_y_continuous(name = sprintf("Species %i",i))+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[1])
      }
      if (Sig)
      {
        plot_datN_sig <- ggplot_build(plotListNormal[[2]][[i]])
        #plot_datF_sig <- ggplot_build(plotListFish[[2]][[i]])
        
        plotSigStore[[i]] <-ggplot(plot_datN_sig$data[[1]])+
          geom_point(aes(x=x,y=y))+
          #geom_point(data = plot_datF_sig$data[[1]], aes(x=x,y=y), shape = 3, size = 3, color = "red") +
          geom_line(data = plot_datN_sig$data[[2]],aes(x=x,y=y,group=group))+
          #geom_line(data = plot_datF_sig$data[[2]], aes(x=x,y=y,group=group), linetype = "dashed", color = "red" ) +
          scale_x_continuous(name = NULL, limits = c(-0.25,0.25))+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[2])
      }
      if (PPMR)
      {
        plot_datN_ppm <- ggplot_build(plotListNormal[[3]][[i]])
        # plot_datF_ppm <- ggplot_build(plotListFish[[3]][[i]])
        
        plotPPMStore[[i]] <-ggplot(plot_datN_ppm$data[[1]])+
          geom_point(aes(x=x,y=y))+
          #geom_point(data = plot_datF_ppm$data[[1]], aes(x=x,y=y), shape = 3, size = 3, color = "red") +
          scale_x_continuous(name = NULL, limits = c(-0.25,0.25))+
          scale_y_continuous(name = NULL, labels = NULL)+
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
                legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
          ggtitle(title[3])
      }
    }
  }
  # do the multiplot
  
  
  temp = NULL
  for (i in SpIdx)
    if (!is.null(plotMatStore[[i]]$data[1][[1]])) # system D, if it's null it means that there is no data and the plot is going to give an error
      temp = c(temp, i)
  SpIdx = temp
  
  cat("final Spidx is set to:")
  print(SpIdx)
  
  path_to_png = paste(directory,"/fitnessMulti",whatTime,".png",sep="")
  
  if (PPMR == F & Sig == F)
  {
    gridID = matrix(c(1,rep(2,3),rep(3,5),rep(1,2),rep(2,4),rep(3,3)), nrow = 2, ncol = 9,byrow = T, dimnames = list(c("nRow","nCol"),seq(1,9))) # grid dimension depending on sp number
    rowPos = rep(1:gridID[1,length(SpIdx)],gridID[2,length(SpIdx)]) # row position depending on the number of plots
    colPos = NULL
    for (i in 1:gridID[2,length(SpIdx)]) colPos = c(colPos,rep(i,gridID[1,length(SpIdx)])) # col position depending on the number of plots
    
    noRow = rowPos[length(rowPos)] # last value to determine number of row/col
    noCol = colPos[length(colPos)]
    
    png(filename=path_to_png, width = 10*noCol, height = 8*noRow, units = "cm",res = 1000)
    
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=noRow, ncol=noCol, 
                                               widths = unit(rep(10,noCol), "cm"), 
                                               heights = unit(rep(8,noRow), "cm"))))
    
    for (i in 1:length(SpIdx))
    {print(i)
      print(plotMatStore[[SpIdx[i]]], vp = viewport(layout.pos.row = rowPos[i], layout.pos.col = colPos[i])) 
    }
    
  } else {
    png(filename=path_to_png,width = 30, height = 45, units = "cm",res = 1000)
    
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=9, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(5,9), "cm"))))
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)) 
      print(plotSigStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 2)) 
      print(plotPPMStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 3)) 
    }
  }
  dev.off()
  #if (save_it) saveRDS(object = plotMatStore, file = paste(directory,"/fitnessMulti",whatTime,".rds",sep=""))
}


# function that plot the fitness at different time steps per species
plotFitnessTime <- function(folder, returnData = F, SpIdx = NULL, Mat = T, Sig = T, PPMR = T, comments = T, whatTime = NULL, window = NULL)
{
  if (is.null(window)) window = c(-1,1)
  if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
  
  if (is.null(whatTime)) whatTime = c(1500.1,1999.1)
  
  listFitnessPlots <-list() #stack the plots at different times
  for (i in whatTime)
    listFitnessPlots[[i]] <- plotFitnessMulti(folder = folder,Mat = Mat, PPMR = PPMR, Sig = Sig, returnData = T, whatTime = i, SpIdx = SpIdx)
  listFitnessPlots <- listFitnessPlots[lapply(listFitnessPlots,length)>0]
  
  #SpIdx <- seq(1,9)
  speciesList <- list()
  for (j in SpIdx) # for every species
  {
    speciesData <- NULL
    for (i in 1:length(listFitnessPlots)) # at each time
    {
      if (!is.null(listFitnessPlots[[i]][[1]][[j]]))
      {
        a <- ggplot_build(listFitnessPlots[[i]][[1]][[j]]) # take the plot data
        a$data[[1]]$group <- whatTime[i] # change the group to the right time
        speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
      }
    }
    speciesList[[j]] <- speciesData #this is a list of all the species at each time
  }
  
  # plot time
  plotMatStore = list()
  plotSigStore = list()
  plotPPMStore = list()
  
  temp = NULL
  for (i in SpIdx) if (!is.null(speciesList[[i]])) temp = c(temp,i) #update SpIdx if species are extinct at all times
  SpIdx = temp
  for (i in SpIdx) # do the plots for each species and traits
  {
    if (comments) cat(sprintf("plot %i\n",i))
    #if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
    
    if (Mat){
      plotMatStore[[i]] <-ggplot(speciesList[[i]])+
        geom_point(aes(x=x,y=y,colour=as.factor(group))) +
        scale_x_continuous(name = NULL, limits = window)+
        scale_y_continuous(name = sprintf("Species %i",i), limits = c(0,NA))+
        scale_color_grey(name = "Species")+ # grey
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"), legend.position="none", text = element_text(size=20),
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(NULL)
    }
    if (Sig)
    {
      plotSigStore[[i]] <-ggplot(speciesList[[i]])+
        geom_point(aes(x=x,y=y,colour=as.factor(group))) +
        geom_line(data =trendList[[i]], aes(x=x,y=y,group = run))+
        scale_x_continuous(name = NULL, limits = window)+
        scale_y_continuous(name = NULL, labels = NULL)+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(title[2])
    }
    if(PPMR)
    {
      plotPPMStore[[i]] <-ggplot(speciesList[[i]])+
        geom_point(aes(x=x,y=y,colour=as.factor(group))) +
        scale_x_continuous(name = NULL, limits = window)+
        scale_y_continuous(name = NULL, labels = NULL)+
        theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),# legend.position="none", 
              legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        ggtitle(title[3])
    }
  }
  
  if (comments) print("plot done, starting the multiplot")
  # do the multiplot
  path_to_png = paste(folder,"/fitnessTime.png",sep="")
  
  if (PPMR == F & Sig == F)
  {
    png(filename=path_to_png,width = 30, height = 30, units = "cm",res = 500)
    rowPos = rep(c(1,2,3),3)
    colPos = c(rep(1,3),rep(2,3),rep(3,3))
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=3, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(10,3), "cm"))))
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = rowPos[i], layout.pos.col = colPos[i])) 
    }
    
  } else {
    
    png(filename=path_to_png,width = 30, height = 45, units = "cm",res = 500)
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=9, ncol=3, 
                                               widths = unit(rep(10,3), "cm"), 
                                               heights = unit(rep(5,9), "cm"))))
    for (i in SpIdx)
    {
      print(plotMatStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)) 
      print(plotSigStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 2)) 
      print(plotPPMStore[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 3)) 
    }
  }
  dev.off()
  
  if (returnData)
  {
    
    
    return(speciesList)
  }
  
}


# 4 plots that describe the dynamics of a specific cohort
plotCohort <- function(object, dt = 0.1, T = 20/dt, t_start = object@params@species_params$timeMax[1]*dt - T,save_it = T, comments=T, path_to_png = NULL)
{
  
  if (is.null(path_to_png)) path_to_png <- getwd()
  
  no_sp = length(unique(object@params@species_params$species))
  effort = sim@effort[1]
  sex_ratio = 0.5
  # if (is.null(path_to_save)) path_to_save = paste(getwd(),"/temporary",sep="")
  # ifelse(!dir.exists(file.path(path_to_save)), dir.create(file.path(path_to_save),recursive = T), FALSE)
  
  # TRACK THE DENSITY OF A COHORT THROUGH TIME 
  
  # LOOK AT GROWTH AND SURVIVAL
  cohortW = array(0, c(no_sp, T+1)); # row vector for following cohort weight
  cohortS = array(0, c(no_sp, T+1)); # vector for cohort survival
  cohortR = array(0, c(no_sp, T+1)); # vector for cohort spawning
  cohortR_sol = array(0, c(no_sp, T+1)); # vector for cohort spawn at size
  
  # NEWBORNS OVER LIFETIME
  cohortW[,1] = object@params@w[1]; # log weight initially (newborn)
  cohortS[,1] = object@n[t_start,,1]; # initial population in spectrum
  
  for (q in seq(1,no_sp)){ 
    for (t in seq(1,T)){ # within time period you're interested in
      cohortWprev = max(which(cohortW[q,t] - object@params@w >= 0)) # weight bin of cohort from last time step 
      growth = getEGrowth(object@params,n = object@n[t_start+t-1,,],n_pp = object@n_pp[t_start+t-1,])
      cohortW[q,t+1] = cohortW[q,t]+dt*growth[q,cohortWprev] # using growth rate in that bin to update to cohortW(t-t_start+1)
      z = getZ(object = object@params, n = object@n[t_start+t-1,,],n_pp = object@n_pp[t_start+t-1,], effort = effort)
      cohortS[q,t+1] = cohortS[q,t]*exp(-dt*z[q,cohortWprev]) # updating amount surviving using death rate
      e_spawning <- getESpawning(object = object@params, n = object@n[t_start+t-1,,]*cohortS[,t]/cohortS[,1]
                                 ,n_pp = object@n_pp[t_start+t-1,])
      e_spawning_pop <- apply((e_spawning*object@n[t_start+t-1,,]*cohortS[,t]/cohortS[,1]),1,"*",object@params@dw)
      rdi <- sex_ratio*(e_spawning_pop * object@params@species_params$erepro)/object@params@w[object@params@species_params$w_min_idx] 
      cohortR[q,t+1] = cohortR[q,t] + dt*rdi[cohortWprev,q]  #[q,cohortWprev]
      cohortR_sol[q,t+1] = dt*rdi[cohortWprev,q] 
    }
  }
  
  # plots
  min_value <- 1e-300
  title = paste("Cohort",t_start)
  
  dimnames(cohortW) <- list(seq(1,no_sp),seq(1,T+1))
  names(dimnames(cohortW)) <- list("species","time")
  plot_datG <- melt(cohortW)
  
  pG <- ggplot(plot_datG) +
    geom_line(aes(x = time, y = value, color = as.factor(species))) +
    scale_x_continuous(name = "Time (years)") +
    scale_y_continuous(name = "Size (g)", trans = "log10")+
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
    ggtitle(title)
  
  dimnames(cohortS) <- list(seq(1,no_sp),seq(1,T+1))
  names(dimnames(cohortS)) <- list("species","time")
  plot_datS <- melt(cohortS)
  plot_datS <- plot_datS[plot_datS$value >= min_value,]
  
  
  pS <- ggplot(plot_datS) +
    geom_line(aes(x = time, y = value, color = as.factor(species))) +
    scale_x_continuous(name = "Time (years)") +
    scale_y_continuous(name = "Biomass remaining", trans = "log10")+
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
    ggtitle(NULL)
  
  dimnames(cohortR) <- list(seq(1,no_sp),seq(1,T+1))
  names(dimnames(cohortR)) <- list("species","time")
  plot_datR <- melt(cohortR)
  plot_datR <- plot_datR[plot_datR$value >= min_value,]
  
  
  pR <- ggplot(plot_datR) +
    geom_line(aes(x = time, y = value, color = as.factor(species))) +
    scale_x_continuous(name = "Time (years)") +
    scale_y_continuous(name = "Total spawn", trans = "log10") +
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
    ggtitle(NULL)
  dimnames(cohortR_sol) <- list(seq(1,no_sp),seq(1,T+1))
  names(dimnames(cohortR_sol)) <- list("species","time")
  plot_datR_sol <- melt(cohortR_sol)
  plot_datR_sol <- plot_datR_sol[plot_datR_sol$value >= min_value,]
  
  pR_sol <- ggplot(plot_datR_sol) +
    geom_line(aes(x = time, y = value, color = as.factor(species))) +
    scale_x_continuous(name = "Time (years)") +
    scale_y_continuous(name = "Spawn output", trans = "log10")+
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
    ggtitle(NULL)
  
  if (save_it)
  {
    png(filename=paste(path_to_png,"/",title,".png",sep=""),width = 20, height = 20, units = "cm",res = 500)
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(nrow=2, ncol=2, 
                                               widths = unit(rep(10,2), "cm"), 
                                               heights = unit(rep(10,2), "cm"))))
    
    print(pG, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) 
    print(pS, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) 
    print(pR, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(pR_sol, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) 
    
    dev.off()
  }
  
  #return(list(pG,pS,pR,pR_sol))
  
}

# plot the number of alive phenotype per species and make distincttion between scenarios
plotNoPhen <- function(folder, comments = T, print_it = T, returnData = F, SpIdx = NULL, dt = 0.1)
{
  
  dataList <- list()
  for (dirContentIdx in dir(paste(folder,"/init",sep=""))) # for each sim in the init folder
  {
    if (file.exists(paste(folder,"/init/",dirContentIdx,"/run.Rdata",sep = "")) && file.exists(paste(folder,"/normal/",dirContentIdx,"/run.Rdata",sep = "")) && file.exists(paste(folder,"/fisheries/",dirContentIdx,"/run.Rdata",sep = ""))) # check that the sim has following scenarios as well
    {
      simList <- list(get(load(paste(folder,"/init/",dirContentIdx,"/run.Rdata",sep=""))),get(load(paste(folder,"/normal/",dirContentIdx,"/run.Rdata",sep=""))),get(load(paste(folder,"/fisheries/",dirContentIdx,"/run.Rdata",sep="")))) # load the 3 scenarios (init/normal/fish) if it's the case
      
      # create a dataframe counting the number of phenotypes alive at each time steps per species. 
      #We will start by a dataframe per species and then melt
      if (is.null(SpIdx)) SpIdx = sort(unique(simList[[1]]@params@species_params$species))
      
      # stitch together the params of the 3 scenarios
      SummaryParamsN = simList[[1]]@params # structure and basics params
      SummaryParamsF = SummaryParamsN
      # fishing differenciation
      simList[[1]]@params@species_params$Fished <- FALSE
      simList[[2]]@params@species_params$Fished <- FALSE
      simList[[3]]@params@species_params$Fished <- TRUE
      # timeMax?
      timeMax = simList[[1]]@params@species_params$timeMax[1] + simList[[2]]@params@species_params$timeMax[1]
      # create 2 dataframes
      SummaryParamsN@species_params = rbind(simList[[1]]@params@species_params,simList[[2]]@params@species_params) # get the sp ID from everything
      SummaryParamsN@species_params$timeMax = timeMax # update the timemax 
      a <- SummaryParamsN@species_params
      a <- a[order(a$ecotype, a$extinct, decreasing=TRUE),] # weird 3 lines to get rid of duplicates and keep the ones with the extinction value
      a <- a[!duplicated(a$ecotype),]
      SummaryParamsN@species_params = a[order(a$pop,a$ecotype),]
      
      SummaryParamsF@species_params = rbind(simList[[1]]@params@species_params,simList[[3]]@params@species_params) # get the sp ID from everything
      SummaryParamsF@species_params$timeMax = timeMax # update the timemax 
      a <- SummaryParamsF@species_params
      a <- a[order(a$ecotype, a$extinct, decreasing=TRUE),] # weird 3 lines to get rid of duplicates and keep the ones with the extinction value
      a <- a[!duplicated(a$ecotype),]
      SummaryParamsF@species_params = a[order(a$pop,a$ecotype),]
      
      # so here I should have 2 dataframes that sums up the phenotypes ID and distinghuished if its form fished or unfished scenario
      scenarioList <- list()
      for (scenario in c("normal","fished")){
        #short cut the df
        if (scenario == "normal") 
        {
          SumPar = data.frame(SummaryParamsN@species_params$species,SummaryParamsN@species_params$pop,SummaryParamsN@species_params$extinct)
          listPosition = 1 } else {
            SumPar = data.frame(SummaryParamsF@species_params$species,SummaryParamsF@species_params$pop,SummaryParamsF@species_params$extinct)
            listPosition = 2
          }
        colnames(SumPar) = c("species","pop","exit")
        SumPar$pop[which(SumPar$pop == 0)] = 1 # initial species pop at 1, not 0
        for (i in 1:dim(SumPar)[1]) if (SumPar$exit[i] == 0) SumPar$exit[i] = timeMax # the not extinct ones get the end sim as extinction value
        
        DemList <- list()
        for (x in SpIdx) # for every species
        {
          DemCount = matrix(data = 0, nrow = (timeMax*dt), ncol =2) # create a matrix that count pop and extinction for each time step *dt
          DemCount[1,1] = 1 # fill the first row
          SpSumPar <- SumPar[which(SumPar$species == x),] # select the right species df
          
          if(dim(SpSumPar)[1] != 1)
          {
            for (j in 2:dim(SpSumPar)[1]) # along the df
              DemCount[ceiling(SpSumPar$pop[j]*dt),1] = DemCount[ceiling(SpSumPar$pop[j-1]*dt),1] + 1 # total number of phenotypes at that time
            for (j in 1:dim(SpSumPar)[1]) # along the df 
              if (SpSumPar$exit[j] != timeMax) # do not take into account extinction at the last step (because its not)
                DemCount[ceiling(SpSumPar$exit[j]*dt),2] = DemCount[ceiling(SpSumPar$exit[j]*dt),2] -1 # an extinction happened at that time
              
              # I have a matrix with holes, need to fill them
              ExCount = 0 
              for (i in 2:dim(DemCount)[1])
              {
                if (DemCount[i,1] == 0) DemCount[i,1] = DemCount[i-1,1]
                if (DemCount[i,2] != 0)
                {
                  DemCount[i,2] = ExCount + abs(DemCount[i,2])
                  ExCount = abs(DemCount[i,2]) } else  {
                    DemCount[i,2] =  ExCount 
                  }
              }
              DemCount <- as.data.frame(DemCount)
              colnames(DemCount) <- c("pop","exit")
              DemCount$time <- as.numeric(rownames(DemCount))
              DemCount$species <- x
              DemList[[x]] <- DemCount
          }
        }
        scenarioList[[listPosition]] <- do.call(rbind,lapply(DemList,function(x)x)) # store the list as only one dataframe (with the species identity)
      }
      dataList[[dirContentIdx]] <- scenarioList
    }
  }
  
  # do an average across simulation
  plot_dat <- data.frame(dataList[[1]][[1]]$time,dataList[[1]][[1]]$species,
                         apply(do.call(cbind,lapply(dataList, function(x) x[[1]]$pop)),1,mean),
                         apply(do.call(cbind,lapply(dataList, function(x) x[[1]]$exit)),1,mean),
                         apply(do.call(cbind,lapply(dataList, function(x) x[[2]]$pop)),1,mean),
                         apply(do.call(cbind,lapply(dataList, function(x) x[[2]]$exit)),1,mean))
  colnames(plot_dat) <- c("time","species","popN","exitN","popF","exitF")
  
  # PLOTTING TIME
  # color gradient
  colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
  colGrad <- colfunc(length(SpIdx))
  names(colGrad) <- SpIdx
  if (length(SpIdx) == 9) colLine <- c("solid","dashed","solid","dashed","solid","longdash","solid","longdash","solid") else colLine = rep("solid",3) # 3,6 and 8 are dashed
  
  
  
  p <- ggplot(plot_dat) +
    geom_line(aes(x = time, y = popN-exitN, color = as.factor(species)))+
    geom_line(aes(x = time, y = popF-exitF,color = as.factor(species)), linetype = "dashed")+
    scale_x_continuous(name = "Time (yr)") +
    scale_y_continuous(name = "Number of phenotypes")+
    scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    #scale_linetype_manual(name = "Fished") +
    theme(legend.title=element_blank(),
          legend.position=c(0.19,0.95),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    ggtitle("Variation of phenotype's number throughout the simulation")
  
  if(returnData) return(plot_dat) else return(p)
  
}
