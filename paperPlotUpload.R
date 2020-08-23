### Script grouping all the figures of Interacting forces of predation and fishing affect species maturation size
### Figures do not need raw data but respective .rds of processed data (lighter)
require(tidyverse)
require(gridExtra)
require(grid)
# devtools::install_github("zeehio/facetscales")
require(facetscales)

### initialise this for all plots ###
no_sp = 9
SpIdx <- seq(1,no_sp)
colfunc <- colorRampPalette(c("firebrick3","darkturquoise", "orange"))
colGrad <- colfunc(length(SpIdx))
names(colGrad) <- SpIdx
colLine <- c("solid","dashed")
names(colLine) <- c("un-fished","fished")
fishStart <- 3000
# for multiplots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

### Figure 2 - Biomass and trait variations ----------
# scaling for easier comparison
windowTrait = c(-100,160)

### PANNEL a) - Biomass variation without predation

plot_dat <- readRDS(file = "data/BiomassNoPred.rds")

# delete first time step of every species
plot_dat <- lapply(plot_dat,FUN = function(x) x[-c(seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp),
                                                   (seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp)+1),
                                                   (seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp)+2),
                                                   (seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp)+3)),])

# prepare the data for the ribbon (stadndard deviation)
g1 <- ggplot(plot_dat[[1]])+
  geom_smooth(aes(x = time, y = meanSp-sdSp, group = as.factor(sp))) +
  geom_smooth(aes(x= time, y = meanSp+sdSp, group = as.factor(sp))) 
gg1 <- ggplot_build(g1)
dfRibbon1 <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y, group = gg1$data[[1]]$group) #and extract the smooth data for the ribbon

g2 <- ggplot(plot_dat[[2]])+
  geom_smooth(aes(x = time, y = meanSp-sdSp, group = as.factor(sp))) +
  geom_smooth(aes(x= time, y = meanSp+sdSp, group = as.factor(sp))) 
gg2 <- ggplot_build(g2)
dfRibbon2 <- data.frame(x = gg2$data[[1]]$x, ymin = gg2$data[[1]]$y, ymax = gg2$data[[2]]$y, group = gg2$data[[1]]$group) #and extract the smooth data for the ribbon
dfRibbon2 <- dfRibbon2[which(dfRibbon2$x >= fishStart),]

#improving main data
plot_dat[[1]]$fisheries <- "un-fished"
plot_dat[[2]]$fisheries <- "fished"
plot_dat <- rbind(plot_dat[[1]],plot_dat[[2]])
colnames(plot_dat)[2] <- "species"
plot_dat$size <- NA
for(i in 1:dim(plot_dat)[1])
{
  if(plot_dat$species[i]  <=3) plot_dat$size[i] = "small"
  else if(plot_dat$species[i] >= 7) plot_dat$size[i] = "large"
  else plot_dat$size[i] = "medium"
}
plot_dat$species <- as.factor(plot_dat$species)

#improving ribbon data
dfRibbon1$size <- NA
for(i in 1:dim(dfRibbon1)[1])
{
  if(dfRibbon1$group[i]  <=3) dfRibbon1$size[i] = "small"
  else if(dfRibbon1$group[i] >= 7) dfRibbon1$size[i] = "large"
  else dfRibbon1$size[i] = "medium"
}

dfRibbon2$size <- NA
for(i in 1:dim(dfRibbon2)[1])
{
  if(dfRibbon2$group[i]  <=3) dfRibbon2$size[i] = "small"
  else if(dfRibbon2$group[i] >= 7) dfRibbon2$size[i] = "large"
  else dfRibbon2$size[i] = "medium"
}

p1 <- ggplot(plot_dat)+
  geom_line(aes(x = time, y = meanSp, color = species, linetype = fisheries))+
  facet_grid(size~.) +
  geom_ribbon(data = dfRibbon1, aes(x = x, ymin = ymin, ymax = ymax, group = group), fill = "grey", alpha = 0.2)+
  # geom_line(data = plot_dat[[2]],aes(x = time, y = meanSp, color = as.factor(sp)),linetype = "dashed")+
  geom_ribbon(data = dfRibbon2, aes(x = x, ymin = ymin, ymax = ymax, group = group), fill = "red", alpha = 0.06)+
  scale_x_continuous(name = NULL, limits = c(NA,NA))+
  scale_y_continuous(name = expression(paste("Biomass density in g.m"^"-3")), limits = c(NA,NA), trans = "log10", breaks = seq(0.6,2,0.2))+#, breaks = seq(0.1,10,0.2))+
  geom_vline(xintercept = fishStart, linetype = "dashed", alpha = .5, size = .5) +
  scale_color_manual(name = "Species", values = colGrad)+
  scale_linetype_manual(name = "Fisheries", values = colLine)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow=1)) +
  ggtitle("(a) Without predation")

## PANNEL b) - Trait variation without predation 

traitDF <- readRDS(file = "data/TraitNoPred.rds")

# converts to %
traitDF[[1]]$percentMean <- traitDF[[1]]$percentMean*100
traitDF[[2]]$percentMean <- traitDF[[2]]$percentMean*100
traitDF[[1]]$sd <- traitDF[[1]]$sd*100
traitDF[[2]]$sd <- traitDF[[2]]$sd*100

# do the trait average in the traitDF and sd
plot_dat <- list()
for (type in c(1,2))
{
  a<- traitDF[[type]]
  a$group <- sapply(a$group, function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  a<-a[order(a$time),]
  
  myMat <- matrix(NA, nrow = length(unique(a$time))*length(unique(a$group)), ncol = 4, dimnames = list(NULL, c("time","percentMean","sd","group")))
  counter = 1
  for (no_sp in unique(a$group))
  {
    tempSp <- a[which(a$group == no_sp),]
    for (x in unique(tempSp$time))
    {
      temp_dat <- tempSp[which(tempSp$time == x),]
      meanX <- mean(temp_dat$percentMean)
      sdX <- sd(temp_dat$percentMean)
      myMat[counter,] <- c(x,meanX,sdX,no_sp)
      counter = counter +1
    }
  }
  plot_dat[[type]] <- as.data.frame(myMat) # this a 3 col data frame with the averaged trait values of the species across runs
}

temp <- data.frame("time" = rep(0,no_sp), "percentMean" = rep(0,no_sp), "sd" = rep(0,no_sp),"group" = seq(1,no_sp))
plot_dat[[1]] <- rbind(temp,plot_dat[[1]])

# force smoothing to pass by 0
pTemp <- ggplot(plot_dat[[1]])+ # do the smoothing
  stat_smooth(aes(x=time,y=percentMean, group = group, color = as.factor(group), linetype = "un-fished"), method = "loess", span = 0.15, se = F, size = 0.5) +
  stat_smooth(data = plot_dat[[2]], aes(x=time,y=percentMean, group = group, color = as.factor(group) , linetype = "fished"),method = "loess", span = 0.15, se = F, size = 0.5)

gg3 <- ggplot_build(pTemp)
# change value for un fished
a <- gg3[[1]][[1]]
a$y[which(a$x == 0)] <- 0
gg3[[1]][[1]] <- a
# and fished
b <- gg3[[1]][[2]]
b$y[which(a$x == 0)] <- 0
gg3[[1]][[2]] <- b

p2 <- ggplot(gg3[[1]][[1]])+
  geom_line(aes(x=x,y=y, color = as.factor(group), linetype = "un-fished"), size = 0.5) +
  geom_line(data = gg3[[1]][[2]], aes(x=x,y=y, color = as.factor(group), linetype = "fished"), size = 0.5) +
  scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
  scale_y_continuous(name = "Trait difference in %", limits = windowTrait, breaks = seq(-100,150,50))+
  geom_vline(xintercept = fishStart, linetype = "dashed") +
  scale_color_manual(name = "Species", values = colGrad)+
  scale_linetype_manual(name = "Fisheries", values = colLine)+
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="none",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow=2),
         linetype = guide_legend(order = 2,override.aes = list(colour = "black")))+
  ggtitle("(b)")

# PANNEL c) - Biomass variation with predation
plot_dat <- readRDS(file = "data/BiomassPred.rds")


# delete first 4 time step of every species (for smoothing)
plot_dat <- lapply(plot_dat,FUN = function(x) x[-c(seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp),
                                                   (seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp)+1),
                                                   (seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp)+2),
                                                   (seq(1,dim(x)[1],dim(plot_dat[[1]])[[1]]/no_sp)+3)),])

# prepare the data for the ribbon (stadndard deviation)
g1 <- ggplot(plot_dat[[1]])+
  geom_smooth(aes(x = time, y = meanSp-sdSp, group = as.factor(sp))) +
  geom_smooth(aes(x= time, y = meanSp+sdSp, group = as.factor(sp))) 
gg1 <- ggplot_build(g1)
dfRibbon1 <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y, group = gg1$data[[1]]$group) #and extract the smooth data for the ribbon

g2 <- ggplot(plot_dat[[2]])+
  geom_smooth(aes(x = time, y = meanSp-sdSp, group = as.factor(sp))) +
  geom_smooth(aes(x= time, y = meanSp+sdSp, group = as.factor(sp))) 
gg2 <- ggplot_build(g2)
dfRibbon2 <- data.frame(x = gg2$data[[1]]$x, ymin = gg2$data[[1]]$y, ymax = gg2$data[[2]]$y, group = gg2$data[[1]]$group) #and extract the smooth data for the ribbon
dfRibbon2 <- dfRibbon2[which(dfRibbon2$x >= fishStart),]

# polishing data
plot_dat[[1]]$fisheries <- "un-fished"
plot_dat[[2]]$fisheries <- "fished"
plot_dat <- rbind(plot_dat[[1]],plot_dat[[2]])
colnames(plot_dat)[2] <- "species"

plot_dat$size <- NA
for(i in 1:dim(plot_dat)[1])
{
  if(plot_dat$species[i]  <=3) plot_dat$size[i] = "small"
  else if(plot_dat$species[i] >= 7) plot_dat$size[i] = "large"
  else plot_dat$size[i] = "medium"
}
plot_dat$species <- as.factor(plot_dat$species)

# polishing ribbon
dfRibbon1$size <- NA
for(i in 1:dim(dfRibbon1)[1])
{
  if(dfRibbon1$group[i]  <=3) dfRibbon1$size[i] = "small"
  else if(dfRibbon1$group[i] >= 7) dfRibbon1$size[i] = "large"
  else dfRibbon1$size[i] = "medium"
}

dfRibbon2$size <- NA
for(i in 1:dim(dfRibbon2)[1])
{
  if(dfRibbon2$group[i]  <=3) dfRibbon2$size[i] = "small"
  else if(dfRibbon2$group[i] >= 7) dfRibbon2$size[i] = "large"
  else dfRibbon2$size[i] = "medium"
}

p3 <- ggplot(plot_dat)+
  geom_line(aes(x = time, y = meanSp, color = species, linetype = fisheries))+
  facet_grid(size~.) +
  geom_ribbon(data = dfRibbon1, aes(x = x, ymin = ymin, ymax = ymax, group = group), fill = "grey", alpha = 0.2)+
  geom_ribbon(data = dfRibbon2, aes(x = x, ymin = ymin, ymax = ymax, group = group), fill = "red", alpha = 0.06)+
  scale_x_continuous(name = NULL, limits = c(NA,NA))+
  scale_y_continuous(name = NULL, limits = window, trans = "log10", breaks = seq(0.01,0.1,0.02))+
  geom_vline(xintercept = fishStart, linetype = "dashed", alpha = .5, size = .5) +
  scale_color_manual(name = "Species", values = colGrad)+
  scale_linetype_manual(name = "Fisheries", values = colLine)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow=1)) +
  ggtitle("(c) With predation")

## PANNEL d) - Trait variation with predation 

traitDF <- readRDS(file = "data/TraitPred.rds")
# converts to %
traitDF[[1]]$percentMean <- traitDF[[1]]$percentMean*100
traitDF[[2]]$percentMean <- traitDF[[2]]$percentMean*100
traitDF[[1]]$sd <- traitDF[[1]]$sd*100
traitDF[[2]]$sd <- traitDF[[2]]$sd*100

# do the trait average in the traitDF and sd
plot_dat <- list()
for (type in c(1,2))
{
  a<- traitDF[[type]]
  a$group <- sapply(a$group, function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  a<-a[order(a$time),]
  
  myMat <- matrix(NA, nrow = length(unique(a$time))*length(unique(a$group)), ncol = 4, dimnames = list(NULL, c("time","percentMean","sd","group")))
  counter = 1
  for (no_sp in unique(a$group))
  {
    tempSp <- a[which(a$group == no_sp),]
    for (x in unique(tempSp$time))
    {
      temp_dat <- tempSp[which(tempSp$time == x),]
      meanX <- mean(temp_dat$percentMean)
      sdX <- sd(temp_dat$percentMean)
      myMat[counter,] <- c(x,meanX,sdX,no_sp)
      counter = counter +1
    }
  }
  plot_dat[[type]] <- as.data.frame(myMat) # this a 3 col data frame with the averaged trait values of the species across runs
}

temp <- data.frame("time" = rep(0,no_sp), "percentMean" = rep(0,no_sp), "sd" = rep(0,no_sp),"group" = seq(1,no_sp))
plot_dat[[1]] <- rbind(temp,plot_dat[[1]])

# force smoothing to pass by 0
pTemp <- ggplot(plot_dat[[1]])+ # do the smoothing
  stat_smooth(aes(x=time,y=percentMean, group = group, color = as.factor(group), linetype = "un-fished"), method = "loess", span = 0.15, se = F, size = 0.5) +
  stat_smooth(data = plot_dat[[2]], aes(x=time,y=percentMean, group = group, color = as.factor(group) , linetype = "fished"),method = "loess", span = 0.15, se = F, size = 0.5)

gg3 <- ggplot_build(pTemp)
# change value for un fished
a <- gg3[[1]][[1]]
a$y[which(a$x == 0)] <- 0
gg3[[1]][[1]] <- a
# and fished
b <- gg3[[1]][[2]]
b$y[which(a$x == 0)] <- 0
gg3[[1]][[2]] <- b


p4 <- ggplot(gg3[[1]][[1]])+
  geom_line(aes(x=x,y=y, color = as.factor(group), linetype = "un-fished"), size = 0.5) +
  geom_line(data = gg3[[1]][[2]], aes(x=x,y=y, color = as.factor(group), linetype = "fished"), size = 0.5) +
  scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
  scale_y_continuous(name = NULL, limits = windowTrait, breaks = seq(-100,150,50))+
  geom_vline(xintercept = fishStart, linetype = "dashed") +
  scale_color_manual(name = "Species", values = colGrad)+
  scale_linetype_manual(name = "Fisheries", values = colLine)+
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow=1),
         linetype = guide_legend(order = 2,override.aes = list(colour = "black")))+
  ggtitle("(d)")

### Multiplot
mylegend<-g_legend(p1)


p5 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p3 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               p4 + theme(legend.position="none"),
                               nrow=2),
                    mylegend, nrow=2,heights=c(9.5,0.5))

ggsave(p5,filename = "figures/Fig2.eps",width = 20,height = 20,units = "cm", device=cairo_ps)


### Figure 3 - Statistic analysis pannel ------------

edf <- readRDS("data/anovaDat.rds")

#backtransform to orginal eta scale
edf$fit <-exp(edf$fit)
edf$lower <-exp(edf$lower)
edf$upper <-exp(edf$upper)


## plot effects
edf$etaMean <-edf$fit

# formating for nice plot
edf$predation <- as.character(edf$predation)
edf$fisheries <- as.character(edf$fisheries)
edf$predation[which(edf$predation == 0)] <- "a) without predation"
edf$predation[which(edf$predation == 1)] <- "b) with predation"
edf$fisheries[which(edf$fisheries == 0)] <- "un-fished"
edf$fisheries[which(edf$fisheries == 1)] <- "fished"


fie <- readRDS(file = "data/anovaWeighted.rds") #weighted species mean (weighted by phenotype abundance) at year 6000
fie$run <- fie$simulation
fie$species <- as.factor(fie$species)
fie$predation[which(fie$predation == 0)] <- "a) without predation"
fie$predation[which(fie$predation == 1)] <- "b) with predation"
fie$fisheries[which(fie$fisheries == 0)] <- "un-fished"
fie$fisheries[which(fie$fisheries == 1)] <- "fished"


p <- ggplot(edf,aes(x=species,y=etaMean))  +
  geom_point(data=fie,color="dark grey",size=0.5,alpha=0.8)  +   
  facet_wrap(~predation)+
  geom_pointrange(data=edf,aes(color=fisheries,ymin=lower,ymax=upper),position=position_dodge(width=.1),alpha=0.8) +
  scale_color_manual(values = c("red","black")) +
  xlab("Species") + ylab("Eta") + 
  geom_hline(yintercept=0.25, linetype="dashed", color="grey")  +
  ylab(bquote(eta))  +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(),#element_text(c("no predation","predation")),
        legend.key = element_rect(fill = "white"))

ggsave(p, file = "figures/Fig3.eps", units = "cm", width = 20, height = 15, device=cairo_ps)

### Figure 4 - Fitness -----

iTime = 1
fitnessNorm <- readRDS("data/fitnessNormal.rds")
fitnessFish <- readRDS("data/fitnessFisheries.rds")
fitnessNormNoI <- readRDS("data/fitnessNormNoInter.rds")
fitnessFishNoI <- readRDS("data/fitnessFishNoInter.rds")

myDataN <- fitnessNorm[,c("trait","species","sim",iTime)]
colnames(myDataN) <- c("trait","species","sim","fitness")
myDataN <- myDataN[-which(myDataN$fitness == 0),]
myDataN$scenario <- "un-fished"
myDataN$interaction <- "predation"

myDataF <- fitnessFish[,c("trait","species","sim",iTime)]
colnames(myDataF) <- c("trait","species","sim","fitness")
myDataF <- myDataF[-which(myDataF$fitness == 0),]
myDataF$scenario <- "fished"
myDataF$interaction <- "predation"

myDataNnoI <- fitnessNormNoI[,c("trait","species","sim",iTime)]
colnames(myDataNnoI) <- c("trait","species","sim","fitness")
myDataNnoI <- myDataNnoI[-which(myDataNnoI$fitness == 0),]
myDataNnoI$scenario <- "un-fished"
myDataNnoI$interaction <- "no predation"

myDataFnoI <- fitnessFishNoI[,c("trait","species","sim",iTime)]
colnames(myDataFnoI) <- c("trait","species","sim","fitness")
myDataFnoI <- myDataFnoI[-which(myDataFnoI$fitness == 0),]
myDataFnoI$scenario <- "fished"
myDataFnoI$interaction <- "no predation"

plot_dat <- rbind(myDataN,myDataF,myDataNnoI,myDataFnoI)
plot_dat$interaction_f = factor(plot_dat$interaction, levels=c("no predation","predation")) # to define column order

# plots without facet for paper format
stripNames <- c("no predation" = "a) Without predation", "predation" = "b) With predation")
p1 <- ggplot(plot_dat %>% filter(species == 1)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 1") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free", labeller = as_labeller(stripNames)) +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(),#element_text(c("no predation","predation")),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p2 <- ggplot(plot_dat %>% filter(species == 2)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 2") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)


p3 <- ggplot(plot_dat %>% filter(species == 3)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 3") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p4 <- ggplot(plot_dat %>% filter(species == 4)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 4") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p5 <- ggplot(plot_dat %>% filter(species == 5)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 5") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p6 <- ggplot(plot_dat %>% filter(species == 6)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 6") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p7 <- ggplot(plot_dat %>% filter(species == 7)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 7") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free",labeller = as_labeller(stripNames)) +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p8 <- ggplot(plot_dat %>% filter(species == 8)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 8") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p9 <- ggplot(plot_dat %>% filter(species == 9)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 9") +
  scale_x_continuous(name = "Maturation size in g") +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

mylegend<-g_legend(p1)

p10 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                p3 + theme(legend.position="none"),
                                p4 + theme(legend.position="none"),
                                p5 + theme(legend.position="none"),
                                p6 + theme(legend.position="none"),
                                p7 + theme(legend.position="none"),
                                p8 + theme(legend.position="none"),
                                p9 + theme(legend.position="none"),
                                nrow=9),
                    mylegend, nrow=2,heights=c(9.5,.5),
                    left=textGrob("Fitness", rot = 90, vjust = 1))

ggsave(p10,filename = "figures/Fig4.eps",width = 15,height = 45,units = "cm",device=cairo_ps)



### Figure 5 - Effort sensitivity -----------

plot_dat <- readRDS("data/traitEffort.rds")

# make a pannel with small, medium, large
plot_dat$size <- NA
for(i in 1:dim(plot_dat)[1])
{
  if(plot_dat$group[i]  <=3) plot_dat$size[i] = "small"
  else if(plot_dat$group[i] >= 7) plot_dat$size[i] = "large"
  else plot_dat$size[i] = "medium"
}

plot_dat$group <- as.factor(plot_dat$group)
colnames(plot_dat)[4] <- "species"
plot_dat$time <- NULL

plot_dat$percentMean <- plot_dat$percentMean *100
plot_dat$sd <- plot_dat$sd *100

colGrad <- colfunc(length(SpIdx))
names(colGrad) <- SpIdx

p1 <-ggplot(plot_dat)+
  geom_pointrange(aes(x = effort, y = percentMean,ymin = percentMean - sd, ymax = percentMean + sd, group = species, color = species)) +
  geom_line(aes(x = effort, y = percentMean, group = species, color = species)) +
  scale_x_continuous(name = "Fishing mortality", breaks = seq(0,1,.1))+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  facet_grid(size~., scales = "free") +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="right",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p1, file = "figures/Fig5.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure 6 - winf/wmat --------------
myDat <- readRDS("data/maxSize.rds")

# Averaging across species/sim/scenario
plot_dat <- NULL
for(iSpecies in unique(myDat$species))
{
  for(iScenario in unique(myDat$scenario))
  {
    for(iFish in unique(myDat$fisheries))
    {  
      temp <- filter(myDat, species == iSpecies & scenario == iScenario & fisheries == iFish)
      wmatMean <- mean(temp$maturationSize)
      wmatSd <- sd(temp$maturationSize)
      wmaxMean <- mean(temp$maxSizeTheoretical)
      wmaxSd <- sd(temp$maxSizeTheoretical) 
      plot_dat <- rbind(plot_dat,c(iSpecies,iScenario,iFish,wmatMean,wmatSd,wmaxMean,wmaxSd))
    }
  }
}
# formating the data
plot_dat <- as.data.frame(plot_dat)
colnames(plot_dat) <- c("species","scenario","fisheries","wmatMean","wmatSd","wmaxMean","wmaxSd")
plot_dat$wmatMean <- as.numeric(as.character(plot_dat$wmatMean))
plot_dat$wmatSd <- as.numeric(as.character(plot_dat$wmatSd))
plot_dat$wmaxMean <- as.numeric(as.character(plot_dat$wmaxMean))
plot_dat$wmaxSd <- as.numeric(as.character(plot_dat$wmaxSd))

plot_dat$scenario <- as.character(plot_dat$scenario)
plot_dat$scenario[which(plot_dat$scenario == "Sim9")] <- "predation"
plot_dat$scenario[which(plot_dat$scenario == "SimNoInter9")] <- "no predation"

plot_dat <- filter(plot_dat, fisheries != "init")
plot_dat$fisheries <- as.character(plot_dat$fisheries)
plot_dat$fisheries[which(plot_dat$fisheries == "normal")] <- "un-fished"
plot_dat$fisheries[which(plot_dat$fisheries == "fisheries")] <- "fished"

## add empirical data to the plot
empDat <- readxl::read_excel("data/OlssonGislasonDATA.xlsx")
empDat <- as.data.frame(empDat)

empDat2 <- readxl::read_excel("data/marineSizeData.xlsx")
empDat2 <- as.data.frame(empDat2)
colnames(empDat2)[1] <- "species"
empDat2$...5 <- NULL
colnames(empDat2)[2] <- "Lmax"
empDat2$wmax <- empDat2$a*empDat2$Lmax^(empDat2$b)
empDat2$wmat <- empDat2$a*empDat2$Lmat^(empDat2$b)

empDat3 <- readxl::read_excel("data/smallFish.xlsx")
empDat3 <- as.data.frame(empDat3)
# merge empirical data
empDat1 <- empDat[,c("Wmat","Winf")]
empDat2 <- empDat2[,c("wmat","wmax")]
colnames(empDat2) <- colnames(empDat1)
colnames(empDat3) <- colnames(empDat1)
empDat <- rbind(empDat1,empDat2,empDat3)

p <- ggplot(filter(plot_dat)) +  
  geom_pointrange(aes(x = wmaxMean, y = wmatMean,ymin = wmatMean - wmatSd, ymax = wmatMean + wmatSd, color = fisheries, shape = scenario), alpha = .8, size = 1) +
  geom_errorbarh(aes(y = wmatMean,xmin = wmaxMean - wmaxSd, xmax = wmaxMean + wmaxSd, height = 0, color = fisheries)) +
  geom_point(data = empDat, aes(x = Winf, y = Wmat), color = "grey", shape = "*", size = 10) +
  scale_x_continuous(trans = "log10", name = "Asymptotic size in g") +
  scale_y_continuous(trans = "log10", name = "Maturation size in g") +
  scale_color_manual(values = c("red","black")) +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p, file = "figures/Fig6.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S1 - eta sensitivity --------------
plot_dat <- readRDS("data/etaSensitivity.rds")
plot_dat <- filter(plot_dat, time == max(plot_dat$time))
plot_dat$time <- NULL

# conversion eta to % change
plot_dat$trait <- as.numeric(as.character(plot_dat$trait))
plot_dat$mean <- (plot_dat$mean/plot_dat$trait-1)*100
plot_dat$sd <- (plot_dat$sd/plot_dat$trait) *100


scales_y <- list(
  large = scale_y_continuous(name = "Trait difference in %",limits = c(-100,100)),
  medium = scale_y_continuous(name = "Trait difference in %",limits = c(-100,200)),
  small = scale_y_continuous(name = "Trait difference in %",limits = c(-100,100)))


p2 <- ggplot(plot_dat) +
  geom_pointrange(aes(x = as.factor(trait), y = mean, ymax = mean+sd, ymin = mean-sd, group = species, color = species),size = .3) +
  facet_grid_sc(size~., scales = list(y = scales_y)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = "Eta initial value")+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  guides(colour = guide_legend(nrow = 1)) +
  ggtitle(NULL)

ggsave(p2,filename = "figures/FigS1.eps",width = 20,height = 20,units = "cm",device=cairo_ps)

### Figure S2 - zeta sensitivity ----------------------
plot_dat <- readRDS("data/mAmplitudeSensitivity.rds")
plot_dat <- filter(plot_dat, time == max(plot_dat$time))
plot_dat$time <- NULL

plot_dat$size <- NA
for(i in 1:dim(plot_dat)[1])
{
  if(plot_dat$species[i]  <=3) plot_dat$size[i] = "small"
  else if(plot_dat$species[i] >= 7) plot_dat$size[i] = "large"
  else plot_dat$size[i] = "medium"
}

plot_dat$species <- as.factor(plot_dat$species)

# conversion eta to % change
plot_dat$mean <- (plot_dat$mean/0.25-1)*100
plot_dat$sd <- (plot_dat$sd/0.25) *100

plot_dat$type <- as.character(plot_dat$type)
plot_dat$type[which(plot_dat$type == "un-fished")] <- "un-fished: predation"
plot_dat$type[which(plot_dat$type == "fished")] <- "fished: predation"

p <- ggplot(plot_dat) +
  geom_pointrange(aes(x = trait, y = mean, ymax = mean+sd, ymin = mean-sd, group = species, color = species),size = .3) +
  facet_grid(size~type) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = "trait Magnitude")+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position="right",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p, file = "figures/FigS2.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S3 - chi sensitivity ---------------
plot_dat <- readRDS("data/mAmplitudeSensitivity.rds")
plot_dat <- filter(plot_dat, time == max(plot_dat$time))
plot_dat$time <- NULL

plot_dat$size <- NA
for(i in 1:dim(plot_dat)[1])
{
  if(plot_dat$species[i]  <=3) plot_dat$size[i] = "small"
  else if(plot_dat$species[i] >= 7) plot_dat$size[i] = "large"
  else plot_dat$size[i] = "medium"
}

plot_dat$species <- as.factor(plot_dat$species)
# conversion eta to % change
plot_dat$mean <- (plot_dat$mean/0.25-1)*100
plot_dat$sd <- (plot_dat$sd/0.25) *100

plot_dat$type <- as.character(plot_dat$type)
plot_dat$type[which(plot_dat$type == "un-fished")] <- "un-fished: predation"
plot_dat$type[which(plot_dat$type == "fished")] <- "fished: predation"

p <- ggplot(plot_dat) +
  geom_pointrange(aes(x = trait, y = mean, ymax = mean+sd, ymin = mean-sd, group = species, color = species),size = .3) +
  facet_grid(size~type) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = expression(paste("1000.",chi,sep="")))+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="right",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p, file = "figures/FigS3.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S4 - initial phen fraction abundance sensitivity ------------

plot_dat <- readRDS("data/initPhenSensitivity.rds")
plot_dat <- filter(plot_dat, time == max(plot_dat$time))
plot_dat$time <- NULL

# conversion eta to % change
plot_dat$mean <- (plot_dat$mean/0.25-1)*100
plot_dat$sd <- (plot_dat$sd/0.25) *100

p <- ggplot(plot_dat) +
  geom_pointrange(aes(x = trait, y = mean, ymax = mean+sd, ymin = mean-sd, group = species, color = species),size = .3) +
  facet_grid(size~.) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = "Initial phenotype abundance fraction")+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="right",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p, file = "figures/FigS4.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S5 - sigma sensitivity --------------------

plot_dat <- readRDS("data/sigmaSensitivity.rds")
plot_dat <- filter(plot_dat, time == max(plot_dat$time))
plot_dat$time <- NULL

# conversion eta to % change
plot_dat$mean <- (plot_dat$mean/0.25-1)*100
plot_dat$sd <- (plot_dat$sd/0.25) *100

plot_dat$type <- as.character(plot_dat$type)
plot_dat$type[which(plot_dat$type == "un-fished")] <- "un-fished: predation"
plot_dat$type[which(plot_dat$type == "fished")] <- "fished: predation"

p <- ggplot(plot_dat) +
  geom_pointrange(aes(x = trait, y = mean, ymax = mean+sd, ymin = mean-sd, group = species, color = species),size = .3) +
  facet_grid(size~type) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = expression(sigma))+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="right",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p, file = "figures/FigS5.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S6 - beta sensitivity ----------------

plot_dat <- readRDS("data/betaSensitivity.rds")
plot_dat <- filter(plot_dat, time == max(plot_dat$time))
plot_dat$time <- NULL

# conversion eta to % change
plot_dat$mean <- (plot_dat$mean/0.25-1)*100
plot_dat$sd <- (plot_dat$sd/0.25) *100

plot_dat$type <- as.character(plot_dat$type)
plot_dat$type[which(plot_dat$type == "un-fished")] <- "un-fished: predation"
plot_dat$type[which(plot_dat$type == "fished")] <- "fished: predation"

p <- ggplot(plot_dat) +
  geom_pointrange(aes(x = trait, y = mean, ymax = mean+sd, ymin = mean-sd, group = species, color = species),size = .3) +
  facet_grid(size~type) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = expression(beta))+
  scale_y_continuous(name = "Trait difference in %") +
  scale_color_manual(name = "species",values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"), 
        legend.position="right",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

ggsave(p, file = "figures/FigS6.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S7|8 - number of simulations un-fished|fished -----------

plot_dat <- readRDS("data/simNumberDat.rds")

# add grouping per size
plot_dat$size <- NA
for(i in 1:dim(plot_dat)[1])
  if(plot_dat$group[i] < 4) plot_dat$size[i] <- "small" else if (plot_dat$group[i] > 6) plot_dat$size[i] <- "large" else plot_dat$size[i] <- "medium"

plot_dat$group <- as.factor(plot_dat$group)
plot_dat$simulations <- as.factor(plot_dat$simulations)

# change y to proportional change
plot_dat$percentMean <- plot_dat$percentMean/100
plot_dat$sd<- plot_dat$sd/100


p1 <- ggplot(filter(plot_dat, fisheries == "un-fished"))+
  geom_line(aes(x=simulations, y= percentMean, color = group, group = group)) +
  geom_point(aes(x=simulations, y= percentMean, color = group, group = group)) +
  geom_errorbar(aes(x = simulations, y = percentMean,ymin = percentMean - sd, ymax = percentMean + sd, group = group, color = group), width = .1) +
  facet_grid(size~., scales = "free")+
  scale_x_discrete(name = "Simulation number")+
  scale_y_continuous(name = "Proportional change")+
  scale_color_manual(name = "Species", values = colGrad)+
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow=1))+
  ggtitle(NULL)

ggsave(p1, file = "figures/FigS7.eps", units = "cm", width = 20, height = 20,device=cairo_ps)



p2 <- ggplot(filter(plot_dat, fisheries == "fished"))+
  geom_line(aes(x=simulations, y= percentMean, color = group, group = group)) +
  geom_point(aes(x=simulations, y= percentMean, color = group, group = group)) +
  geom_errorbar(aes(x = simulations, y = percentMean,ymin = percentMean - sd, ymax = percentMean + sd, group = group, color = group), width = .1) +
  facet_grid(size~., scales = "free")+
  scale_x_discrete(name = "Simulation number")+
  scale_y_continuous(name = "Proportional change")+
  scale_color_manual(name = "Species", values = colGrad)+
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow=1))+
  ggtitle(NULL)

ggsave(p2, file = "figures/FigS8.eps", units = "cm", width = 20, height = 20,device=cairo_ps)
### Figure S9 - Fitness end of simulation --------------
iTime = 2500
fitnessNorm <- readRDS("data/fitnessNormal.rds")
fitnessFish <- readRDS("data/fitnessFisheries.rds")
fitnessNormNoI <- readRDS("data/fitnessNormNoInter.rds")
fitnessFishNoI <- readRDS("data/fitnessFishNoInter.rds")

myDataN <- fitnessNorm[,c("trait","species","sim",iTime)]
colnames(myDataN) <- c("trait","species","sim","fitness")
myDataN <- myDataN[-which(myDataN$fitness == 0),]
myDataN$scenario <- "un-fished"
myDataN$interaction <- "predation"

myDataF <- fitnessFish[,c("trait","species","sim",iTime)]
colnames(myDataF) <- c("trait","species","sim","fitness")
myDataF <- myDataF[-which(myDataF$fitness == 0),]
myDataF$scenario <- "fished"
myDataF$interaction <- "predation"

myDataNnoI <- fitnessNormNoI[,c("trait","species","sim",iTime)]
colnames(myDataNnoI) <- c("trait","species","sim","fitness")
myDataNnoI <- myDataNnoI[-which(myDataNnoI$fitness == 0),]
myDataNnoI$scenario <- "un-fished"
myDataNnoI$interaction <- "no predation"

myDataFnoI <- fitnessFishNoI[,c("trait","species","sim",iTime)]
colnames(myDataFnoI) <- c("trait","species","sim","fitness")
myDataFnoI <- myDataFnoI[-which(myDataFnoI$fitness == 0),]
myDataFnoI$scenario <- "fished"
myDataFnoI$interaction <- "no predation"

plot_dat <- rbind(myDataN,myDataF,myDataNnoI,myDataFnoI)
plot_dat$interaction_f = factor(plot_dat$interaction, levels=c("no predation","predation")) # to define column order

# plots without facet for paper format
stripNames <- c("no predation" = "a) Without predation", "predation" = "b) With predation")
p1 <- ggplot(plot_dat %>% filter(species == 1)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 1") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free", labeller = as_labeller(stripNames)) +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(),#element_text(c("no predation","predation")),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p2 <- ggplot(plot_dat %>% filter(species == 2)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 2") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)


p3 <- ggplot(plot_dat %>% filter(species == 3)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 3") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p4 <- ggplot(plot_dat %>% filter(species == 4)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 4") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p5 <- ggplot(plot_dat %>% filter(species == 5)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 5") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p6 <- ggplot(plot_dat %>% filter(species == 6)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 6") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p7 <- ggplot(plot_dat %>% filter(species == 7)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 7") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free",labeller = as_labeller(stripNames)) +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p8 <- ggplot(plot_dat %>% filter(species == 8)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 8") +
  scale_x_continuous(name = NULL) +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

p9 <- ggplot(plot_dat %>% filter(species == 9)) +
  geom_point(aes(x=trait,y=fitness, color = scenario), size = 0.5) +
  scale_y_continuous(trans = "log10", name = "Species 9") +
  scale_x_continuous(name = "Maturation size in g") +
  scale_color_manual(values = c("red","black")) +
  facet_grid(.~interaction_f, scales = "free") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.minor = element_line(colour = "grey92"),
        # legend.justification=c(1,1),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

mylegend<-g_legend(p1)

p10 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                p3 + theme(legend.position="none"),
                                p4 + theme(legend.position="none"),
                                p5 + theme(legend.position="none"),
                                p6 + theme(legend.position="none"),
                                p7 + theme(legend.position="none"),
                                p8 + theme(legend.position="none"),
                                p9 + theme(legend.position="none"),
                                nrow=9),
                    mylegend, nrow=2,heights=c(9.5,.5),
                    left=textGrob("Fitness", rot = 90, vjust = 1))

ggsave(p10,filename = "figures/FigS9.eps",width = 15,height = 45,units = "cm",device=cairo_ps)

### Figure S10 - phenotype's number through time -----------

# making data just right
plot_datPred <- readRDS(file = "data/PhenDataPred.rds")
plot_datPred <- plot_datPred[[1]]
plot_datPred$scenario <- "predation"
temp <- plot_datPred
plot_datPred <- plot_datPred[,-4]
plot_datPred$fisheries <- "un-fished"
colnames(plot_datPred)[3] <- "value"
temp$fisheries <- "fished"
temp <- temp[,-3]
colnames(temp)[3] <- "value"
plot_datPred <- rbind(plot_datPred,temp)

plot_datNoPred <- readRDS(file = "data/PhenDataNoPred.rds")
plot_datNoPred <- plot_datNoPred[[1]]
plot_datNoPred$scenario <- "no predation"
temp <- plot_datNoPred
plot_datNoPred <- plot_datNoPred[,-4]
plot_datNoPred$fisheries <- "un-fished"
colnames(plot_datNoPred)[3] <- "value"
temp$fisheries <- "fished"
temp <- temp[,-3]
colnames(temp)[3] <- "value"
plot_datNoPred <- rbind(plot_datNoPred,temp)

plot_dat <- rbind(plot_datPred,plot_datNoPred) 
plot_dat$species <- as.factor(plot_dat$species)

p <- ggplot(plot_dat) +
  stat_smooth(aes(x = time, y = value, color = species, linetype = fisheries), method = "loess", span = 0.15, se = F, size = 0.5)+
  facet_grid(scenario~.) +
  scale_x_continuous(name = "Time in years")+
  scale_y_continuous(name = "Number of phenotypes") +
  scale_color_manual(name = "Species", values = colGrad)+
  geom_vline(xintercept = 3000, linetype = "dashed") +
  scale_linetype_manual(name = "Fisheries", values = colLine)+
  theme(panel.background = element_rect(fill = "white", color = "black"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"))+ 
  guides(color = guide_legend(nrow=1)) +
  ggtitle(NULL)

ggsave(p, file = "figures/FigS10.eps", units = "cm", width = 20, height = 20,device=cairo_ps)

### Figure S11 - Feeding and mortality ---------------
plot_dat <- readRDS("data/foodData.rds")
w_inf <- c(10, 31, 100, 316, 1000, 3162, 10000, 31622, 100000)
vlines <- data.frame(xint = w_inf, grp = SpIdx)

p1 <- ggplot(plot_dat[[1]]) + 
  geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
  geom_hline(yintercept = plot_dat[[2]][1], linetype = "dashed", color = "red") +
  scale_x_log10(name = "Size") + 
  scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
  geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") +
  scale_color_manual(name = "Species", values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow = 1))+
  ggtitle("a)")

plot_dat2 <- readRDS("data/mortData.rds")

p2 <- ggplot(plot_dat2) + 
  geom_line(aes(x=w, y = value, colour = as.factor(Species))) + 
  scale_x_log10(name = "Size") + 
  scale_y_continuous(name = "Predation mortality")+
  geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") +
  scale_color_manual(name = "Species", values = colGrad)+
  theme(legend.title=element_text(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position="bottom",
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(nrow = 1))+
  ggtitle("b)")

mylegend<-g_legend(p1)

p3<- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2+ theme(legend.position="none"),
                                nrow=2),
                    mylegend, nrow=2,heights=c(9.5,.5))

ggsave(p3,filename = "figures/FigS11.eps",width = 20,height = 20,units = "cm",device=cairo_ps)
