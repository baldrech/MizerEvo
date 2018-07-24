# cohort plot -------------
fishStart = 3000
folder = paste(getwd(),"/SimNoInter92",sep="")
cohortSpan <- c(1000,2800,3200,5500)
myData <- plotCohortLong(folder = folder, cohortSpan = cohortSpan, t_steps = 30, returnData = T )
saveRDS(myData,file = "fitnessNoInter92Long.rds")
myData <- readRDS("fitnessNoInter92Long.rds")



### multiplot (species per time) ###

fishStart = 3000
folder = paste(getwd(),"/SimDefault",sep="")
no_sp = 9
plotSpList <- vector("list",no_sp)
counterSp = 0
for(iSpecies in seq(1:no_sp)) #c(1,6,9)) #
{
  counterSp = counterSp +1
  myMat <- as.data.frame(myData[[1]][[iSpecies]])
  myMat2 <- as.data.frame(myData[[2]][[iSpecies]])
  
  # then process
  myMat$species <- NULL
  myMat$group <- NULL
  myMat[myMat==0] <- NA
  
  # then process
  temp <- myMat2[,which(as.numeric(colnames(myMat2)) >= fishStart)]
  myMat2 <- cbind(temp,myMat2$w_mat)
  myMat2[myMat2==0] <- NA
  
  # change slightly the w_mat of the original species so they do not overlap across sim
  #myMat$w_mat[myMat$w_mat == myMat$w_mat[1]] <- sapply(myMat$w_mat[myMat$w_mat == myMat$w_mat[1]],function(x) x*rnorm(1,1,0.00001))
  
  plot_dat <- melt(myMat, id = "w_mat")
  plot_dat$variable <- as.numeric(as.character(plot_dat$variable))
  
  colnames(plot_dat) <- c("w_mat","Time","Fitness")
  
  plot_dat2 <- melt(myMat2, id = "myMat2$w_mat")
  plot_dat2$variable <- as.numeric(as.character(plot_dat2$variable))
  
  colnames(plot_dat2) <- c("w_mat","Time","Fitness")
  
  whichTime <- unique(plot_dat$Time)#[c(1,3,5,6,9,10,11)]
  plotTimeList <- vector("list",length(whichTime))
  counter = 0
  for(iTime in whichTime)
  {
    counter <- counter +1
    plot_dat_time <- plot_dat[which(plot_dat$Time == iTime),]
    plot_dat_time2 <- plot_dat2[which(plot_dat2$Time == iTime),]
    
    if (counter == 1) yname = paste("Species",iSpecies) else yname = NULL
    if (counterSp == no_sp) xname = paste(iTime,"years") else xname = NULL
    
    
    plotTimeList[[counter]] <- ggplot(plot_dat_time) +
      geom_point(aes(x=w_mat,y=Fitness)) +#,group = Time, color = Time), size = 1) +
      geom_point(data = plot_dat_time2, aes(x=w_mat,y=Fitness), color = "red") +
      scale_y_continuous(trans = "log10", name = yname) +
      scale_x_continuous(name = xname) +
      scale_colour_gradient2(low = "blue", mid = "orange", high = "black", midpoint = fishStart) +
      theme(legend.title=element_text(),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.justification=c(1,1),
            legend.position = "right",
            legend.key = element_rect(fill = "white"))+
      ggtitle(NULL)
    
  }
  plotSpList[[counterSp]] <- plotTimeList
}

png(filename=paste(folder,"/Fitness.png",sep=""), width = (8*(length(plotSpList[[1]]))), height = (8*(length(plotSpList))), units = "cm",res = 600)

grid.newpage() 
pushViewport(viewport(layout = grid.layout(nrow=length(plotSpList), ncol=length(plotSpList[[1]]), 
                                           widths = unit(8, "cm"), 
                                           heights = unit(8, "cm"))))

for(i in 1:length(plotSpList))
  for(j in 1:length(plotSpList[[i]]))
    print(plotSpList[[i]][[j]], vp = viewport(layout.pos.row = i, layout.pos.col = j))

dev.off()



# paste data together
myData1 <- readRDS("fitnessMat5.rds")
myData2 <- readRDS("fitnessMat6.rds")
myData3 <- readRDS("fitnessMat7.rds")

fitnessSpList <- vector("list",9)
fitnessList <- list(fitnessSpList,fitnessSpList)

for(iList in 1:length(myData1))
{
  for(jList in 1:length(myData1[[iList]]))
  {
    
    a1 <- myData1[[iList]][[jList]]
    a2 <- myData2[[iList]][[jList]]
    a3 <- myData3[[iList]][[jList]]
    
    spNames <- unique(c(rownames(a1),rownames(a2),rownames(a3)))
    timeNames <- sort(unique(c(colnames(a1),colnames(a2),colnames(a3))))
    timeNames <- c("500",timeNames[-10]) # move the 500 in the proper position
    
    fitnessMat <- matrix(data=0,ncol = length(timeNames),nrow = length(spNames),dimnames = list(spNames,timeNames))
    
    for (irow in spNames)
    {
      for (icol in timeNames) 
      {
        
        if (sum(sum(rownames(a1) %in% irow), sum(colnames(a1) %in% icol)) == 2)
          fitOpt = 1
        else if (sum(sum(rownames(a2) %in% irow), sum(colnames(a2) %in% icol)) == 2)
          fitOpt = 2
        else if (sum(sum(rownames(a3) %in% irow), sum(colnames(a3) %in% icol)) == 2)
          fitOpt = 3
        else
          fitOpt = 4
        
        
        switch(fitOpt,
               "1" = {fitnessMat[irow,icol] <- a1[irow,icol]},
               "2" = {fitnessMat[irow,icol] <- a2[irow,icol]},
               "3" = {fitnessMat[irow,icol] <- a3[irow,icol]},
               {})
        
        #cat(sprintf("irow = %s, icol = %s, fitOpt = %g and fitness = %g\n",irow,icol,fitOpt,fitnessMat[irow,icol]))
      }
    }
    
    fitnessList[[iList]][[jList]] <- fitnessMat
    
    
    
    
    
    
    
    
  }
}
saveRDS(fitnessList,"fitnessData.rds")
fitnessList <- readRDS("fitnessData.rds")
myData <- fitnessList



# I don't know
a <- get(load("SimDefault/init/run1/run.Rdata"))
b <- plotCohort(a,t_steps = 1, cohortSpan = c(1000,2000),returnData = T)
b$group <- 3
c <- list(b)
c[[1]]@group <- 3
# Let's play with the data
# myData is a list of normal fitness (1) et fisheries fitness(2)
# inside each list are 9 matrix, one per species, so I lost my multiple sims on the way


#plot_dat <- myData[[1]][[1]] # fitness of all the phenotypes across sims of species 1

fitness_dat <- as.data.frame(myData[[1]][[9]])
fitness_dat$species <- NULL

fitness_dat[fitness_dat==0] <- NA

# change slightly the w_mat of the original species so they do not overlap across sim
fitness_dat$w_mat[fitness_dat$w_mat == fitness_dat$w_mat[1]] <- sapply(fitness_dat$w_mat[fitness_dat$w_mat == fitness_dat$w_mat[1]],function(x) x*rnorm(1,1,0.00001))

plot_dat <- melt(fitness_dat, id = "w_mat")
plot_dat$variable <- as.numeric(as.character(plot_dat$variable))

ggplot(plot_dat) +
  geom_line(aes(x=as.numeric(variable),y=value,group = w_mat, color = w_mat)) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(name = "Time") +
  scale_colour_gradient2(low = "seagreen4", mid = "turquoise2", high = "dodgerblue4", midpoint = fitness_dat$w_mat[1]) +
  #geom_vline(xintercept = 3000, linetype = "dashed") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),
        legend.justification=c(1,1),
        legend.position = "right",
        legend.key = element_rect(fill = "white"))+
  ggtitle("Species 1")

ggsave("fitnessLine.png")
# fitness plot
# one time step and one species, plot relative fitness (to compare runs) 3 time steps



myData <- readRDS("fitnessMat5.rds")
# need to normalise fitness per runs
fitness_dat <- as.data.frame(myData[[1]][[1]])

myMat <- NULL
for (i in unique(fitness_dat$group))
{
  fitness_temp <- fitness_dat[which(fitness_dat$group == i),]
  myMat <- rbind(myMat,apply(fitness_temp,2,function(x) x/max(x)))
}
myMat <- as.data.frame(myMat)
# then process
myMat$species <- NULL
myMat$group <- NULL
myMat$w_mat <- fitness_dat$w_mat

myMat[myMat==0] <- NA

# change slightly the w_mat of the original species so they do not overlap across sim
myMat$w_mat[myMat$w_mat == myMat$w_mat[1]] <- sapply(myMat$w_mat[myMat$w_mat == myMat$w_mat[1]],function(x) x*rnorm(1,1,0.00001))

plot_dat <- melt(myMat, id = "w_mat")
plot_dat$variable <- as.numeric(as.character(plot_dat$variable))

ggplot(plot_dat) +
  #geom_point(aes(x=as.numeric(variable),y=value,group = w_mat, color = w_mat)) +
  geom_point(aes(x=w_mat,y=value,group = variable, color =variable)) +
  #scale_y_continuous(trans = "log10") +
  scale_x_continuous(name = "Time") +
  #scale_colour_gradient2(low = "seagreen4", mid = "turquoise2", high = "dodgerblue4", midpoint = myMat$w_mat[1]) +
  scale_colour_gradient2(low = "blue", mid = "orange", high = "black", midpoint = myMat$w_mat[1]) +
  #geom_vline(xintercept = 3000, linetype = "dashed") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),
        legend.justification=c(1,1),
        legend.position = "right",
        legend.key = element_rect(fill = "white"))+
  ggtitle("Species 1")
ggsave("fitnessDot.png")

myData <- readRDS("fitnessMat6.rds")

# need to normalise fitness per runs
sp=9
fitness_dat <- as.data.frame(myData[[1]][[sp]])
fitness_dat2 <- as.data.frame(myData[[2]][[sp]])

myMat <- NULL
for (i in unique(fitness_dat$group))
{
  fitness_temp <- fitness_dat[which(fitness_dat$group == i),]
  myMat <- rbind(myMat,apply(fitness_temp,2,function(x) x/sum(x)))
}
myMat <- as.data.frame(myMat)

myMat2 <- NULL
for (i in unique(fitness_dat2$group))
{
  fitness_temp <- fitness_dat[which(fitness_dat2$group == i),]
  myMat2 <- rbind(myMat2,apply(fitness_temp,2,function(x) x/sum(x)))
}
myMat2 <- as.data.frame(myMat2)

# then process
myMat$species <- NULL
myMat$group <- NULL
myMat$w_mat <- fitness_dat$w_mat
myMat[myMat==0] <- NA

# then process
temp <- myMat2[,which(as.numeric(colnames(myMat2)) >= fishStart)]
myMat2 <- cbind(temp,fitness_dat2$w_mat)
myMat2[myMat2==0] <- NA

# change slightly the w_mat of the original species so they do not overlap across sim
#myMat$w_mat[myMat$w_mat == myMat$w_mat[1]] <- sapply(myMat$w_mat[myMat$w_mat == myMat$w_mat[1]],function(x) x*rnorm(1,1,0.00001))

plot_dat <- melt(myMat, id = "w_mat")
plot_dat$variable <- as.numeric(as.character(plot_dat$variable))

colnames(plot_dat) <- c("w_mat","Time","Fitness")

plot_dat2 <- melt(myMat2, id = "fitness_dat2$w_mat")
plot_dat2$variable <- as.numeric(as.character(plot_dat2$variable))

colnames(plot_dat2) <- c("w_mat","Time","Fitness")



# or not
sp=9
myMat <- as.data.frame(myData[[1]][[sp]])
myMat2 <- as.data.frame(myData[[2]][[sp]])

# then process
myMat$species <- NULL
myMat$group <- NULL
myMat[myMat==0] <- NA

# then process
temp <- myMat2[,which(as.numeric(colnames(myMat2)) >= fishStart)]
myMat2 <- cbind(temp,myMat2$w_mat)
myMat2[myMat2==0] <- NA

# change slightly the w_mat of the original species so they do not overlap across sim
#myMat$w_mat[myMat$w_mat == myMat$w_mat[1]] <- sapply(myMat$w_mat[myMat$w_mat == myMat$w_mat[1]],function(x) x*rnorm(1,1,0.00001))

plot_dat <- melt(myMat, id = "w_mat")
plot_dat$variable <- as.numeric(as.character(plot_dat$variable))

colnames(plot_dat) <- c("w_mat","Time","Fitness")

plot_dat2 <- melt(myMat2, id = "myMat2$w_mat")
plot_dat2$variable <- as.numeric(as.character(plot_dat2$variable))

colnames(plot_dat2) <- c("w_mat","Time","Fitness")


ggplot(plot_dat) +
  #geom_point(aes(x=as.numeric(variable),y=value,group = w_mat, color = w_mat)) +
  geom_point(aes(x=w_mat,y=Fitness,group = Time, color = Time), size = 1) +
  #geom_point(data = plot_dat2, aes(x=w_mat,y=Fitness,group = Time, color = Time), shape = "+", size = 3) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(name = "Maturation size") +
  #scale_colour_gradient2(low = "seagreen4", mid = "turquoise2", high = "dodgerblue4", midpoint = myMat$w_mat[1]) +
  scale_colour_gradient2(low = "blue", mid = "orange", high = "black", midpoint = fishStart) +
  #geom_vline(xintercept = 3000, linetype = "dashed") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),
        legend.justification=c(1,1),
        legend.position = "right",
        legend.key = element_rect(fill = "white"))+
  ggtitle("Species 1")
ggsave("fitnessTime.png")

# old fitness plot

when = 4
plot_dat_time <- plot_dat[which(plot_dat$Time == unique(plot_dat$Time)[when]),]
plot_dat_time2 <- plot_dat2[which(plot_dat2$Time == unique(plot_dat$Time)[when]),]

ggplot(plot_dat_time) +
  #geom_point(aes(x=as.numeric(variable),y=value,group = w_mat, color = w_mat)) +
  geom_point(aes(x=w_mat,y=Fitness)) +#,group = Time, color = Time), size = 1) +
  geom_point(data = plot_dat_time2, aes(x=w_mat,y=Fitness), color = "red") +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(name = "Maturation size") +
  #scale_colour_gradient2(low = "seagreen4", mid = "turquoise2", high = "dodgerblue4", midpoint = myMat$w_mat[1]) +
  scale_colour_gradient2(low = "blue", mid = "orange", high = "black", midpoint = fishStart) +
  #geom_vline(xintercept = 3000, linetype = "dashed") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),
        legend.justification=c(1,1),
        legend.position = "right",
        legend.key = element_rect(fill = "white"))+
  ggtitle("Species 9")
ggsave("fitness5500yr.png")
