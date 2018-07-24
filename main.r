# This is the main script where you run the simulation (or experiment new things)

#setting things up -----------------------
rm(list = ls())
require(scales)
require(ggplot2)#because always need these two
require(reshape2)
require(plyr)# for aaply
require(grid)# for grid.newpage (plotSummary)
require(abind) # to use abind (bind of arrays)
require(rmarkdown)
require(RColorBrewer)
require(tictoc)
require(limSolve)
#require(tidyverse)

source("MizerParams-class.r") #to get the Constructor
source("selectivity_funcs.r") #to get the knife_edge function
source("methods.r") #I'm doing my own methods then!
source("summaryFunction.r") #to have all the GetSomething functions
source("plotFunction.r") #to draw the plots
source("TBM1.r") # the model from mizer (more like a set up)
source("model.r") # my model 
source("utility.r") # helpful functions

# # little script to check sim content ----------------
# a<- get(load("eta5/fisheries/run1/run.Rdata"))
# a@params@species_params$knife_edge_size
# a@params@interaction
# a@params@species_params$r_max

# multi species simulations ----------
file_name = "/Sim9"
noInter = F
# PARAMETERS
# physio
no_sp = 9
min_w_inf <- 10
max_w_inf <- 1e5
RMAX = T
w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp) # for fisheries gear

# varying param
# parameters worth checking: h, ks, z0pre, sigma, beta, f0, erepro, w_pp_cutoff
# defaults
h = 85
ks = 4
z0pre = 2
sigma = 1
beta = 100
f0 = 0.5
erepro = 1
w_pp_cutoff = 1
interaction = 0.5
overlap = 0.5
eta = 0.25
mAmplitude = 0.2
mu=1
kappa = 0.05

if(noInter)
{
  w_pp_cutoff = 1e5
  interaction = 0
  overlap = 0
  kappa = 0.5
}

# fisheries
gear_names <- rep("FishingStuff", no_sp)
knife_edges <- w_inf * eta 


# other
t_max = 50
no_run = 60
no_sim = 10
i_start = 1
# or
simulationVec <- c(10)
# initialisation phase (4000 yr)
#for (i in i_start:no_sim)
for (i in simulationVec)
{

  
  # switch(i,
  #        "1" = {mu = 0.01},
  #        "2" = {mu = 0.1},
  #        "3" = {mu = 0.5},
  #        "4" = {mu = 1},
  #        "5" = {mu = 1.5},
  #        "6" = {mu = 3},
  #        "7" = {mu = 5},
  #        {})
  

  
  
  tic()
  cat(sprintf("Simulation number %g\n",i))
  path_to_save = paste(getwd(),file_name,"/init/run", i, sep = "")
  
  sim <- myModel(no_sp = no_sp, eta = eta, t_max = t_max, no_run = no_run, min_w_inf = min_w_inf,extinct = T, 
                 max_w_inf = max_w_inf, RMAX = RMAX, 
                 ken = F, initTime = 1, initPool = 9, ks = ks, z0pre = z0pre, f0 = f0, overlap = overlap, sigma = sigma, beta = beta, w_pp_cutoff = w_pp_cutoff,
                 kappa = kappa,
                 OptMutant = "M5", mAmplitude = mAmplitude, mu= mu,
  effort = 0, #knife_edge_size = knife_edges, gear_names = gear_names, 
  save_it = T, path_to_save = path_to_save,
  print_it = T, normalFeeding = F, Traits = "eta")
  
  #rm(sim)
  for (j in 1:20) gc()
  toc()
}
# simulation after initialisation

folder <- paste(getwd(),file_name,sep="")
initFolder <- paste(folder,"/init",sep="")
dirContent <- dir(initFolder)[1:11]
#dirContent <- "run4"
no_run = 60
i_start = 2
# NO fisheries
for (i in i_start:length(dirContent))
{
  # switch(i,
  #        "1" = {mu = 0.01},
  #        "2" = {mu = 0.1},
  #        "3" = {mu = 0.5},
  #        "4" = {mu = 1},
  #        "5" = {mu = 1.5},
  #        "6" = {mu = 3},
  #        "7" = {mu = 5},
  #        {})
  if (file.exists(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = "")))
    path_to_save <- paste(folder,"/normal/",dirContent[i],sep = "")
    
    cat(sprintf("Using %s\n",i))
    
    output <- myModel(no_sp = no_sp, eta = eta, t_max = t_max, no_run = no_run, min_w_inf = min_w_inf,extinct = T, 
                   max_w_inf = max_w_inf, RMAX = RMAX, 
                   ken = F, initTime = 1, initPool = 9, ks = ks, z0pre = z0pre, f0 = f0, overlap = overlap, sigma = sigma, beta = beta, w_pp_cutoff = w_pp_cutoff,
                   kappa = kappa,
                   OptMutant = "M5", mAmplitude = mAmplitude, mu = mu, initCondition = sim,
                   effort = 0, #knife_edge_size = knife_edges, gear_names = gear_names, 
                   save_it = T, path_to_save = path_to_save,
                   print_it = T, normalFeeding = F, Traits = "eta")
    rm(output)
    for (j in 1:20) gc()
  }
}
# Fisheries
dirContent <- "run6"
for (i in i_start:length(dirContent))
{
  # switch(i,
  #        "1" = {mu = 0.01},
  #        "2" = {mu = 0.1},
  #        "3" = {mu = 0.5},
  #        "4" = {mu = 1},
  #        "5" = {mu = 1.5},
  #        "6" = {mu = 3},
  #        "7" = {mu = 5},
  #        {})
  
  if (file.exists(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = "")))
    path_to_save <- paste(folder,"/fisheries/",dirContent[i],sep = "")
    
    cat(sprintf("Using %s\n",i))
    
    output <- myModel(no_sp = no_sp, eta = eta, t_max = t_max, no_run = no_run, min_w_inf = min_w_inf,extinct = T, 
                   max_w_inf = max_w_inf, RMAX = RMAX, 
                   ken = F, initTime = 1, initPool = 9, ks = ks, z0pre = z0pre, f0 = f0, overlap = overlap, sigma = sigma, beta = beta, w_pp_cutoff = w_pp_cutoff,
                   kappa = kappa,
                   OptMutant = "M5", mAmplitude = mAmplitude, mu= mu, initCondition = sim,
                   effort = 0.8, knife_edge_size = knife_edges, gear_names = gear_names, 
                   save_it = T, path_to_save = path_to_save,
                   print_it = T, normalFeeding = F, Traits = "eta")
    rm(output)
    for (j in 1:20) gc()
  }
}

for (i in i_start:length(dirContent))
{
  for (effort in c(seq(0.1,0.7,0.1),0.9,1))
{

  
  if (file.exists(paste(initFolder,"/",dirContent[1],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[1],"/run.Rdata",sep = "")))
    path_to_save <- paste(folder,"/fisheries/effort",effort,"/",dirContent[i],sep = "")

    cat(sprintf("Using %s\n",i))
    
    output <- myModel(no_sp = no_sp, eta = eta, t_max = t_max, no_run = no_run, min_w_inf = min_w_inf,extinct = T, 
                      max_w_inf = max_w_inf, RMAX = RMAX, 
                      ken = F, initTime = 1, initPool = 9, ks = ks, z0pre = z0pre, f0 = f0, overlap = overlap, sigma = sigma, beta = beta, w_pp_cutoff = w_pp_cutoff, kappa = kappa,
                      OptMutant = "M5", mAmplitude = mAmplitude, mu= mu, initCondition = sim,
                      effort = effort, knife_edge_size = knife_edges, gear_names = gear_names, 
                      save_it = T, path_to_save = path_to_save,
                      print_it = T, normalFeeding = F, Traits = "eta")
    rm(output)
    for (j in 1:20) gc()
  }
}
}


#with parallel / need to update the function--------------------------- 
rm(list = ls())
library(parallel)
library(ggplot2)#because always need these two
library(reshape2)
library(plyr)# for aaply
library(grid)# for grid.newpage (plotSummary)
library(abind) # to use abind (bind of arrays)
library(rmarkdown)
library(RColorBrewer)
library(tictoc)

source("MizerParams-class.r") #to get the Constructor
source("selectivity_funcs.r") #to get the knife_edge function
source("methods.r") #I'm doing my own methods then!
source("summaryFunction.r") #to have all the GetSomething functions
source("plotFunction.r") #to draw the plots
source("TBM1.r") # the model from mizer (more like a set up)
source("model.r") # my model 
source("utility.r") 

#(optional) record start time, for timing
ptm=proc.time()
tic()

#unsure  what this setting does
options(warn=-1) #?

#Adjust this for num of targeted cpu/cores
# e.g. Numcores = detectCores()-1

where = paste(getwd(),"/parallel",sep="")

numcores=4
cl <- makeForkCluster(getOption("cl.cores", numcores), outfile = "")
sim <- clusterApplyLB(cl
                      ,x=1:numcores
                      ,fun=multiRun
                      ,no_sp = 9
                      ,t_max = 50
                      ,mu = 5
                      ,no_run = 80
                      ,min_w_inf = 10
                      ,max_w_inf = 10e5
                      ,effort = 0
)
stopCluster(cl)

## Option 1: future package (and safely)
library(future)
plan(multiprocess)
## optionally, safely
safe_multiRun <- purrr::safely(multiRun)
sim <- future::future_lapply(1:numcores, safe_multiRun, no_sp , )

library(purrr)
## Option 2: purrr package
safe_multiRun <- purrr::safely(multiRun)
sim <- purrr::map(1:numcores, safe_multiRun, no_sp = 9, ...)

#(optional) compare end with start time, for timing


# saving
for (i in 1:length(sim)) 
{
  path_to_save = paste(where,"/run",i,sep="")
  ifelse(!dir.exists(file.path(path_to_save)), dir.create(file.path(path_to_save),recursive = T), FALSE)
  saveRDS(file = paste(path_to_save,"/run.RDS",sep=""),object = sim[[i]])
}

print((proc.time()-ptm)/60.0)
toc()






# working on that right now -------------------


rownames(interactionBeta) <- c("1","2","3","4")
colnames(interactionBeta) <- c("1","2","3","4")
rownames(interactionAlpha) <- c("1","2","3","4","5")
colnames(interactionAlpha) <- c("1","2","3","4","5")

interactionAlpha<-interactionAlpha[-3,-3]


which(rownames(interactionAlpha) != rownames(interactionBeta))

a <- rownames(interactionAlpha)
b <- rownames(interactionBeta)

 c <- which(!(a %in% b))

interactionSave <- rbind(interactionBeta,interactionAlpha[c,])
interactionSave <- cbind(interactionSave,interactionAlpha[,c])

# investigate this fucking growth
object <- get(load("ParamChap1/init/run4/run.Rdata"))
# Plot realised intake versus maximum intake of small and large individuals to see what is causing decrease in growth at large sizes. 
# Is it different among large and small species? 
# Is it food limitation or is metabolism too high? 
# If this is food limitation that affects growth, would changing PPMR or feeding kernel improve growth? 
# And how much do you need to change it for some substantial effect to happen?

#look at mortality
plotScythe(object)



# trait value picking
# Trait = eta
mAmplitude = 0.05
eta = 0.5
sd = as.numeric(mAmplitude *  eta) # standard deviation
#x <- eta + rnorm(1, 0, sd) # change a bit eta

df = NULL
for (i in 1:10000)  df <- c(df,( rnorm(1, 0, sd)))

summary(df)
plot(df)
plot(density(df))

x <- seq(-0.5,0.5, length=500)
y<-dnorm(x,mean=0, sd=0.025)
plot(x,y, type="l")

# Plots of every kind of output + for loop to see the variation of one parameter ------------------

res = 1000 # figure resolution
subdir = "/weighted" # where to store the plots
parameter = "none" # name of parameter varying for plot title
t_max = 100
no_sp = 4
mu = 5
for (i in c(1 %o% 10^(-5:0)))
{
  output <- myModel(no_sp = no_sp, t_max = t_max, mu = mu, OptMutant = "yo", no_run = 1,
                    min_w_inf = 10, ks=2, max_w_inf = 10000, 
                    #param = sim@params, # option to give param of another sim to have mutant relations
                    effort = 0, data = TRUE)
  
  # when data = true, mutation do not work but I get the values of lots of function at each time step of the simulation
  # do that for short runs
  # sort the output
  energy = output[[1]]
  rd = output[[2]]
  eggs = output[[3]]
  sim = output[[4]]
  food = output[[5]]
  m2 = output[[6]]
  z = output[[7]]
  m2_background = output[[8]]
  phi_fish = output[[9]]
  phi_pltk = output[[10]]
  end = dim(energy)[1]
  
  # thing to fix: if I give parameters to the sim, it won't have the right name (sp name instead of ecotype)
  # if there are no mutants I guess its fine
  dimnames(sim@n)$sp = sim@params@species_params$ecotype
  
  ifelse(!dir.exists(file.path(dir, subdir)), dir.create(file.path(dir, subdir)), FALSE) #create the file if it does not exists
  dir.create(file.path(dir,subdir,"/reproduction")) # create tree files to ease comparison
  dir.create(file.path(dir,subdir,"/growth"))
  dir.create(file.path(dir,subdir,"/mortality"))
  dir.create(file.path(dir,subdir,"/spawn"))
  dir.create(file.path(dir,subdir,"/RDD"))
  dir.create(file.path(dir,subdir,"/feeding"))
  
  # plots ----------------
  plotDynamics(sim)
  setwd(paste(dir,subdir, sep = "")) #to have the figures in the right directory
  mytitle = paste("biomass_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  plotSS(sim)
  setwd(paste(dir,subdir, sep = "")) #to have the figures in the right directory
  mytitle = paste("sizespectrum_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # RDI
  rdi <- rd[,,1]
  RDI <- melt(rdi)
  print(ggplot(RDI) +
          geom_line(aes(x=Time,y=value ,colour = as.factor(Species))) +
          scale_x_continuous(name = "Time") + 
          scale_y_log10(name = "Energy") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Reproduction Density Independent"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("rdi_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # RDD
  rdd <- rd[,,2]
  RDD <- melt(rdd)
  print(ggplot(RDD) +
          geom_line(aes(x=Time,y=value ,colour = as.factor(Species))) +
          scale_x_continuous(name = "Time") + 
          scale_y_log10(name = "Energy") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Reproduction Density Dependent"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("rdd_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # ratio RDD/RDI
  ratio <- rdd/rdi
  RAT <- melt(ratio)
  print(ggplot(RAT) +
          geom_line(aes(x=Time,y=value ,colour = as.factor(Species))) +
          scale_x_continuous(name = "Time") + 
          scale_y_log10(name = "Ratio", breaks = c(1 %o% 10^(-5:-2)) ) +
          scale_colour_discrete(name = "Species") +
          ggtitle("RDD/RDI"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("RddRdi_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 

  # e
  # energy after metabolism, for the moment equal between every species
  e <- energy[,,,1]
  etot = apply(e, c(1,2), sum)
  E <- melt(etot)
  ggplot(E) +
    geom_line(aes(x=Time,y=value, colour = as.factor(Species))) +  # the as.factor convert to discrete as linetype doesnt work with continuous value
    scale_x_continuous(name = "Time") + 
    scale_y_continuous(name = "Energy")+
    scale_colour_discrete(name = "Species") +
    ggtitle("Total energy available after metabolism")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("energy_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # energy by weight by sp at simulation end
  eSP = e[end,,]
  ESP <- melt(eSP)
  ggplot(ESP) +
    geom_line(aes(x=Size,y=value,colour = as.factor(Species))) +
    scale_x_log10(name = "Weight") + 
    scale_y_continuous(name = "Energy")+
    scale_colour_discrete(name = "Species") +
    ggtitle("Energy available after metabolism by weight")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("energy_size_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  #growth
  # energy for through time
  g <- energy[,,,3]
  gtot = apply(g, c(1,2), sum)
  G <- melt(gtot)
  ggplot(G) +
    geom_line(aes(x=Time,y=value, colour = as.factor(Species)))+
    scale_x_continuous(name = "Time") + 
    scale_y_continuous(name = "Energy") +
    scale_colour_discrete(name = "Species") +
    ggtitle("Energy available for growth")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("growth_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # plot of energy by weight by sp at simulation end
  gSP = g[end,,]
  GSP <- melt(gSP)
  ggplot(GSP) +
    geom_line(aes(x=Size,y=value,colour = as.factor(Species)))+
    scale_x_log10(name = "Weight") + 
    scale_y_continuous(name = "Energy") +
    scale_colour_discrete(name = "Species") +
    ggtitle("Energy available for growth by weight")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("growth_size_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # reproduction
  # energy through time
  s <- energy[,,,2]
  
  stot = apply(s, c(1,2), sum)
  S <- melt(stot)
  print(ggplot(S) +
          geom_line(aes(x=Time,y=value, colour = as.factor(Species)))+
          scale_x_continuous(name = "Time") +
          scale_y_continuous(name = "Energy") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Energy available for reproduction"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("reproduction_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # energy by weight by sp at simulation end
  sSP = s[end,,]
  stot = apply(s, c(1,2), sum)
  SSP <- melt(sSP)
  print(ggplot(SSP) +
          geom_line(aes(x=Size,y=value,colour = as.factor(Species))) +
          scale_x_log10(name = "Weight") +
          scale_y_continuous(name = "Energy") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Energy available for reproduction by weight"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("reproduction_size_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # energy by weight by sp at simulation end and weighted by n 
  sSP = s[end,,]
  sSPN = sSP * sim@n[dim(sim@n)[1],,]
  SSPN <- melt(sSPN)
  print(ggplot(SSPN) +
          geom_line(aes(x=Size,y=value,colour = as.factor(Species))) +
          scale_x_log10(name = "Weight") +
          scale_y_log10(name = "Eggs in g/m3") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Real reproduction"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("weighted_reproduction_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # plot of number of eggs by sp by size at end sim
  EGG = melt(eggs)
  print(ggplot(EGG) +
          geom_line(aes(x=Time,y=value,colour = as.factor(Species))) +
          scale_x_continuous(name = "TIme") + 
          scale_y_log10(name = "Eggs in g/m3") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Boudarie condition"))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("spawn_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # feeding
  # throught time
  feeding <- energy [,,,4]
  ftot = apply(feeding, c(1,2), sum)
  FEED <- melt(ftot)
  ggplot(FEED) +
    geom_line(aes(x=Time,y=value, colour = as.factor(Species)))+
    scale_x_continuous(name = "Time") + 
    scale_y_continuous(name = "Energy") +
    scale_colour_discrete(name = "Species") +
    ggtitle("Energy issue from feeding")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("feeding_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # energy by weight by sp at simulation end
  fSP = feeding[end,,]
  FSP <- melt(fSP)
  ggplot(FSP) +
    geom_line(aes(x=Size,y=value,colour = as.factor(Species)))+
    scale_x_log10(name = "Weight") + 
    scale_y_continuous(name = "Energy") +
    scale_colour_discrete(name = "Species") +
    ggtitle("Energy issue from feeding by weight")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("feeding_size_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # Phi
  a = phi_fish[end,1,]
  A = melt(a)
  b= phi_pltk[end,1,]
  B = melt(b)
  feeding = energy[end,1,,4] # feeding level of one sp as they have the same profile   
  Fe = melt(feeding)
  S = melt(sim@params@search_vol[1,]) # search volume
  
  # plot of phi and others
  ggplot()+
    geom_line(data = A,aes(x = as.numeric(rownames(A)), y = value, color = "Phi fish")) +
    geom_line(data = B, aes(x = as.numeric(rownames(B)), y = value, color = "Phi plankton")) +
    # geom_line(data = Fe, aes(x = as.numeric(rownames(Fe)), y = value, color = "Feeding level")) +
    #geom_line(data = S, aes(x = as.numeric(rownames(S)), y = value, color = "Search Volume")) +
    scale_x_log10(name = "Predator size",breaks = c(1 %o% 10^(-10:5)))+
    scale_y_continuous(name = "value of phy prey")+
    ggtitle("Relative proportion of food eaten between plankton and fish")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("phi_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # Mortality
  # predation mortality
  a = m2[end,,]
  A = melt(a)
  ggplot(A) + 
    geom_line(aes(x = PreySize, y = value, color = as.factor(PreySp)))+
    scale_x_log10()+
    scale_y_continuous( limits = c(0,30))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("PredMort_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # total mortality
  a = z[end,,]
  A = melt(a)
  ggplot(A) + 
    geom_line(aes(x = PreySize, y = value, color = as.factor(PreySp)))+
    scale_x_log10()+
    scale_y_continuous( limits = c(0,30))
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("TotMort_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  #mortality on plankton
  a = m2_background[end,]
  A = melt(a)
  A = cbind(A,rownames(A))
  colnames(A) = c("value", "size")
  ggplot(A) + 
    geom_line(aes(x = as.numeric(size), y = value, group = 1))+
    scale_y_log10() +
  scale_x_log10()
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste("PlktMort_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  #weighted plots -------------------
  
  for (j in seq(t_max,t_max*10,t_max))
  {
  time = j
  
  if (time == 1000) time = 992 # I know my sim is weird (last step is 992)
  # reproduction by weight by sp at simulation end and weighted by n 
  
  s <- energy[,,,2]
  sSP = s[time,,]
  sSPN = sSP * sim@n[time,,]
  SSPN <- melt(sSPN)
  name = paste("Real reproduction at time ",time, sep ="")
  print(ggplot(SSPN) +
          geom_line(aes(x=Size,y=value,colour = as.factor(Species))) +
          scale_x_log10(name = "Size") +
          scale_y_log10(name = "Eggs in g/m3") +
          scale_colour_discrete(name = "Species") +
          ggtitle(name))
  
  setwd(paste(dir,subdir,"/reproduction", sep = ""))
  mytitle = paste("weighted_reproduction_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  
  # growth by weight by sp at simulation end and weighted by n 
  g <- energy[,,,3]
  gSP = g[time,,]
  gSPN = gSP * sim@n[time,,]
  GSPN <- melt(gSPN)
  name = paste("Real growth at time ",time, sep ="")
  ggplot(GSPN) +
    geom_line(aes(x=Size,y=value,colour = as.factor(Species)))+
    scale_x_log10(name = "Size") + 
    scale_y_log10(name = "Energy") +
    scale_colour_discrete(name = "Species") +
    ggtitle(name)
  
  setwd(paste(dir,subdir,"/growth", sep = ""))
  mytitle = paste("growth_size_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # energy by weight by sp at simulation end
  feeding <- energy[,,,4]
  fSP = feeding[time,,]
  fSPN = fSP * sim@n[time,,]
  FSPN <- melt(fSPN)
  name = paste("Energy issue from feeding weighted by abundance of species at time ",time, sep ="")
  ggplot(FSPN) +
    geom_line(aes(x=Size,y=value,colour = as.factor(Species)))+
    scale_x_log10(name = "Size") + 
    scale_y_log10(name = "Feeding level") +
    scale_colour_discrete(name = "Species") +
    ggtitle(name)
  
  setwd(paste(dir,subdir,"/feeding", sep = ""))
  mytitle = paste("feeding_size_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # # predation rate (to set up)
  # pred <- food[,,,4]
  # fSP = feeding[time,,]
  # fSPN = fSP * sim@n[time,,]
  # FSPN <- melt(fSPN)
  # name = paste("Energy issue from feeding weighted by abundance of species at time ",time, sep ="")
  # ggplot(FSPN) +
  #   geom_line(aes(x=Size,y=value,colour = as.factor(Species)))+
  #   scale_x_log10(name = "Size") + 
  #   scale_y_log10(name = "Feeding level") +
  #   scale_colour_discrete(name = "Species") +
  #   ggtitle(name)
  # 
  # setwd(paste(dir,subdir, sep = ""))
  # mytitle = paste("feeding_size_", parameter, "_",i,".png", sep = "")
  # dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # total mortality
  mortality = z[time,,]
  mN = mortality * sim@n[time,,]
  MN = melt(mN)
  name = paste("Total mortality at time ",time, sep ="")
  ggplot(MN) + 
    geom_line(aes(x = PreySize, y = value, color = as.factor(PreySp)))+
    scale_x_log10()+
    scale_y_log10()+
    scale_colour_discrete(name = "Species") +
    ggtitle(name)
  
  setwd(paste(dir,subdir,"/mortality", sep = ""))
  mytitle = paste("TotMort_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  
  # egg number
  EGG = melt(eggs)
  print(ggplot(EGG) +
          geom_line(aes(x=Time,y=value,colour = as.factor(Species))) +
          scale_x_continuous(name = "TIme") + 
          scale_y_log10(name = "Eggs in g/m3") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Boudarie condition"))
  
  setwd(paste(dir,subdir,"/spawn", sep = ""))
  mytitle = paste("spawn_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  # RDD
  rdd <- rd[,,2]
  RDD <- melt(rdd)
  print(ggplot(RDD) +
          geom_line(aes(x=Time,y=value ,colour = as.factor(Species))) +
          scale_x_continuous(name = "Time") + 
          scale_y_log10(name = "Energy") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Reproduction Density Dependent"))
  
  setwd(paste(dir,subdir,"RDD", sep = ""))
  mytitle = paste("rdd_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  
  # RDI
  rdi <- rd[,,1]
  RDI <- melt(rdi)
  RDI <- RDI[RDI$value >= min_value,]
  print(ggplot(RDI) +
          geom_line(aes(x=Time,y=value ,colour = as.factor(Species))) +
          scale_x_continuous(name = "Time") + 
          scale_y_log10(name = "Energy") +
          scale_colour_discrete(name = "Species") +
          ggtitle("Reproduction Density Independent"))
  
  setwd(paste(dir,subdir,"/RDI", sep = ""))
  mytitle = paste("rdi_",time,"_", parameter, "_",i,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  
  }
}

# predation traits analyses / detail of the predation equation here (not updated though)--------------

# phi prey
# n_eff_prey is the total prey abundance by size exposed to each predator
# (prey not broken into species - here we are just working out how much a predator eats - not which species are being eaten - that is in the mortality calculation
n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*") 
# Quick reference to just the fish part of the size spectrum
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
# predKernal is predator x predator size x prey size
# So multiply 3rd dimension of predKernal by the prey abundance
# Then sum over 3rd dimension to get total eaten by each predator by predator size
phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*"),dims=2)
# Eating the background
phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*"),dims=2)
return(phi_prey_species+phi_prey_background)

#feeding level
encount <- object@search_vol * phi_prey
# calculate feeding level
f <- encount/(encount + object@intake_max)
return(f)

#pred rate
n_total_in_size_bins <- sweep(n, 2, object@dw, '*')
pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*")
return(pred_rate)


#pred kernel
res@pred_kernel[] <- object$beta
res@pred_kernel <- exp(-0.5*sweep(log(sweep(sweep(res@pred_kernel,3,res@w_full,"*")^-1,2,res@w,"*")),1,object$sigma,"/")^2)
res@pred_kernel <- sweep(res@pred_kernel,c(2,3),combn(res@w_full,1,function(x,w)x<w,w=res@w),"*") # find out the untrues and then multiply




# trait study --------------
# draw plots that show the growth rate with different trait varying


# need some n values to get the rest
sim <- myModel(no_sp = 9, t_max = 50, mu = 5, OptMutant = "yo", RMAX = TRUE, hartvig = TRUE)
endList <- length(sim) # shortcut to have ref to the last simulation which has the right dim, names, ...
PSim <- sim[[endList]] # if I want to look at params and such I'm taking the last sim
PSim@params@species_params
plotDynamics(PSim)
end = dim(PSim@n)[1]

# and some parameters
eta = 0.25
z0pre = 0.84
n = 0.75 # exponent of maximum intake (scaling of intake)
q = 0.8 # exponent of search volume
kappa = 0.005 # ressource spectrum carrying capacity
lambda = 2+q-n # exponent of the background spectrum.
h = 85 # factor of maximum intake
f0 = 0.6 # average feeding level of the community/feeding level of small individuals feeding on background

# Asymptotic size

min_w_inf = 10
max_w_inf = 10e5
w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=1000) # asymptotic mass of the species
w_mat <- w_inf * eta
z0 <- z0pre * w_inf^(n-1)

size = data.frame(w_inf,w_mat,z0)
ggplot(size) +
  geom_line(aes(x = w_inf, y = w_mat, color = "Maturation size")) +
  geom_line(aes(x = w_inf, y = z0, color = "Background mortality")) +
  scale_x_log10(name = "Asymptotic size") +
  scale_y_log10(name = "Size") +
  ggtitle("Effect of varition of asymptotic size")

# w_mat is only used in psi (allocation reproduction)
# w_inf is used for h and I dont know what that is

# PPMR 
beta = 100 # preferred predator-prey weight ratio
sigma = 1.3 # width of selection function

beta = seq (10,200,10)
alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2)
gamma <- h * f0 / (alpha_e * kappa * (1-f0))

PPMR <- data.frame(beta,gamma)
ggplot(PPMR)+
  geom_line(aes(x = beta, y = gamma))+
  ggtitle("Gamma function of beta")

# impact of beta variation
beta_min = 10
beta_max = 200
dBeta = 10
results = list()
for (i in seq (beta_min,beta_max,dBeta))
{
  sim <- myModel(no_sp = 9, t_max = 50, OptMutant = "yo", RMAX = TRUE, min_w_inf = 10, max_w_inf = 10000, beta = i, extinct = FALSE, hartvig = TRUE)
  sim <- sim[[endList]] 
  
a = getPhiPrey(object = sim@params, n=sim@n[end,,], n_pp = sim@n_pp[end,])
b = getFeedingLevel(object = sim@params, n=sim@n[end,,], n_pp = sim@n_pp[end,], phi_prey = a)
betaPred = cbind(a[2,],b[2,])
name <- paste('beta',i,sep='')
results[[name]] = betaPred
}

#plots
# beta by size
res = 600
for (i in seq (beta_min,beta_max,dBeta)) 
{
  name <- paste('beta',i,sep='')
pred = as.data.frame(results[[name]])
print(
  ggplot(pred) +
geom_line(aes(x = as.numeric(rownames(pred)), y=V2, colour = "Feeding level"), group = 1)+
  geom_line(aes(x = as.numeric(rownames(pred)), y=V1,colour = "Phi prey"), group = 1)+
  scale_x_log10(name = "Size")+
  scale_y_continuous(name = "Function output",limits = c(0,0.7))+
  ggtitle(name)
)
setwd(paste(dir,"/Traits/Beta", sep = ""))
mytitle = paste(name,".png", sep = "")
dev.print(png, mytitle, width = res, height = res) 
}

# global impact of beta on feeding and phi when summing weights
bigBeta = matrix(data = NA, nrow = length(seq (beta_min,beta_max,dBeta)), ncol = 2, dimnames = list(c(seq (beta_min,beta_max,dBeta)), c("Phi","Feed")))
for (i in seq (beta_min,beta_max,dBeta)) 
{
  name <- paste('beta',i,sep='')
  bigBeta[i/dBeta,] = colSums(results[[name]]) 
}

bigBeta = as.data.frame(bigBeta)

  ggplot(bigBeta) +
    geom_line(aes(x = as.numeric(rownames(bigBeta)), y=Feed, colour = "Feeding level"), group = 1)+
    geom_line(aes(x = as.numeric(rownames(bigBeta)), y=Phi,colour = "Phi prey"), group = 1)+
    scale_x_continuous(name = "Beta value")+
    scale_y_continuous(name = "Function output")+
    ggtitle("Impact of beta")

setwd(paste(dir,"/Traits/Beta", sep = ""))
mytitle = paste("betaVar",".png", sep = "")
dev.print(png, mytitle, width = res, height = res) 


# sigma
sigma_min = 0.1
sigma_max = 2
dSigma = 0.1
results = list()
for (i in seq (sigma_min,sigma_max,dSigma))
{
  sim <- myModel(no_sp = 9, t_max = 50, OptMutant = "yo", RMAX = TRUE, min_w_inf = 10, max_w_inf = 10000, sigma = i, extinct = FALSE, hartvig = TRUE)
  
  sim <- sim[[endList]] # if I want to look at params and such I'm taking the last sim
  
  a = getPhiPrey(object = sim@params, n=sim@n[end,,], n_pp = sim@n_pp[end,])
  b = getFeedingLevel(object = sim@params, n=sim@n[end,,], n_pp = sim@n_pp[end,], phi_prey = a)
  sigmaPred = cbind(a[2,],b[2,])
  name <- paste('sigma',i,sep='')
  results[[name]] = sigmaPred
}

#plots
# sigma by size
res = 600
for (i in seq (sigma_min,sigma_max,dSigma)) 
{
  name <- paste('sigma',i,sep='')
  pred = as.data.frame(results[[name]])
  print(
    ggplot(pred) +
      geom_line(aes(x = as.numeric(rownames(pred)), y=V2, colour = "Feeding level"), group = 1)+
      geom_line(aes(x = as.numeric(rownames(pred)), y=V1,colour = "Phi prey"), group = 1)+
      scale_x_log10(name = "Size")+
      scale_y_continuous(name = "Function output",limits = c(0,0.8))+
      ggtitle(name)
  )
  setwd(paste(dir,"/Traits/Sigma", sep = "")) #to have the figures in the right directory
  mytitle = paste(name,".png", sep = "")
  dev.print(png, mytitle, width = res, height = res) 
}

# energy by sigma
bigSigma = matrix(data = NA, nrow = length(seq (sigma_min,sigma_max,dSigma)), ncol = 2, dimnames = list(c(seq (sigma_min,sigma_max,dSigma)), c("Phi","Feed")))
for (i in seq (sigma_min,sigma_max,dSigma)) 
{
  name <- paste('sigma',i,sep='')
  idx = i/dSigma
  bigSigma[idx,] = colSums(results[[name]]) 
}

bigSigma = as.data.frame(bigSigma)

ggplot(bigSigma) +
  geom_line(aes(x = as.numeric(rownames(bigSigma)), y=Feed, colour = "Feeding level"), group = 1)+
  geom_line(aes(x = as.numeric(rownames(bigSigma)), y=Phi,colour = "Phi prey"), group = 1)+
  scale_x_continuous(name = "Sigma value")+
  scale_y_continuous(name = "Function output")+
  ggtitle("Impact of sigma")

setwd(paste(dir,"/Traits/Sigma", sep = "")) #to have the figures in the right directory
mytitle = paste("sigmaVar",".png", sep = "")
dev.print(png, mytitle, width = res, height = res) 
    


# Other parameters variation I dont remember what that is------------------------------
#psi

psi = PSim@params@psi
#psi = as.data.frame(psi)
PSI = melt(psi)
ggplot(data = PSI, aes(x = w, y = value), group = sp) +
geom_point()  +
scale_x_continuous(breaks = c(1 %o% 10^(-3:5)))
  
res@psi[] <- unlist(tapply(res@w,1:length(res@w),function(wx,w_inf,w_mat,n)
  {
  ((1 + (wx/(w_mat))^-10)^-1) * (wx/w_inf)^(1-n)
  }
  ,w_inf=object$w_inf,w_mat=object$w_mat,n=n))

# metabolsim maintenance

es@std_metab[] <-  unlist(tapply(res@w,1:length(res@w),function(wx,ks,p)
  
  ks * wx^p
  
  , ks=object$ks,p=p))

ks = 4
p = 0.75
size = as.numeric(dimnames(PSim@n)$w)
metabolism = ks*size^p
ratio = metabolism/size
truc = cbind(size,metabolism,ratio)
truc = as.data.frame(truc)
ggplot(truc)+
  geom_line(aes(x=size,y=metabolism)) +
  geom_line(aes(x=size, y=ratio))+
scale_x_log10() +
  scale_y_log10()



# plot biomass sum by family, to see if the relative biomass difference between the species change when I introduce new ecotypes 
# I could add stars when a new ecotype appear on the graph to do that need to do a geom_point (data = , aes ...)

truc = getBiomass(sim)
dimnames(truc)$sp <- sim@params@species_params$species
truc <- as.data.frame(truc)
Struc <- sapply(unique(names(truc)[duplicated(names(truc))]), 
         function(x) Reduce("+", truc[ , grep(x, names(truc))]) ) # magic thing that sum col with same names
names(dimnames(Struc)) <- list("Time","Species")
TRUC = melt(Struc)
ggplot(TRUC)+
  geom_line(aes(x = Time, y = value, colour = as.factor(Species)))


# egg interference
I = exp(-(log(mi/mj)^2)/2*sigma^2)

f= I *sum(vol search rate * n * dw ) of w

n_total_in_size_bins <- sweep(n, 2, object@dw, '*')
object@search_vol*n_total_in_size_bins


# plot function of egg reduction
no_sp = 10
sim <- myModel(no_sp = no_sp, t_max = 20, mu = 0, OptMutant = "yo", RMAX = TRUE, cannibalism = 1, r_mult = 1e0, erepro = 0.001, p =0.75, ks=4, extinct = FALSE, k0 = 25)
endList <- length(sim) # shortcut to have ref to the last simulation which has the right dim, names, ...
PSim <- sim[[endList]] # if I want to look at params and such I'm taking the last sim
r_max = PSim@params@species_params$r_max
# plotFeedingLevel(PSim)
# plotM2(PSim)
# plotDynamics(PSim)
# plotSS(PSim) 
rdi = seq(0,1,0.001) #fake rdi to get values
# mizer egg production

a = matrix(nrow = length(rdi), ncol = no_sp, dimnames = list(as.character(rdi),as.character(c(1:no_sp))))
names(dimnames(a)) = list("RDI","Species")

for(i in 1:dim(a)[1])
{
  for (j in 1:dim(a)[2])
  {
    a[i,j] = r_max[j] * rdi[i] / (r_max[j]+rdi[i]) 
  }
}
# a is the matrix showing the recruitment (RDI processed by rmax) in function of the rdi
MEgg = melt(a)
ggplot(MEgg) +
  geom_line(aes(x = RDI, y = value, colour = as.factor(Species))) +
  scale_y_continuous(name = "Recruitement", limits = c(0,0.08)) +
  geom_abline(intercept = 0, slope = 1)

b = sweep(a,2,r_max,"/") # a divided by rmax for the graph

MEggR = melt(b)
ggplot(MEggR) +
  geom_line(aes(x = RDI, y = value, colour = as.factor(Species))) +
  scale_y_continuous(name = "Recruitement/Rmax", limits = c(0,2.5)) +
  scale_x_continuous(name = "RDI") +
  geom_abline(intercept = 0, slope = 1)

# what changes when I poke parameters ? nothing




#test starvation mortality
  
    
    




   # feeding plots----------------
   #behvior of gamma with the traits
   # gamma study
   n = 0.75 # exponent of maximum intake (scaling of intake)
   p = 0.75 # exponent of standard metabolism
   q = 0.8 # exponent of search volume
   lambda = 2+q-n # exponent of the background spectrum.
   h = 85 # factor of maximum intake
   beta = 100 # preferred predator-prey weight ratio
   sigma = 1.3 # width of selection function
   f0 = 0.6 # average feeding level of the community/feeding level of small individuals feeding on background
   kappa = 0.008 # ressource spectrum carrying capacity
   
   #plots
   # beta
   beta = seq(50,150,1)
   sigma = 1.3
   gamma <- h * f0 / ((sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2)) * kappa * (1-f0))
   betaM = matrix(data = cbind(beta,gamma), nrow = length(beta), ncol = 2, dimnames = list(NULL,c("beta","gamma")))
   betaDF = as.data.frame(betaM)
   ggplot(betaDF)+
     geom_line(aes(x = beta, y =gamma))+
     scale_x_continuous(name = "beta (PPMR)")+
     scale_y_continuous(name = "gamma (factor for search volume)")
   
   setwd(paste(dir,subdir, sep = ""))
   mytitle = "beta_gamma.png"
   dev.print(png, mytitle, width = res, height = 0.6*res) 
   
   #sigma
   sigma = seq(0.1,2.5,0.025)
   beta = 100
   gamma <- h * f0 / ((sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2)) * kappa * (1-f0))
   sigmaM = matrix(data = cbind(sigma,gamma), nrow = length(sigma), ncol = 2, dimnames = list(NULL,c("beta","gamma")))
   sigmaDF = as.data.frame(sigmaM)
   ggplot(sigmaDF)+
     geom_line(aes(x = sigma, y =gamma))+
     scale_y_log10(breaks = c(1000,5000,10000,50000))+
   scale_x_continuous(name = "sigma (diet breadth)")+
     scale_y_continuous(name = "gamma (factor for search volume)")
   
   setwd(paste(dir,subdir, sep = ""))
   mytitle = "sigma_gamma.png"
   dev.print(png, mytitle, width = res, height = 0.6*res) 
   
   # building a matrix to plot a surface
   beta = seq(10,200,1)
   sigma = seq(0.5,2.5,0.025)
   mat = matrix(data = NA, nrow = length(sigma),ncol = length(beta))
   for(i in 1:dim(mat)[1])
     for (j in 1:dim(mat)[2])
       mat[i,j] =  h * f0 / ((sqrt(2*pi) * sigma[i] * beta[j]^(lambda-2) * exp((lambda-2)^2 * sigma[i]^2 / 2)) * kappa * (1-f0))
  
   
   dimnames(mat) = list(as.character(sigma),as.character(beta))
  data = melt(mat)
  colnames(data) = c("sigma","beta", "gamma")
  data$logG = log10(data$gamma) # check log values
  
  # ggplot(data)+
  #   geom_raster(aes(x = sigma, y = beta, fill = logG))+
  # scale_fill_gradient(low = "white",high = "black")
  
  ggplot(data)+
    geom_raster(aes(x = sigma, y = beta, fill = gamma))+
    scale_fill_gradient(low = "white",high = "black")+
    scale_x_continuous(name = "sigma (diet breadth)")+
    scale_y_continuous(name = "beta (PPMR)")
  
  setwd(paste(dir,subdir, sep = ""))
  mytitle = "sigma_beta.png"
  dev.print(png, mytitle, width = res, height = 0.6*res) 
 
# the result is shit  

  # relationship between traits ------------------
  
  # beta/sigma ratio
  for (i in SpIdx)
  {
    # empty matrix of ecotype of species i by time
    A = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:SumPar$timeMax[1])))
    # fill the matrix with ones when the ecotype exists
    for (x in 1:nrow(A)) # I'm sure I can do an apply but don't know how
    {
      for (j in 1:ncol(A))
      {
        if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) A[x,j] = 1
      }}
    # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
    # change the ones by the trait value of the ecotype
    BetaA = A * TT[TT$Lineage == i,]$PPMR
    SigmaA = A * TT[TT$Lineage == i,]$Diet_breadth
    # calculate mean trait value at each time step
    no_trait = apply(A,2,sum) # this vector is the number of traits present at each time step
    BetaSum = apply(BetaA,2,sum) # this vector is the sum of the traits value at each time step
    SigmaSum = apply(SigmaA,2,sum) # this vector is the sum of the traits value at each time step
    BetaMean = BetaSum/no_trait # this is the mean trait value at each time step
    SigmaMean = SigmaSum/no_trait # this is the mean trait value at each time step
    
    # Matrix with all traits combination and at what time they go extinct
    TTi = TT[TT$Lineage == i,]
    comb=data.frame(TTi$Ecotype,TTi$Apparition,TTi$Extinction,TTi$PPMR,TTi$Diet_breadth)
    
    # plot of extinction of combinations
    # title = paste("Combination of PPMR and diet breath value of species ",i, sep = "")
    # print(
    #   ggplot(comb) +
    #     geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Extinction)) +
    #     scale_x_continuous(name = "PPMR") +
    #     scale_y_continuous(name = "Diet breath") +
    #     scale_color_continuous(name = "Extinction time in year / dt") +
    #     ggtitle(title)
    # )
    # name = paste("Extinction BetaSigma of species",i, sep = "")
    # setwd(paste(dir,subdir, sep = ""))
    # mytitle = paste(name,".png", sep = "")
    # dev.print(png, mytitle, width = res, height = 2/3* res)
    # 
    # #plot of apparition of combinations
    # print(
    #   ggplot(comb) +
    #     geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Apparition)) +
    #     scale_x_continuous(name = "PPMR") +
    #     scale_y_continuous(name = "Diet breath") +
    #     scale_color_continuous(name = "Apparition time in year / dt") +
    #     ggtitle(title)
    # )
    # name = paste("Apparition BetaSigma of species",i, sep = "")
    # setwd(paste(dir,subdir, sep = ""))
    # mytitle = paste(name,".png", sep = "")
    # dev.print(png, mytitle, width = res, height = 2/3* res)
    
    
    # Mean combination values throughout sim
    stat = data.frame(BetaMean,SigmaMean) 
    #get rid of duplicates
    stat = stat[!duplicated(stat),]
    # need to work on the time because of fucked up legend
    stat = cbind(stat,rownames(stat))
    dimnames(stat)[[2]] = list("BetaMean", "SigmaMean", "Time")
    stat$Time <- as.numeric(substr(stat$Time,1,5)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5

    print(
      ggplot(stat) +
        geom_point(data = stat, aes(x=BetaMean,y=SigmaMean, color = Time)) +
        scale_x_continuous(name = "PPMR") +
        scale_y_continuous(name = "Diet breath") +
        scale_color_continuous(name = "Time in year / dt", low = "blue", high = "red")
    )
    name = paste("MeanBS_SP",i, sep = "")
    setwd(paste(dir,subdir, sep = ""))
    mytitle = paste(name,".png", sep = "")
    dev.print(png, mytitle, width = res, height = 2/3* res)
  }
  

  
  
  
   
  
  #rmax------------
  # some rmax plots 
  #Plot rmax per species and rdd per ecotypes
  #rdd is when rmax is applied
  rdd = rd[,,2] #rdd per ecotype per size
  rdi = rd[,,1]
  rmax = sim@params@species_params$r_max # rmax per species
  
  RDD = apply(rdd,c(1,2),sum)# sum through sizes
  RDI = apply(rdi,c(1,2),sum)
  
  # I need to decide for a specific time and which species I want to look at
  # And I need to run a small simulation to get some ecotypes and then run it with the data function
  
  dimnames(RDD)$Species = sim@params@species_params$species # ecotypes have the same name
  dimnames(RDI)$Species = sim@params@species_params$species
  
  # choose a species
  i = 5
  ecoName = sim@params@species_params[sim@params@species_params$species == i,]$ecotype # name of the ecotypes in the species
  
  RDDSp = RDD[,which(colnames(RDD)==i)]
  colnames(RDDSp) = ecoName #values are isolated from other species and have the right names
  
  RDISp = RDI[,which(colnames(RDI)==i)]
  colnames(RDISp) = ecoName
  
  egg = apply(RDISp,2, function (x) x/rmax[i]) # x axis (rdi/rmax)
  reproduction = apply(RDDSp,2, function (x) x/rmax[i]) # y axis ( rdd/rmax)
  
  EGG = melt(egg)
  REPRO = melt(reproduction)
  
  graphdata = EGG
  graphdata$repro = REPRO$value # make only one dataframe
  
  #add the total spawn from the species
  RDItot = apply(RDISp,1,sum)
  RDDtot = apply(RDDSp,1,sum)
  
  # rmax with real data
  ggplot(graphdata)+
    geom_line(aes(x = value, y = repro, color = as.factor(Species)))+
    geom_hline(yintercept = 1)# rmax
  
  
  #example with fake values
  SumPar = sim@params@species_params
  truc = sim@n
  # put 0 in sim@n when w < w_mat
  for (i in 1:dim(truc)[1]) # for each time step
  {
    for (j in 1:dim(truc)[2]) # for each ecotypes
    {
      w_lim = sim@params@species_params$w_mat[j] # get the maturation size of the ecotype
      S <- numeric(length(sim@params@w))
      S[sapply(w_lim, function(i) which.min(abs(i - sim@params@w)))] <- 1 # find what w bin is the closest of the maturation size
      NoW_mat = which(S == 1) # what is the col number of that size bin
      truc[i,j,1:NoW_mat-1] <-0 # everything under this value become 0
    }
  }
  abundanceM = apply(truc, c(1,2),sum) # sum the abundance left 
  
  #2 normalisation per species 
  colnames(abundanceM) = sim@params@species_params$species
  abundanceNormal = matrix(0,nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])
  
  # I am getting rid of the species which went instinct at the begining
  SpIdx = NULL
  for (i in unique(sim@params@species_params$species))
    if (sum(abundanceM[,i]) != 0 & dim(SumPar[SumPar$species == i,])[1] != 1)
      SpIdx = c(SpIdx,i)
  
  #  I also need to get rid of the species that went extinct without having mutants (no trait variation)
  
  
  for (i in SpIdx)
  {
    abundanceSp = abundanceM # save to manip
    abundanceSp[,which(colnames(abundanceM) != i)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
    abundanceSp = sweep(abundanceSp,1,apply(abundanceSp,1,sum),"/") # normalise
    abundanceSp[is.nan(abundanceSp)] <-0
    abundanceNormal = abundanceNormal + abundanceSp # I just need to add them up to get the final matrix
  }
  colnames(abundanceNormal) = SumPar$ecotype
  LastAb = abundanceNormal[dim(abundanceNormal)[1],]# normalised abundance of the ecotypes at the last simulation step
  
  rmaxN = sim@params@species_params$r_max
  rmaxN = rmaxN * LastAb # normlised rmax following the ecotypes abundance
  # I have the rmax value for the last abundance step
  # fake rdi values
  rdi = seq(0,0.0001,0.00000001)
  RDI = matrix(rdi , length(rdi) , length(rmaxN) )
  colnames(RDI) = SumPar$ecotype
  
  RDD1 = sweep(RDI,2,-rmaxN) # rdi+rmax
  RDD2 = sweep(RDI,2,rmaxN,"*")  # rdi*rmax
  
  RDD = RDD2/RDD1 # reproduction with rmax applied (rmax is different for each ecotype)
  
  
  egg = sweep(RDI,2,rmaxN,"/") # x axis (rdi/rmax)
  reproduction = sweep(RDD,2,rmaxN,"/") # y axis ( rdd/rmax)
  
  #get rid of Nan / Inf
  i = 1
  while (i <= dim(egg)[2])
  {
    if (is.nan(egg[1, i]))
    {
      egg = egg[, -i]
      reproduction = reproduction[, -i]
    }
    i = i + 1
  }
  
  
  
  EGG = melt(egg)
  REPRO = melt(reproduction)
  
  graphdata = EGG
  dimnames(graphdata)[[2]] = list("col","species","rdi")
  graphdata$repro = REPRO$value # make only one dataframe
  graphdata$bloodline = sapply(graphdata[,2], function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  
  #shows all ecotypes
  ggplot(graphdata)+
    geom_point(aes(x = rdi, y = repro, color = as.factor(species)))+
    scale_x_log10()+
    geom_hline(yintercept = 1)# rmax
  
  #same color within species
  ggplot(graphdata)+
    geom_point(aes(x = rdi, y = repro, color = as.factor(bloodline), group = species))+
    scale_x_log10()+
    geom_hline(yintercept = 1)# rmax
  
  #only one specific species
  # select the species
  i = 3
  graphdataSp = graphdata[which(graphdata$bloodline==i),]
  
  ggplot(graphdataSp)+
    geom_point(aes(x = rdi, y = repro, color = as.factor(species)))+
    scale_x_log10()+
    geom_hline(yintercept = 1)
  
  
  #fisheries scenarios---------------------
  source("TBM1.r") # the model from mizer (more like a set up)
  source("model.r") # my model 
  ##   scenario 1 biggest species fished, selectivity above maturation size
  
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #dividing species between the gears (fished and non-fished)
 # other_gears <- w_inf >= 10000
  gear_names <- rep("None", no_sp)
  #gear_names[other_gears] <- "FishingStuff"
  
  #setting up knife edge
  knife_edges <- w_inf * 0.35 # slightly above maturation size
  # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
  #knife_edges[1:6] <-1e6
  knife_edges <- 1000

  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 20,
                    kappa = 1, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 90, 
                    effort = 0.2, knife_edge_size = knife_edges, gear_names = gear_names)
  sim = processing(output, plot = T, where =  paste(dir,"/scenario1",sep=""))
  gc()

  ##   scenario 2 biggest species fished, selectivity at maturation size
  
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #dividing species between the gears (fished and non-fished)
  other_gears <- w_inf >= 10000
  gear_names <- rep("None", no_sp)
  gear_names[other_gears] <- "FishingStuff"
  
  #setting up knife edge
  knife_edges <- w_inf * 0.25 # at maturation size
  # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
  knife_edges[1:6] <-1e6
  
  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 40,
                    kappa = 0.05, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 95, 
                    effort = 0.4, knife_edge_size = knife_edges, gear_names = gear_names)
  sim = processing(output, plot = F, where =  paste(dir,"/scenario2",sep=""))
gc()
  
  ##   scenario 3 small species fished, selectivity above maturation size
  
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #dividing species between the gears (fished and non-fished)
  # 100 to 1000
  other_gears <- w_inf <= 1000 & w_inf >=100
  gear_names <- rep("None", no_sp)
  gear_names[other_gears] <- "FishingStuff"
  
  #setting up knife edge
  knife_edges <- w_inf * 0.35 # above maturation size
  # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
  knife_edges[1:2] <-1e6
  knife_edges[6:9] <-1e6
  
  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 20,
                    kappa = 1, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 90, 
                    effort = 0.4, knife_edge_size = knife_edges, gear_names = gear_names)
  sim = processing(output, plot = T, where =  paste(dir,"/scenario3",sep=""))

  
  ##   scenario 4 small species fished, selectivity at maturation size
  
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #dividing species between the gears (fished and non-fished)
  # 100 to 1000
  other_gears <- w_inf <= 1000 & w_inf >=100
  gear_names <- rep("None", no_sp)
  gear_names[other_gears] <- "FishingStuff"
  
  #setting up knife edge
  knife_edges <- w_inf * 0.25 # at maturation size
  # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
  knife_edges[1:2] <-1e6
  knife_edges[6:9] <-1e6
  
  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 20,
                    kappa = 1, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 90, 
                    effort = 0.4, knife_edge_size = knife_edges, gear_names = gear_names)
  sim = processing(output, plot = T, where =  paste(dir,"/scenario4",sep=""))

  
  ##   scenario 5 everyone fished, selectivity above maturation size
  
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #dividing species between the gears (fished and non-fished)
  
  other_gears <-  w_inf >=50
  gear_names <- rep("None", no_sp)
  gear_names[other_gears] <- "FishingStuff"
  
  
  #setting up knife edge
  knife_edges <- w_inf * 0.35 # at maturation size
  # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
  knife_edges[1:2] <-1e6
  
  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 20,
                    kappa = 1, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 90, 
                    effort = 0.4, knife_edge_size = knife_edges, gear_names = gear_names)
  sim = processing(output, plot = T, where =  paste(dir,"/scenario5",sep=""))
  
  ##   scenario 6 everyone fished, selectivity at maturation size
  
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #dividing species between the gears (fished and non-fished)
  
  other_gears <-  w_inf >=50
  gear_names <- rep("None", no_sp)
  gear_names[other_gears] <- "FishingStuff"
  
  
  #setting up knife edge
  knife_edges <- w_inf * 0.25 # at maturation size
  # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
  knife_edges[1:2] <-1e6
  
  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 20,
                    kappa = 1, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 90, 
                    effort = 0.4, knife_edge_size = knife_edges, gear_names = gear_names)
  sim = processing(output, plot = T, where =  paste(dir,"/scenario6",sep=""))
  

  # does not happen anymore -----------------------------
  #checking for errors
  param = output[[2]]@species_params
  paramS = data.frame(param$species,param$ecotype,param$pop,param$extinct,param$run,param$error)
  
  #in case of umbrella
  output[[length(output)]] = NULL # delete last half sim
  param = NULL
  for (i in 1:length(output)) param = rbind(param,output[[i]]@params@species_params) # create the dataframe for species
  param <- param[order(param$ecotype, param$extinct, decreasing=TRUE),] 
  param <- param[!duplicated(param$ecotype),]
  SummaryParams = param[order(param$pop,param$ecotype),]
  FinalParam <- MizerParams(SummaryParams, min_w =0.001, max_w=10000 * 1.1, no_w = 100, min_w_pp = 1e-10, w_pp_cutoff = 0.5, n = 0.75, p=0.75, q=0.8, r_pp=4, kappa=0.1, lambda = 2.05) #create the mizer param from the dataframe
  result=list(output,FinalParam) # put it in the right disposition
  rm(output)
  sim = processing(result,plot = T, where = "/umbrella")
  