# This is the main script where you run the simulation (or experiment new things)
# to run the model, first set up your working directory in "dir"
# then you can use the function myModel to make the projection
# the basic options are to set up the number of species and the time of the projection (default of 10 and 100 respectively)
# the advanced options are mu (default = 2) the mutation rate, extinction (default = TRUE), 
# RMAX (default = FALSE) to enable/disable the Rmax assumption of Mizer (limitation of eggs number)
# There is one option called data (default = FALSE) that gives back other things instead of the simulation if I'm checking stuff (the energy allocation for example)
# After running the model you have to process the output with the "stitch" function and you can plot with "plotDynamics"

# For mutations: 4 options
# optMutant = M1: mutation are egg density dependent
# optMutant = M2: mutation is a commom rate between species
# optMutant = M3: user choose the time of apparition of the mutant, can choose specific mutant as well
# optMutant = M4: mutants can appear at the same time from different species
# optMutant = stuff: Default option, disable mutations

# For now just run the section "setting up little model" to get the simulation running!


#setting things up -----------------------
rm(list = ls())
dir <-"/data/home/romainf/romain"
setwd(dir)
library(ggplot2)#because always need these two
library(reshape2)
library(plyr)# for aaply
library(grid)# for grid.newpage (plotSummary)
library(abind) # to use abind (bind of arrays)
#library(rmarkdown)
library(RColorBrewer)
library(tictoc)

source("MizerParams-class.r") #to get the Constructor
source("selectivity_funcs.r") #to get the knife_edge function
source("methods.r") #I'm doing my own methods then!
source("summaryFunction.r") #to have all the GetSomething functions
source("plotFunction.r") #to draw the plots
source("TBM1.r") # the model from mizer (more like a set up)
source("model.r") # my model 
source("utility.r") # helpful functions

# parametrisation that works for 3sp en deterministic
file_name = "/eta5"

#asymptotic size
no_sp = 9
min_w_inf <- 10
max_w_inf <- 1e5
w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)

#fisheries
gear_names <- rep("FishingStuff", no_sp)
knife_edges <- w_inf * 0.25 

#other
t_max = 50
no_run = 80
effort = 0

for (i in 1:5)
{
  tic()
  cat(sprintf("Run number %g\n",i))
  path_to_save = paste(getwd(),file_name,"/init/run", i, sep = "")
  


sim <- myModel(no_sp = no_sp, t_max = t_max, no_run = no_run, OptMutant = "M2",w_pp_cutoff = 1e3, erepro = 1,
               min_w_inf = 10, max_w_inf = max_w_inf, RMAX = T, cannibalism = 0.5, #r_mult = 1e2,
               hartvig = T, f0 = 0.5, eta = 0.5,initTime = 1,
               effort = effort, #knife_edge_size = knife_edges, gear_names = gear_names, 
               save_it = T, path_to_save = path_to_save,
               print_it = T, normalFeeding = F, Traits = "eta")

rm(sim) # clean behind
for (j in 1:20) gc()
toc()
}

# starting from previous abundance

folder <- paste(getwd(),file_name,sep="")
initFolder <- paste(folder,"/init",sep="")
dirContent <- dir(initFolder)
no_run = 40
#sim normal
for (i in 1:length(dirContent))
{
  if (file.exists(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = "")))
    path_to_save <- paste(folder,"/normal/",dirContent[i],sep = "")
    
    output <- myModel(no_sp = no_sp, t_max = t_max, no_run = no_run, OptMutant = "M2",w_pp_cutoff = 1e3, erepro = 1,
                   min_w_inf = 10, max_w_inf = max_w_inf, RMAX = T, cannibalism = 0.5, #r_mult = 1e2,
                   hartvig = T, f0 = 0.5, eta = 0.5,initTime = 1, initCondition = sim,
                   effort = effort, #knife_edge_size = knife_edges, gear_names = gear_names, 
                   save_it = T, path_to_save = path_to_save,
                   print_it = T, normalFeeding = F, Traits = "eta")
    for (j in 1:20) gc()
  }
}
#sim fish
for (i in 1:length(dirContent))
{
  if (file.exists(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = "")))
    path_to_save <- paste(folder,"/fisheries/",dirContent[i],sep = "")
    
    output <- myModel(no_sp = no_sp, t_max = t_max, no_run = no_run, OptMutant = "M2",w_pp_cutoff = 1e3, erepro = 1,
                      min_w_inf = 10, max_w_inf = max_w_inf, RMAX = T, cannibalism = 0.5, #r_mult = 1e2,
                      hartvig = T, f0 = 0.5, eta = 0.5,initTime = 1, initCondition = sim,
                      effort = 0.8, knife_edge_size = knife_edges, gear_names = gear_names, 
                      save_it = T, path_to_save = path_to_save,
                      print_it = T, normalFeeding = F, Traits = "eta")
    for (j in 1:20) gc()
  }
}

file_name = "/eta5"
no_sp = 9
plotTraitOverlap(directory = paste(getwd(),file_name,sep=""),PPMR = F,Sig = F, SpIdx = seq(1,no_sp), init = T)
plotFitnessMultiOverlap(directory = paste(getwd(),file_name,sep=""),PPMR = F,Sig = F, SpIdx = seq(1,no_sp))
plotDynamicsMulti(folder = paste(getwd(),file_name,sep=""))

a <- PlotNoSp(sim)
ggsave(filename = "eta3.png",plot = a)

# one sp run
for (i in 2:5)
{
  path_to_save = paste(getwd(),"/oneSp/2/fisheries/run", i, sep = "")
  #asymptotic size
  no_sp = 2
  size = 1e2
  min_w_inf <- 10
  max_w_inf <- 1e2
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #fisheries
  gear_names <- rep("FishingStuff", no_sp)
  knife_edges <- w_inf * 0.25 
  
  output <- myModel(no_sp = no_sp, t_max = 100, no_run = 20, w_pp_cutoff = 1e4, kappa = 0.00005,
                    min_w_inf = min_w_inf, max_w_inf = max_w_inf,
                    effort = 0.8, knife_edge_size = knife_edges, gear_names = gear_names, 
                    print_it = T, normalFeeding = F, hartvig = T, Traits = F)
  
  sim = processing(output, plot = F, where =  path_to_save,save_it = F) #handle data and save not optimised output
  gc()
  simOpt = superOpt(sim) #optimise output
  save(simOpt,file = paste(path_to_save,"/run",".Rdata",sep="")) #save it
  rm(sim,output) # clean behind
  for (j in 1:20) gc()
}

plotTraitOverlap(directory = paste(dir,"/eta3",sep=""),SpIdx = seq(1,3),PPMR = F,Sig = F)


TotAnalysis(folder = paste(getwd(),"/biomassCompareNN1T/normal",sep=""))

dualAnalysis(folder = paste(getwd(),"/biomassCompareNN1T",sep=""),SpIdx = seq(2,9,1), biomass = F, PPMR = F, Sig = F, trait = T, fitness = F)

dualAnalysis(folder = paste(getwd(),"/biomCompareNN3T",sep=""),SpIdx = seq(2,9,1))

# multiple runs
file_name = "/eta2"

for (i in 1:1)
{
  tic()
  cat(sprintf("Run number %g\n",i))
  path_to_save = paste(getwd(),file_name,"/init/run", i, sep = "")
  #asymptotic size
  no_sp = 9
  min_w_inf <- 10
  max_w_inf <- 1e5
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  
  #fisheries
  gear_names <- rep("FishingStuff", no_sp)
  knife_edges <- w_inf * 0.25 
  
  output <- myModel(no_sp = no_sp, t_max = 50, no_run = 80, z0pre = 0, kappa = 0.01,
                    min_w_inf = min_w_inf, max_w_inf = max_w_inf, # cannibalism = 1,
                    effort = 0, #knife_edge_size = knife_edges, gear_names = gear_names, 
                    save_it = T, path_to_save = path_to_save, hartvig = T,
                    print_it = T, normalFeeding = F, Traits = "eta")
  
  #rm(output) # clean behind
  for (j in 1:20) gc()
  toc()
}
plotDynamics(output)
# starting from previous abundance

folder <- paste(getwd(),file_name,sep="")
initFolder <- paste(folder,"/init",sep="")
dirContent <- dir(initFolder)
#asymptotic size
no_sp = 9
min_w_inf <- 10
max_w_inf <- 1e5
w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)

#fisheries
gear_names <- rep("FishingStuff", no_sp)
knife_edges <- w_inf * 0.25 

#sim normal
for (i in 1:length(dirContent))
{
  if (file.exists(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = "")))
path_to_save <- paste(folder,"/normal/",dirContent[i],sep = "")
output <- myModel(no_sp = no_sp, t_max = 50, no_run = 40, z0pre = 0, #w_pp_cutoff = 1e4, kappa = 0.00005,
                  min_w_inf = min_w_inf, max_w_inf = max_w_inf, # cannibalism = 1,
                  effort = 0, #knife_edge_size = knife_edges, gear_names = gear_names, 
                  save_it = T, path_to_save = path_to_save, initCondition = sim,
                  print_it = T, normalFeeding = F, Traits = "eta")
for (j in 1:20) gc()
}
}
#sim fish
for (i in 1:length(dirContent))
{
  if (file.exists(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = ""))) 
  {
    sim <- get(load(paste(initFolder,"/",dirContent[i],"/run.Rdata",sep = "")))
    path_to_save <- paste(folder,"/fisheries/",dirContent[i],sep = "")
    output <- myModel(no_sp = no_sp, t_max = 50, no_run = 40,  z0pre = 0, #w_pp_cutoff = 1e4, kappa = 0.00005,
                      min_w_inf = min_w_inf, max_w_inf = max_w_inf, # cannibalism = 1,
                      effort = 0.8, knife_edge_size = knife_edges, gear_names = gear_names, 
                      save_it = T, path_to_save = path_to_save, initCondition = sim,
                      print_it = T, normalFeeding = F, Traits = "eta")
    for (j in 1:20) gc()
  }
}

file_name = "/eta4"
plotTraitOverlap(directory = paste(getwd(),file_name,sep=""),PPMR = F,Sig = F)
plotFitnessMultiOverlap(directory = paste(getwd(),file_name,sep=""),PPMR = F,Sig = F)
plotDynamicsMulti(folder = paste(getwd(),file_name,sep=""))

#with parallel
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

dir <-"/mnt/home/romain/mystuff"
setwd(dir)
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

# deterministic runs and fitness calcul

# get a deterministic run to try out

for (i in c(0.35,0.5))
{
dt = 0.1
no_sp = 3
t_max = 100
no_run = 5
effort = 0
sex_ratio = 0.5
laststep = t_max * no_run

sim <- myModel(no_sp = no_sp, t_max = t_max, no_run = no_run, OptMutant = "No",w_pp_cutoff = 1e3, erepro = 1,
               min_w_inf = 10, max_w_inf = 1e3, RMAX = F, cannibalism = 0.5, #r_mult = 1e2,
               hartvig = T, f0 = 0.5, eta = 0.35,initTime = 1,
               effort = effort, #knife_edge_size = knife_edges, gear_names = gear_names, 
               save_it = F, path_to_save = path_to_save,
               print_it = T, normalFeeding = F, Traits = "eta")


# plotCohort(sim,path_to_png = paste(getwd(),sep=""), t_start = 150)
# plotBiomass(sim,end_time = 250)
# fitness calcul
# I want the spawn through life -> last value of total cohortR
cohortSpan <- seq(100,1500,100)
no_sp = dim(sim@params@species_params)[1]
fitness <- array(0,c(no_sp, length(cohortSpan)), dimnames = list(dimnames(sim@n)$species,cohortSpan)) #collect the total spawn per time (start of cohort) per species
names(dimnames(fitness)) <- list("species","cohort")
for (cohortGen in cohortSpan)
{
  t_start = cohortGen
  # get rid of the non-existing species at that time
  SpIdx <- which(sim@n[t_start,,1]>0)
  # sim@n <- simulation@n[,SpIdx,]
  # no_sp = dim(sim@n)[2]
  
  cat(sprintf("cohort is %g\n",cohortGen))
  T = 5/dt; # number of time steps you want to follow cohort for
  #t_start = laststep - T; # setting the start time for tracking the cohort
  t_start = cohortGen
  cohortW = array(0, c(no_sp, T+1)); # row vector for following cohort weight
  cohortS = array(0, c(no_sp, T+1)); # vector for cohort survival
  cohortR = array(0, c(no_sp, T+1)); # vector for cohort spawning
  cohortR_sol = array(0, c(no_sp, T+1)); # vector for cohort spawn at size
  
  # NEWBORNS OVER LIFETIME
  cohortW[,1] = sim@params@w[1]; # log weight initially (newborn)
  cohortS[,1] = sim@n[t_start,,1]; # initial population in spectrum
  

  
  for (q in SpIdx){ 
    cat(sprintf("q is %g\n",q))
    for (t in seq(1,T)){ # within time period you're interested in
      #cat(sprintf("t is %g\n",t))
      
      cohortWprev = max(which(cohortW[q,t] - sim@params@w >= 0)) # weight bin of cohort from last time step 
      growth = getEGrowth(sim@params,n = sim@n[t_start+t-1,,],n_pp = sim@n_pp[t_start+t-1,])
      cohortW[q,t+1] = cohortW[q,t]+dt*growth[q,cohortWprev] # using growth rate in that bin to update to cohortW(t-t_start+1)
      z = getZ(object = sim@params, n = sim@n[t_start+t-1,,],n_pp = sim@n_pp[t_start+t-1,], effort = effort)
      cohortS[q,t+1] = cohortS[q,t]*exp(-dt*z[q,cohortWprev]) # updating amount surviving using death rate
      
      # need to prepare n for no NAN
      n = sim@n[t_start+t-1,,]*cohortS[,t]/cohortS[,1]
      n[!is.finite(n)] <- 1e-30
      
      e_spawning <- getESpawning(object = sim@params, n = n,n_pp = sim@n_pp[t_start+t-1,])
      e_spawning_pop <- apply((e_spawning*n),1,"*",sim@params@dw)
      rdi <- sex_ratio*(e_spawning_pop * sim@params@species_params$erepro)/sim@params@w[sim@params@species_params$w_min_idx] # need to get RDI that way so it is not summed up
      cohortR[q,t+1] = cohortR[q,t] + dt*rdi[cohortWprev,q]  #[q,cohortWprev]
      #cohortR_sol[q,t+1] = dt*rdi[cohortWprev,q] # do not sum the spawn so it is the spawn at time
    }
    fitness[q,which(cohortGen==cohortSpan)] = cohortR[q,T]
  }
}

# need to take the average fitness
#dimnames(fitness)$species <- do.call(cbind,lapply(as.numeric(strsplit(as.character(dimnames(fitness)$species), "")), function (x) x[1])) # the name of the phenotypes become their first digit
dimnames(fitness)$species <- sim@params@species_params$species #or I can do that ^^
SpIdx <- unique(sim@params@species_params$species)

newFit <- array(0,c(length(SpIdx), length(cohortSpan)), dimnames = list(SpIdx,cohortSpan))
for (i in SpIdx)
{
lol <- fitness[which(dimnames(fitness)$species == as.character(SpIdx[i])),]
lol[lol==0] <- NA
newFit[i,] <- apply(lol,2,mean,na.rm=T)
}

newFit[!is.finite(newFit)] <- 0

cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind


#plot_datF <- melt(fitness)
plot_datF <- melt(newFit)
colnames(plot_datF) <- list("species","cohort","fitness")

pF <- ggplot(plot_datF) +
  geom_line(aes(x=cohort,y=fitness,color = as.factor(species)),size =2) +
  scale_x_continuous(name = "Cohorts") +
  scale_y_continuous(name = "Fitness", trans = "log10")+
  scale_colour_manual(values=cbPalette)+ # colorblind
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
        legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

#pB <- plotBiomass(sim,start_time = 10.1,end_time = 300.1)
pB <- plotDynamics(sim)

png(filename = paste("fitness",i,".png",sep=""),width = 20,height = 20,units = "cm",res = 500)
grid.newpage() 
pushViewport(viewport(layout = grid.layout(nrow=2, ncol=1, 
                                           widths = unit(20, "cm"), 
                                           heights = unit(rep(10,2), "cm"))))

print(pF, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) 
print(pB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) 

dev.off()

}


# working on that right now -------------------

# plot fitness over time
where = paste(getwd(),"/biomassCompareNN1T/normal",sep="")
where = paste(getwd(),"/biomassCompareNN1T/fisheries",sep="")
plotFitnessMulti(folder = where, PPMR = F, Sig = F )

  plotFitnessTime <- function(folder, returnData = F, SpIdx = NULL, Mat = T, Sig = T, PPMR = T, comments = T, whatTime = NULL)
  {
    window = c(-0.25,0.25)
    if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))
    
    listFitnessPlots <-list() #stack the plots at different times
    for (i in whatTime)
      listFitnessPlots[[i]] <- plotFitnessMulti(folder = folder,Mat = Mat, PPMR = PPMR, Sig = Sig, returnData = T, whatTime = i, SpIdx = SpIdx)
    listFitnessPlots <- listFitnessPlots[lapply(listFitnessPlots,length)>0]
    
    SpIdx <- seq(1,9)
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
          scale_y_continuous(name = sprintf("Species %i",i))+
          scale_color_grey(name = "Species")+ # grey
          theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
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

  
time_range = c(seq(500.1,3900.1,500),3999.1)
where = paste(getwd(),"/biomassCompareNN1T/normal",sep="")
plotFitnessTime(folder = where, whatTime = time_range, PPMR = F, Sig = F)
where = paste(getwd(),"/biomassCompareNN1T/fisheries",sep="")
plotFitnessTime(folder = where, whatTime = time_range, PPMR = F, Sig = F)

# overlapping of traits

plotTraitOverlap <- function(directory, SpIdx = NULL, comments = T)
{
 
  if (is.null(SpIdx)) SpIdx = seq(1,9)
  
  # NO FISHERIES PART
  # need to load the data first
  normalList <- bunchLoad(folder = paste(directory,"/normal",sep=""))
  # get the plots
  normalTraitsList <- list()
  for (i in 1:length(normalList))
  normalTraitsList[[i]] <- plotTraitsMulti(object = normalList[[i]],PPMR = F,Sig = F,returnData = T)

speciesListN <- list()
for (j in SpIdx) # for every species
{
  speciesData <- NULL
  for (i in 1:length(normalTraitsList)) # at each time
  {
    if (!is.null(normalTraitsList[[i]][[1]][[j]]))
    {
      a <- ggplot_build(normalTraitsList[[i]][[1]][[j]]) # take the plot data
      a$data[[1]]$group <- i # change the group to the run number
      speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
    }
  }
  speciesListN[[j]] <- speciesData #this is a list of all the species at each time
}

# FISHERIES PART
# need to load the data first
fishList <- bunchLoad(folder = paste(directory,"/fisheries",sep=""))
# get the plots
fishTraitsList <- list()
for (i in 1:length(fishList))
  fishTraitsList[[i]] <- plotTraitsMulti(object = fishList[[i]],PPMR = F,Sig = F,returnData = T)

speciesListF <- list()
for (j in SpIdx) # for every species
{
  speciesData <- NULL
  for (i in 1:length(fishTraitsList)) # at each time
  {
    if (!is.null(fishTraitsList[[i]][[1]][[j]]))
    {
      a <- ggplot_build(fishTraitsList[[i]][[1]][[j]]) # take the plot data
      a$data[[1]]$group <- i # change the group to the run number
      speciesData <- rbind(speciesData,a$data[[1]]) #bind the same species at different time
    }
  }
  speciesListF[[j]] <- speciesData #this is a list of all the species at each time
}
  
# plot time
plotMatStore = list()
plotSigStore = list()
plotPPMStore = list()

window = c(-0.2,0.2)
if (comments) cat(sprintf("windows is set from %g to %g\n",window[1], window[2]))

temp = NULL
for (i in SpIdx) if (!is.null(speciesListN[[i]])) temp = c(temp,i) #update SpIdx if species are extinct at all times
SpIdx = temp

for (i in SpIdx) # do the plots for each species and traits
{
  if (comments) cat(sprintf("plot %i\n",i))
  #if (i == SpIdx[1]) title = c("a)","b)","c)") else title = NULL
  
  if (Mat){
    plotMatStore[[i]] <-ggplot(speciesListN[[i]])+
      geom_line(aes(x=x,y=y, group = group)) +#colour=as.factor(group))) +
      geom_line(data = speciesListF[[i]], aes(x=x,y=y, group = group), linetype = "dashed") +
      scale_x_continuous(name = NULL)+
      scale_y_continuous(name = sprintf("Species %i",i), limits = window)+
      scale_color_grey(name = "Species")+ # grey
      theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
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
path_to_png = paste(directory,"/Trait.png",sep="")

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
  
}

plotTraitOverlap(directory = folder, PPMR = F, Sig = F)
  

#fishereis mortality
source("methods.r") #I'm doing my own methods then!

p <- plotFM(sim)

plotFM <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
  survivant = object@params@species_params[object@params@species_params$extinct ==0,]$ecotype # get the non-extinct phen
  
   object@n = object@n[,as.character(survivant),] #get rid of the species extinct in the abundance data
  
  
  f_time <- getFMort(object, time_range=time_range, drop=FALSE)#, ...)
  f <- apply(f_time, c(2,3), mean)
  plot_dat <- data.frame(value = c(f), Species = dimnames(f)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
  #p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Total fishing mortality", lim=c(0,max(plot_dat$value)))

    p <- ggplot(plot_dat) + 
      geom_line(aes(x=w, y = value, group = Species)) +
      scale_x_continuous(name = "Size", trans="log10") + 
      scale_y_continuous(name = "Total fishing mortality", lim=c(0,max(plot_dat$value))) +
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
          legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
      ggtitle("Fisheries mortality")
  

}

ggsave("mort.png", plot = p)


#' @rdname getFMort-methods
#' @aliases getFMort,MizerParams,matrix-method
setMethod('getFMort', signature(object='MizerParams', effort='matrix'),
          function(object, effort, ...){
            fMortGear <- getFMortGear(object, effort, ...)
            fMort <- apply(fMortGear, c(1,3,4), sum)
            return(fMort)
          })

#' @rdname getFMort-methods
#' @aliases getFMort,MizerSim,missing-method
setMethod('getFMort', signature(object='MizerSim', effort='missing'),
          function(object, effort, time_range=dimnames(object@effort)$time, drop=TRUE, ...){
            
            
            time_elements <- get_time_elements(object,time_range, slot="effort")
            
            fMort <- getFMort(object@params, object@effort, ...)
            
            return(fMort[time_elements,,,drop=drop])
          })


get_time_elements <- function(sim,time_range,slot_name="n"){
  if (!(slot_name %in% c("n","effort")))
    stop("'slot_name' argument should be 'n' or 'effort'")
  if (!is(sim,"MizerSim"))
    stop("First argument to get_time_elements function must be of class MizerSim")
  time_range <- range(as.numeric(time_range))
  # Check that time range is even in object
  sim_time_range <- range(as.numeric(dimnames(slot(sim,slot_name))$time))
  if ((time_range[1] < sim_time_range[1]) | (time_range[2] > sim_time_range[2]))
    stop("Time range is outside the time range of the modell")
  time_elements <- (as.numeric(dimnames(slot(sim,slot_name))$time) >= time_range[1]) & (as.numeric(dimnames(slot(sim,slot_name))$time) <= time_range[2])
  names(time_elements) <- dimnames(slot(sim,slot_name))$time
  return(time_elements)
}
# number of phenotypes per species through time --------------------
ID = sim@params@species_params
# I will calculate the number of apparition and number of extinction per species and make the difference to have the number of phen alive at time
# I need to calculate at each time step
dt = 0.1
for (i in unique(ID$species))
{
  a = matrix(data = NA, nrow = (ID$timeMax[1]), ncol = 3, dimnames = list(seq(1:(ID$timeMax[1])), c("pop","ext","current")))
  
  for (j  in 1:dim(a)[1])
  {
    if (any(rownames(a)[j] == ID$pop))
  }}

PlotNoSp <- function(object, print_it = T, species = F){
  SumPar = cbind(object@params@species_params$pop,object@params@species_params$extinct) # weird things happen without the as.numeric
  colnames(SumPar) = c("Apparition","Extinction")
  rownames(SumPar) = object@params@species_params$ecotype
  SumPar = SumPar[order(SumPar[,1],decreasing=FALSE),]
  SumPar = as.data.frame(SumPar) # little dataframe with species apparition and extinction
  for (i in 1:dim(SumPar)[1]) if (SumPar$Extinction[i] == 0) SumPar$Extinction[i] = object@params@species_params$timeMax[1] # the not extinct ones get the end sim as extinction value
  
  avg = matrix(data = 0, nrow = object@params@species_params$timeMax[1], ncol =4)
  avg[,1] = seq(1:object@params@species_params$timeMax[1])*0.1
  ApCount = length(which(SumPar$Apparition == 0)) # number of starting species
  ExCount = 0
  SumParEx = SumPar[order(SumPar$Extinction,decreasing=F),] #need this order for extinction in loop
  
  if (species)
  {
    # work to do here
    
  }
  
  
  
  for (i in 1:object@params@species_params$timeMax[1]) # for each time step
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
    geom_smooth(aes(x = time, y = number, color = "green")) +
    geom_smooth(aes(x=time, y = apparition, color = "blue"))+
    geom_smooth(aes(x=time, y = extinction, color = "red"))+
    scale_x_continuous(name = "Time (yr)") +
    scale_y_continuous(name = "Number of phenotypes")+
    scale_colour_discrete(labels = c("Total","Alive","Extinct"))+
    theme(legend.title=element_blank(),
          legend.position=c(0.19,0.95),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    ggtitle("Variation of phenotype's number throughout the simulation")
  
  if (print_it)  print(p)
  
  return(p)
}




## Fisheries scenario with all species fished, selectivity at maturation size -------------------------

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
knife_edges <- w_inf * 0.35 # at maturation size
# fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
#knife_edges[1:6] <-1e6

output <- myModel(no_sp = no_sp, t_max = 100, no_run = 40,
                  kappa = 0.05, min_w_inf = min_w_inf, max_w_inf = max_w_inf, h = 95, 
                  effort = 0.2, knife_edge_size = knife_edges, gear_names = gear_names, print_it = T)
sim = processing(output, plot = F, save_it = T)
#sim = processing(output, plot = F, where =  paste(dir,"/scenario2",sep=""))

gc()
#parallel -----------------------------------------------------
# with foreach
library(ggplot2)#because always need these two
library(reshape2)
library(plyr)# for aaply
library(grid)# for grid.newpage (plotSummary)
library(abind) # to use abind (bind of arrays)
library(rmarkdown)
library(RColorBrewer)

library(foreach)
library(iterators)
library(snow)
library(doSNOW)

cl<-makeCluster(5) # number of CPU cores
registerDoSNOW(cl)

op <- foreach(i=1:5, .packages = c("ggplot2","plyr","reshape2","grid","abind","RColorBrewer")) %dopar% {
  
  #setting things up -----------------------
  
  dir <-"/mnt/home/romain/mystuff"
  setwd(dir)
  source("MizerParams-class.r") #to get the Constructor
  source("selectivity_funcs.r") #to get the knife_edge function
  source("methods.r") #I'm doing my own methods then!
  source("summaryFunction.r") #to have all the GetSomething functions
  source("plotFunction.r") #to draw the plots
  source("TBM1.r") # the model from mizer (more like a set up)
  source("model.r") # my model 
  source("utility.r") # helpful functions
  # Run the model and plots------------------------------
  
  
  output <- myModel(no_sp = 9, t_max = 50, mu = 5, OptMutant = "M2", no_run = 2,
                    kappa = 0.05, ks=2, min_w_inf = 10, max_w_inf = 100000, h= 95,
                    #rm =1,
                    effort = 0,print_it = F)
  # path_to_save = paste(dir,"40run", i, sep = "")
  sim = processing(output, plot = F, save_it = F)
  
}

#save the output
dir <-"/mnt/home/romain/mystuff"
source("utility.r")
for (i in 1:length(op))
  {
  sim = superOpt(op[[i]])
  path_to_save = paste(dir,"/40run", i, sep = "")
  ifelse(!dir.exists(file.path(path_to_save)), dir.create(file.path(path_to_save)), FALSE) #create the file if it does not exists
  
  save(sim,file = paste(path_to_save,"/run",".Rdata",sep=""))
}

#with parallel
library(parallel)
library(ggplot2)#because always need these two
library(reshape2)
library(plyr)# for aaply
library(grid)# for grid.newpage (plotSummary)
library(abind) # to use abind (bind of arrays)
library(rmarkdown)
library(RColorBrewer)

dir <-"/mnt/home/romain/mystuff"
setwd(dir)
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

#unsure  what this setting does
options(warn=-1) #?

#Adjust this for num of targeted cpu/cores
# e.g. Numcores = detectCores()-1

numcores=5
cl <- makeForkCluster(getOption("cl.cores", numcores))
clusterApplyLB(cl
               ,x=1:5
               ,fun=multiRun
               ,no_sp = 9
               ,t_max = 100
               ,mu = 5
               ,no_run = 40
               ,min_w_inf = 10
               ,max_w_inf = 10e5
               ,effort = 0.2
               ,fisheries = T)
stopCluster(cl)


#(optional) compare end with start time, for timiing
print((proc.time()-ptm)/60.0)


# plot at different time to compare ------------------
res = 800

plotDynamics(sim, species = T)
setwd(paste(dir,"feedStat", sep = ""))
mytitle = "biomass.png"
dev.print(png, mytitle, width = res, height = 0.6*res) 

plotSS(sim)
setwd(paste(dir,"feedStat", sep = ""))
mytitle = "sizespectrum.png"
dev.print(png, mytitle, width = res, height = 0.6*res) 


when = c(100,500,1000,1500)
for(i in when)
{
  plotFood(sim, time_range = i)
  setwd(paste(dir,"feedStat", sep = ""))
  mytitle = paste("feeding t= ",i,".png",sep="")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  plotUdead(sim, time_range = i)
  setwd(paste(dir,"feedStat", sep = ""))
  mytitle = paste("mortality t= ",i,".png",sep="")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
  
  plotGrowth(sim, time_range = i)
  setwd(paste(dir,"feedStat", sep = ""))
  mytitle = paste("growth t= ",i,".png",sep="")
  dev.print(png, mytitle, width = res, height = 0.6*res) 
}

# automatisation of the results --------------------


for (i in 1:30)
{
  
  dir <-"/home/romain/mystuff"
  setwd(dir)
  library(ggplot2)#because always need these two
  library(reshape2)
  library(plyr)# for aaply
  library(grid)# for grid.newpage (plotSummary)

  library(abind) # to use abind (bind of arrays)
  library(rmarkdown)
  library(RColorBrewer)

  source("MizerParams-class.r") #to get the Constructor
  source("selectivity_funcs.r") #to get the knife_edge function
  source("methods.r") #I'm doing my own methods then!
  source("summaryFunction.r") #to have all the GetSomething functions
  source("plotFunction.r") #to draw the plots
  source("TBM1.r") # the model from mizer (more like a set up)
  source("model.r") # my model 
  source("utility.r") # helpful functions
  gc()
  
  output <- myModel(no_sp = 9, t_max = 100, mu = 5, OptMutant = "M2", no_run = 50,
                    kappa = 0.025, ks=2, min_w_inf = 10, max_w_inf = 100000, h= 85,
                    #rm =1,
                    effort = 0)
  
  here = paste(dir,"/","50run",i,sep="")
 
  processing(output, plot = T, where = here)
  rm(list = ls())
}


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
  
    
    




 # paper plots ----------------------------
 # need some data to play with
 res = 650 # figure resolution
 subdir = "/paper" # where to store the plots
 t_max = 100
 no_sp = 9
 mu = 5
 no_run = 4
 dt = 0.1
 # first run with mutants
 output <- myModel(no_sp = no_sp, t_max = t_max, mu = mu, OptMutant = "M2", RMAX = TRUE, kappa = 0.008, cannibalism = 1, r_mult = 1e0,
                   erepro = 0.1, min_w_inf = 10, p =0.75, ks=2, max_w_inf = 10000, extinct = TRUE, Traits = TRUE, 
                   hartvig = TRUE, no_run = no_run, effort = 0)#,rm =1)
 #checking for errors
 param = output[[2]]@species_params
 sum(param$error)
 paramS = data.frame(param$species,param$ecotype,param$pop,param$extinct,param$run,param$error)
 
 simMute = processing(output, plot = F)
 

# then 1 run to get the data of everything, mutations not allowed
output <- myModel(no_sp = 9, t_max = 50, mu = 5, OptMutant = "yo",  
                     #min_w_inf = 10, 
                      #rm = i,
                     #max_w_inf = 10, 
                      no_run = 1, effort = 0, k0 = 25, r_pp = 4, kappa = 0.008, 
                     data = TRUE, param = simMute@params,
                     ks = 2)
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
   
   ifelse(!dir.exists(file.path(dir, subdir)), dir.create(file.path(dir, subdir)), FALSE) #create the file if it does not exists
   
   # dynamics plots----------------------------
   # Biomass
   biomass <- getBiomass(sim) # n * w * dw and sum by species
   Sbm <- melt(biomass) # melt for ggplot
   names(Sbm) = list("time","sp","value") #,"bloodline") # just in case
   min_value <- 1e-300
   Sbm <- Sbm[Sbm$value >= min_value,]
   colourCount = length(unique(Sbm$sp))
   getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used
   Sbm$bloodline = sapply(Sbm[,2], function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
   
  # Abm = Sbm[Sbm$bloodline == 1,] to plot a specific sepcies
   
       ggplot(Sbm) +
       geom_line(aes(x = time, y = value, colour = as.factor(bloodline), group = sp)) +
       scale_y_log10(name = "Biomass in g.m^-3", limits = c(1e-13, NA), breaks = c(1 %o% 10^(-11:-1))) +
       scale_x_continuous(name = "Time in years") +
       labs(color='Species') +
       theme(legend.key = element_rect(fill = "white"))+
       theme_bw()+
       scale_color_grey()

   setwd(paste(dir,subdir, sep = "")) #to have the figures in the right directory
   mytitle = paste("biomass",".png", sep = "")
   dev.print(png, mytitle, width = res, height = 0.6*res) 
   
# size spectrum
   time_range = max(as.numeric(dimnames(sim@n)$time))
   min_w = min(sim@params@w) / 100
   time_elements <- get_time_elements(sim, time_range)
   spec_n <- apply(sim@n[time_elements, , , drop = FALSE], c(2, 3), mean)
   background_n <- apply(sim@n_pp[89.2:99.2, , drop = FALSE], 2, mean)
   spec_n <- sweep(spec_n, 2, sim@params@w, "*")
   background_n <- background_n * sim@params@w_full
   y_axis_name = "Biomass"
    
     # Make data.frame for plot
     plot_datSP <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(sim@params@w, each=nrow(sim@params@species_params)), bloodline = sim@params@species_params$species)
     plot_datPkt <- data.frame(value = c(background_n), Species = "Background", w = sim@params@w_full)
     # lop off 0s in background and apply min_w
     plot_datSP <- plot_datSP[(plot_datSP$value > 0) & (plot_datSP$w >= min_w),]
     plot_datPkt <- plot_datPkt[(plot_datPkt$value > 0) & (plot_datPkt$w >= min_w),]
     getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used
     
     
  ggplot(plot_datSP) + 
       geom_line(aes(x=w, y = value, colour = as.factor(bloodline), group = Species)) + 
       geom_line(data = plot_datPkt, aes(x = w, y = value, colour = Species), size = 1.5) +
       scale_x_log10(name = "Size in g", breaks = c(1 %o% 10^(-6:5)))+
       scale_y_log10(name = "Abundance density in individuals.m^-3", limits = c(1e-27,1e4)) +
       #scale_color_discrete(name = "Species")+
       theme(panel.background = element_blank(),legend.key = element_rect(fill = "white"))+
       theme_bw()+
       scale_color_grey(name = "Species")
   
     
   setwd(paste(dir,subdir, sep = "")) #to have the figures in the right directory
   mytitle = paste("sizespectrum_",".png", sep = "")
   dev.print(png, mytitle, width = res, height = 0.6*res) 
   
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

  # traits dynamics ----------------
  # Traits plots
  # little initialisation to plot the trait values
  SumPar = sim@params@species_params
  TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),SumPar$pop,SumPar$extinct,SumPar$w_mat,SumPar$beta,SumPar$sigma) # weird things happen without the as.numeric
  colnames(TT) = c("Lineage","Ecotype","Apparition","Extinction","Maturation_size","PPMR","Diet_breadth")
  rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  TT = as.data.frame(TT)
  
  for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = SumPar$timeMax[i]
  
  
  # I'm going to do a weighted mean so I need to multiply each trait by their abundance proportion to the other (sum of the weights equal 1)
  # 1 matrix of summed abundance of mature ind at each time step
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
  #colnames(abundanceNormal) = sim@params@species_params$ecotype
  
  # Maturation size---------------
  #3 Multiply by their trait value
  abundanceT = sweep(abundanceNormal,2,sim@params@species_params$w_mat,"*")
  
  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(sim@params@species_params$species)), dimnames = list(rownames(abundanceT),unique(sim@params@species_params$species)))
  names(dimnames(TotMean)) = list("time","species")
  
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
    TotMean[,i] = AMean
  }
  
  # Calculate variance and standard deviation
  # it is the sum of the difference between value and mean squared and multiplied by the weight
  # I think I have to do it for each species from there
  
  statMS = list() # list with the stats of all species
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = sim@params@species_params$w_mat[sim@params@species_params$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(rownames(abundanceT),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","variance","sd"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
      stat[j,3] = variance
      stat[j,4] = sqrt(variance) # calculate the standard deviation
    }
    statMS[[i]] = as.data.frame(stat) # put in the list
  }
  
  plotMatStore <- list()
  
  for (i in SpIdx)
  {
    # data for the trait values
    # empty matrix of ecotype of species i by time
    a = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:SumPar$timeMax[1])))
    # fill the matrix with ones when the ecotype exists
    for (x in 1:nrow(a)) # I'm sure I can do an apply but don't know how
    {
      for (j in 1:ncol(a))
      {
        if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) a[x,j] = 1
      }}
    # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
    # calculate the mean trait value + standard deviation
    # change the ones by the trait value of the ecotype
    b = a * TT[TT$Lineage == i,]$Maturation_size
    
    # calculate mean trait value at each time step (normal mean)
    # c = apply(a,2,sum) # this vector is the number of traits present at each time step
    # d = apply(b,2,sum) # this vector is the sum of the traits value at each time step
    # e = d/c # this is the mean trait value at each time step
    # f = apply((sweep(b,2,e,"-")*a)^2,2,sum)/c # this is the variance at each time step
    # g = sqrt(f) # this is the standard population deviation at each time step
    # 
    # stat = data.frame(e,g,e-g,e+g)
    # colnames(stat) = list("mean","variance","low","up")
    
    #X = TT$Maturation_size[TT$Ecotype==i] # value around which to do the breaks
    
    stat = statMS[[i]] # take the stats of the right species
    stat$mean <- as.numeric(substr(stat$mean,1,20)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5
    stat$sd <- as.numeric(substr(stat$sd,1,20)) # otherwise I get errors from ggplot
    stat$up = stat$mean + stat$sd
    stat$low = stat$mean - stat$sd
    
    family = TT[TT$Lineage == i,3:5]
    Fam = melt(family,"Maturation_size") 
    name = paste("Species",i, sep = " ")
    
    if (i == 1) # top of the multi plot
    {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) +
        geom_line(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(breaks = NULL, name = NULL) +
        scale_y_continuous(name = name) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle("Maturation size (g)")
      
    }
    
    else if (i == 9) # bottom of multiplots
    
    {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) +
        geom_line(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(name = "Time (yr)") +
        scale_y_continuous(name = name) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
      
    }
    
    else {
    p = ggplot() +
      geom_point(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) +
      geom_line(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) + 
      geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
      geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
      geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
      scale_x_continuous(breaks = NULL, name = NULL) +
      scale_y_continuous(name = name) +
      scale_colour_discrete(labels = c("mean","standard deviation"))+
      theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
      guides(color=guide_legend(override.aes=list(fill=NA)))
    }
    # print(p)
    # setwd(paste(dir,subdir, sep = ""))
    # mytitle = paste(name,".png", sep = "")
    # dev.print(png, mytitle, width = res, height = 2/3* res) 
    plotMatStore[[i]] = p
  }
  
  # PPMR --------------------------
  #3 Multiply by their trait value
  abundanceT = sweep(abundanceNormal,2,sim@params@species_params$beta,"*")
  
  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(sim@params@species_params$species)), dimnames = list(rownames(abundanceT),unique(sim@params@species_params$species)))
  names(dimnames(TotMean)) = list("time","species")
  
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean
    TotMean[,i] = AMean
  }
  
  # Calculate variance and standard deviation
  # it is the sum of the difference between value and mean squared and multiplied by the weight
  # I think I have to do it for each species from there
  
  statPPMR = list() # list with the stats of all species
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = sim@params@species_params$beta[sim@params@species_params$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(rownames(abundanceT),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","variance","sd"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
      stat[j,3] = variance
      stat[j,4] = sqrt(variance) # calculate the standard deviation
    }
    statPPMR[[i]] = as.data.frame(stat) # put in the list
  }
  
  plotPPMRStore <- list()
  
  for (i in SpIdx)
  {
    # empty matrix of ecotype of species i by time
    a = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
    # fill the matrix with ones when the ecotype exists
    for (x in 1:nrow(a)) # I'm sure I can do an apply but don't know how
    {
      for (j in 1:ncol(a))
      {
        if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) a[x,j] = 1
      }}
    # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
    # change the ones by the trait value of the ecotype
    b = a * TT[TT$Lineage == i,]$PPMR
    # calculate mean trait value at each time step
    # c = apply(a,2,sum) # this vector is the number of traits present at each time step
    # d = apply(b,2,sum) # this vector is the sum of the traits value at each time step
    # e = d/c # this is the mean trait value at each time step
    # f = apply((sweep(b,2,e,"-")*a)^2,2,sum)/c # this is the variance at each time step
    # g = sqrt(f) # this is the standard population deviation at each time step
    # 
    # stat = data.frame(e,g,e-g,e+g)
    # colnames(stat) = list("mean","variance","low","up")
    
    stat = statPPMR[[i]] # take the stats of the right species
    stat$mean <- as.numeric(substr(stat$mean,1,10)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5
    stat$sd <- as.numeric(substr(stat$sd,1,10)) # otherwise I get errors from ggplot
    stat$up = stat$mean + stat$sd
    stat$low = stat$mean - stat$sd
    
    family = TT[TT$Lineage == i,c(3,4,6)]
    Fam = melt(family,"PPMR") 
    name = paste("PPMR_SP",i, sep = "")
   
    if (i == 1) # top of the multi plot
    {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = PPMR, group = PPMR)) +
        geom_line(data = Fam, aes(x = value, y = PPMR, group = PPMR)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(breaks = NULL, name = NULL) +
        scale_y_continuous(limits = c(80,120), name = NULL) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle("Preferred PPMR")
      
    }
    
    else if (i == 9) # bottom of multiplots
      
    {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = PPMR, group = PPMR)) +
        geom_line(data = Fam, aes(x = value, y = PPMR, group = PPMR)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(name = "Time (yr)") +
        scale_y_continuous(limits = c(80,120), name = NULL) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
      
    }
    
    else {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = PPMR, group = PPMR)) +
        geom_line(data = Fam, aes(x = value, y = PPMR, group = PPMR)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(breaks = NULL, name = NULL) +
        scale_y_continuous(limits = c(80,120), name = NULL) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
    }
    # print(p)
    # setwd(paste(dir,subdir, sep = ""))
    # mytitle = paste(name,".png", sep = "")
    # dev.print(png, mytitle, width = res , height = 2/3* res) 
    
    plotPPMRStore[[i]] <- p
    
  }
  
  # Diet breath -------------------------
  #3 Multiply by their trait value
  abundanceT = sweep(abundanceNormal,2,sim@params@species_params$sigma,"*")
  
  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = length(unique(sim@params@species_params$species)), dimnames = list(rownames(abundanceT),unique(sim@params@species_params$species)))
  names(dimnames(TotMean)) = list("time","species")
  
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean
    TotMean[,i] = AMean
  }
  
  # Calculate variance and standard deviation
  # it is the sum of the difference between value and mean squared and multiplied by the weight
  # I think I have to do it for each species from there
  
  statDB = list() # list with the stats of all species
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = sim@params@species_params$sigma[sim@params@species_params$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(rownames(abundanceT),meanSp,0,0), ncol = 4,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","variance","sd"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,]) # calculate the variance
      stat[j,3] = variance
      stat[j,4] = sqrt(variance) # calculate the standard deviation
    }
    statDB[[i]] = as.data.frame(stat) # put in the list
  }
  
  plotDBStore <- list()
  
  
  for (i in SpIdx)
  {
    # empty matrix of ecotype of species i by time
    a = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
    # fill the matrix with ones when the ecotype exists
    for (x in 1:nrow(a)) # I'm sure I can do an apply but don't know how
    {
      for (j in 1:ncol(a))
      {
        if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) a[x,j] = 1
      }}
    # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
    # change the ones by the trait value of the ecotype
    b = a * TT[TT$Lineage == i,]$Diet_breadth
    # calculate mean trait value at each time step
    # c = apply(a,2,sum) # this vector is the number of traits present at each time step
    # d = apply(b,2,sum) # this vector is the sum of the traits value at each time step
    # e = d/c # this is the mean trait value at each time step
    # f = apply((sweep(b,2,e,"-")*a)^2,2,sum)/c # this is the variance at each time step
    # g = sqrt(f) # this is the standard population deviation at each time step
    # 
    # stat = data.frame(e,g,e-g,e+g)
    # colnames(stat) = list("mean","variance","low","up")
    
    stat = statDB[[i]] # take the stats of the right species
    stat$mean <- as.numeric(substr(stat$mean,1,10)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5
    stat$sd <- as.numeric(substr(stat$sd,1,10)) # otherwise I get errors from ggplot
    stat$up = stat$mean + stat$sd
    stat$low = stat$mean - stat$sd
    
    family = TT[TT$Lineage == i,c(3,4,7)]
    Fam = melt(family,"Diet_breadth") 
    #name = paste("DB_SP",i, sep = "")
    
    
    if (i == 1) # top of the multi plot
    {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) +
        geom_line(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red"), method = "auto") +
       geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue"), method = "auto") +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue"), method = "auto") +
        scale_x_continuous(breaks = NULL, name = NULL) +
        scale_y_continuous(limits = c(0.75, 1.4), name = NULL) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend
        guides(color=guide_legend(override.aes=list(fill=NA)))+
        ggtitle("Diet Breadth")
      
    }
    
    else if (i == 9) # bottom of multiplots
      
    {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) +
        geom_line(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(name = "Time (yr)") +
        scale_y_continuous(limits = c(0.75, 1.4), name = NULL) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
      
    }
    
    else {
      p = ggplot() +
        geom_point(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) +
        geom_line(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) + 
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
        geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
        scale_x_continuous(breaks = NULL, name = NULL) +
        scale_y_continuous(limits = c(0.75, 1.4), name = NULL) +
        scale_colour_discrete(labels = c("mean","standard deviation"))+
        theme(legend.title=element_blank(),panel.background = element_blank(), legend.position="none", legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
    }
    
    
    # print(p)
    # setwd(paste(dir,subdir, sep = ""))
    # mytitle = paste(name,".png", sep = "")
    # dev.print(png, mytitle, width = res, height = 2/3* res)
    
    plotDBStore[[i]] <- p
    
  }
  
  # multi plots ----------------------------------
  
  multiplot(plotMatStore[[1]], plotMatStore[[3]],plotMatStore[[4]],plotMatStore[[5]],plotMatStore[[9]], 
            plotPPMRStore[[1]], plotPPMRStore[[3]],plotPPMRStore[[4]], plotPPMRStore[[5]],plotPPMRStore[[9]],
            plotDBStore[[1]],plotDBStore[[3]],plotDBStore[[4]], plotDBStore[[5]],plotDBStore[[9]], 
            cols=3)
  
res = 800
  mytitle = paste("trait",".png", sep = "")
  dev.print(png, mytitle, width = res, height = 2/3* res)
  
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
  
  
  #fisheries run---------------------
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
  