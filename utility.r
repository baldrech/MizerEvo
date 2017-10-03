# useful function


# give it an object and it will return the biomass data ready to be plotted (for species and phenotypes)
biom <-function(object,phenotype=T)
{
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
  plotB <- function(x)
  {
    Biom <- melt(x) # melt for ggplot
    names(Biom) = list("time","sp","value")
    # Due to log10, need to set a minimum value, seems like a feature in ggplot
    min_value <- 1e-300
    Biom <- Biom[Biom$value >= min_value,]
    # take the first digit of the species column and put it in a new column
    Biom$bloodline = sapply(Biom[,2], function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
    return(Biom)
  }
  BiomPhen = NULL
  if (phenotype) BiomPhen <- plotB(biomass)
  BiomSp <- plotB(biomassTot)
  
  return(list(BiomSp,BiomPhen))
}

# function that takes the output of the model and make it usable for plotting and shit
finalTouch <- function(result, dt = 0.1, print_it = T)
{
  ## processing data
  if (print_it) cat(sprintf("Data handling\n"))
  # a result will be a list of simulations (number of runs) and the ID card of the ecosystem, 
  
  gc()
  sim = result[[1]] # get the dat
  SummaryParams = result[[2]] # get the params
  rm(result) # get space back
  sim <- sim[lapply(sim,length)>0] # if a sim is empty
  endlist = length(sim) # updating the number
  template = sim[[endlist]] # to keep a template of mizer object somewhere
  gc()
  
  # stitiching the sims together
  Dtime = dim(sim[[1]]@n)[1] # need these dimensions for making new arrays but separatly as I'm going to work on sp dim
  Dsp = length(SummaryParams@species_params$ecotype)
  Dw = dim(sim[[1]]@n)[3]
  # put all the sim at the same dimension
  for (i in 1:endlist) 
  {
    # print(i)
    for (j in 1:Dsp)
    {
      gc()
      if (is.na(dimnames(sim[[i]]@n)$sp[j])) # if I need to add species at the end of the array
      {
        sim[[i]]@n <- abind(sim[[i]]@n,array(dim = c(Dtime,Dw)),along = 2)
        names(dimnames(sim[[i]]@n)) <- list("time","sp","w") # the abind make me loose the info
      }
      else if (dimnames(sim[[i]]@n)$sp[j] != SummaryParams@species_params$ecotype[j]) # if I need to add in the middle of the array
      {
        rmbr = NULL # when I abind only one column of the array (happens when there is an extinction right before the last ecotype), I lose its name, leading to bugs
        
        if (j == dim(sim[[i]]@n)[2]) rmbr = dimnames(sim[[i]]@n)$sp[j] # with this I know when it happens
        
        sim[[i]]@n <- abind(sim[[i]]@n[,1:j-1,],array(dim = c(Dtime,Dw)),sim[[i]]@n[,j:dim(sim[[i]]@n)[2],], along = 2) # if sp 1 goes extinct I will have a bug, or not it seems
        names(dimnames(sim[[i]]@n)) <- list("time","sp","w") # the abind makes me loose the info
        
        if (is.null(rmbr)==F) dimnames(sim[[i]]@n)$sp[j+1] = rmbr # and I add back the name
        
      }
    }
    dimnames(sim[[i]]@n)$sp = SummaryParams@species_params$ecotype # put back the right names
    sim[[i]]@n = sim[[i]]@n[-1,,] # in my tentative to include fisheries, I am deleting the first time step (which is 0 now) to have the same length as the effort matrix. I'm running out of time, Ill fix that later I hope
  }
  # getting rid of the NAs
  gc()
  for (j in 1:endlist) sim[[j]]@n[is.na( sim[[j]]@n)] <- 0
  gc()
  # now add them in a row
  biom <- do.call(abind, c(lapply(sim, function(isim) isim@n),along = 1))
  gc()
  
  names(dimnames(biom)) = list("time","species","size")
  dimnames(biom)$time = seq(1, SummaryParams@species_params$timeMax[1]) * dt # meh
  gc()
  
  # I have to do the phyto aussi
  phyto = sim[[1]]@n_pp[-1, ] #get rid of time 0
  for (i in 2:endlist) phyto = abind(phyto, sim[[i]]@n_pp[-1, ], along = 1)
  
  if (print_it) cat(sprintf("Biomass handeled successfully, now starting the effort.\n"))
  
  # taking care of the effort
  # for now it is constant in time so we just need to multiply effort by no_run
  effort = sim[[1]]@effort
  for(i in 2:endlist) effort = rbind(effort,sim[[i]]@effort)
  names(dimnames(effort)) = list("time","effort")
  dimnames(effort)$time = seq(1, SummaryParams@species_params$timeMax[1]) * dt # i don't know whats wrong with the else
  
  # I need to assemble the mizer object now
  sim = template
  sim@params=SummaryParams
  sim@n = biom
  sim@effort = effort
  sim@n_pp = phyto
  
  rm(list = "template","biom","phyto","effort")
  gc()
  
  return(sim)  
}

# #function that prep files for download
# DLfig <- function(folder)
# {
#   dirContent = dir(folder)
#   runFolder = dirContent
#   for (i in 1:length(dirContent))
#     if (runFolder[i] == "fig" || runFolder[i] == "results")
#       runFolder 
#   
#   
#   
#   dir.create(paste(folder,"/fig",sep=""))
#   
#   
# }

# function that load runs from a folder and put them in a list for further uses
bunchLoad <- function(folder)
{
  dirContent <- dir(folder)
  listOfSim = list()
  for(i in 1:length(dirContent))
  {
    if (file.exists(paste(folder,"/",dirContent[i],"/run.Rdata",sep = "")))
    {
      sim <- get(load(paste(folder,"/",dirContent[i],"/run.Rdata",sep="")))
      listOfSim[[i]] = sim
    } else if (file.exists(paste(folder,"/",dirContent[i],"/defaultRun.Rdata",sep = "")))
    {
      sim <- get(load(paste(folder,"/",dirContent[i],"/defaultRun.Rdata",sep="")))
      sim = superOpt(sim)
      listOfSim[[i]] = sim
    } 
    listOfSim <- listOfSim[lapply(listOfSim,length)>0]
  }
  return(listOfSim)
}

#function that takes a bunch of run in different folder and plots for each run in their respective folder and also a commom folder of average plots
TotAnalysis <- function(folder)
{
  dirContent <- dir(folder)
  bigSim = list()
  for(i in 1:length(dirContent))
  {
    if (file.exists(paste(folder,"/",dirContent[i],"/run.Rdata",sep = "")))
    {
      sim <- get(load(paste(folder,"/",dirContent[i],"/run.Rdata",sep="")))
    } else {# (file.exists(paste(folder,"/",dirContent[i],"/defaultRun.Rdata",sep = ""))) 
      
      sim <- get(load(paste(folder,"/",dirContent[i],"/defaultRun.Rdata",sep="")))
      sim = superOpt(sim)
    } 
    
    bigSim[[i]] = sim
    
    where = paste(folder,"/",dirContent[i],sep = "")
    
    # p = plotDynamics(sim)
    # mytitle = paste("BiomassSp",".png", sep = "")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    # 
    # p = plotSS(sim, time_range = ((sim@params@species_params$timeMax[1]*0.1-500) : (sim@params@species_params$timeMax[1]*0.1-1)) )
    # mytitle = paste("SizeSpectrum",".png", sep = "")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    # p = PlotNoSp(sim)
    # mytitle = paste("noSp",".png",sep="")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    # p= plotBarSp(sim)
    # mytitle = paste("barPhen",".png",sep="")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    # 
    # p = plotGrowthCurve(sim, generation = 50)
    # mytitle = paste("growthCurve",".png",sep="")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    # 
    # path_to_png = paste(where,"/multi",".png",sep="")
    # png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
    # plotTraitsMulti(sim)
    # dev.off()
    
    # path_to_png = paste(where,"/performance",".png",sep="")
    # png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
    # plotPerformance(sim)
    # dev.off()
    
    # if (rawRun)
    # {
    # p = plotFood(sim)
    # mytitle = paste("food",".png",sep="")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    #
    # p = plotUdead(sim)
    # mytitle = paste("mortality",".png",sep="")
    # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    # }
    rm(sim)
    for (i in 1:20) gc()
    
  }
  
  sim = superStich(bigSim)
  for (i in 1:20) gc()
  
  where = paste(folder,"/results",sep="")
  
  ifelse(!dir.exists(file.path(where)), dir.create(file.path(where)), FALSE) #create the file if it does not exists
  
  save (bigSim,file = paste(where,"/sim.Rdata",sep=""))
  
  # p = plotDynamics(sim)
  # mytitle = paste("BiomassSp",".png", sep = "")
  # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  # 
  # p = plotSS(sim, time_range = ((sim@params@species_params$timeMax[1]*0.1-500) : (sim@params@species_params$timeMax[1]*0.1-1)) )
  # mytitle = paste("SizeSpectrum",".png", sep = "")
  # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  # p = PlotNoSp(bigSim)
  # mytitle = paste("noSp",".png",sep="")
  # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  # p= plotBarSp(bigSim)
  # mytitle = paste("barPhen",".png",sep="")
  # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  # 
  # p = plotGrowthCurve(sim, generation = 50)
  # mytitle = paste("growthCurve",".png",sep="")
  # ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  path_to_png = paste(where,"/multi",".png",sep="")
  png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
  plotTraitsMulti(sim)
  dev.off()
  
  # path_to_png = paste(where,"/performance",".png",sep="")
  # png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
  # plotPerformance(sim)
  # dev.off()
  
}

# load a big run, optimise it and save it in the same folder
dir.Opt <- function(folder)
{
  dirContent <- dir(folder)
  for(i in 1:length(dirContent))
  {
    
    load(paste(folder,"/",dirContent[i],"/defaultRun.Rdata",sep=""))
    sim = superOpt(sim)
    save(sim, file = paste(folder,"/",dirContent[i],"/Run.Rdata",sep="") )
    rm(sim)
  }
}

# take a bunch of sims and get average plots from them
AvgAnalysis <- function(folder,where)
{
  dirContent <- dir(folder)
  bigSim = list()
  for(i in 1:length(dirContent))
  {
    if (file.exists(paste(folder,"/",dirContent[i],"/run.Rdata",sep = "")))
    {
      sim <- get(load(paste(folder,"/",dirContent[i],"/run.Rdata",sep="")))
    } else {
      sim <- get(load(paste(folder,"/",dirContent[i],"/defaultRun.Rdata",sep="")))
      sim = superOpt(sim)
    }
    bigSim[[i]] = sim
    rm(sim)
  }
  sim = superStich(bigSim)
  for (i in 1:20) gc()
  
  ifelse(!dir.exists(file.path(paste(getwd(),"/",where,sep=""))), dir.create(file.path(paste(getwd(),"/",where,sep=""))), FALSE) #create the file if it does not exists
  
  save (bigSim,file = paste(where,"/sim.Rdata",sep=""))
  
  p = plotDynamics(sim)
  mytitle = paste("BiomassSp",".png", sep = "")
  ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  p = plotSS(sim, time_range = ((sim@params@species_params$timeMax[1]*0.1-500) : (sim@params@species_params$timeMax[1]*0.1-1)) )
  mytitle = paste("SizeSpectrum",".png", sep = "")
  ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  p = PlotNoSp(bigSim)
  mytitle = paste("noSp",".png",sep="")
  ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  p= plotBarSp(bigSim)
  mytitle = paste("barPhen",".png",sep="")
  ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  p = plotGrowthCurve(sim, generation = 50)
  mytitle = paste("growthCurve",".png",sep="")
  ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
  
  path_to_png = paste(where,"/multi",".png",sep="")
  png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
  plotTraitsMulti(sim)
  dev.off()
  
  path_to_png = paste(where,"/performance",".png",sep="")
  png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
  plotPerformance(sim)
  dev.off()
  
  
}

# to average multiple sims
superStich <- function(listOfSim)
{
  # params needed
  
  no_w = 100 
  min_w = 0.001 
  max_w = 1e5 * 1.1
  min_w_pp = 1e-10
  w_pp_cutoff = 0.5
  n = 0.75 
  p = 0.75 
  q = 0.8
  r_pp = 4
  kappa = 0.05
  lambda = 2+q-n
  
  EcoName = do.call(c,lapply(listOfSim, function(isim) isim@params@species_params$ecotype))
  
  while (sum(duplicated(EcoName))>0) # while there are some duplicated name in the vector
  {
    for (i in EcoName[duplicated(EcoName)]) # for every duplicated names
    {
      
      x <- as.character(EcoName[which(EcoName == i)[2]]) # take the position in the ecotype vect of the duplicated name and extract its name as character
      first = as.numeric(unlist(strsplit(x, "")))[1] # extract the first number (meaning species number that we want to keep)
      EcoName[which(EcoName == i)[2]] = as.numeric(paste(first,sample(x = seq(1:100000),1),sep="")) #produce a new name but keep the first digit (species identity) 
    }
  }
  
  # Now that I have ecotypes name without duplicate, I can stitch the sims
  
  a = do.call(rbind, lapply(listOfSim, function(isim) isim@params@species_params))
  a$ecotype = EcoName
  
  trait_params <- MizerParams(a, min_w = min_w, max_w=max_w, no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_pp_cutoff, n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda)
  
  
  biom <- do.call(abind, c(lapply(listOfSim, function(isim) isim@n),along = 2))
  dimnames(biom)[[2]] = EcoName
  names(dimnames(biom)) = list("time","species","w")
  
  # need to balance the abundance of biom -> let's mean everything
  biom = biom/length(listOfSim)
  
  sim = listOfSim[[1]]
  sim@params = trait_params
  sim@n = biom
  
  return(sim)
}


# to run in parallel, everything end up in one function for the apply
multiRun <- function(no_sim, no_sp, t_max, mu, no_run,min_w_inf, max_w_inf, effort, print_it = F,
                     fisheries = F,save_it = F, path_to_save = NULL, normalFeeding = F, hartvig = T, Traits = "eta" )
{
  
  if (fisheries)
  {
    #asymptotic size
    w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
    
    #dividing species between the gears (fished and non-fished)
    other_gears <- w_inf >= 10000
    gear_names <- rep("None", no_sp)
    gear_names[other_gears] <- "FishingStuff"
    
    #setting up knife edge
    knife_edges <- w_inf * 0.35 # at maturation size
    # fisheries are primitve now, so I need to set the knife edge a lot above the size of the species that I do not want to fish
    #knife_edges[1:6] <-1e6
    
    output <- myModel(no_sp = no_sp, t_max = t_max, no_run = no_run,
                      mu = mu, max_w_inf = max_w_inf, 
                      effort = effort, knife_edge_size = knife_edges, gear_names = gear_names,
                      save_it = save_it, path_to_save = path_to_save,
                      print_it = print_it, normalFeeding = normalFeeding, hartvig = hartvig, Traits = Traits)
    
    
    
  } else {
    
    output <- myModel(no_sp = no_sp, t_max = t_max, no_run = no_run,
                      mu = mu, max_w_inf = max_w_inf, 
                      effort = 0,
                      save_it = save_it, path_to_save = path_to_save,
                      print_it = print_it, normalFeeding = normalFeeding, hartvig = hartvig, Traits = Traits)
  }
  
  
}


# to lighten the size of the sims
superOpt <- function(object, opt = 10) # optimise the run to make them lighter
{
  object@n = object@n[seq(1,dim(object@n)[1],opt),,]
  object@n_pp = object@n_pp[seq(1,dim(object@n_pp)[1],opt),]
  object@effort = as.matrix(object@effort[seq(1,dim(object@effort)[1],opt),])
  return(object)
}

stitch <- function(sim) # this function takes all the sim output from model and paste together the data to make one n file (for plotting and stuff)
{
  
  endList <- length(sim) # shortcut to have ref of the last simulation which has the right dim, names, ...
  Dtime = dim(sim[[endList]]@n)[1] # need these dimentions for making new arrays but separatly as I'm going to work on sp dim
  Dsp = dim(sim[[endList]]@n)[2]
  Dw = dim(sim[[endList]]@n)[3]
  
  # if no mutants
  
  if (endList == 1) return(list(sim[[1]]@n, sim[[1]]@n_pp)) # if there is only one object in the list, it is the mizer object without mutants
  
  # from there I can add empty columns to every other array and add them 
  for (i in 1:(endList-1))# check every array except last one because no need
  {
    sim[[i]]@n <-abind(sim[[i]]@n, array(dim = c(Dtime, Dsp - dim(sim[[i]]@n)[2],Dw)), along = 2) # I add to each sim the missing columns pre-mutations
    sim[[i]]@n_pp <-abind(sim[[i]]@n_pp, array(dim = c(Dtime, Dsp - dim(sim[[i]]@n)[2])), along = 2)
  }
  # Now every array in the output has the dimemsions of the final array, but with NAs everywhere
  
  # Now to add them up
  # getting rid of the NAs
  for (j in 1:endList) sim[[j]]@n[is.na( sim[[j]]@n)] <- 0
  for (j in 1:endList) sim[[j]]@n_pp[is.na( sim[[j]]@n_pp)] <- 0
  
  # addition
  biomass = 0
  for (j in 1:endList) biomass = biomass + sim[[j]]@n
  colnames(biomass) <- colnames(sim[[endList]]@n)
  names(dimnames(biomass)) <- list("time","sp","w")
  
  biomass_pp = 0
  for (j in 1:endList) biomass_pp = biomass_pp + sim[[j]]@n_pp
  colnames(biomass_pp) <- colnames(sim[[endList]]@n_pp)
  names(dimnames(biomass_pp)) <- list("time","sp")
  
  return(list(biomass,biomass_pp))
}

# to plot only one species family Not sure its working anymore, but function included in plots directly
lineage = function(sim,thread)
{
  noSp = length(sim@params@species_params$species) # number of species after the simulation
  family = sim@params@species_params[sim@params@species_params$species == thread, ] # create a smaller table with only one lineage
  idxFamily = rownames(family)
  tree = sim@n[,idxFamily,] # get their biomass from the big table
  sim@n = tree
  plotDynamics(sim)
}

# to process the mizer output and draw plots, works with the old sim output
processing <- function(result, 
                       plot = FALSE, 
                       where = NULL,
                       who = "defaultRun",
                       dt = 0.1, # I need dt at some point and I don't know how to get it there
                       optimisation = F,
                       save_it = T
)
{
  # setting up workspace things
  if (is.null(where)) where = paste(getwd(),"/temporary",sep="")
  ifelse(!dir.exists(file.path(where)), dir.create(file.path(where),recursive = T), FALSE) #create the file if it does not exists
  
  ## processing data
  # a result will be a list of simulations (number of runs) and the ID card of the ecosystem, 
  # if not it means that I am inputing directly an already processed mizer object, ready for the plots
  if (class(result) == "list")
  {
    if (class(result[[length(result)]]) == "list") #if this is a list it means I got umbrella
    {
      print("umbrella time")
      result[[length(result)]] = NULL # delete last half sim
      param = NULL
      for (i in 1:length(result)) param = rbind(param,result[[i]]@params@species_params) # create the dataframe for species
      param <- param[order(param$ecotype, param$extinct, decreasing=TRUE),] 
      param <- param[!duplicated(param$ecotype),]
      SummaryParams = param[order(param$pop,param$ecotype),]
      FinalParam <- MizerParams(SummaryParams, min_w =0.001, max_w=10000 * 1.1, no_w = 100, min_w_pp = 1e-10, w_pp_cutoff = 0.5, n = 0.75, p=0.75, q=0.8, r_pp=4, kappa=1, lambda = 2.05) #create the mizer param from the dataframe
      temp=list(result,FinalParam) # put it in the right disposition
      result = temp
      # something is wrong with it, don't have the time or interest to debug, just return nothing, it's not like it's useful
      return()
    }
    
    gc()
    sim = result[[1]]
    SummaryParams = result[[2]]
    #endlist = length(sim) # shortcut to have the number of sim
    rm(result)
    # lol = list() # here I get rid of the empty spots of my list. It will happen if I started a simulation from a previous one as initial conditions 
    # for (i in 1:endlist) if (!is.null(sim[[i]])) lol = c(lol, sim[[i]])
    # sim = lol
    # rm(lol)
    sim <- sim[lapply(sim,length)>0]
    endlist = length(sim) # updating the number
    template = sim[[endlist]] # to keep a template of mizer object somewhere
    gc()
    # if (optimisation)
    # {
    #   keep = seq(1,SummaryParams@species_params$timeMax[1]*0.1) #times to keep (too long but doesnt matter)
    #   for (i in 1:endlist)
    #   {
    #     timeKeep =  dimnames(sim[[i]]@n)$time %in% keep
    #    sim[[i]]@n =  sim[[i]]@n[which(timeKeep),,]
    #    sim[[i]]@n_pp =  sim[[i]]@n_pp[which(timeKeep),]
    #    sim[[i]]@effort =  as.matrix(sim[[i]]@effort[which(timeKeep[-1]),]) # effort is being a pain again as it's missing the time 0
    #   }
    # }   
    
    
    if(optimisation)
    {
      for (i in 1:endlist)
      {
        sim[[i]]@n = sim[[i]]@n[seq(1,dim(sim[[i]]@n)[1]-1,2),,]
        sim[[i]]@n_pp = sim[[i]]@n_pp[seq(1,dim(sim[[i]]@n)[1]-1,2),]
        #sim[[i]]@effort = sim[[i]]@effort[seq(1,dim(sim[[i]]@n)[1],2),]
      }
    }
    
    # stitiching the sims together
    Dtime = dim(sim[[1]]@n)[1] # need these dimensions for making new arrays but separatly as I'm going to work on sp dim
    Dsp = length(SummaryParams@species_params$ecotype)
    Dw = dim(sim[[1]]@n)[3]
    # put all the sim at the same dimension
    for (i in 1:endlist) 
    {
      # print(i)
      for (j in 1:Dsp)
      {
        gc()
        if (is.na(dimnames(sim[[i]]@n)$sp[j])) # if I need to add species at the end of the array
        {
          sim[[i]]@n <- abind(sim[[i]]@n,array(dim = c(Dtime,Dw)),along = 2)
          names(dimnames(sim[[i]]@n)) <- list("time","sp","w") # the abind make me loose the info
        }
        else if (dimnames(sim[[i]]@n)$sp[j] != SummaryParams@species_params$ecotype[j]) # if I need to add in the middle of the array
        {
          rmbr = NULL # when I abind only one column of the array (happens when there is an extinction right before the last ecotype), I lose its name, leading to bugs
          
          if (j == dim(sim[[i]]@n)[2]) rmbr = dimnames(sim[[i]]@n)$sp[j] # with this I know when it happens
          
          sim[[i]]@n <- abind(sim[[i]]@n[,1:j-1,],array(dim = c(Dtime,Dw)),sim[[i]]@n[,j:dim(sim[[i]]@n)[2],], along = 2) # if sp 1 goes extinct I will have a bug, or not it seems
          names(dimnames(sim[[i]]@n)) <- list("time","sp","w") # the abind makes me loose the info
          
          if (is.null(rmbr)==F) dimnames(sim[[i]]@n)$sp[j+1] = rmbr # and I add back the name
          
        }
      }
      dimnames(sim[[i]]@n)$sp = SummaryParams@species_params$ecotype # put back the right names
      if (optimisation == F) sim[[i]]@n = sim[[i]]@n[-1,,] # in my tentative to include fisheries, I am deleting the first time step (which is 0 now) to have the same length as the effort matrix. I'm running out of time, Ill fix that later I hope
    }
    # getting rid of the NAs
    gc()
    
    for (j in 1:endlist) sim[[j]]@n[is.na( sim[[j]]@n)] <- 0
    gc()
    biom <- do.call(abind, c(lapply(sim, function(isim) isim@n),along = 1))
    
    gc()
    
    
    
    # now add them in a row
    # biom = sim[[1]]@n
    # for (i in 2:endlist) biom = abind(biom,sim[[i]]@n,along = 1)
    
    
    
    names(dimnames(biom)) = list("time","species","size")
    
    if (optimisation) dimnames(biom)$time = seq(0, (SummaryParams@species_params$timeMax[1]-1) *dt,0.2) # skip every other time step and reduce to years (not year/dt)
    if (optimisation == F) dimnames(biom)$time = seq(1, SummaryParams@species_params$timeMax[1]) * dt # i don't know whats wrong with the else
    gc()
    
    # I have to do the phyto aussi
    
    
    # phyto = template@n_pp
    # phyto = phyto[seq(1,dim(phyto)[1]-1,2),]
    # phyta = phyto
    # for (i in 2:endlist) phyto = abind(phyto,phyta,along = 1)
    
    if (optimisation)
    {
      phyto = sim[[1]]@n_pp
      for (i in 2:endlist) phyto = abind(phyto,sim[[i]]@n_pp,along = 1)
    }
    if (optimisation == F)
    {
      phyto = sim[[1]]@n_pp[-1,] #get rid of time 0
      for (i in 2:endlist) phyto = abind(phyto,sim[[i]]@n_pp[-1,],along = 1)
    }
    cat(sprintf("Biomass handeled successfully, now starting the effort.\n"))
    
    # taking care of the effort
    # for now it is constant in time so we just need to multiply effort by no_run
    effort = sim[[1]]@effort
    for(i in 2:endlist) effort = rbind(effort,sim[[i]]@effort)
    names(dimnames(effort)) = list("time","effort")
    if (optimisation) dimnames(effort)$time = seq(1, SummaryParams@species_params$timeMax[1] * dt)
    if (optimisation == F) dimnames(effort)$time = seq(1, SummaryParams@species_params$timeMax[1]) * dt # i don't know whats wrong with the else
    
    # I need to assemble the mizer object now
    sim = template
    sim@params=SummaryParams
    sim@n = biom
    sim@effort = effort
    sim@n_pp = phyto
    
    rm(list = "template","biom","phyto","effort")
    gc()
    
    #saving
    path_to_save <- paste(where,"/",who,".Rdata", sep="")
    if (save_it) save(sim,file = path_to_save)
  }
  
  else
  {
    sim = result
    rm(result)
    #saving
    path_to_save <- paste(where,"/",who,".Rdata", sep="")
    if (save_it) save(sim,file = path_to_save)
  }
  
  
  if (plot == TRUE)
  {
    cat(sprintf("Data processed, starting the plots"))
    
    SumPar = sim@params@species_params # little shortcut
    
    
    # biomass variation
    
    p = plotDynamics(sim, species = T)
    
    mytitle = paste("BiomassSp",".png", sep = "")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    # size spectrum
    
    p = plotSS(sim, time_range = ((sim@params@species_params$timeMax[1]*dt-500) : (sim@params@species_params$timeMax[1]*dt-1)) )
    mytitle = paste("SizeSpectrum",".png", sep = "")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    # traits plots
    p = plotTraits(sim)
    
    for (i in 1:length(p))
    {
      for (j in 1:length(p[[1]]))
      {
        y = as.character(i)
        switch(y,
               "1" = {mytitle = paste("Maturation size of specie ",j,".png", sep = "")},
               "2" = {mytitle = paste("PPMR of specie ",j,".png", sep = "")},
               "3" = {mytitle = paste("Diet breadth of specie ",j,".png", sep = "")},
               {})
        ggsave(mytitle, plot = p[[i]][[j]], path = where, width = 20, height = 20, units = "cm")
      }
    }
    
    # # beta/sigma ratio---------------
    # for (i in SpIdx)
    # {
    #   # empty matrix of ecotype of species i by time
    #   A = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
    #   # fill the matrix with ones when the ecotype exists
    #   for (x in 1:nrow(A)) # I'm sure I can do an apply but don't know how
    #   {
    #     for (j in 1:ncol(A))
    #     {
    #       if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) A[x,j] = 1
    #     }}
    #   # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
    #   # change the ones by the trait value of the ecotype
    #   BetaA = A * TT[TT$Lineage == i,]$PPMR
    #   SigmaA = A * TT[TT$Lineage == i,]$Diet_breadth
    #   # calculate mean trait value at each time step
    #   no_trait = apply(A,2,sum) # this vector is the number of traits present at each time step
    #   BetaSum = apply(BetaA,2,sum) # this vector is the sum of the traits value at each time step
    #   SigmaSum = apply(SigmaA,2,sum) # this vector is the sum of the traits value at each time step
    #   BetaMean = BetaSum/no_trait # this is the mean trait value at each time step
    #   SigmaMean = SigmaSum/no_trait # this is the mean trait value at each time step
    #   
    #   # Matrix with all traits combination and at what time they go extinct
    #   TTi = TT[TT$Lineage == i,]
    #   comb=data.frame(TTi$Ecotype,TTi$Apparition,TTi$Extinction,TTi$PPMR,TTi$Diet_breadth)
    #   
    #   # plot of extinction of combinations
    #   title = paste("Combination of PPMR and diet breath value of species ",i, sep = "")
    #   print(
    #     ggplot(comb) +
    #       geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Extinction)) +
    #       scale_x_continuous(name = "PPMR") +
    #       scale_y_continuous(name = "Diet breath") +
    #       scale_color_continuous(name = "Extinction time in year / dt") +
    #       ggtitle(title)
    #   )
    #   name = paste("Extinction BetaSigma of species",i, sep = "")
    #   setwd(where)
    #   mytitle = paste(name,".png", sep = "")
    #   dev.print(png, mytitle, width = res, height = 2/3* res)
    #   
    #   #plot of apparition of combinations
    #   print(
    #     ggplot(comb) +
    #       geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Apparition)) +
    #       scale_x_continuous(name = "PPMR") +
    #       scale_y_continuous(name = "Diet breath") +
    #       scale_color_continuous(name = "Apparition time in year / dt") +
    #       ggtitle(title)
    #   )
    #   name = paste("Apparition BetaSigma of species",i, sep = "")
    #   setwd(where)
    #   mytitle = paste(name,".png", sep = "")
    #   dev.print(png, mytitle, width = res, height = 2/3* res)
    #   
    #   
    #   # Mean combination values throughout sim
    #   stat = data.frame(BetaMean,SigmaMean) 
    #   #get rid of duplicates
    #   stat = stat[!duplicated(stat),]
    #   # need to work on the time because of fucked up legend
    #   stat = cbind(stat,rownames(stat))
    #   dimnames(stat)[[2]] = list("BetaMean", "SigmaMean", "Time")
    #   stat$Time <- as.numeric(substr(stat$Time,1,5)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5
    #   
    #   title = paste("Mean trait's combination of PPMR and diet breath value of species ",i, sep = "")
    #   print(
    #     ggplot(stat) +
    #       geom_point(data = stat, aes(x=BetaMean,y=SigmaMean, color = Time)) +
    #       scale_x_continuous(name = "PPMR") +
    #       scale_y_continuous(name = "Diet breath") +
    #       scale_color_continuous(name = "Time in year / dt") +
    #       ggtitle(title)
    #   )
    #   name = paste("Mean BetaSigma of species",i, sep = "")
    #   setwd(where)
    #   mytitle = paste(name,".png", sep = "")
    #   dev.print(png, mytitle, width = res, height = 2/3* res)
    # }
    # 
    # biology plots -------------
    
    # when = seq(500,sim@params@species_params$timeMax[1]*0.1,1000)
    # for(i in when)
    # {
    #   p = plotFood(sim, time_range = i)
    #   mytitle = paste("feeding t= ",i,".png",sep="")
    #   ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    #   
    #   p = plotUdead(sim, time_range = i)
    #   mytitle = paste("mortality t= ",i,".png",sep="")
    #   ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    #   
    #   p = plotGrowth(sim, time_range = i)
    #   mytitle = paste("growth t= ",i,".png",sep="")
    #   ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    # }
    
    p = plotFood(sim)
    mytitle = paste("food",".png",sep="")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    p = plotUdead(sim)
    mytitle = paste("mortality",".png",sep="")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    p = PlotNoSp(sim)
    mytitle = paste("noSp",".png",sep="")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm")
    
    p = plotGrowthSpeed(sim)
    mytitle = paste("growthSpeed",".png",sep="")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm") 
    
    p = plotGrowthCurve(sim, generation = 50)
    mytitle = paste("growthCurve",".png",sep="")
    ggsave(mytitle, plot = p, path = where, width = 20, height = 20, units = "cm") 
    
    path_to_png = paste(where,"/multi",".png",sep="")
    png(filename=path_to_png,width = 20, height = 20, units = "cm",res = 400)
    plotTraitsMulti(sim)
    dev.off()
    
  }
  
  
  
  
  
  
  
  
  
  return(sim)
}


# function that draw the bloodline of an ecotype by plotting its trait combination (could do 3D) Need work
ancestor <- function ()
{
  SumPar = sim@params@species_params
  TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),SumPar$pop,SumPar$extinct,SumPar$w_mat,SumPar$beta,SumPar$sigma) # weird things happen without the as.numeric
  colnames(TT) = c("Lineage","Ecotype","Apparition","Extinction","Maturation_size","PPMR","Diet_breadth")
  rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  TT = as.data.frame(TT)
  
  for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = SumPar$timeMax[1]
  res = 1000
  
  i = 1 # sp
  # empty matrix of ecotype of species i by time
  A = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
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
  
  
  tree = NULL
  for (j in 2:dim(comb)[1]) # for each fucking ecotype (except first)
  {
    x = comb[j, "TTi.Ecotype"]
    bloodLine = x
    while (x != i) # I determine from whom they descend
    {
      x = x %/% 10
      bloodLine = c(bloodLine, x)
    }
    ancestor = NULL
    for (l in 1:dim(comb)[1]) # and I put all these people in one matrix
    {
      for (k in 1:length(bloodLine))
      {
        if (comb$TTi.Ecotype[l] == bloodLine[k])
          ancestor = rbind(ancestor, comb[l, ])
      }
    }
    ancestor = cbind(ancestor, j)
    tree = rbind(tree,ancestor)
  }
  
  print(
    ggplot(comb) +
      geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Apparition, size = 1) ) +
      geom_path(data = tree, aes(x = TTi.PPMR, y = TTi.Diet_breadth, group = j), arrow = arrow(angle = 15, type = "closed", length = unit(0.1,"inches"))) +
      scale_x_continuous(name = "PPMR") +
      scale_y_continuous(name = "Diet breath") +
      scale_color_continuous(name = "Apparition time in year / dt") +
      ggtitle(title)
  )
  
  
  
  # plot of extinction of combinations
  title = paste("Combination of PPMR and diet breath value of species ",i, sep = "")
  print(
    ggplot(comb) +
      geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Extinction)) +
      scale_x_continuous(name = "PPMR") +
      scale_y_continuous(name = "Diet breath") +
      scale_color_continuous(name = "Extinction time in year / dt") +
      ggtitle(title)
  )
  name = paste("Extinction BetaSigma of species",i, sep = "")
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste(name,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 2/3* res)
  
  #plot of apparition of combinations
  
  print(
    ggplot(comb) +
      geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Apparition)) +
      geom_path(data = tree, aes(x = TTi.PPMR, y = TTi.Diet_breadth, group = j)) +
      scale_x_continuous(name = "PPMR") +
      scale_y_continuous(name = "Diet breath") +
      scale_color_continuous(name = "Apparition time in year / dt") +
      ggtitle(title)
  )
  name = paste("Apparition BetaSigma of species",i, sep = "")
  setwd(paste(dir,subdir, sep = ""))
  mytitle = paste(name,".png", sep = "")
  dev.print(png, mytitle, width = res, height = 2/3* res)
  
  
  # Mean combination values throughout sim
  stat = data.frame(BetaMean,SigmaMean) 
  #get rid of duplicates
  stat = stat[!duplicated(stat),]
  # need to work on the time because of fucked up legend
  stat = cbind(stat,rownames(stat))
  dimnames(stat)[[2]] = list("BetaMean", "SigmaMean", "Time")
  stat$Time <- as.numeric(substr(stat$Time,1,5)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5
  
  title = paste("Mean trait's combination of PPMR and diet breath value of species ",i, sep = "")
  print(
    ggplot(stat) +
      geom_point(data = stat, aes(x=BetaMean,y=SigmaMean, color = Time)) +
      scale_x_continuous(name = "PPMR") +
      scale_y_continuous(name = "Diet breath") +
      scale_color_continuous(name = "Time in year / dt") +
      ggtitle(title)
  )
  
  # all present traits at last time step
  # last = data.frame(BetaA[,t_max*no_run],SigmaA[,t_max*no_run])
  # colnames(last) = c("beta","sigma")
  # row_sub = apply(last, 1, function(row) all(row !=0 ))
  # last = last[row_sub,]
  # 
  # 
  # print(
  #   ggplot(last) +
  #     geom_point(aes(x=beta,y=sigma)) +
  #     ggtitle("Last time step")
  # )
  
  
  
}


#function that reduce the number of data if too heavy for plotting
optimisation <- function(sim, keep = 2)
{
  sim@n = sim@n[seq(1,dim(sim@n)[1],keep),,]
  sim@n_pp = sim@n_pp[seq(1,dim(sim@n_pp)[1],keep),]
  sim@effort = sim@effort[seq(1,dim(sim@effort)[1],keep),]
  
  return(sim)
}


# # old trait tree function included in processing but keeping it just in case (has not weighted plots) -----------------
# #TT = cbind(sim@params@species_params$species,sim@params@species_params$pop,sim@params@species_params$extinct,sim@params@species_params$w_mat,sim@params@species_params$beta,sim@params@species_params$sigma)
# TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),SumPar$pop,SumPar$extinct,SumPar$w_mat,SumPar$beta,SumPar$sigma) # weird things happen without the as.numeric
# colnames(TT) = c("Lineage","Ecotype","Apparition","Extinction","Maturation_size","PPMR","Diet_breadth")
# rownames(TT) = rownames(SumPar)
# TT = TT[order(TT[,1],decreasing=FALSE),]
# TT = as.data.frame(TT)
# no_run = 25 # need to specify that if I did more than one succession
# for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = SumPar$timeMax[1]
# res = 1000
# #plot of traits maturation size
# for (i in 1:no_sp)
# {
#   # empty matrix of ecotype of species i by time
#   a = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
#   # fill the matrix with ones when the ecotype exists
#   for (x in 1:nrow(a)) # I'm sure I can do an apply but don't know how
#   {
#     for (j in 1:ncol(a))
#     {
#       if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) a[x,j] = 1
#     }}
#   # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
#   # calculate the mean trait value + standard deviation
#   # change the ones by the trait value of the ecotype
#   b = a * TT[TT$Lineage == i,]$Maturation_size
#   # calculate mean trait value at each time step
#   c = apply(a,2,sum) # this vector is the number of traits present at each time step
#   d = apply(b,2,sum) # this vector is the sum of the traits value at each time step
#   e = d/c # this is the mean trait value at each time step
#   f = apply((sweep(b,2,e,"-")*a)^2,2,sum)/c # this is the variance at each time step
#   g = sqrt(f) # this is the standard population deviation at each time step
#   
#   stat = data.frame(e,g,e-g,e+g)
#   colnames(stat) = list("mean","variance","low","up")
#   
#   X = TT$Maturation_size[TT$Ecotype==i] # value around which to do the breaks
#   
#   family = TT[TT$Lineage == i,3:5]
#   Fam = melt(family,"Maturation_size") 
#   name = paste("Maturation size of species ",i, sep = "")
#   p = ggplot() +
#     geom_point(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) +
#     geom_line(data = Fam, aes(x = value, y = Maturation_size, group = Maturation_size)) + 
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
#     scale_x_continuous(name = "Time") +
#     scale_y_continuous(name = "Trait value") +
#     scale_colour_discrete(labels = c("mean","standard deviation"))+
#     theme(legend.title=element_blank(),panel.background = element_blank(), legend.position=c(0.9,1), legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
#     guides(color=guide_legend(override.aes=list(fill=NA)))+
#     ggtitle(name)
#   print(p)
#   setwd(paste(dir,subdir, sep = ""))
#   mytitle = paste(name,".png", sep = "")
#   dev.print(png, mytitle, width = res, height = 2/3* res) 
# }
# 
# 
# #PPMR
# for (i in 1:no_sp)
# {
#   # empty matrix of ecotype of species i by time
#   a = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
#   # fill the matrix with ones when the ecotype exists
#   for (x in 1:nrow(a)) # I'm sure I can do an apply but don't know how
#   {
#     for (j in 1:ncol(a))
#     {
#       if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) a[x,j] = 1
#     }}
#   # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
#   # change the ones by the trait value of the ecotype
#   b = a * TT[TT$Lineage == i,]$PPMR
#   # calculate mean trait value at each time step
#   c = apply(a,2,sum) # this vector is the number of traits present at each time step
#   d = apply(b,2,sum) # this vector is the sum of the traits value at each time step
#   e = d/c # this is the mean trait value at each time step
#   f = apply((sweep(b,2,e,"-")*a)^2,2,sum)/c # this is the variance at each time step
#   g = sqrt(f) # this is the standard population deviation at each time step
#   
#   stat = data.frame(e,g,e-g,e+g)
#   colnames(stat) = list("mean","variance","low","up")
#   
#   family = TT[TT$Lineage == i,c(3,4,6)]
#   Fam = melt(family,"PPMR") 
#   name = paste("PPMR size of species ",i, sep = "")
#   p = ggplot() +
#     geom_point(data = Fam, aes(x = value, y = PPMR, group = PPMR)) +
#     geom_line(data = Fam, aes(x = value, y = PPMR, group = PPMR)) + 
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
#     scale_x_continuous(name = "Time") +
#     scale_y_continuous(name = "Trait value") +
#     scale_colour_discrete(labels = c("mean","standard deviation"))+
#     theme(legend.title=element_blank(),panel.background = element_blank(), legend.position=c(0.9,1), legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
#     guides(color=guide_legend(override.aes=list(fill=NA)))+
#     ggtitle(name)
#   print(p)
#   setwd(paste(dir,subdir, sep = ""))
#   mytitle = paste(name,".png", sep = "")
#   dev.print(png, mytitle, width = res , height = 2/3* res) 
# }
# 
# # Diet breath
# for (i in 1:no_sp)
# {
#   # empty matrix of ecotype of species i by time
#   a = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
#   # fill the matrix with ones when the ecotype exists
#   for (x in 1:nrow(a)) # I'm sure I can do an apply but don't know how
#   {
#     for (j in 1:ncol(a))
#     {
#       if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) a[x,j] = 1
#     }}
#   # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
#   # change the ones by the trait value of the ecotype
#   b = a * TT[TT$Lineage == i,]$Diet_breadth
#   # calculate mean trait value at each time step
#   c = apply(a,2,sum) # this vector is the number of traits present at each time step
#   d = apply(b,2,sum) # this vector is the sum of the traits value at each time step
#   e = d/c # this is the mean trait value at each time step
#   f = apply((sweep(b,2,e,"-")*a)^2,2,sum)/c # this is the variance at each time step
#   g = sqrt(f) # this is the standard population deviation at each time step
#   
#   stat = data.frame(e,g,e-g,e+g)
#   colnames(stat) = list("mean","variance","low","up")
#   
#   family = TT[TT$Lineage == i,c(3,4,7)]
#   Fam = melt(family,"Diet_breadth") 
#   name = paste("Diet breadth size of species ",i, sep = "")
#   p = ggplot() +
#     geom_point(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) +
#     geom_line(data = Fam, aes(x = value, y = Diet_breadth, group = Diet_breadth)) + 
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = mean, color =" red")) +
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = low, color ="blue")) +
#     geom_smooth(data = stat, aes(x = as.numeric(row.names(stat)), y = up, color ="blue")) +
#     scale_x_continuous(name = "Time") +
#     scale_y_continuous(name = "Trait value") +
#     scale_colour_discrete(labels = c("mean","standard deviation"))+
#     theme(legend.title=element_blank(),panel.background = element_blank(), legend.position=c(0.9,1), legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
#     guides(color=guide_legend(override.aes=list(fill=NA)))+
#     ggtitle(name)
#   print(p)
#   setwd(paste(dir,subdir, sep = ""))
#   mytitle = paste(name,".png", sep = "")
#   dev.print(png, mytitle, width = res, height = 2/3* res)
# }
# 
# # beta/sigma ratio
# for (i in 1:no_sp)
# {
#   # empty matrix of ecotype of species i by time
#   A = matrix(0, ncol = SumPar$timeMax[1], nrow = dim(TT[TT$Lineage == i,])[1], dimnames = list(as.numeric(TT$Ecotype[TT$Lineage == i]), c(1:(SumPar$timeMax[1]))))
#   # fill the matrix with ones when the ecotype exists
#   for (x in 1:nrow(A)) # I'm sure I can do an apply but don't know how
#   {
#     for (j in 1:ncol(A))
#     {
#       if (TT$Apparition[x] <= j & TT$Extinction[x] >= j) A[x,j] = 1
#     }}
#   # a is a matrix of 0 and 1 showing if the ecotype is present or not at time t
#   # change the ones by the trait value of the ecotype
#   BetaA = A * TT[TT$Lineage == i,]$PPMR
#   SigmaA = A * TT[TT$Lineage == i,]$Diet_breadth
#   # calculate mean trait value at each time step
#   no_trait = apply(A,2,sum) # this vector is the number of traits present at each time step
#   BetaSum = apply(BetaA,2,sum) # this vector is the sum of the traits value at each time step
#   SigmaSum = apply(SigmaA,2,sum) # this vector is the sum of the traits value at each time step
#   BetaMean = BetaSum/no_trait # this is the mean trait value at each time step
#   SigmaMean = SigmaSum/no_trait # this is the mean trait value at each time step
#   
#   # Matrix with all traits combination and at what time they go extinct
#   TTi = TT[TT$Lineage == i,]
#   comb=data.frame(TTi$Ecotype,TTi$Apparition,TTi$Extinction,TTi$PPMR,TTi$Diet_breadth)
#   
#   # plot of extinction of combinations
#   title = paste("Combination of PPMR and diet breath value of species ",i, sep = "")
#   print(
#     ggplot(comb) +
#       geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Extinction)) +
#       scale_x_continuous(name = "PPMR") +
#       scale_y_continuous(name = "Diet breath") +
#       scale_color_continuous(name = "Extinction time in year / dt") +
#       ggtitle(title)
#   )
#   name = paste("Extinction BetaSigma of species",i, sep = "")
#   setwd(paste(dir,subdir, sep = ""))
#   mytitle = paste(name,".png", sep = "")
#   dev.print(png, mytitle, width = res, height = 2/3* res)
#   
#   #plot of apparition of combinations
#   print(
#     ggplot(comb) +
#       geom_point(aes(x=TTi.PPMR,y=TTi.Diet_breadth, color = TTi.Apparition)) +
#       scale_x_continuous(name = "PPMR") +
#       scale_y_continuous(name = "Diet breath") +
#       scale_color_continuous(name = "Apparition time in year / dt") +
#       ggtitle(title)
#   )
#   name = paste("Apparition BetaSigma of species",i, sep = "")
#   setwd(paste(dir,subdir, sep = ""))
#   mytitle = paste(name,".png", sep = "")
#   dev.print(png, mytitle, width = res, height = 2/3* res)
#   
#   
#   # Mean combination values throughout sim
#   stat = data.frame(BetaMean,SigmaMean) 
#   #get rid of duplicates
#   stat = stat[!duplicated(stat),]
#   # need to work on the time because of fucked up legend
#   stat = cbind(stat,rownames(stat))
#   dimnames(stat)[[2]] = list("BetaMean", "SigmaMean", "Time")
#   stat$Time <- as.numeric(substr(stat$Time,1,5)) # read the time and delete all conditions on it (like factor) # the 5 means handling number up to 10^5
#   
#   title = paste("Mean trait's combination of PPMR and diet breath value of species ",i, sep = "")
#   print(
#     ggplot(stat) +
#       geom_point(data = stat, aes(x=BetaMean,y=SigmaMean, color = Time)) +
#       scale_x_continuous(name = "PPMR") +
#       scale_y_continuous(name = "Diet breath") +
#       scale_color_continuous(name = "Time in year / dt") +
#       ggtitle(title)
#   )
#   
#   # all present traits at last time step
#   # last = data.frame(BetaA[,t_max*no_run],SigmaA[,t_max*no_run])
#   # colnames(last) = c("beta","sigma")
#   # row_sub = apply(last, 1, function(row) all(row !=0 ))
#   # last = last[row_sub,]
#   # 
#   # 
#   # print(
#   #   ggplot(last) +
#   #     geom_point(aes(x=beta,y=sigma)) +
#   #     ggtitle("Last time step")
#   # )
#   
# }


# plotting a normal distribution ------------------
# x<-seq(-4,4,length=200)
# y<-dnorm(x,mean=0, sd=1)
# 
# df = data.frame(x,y)
# 
# ggplot(df)+
#   geom_line(aes(x=x, y=y), size = 2) +
#   theme(panel.background = element_blank()) +
#   geom_vline(xintercept = 0, linetype = "dashed", size = 2) +
#   scale_x_continuous(name = element_blank(), breaks = NULL) +
#   scale_y_continuous( name = element_blank(), breaks = NULL) +
#   theme(panel.background = element_rect(linetype = "solid"))

# one line code to get the size of the objects in the workspace
#for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }