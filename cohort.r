# cohort plot -------------
# I want the spawn through life -> last value of total cohortR
plotCohort <- function(object, dt = 0.1, t_steps = 5, iSpecies = 1, effort = 0, cohortSpan = seq(max(dim(object@n)[1])-30,max(dim(object@n)[1]),2),
                       print_it = T, returnData = F, save_it = F, nameSave = paste("CohortSpecies",iSpecies,".png",sep=""))
{
  # setting up some parameters
sex_ratio = 0.5
T = t_steps/dt; # number of time steps you want to follow cohort for
no_sp <- dim(object@params@species_params)[1]
#PhenIdx <- which(object@params@species_params$species == iSpecies) # select who to track, not the name of the species but their position in the dataframe
PhenIdx <- seq(1,no_sp)
PhenName<- object@params@species_params$ecotype[PhenIdx] # this is their name
no_Phen = length(PhenIdx)
fitness <- array(0,c(no_Phen, length(cohortSpan)), dimnames = list(PhenIdx,cohortSpan)) #collect the total spawn per time (start of cohort) per species
names(dimnames(fitness)) <- list("species","cohort")

for (t_start in cohortSpan)
{
  cat(sprintf("Cohort number %g\n",t_start))
  
  # Initialising matrixes
  cohortW = array(0, c(no_sp, T+1)); # row vector for following cohort weight
  cohortS = array(0, c(no_sp, T+1)); # vector for cohort survival
  cohortR = array(0, c(no_sp, T+1)); # vector for cohort spawning
  cohortR_sol = array(0, c(no_sp, T+1)); # vector for cohort spawn at size
  cohortW[,1] = object@params@w[1]; # log weight initially (newborn)
  cohortS[,1] = object@n[t_start,,1]; # initial population in spectrum

    for (t in seq(1,T)){ # within time period you're interested in
      # vector of the previous size bin for every phenotypes
      cohortWprev = unlist(lapply(lapply(cohortW[PhenIdx,t], FUN = function(x) x-object@params@w), FUN = function(x) max(which(x>= 0)))) # yolo
      # growth matrix
      growth = getEGrowth(object@params,n = object@n[t_start+t-1,,],n_pp = object@n_pp[t_start+t-1,])
      # update the new size bin with the growth
      cohortW[PhenIdx,t+1] = cohortW[PhenIdx,t]+dt*diag(growth[PhenIdx,cohortWprev])
      # mortality matrix
      z = getZ(object = object@params, n = object@n[t_start+t-1,,],n_pp = object@n_pp[t_start+t-1,], effort = effort)
      # update the amount surviving the time-step
      cohortS[PhenIdx,t+1] = cohortS[PhenIdx,t]*exp(-dt*diag(z[PhenIdx,cohortWprev]))
      # need to take global n to have the right amount of resources available but need to take the right fraction at the end for the fitness, as not all the individuals reproducing are part of the cohort.
      # need to prepare n for no NAN, I just want the n of the specific cohort so I extract the right fraction
      n = object@n[t_start+t-1,,]#*cohortS[,t]/cohortS[,1] # this is my biomass x fraction
      #n[!is.finite(n)] <- 0
      # get the rdi manually to have it spread over size bins
      e_spawning <- getESpawning(object = object@params, n = n,n_pp = object@n_pp[t_start+t-1,])
      #print(e_spawning[PhenIdx,70:100])
      e_spawning_pop <- apply((e_spawning*n),1,"*",object@params@dw)
      #print(e_spawning_pop[70:100,PhenIdx])
      rdi <- sex_ratio*(e_spawning_pop * object@params@species_params$erepro)/object@params@w[object@params@species_params$w_min_idx] # global rdi
      #rdi <- rdi*cohortS[,t]/cohortS[,1]
      #print(rdi[70:100,PhenIdx])
      # get the total abundance in each species/size combination from PhenIdx and cohortWprev
      cohortN <- NULL
      for (i in 1:dim(n)[1]) cohortN <- c(cohortN,n[PhenIdx[i],cohortWprev[i]])
      # get the proportion of abundance from the followed cohort
      cohortF <- cohortS[,t]/cohortN
      cohortF[!is.finite(cohortF)] <- 0
      # update the total spawn for fitness
      cohortR[PhenIdx,t+1] = cohortR[PhenIdx,t] + dt*diag(rdi[cohortWprev,PhenIdx])*cohortF
      #cohortR_sol[q,t+1] = dt*rdi[cohortWprev,q] # do not sum the spawn so it is the spawn at time
      #print(cohortR[PhenIdx,1:t+1])
    }
      fitness[which(dimnames(fitness)[[1]] == PhenIdx),which(t_start==cohortSpan)] = cohortR[PhenIdx,T] # fitness is the total spawn within the time period
}
# rownames(fitness) <- PhenName # put the right name
# rownames(fitness) <- round(object@params@species_params$w_mat[PhenIdx],2) # put the maturation size instead of the names
# colnames(fitness) <- cohortSpan

rownames(fitness) <- object@params@species_params$species[PhenIdx] # this is their name
fitness <- fitness[which(rownames(fitness) == iSpecies),]
rownames(fitness) <- round(object@params@species_params$w_mat[which(object@params@species_params$species == iSpecies)],2) # put the maturation size instead of the names

fitness <- fitness[!rowSums(fitness) == 0,] # get rid of phenotypes not appeared yet

plot_dat<- melt(fitness)

p <- ggplot(plot_dat) +
  geom_line(aes(x=cohort,y=value,color = as.factor(species)),size =1) +
  scale_x_continuous(name = "Cohorts", breaks = as.numeric(colnames(fitness))) +
  scale_y_continuous(name = "Fitness", trans = "log10")+
  labs(color = "Phenotypes") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),
        legend.justification=c(1,1),
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

if(save_it) ggsave(plot = p, filename = nameSave)

if (returnData) return(fitness) else if(print_it) return(p)
}



# upgrade of cohort function to run on folders

folder = paste(getwd(),"/SimDefault",sep="")
cohortSpan <- seq(1000,5000,1000)

plotCohortLong <- function(folder, dt = 0.1, t_steps = 5, SpIdx = NULL, iSpecies = NULL, cohortSpan = NULL, print_it = T, returnData = F, save_it = F, comments = T, what = c("normal","fisheries"))
{
  path_to_png = paste(folder,"/FitnessCompare.png",sep="")
  
  plotList <- list() # to store the plots
  plot_datList <- list() #to store the data
  multiSimInit <- bunchLoad(folder = paste(folder,"/init",sep=""))
  
  if(comments) cat("Init sims loaded\n")
  
  for (column in what) # are we plotting normal or fisheries sims ?
  {
    if (column == "normal") 
    {
      listPosition = 1
      multiSim <- bunchLoad(folder = paste(folder,"/normal",sep=""))
      title = c("(a) Without fisheries","(b)","(c)")
      effortInit = 0
      cat("Normal runs\n")
    } else if (column == "fisheries") 
    {
      listPosition = 2
      multiSim <- bunchLoad(folder = paste(folder,"/fisheries",sep=""))
      title = c("(d) With fisheries","(e)","(f)")
      effortInit = multiSim[[1]]@effort[1,1]
      cat("Fishery runs\n")
    } else stop("I don't know what to plot")
    
    
    template = multiSimInit[[1]]
    longSimList <- list()
    timeMax = multiSimInit[[1]]@params@species_params$timeMax[1] + multiSim[[1]]@params@species_params$timeMax[1]
    #prepare the data 
    for (x in 1:length(multiSimInit))
    {
      if(comments)  cat(sprintf("Using run %i\n",x))
      
      # get the right species params
      SummaryParams = template@params # structure and basics params
      SummaryParams@species_params = rbind(multiSimInit[[x]]@params@species_params,multiSim[[x]]@params@species_params) # get the sp ID from both sim
      SummaryParams@species_params$timeMax = timeMax # update the timemax
      a <- SummaryParams@species_params
      a <- a[order(a$ecotype, a$extinct, decreasing=TRUE),] # weird 3 lines to get rid of duplicates and keep the ones with the extinction value
      a <- a[!duplicated(a$ecotype),]
      SummaryParams@species_params = a[order(a$pop,a$ecotype),]
      
      #if (is.null(SpIdx)) SpIdx = sort(unique(SummaryParams@species_params$species))
      
      if (comments) cat("Data handling\n")
      
      result = list(list(multiSimInit[[x]],multiSim[[x]]),SummaryParams) # cannot use finalTOuch exactly as the data as been divided by 10 (dimnames issues)
      
      gc()
      sim = result[[1]] # get the dat
      SummaryParams = result[[2]] # get the params
      rm(result) # get space back
      sim <- sim[lapply(sim, length) > 0] # if a sim is empty
      template = sim[[length(sim)]] # to keep a template of mizer object somewhere
      gc()
      
      # stitiching the sims together
      Dtime = SummaryParams@species_params$timeMax[1] * dt
      Dsp = length(SummaryParams@species_params$ecotype)
      Dw = dim(sim[[1]]@n)[3]
      # put all the sim at the same dimension
      biomList <- list()
      for (i in 1:length(sim)) # for each sim
      {
        biom <- array(data = 0, dim = c(dim(sim[[i]]@n)[1], Dsp, Dw), dimnames = list(dimnames(sim[[i]]@n)$time, SummaryParams@species_params$ecotype, SummaryParams@w)) # create an array of the right dimension
        names(dimnames(biom)) = c("time", "species", "size")
        for (j in dimnames(sim[[i]]@n)$sp) # fill it when necessary
          biom[, which(dimnames(biom)$species == j), ] = sim[[i]]@n[, which(dimnames(sim[[i]]@n)$sp == j), ]
        
        biomList[[i]] <- biom[-1,,] # store it
      }
      biom <- do.call("abind", list(biomList, along = 1)) # abind the list
      names(dimnames(biom)) = list("time", "species", "size")
      dimnames(biom)$time = seq(1, SummaryParams@species_params$timeMax[1]*dt)[-c(length(seq(1, SummaryParams@species_params$timeMax[1]*dt))-1,length(seq(1, SummaryParams@species_params$timeMax[1]*dt)))] # system D
      
      # I have to do the phyto aussi
      phyto <- do.call(abind, c(lapply(sim, function(isim)
        isim@n_pp), along = 1))
      
      # taking care of the effort
      effort <- do.call(rbind, lapply(sim, function(isim)
        isim@effort))
      names(dimnames(effort)) = list("time", "effort")
      dimnames(effort)$time = seq(1, SummaryParams@species_params$timeMax[1]*dt)
      
      # reconstruct the mizer object
      sim = template
      sim@n = biom
      sim@effort = effort
      sim@n_pp = phyto
      
      #need to reconstruct the mizer param object so need to get all the original parameters somehow
      
      sim@params <- MizerParams(SummaryParams@species_params, 
                           max_w=1.1*SummaryParams@species_params[length(unique(SummaryParams@species_params$species)),"w_inf"], #w_inf of biggest species at the start
                           no_w = length(SummaryParams@w), 
                           min_w_pp = SummaryParams@w_full[1], 
                           w_pp_cutoff = round(as.numeric(names(SummaryParams@cc_pp[SummaryParams@cc_pp == 0][1]))),# first size where the background is 0
                           n = 0.75, p = 0.75, q = 0.8, # I never touch these
                           r_pp=4, kappa=0.05, #lambda = (2+q-n), # nor these
                           normalFeeding = F, 
                           tau = 7,
                           interaction = matrix(data = SummaryParams@interaction[1,1],nrow = dim(SummaryParams@species_params)[1],ncol = dim(SummaryParams@species_params)[1],dimnames = list(SummaryParams@species_params$ecotype,SummaryParams@species_params$ecotype)))

      rm(list = "biom", "phyto", "effort")
      gc()
      longSimList[[x]] <- sim
    }
    
    if (comments) cat("Data ready\n")  
    # get the initial stuff
    SpIdx = sort(unique(longSimList[[1]]@params@species_params$species)) # determine spidx if not already given by user
    no_sp = length(SpIdx) # get the species number from this
    if (comments) cat("Initialisation done\n")
    
    # color gradient
    colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
    colGrad <- colfunc(no_sp)
    names(colGrad) <- seq(1,no_sp)
    
    #Fitness
    # get the plots
    FitList <- list()
    if(is.null(cohortSpan)) cohortSpan <- seq(min(dim(longSimList[[1]]@n)[1]),max(dim(longSimList[[1]]@n)[1])-1,500)
    for (i in 1:length(longSimList))
    {
      if (comments) cat(sprintf("Using run %g\n",i))
      FitList[[i]] <- plotCohort(object = longSimList[[i]], returnData = T, iSpecies = iSpecies, t_steps = t_steps, effort = effortInit, cohortSpan = cohortSpan)
    }
    if (comments) cat("Plots loaded")

    # what kind of plot do I want?
    # average over species over sim
    # therefore 9 plots? 1 per species
    
    # need to divide the matrix per species, add their w_mat and then regroup sim per species
    phenList <- list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) # need to do something neater
    for (i in 1:length(FitList))
      for (jSpecies in sort(unique(FitList[[i]]$species))) 
        phenList[[jSpecies]] <- rbind(phenList[[jSpecies]], FitList[[i]][which(FitList[[i]]$species == jSpecies),])

    plot_datList[[listPosition]] <- phenList
   } 
      

    # p <- ggplot(plot_datTrait)+
    #   geom_smooth(aes(x=time,y=percentMean, group = group, color = as.factor(species))) +
    #   scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
    #   scale_y_continuous(name = "Trait relative proportion to initial value", limits = window)+
    #   geom_vline(xintercept = 4000, linetype = "dashed") +
    #   scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    #   theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
    #         panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom", 
    #         #legend.justification=c(1,1),
    #         legend.key = element_rect(fill = "white"))+ 
    #   guides(color = guide_legend(nrow=1)) +
    #   ggtitle(title[1])

  if(returnData) return(plot_datList) else{}

}

cohortSpan <- seq(100,5900,250)
myData <- plotCohortLong(folder = paste(getwd(),"/SimDefault",sep=""), cohortSpan = cohortSpan, returnData = T )
saveRDS(myData,file = "fitnessMat.rds")
myData <- readRDS("fitnessMat.rds")

# Let's play with the data
# myData is a list of normal fitness (1) et fisheries fitness(2)
# inside each list are 9 matrix, one per species, so I lost my multiple sims on the way


#plot_dat <- myData[[1]][[1]] # fitness of all the phenotypes across sims of species 1

plot_dat <- as.data.frame(myData[[2]][[4]])
plot_dat$species <- NULL

b <- plot_dat[!rowSums(plot_dat) == plot_dat$w_mat,]
b[b==0] <- NA

a <- melt(b, id = "w_mat")

a$variable <- as.numeric(as.character(a$variable))

colfunc <- colorRampPalette(c("black", "orange"))
colGrad <- colfunc(length(unique(a$w_mat)))

ggplot(a) +
  geom_line(aes(x=as.numeric(variable),y=value,group = as.factor(w_mat), color = as.factor(w_mat))) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(name = "Time") +
  scale_color_manual(name = "w_mat", values = colGrad)+
  geom_vline(xintercept = 3000, linetype = "dashed") +
  theme(legend.title=element_text(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),
        legend.justification=c(1,1),
        legend.position = "none",
        legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)


