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


#debugging
# fitness2 <- fitness #save
# rownames(fitness2) <- PhenName # put the right name
# PhenIdx <- which(object@params@species_params$species == 9)
# SpName<- object@params@species_params$species[PhenIdx] # this is their name
# rownames(fitness2) <- SpName # put the right name
# 
# a <- fitness2[which(rownames(fitness2) == 9),]
# b <- getEReproAndGrowth(object = object@params, n = n,n_pp = object@n_pp[t_start+t-1,])
# c <- getE(object = object@params, n = n,n_pp = object@n_pp[t_start+t-1,])
# 
# 
# energyMat <- NULL
# time_range <- cohortSpan +195.1
# for (i in time_range)
# {
#   #time_range = max(as.numeric(dimnames(object@n)$time))
#   time_elements <- get_time_elements(sim,i)
#   myData <- aaply(which(time_elements), 1, function(x){
#     # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
#     
#     n <- array(sim@n[x,,],dim=dim(sim@n)[2:3])
#     dimnames(n) <- dimnames(sim@n)[2:3]
#     myData <- getE(sim@params, n=n, n_pp = sim@n_pp[x,])
#     return(myData)})
#   
#   
#   dimnames(myData)$sp = sim@params@species_params$species
#   #SpIdx = sort(unique(sim@params@species_params$species)) # get the species names
#   SpIdx = 9
#   growth_sp = matrix(data = NA, ncol = dim(myData)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(myData)$w)) # prepare the new sim
#   names(dimnames(growth_sp))=list("species","size")
#   
#   for (i in SpIdx)
#   {
#     temp = myData # save to manip
#     temp[which(rownames(myData) != i), ] = 0 # keep the ecotypes from the species only
#     temp = apply(temp, 2, sum)
#     temp = temp / length(which(rownames(myData)==i)) # do the mean (in 2 steps)
#     growth_sp[which(rownames(growth_sp)==i), ] = temp
#   }
#   energyMat = rbind(energyMat,growth_sp)
# }
# 
# # Calculate energy when intake is maximum (feeding level = 1)
# getEmax <-  function(sim, n, n_pp, feeding_level){
#   e <- sweep(sim@intake_max,1,sim@species_params$alpha,"*")
#   e[e<0] <- 0
#   return(e)}
# 
# time_elements <- get_time_elements(sim,100.1)
# myData <- aaply(which(time_elements), 1, function(x){
#   n <- array(sim@n[x,,],dim=dim(sim@n)[2:3])
#   dimnames(n) <- dimnames(sim@n)[2:3]
#   myData <- getEmax(sim@params, n=n, n_pp = sim@n_pp[x,])
#   return(myData)})
# 
# growthMax = myData[1,]
# growthMaxDf <- data.frame(size = as.numeric(sim@params@w), value = growthMax)
# 
# rownames(energyMat) <- time_range
# plot_dat <- melt(energyMat)
# colnames(plot_dat) <- c("time","size","value")
# 
# colfunc <- colorRampPalette(c("black", "orange"))
# colGrad <- colfunc(length(unique(plot_dat$time)))
# 
# # metabolism
# ks <- sim@params@species_params$ks[1]
# p = 0.75
# metab <-  unlist(tapply(sim@params@w,1:length(sim@params@w),function(wx,ks,p)ks * wx^p, ks=ks,p=p))
# metabDF <- data.frame(size = as.numeric(sim@params@w), value = metab)
# 
# p <- ggplot(plot_dat) + 
#   geom_line(aes(x=size, y = value, colour = as.factor(time))) + 
#   geom_line(data = metabDF, aes(x = size,y = value), color = "red") +    
#   geom_line(data = growthMaxDf, aes(x = size,y = value), color = "green") + 
#   geom_vline(xintercept = sim@params@species_params$w_mat[9], linetype = "dashed") +
#   scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
#   scale_y_continuous(name = "Value", trans ="log10")+
#   scale_color_manual(name = "Time", values = colGrad)+
#   theme(legend.title=element_blank(),
#         legend.justification=c(1,1),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white", color = "black"),
#         panel.grid.minor = element_line(colour = "grey92"))+
#   ggtitle("Energy before metabolism")
# 
# p
