# cohort plot -------------
# I want the spawn through life -> last value of total cohortR
plotCohort <- function(object, dt = 0.1, t_steps = 5, iSpecies = 1, effort = 0, cohortSpan = seq(max(dim(object@n)[1])-30,max(dim(object@n)[1]),2),
                       print_it = T, returnData = F, save_it = F, nameSave = paste("CohortSpecies",iSpecies,".png",sep=""))
{
  # setting up some parameters
sex_ratio = 0.5
T = t_steps/dt; # number of time steps you want to follow cohort for
no_sp <- dim(sim@params@species_params)[1]
PhenIdx <- which(sim@params@species_params$species == iSpecies) # select who to track, not the name of the species but their position in the dataframe
PhenName<- sim@params@species_params$ecotype[PhenIdx] # this is their name
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
  cohortW[,1] = sim@params@w[1]; # log weight initially (newborn)
  cohortS[,1] = sim@n[t_start,,1]; # initial population in spectrum

    for (t in seq(1,T)){ # within time period you're interested in
      # vector of the previous size bin for every phenotypes
      cohortWprev = unlist(lapply(lapply(cohortW[PhenIdx,t], FUN = function(x) x-sim@params@w), FUN = function(x) max(which(x>= 0)))) # yolo
      # growth matrix
      growth = getEGrowth(sim@params,n = sim@n[t_start+t-1,,],n_pp = sim@n_pp[t_start+t-1,])
      # update the new size bin with the growth
      cohortW[PhenIdx,t+1] = cohortW[PhenIdx,t]+dt*diag(growth[PhenIdx,cohortWprev])
      # mortality matrix
      z = getZ(object = sim@params, n = sim@n[t_start+t-1,,],n_pp = sim@n_pp[t_start+t-1,], effort = effort)
      # update the amount surviving the time-step
      cohortS[PhenIdx,t+1] = cohortS[PhenIdx,t]*exp(-dt*diag(z[PhenIdx,cohortWprev]))
      # need to prepare n for no NAN, I just want the n of the specific cohort so I extract the right fraction
      n = sim@n[t_start+t-1,,]*cohortS[,t]/cohortS[,1]
      n[!is.finite(n)] <- 0
      # get the rdi manually to have it spread over size bins
      e_spawning <- getESpawning(object = sim@params, n = n,n_pp = sim@n_pp[t_start+t-1,])
      e_spawning_pop <- apply((e_spawning*n),1,"*",sim@params@dw)
      rdi <- sex_ratio*(e_spawning_pop * sim@params@species_params$erepro)/sim@params@w[sim@params@species_params$w_min_idx]
      # update the total spawn for fitness
      cohortR[PhenIdx,t+1] = cohortR[PhenIdx,t] + dt*diag(rdi[cohortWprev,PhenIdx])
      #cohortR_sol[q,t+1] = dt*rdi[cohortWprev,q] # do not sum the spawn so it is the spawn at time
    }
      fitness[which(dimnames(fitness)[[1]] == PhenIdx),which(t_start==cohortSpan)] = cohortR[PhenIdx,T] # fitness is the total spawn within the time period
}
rownames(fitness) <- PhenName # put the right name
rownames(fitness) <- round(sim@params@species_params$w_mat[PhenIdx],2) # put the maturation size instead of the names
colnames(fitness) <- cohortSpan
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

