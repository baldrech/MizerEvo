# cohort plot -------------
# I want the spawn through life -> last value of total cohortR

sim <- get(load("eta4and5/9sp/normal/run10/run.Rdata"))
sim <- get(load("eta4and5/3sp/normal/run5/run.Rdata"))
dt=0.1
effort=0
sex_ratio = 0.5

cohortSpan <- seq(1,250,1)
no_sp = dim(sim@params@species_params)[1]
PhenIdx <- which(sim@params@species_params$species == 6)[1] # select who to track
no_Phen = length(PhenIdx)
fitness <- array(0,c(no_Phen, length(cohortSpan)), dimnames = list(PhenIdx,cohortSpan)) #collect the total spawn per time (start of cohort) per species
names(dimnames(fitness)) <- list("species","cohort")
for (cohortGen in cohortSpan)
{
  t_start = cohortGen
  # get rid of the non-existing species at that time
  # SpIdx <- which(sim@n[t_start,PhenIdx,1]>0)
  # a <- as.numeric(names(SpIdx))
  # SpIdx <- a
  # sim@n <- simulation@n[,SpIdx,]
  # no_sp = dim(sim@n)[2]
  
  cat(sprintf("cohort is %g\n",cohortGen))
  T = 5/dt; # number of time steps you want to follow cohort for
  #t_start = laststep - T; # setting the start time for tracking the cohort
  
  cohortW = array(0, c(no_sp, T+1)); # row vector for following cohort weight
  cohortS = array(0, c(no_sp, T+1)); # vector for cohort survival
  cohortR = array(0, c(no_sp, T+1)); # vector for cohort spawning
  cohortR_sol = array(0, c(no_sp, T+1)); # vector for cohort spawn at size
  
  # NEWBORNS OVER LIFETIME
  cohortW[,1] = sim@params@w[1]; # log weight initially (newborn)
  cohortS[,1] = sim@n[t_start,,1]; # initial population in spectrum
  
  
  
  for (q in PhenIdx){ 
    cat(sprintf("q is %g\n",q))
    for (t in seq(1,T)){ # within time period you're interested in
      #cat(sprintf("t is %g\n",t))
      
      cohortWprev = max(which(cohortW[q,t] - sim@params@w >= 0)) # weight bin of cohort from last time step 
      #print(cohortWprev)
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
    fitness[which(dimnames(fitness)[[1]] == q),which(cohortGen==cohortSpan)] = cohortR[q,T]
  }
}

# need to take the average fitness
#dimnames(fitness)$species <- do.call(cbind,lapply(as.numeric(strsplit(as.character(dimnames(fitness)$species), "")), function (x) x[1])) # the name of the phenotypes become their first digit
# dimnames(fitness)$species <- sim@params@species_params$species #or I can do that ^^
# SpIdx <- unique(sim@params@species_params$species)
# 
# newFit <- array(0,c(length(SpIdx), length(cohortSpan)), dimnames = list(SpIdx,cohortSpan))
# for (i in SpIdx)
# {
#   lol <- fitness[which(dimnames(fitness)$species == as.character(SpIdx[i])),]
#   lol[lol==0] <- NA
#   newFit[i,] <- apply(lol,2,mean,na.rm=T)
# }
# 
# newFit[!is.finite(newFit)] <- 0

#cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind


#plot_datF <- melt(fitness)
# plot_datF <- melt(newFit)
# colnames(plot_datF) <- list("species","cohort","fitness")


plot_dat<- melt(fitness)

pF <- ggplot(plot_dat) +
  geom_line(aes(x=cohort,y=value,color = as.factor(species)),size =2) +
  scale_x_continuous(name = "Cohorts") +
  scale_y_continuous(name = "Fitness", trans = "log10")+
  #scale_colour_manual(values=cbPalette)+ # colorblind
  theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"),legend.position="none", 
        legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+
  ggtitle(NULL)

#pB <- plotBiomass(sim,start_time = 10.1,end_time = 300.1)
pB <- plotDynamics(sim)


