#give a mizer param object to that thing
set_TBM <- function(no_sp = 10, # number of species #param described in Andersen & Pedersen 2010
                    min_w_inf = 10, # minimum weight of sp
                    max_w_inf = 1e5, # maximum weight of sp
                    no_w = 100, # number of size bins community spectrum
                    min_w = 0.001, #min size bin of community spectrum/The smallest size of the species community size spectrum
                    max_w = max_w_inf * 1.1, #max size bin of both spectrum
                    min_w_pp = 1e-10, #min size bin of background size spectrum
                    no_w_pp = round(no_w)*0.3, # number of size bins background spectrum
                    w_pp_cutoff = 0.5, # cut of size of the background spectrum
                    k0 = 50, # recruitment adjustment parameter
                    n = 0.75, # exponent of maximum intake (scaling of intake)
                    p = 0.75, # exponent of standard metabolism
                    q = 0.8, # exponent of search volume
                    eta = 0.25, # size at maturation relative to Mg (mass in grams ?)
                    r_pp = 4, # growth rate of resource spectrum (primary production)
                    kappa = 0.05, # ressource spectrum carrying capacity
                    lambda = 2+q-n, # exponent of the background spectrum.
                    alpha = 0.6, # assimilation efficiency
                    ks = 2, # factor for standard metabolism
                    z0pre = 0.84, # background mortality factor
                    h = 85, # factor of maximum intake
                    beta = 100, # preferred predator-prey weight ratio
                    sigma = 1.3, # width of selection function
                    f0 = 0.6, # average feeding level of the community/feeding level of small individuals feeding on background
                    knife_edge_size = 1000, #knife edge position
                    gear_names = "knife_edge_gear",
                    r_mult = 1e0, #rmax multiplier to try things
                    cannibalism = 1, # to tweak cannibalism in the interaction matrix
                    erepro = 0.1, # reproduction efficiency
                    rm = NULL, # rmax if want to set up constant
                    s_max = 1000, # time max of the simulation
                    normalFeeding = T, # if wants to normalise the feeding
                    tau = 10, # exponent in psi function
                    interaction = NULL,
                    ...){

  # Calculate gamma using equation 2.1 in Andersen & Pedersen 2010
  alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2) # see A&P 2009
  gamma <- h * f0 / (alpha_e * kappa * (1-f0)) # see A&P 2009 / volumetric search rate
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp) # asymptotic mass of the species
  w_mat <- w_inf * eta # maturation mass / mass at first maturity
  
  #check if things are ok
  # cat(sprintf("beta = %g \n",beta))
  # cat(sprintf("sigma = %g \n",sigma))
  # cat(sprintf("lambda = %g \n",lambda))
  # cat(sprintf("alpha e = %g \n",alpha_e))
  # cat(sprintf("h = %g \n",h))
  # cat(sprintf("f0 = %g \n",f0))
  # cat(sprintf("kappa = %g \n",kappa))
  # cat(sprintf("gamma = %g \n", gamma))
  
  # Check if gears ok
  if (length(knife_edge_size) > no_sp){
    stop("There cannot be more gears than species in the model")
  }
  if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)){
    warning("Number of gears is less than number of species so gear information is being recycled. Is this what you want?")
  }
  if ((length(gear_names) != 1) & (length(gear_names) != no_sp)){
    stop("Length of gear_names argument must equal the number of species.")
  }
  
  # Make the species parameters data.frame
  trait_params_df <- data.frame(
    species = 1:no_sp,
    w_inf = w_inf,
    w_mat = w_mat,
    h = h, # max food intake
    gamma = gamma, # vol. search rate,
    ks = ks,# standard metabolism coefficient,
    beta = beta,
    sigma = sigma,
    eta = eta,
    z0 = z0pre * w_inf^(n-1), # background mortality
    alpha = alpha,
    w_min = min_w,
    sel_func = "knife_edge",
    knife_edge_size = knife_edge_size,
    gear = gear_names,
    erepro = erepro, # not used but included out of necessity
    extinct = FALSE,
    cannibalism = cannibalism,
    pop = 0, # to get the time of apparition
    run = 1,
    ecotype = 1:no_sp,
    error = 0 # to trace errors
  )
  # Make the MizerParams
  # MizerParams is in MizerParams-class. Use Source or something because it's freaking long.
  trait_params <- MizerParams(trait_params_df, max_w=max_w, no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_pp_cutoff, n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda, normalFeeding = normalFeeding, tau = tau, interaction = interaction) 
  # Sort out maximum recruitment - see A&P 2009
  # Get max flux at recruitment boundary, R_max
  # R -> | -> g0 N0
  # R is egg flux, in numbers per time
  # Actual flux at recruitment boundary = RDD = NDD * g0 (where g0 is growth rate)
  # So in our BH SRR we need R_max comparable to RDI (to get RDD)
  # R_max = N0_max * g0 (g0 is the average growth rate of smallest size, i.e. at f0 = 0.5)
  # N0 given by Appendix A of A&P 2010 - see Ken's email 12/08/13
  # Taken from Ken's code 12/08/13 - equation in paper is wrong!

  if (is.null(rm))
  {
    alpha_p <- f0 * h * beta^(2 * n - q - 1) * exp((2 * n * (q - 1) - q^2 + 1) * sigma^2 / 2)
    alpha_rec <- alpha_p / (alpha * h * f0 - ks)
    # Calculating dw using Ken's code - see Ken's email 12/08/13
    tmpA <- w_inf[1]
    tmpB <- (log10(w_inf[length(w_inf)]) - log10(w_inf[1])) / (no_sp - 1) # Difference between logged w_infs, fine
    
    #if (length(no_sp) == 1 ) dw_winf <-  tmpA *10
    if (no_sp == 1 ) dw_winf <-  tmpA *10
    
    else  dw_winf <- tmpB * tmpA *10^(tmpB*((1:no_sp)-1)) # ?
    
    N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
    # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
    g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
    r_max <- N0_max * g0 * r_mult
    
    trait_params@species_params$r_max <- r_max
  }
  else trait_params@species_params$r_max <- rm
  
  # addition of the maximum time of the simulation to have it somewhere in the mizer object
  trait_params@species_params$timeMax <- s_max
  
  return(trait_params)
}
# what's in it?
#param <- set_TBM(no_sp = 3)
#slotNames(param)
# w: start value of each size bins of the community spectrum
# dw: lenght of the size bin (w1+dw=w2)
#dw_full: no idea, might a limit size for the bins
#psi: allocation to reproduction #psi to std_metabolism are fucntion applied on w, therefore one value for each sp of each size bin
#intake_max:
#search_vol:
#activity:
#std_metab:
#pred_kernel: huge stuff
#rr_pp: background size spectrum growth rate (weight specific)
#cc_pp: background size spectrum ressource carrying capacity
#species_params: summary of the sp params. you can get the type of gear use to catch them (only one for now)
#interaction: does the sp interact or not (I guess)
#srr: Beverton Holt esque relationship (function) didn't know what that do
#selectivity: does the gear catch or not depending on the sp/w
#catchability: eeeh I dunno

#summary(param)
#use param@something to explore slots

# now the projection in time
#there are no fisheries for the moment, issue in the fMortGear

#for testing purpose
# effort = 0
# t_max = 100
# t_save = 1
# dt = 0.1
# initial_n=get_initial_n(param)
# initial_n_pp=param@cc_pp

project <-  function(object, effort=0,  t_max = 100, t_save=0.1, dt=0.1, initial_n=get_initial_n(object), initial_n_pp=object@cc_pp, 
                     mu = 2, i_stop = NULL, resident = NULL, data = FALSE, extinct = TRUE, RMAX = TRUE, OptMutant = "M1", M3List = NULL ,
                     checkpoint, print_it, predMort = NULL, ...){
  
  
  umbrella = FALSE # parameter that says if there are still things alive
  
  #first, let's convert the effort to the good dim/class
  if(class(effort) == "numeric"){
    if (!all((t_max %% dt) == 0)) # %% is the remainder
      stop("t_max must be divisible by dt with no remainder")
    
    no_gears <- dim(object@catchability)[1] #number of gears
    if ((length(effort)>1) & (length(effort) != no_gears))
      stop("Effort vector must be the same length as the number of fishing gears\n")
    
    # If more than 1 gear need to check that gear names match
    gear_names <- dimnames(object@catchability)[[1]]
    effort_gear_names <- names(effort)
    
    if (length(effort) == 1 & is.null(effort_gear_names)){
      effort_gear_names <- gear_names
    }
    
    if(!all(gear_names %in% effort_gear_names)){
      gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort vector.", sep="")
      stop(gear_names_error_message)
    }
    
    # Set up the effort array transposed so we can use the recycling rules
    #time_dimnames <- signif(seq(from=1,to=t_max,by=dt),3)
    time_dimnames <- signif(seq(from=1*dt,to=(t_max/dt)*dt,by=dt),3) #if I keep the previous one, I' missing the first 8 time step, which should be the inverse situation # deleted the -8 and it did not change anything (fingers crossed)
    # that's super weird
    
    effort_array <- t(array(effort, dim=c(no_gears,length(time_dimnames)), dimnames=list(gear=effort_gear_names,time=time_dimnames)))
    effort <- effort_array
  }
  # now we have the effort check an in the array format
  # print(class(effort))
  # print(dim(effort))
  # print(dimnames(effort))
  
  
  validObject(object)
  
  # Check that number and names of gears in effort array is same as in MizerParams object
  no_gears <- dim(object@catchability)[1]
  if(dim(effort)[2] != no_gears){
    no_gears_error_message <- paste("The number of gears in the effort array (length of the second dimension = ", dim(effort)[2], ") does not equal the number of gears in the MizerParams object (", no_gears, ").", sep="")
    stop(no_gears_error_message)
  }
  gear_names <- dimnames(object@catchability)[[1]]
  if(!all(gear_names %in% dimnames(effort)[[2]])){
    gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort array.", sep="")
    stop(gear_names_error_message)
  }
  # Sort effort array to match order in MizerParams
  effort <- effort[,gear_names, drop=FALSE]
  
  # Blow up time dimension of effort array
  # i.e. effort might have been passed in using time steps of 1, but actual dt = 0.1, so need to blow up
  if (is.null(dimnames(effort)[[1]])){
    stop("The time dimname of the effort argument must be numeric.")
  }
  if (any(is.na(as.numeric(dimnames(effort)[[1]])))){
    stop("The time dimname of the effort argument must be numeric.")
  }
  time_effort <- as.numeric(dimnames(effort)[[1]])
  
  ##_max <- time_effort[length(time_effort)] # commenting that because bugs
  
  # Blow up effort so that rows are dt spaced
  time_effort_dt <- seq(from = time_effort[1], to = t_max, by = dt)
  
  effort_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(effort)[2]), dimnames=list(time = time_effort_dt, dimnames(effort)[[2]])))
  
  for (i in 1:length(time_effort)){
    effort_dt[,time_effort_dt >= time_effort[i]] <- effort[i,]
  }
  effort_dt <- t(effort_dt)
  
  #now the effort is done, let's do something interesting
  
  # Make the MizerSim object with the right size
  # We only save every t_save steps, default is 1
  if (!all((t_save %% dt)  == 0))
    stop("t_save must be divisible by dt with no remainder")
  t_dimnames_index <- as.integer(seq(from = 1+ ((t_save-1) / dt), to = length(time_effort_dt), by = t_save/dt)) #create a vector of all the times step where there is a save (every 10 dt)
  t_dimnames_index <- t_dimnames_index[t_dimnames_index>0] #get rid of non positive if so
  t_dimnames <- time_effort_dt[t_dimnames_index]
  sim <- MizerSim(object, t_dimnames = t_dimnames) #build the object, pretty much empty at this stade but has the right dimensions to be filled
  # I dont know why but the mizer object is created 1 time step longer (add time step 0) and its a fucking pain in the ass
  # Fill up the effort array
  sim@effort[] <- effort_dt[t_dimnames_index,]
  
  # # Set initial population
  #   sim@n[1, , ] <- initial_n
  #   sim@n_pp[1, ] <- initial_n_pp
  
  # Handy things
  no_sp <- nrow(sim@params@species_params)
  no_w <- length(sim@params@w)
  no_w_pp <- length(sim@params@w_full)
  idx <- 2:no_w
  
  # If no w_min_idx column in species_params, add one
  if (!("w_min_idx" %in% names(sim@params@species_params)))
    sim@params@species_params$w_min_idx <- 1
  
  # Hacky shortcut to access the correct element of a 2D array using 1D notation
  # this thing get you the first value of each species in whatever function you want
  w_min_idx_array_ref <- (sim@params@species_params$w_min_idx-1) * no_sp + (1:no_sp)
  
  # sex ratio - DO SOMETHING LATER WITH THIS
  sex_ratio <- 0.5
  
  # Matrices for solver
  # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
  # A <- matrix(0,nrow=no_sp,ncol=no_w)
  # B <- matrix(0,nrow=no_sp,ncol=no_w)
  # S <- matrix(0,nrow=no_sp,ncol=no_w)
  
  # new version
  A = matrix(data=0, nrow=no_sp, ncol=no_w-1)
  B = matrix(data=0, nrow=no_sp, ncol=no_w)
  C = rep(0,no_w-1)
  S = matrix(data=0, nrow=no_sp, ncol=no_w)
  
  if(dim(sim@n)[2] == 1) dimnames(sim@n)$sp = 1
  else   dimnames(sim@n)$sp = rownames(initial_n) # the object created by mizer doesnt keep in memomry my mutant names, so Im putting them here
  
  # Set initial population
  if (missing(i_stop) == TRUE)
  {
    sim@n[1, , ] <- initial_n
    sim@n_pp[1, ] <- initial_n_pp
    # initialise n and nPP (pp is background)
    # We want the first time step only but cannot use drop as there may only be a single species
    n <- array(sim@n[1,,],dim=dim(sim@n)[2:3]) #take the first line of sim (for each weight) and put it in the matrix of right dimension (= sim@n at t=1)
    dimnames(n) <- dimnames(sim@n)[2:3] # now it has the weights as names
    n_pp <- sim@n_pp[1,] # no need for an array, there is only one line
    t_steps <- dim(effort_dt)[1] #time stpes = max number of dt (not only the 100 saved)
    init = 1 # for the for loop
  }
  
  else # this loop allow to continue the simulation where it stopped previously (on a time step point of view), if it has
  {
    # i_stop is not the real time step that the user enter, i_stop = i_step/dt
    # I'm going to start at the next i_step then
    t_init = i_stop# round up / 
    
    # Set initial population
    
    dimnames(sim@n)[[2]] <- rownames(initial_n) # updating the names accordingly (could do that during the object creation)
    sim@n[t_init+1, , ] <- initial_n #+1 is there because the array start at 0 and not 1 # name bug here
    sim@n_pp[t_init+1, ] <- initial_n_pp
    # initialise n and nPP (pp is background)
    # We want the first time step only but cannot use drop as there may only be a single species
    n <- array(sim@n[t_init+1,,],dim=dim(sim@n)[2:3]) #take the first line of sim (for each weight) and put it in the matrix of right dimension (= sim@n at t=1)
    dimnames(n) <- dimnames(sim@n)[2:3] # now it has the weights as names
    n_pp <- sim@n_pp[t_init+1,] # no need for an array, there is only one line
    t_steps <- dim(effort_dt)[1] #time steps = max number of dt (not only the 100 saved)
    init = i_stop+1 # I need to start right after the last save (t_init), but  pass in small steps
  }
  # the sim is fully initialised now, time to move forwards           
  #time projection
  
  if (data == TRUE){
    # arrays that gets the details of energy allocation (only works without mutants) c(as.character(seq(1:no_sp)))
    energy <- array(dim = c(t_steps,no_sp,no_w,4), dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp,dimnames(sim@n)$w,c("reproAndGrowth", "spawning", "growth", "feeding")))
    names(dimnames(energy)) <- list("Time","Species","Size","Energy")
    rd <- array(dim = c(t_steps,no_sp,2), dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp, c("RDI", "RDD")))
    names(dimnames(rd)) <- list("Time","Species","Energy")
    eggs <- array(dim = c(t_steps,no_sp), dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp))
    names(dimnames(eggs)) <- list("Time","Species")
    food <-array(dim = c(t_steps,no_sp,no_w,no_w_pp), dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp,dimnames(sim@n)$w,dimnames(sim@n_pp)$w))
    names(dimnames(food)) <- list("Time","Species","PredSize","PreySize")
    death <-array(dim = c(t_steps,no_sp,no_w), dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp,dimnames(sim@n)$w ))
    names(dimnames(death)) <- list("Time","PreySp","PreySize")
    Tdeath <-array(dim = c(t_steps,no_sp,no_w), dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp,dimnames(sim@n)$w ))
    names(dimnames(Tdeath)) <- list("Time","PreySp","PreySize")
    Pdeath <- array(dim = c(t_steps,no_w_pp),dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n_pp)$w))
    names(dimnames(Pdeath)) <- list("Time","PreySize")
    trouveF <- array(dim = c(t_steps,no_sp,no_w),  dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp,dimnames(sim@n)$w ))
    names(dimnames(trouveF)) <- list("Time","PredSp","PredSize")
    trouveB <- array(dim = c(t_steps,no_sp,no_w),  dimnames = list(c(as.character(seq(1:t_steps))),dimnames(sim@n)$sp,dimnames(sim@n)$w ))
    names(dimnames(trouveB)) <- list("Time","PredSp","PredSize")
  }
  
  for (i_time in init:t_steps)
  {
    # if (i_time %% check_point == 0) return() # stop the simulation every 500 loop to reduce the size of the arrays and clean the extinct species
    # Do it piece by piece to save repeatedly calling methods, functions found in porject_methods.r
    phi_prey <- getPhiPrey(sim@params, n=n, n_pp=n_pp,opt = T)
    
    feeding_level <- getFeedingLevel(sim@params, n=n, n_pp=n_pp, phi_prey=phi_prey)
    pred_rate <- getPredRate(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
    
    if (!is.null(predMort)) m2 = predMort else m2 <- getM2(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate) # for cst mortality rate (I think)
    # print(class(m2))
    # print(dim(m2))
    # print(dimnames(m2))
    # print(m2)
    m2_background <- getM2Background(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate)
    #print(effort_dt[i_time,])
    z <- getZ(sim@params, n=n, n_pp=n_pp, effort=effort_dt[i_time,], m2=m2) #total mortality 
    e <- getEReproAndGrowth(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
    e_spawning <- getESpawning(sim@params, n=n, n_pp=n_pp, e=e)
    e_growth <- getEGrowth(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, e=e)
    rdi <- getRDI(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, sex_ratio=sex_ratio)
    rdd <- getRDD(sim@params, n=n, n_pp=n_pp, rdi=rdi, sex_ratio=sex_ratio)
    
    # Iterate species one time step forward:
    # A[,idx] <- sweep(-e_growth[,idx-1,drop=FALSE]*dt, 2, sim@params@dw[idx], "/") 
    # # idx start at 2, the -1 makes it include the first column # the "-" makes all the value negative and *dt reduce accordingly to one time step
    # # the operation takes the first column of e_growth and divide it by the first column of sim@params@dw (which is in reality the secon one because of the idx-1)
    # #the result is a negative value placed in the second column of A (dw is small so dividing by it makes a big number)
    # 
    # B[,idx] <- 1 + sweep(e_growth[,idx,drop=FALSE]*dt,2,sim@params@dw[idx],"/") + z[,idx,drop=FALSE]*dt
    # # in this one, with start with the second column of e_growth, divided by the same of sim@params@dw
    # # why m2 pos? # I think it's the sum of everything that leaves w, hence the growth that goes in the next and predation
    # 
    # S[,idx] <- n[,idx,drop=FALSE]
    # # Boundary condition upstream end (recruitment), add values to the first column that stayed empty
    # 
    # B[w_min_idx_array_ref] <- 1+e_growth[w_min_idx_array_ref]*dt/sim@params@dw[sim@params@species_params$w_min_idx]+z[w_min_idx_array_ref]*dt
    # 
    # # Update first size group of n
    # #actual value + density dependent reproduction by dt / dw (bins size) / first column of B
    # if (RMAX == FALSE)  n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdi*dt/sim@params@dw[sim@params@species_params$w_min_idx]) / B[w_min_idx_array_ref] #changed rdd to rdi
    # else  n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdd*dt/sim@params@dw[sim@params@species_params$w_min_idx]) / B[w_min_idx_array_ref]
    # # print("rmax")
    # # print(sim@params@species_params$r_max)
    # # 
    # # print("rdi")
    # # print(rdi)
    # # print("rdd")
    # # print(rdd)
    # # 
    # # print("dt machin")
    # # print(dt/sim@params@dw[sim@params@species_params$w_min_idx])
    # # 
    # # print("B")
    # # print(B[w_min_idx_array_ref])
    # # Invert matrix
    # for (i in 1:no_sp)
    #   for (j in (sim@params@species_params$w_min_idx[i]+1):no_w) # the 2 loops sweep all the matrix n, change the whole n array at one time step, start at 2
    #     n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
    
    # new version
    #
    # Set up matrix:
    #
    for (iSpecies in 1:no_sp) {
      A[iSpecies, ] <- -e_growth[iSpecies, 1:(no_w-1)]*dt/sim@params@dw[2:no_w]
      B[iSpecies, ] <- 1 + dt*(e_growth[iSpecies,]/sim@params@dw + z[iSpecies,])
      S[iSpecies, ] <- n[iSpecies, ]
      
      if (RMAX) S[iSpecies, 1] <- n[iSpecies, 1] + rdd[iSpecies]*dt/sim@params@dw[1]
      else S[iSpecies, 1] <- n[iSpecies, 1] + rdi[iSpecies]*dt/sim@params@dw[1]
    }
    
    #
    # Invert matrix:
    #
    for (iSpecies in 1:no_sp) {
      n[iSpecies,] <- Solve.tridiag(A[iSpecies,],B[iSpecies,],C,S[iSpecies,])
    }
    

    
    
    
    
    
    

    
    # extinction part
    if (extinct == TRUE)
    {
      extinction = 1e-30
      # remove all rows with non-finite values
      n[!rowSums(!is.finite(n)), ]
      # replace all non-finite values with 1e-30 (not 0 but lower than extinction threshold)
      n[!is.finite(n)] <- 1e-30
      
      for (i in 1:no_sp)
      {
        if (sum(n[i,]) < extinction &
            0 < sum(n[i,]))
          # if species abundance under extinction threshold but not already extinct, kill it
        {
          n[i,] = 0
          # find the name of the species going extinct
          toto = which(sim@params@species_params$ecotype == rownames(n)[i])
          
          if (sim@params@species_params$extinct[toto] == FALSE)
            # security for bugs
            if (print_it) cat(
              sprintf(
                "Extinction of species %s at time %s\n",
                sim@params@species_params$ecotype[toto],
                i_time
              )
            )
          
          else if (sim@params@species_params$extinct[toto] != FALSE)
          {
            if (print_it)  cat(
              sprintf(
                "Species %s at time %s is a zombie\n",
                sim@params@species_params$ecotype[toto],
                i_time
              )
            ) # to check if they come back from the dead
            sim@params@species_params$erro[toto] = 1 # if this happen it will be noted by a 1 in the sp ID
          }
          
          sim@params@species_params$extinct[toto] <-
            i_time + (checkpoint - 1) * t_max / dt # update the extinction status
          #print(sim@params@species_params)
          if (sim@params@species_params$extinct[toto] < sim@params@species_params$pop[toto])
            sim@params@species_params$error[toto] = 2 # if this happen it will be noted by a 2 in the sp ID
        }
      }
      if (dim(sim@params@species_params[sim@params@species_params$extinct != FALSE,])[1] == dim(sim@params@species_params)[1])
        umbrella = TRUE # if this is true, evrything is dead
    }
    
    #why the -A ? why not directly construct A the right way?
    #here, B = S -A*n(-1) +e_growth + m2, it's like the total of everything and we get the portion that stay in the bin ?
    
    # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
    tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
    n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp+m2_background)*dt)
    
    # time to save!, is i_time in t_dimnames?
    store <- t_dimnames_index %in% i_time # test if t is in i
    
    
    if (any(store))
    {
      sim@n[which(store)+1,,] <- n # 'which' tells how many true are in store, indicate the time step where to store the n
      sim@n_pp[which(store)+1,] <- n_pp
    }
    
    if (umbrella == TRUE) # in that case nothing is left and the simulation stop
    {
      if (print_it) cat(sprintf("Life has left your simulation, game over.\nSimulation stopped at time %s.\n", i_time))
      if (data == TRUE) return(list(energy,rd,eggs,sim,food,death,Tdeath,Pdeath,trouveF,trouveB)) # when I want to run the sim but get something else from it (like any other data)
      
      else return(list(sim,umbrella)) #I just want something size 2
    }
    
    # getting the energy allocation data
    if (data == TRUE)
    {
      energy[i_time, ,, ] =  cbind(e, e_spawning, e_growth,feeding_level)
      rd[i_time,,] = cbind(rdi,rdd)
      eggs[i_time,] = n[w_min_idx_array_ref]
      food[i_time,,,] = pred_rate
      death[i_time,,] = m2
      Tdeath[i_time,,] = z
      Pdeath[i_time,] = m2_background
      trouveF[i_time,,] = phi[[1]]
      trouveB[i_time,,] = phi[[2]]
    }
    
    
    # MUTANT TIME
    mute = FALSE
    multiple = FALSE
    switch(OptMutant,
           
           M2 = { # default mutation rate, with one mutant max per time step, randomly drawn from every species (not phenoytpes)
             if (mu >= sample(1:1000, 1))
             {
               residentPool = sim@params@species_params[sim@params@species_params$extinct == FALSE,] # only keep the available residents (the one not extinct)
               # to block exponential evolution of species, I'm first picking a lineage randomly and then an ecotype in this lineage
               lineagePool = unique(residentPool$species)

               if (length(lineagePool) == 1) lineage = lineagePool else lineage = sample(lineagePool, 1)
               
               residentPool=residentPool[residentPool$species == lineage,]
               resident <- sample(1:nrow(residentPool), 1) # this is the rownumber of the selected resident
               resident <- residentPool$ecotype[resident] # this is his name now
                mute = TRUE
             }
           },
           M3 = { # if the user define a specific time to mutate, mainly for debugging
             for (i in 1:length(M3List[[1]]))
             {
               if (M3List[[1]][i] == i_time)
               {
                 # old version
                 #   residentPool = sim@params@species_params[sim@params@species_params$extinct == FALSE,] # only keep the available residents (the one not extinct)
                 # resident <- sample(1:nrow(residentPool), 1) # this is the rownumber of the selected resident
                 # resident <- rownames(residentPool)[resident] # this is his name now
                 #new version
                 residentPool = sim@params@species_params[sim@params@species_params$extinct == FALSE,] # only keep the available residents (the one not extinct)
                 #print(residentPool)
                 # to block exponential evolution of species, I'm first picking a lineage randomly and then an ecotype in this lineage
                 lineagePool = unique(residentPool$species)
                 #print(lineagePool)
                 
                 if (length(lineagePool) == 1) lineage = lineagePool
                 
                 else lineage = sample(lineagePool, 1)
                 
                 #print(lineage)
                 residentPool=residentPool[residentPool$species == lineage,]
                 #print(residentPool)
                 resident <- sample(1:nrow(residentPool), 1) # this is the rownumber of the selected resident
                 #print(resident)
                 
                 resident <- residentPool$ecotype[resident]
                 mute = TRUE
               }}
           },
           M4 = { # multiple residents at one time, not sure if it still works
             residentPool = sim@params@species_params[sim@params@species_params$extinct == FALSE,]
             resident = NULL
             for (i in 1:nrow(residentPool))
             {
               if (mu >= sample(1:1000, 1))
               {
                 resident <- c(resident, rownames(residentPool)[i]) # this is his name now
                 mute = TRUE
               }
               if (length(resident) >1) multiple = TRUE
             }
           },
           M5 = { # pick only one mutant but give a chance to every species
             residentPool = sim@params@species_params[sim@params@species_params$extinct == FALSE,] # only keep the available residents (the one not extinct)
             speciesPool = unique(residentPool$species) # which species are available to produce new phenotypes
             challengers <- NULL
             for (iSpecies in speciesPool) # do the picking for every species
             {
               if (mu >= sample(1:1000, 1)) # if mutant happens
               {
                 resident <- sample(residentPool[residentPool$species == iSpecies,]$ecotype, 1) # get the name of one phenotype in the selected species
                 mute = TRUE
                 challengers <- c(challengers,resident)
               } 
             }
             if (length(challengers) >1){
             cat(sprintf("Possible new phenotypes\n"))
             print(challengers)
             resident <- sample(challengers,1) # select only one to mutate (I know I'm lazy)
             } else if (length(challengers == 1)) resident <- challengers
             },
           {})
    
    
    
    if (mute == TRUE & i_time!=t_steps) { # mutation rate egg dependent # if  Iget a mutant on last time step I get bugs because the sim restart at last +1 time step
      # save the data
      sim_stop = sim
      # I need to get rid of the first or last line for no overlapping
      # it will be the first (initial conditions of this sim), but only after the first mutation
      # t_init +1 is the line number of the initialisation
      if (missing(i_stop) == FALSE) sim_stop@n[t_init+1,,]<- NA # if it's not the first run, delete the initialisation (first line where the mutant is introduce, easier for pasting later)
      
      i_stop = i_time  # to conserve the time of the projection to restart later
      stopList <- list(sim_stop, i_stop, resident,n,n_pp)
      names(stopList) <- c("data", "i_stop", "resident","n","n_pp")
      
      if(length(stopList)!=5) cat(sprintf("error in stop list, length is %i\n",length(stopList)))
      # now I need to leave the projection and keep resident, i_stop and sim_stop
      if (multiple == FALSE)
      {
        if (print_it) cat(sprintf(
          "A mutant from species %s has appeared at time %s\n",
          resident,
          i_time
        ))}
      
      else 
      {
        if (print_it) cat(sprintf(
          "Mutants from species %s have appeared at time %s\n",
          resident,
          i_time
        ))}
      
      
      
      return(stopList)
      
      
    }  
    
    
    
  }
  # and end
  if (missing(i_stop) == FALSE) sim@n[t_init+1,,]<- NA # need to get rid of the initialisation for the last run before exiting
  # I'm keeping the if to not have this enable during the initialisation phase
  # I'm assuming that I have at least one mutant per run
  
  if (data == TRUE) return(list(energy,rd,eggs,sim,food,death,Tdeath,Pdeath,trouveF,trouveB)) # when I want to run the sim but get something else from it (like any other data)
  
  else return(sim)
}



