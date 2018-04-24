#function used for the time projection

#----------------------------------------------------------------
# Class specification and constructors for the simulation class
#----------------------------------------------------------------

#check if there is an error somewhere in the object (mizer param)
valid_MizerSim <- function(object){
  errors <- character()
  validObject(object@params)
  # array dimensions
  if(length(dim(object@n)) != 3){
    msg <- "n slot must have three dimensions"
    errors <- c(errors, msg)
  }
  if(length(dim(object@effort)) != 2){
    msg <- "effort slot must have two dimensions"
    errors <- c(errors, msg)
  }
  if(length(dim(object@n_pp)) != 2){
    msg <- "n_pp slot must have two dimensions"
    errors <- c(errors, msg)
  }
  # Check time dimension is good - size, dim name, and names
  if(!all(c(dim(object@n)[1]-1,dim(object@n_pp)[1]-1) == dim(object@effort)[1])){
    msg <- "First dimension of effort, n and n_pp slots must be the same length (n and n_pp are + 1 longer)"
    errors <- c(errors, msg)
  }
  if(!all(c(names(dimnames(object@n))[1], names(dimnames(object@n_pp))[1], names(dimnames(object@effort))[1]) == "time")){
    msg <- "First dimension of effort, n and n_pp slots must be called 'time'"
    errors <- c(errors, msg)
  }
  if(!all(c(names(dimnames(object@n))[1], names(dimnames(object@n_pp))[1]) == 
          names(dimnames(object@effort))[1])){
    msg <- "First dimension of effort, n and n_pp slots must have the same names"
    errors <- c(errors, msg)
  }
  # species dimension of n
  if(dim(object@n)[2] != dim(object@params@psi)[1]){
    msg <- "Second dimension of n slot must have same length as the species names in the params slot"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@n))[2] != "sp"){
    msg <- "Second dimension of n slot must be called 'sp'"
    errors <- c(errors, msg)
  }
  if(!all(names(dimnames(object@n))[2] == names(dimnames(object@params@psi))[1])){
    msg <- "Second dimension of n slot must have same species names as in the params slot"
    errors <- c(errors, msg)
  }
  # w dimension of n
  if(dim(object@n)[3] != length(object@params@w)){
    msg <- "Third dimension of n slot must have same length as w in the params slot"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@n))[3] != "w"){
    msg <- "Third dimension of n slot must be called 'w'"
    errors <- c(errors, msg)
  }
  if(!all(names(dimnames(object@n))[3] == names(dimnames(object@params@psi))[2])){
    msg <- "Third dimension of n slot must have same size names as in the params slot"
    errors <- c(errors, msg)
  }
  # w dimension of n_pp
  if(dim(object@n_pp)[2] != length(object@params@w_full)){
    msg <- "Second dimension of n_pp slot must have same length as w_full in the params slot"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@n_pp))[2] != "w"){
    msg <- "Second dimension of n_pp slot must be called 'w'"
    errors <- c(errors, msg)
  }
  if(!all(dimnames(object@n_pp)$w == names(object@params@rr_pp))){
    msg <- "Second dimension of n_pp slot must have same size names as rr_pp in the params slot"
    errors <- c(errors, msg)
  }
  # gear dimension of effort
  if(dim(object@effort)[2] != dim(object@params@catchability)[1]){
    msg <- "Second dimension of effort slot must have same number of gears as in the params slot"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@effort))[2] != "gear"){
    msg <- "Second dimension of effort slot must be called 'gear'"
    errors <- c(errors, msg)
  }
  if(!all(names(dimnames(object@effort))[2] == names(dimnames(object@params@catchability)[1]))){
    msg <- "Second dimension of effort slot must have same gear names as in the params slot"
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}

#set a new class, the simulation class
setClass("MizerSim",
         representation(
           params = "MizerParams",
           n = "array",
           effort = "array",
           n_pp = "array"
         ),
         prototype = prototype(
           params = new("MizerParams"),
           n = array(NA,dim=c(1,1,1), dimnames = list(time = NULL, sp = NULL, w = NULL)),
           effort = array(NA,dim=c(1,1), dimnames = list(time = NULL, gear = NULL)),
           n_pp = array(NA,dim=c(1,1), dimnames = list(time = NULL, w = NULL))
         ),
         validity = valid_MizerSim
)
setValidity("MizerSim", valid_MizerSim)
remove(valid_MizerSim)

# the constructor

MizerSim <- function(object, t_dimnames = NA, t_max = 100, t_save=1, ...){
  # If the dimnames for the time dimension not passed in, calculate them from t_max and t_save 
  if (any(is.na(t_dimnames))){
    if((t_max %% t_save) != 0)
      stop("t_max must be divisible by t_save with no remainder")
    t_dimnames <- seq(from = t_save, to = t_max, by = t_save)
  }
  if (is.character(t_dimnames)){
    stop("The t_dimnames argument must be numeric.")
  }
  no_sp <- nrow(object@species_params)
  species_names <- dimnames(object@psi)$sp
  no_w <- length(object@w)
  w_names <- dimnames(object@psi)$w
  t_dimnames_n <- c(t_dimnames[1] - (t_dimnames[2]-t_dimnames[1]),t_dimnames) # N is 1 bigger because it holds the initial population
  t_dim_n <- length(t_dimnames_n) 
  t_dim_effort <- length(t_dimnames)
  array_n <- array(NA, dim = c(t_dim_n, no_sp, no_w), dimnames = list(time = t_dimnames_n, sp = species_names, w = w_names))
  
  no_gears <- dim(object@selectivity)[1]
  gear_names <- dimnames(object@selectivity)$gear
  array_effort <- array(NA, dim = c(t_dim_effort, no_gears), dimnames = list(time = t_dimnames, gear = gear_names))
  
  no_w_full <- length(object@w_full)
  w_full_names <- names(object@rr_pp)
  array_n_pp <- array(NA, dim = c(t_dim_n, no_w_full), dimnames = list(time=t_dimnames_n, w = w_full_names))
  
  sim <- new('MizerSim',
             n = array_n, 
             effort = array_effort,
             n_pp = array_n_pp,
             params = object)
  return(sim)
}



#----------------------------------------------------------------
# Phi prey
#----------------------------------------------------------------
#' getPhiPrey method for the size based model
#'
#' Calculates the amount of food exposed to each predator by predator size.
#' This method is used by the \code{\link{project}} method for performing simulations.
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A two dimensional array (predator species x predator size) 
#' @seealso \code{\link{project}}
#' @export
#' @docType methods
#' @rdname getPhiPrey-methods
#' @aliases getPhiPrey-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPhiPrey(params,n,n_pp)
#' }

setGeneric('getPhiPrey', function(object, n, n_pp, opt = T,...)
  standardGeneric('getPhiPrey'))

#' @rdname getPhiPrey-methods
#' @aliases getPhiPrey,MizerParams,matrix,numeric-method
setMethod('getPhiPrey', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
          function(object, n, n_pp, opt = F, ...){
            #        cat("In getPhiPrey\n")
            # Check n dims
            if(dim(n)[1] != dim(object@interaction)[1])
              stop("n does not have the right number of species (first dimension)")
            if(dim(n)[2] != length(object@w))
              stop("n does not have the right number of size groups (second dimension)")
            if(length(n_pp) != length(object@w_full))
              stop("n_pp does not have the right number of size groups")
            
            if (opt){
              n_eff_prey <- object@interaction[1,] %*% n * object@w *object@dw
              idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
              
              phi_prey_species_1 <- rowSums(sweep(object@pred_kernel[1,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*"),dims=2)
              phi_prey_species <- matrix(rep(phi_prey_species_1,each=dim(n)[1]),nrow=dim(n)[1], dimnames = dimnames(n))
              names(dimnames(phi_prey_species)) <- c("sp","w_pred")
              
              phi_prey_background_1 <- rowSums(sweep(object@pred_kernel[1,,],2,object@dw_full*object@w_full*n_pp,"*"))
              phi_prey_background <- matrix(rep(phi_prey_background_1,each=dim(n)[1]),nrow=dim(n)[1], dimnames = dimnames(n))
              names(dimnames(phi_prey_background)) <- c("sp","w_pred")
            } else{
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
            }
            return(phi_prey_species+phi_prey_background)
          })

# my own phi prey function to get the details of fish and background
setGeneric('getPhi', function(object, n, n_pp,...)
  standardGeneric('getPhi'))

#' @rdname getPhiPrey-methods
#' @aliases getPhiPrey,MizerParams,matrix,numeric-method
setMethod('getPhi', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
          function(object, n, n_pp, ...){
            #        cat("In getPhiPrey\n")
            # Check n dims
            if(dim(n)[1] != dim(object@interaction)[1])
              stop("n does not have the right number of species (first dimension)")
            if(dim(n)[2] != length(object@w))
              stop("n does not have the right number of size groups (second dimension)")
            if(length(n_pp) != length(object@w_full))
              stop("n_pp does not have the right number of size groups")
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
            #return(phi_prey_species+phi_prey_background)
            return(list(phi_prey_species,phi_prey_background))
          })
#----------------------------------------------------------------
# Feeding Level
#----------------------------------------------------------------
# The amount of food consumed by a predator, by each predator size

#' getFeedingLevel method for the size based model
#'
#' Calculates the amount of food consumed by a predator by predator size based on food availability, search volume and maximum intake.
#' This method is used by the \code{\link{project}} method for performing simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object
#' @param n A matrix of species abundance (species x size). Only used if \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the background abundance by size. Only used if \code{object} argument is of type \code{MizerParams}.
#' @param phi_prey The PhiPrey matrix (optional) of dimension no. species x no. size bins. If not passed in, it is calculated internally using the \code{getPhiPrey()} method. Only used if \code{object} argument is of type \code{MizerParams}.
#' @param time_range Subset the returned fishing mortalities by time. The time range is either a vector of values, a vector of min and max time, or a single value. Default is the whole time range. Only used if the \code{object} argument is of type \code{MizerSim}.
#' @param drop should extra dimensions of length 1 in the output be dropped, simplifying the output. Defaults to TRUE  
#'
#' @note
#' If a \code{MizerParams} object is passed in, the method returns a two dimensional array (predator species x predator size) based on the abundances also passed in.
#' If a \code{MizerSim} object is passed in, the method returns a three dimensional array (time step x predator species x predator size) with the feeding level calculated at every time step in the simulation.
#' @seealso \code{\link{getPhiPrey}}
#' @export
#' @docType methods
#' @rdname getFeedingLevel-methods
#' @aliases getFeedingLevel-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' fl <- getFeedingLevel(params,n,n_pp)
#' # Get the feeding level at all saved time steps
#' fl <- getFeedingLevel(sim)
#' # Get the feeding level for time 15 - 20
#' fl <- getFeedingLevel(sim, time_range = c(15,20))
#' }
setGeneric('getFeedingLevel', function(object, n, n_pp, phi_prey, ...)
  standardGeneric('getFeedingLevel'))

#' @rdname getFeedingLevel-methods
#' @aliases getFeedingLevel,MizerParams,matrix,numeric,matrix-method
setMethod('getFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric', phi_prey='matrix'),
          function(object, n, n_pp, phi_prey, ...){
            # Check dims of phi_prey
            if (!all(dim(phi_prey) == c(nrow(object@species_params),length(object@w)))){
              stop("phi_prey argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            # encountered food = available food * search volume
            encount <- object@search_vol * phi_prey
            # calculate feeding level
            f <- encount/(encount + object@intake_max)
            return(f)
          })
#' @rdname getFeedingLevel-methods
#' @aliases getFeedingLevel,MizerParams,matrix,numeric,missing-method
setMethod('getFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric', phi_prey='missing'),
          function(object, n, n_pp, ...){
            phi_prey <- getPhiPrey(object, n=n, n_pp=n_pp)
            # encountered food = available food * search volume
            #encount <- object@search_vol * phi_prey
            # calculate feeding level
            #f <- encount/(encount + object@intake_max)
            f <- getFeedingLevel(object=object, n=n, n_pp=n_pp, phi_prey=phi_prey)
            return(f)
          })

#' @rdname getFeedingLevel-methods
#' @aliases getFeedingLevel,MizerSim,missing,missing,missing-method
setMethod('getFeedingLevel', signature(object='MizerSim', n = 'missing', n_pp='missing', phi_prey='missing'),
          function(object, time_range=dimnames(object@n)$time, drop=FALSE, ...){
            time_elements <- get_time_elements(object,time_range)
            feed_time <- aaply(which(time_elements), 1, function(x){
              # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
              n <- array(object@n[x,,],dim=dim(object@n)[2:3])
              dimnames(n) <- dimnames(object@n)[2:3]
              feed <- getFeedingLevel(object@params, n=n, n_pp = object@n_pp[x,])
              return(feed)}, .drop=drop)
            return(feed_time)
          })


#----------------------
# Critical feeding level
#------------------------

setGeneric('getCFeedingLevel', function(object, n, n_pp, phi_prey, ...)
  standardGeneric('getCFeedingLevel'))

#' @rdname getCFeedingLevel-methods
#' @aliases getFeedingLevel,MizerParams,matrix,numeric,matrix-method
setMethod('getCFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric', phi_prey='matrix'),
          function(object, n, n_pp, phi_prey, ...){
            # Check dims of phi_prey
            if (!all(dim(phi_prey) == c(nrow(object@species_params),length(object@w)))){
              stop("phi_prey argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            f = sweep(object@std_metab,2,(object@intake_max * object@species_params$alpha)[1,],"/")
            return(f)
          })

#' @rdname getCFeedingLevel-methods
#' @aliases getFeedingLevel,MizerParams,matrix,numeric,missing-method
setMethod('getCFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric', phi_prey='missing'),
          function(object, n, n_pp, ...){
            phi_prey <- getPhiPrey(object, n=n, n_pp=n_pp)
            # encountered food = available food * search volume
            #encount <- object@search_vol * phi_prey
            # calculate feeding level
            #f <- encount/(encount + object@intake_max)
            f <- getCFeedingLevel(object=object, n=n, n_pp=n_pp, phi_prey=phi_prey)
            return(f)
          })

#' @rdname getCFeedingLevel-methods
#' @aliases getFeedingLevel,MizerSim,missing,missing,missing-method
setMethod('getCFeedingLevel', signature(object='MizerSim', n = 'missing', n_pp='missing', phi_prey='missing'),
          function(object, time_range=dimnames(object@n)$time, drop=FALSE, ...){
        
            time_elements <- get_time_elements(object,time_range)
            feed_time <- aaply(which(time_elements), 1, function(x){
              # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
              n <- array(object@n[x,,],dim=dim(object@n)[2:3])
              dimnames(n) <- dimnames(object@n)[2:3]
              feed <- getCFeedingLevel(object@params, n=n, n_pp = object@n_pp[x,])
              return(feed)}, .drop=drop)
            return(feed_time)
          })

#----------------------------------------------------------------
# Predation rate
#----------------------------------------------------------------
#' getPredRate method for the size based model
#'
#' Calculates the predation rate of each predator species at size on prey size.
#' This method is used by the \code{\link{project}} method for performing simulations. In the simulations, it is combined with the interaction matrix (see \code{\link{MizerParams}}) to calculate the realised predation mortality (see \code{\link{getM2}}).
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size no. species x no. size bins. If not supplied, is calculated internally using the \code{getFeedingLevel()} method.
#'
#' @return A three dimensional array (predator species x predator size x prey size) 
#' @export
#' @seealso \code{\link{project}}, \code{\link{getM2}}, \code{\link{getFeedingLevel}} and \code{\link{MizerParams}}
#' @docType methods
#' @rdname getPredRate-methods
#' @aliases getPredRate-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredRate(params,n,n_pp)
#' }
setGeneric('getPredRate', function(object, n, n_pp, feeding_level,...)
  standardGeneric('getPredRate'))

#' @rdname getPredRate-methods
#' @aliases getPredRate,MizerParams,matrix,numeric,matrix-method
setMethod('getPredRate', signature(object='MizerParams', n = 'matrix', n_pp='numeric', feeding_level = 'matrix'),
          function(object, n, n_pp, feeding_level, ...){
            if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
              stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            n_total_in_size_bins <- sweep(n, 2, object@dw, '*')
            pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*")
            return(pred_rate)
          })

#' @rdname getPredRate-methods
#' @aliases getPredRate,MizerParams,matrix,numeric,missing-method
setMethod('getPredRate', signature(object='MizerParams', n = 'matrix', n_pp='numeric', feeding_level = 'missing'),
          function(object, n, n_pp, ...){
            n_total_in_size_bins <- sweep(n, 2, object@dw, '*')
            feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
            #pred_rate <- sweep(object@pred_kernel,c(1,2),(1-f)*object@search_vol*n_total_in_size_bins,"*")
            pred_rate <- getPredRate(object=object, n=n, n_pp=n_pp, feeding_level = feeding_level)
            return(pred_rate)
          })


#----------------------------------------------------------------
# M2, mortality by predation
#----------------------------------------------------------------
# This uses the predation rate which is also used in M2background
# Too much overlap? Inefficient? Same thing is calculated twice

#' getM2 method for the size based model
#'
#' Calculates the total predation mortality on each prey species by prey size.
#' This method is used by the \code{\link{project}} method for performing simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object.
#' @param n A matrix of species abundance (species x size). Only used if \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the background abundance by size. Only used if \code{object} argument is of type \code{MizerParams}.
#' @param pred_rate An array of predation rates of dimension no. sp x no. community size bins x no. of size bins in whole spectra (i.e. community + background, the w_full slot). The array is optional. If it is not provided it is calculated by the \code{getPredRate()} method.
#' @param time_range Subset the returned fishing mortalities by time. The time range is either a vector of values, a vector of min and max time, or a single value. Default is the whole time range. Only used if the \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should dimensions of length 1 in the output be dropped, simplifying the output. Defaults to TRUE  
#'
#' @note
#' If a \code{MizerParams} object is passed in, the method returns a two dimensional array (prey species x prey size) based on the abundances also passed in.
#' If a \code{MizerSim} object is passed in, the method returns a three dimensional array (time step x prey species x prey size) with the predation mortality calculated at every time step in the simulation.
#' @seealso \code{\link{getPredRate}} and \code{\link{project}}.
#' @export
#' @docType methods
#' @rdname getM2-methods
#' @aliases getM2-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get M2 at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2(params,n,n_pp)
#' # Get M2 at all saved time steps
#' getM2(sim)
#' # Get M2 over the time 15 - 20
#' getM2(sim, time_range = c(15,20))
#' }
setGeneric('getM2', function(object, n, n_pp, pred_rate,...)
  standardGeneric('getM2'))

#' @rdname getM2-methods
#' @aliases getM2,MizerParams,matrix,numeric,array-method
setMethod('getM2', signature(object='MizerParams', n = 'matrix', n_pp='numeric', pred_rate = 'array'),
          function(object, n, n_pp, pred_rate, ...){
            if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
              stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
            }
            # get the element numbers that are just species
            idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
            # Interaction is predator x prey so need to transpose so it is prey x pred
            # Sum pred_kernel over predator sizes to give total predation rate of each predator on each prey size
            m2 <- t(object@interaction) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]
            return(m2)
          })
#' @rdname getM2-methods
#' @aliases getM2,MizerParams,matrix,numeric,missing-method
setMethod('getM2', signature(object='MizerParams', n = 'matrix', n_pp='numeric', pred_rate = 'missing'),
          function(object, n, n_pp, ...){
            pred_rate <- getPredRate(object,n=n,n_pp=n_pp)
            m2 <- getM2(object,n=n,n_pp=n_pp, pred_rate=pred_rate)
            return(m2)
          })
#' @rdname getM2-methods
#' @aliases getM2,MizerSim,missing,missing,missing-method
setMethod('getM2', signature(object='MizerSim', n = 'missing', n_pp='missing', pred_rate = 'missing'),
          function(object, time_range=dimnames(object@n)$time, drop=TRUE, ...){
            time_elements <- get_time_elements(object,time_range)
            m2_time <- aaply(which(time_elements), 1, function(x){
              n <- array(object@n[x,,],dim=dim(object@n)[2:3])
              dimnames(n) <- dimnames(object@n)[2:3]
              m2 <- getM2(object@params, n=n, n_pp = object@n_pp[x,])
              return(m2)}, .drop=drop)
            return(m2_time)
          })


#----------------------------------------------------------------
# Starvation mortality
#----------------------------------------------------------------

setGeneric('getSmort', function(object, n, n_pp, e,...)
  standardGeneric('getSmort'))

#' @rdname getSmort-methods
#' @aliases getSmort,MizerParams,matrix,numeric,array-method
setMethod('getSmort', signature(object='MizerParams', n = 'matrix', n_pp='numeric', e = 'matrix'),
          function(object, n, n_pp, e, ...){
            mu_S <- e - object@std_metab - object@activity
            mu_S[mu_S<0]<- t(apply(mu_S,1,function(x,y = object@dw){x[which(x<0)]/y[which(x<0)]}))
            mu_S[mu_S>0] <- 0
            return(mu_S)
          })
#' @rdname getSmort-methods
#' @aliases getSmort,MizerParams,matrix,numeric,missing-method
setMethod('getSmort', signature(object='MizerParams', n = 'matrix', n_pp='numeric', e = 'missing'),
          function(object, n, n_pp, ...){
            e <- getE(object,n=n,n_pp=n_pp)
            mu_S <- getSmort(object,n=n,n_pp=n_pp, e=e)
            return(mu_S)
          })
#' @rdname getSmort-methods
#' @aliases getSmort,MizerSim,missing,missing,missing-method
setMethod('getSmort', signature(object='MizerSim', n = 'missing', n_pp='missing', e = 'missing'),
          function(object, time_range=dimnames(object@n)$time, drop=TRUE, ...){
            time_elements <- get_time_elements(object,time_range)
            e_time <- aaply(which(time_elements), 1, function(x){
              n <- array(object@n[x,,],dim=dim(object@n)[2:3])
              dimnames(n) <- dimnames(object@n)[2:3]
              e <- getE(object@params, n=n, n_pp = object@n_pp[x,])
              return(e)}, .drop=drop)
            return(e_time)
          })


#----------------------------------------------------------------
# Senescence mortality
#----------------------------------------------------------------

setGeneric('getOmort', function(object,...)
  standardGeneric('getOmort'))

#' @rdname getOmort-methods
#' @aliases getOmort,MizerParams,matrix,numeric,array-method
setMethod('getOmort', signature(object='MizerParams'),
          function(object, ...){
              # exponential parameters
              m_min = 0.1
              m_max = 1
              m1 = 1
              
              Omort <- matrix(0,ncol=length(object@w),nrow=dim(object@species_params)[1],dimnames = list(object@species_params$ecotype,object@w))
              for (iSpecies in 1:dim(Omort)[1])
              {
                size = (object@w[which(object@w >= object@species_params$w_mat[iSpecies])[1]:which(object@w >= object@species_params$w_inf[iSpecies])[1]])
                sizeN = size/size[length(size)]
                m_rate <- m_min + (m_max-m_min)* exp(-m1* 1/sizeN)
                a <- length(c(rep(0,(which(size[1]==object@w)-1)), m_rate, rep(0,length(object@w)-(which(size[length(size)]==object@w)))))
                if (a != 100) print(a)
                Omort[iSpecies,] <- c(rep(0,(which(size[1]==object@w)-1)), m_rate, rep(0,length(object@w)-(which(size[length(size)]==object@w))))
              }
            return(Omort)
          })

#' @rdname getOmort-methods
#' @aliases getOmort,MizerSim,missing,missing,missing-method
setMethod('getOmort', signature(object='MizerSim'),
          function(object, ...){

            Omort <- getOmort(object@params)
            return(Omort)
          })

#----------------------------------------------------------------
# M2 background
#----------------------------------------------------------------
#'
#' Calculates the predation mortality on the background spectrum by prey size. Used by the \code{project} method for running size based simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param pred_rate An array of predation rates of dimension no. sp x no. community size bins x no. of size bins in whole spectra (i.e. community + background, the w_full slot). The array is optional. If it is not provided it is calculated by the \code{getPredRate()} method.
#'
#' @return A vector of predation mortalities by background prey size.
#' @seealso \code{\link{project}} and \code{\link{getM2}}.
#' @export
#' @docType methods
#' @rdname getM2Background-methods
#' @aliases getM2Background-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get M2 of the background spectrum at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2Background(params,n,n_pp)
#' }
setGeneric('getM2Background', function(object, n, n_pp, pred_rate,...)
  standardGeneric('getM2Background'))

#' @rdname getM2Background-methods
#' @aliases getM2Background,MizerParams,matrix,numeric,array-method
setMethod('getM2Background', signature(object='MizerParams', n = 'matrix', n_pp='numeric', pred_rate='array'),
          function(object, n, n_pp, pred_rate, ...){
            if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
              stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
            }
            M2background <- colSums(pred_rate,dims=2)
            return(M2background)
          })

#' @rdname getM2Background-methods
#' @aliases getM2Background,MizerParams,matrix,numeric,missing-method
setMethod('getM2Background', signature(object='MizerParams', n = 'matrix', n_pp='numeric',  pred_rate='missing'),
          function(object, n, n_pp, ...){
            pred_rate <- getPredRate(object,n=n,n_pp=n_pp)
            M2background <- getM2Background(object, n=n, n_pp=n_pp, pred_rate=pred_rate)
            return(M2background)
          })

#----------------------------------------------------------------
# Fisheries mortality gear
#----------------------------------------------------------------
# getFMortGear
#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality by gear, species and size at each time step in the \code{effort} argument. Used by the \code{project} method to perform simulations.
#'
#' @param object A \code{MizerParams} object or a \code{MizerSim} object.
#' @param effort The effort of each fishing gear. Only needed if the object argument is of class \code{MizerParams}. See notes below. 
#' @param time_range Subset the returned fishing mortalities by time. The time range is either a vector of values, a vector of min and max time, or a single value. Default is the whole time range. Only used if the \code{object} argument is of type \code{MizerSim}.
#'
#' @return An array. If the effort argument has a time dimension, or a \code{MizerSim} is passed in, the output array has four dimensions (time x gear x species x size). If the effort argument does not have a time dimension (i.e. it is a vector or a single numeric), the output array has three dimensions (gear x species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The \code{effort} argument is only used if a \code{MizerParams} object is passed in. The \code{effort} argument can be a two dimensional array (time x gear), a vector of length equal to the number of gears (each gear has a different effort that is constant in time), or a single numeric value (each gear has the same effort that is constant in time). The order of gears in the \code{effort} argument must be the same the same as in the \code{MizerParams} object.
#'
#' If the object argument is of class \code{MizerSim} then the effort slot of the \code{MizerSim} object is used and the \code{effort} argument is not used.
#' @export
#' @docType methods
#' @rdname getFMortGear-methods
#' @aliases getFMortGear-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Get the fishing mortality when effort is constant
#' # for all gears and time:
#' getFMortGear(params, effort = 1)
#' # Get the fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMortGear(params, effort = c(0.5,1,1.5,0.75))
#' # Get the fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMortGear(params, effort=effort)
#' # Get the fishing mortality using the effort already held in a MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMortGear(sim)
#' getFMortGear(sim, time_range=c(10,20))
#' }
setGeneric('getFMortGear', function(object, effort, ...)
  standardGeneric('getFMortGear'))

#' @rdname getFMortGear-methods
#' @aliases getFMortGear,MizerParams,numeric-method
# Effort is a single value or a numeric vector.
# Effort has no time time dimension
setMethod('getFMortGear', signature(object='MizerParams', effort = 'numeric'),
          function(object, effort, ...){
            no_gear <- dim(object@catchability)[1]
            # If a single value, just repeat it for all gears
            if(length(effort) == 1)
              effort <- rep(effort, no_gear)
            if (length(effort) != no_gear)
              stop("Effort must be a single value or a vector as long as the number of gears\n")
            # turn to array and call next method
            effort <- array(effort,dim=c(1,no_gear))
            fmort_gear <- getFMortGear(object,effort)
            # fmort_gear is 4D, and first D is time with length 1
            # Drop time dimension - bit annoying because we want to keep the other dims even if they have length 1
            out <- array(fmort_gear, dim=dim(fmort_gear)[2:4])
            dimnames(out) <- dimnames(fmort_gear)[2:4]
            return(out)
          }
)

#' @rdname getFMortGear-methods
#' @aliases getFMortGear,MizerParams,matrix-method
# Always returns a 4D array: time x gear x species x size
setMethod('getFMortGear', signature(object='MizerParams', effort = 'matrix'),
          function(object, effort, ...){
            no_gear <- dim(object@catchability)[1]
            if (dim(effort)[2] != no_gear)
              stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
            # F = sel * q * effort
            sel_q <- sweep(object@selectivity, c(1,2), object@catchability, "*")
            # Kinda nasty! ends up with 4D array 
            fmort_gear <- aaply(effort, 1, function(x,sel_q) sweep(sel_q, c(1), x, "*"), sel_q=sel_q, .drop=FALSE)
            return(fmort_gear)
          }
)

# Returns the fishing mortality: time * gear * species * size
#' @rdname getFMortGear-methods
#' @aliases getFMortGear,MizerSim,missing-method
setMethod('getFMortGear', signature(object='MizerSim', effort='missing'),
          function(object,effort, time_range=dimnames(object@effort)$time, ...){
            time_elements <- get_time_elements(object,time_range, slot="effort")
            f_mort_gear <- getFMortGear(object@params, object@effort, ...)
            return(f_mort_gear[time_elements,,,,drop=FALSE])
          })


#----------------------------------------------------------------
# Fisheries mortality
#----------------------------------------------------------------
# Total fishing mortality from all gears
# species x size and maybe also by time if effort is time based

#' Get the total fishing mortality from all fishing gears by time, species and size
#'
#' Calculates the fishing mortality from all gears by species and size at each time step in the \code{effort} argument.
#' The total fishing mortality is just the sum of the fishing mortalities imposed by each gear.
#'
#' @param object A \code{MizerParams} object or a \code{MizerSim} object
#' @param effort The effort of each fishing gear. Only needed if the object argument is of class \code{MizerParams}. See notes below. 
#' @param time_range Subset the returned fishing mortalities by time. The time range is either a vector of values, a vector of min and max time, or a single value. Default is the whole time range. Only used if the \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should dimensions of length 1 be dropped, e.g. if your community only has one species it might make presentation of results easier. Default is TRUE
#'
#' @return An array. If the effort argument has a time dimension, or object is of class \code{MizerSim}, the output array has three dimensions (time x species x size). If the effort argument does not have a time dimension, the output array has two dimensions (species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The \code{effort} argument is only used if a \code{MizerParams} object is passed in. The \code{effort} argument can be a two dimensional array (time x gear), a vector of length equal to the number of gears (each gear has a different effort that is constant in time), or a single numeric value (each gear has the same effort that is constant in time). The order of gears in the \code{effort} argument must be the same the same as in the \code{MizerParams} object.
#'
#' If the object argument is of class \code{MizerSim} then the effort slot of the \code{MizerSim} object is used and the \code{effort} argument is not used.
#' @export
#' @docType methods
#' @rdname getFMort-methods
#' @aliases getFMort-method
#' @seealso \code{getFMortGear}, \code{project}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Get the total fishing mortality when effort is constant for all gears and time:
#' getFMort(params, effort = 1)
#' # Get the total fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMort(params, effort = c(0.5,1,1.5,0.75))
#' # Get the total fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMort(params, effort=effort)
#' # Get the total fishing mortality using the effort already held in a MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMort(sim)
#' getFMort(sim, time_range = c(10,20))
#' }
setGeneric('getFMort', function(object, effort, ...)
  standardGeneric('getFMort'))

#' @rdname getFMort-methods
#' @aliases getFMort,MizerParams,numeric-method
setMethod('getFMort', signature(object='MizerParams', effort='numeric'),
          function(object, effort, ...){
            fMortGear <- getFMortGear(object, effort, ...)
            fMort <- apply(fMortGear, c(2,3), sum)
            return(fMort)
          })

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
          function(object, effort, time_range=dimnames(object@effort)[[1]], drop=TRUE, ...){

 print(time_range)
            time_elements <- get_time_elements(object,time_range, slot="effort")
  
            fMort <- getFMort(object@params, object@effort, ...)
   
            return(fMort[time_elements,,,drop=drop])
          })


#----------------------------------------------------------------
# Total mortality Z
#----------------------------------------------------------------
#' getZ method for the size based model
#'
#' Calculates the total mortality on each species by size from predation mortality (M2), background mortality (M) and fishing mortality for a single time step.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param effort A numeric vector of the effort by gear or a single numeric effort value which is used for all gears.
#' @param m2 A two dimensional array of predation mortality (optional). Has dimensions no. sp x no. size bins in the community. If not supplied is calculated using the \code{getM2()} method.
#'
#' @return A two dimensional array (prey species x prey size). 
#' @export
#' @seealso \code{\link{getM2}}, \code{\link{getFMort}}
#' @docType methods
#' @rdname getZ-methods
#' @aliases getZ-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getZ(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
setGeneric('getZ', function(object, n, n_pp, effort, m2, ...)
  standardGeneric('getZ'))

#' @rdname getZ-methods
#' @aliases getZ,MizerParams,matrix,numeric,numeric,matrix-method
setMethod('getZ', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', effort='numeric', m2 = 'matrix'),
          function(object, n, n_pp, effort, m2){
            if (!all(dim(m2) == c(nrow(object@species_params),length(object@w)))){
              stop("m2 argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            f_mort <- getFMort(object, effort = effort)
            mu_S <- getSmort(object, n=n, n_pp = n_pp)
            mu_O <- getOmort(object)
            z = sweep(m2- mu_S + mu_O + f_mort,1,object@species_params$z0,"+")
            #z = sweep(m2- mu_S + f_mort,1,object@species_params$z0,"+")
            return(z)
          })

#' @rdname getZ-methods
#' @aliases getZ,MizerParams,matrix,numeric,numeric,missing-method
setMethod('getZ', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', effort='numeric', m2 = 'missing'),
          function(object, n, n_pp, effort){
            m2 <- getM2(object, n=n, n_pp=n_pp)
            z <- getZ(object, n=n, n_pp=n_pp, effort=effort, m2=m2)
            return(z)
          })


#----------------------------------------------------------------
# Energy after metabolism and movement
#----------------------------------------------------------------
#' getEReproAndGrowth method for the size based model
#'
#' Calculates the energy available by species and size for reproduction and growth after metabolism and movement have been accounted for.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size no. species x no. size bins. If not supplied, is calculated internally using the \code{getFeedingLevel()} method.
#'
#' @return A two dimensional array (species x size) 
#' @export
#' @docType methods
#' @rdname getEReproAndGrowth-methods
#' @aliases getEReproAndGrowth-method
#' @seealso \code{\link{project}} and \code{\link{getFeedingLevel}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEReproAndGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getEReproAndGrowth', function(object, n, n_pp, feeding_level, ...)
  standardGeneric('getEReproAndGrowth'))

#' @rdname getEReproAndGrowth-methods
#' @aliases getEReproAndGrowth,MizerParams,matrix,numeric,matrix-method
setMethod('getEReproAndGrowth', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', feeding_level='matrix'),
          function(object, n, n_pp, feeding_level){
            if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
              stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            # assimilated intake
            e <- sweep(feeding_level * object@intake_max,1,object@species_params$alpha,"*")
            # Subtract basal metabolism and activity 
            e <- e - object@std_metab - object@activity
            e[e<0] <- 0 # Do not allow negative growth
            return(e)
          })

#' @rdname getEReproAndGrowth-methods
#' @aliases getEReproAndGrowth,MizerParams,matrix,numeric,missing-method
setMethod('getEReproAndGrowth', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', feeding_level='missing'),
          function(object, n, n_pp){
            feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
            e <- getEReproAndGrowth(object, n=n, n_pp=n_pp, feeding_level=feeding_level)
            return(e)
          })


#my own energy function, without the metabolism

setGeneric('getE', function(object, n, n_pp, feeding_level, ...)
  standardGeneric('getE'))

#' @rdname getE
#' @aliases getE,MizerParams,matrix,numeric,matrix-method
setMethod('getE', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', feeding_level='matrix'),
          function(object, n, n_pp, feeding_level){
            if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
              stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            # assimilated intake
            e <- sweep(feeding_level * object@intake_max,1,object@species_params$alpha,"*")
            e[e<0] <- 0 # Do not allow negative growth
            return(e)
          })

#' @rdname getE
#' @aliases getE,MizerParams,matrix,numeric,missing-method
setMethod('getE', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', feeding_level='missing'),
          function(object, n, n_pp){
            feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
            e <- getE(object, n=n, n_pp=n_pp, feeding_level=feeding_level)
            return(e)
          })

#----------------------------------------------------------------
# Energy left for reproduction
#----------------------------------------------------------------
# assimilated food intake, less metabolism and activity, split between reproduction and growth

#' getESpawning method for the size based model
#'
#' Calculates the energy available by species and size for reproduction after metabolism and movement have been accounted for.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e The energy available for reproduction and growth (optional). A matrix of size no. species x no. size bins. If not supplied, is calculated internally using the \code{getEReproAndGrowth()} method.
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getESpawning-methods
#' @aliases getESpawning-method
#' @seealso \code{\link{project}} and \code{\link{getEReproAndGrowth}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getESpawning(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getESpawning', function(object, n, n_pp, e, ...)
  standardGeneric('getESpawning'))

#' @rdname getESpawning-methods
#' @aliases getESpawning,MizerParams,matrix,numeric,matrix-method
setMethod('getESpawning', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', e = 'matrix'),
          function(object, n, n_pp, e){
            if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
              stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            e_spawning <- object@psi * e 
            return(e_spawning)
          }
)
#' @rdname getESpawning-methods
#' @aliases getESpawning,MizerParams,matrix,numeric,missing-method
setMethod('getESpawning', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', e = 'missing'),
          function(object, n, n_pp){
            e <- getEReproAndGrowth(object,n=n,n_pp=n_pp)
            e_spawning <- getESpawning(object, n=n, n_pp=n_pp, e=e)
            return(e_spawning)
          })

#----------------------------------------------------------------
# Energy left for growth
#----------------------------------------------------------------
#'
#' Calculates the energy available by species and size for growth after metabolism and movement have been accounted for.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e The energy available for reproduction and growth (optional, although if specified, e_spawning must also be specified). A matrix of size no. species x no. size bins. If not supplied, is calculated internally using the \code{getEReproAndGrowth()} method.
#' @param e_spawning The energy available for spawning (optional, although if specified, e must also be specified). A matrix of size no. species x no. size bins. If not supplied, is calculated internally using the \code{getESpawning()} method.
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getEGrowth-methods
#' @aliases getEGrowth-method
#' @seealso \code{\link{project}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getEGrowth', function(object, n, n_pp, e_spawning, e, ...)
  standardGeneric('getEGrowth'))

#' @rdname getEGrowth-methods
#' @aliases getEGrowth,MizerParams,matrix,numeric,matrix,matrix-method
setMethod('getEGrowth', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', e_spawning='matrix', e='matrix'),
          function(object, n, n_pp, e_spawning, e){
            if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
              stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
              stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            # Assimilated intake less activity and metabolism
            # energy for growth is intake - energy for growth
            e_growth <- e - e_spawning
            return(e_growth)
          })
#' @rdname getEGrowth-methods
#' @aliases getEGrowth,MizerParams,matrix,numeric,missing,missing-method
setMethod('getEGrowth', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', e_spawning='missing', e='missing'),
          function(object, n, n_pp){
            # Assimilated intake less activity and metabolism
            e <- getEReproAndGrowth(object,n=n,n_pp=n_pp)
            e_spawning <- getESpawning(object,n=n,n_pp=n_pp)
            # energy for growth is intake - energy for growth
            e_growth <- getEGrowth(object, n=n, n_pp=n_pp, e_spawning=e_spawning, e=e)
            return(e_growth)
          })

#----------------------------------------------------------------
# density independent recruitment (total egg production)
#----------------------------------------------------------------
#' getRDI method for the size based model
#'
#' Calculates the density independent recruitment (total egg production) before density dependence, by species.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e_spawning The energy available for spawning (optional). A matrix of size no. species x no. size bins. If not supplied, is calculated internally using the \code{getESpawning()} method.
#' @param sex_ratio Proportion of the population that is female. Default value is 0.5.
#'
#' @return A numeric vector the length of the number of species 
#' @export
#' @docType methods
#' @rdname getRDI-methods
#' @aliases getRDI-method
#' @seealso \code{\link{project}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the recruitment at a particular time step
#' getRDI(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getRDI', function(object, n, n_pp, e_spawning, ...)
  standardGeneric('getRDI'))

#' @rdname getRDI-methods
#' @aliases getRDI,MizerParams,matrix,numeric,matrix-method
setMethod('getRDI', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', e_spawning='matrix'),
          function(object, n, n_pp, e_spawning, sex_ratio = 0.5){
            if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
              stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            # Should we put this in the class as part of species_params?
            # Index of the smallest size class for each species
            #w0_idx <- as.vector(tapply(object@species_params$w_min,1:length(object@species_params$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=params@w))
            e_spawning_pop <- (e_spawning*n) %*% object@dw
            rdi <- sex_ratio*(e_spawning_pop * object@species_params$erepro)/object@w[object@species_params$w_min_idx] 
            return(rdi)
          })
#' @rdname getRDI-methods
#' @aliases getRDI,MizerParams,matrix,numeric,missing-method
setMethod('getRDI', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', e_spawning='missing'),
          function(object, n, n_pp, sex_ratio = 0.5){
            # Should we put this in the class as part of species_params?
            # Index of the smallest size class for each species
            #w0_idx <- as.vector(tapply(object@species_params$w_min,1:length(object@species_params$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=params@w))
            e_spawning <- getESpawning(object, n=n, n_pp=n_pp)
            rdi <- getRDI(object, n=n, n_pp=n_pp, e_spawning=e_spawning, sex_ratio=sex_ratio)
            return(rdi)
          })

#----------------------------------------------------------------
# density dependent recruitment (total egg production)
#----------------------------------------------------------------
#' getRDD method for the size based model
#'
#' Calculates the density dependent recruitment (total egg production) for each species.
#' This is the flux entering the smallest size class of each species. 
#' The density dependent recruiment is the density independent recruitment after it has been put through the density dependent stock-recruitment relationship function. 
#' This method is used by the \code{project} method for performing simulations.
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param rdi A matrix of density independent recruitment (optional) with dimensions no. sp x 1. If not specified rdi is calculated internally using the \code{getRDI()} method.
#' @param sex_ratio Proportion of the population that is female. Default value is 0.5
#'
#' @return A numeric vector the length of the number of species. 
#' @export
#' @docType methods
#' @rdname getRDD-methods
#' @aliases getRDD-method
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getRDD(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getRDD', function(object, n, n_pp, rdi, ...)
  standardGeneric('getRDD'))

#' @rdname getRDD-methods
#' @aliases getRDD,MizerParams,matrix,numeric,matrix-method
setMethod('getRDD', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', rdi='matrix'),
          function(object, n, n_pp, rdi, sex_ratio = 0.5){
            if (!all(dim(rdi) == c(nrow(object@species_params),1))){
              stop("rdi argument must have dimensions: no. species (",nrow(object@species_params),") x 1")
            }
            rdd <- object@srr(rdi = rdi, species_params = object@species_params)
            return(rdd)
          })

#' @rdname getRDD-methods
#' @aliases getRDD,MizerParams,matrix,numeric,missing-method
setMethod('getRDD', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', rdi='missing'),
          function(object, n, n_pp, sex_ratio = 0.5){
            rdi <- getRDI(object, n=n, n_pp=n_pp, sex_ratio = sex_ratio)
            rdd <- getRDD(object, n=n, n_pp=n_pp, rdi=rdi, sex_ratio=sex_ratio)
            return(rdd)
          })


#----------------------------------------------------------------
# time elements
#----------------------------------------------------------------
# internal function to get the array element references of the time dimension for the time based slots of a MizerSim object
# time_range can be character or numeric
# Necessary to include a slot_name argument because the effort and abundance slots have different time dimensions
get_time_elements <- function(sim,time_range,slot_name="n"){

  if (!(slot_name %in% c("n","effort")))
    stop("'slot_name' argument should be 'n' or 'effort'")
  if (!is(sim,"MizerSim"))
    stop("First argument to get_time_elements function must be of class MizerSim")
  time_range <- range(as.numeric(time_range))

  # Check that time range is even in object

  sim_time_range <- range(as.numeric(dimnames(slot(sim,slot_name))[[1]]))

  if ((time_range[1] < sim_time_range[1]) | (time_range[2] > sim_time_range[2]))
    stop("Time range is outside the time range of the model")
  time_elements <- (as.numeric(dimnames(slot(sim,slot_name))[[1]]) >= time_range[1]) & (as.numeric(dimnames(slot(sim,slot_name))[[1]]) <= time_range[2])
  names(time_elements) <- dimnames(slot(sim,slot_name))[[1]]
  return(time_elements)
}

#----------------------------------------------------------------
# to initialise stuff
#----------------------------------------------------------------

get_initial_n<- function(params, n0_mult = NULL, a = 0.35){
  if (!is(params,"MizerParams"))
    stop("params argument must of type MizerParams")
  no_sp <- nrow(params@species_params)
  no_w <- length(params@w)
  initial_n <- array(NA, dim=c(no_sp,no_w))
  dimnames(initial_n) <- dimnames(params@intake_max)
  # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
  # Reverse calc n and q from intake_max and search_vol slots (could add get_n as method)
  n <- (log(params@intake_max[,1] / params@species_params$h) / log(params@w[1]))[1]
  q <- (log(params@search_vol[,1] / params@species_params$gamma) / log(params@w[1]))[1]
  # Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
  if (is.null(n0_mult)){
    lambda <- 2+q-n
    kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
    n0_mult <- kappa / 1000
  }
  initial_n[] <- unlist(tapply(params@w,1:no_w,function(wx,n0_mult,w_inf,a,n,q)
    n0_mult * w_inf^(2*n-q-2+a) * wx^(-n-a),
    n0_mult=n0_mult, w_inf=params@species_params$w_inf, a=a, n=n, q=q))
  #set densities at w > w_inf to 0
  initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
  # Also any densities at w < w_min set to 0
  initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
  return(initial_n)
}


