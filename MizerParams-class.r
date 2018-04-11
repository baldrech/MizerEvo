# Class outline for sbm base parameters class
# Class has members to store parameters of size based model

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS. finlay.scott@cefas.co.uk

#Naming conventions:
#S4 classes and constructors: AClass
#S4 methods: aMethod
#functions a_function


# Validity function - pretty long...
# Not documented as removed later on
valid_MizerParams <- function(object) {
  
  errors <- character()
  # grab some dims
  length_w <- length(object@w)
  length_w_full <- length(object@w_full)
  
  # Check dw and dw_full are correct length
  if(length(object@dw) != length_w){
    msg <- paste("dw is length ", length(object@dw), " and w is length ", length_w, ". These should be the same length", sep="")
    errors <- c(errors, msg)
  }
  
  if(length(object@dw_full) != length_w_full){
    msg <- paste("dw_full is length ", length(object@dw_full), " and w_full is length ", length_w_full, ". These should be the same length", sep="")
    errors <- c(errors, msg)
  }
  # Check the array dimensions are good
  # 2D arrays
  if(!all(c(
    length(dim(object@psi)),
    length(dim(object@intake_max)),
    length(dim(object@search_vol)),
    length(dim(object@activity)),
    length(dim(object@std_metab)),
    length(dim(object@interaction)),
    length(dim(object@catchability))) == 2)){
    msg <- "psi, intake_max, search_vol, activity, std_metab, interaction and catchability must all be two dimensions"
    errors <- c(errors, msg)
  }
  # 3D arrays
  if(!all(c(
    length(dim(object@pred_kernel)),
    length(dim(object@selectivity))) == 3)){
    msg <- "pred_kernel and selectivity must be three dimensions"
    errors <- c(errors, msg)
  }
  # Check number of species is equal across relevant slots
  if(!all(c(
    dim(object@psi)[1],
    dim(object@intake_max)[1],
    dim(object@search_vol)[1],
    dim(object@std_metab)[1],
    dim(object@activity)[1],
    dim(object@pred_kernel)[1],
    dim(object@selectivity)[2],
    dim(object@catchability)[2],
    dim(object@interaction)[1],
    dim(object@interaction)[2]) == 
    dim(object@species_params)[1])){
    msg <- "The number of species in the model must be consistent across the species_params, psi, intake_max, search_vol, activity, pred_kernel, interaction (dim 1), selectivity, catchability and interaction (dim 2) slots"
    errors <- c(errors, msg)
  }
  # Check number of size groups
  if(!all(c(
    dim(object@psi)[2],
    dim(object@intake_max)[2],
    dim(object@search_vol)[2],
    dim(object@activity)[2],
    dim(object@std_metab)[2],
    dim(object@pred_kernel)[2],
    dim(object@selectivity)[3]) ==
    length_w)){
    msg <- "The number of size bins in the model must be consistent across the w, psi, intake_max, search_vol, activity, pred_kernel (dim 2) and selectivity (dim 3) slots"
    errors <- c(errors, msg)
  }
  # Check number of full spectrum size groups
  if(!isTRUE(all.equal(dim(object@pred_kernel)[3],length(object@w_full)))){
    msg <- "The length of the full size spectrum in the third dimension of the pred_kernel slot must be the same length as the w_full slot"
    errors <- c(errors, msg)
  }
  # Check numbe of gears
  if(!isTRUE(all.equal(dim(object@selectivity)[1], dim(object@catchability)[1]))){
    msg <- "The number of fishing gears must be consistent across the catchability and selectivity (dim 1) slots"
    errors <- c(errors, msg)
  }
  # Check names of dimnames of arrays
  # sp dimension
  if(!all(c(
    names(dimnames(object@psi))[1],
    names(dimnames(object@intake_max))[1],
    names(dimnames(object@search_vol))[1],
    names(dimnames(object@activity))[1],
    names(dimnames(object@std_metab))[1],
    names(dimnames(object@pred_kernel))[1],
    names(dimnames(object@selectivity))[2],
    names(dimnames(object@catchability))[2]) == "sp")){
    msg <- "Name of first dimension of psi, intake_max, search_vol, std_metab, activity and pred_kernel and the second dimension of selectivity and catchability must be 'sp'"
    errors <- c(errors, msg)
  }
  #interaction dimension names
  if(names(dimnames(object@interaction))[1] != "predator"){
    msg <- "The first dimension of interaction must be called 'predator'"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@interaction))[2] != "prey"){
    msg <- "The first dimension of interaction must be called 'prey'"
    errors <- c(errors, msg)
  }
  # w dimension
  if(!all(c(
    names(dimnames(object@psi))[2],
    names(dimnames(object@intake_max))[2],
    names(dimnames(object@search_vol))[2],
    names(dimnames(object@std_metab))[2],
    names(dimnames(object@activity))[2],
    names(dimnames(object@selectivity))[3]) == "w")){
    msg <- "Name of second dimension of psi, intake_max, search_vol, std_metab, activity and third dimension of selectivity must be 'w'"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@pred_kernel))[2] != "w_pred"){
    msg <- "Name of second dimension of pred_kernel must be 'w_pred'"
    errors <- c(errors, msg)
  }
  if(names(dimnames(object@pred_kernel))[3] != "w_prey"){
    msg <- "Name of third dimension of pred_kernel must be 'w_prey'"
    errors <- c(errors, msg)
  }
  if(!all(c(
    names(dimnames(object@selectivity))[1],
    names(dimnames(object@catchability))[1]) == "gear")){
    msg <- "Name of first dimension of selectivity and catchability must be 'gear'"
    errors <- c(errors, msg)
  }
  
  # Check dimnames of species are identical
  # Bit tricky this one as I don't know of a way to compare lots of vectors at the same time. Just use == and the recycling rule
  if(!all(c(
    dimnames(object@psi)[[1]],
    dimnames(object@intake_max)[[1]],
    dimnames(object@search_vol)[[1]],
    dimnames(object@std_metab)[[1]],
    dimnames(object@activity)[[1]],
    dimnames(object@pred_kernel)[[1]],
    dimnames(object@selectivity)[[2]],
    dimnames(object@catchability)[[2]],
    dimnames(object@interaction)[[1]],
    dimnames(object@interaction)[[2]]) ==
    object@species_params$ecotype)){
    msg <- "The species names of species_params, psi, intake_max, search_vol, std_metab, activity, pred_kernel, selectivity, catchability and interaction must all be the same"
    errors <- c(errors, msg)
  }
  # Check dimnames of w
  if(!all(c(
    dimnames(object@psi)[[2]],
    dimnames(object@intake_max)[[2]],
    dimnames(object@search_vol)[[2]],
    dimnames(object@std_metab)[[2]],
    dimnames(object@activity)[[2]],
    dimnames(object@pred_kernel)[[2]]) == 
    dimnames(object@selectivity)[[3]])){
    msg <- "The size names of psi, intake_max, search_vol, std_metab, activity, pred_kernel and selectivity must all be the same"
    errors <- c(errors, msg)
  }
  # Check dimnames of gear
  if(!isTRUE(all.equal(
    dimnames(object@catchability)[[1]],
    dimnames(object@selectivity)[[1]]))){
    msg <- "The gear names of selectivity and catchability must all be the same"
    errors <- c(errors, msg)
  }
  # Check the vector slots
  if(length(object@rr_pp) != length(object@w_full)){
    msg <- "rr_pp must be the same length as w_full"
    errors <- c(errors, msg)
  }
  if(!isTRUE(all.equal(names(object@rr_pp),dimnames(object@pred_kernel)[[3]]))){
    msg <- "Names of rr_pp and third dimension of pred_kernel must be consistent"
    errors <- c(errors, msg)
  }
  if(length(object@cc_pp) != length(object@w_full)){
    msg <- "cc_pp must be the same length as w_full"
    errors <- c(errors, msg)
  }
  if(!isTRUE(all.equal(names(object@cc_pp),dimnames(object@pred_kernel)[[3]]))){
    msg <- "Names of cc_pp and third dimension of pred_kernel must be consistent"
    errors <- c(errors, msg)
  }
  
  # SRR
  # Must have two arguments: rdi amd species_params
  if(!isTRUE(all.equal(names(formals(object@srr)), c("rdi", "species_params")))){
    msg <- "Arguments of srr function must be 'rdi' and 'species_params'"
    errors <- c(errors, msg)
  }
  
  # species_params data.frame must have columns: 
  # species, z0, alpha, eRepro
  species_params_cols <- c("species","z0","alpha","erepro")
  if (!all(species_params_cols %in% names(object@species_params))){
    msg <- "species_params data.frame must have 'species', 'z0', 'alpha' and 'erepro' columms"
    errors <- c(errors,msg)
  }
  # must also have SRR params but sorted out yet
  
  # species_params
  # Column check done in constructor
  # If everything is OK
  if (length(errors) == 0) TRUE else errors
}

# Soundtrack: White Hills - Glitter Glamour Atrocity
#' MizerParams
#'
#' A class to hold the parameters for a size based model. These parameters include the model species, their life history parameters and the size ranges of the model.
#'
#' \code{MizerParam} objects can be created using a range of \code{\link{MizerParams}} constructor methods.
#'
#' Dynamic simulations are performed using the \code{\link{project}} method on objects of this class. 
#'
#' @section Slots:
#' \describe{
#'    \item{\code{w}:}{A numeric vector of size bins used for the community (i.e. fish) part of the model. These are usually spaced on a log10 scale}
#'    \item{\code{dw}:}{The absolute difference between the size bins specified in the w slot. A vector the same length as the w slot. The final value is the same as the second to last value}
#'    \item{\code{w_full}:}{A numeric vector of size bins used for the whole model (i.e. the community and background spectra) . These are usually spaced on a log10 scale}
#'    \item{\code{dw_full}:}{The absolute difference between the size bins specified in the w_full slot. A vector the same length as the w_full slot. The final value is the same as the second to last value}
#'    \item{\code{psi}:}{An array (species x size) that holds the allocation to reproduction for each species at size}
#'    \item{\code{intake_max}:}{An array (species x size) that holds the maximum intake for each species at size}
#'    \item{\code{search_vol}:}{An array (species x size) that holds the search volume for each species at size}
#'    \item{\code{activity}:}{An array (species x size) that holds the activity for each species at size}
#'    \item{\code{std_metab}:}{An array (species x size) that holds the standard metabolism for each species at size}
#'    \item{\code{pred_kernel}:}{An array (species x predator size x prey size) that holds the predation coefficient of each predator at size on each prey size}
#'    \item{\code{rr_pp}:}{A vector the same length as the w_full slot. The size specific growth rate of the background spectrum}
#'    \item{\code{cc_pp}:}{A vector the same length as the w_full slot. The size specific carrying capacity of the background spectrum}
#'    \item{\code{species_params}:}{A data.frame to hold the species specific parameters (see Details)}
#'    \item{\code{interaction}:}{The species specific interaction matrix.}
#'    \item{\code{srr}:}{Function to calculate the realised (density dependent) recruitment. Has two arguments which are rdi and species_params}
#'    \item{\code{selectivity}:}{An array (gear x species x w) that holds the selectivity of each species by gear and species size}
#'    \item{\code{catchability}:}{An array (gear x species) that holds the catchability of each species by each gear}
#'     }
#' @note The \code{MizerParams} class is fairly complex with a large number of slots, many of which are multidimensional arrays. The dimensions of these arrays is strictly enforced so that \code{MizerParam} objects are consistent in terms of number of species and number of size classes.
#'
#' Although it is possible to build a \code{MizerParams} object by hand it is not recommended and several constructors are available.
#'
#' The \code{MizerParams} class does not hold any dynamic information, e.g. abundances or harvest effort through time. These are held in \code{\link{MizerSim}} objects.
#' @name MizerParams-class
#' @rdname MizerParams-class
#' @docType class
#' @seealso \code{\link{project}} \code{\link{MizerSim}}
#' @export
setClass("MizerParams",
         representation(
           w = "numeric",
           dw = "numeric",
           w_full = "numeric",
           dw_full = "numeric",
           psi = "array", 
           intake_max = "array",
           search_vol = "array",
           activity = "array",
           std_metab = "array",
           pred_kernel = "array",
           #z0 = "numeric",
           rr_pp = "numeric",
           cc_pp = "numeric", # was NinPP, carrying capacity of background
           species_params = "data.frame",
           interaction = "array",
           srr  = "function",
           selectivity = "array",
           catchability = "array"
         ),
         prototype = prototype(
           w = NA_real_,
           dw = NA_real_,
           w_full = NA_real_,
           dw_full = NA_real_,
           psi = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
           intake_max = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
           search_vol = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
           activity = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
           std_metab = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
           pred_kernel = array(NA,dim=c(1,1,1), dimnames = list(sp=NULL,w_pred=NULL,w_prey=NULL)),
           #z0 = NA_real_,
           rr_pp = NA_real_,
           cc_pp = NA_real_,
           #speciesParams = data.frame(),
           interaction = array(NA,dim=c(1,1), dimnames=list(predator=NULL, prey=NULL)), # which dimension is prey and which is prey?
           selectivity = array(NA, dim=c(1,1,1), dimnames = list(gear=NULL, sp=NULL, w=NULL)),
           catchability = array(NA, dim=c(1,1), dimnames = list(gear=NULL, sp=NULL))
         ),
         validity = valid_MizerParams
)


# Generic constructor
#' Constructors for objects of \code{MizerParams} class
#'
#' Constructor method for the \code{\link{MizerParams}} class. Provides the simplest way of making a \code{MizerParams} object to be used in a simulation.
#'
#' @param object A data frame of species specific parameter values (see notes below).
#' @param interaction Optional argument to specify the interaction matrix of the species (predator by prey). If missing a default interaction is used where all interactions between species are set to 1. Note that any dimnames of the interaction matrix argument are ignored by the constructor. The dimnames of the interaction matrix in the returned \code{MizerParams} object are taken from the species names in the \code{species_params} slot. This means that the order of the columns and rows of the interaction matrix argument should be the same as the species name in the \code{species_params} slot.
#' @param ... Additional arguments used to specify the dimensions and sizes of the model. These include:
#'
#' \itemize{
#'     \item{\code{min_w} The smallest size of the community spectrum}
#'     \item{\code{max_w} The largest size of the community spectrum. Default value is the largest w_inf in the community x 1.1}
#'     \item{\code{no_w} The number of size bins in the community spectrum}
#'     \item{\code{min_w_pp} The smallest size of the background spectrum}
#'     \item{\code{no_w_pp} The number of the extra size bins in the background spectrum (i.e. the difference between the number of sizes bins in the community spectrum and the full spectrum)}
#'     \item{\code{n} Scaling of the intake. Default value is 2/3} 
#'     \item{\code{p} Scaling of the standard metabolism. Default value is 0.7} 
#'     \item{\code{q} Exponent of the search volume. Default value is 0.8} 
#'     \item{\code{r_pp} Growth rate of the primary productivity. Default value is 10} 
#'     \item{\code{kappa} Carrying capacity of the resource spectrum. Default value is 1e11} 
#'     \item{\code{lambda} Exponent of the resource spectrum. Default value is (2+q-n)} 
#'     \item{\code{w_pp_cutoff} The cut off size of the background spectrum. Default value is 10} 
#'     \item{\code{f0} Average feeding level. Used to calculated \code{h} and \code{gamma} if those are not columns in the species data frame. Also requires \code{k_vb} (the von Bertalanffy K parameter) to be a column in the species data frame. If \code{h} and \code{gamma} are supplied then this argument is ignored. Default is 0.6.}
#'     \item{\code{z0pre} If \code{z0}, the mortality from other sources, is not a column in the species data frame, it is calculated as z0pre * w_inf ^ z0exp. Default value is 0.6.}
#'     \item{\code{z0exp} If \code{z0}, the mortality from other sources, is not a column in the species data frame, it is calculated as z0pre * w_inf ^ z0exp. Default value is n-1.}
#' }
#'
#' @return An object of type \code{MizerParams}
#' @note The only essential argument to the \code{MizerParams} constructor is a data frame which contains the species data. The data frame is arranged species by parameter, so each column of the parameter data frame is a parameter and each row has the parameters for one of the species in the model.
#'
#' There are some essential columns that must be included in the parameter data.frame and that do not have default values.
#' Other columns do have default values, so that if they are not included in the species parameter data frame, they will be automatically added when the \code{MizerParams} object is created. 
#' See the accompanying vignette for details of these columns.
#' 
#' An additional constructor method which takes an integer of the number of species in the model. This is only used in internally to set up a \code{MizerParams} object with the correct dimensions. It is not recommended that this method is used by users.
#' @seealso \code{\link{project}} \code{\link{MizerSim}}
#' @export
#' @docType methods
#' @rdname MizerParams-methods
#' @aliases MizerParams-method
#' @examples
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)



setGeneric('MizerParams', function(object, interaction, ...)
  standardGeneric('MizerParams'))

# Basic constructor with only the number of species as dispatching argument
# Only really used to make MizerParams of the right size and shouldn't be used by user
#' @rdname MizerParams-methods
#' @aliases MizerParams,numeric,missing-method
setMethod('MizerParams', signature(object='numeric', interaction='missing'),
          function(object, min_w = 0.001, max_w = 1000, no_w = 100,  min_w_pp = 1e-10, no_w_pp = round(no_w)*0.3, species_names=1:object, gear_names=species_names){
            #args <- list(...)

            # Some checks
            if (length(species_names) != object)
              stop("species_names must be the same length as the value of object argument")
            
            # Set up grids:
            # Community grid
            w <- 10^(seq(from=log10(min_w),to=log10(max_w),length.out=no_w))
            dw <- diff(w)
            dw[no_w] <- dw[no_w-1] # Set final dw as same as one before
            
            # Set up full grid - background + community
            # ERROR if dw > w, nw must be at least... depends on minw, maxw and nw
            if(w[1] <= dw[1])
              stop("Your size bins are too close together. You should consider increasing the number of bins, or changing the size range")
            w_full <- c(10^seq(from=log10(min_w_pp), to =  log10(w[1]-dw[1]),length.out=no_w_pp),w)
            no_w_full <- length(w_full)
            dw_full <- diff(w_full)
            dw_full[no_w_full] <- dw_full[no_w_full-1]
            
            # Basic arrays for templates
            mat1 <- array(NA, dim=c(object,no_w), dimnames = list(sp=species_names,w=signif(w,3)))
            mat2 <- array(NA, dim=c(object,no_w,no_w_full), dimnames = list(sp=species_names,w_pred=signif(w,3), w_prey=signif(w_full,3)))
            selectivity <- array(0, dim=c(length(gear_names), object, no_w), dimnames=list(gear=gear_names, sp=species_names, w=signif(w,3)))
            catchability <- array(0, dim=c(length(gear_names), object), dimnames = list(gear=gear_names, sp=species_names))
            interaction <- array(1, dim=c(object,object), dimnames = list(predator = species_names, prey = species_names))
            vec1 <- as.numeric(rep(NA, no_w_full))
            names(vec1) <- signif(w_full,3)
            
            
            # Make an empty data.frame for species_params
            # This is just to pass validity check. 
            # The project method uses the columns species z0 alpha erepro
            # so these must be in there
            # There is also a seperate function to check the dataframe that is
            # passed in by users (not used in validity check)
            species_params <- data.frame(species = species_names,
                                         z0 = NA, alpha = NA, erepro = NA)
            
            # Make an empty srr function, just to pass validity check
            srr <- function(rdi, species_params) return(0)
            
            # Make the new object
            # Should Z0, rrPP and ccPP have names (species names etc)?
            res <- new("MizerParams",
                       w = w, dw = dw, w_full = w_full, dw_full = dw_full,
                       psi = mat1, intake_max = mat1, search_vol = mat1, activity = mat1, std_metab = mat1, pred_kernel = mat2,
                       selectivity=selectivity, catchability=catchability,
                       rr_pp = vec1, cc_pp = vec1, species_params = species_params,
                       interaction = interaction, srr = srr) 
            return(res)
          }
)

# Constructor that takes the species_params data.frame and the interaction matrix

#' @rdname MizerParams-methods
#' @aliases MizerParams,data.frame,matrix-method
setMethod('MizerParams', signature(object='data.frame', interaction='matrix'),
          function(object, interaction,  n = 2/3, p = 0.7, q = 0.8, r_pp = 10, kappa = 1e11, lambda = (2+q-n), w_pp_cutoff = 10, max_w = max(object$w_inf)*1.1, f0 = 0.6, z0pre = 0.6, z0exp = n-1, normalFeeding = F, tau = 10, ...){

            # Set default values for column values if missing
            # If no gear_name column in object, then named after species
            if(!("gear" %in% colnames(object)))
              object$gear <- object$species
            no_gear <- length(unique(object$gear))
            # If no k column (activity coefficient) in object, then set to 0
            if(!("k" %in% colnames(object)))
              object$k <- 0
            # If no alpha column in object, then set to 0.6
            # Should this be a column? Or just an argument?
            if(!("alpha" %in% colnames(object)))
              object$alpha <- 0.6
            # If no erepro column in object, then set to 1
            if(!("erepro" %in% colnames(object)))
              object$erepro <- 1
            # If no sel_func column in species_params, set to 'sigmoid_length'
            if(!("sel_func" %in% colnames(object))){
              cat("\tNote: No sel_func column in species data frame. Setting selectivity to be 'knife_edge' for all species.\n")
              object$sel_func <- 'knife_edge'
              # Set default selectivity size
              if(!("knife_edge_size" %in% colnames(object))){
                cat("Note: \tNo knife_edge_size column in species data frame. Setting knife edge selectivity equal to w_mat.\n")
                object$knife_edge_size <- object$w_mat
              }
            }
            # If no catchability column in species_params, set to 1
            if(!("catchability" %in% colnames(object)))
              object$catchability <- 1
            # Sort out h column
            # If not passed in directly, is calculated from f0 and k_vb if they are also passed in
            if(!("h" %in% colnames(object))){
              cat("Note: \tNo h column in species data frame so using f0 and k_vb to calculate it.\n")
              if(!("k_vb" %in% colnames(object))){
                stop("\t\tExcept I can't because there is no k_vb column in the species data frame")
              }
              object$h <- ((3 * object$k_vb) / (object$alpha * f0)) * (object$w_inf ^ (1/3))
            }
            # Sorting out gamma column
            if(!("gamma" %in% colnames(object))){
              cat("Note: \tNo gamma column in species data frame so using f0, h, beta, sigma, lambda and kappa to calculate it.\n")
              ae <- sqrt(2*pi) * object$sigma * object$beta^(lambda-2) * exp((lambda-2)^2 * object$sigma^2 / 2)
              object$gamma <- (object$h / (kappa * ae)) * (f0 / (1 - f0))
            }
            # Sort out z0 column
            if(!("z0" %in% colnames(object))){
              cat("Note: \tNo z0 column in species data frame so using z0 = z0pre * w_inf ^ z0exp.\n")
              object$z0 = z0pre*object$w_inf^z0exp    # background natural mortality
            }
            # Sort out ks column
            if(!("ks" %in% colnames(object))){
              cat("Note: \tNo ks column in species data frame so using ks = h * 0.2.\n")
              object$ks <- object$h * 0.2
            }
            # Check essential columns: species (name), wInf, wMat, h, gamma,  ks, beta, sigma 
            check_species_params_dataframe(object)
            
            no_sp <- nrow(object)
            # Make an empty object of the right dimensions
            res <- MizerParams(no_sp, species_names=object$ecotype, gear_names=unique(object$gear), max_w=max_w,...)
            # If not w_min column in species_params, set to w_min of community
            # Check min_w argument is not > w_min in species_params
            if (!("w_min" %in% colnames(object)))
              object$w_min <- min(res@w)
            if(any(object$w_min < min(res@w)))
              stop("One or more of your w_min values is less than the smallest size of the community spectrum")
            
            # Add w_min_idx column which has the reference index of the size class closest to w_min - this is a short cut for later on and prevents repetition
            object$w_min_idx <- as.vector(tapply(object$w_min,1:length(object$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=res@w))
            
            # Start filling the slots
            res@species_params <- object
            # Check dims of interaction argument - make sure it's right

            if (!isTRUE(all.equal(dim(res@interaction), dim(interaction))))
              stop("interaction matrix is not of the right dimensions. Must be number of species x number of species")
            # Check that all values of interaction matrix are 0 - 1. Issue warning if not
            if(!all((interaction>=0) & (interaction<=1)))
              warning("Values in the interaction matrix should be between 0 and 1")
            # In case user has supplied names to interaction matrix which are wrong order
            for (dim_check in 1:length(dimnames(res@interaction))){
              if (!is.null(dimnames(interaction)[[dim_check]]) & (!(isTRUE(all.equal(dimnames(res@interaction)[[dim_check]],dimnames(interaction)[[dim_check]])))))
                
                warning("Dimnames of interaction matrix do not match the order of species names in the species data.frame. I am now ignoring your dimnames so your interaction matrix may be in the wrong order.")}
            res@interaction[] <- interaction
            # Now fill up the slots using default formulations:
            # psi - allocation to reproduction - from original Setup() function
            res@psi[] <- unlist(tapply(res@w,1:length(res@w),function(wx,w_inf,w_mat,n){
              ((1 + (wx/(w_mat))^-tau)^-1) * (wx/w_inf)^(1-n)},w_inf=object$w_inf,w_mat=object$w_mat,n=n))
            # Set w < 10% of w_mat to 0
            res@psi[unlist(tapply(res@w,1:length(res@w),function(wx,w_mat)wx<(w_mat*0.1)  ,w_mat=object$w_mat))] <- 0
            # Set all w > w_inf to 1 # Check this is right...
            res@psi[unlist(tapply(res@w,1:length(res@w),function(wx,w_inf)(wx/w_inf)>1,w_inf=object$w_inf))] <- 1
            
            res@intake_max[] <- unlist(tapply(res@w,1:length(res@w),function(wx,h,n)h * wx^n,h=object$h,n=n))
            res@search_vol[] <- unlist(tapply(res@w,1:length(res@w),function(wx,gamma,q)gamma * wx^q, gamma=object$gamma, q=q))
            res@activity[] <-  unlist(tapply(res@w,1:length(res@w),function(wx,k)k * wx,k=object$k))
            res@std_metab[] <-  unlist(tapply(res@w,1:length(res@w),function(wx,ks,p)ks * wx^p, ks=object$ks,p=p))
            # Could maybe improve this. Pretty ugly at the moment
            res@pred_kernel[] <- object$beta
            res@pred_kernel <- exp(-0.5*sweep(log(sweep(sweep(res@pred_kernel,3,res@w_full,"*")^-1,2,res@w,"*")),1,object$sigma,"/")^2)
            res@pred_kernel <- sweep(res@pred_kernel,c(2,3),combn(res@w_full,1,function(x,w)x<w,w=res@w),"*") # find out the untrues and then multiply
            # if the predation kernel is interpreted as a probability of eating an encountered prey then it should be normalised so that the predation kernel sums to 1 (see Law et al. 2009 Ecology)
            if (normalFeeding) res@pred_kernel <-sweep(res@pred_kernel,c(1,2),apply(res@pred_kernel,c(1,2),sum),"/")
            
            
            # Background spectrum
            res@rr_pp[] <- r_pp * res@w_full^(n-1) #weight specific plankton growth rate ##
            res@cc_pp[] <- kappa*res@w_full^(-lambda) # the resource carrying capacity - one for each mp and m (130 of them)
            res@cc_pp[res@w_full>w_pp_cutoff] <- 0      #set density of sizes < plankton cutoff size
            # Set the SRR to be a Beverton Holt esque relationship
            # Can add more functional forms or user specifies own
            res@srr <- function(rdi, species_params)
            {
              rdiNormal = vector(mode = "numeric", length = length(rdi))
              names(rdi) <- species_params$species
              for (i in unique(species_params$species))
              {
                rdiSp = rdi # save to manip
                rdiSp[which(names(rdi) != i)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
                
                for (i in 1:length(rdiSp))
                  # in case of NA
                  if (is.na(rdiSp[i]) == TRUE)
                    rdiSp[i] = 1e-30
                
                if (sum(rdiSp) != 0)
                  rdiNormal = rdiNormal + rdiSp / sum(rdiSp)
              }
              r_maxN = species_params$r_max * rdiNormal
              
              for (i in 1:length(r_maxN))
                # do not want to divide by 0
                if (r_maxN[i] == 0)
                  r_maxN[i] = species_params$r_max[i]
              
              return(r_maxN * rdi / (r_maxN + rdi))
              #return(species_params$r_max * rdi / (species_params$r_max+rdi))
            }
    
            # Set fishing parameters: selectivity and catchability
            # At the moment, each species is only caught by 1 gear so in species_params there are the columns: gear_name and sel_func
            # BEWARE! This routine assumes that each species has only one gear operating on it
            # So we can just go row by row through the species parameters
            # However, I really hope we can do something better soon
            for (g in 1:nrow(object)){
              # Do selectivity first
              # get args
              # These as.characters are annoying - but factors everywhere
              arg <- names(formals(as.character(object[g,'sel_func'])))
              # lop off w as that is always the first argument of the selectivity functions
              arg <- arg[!(arg %in% "w")]
              if(!all(arg %in% colnames(object)))
                stop("All of the arguments needed for the selectivity function are not in the parameter dataframe")
              # Check that there is only one column in object with the same name
              # Check that column of arguments exists
              par <- c(w=list(res@w),as.list(object[g,arg]))
              sel <- do.call(as.character(object[g,'sel_func']), args=par)
              # Dump Sel in the right place
              res@selectivity[as.character(object[g,'gear']), g, ] <- sel
              # Now do catchability
              res@catchability[as.character(object[g,'gear']), g] <- object[g,"catchability"]
            }
            
            # Remove catchabiliy from species data.frame, now stored in slot
            #params@species_params[,names(params@species_params) != "catchability"]
            res@species_params <- res@species_params[,-which(names(res@species_params)=="catchability")]
            return(res)
          }
)

# If interaction is missing, make one of the right size and fill with 1s
#' @rdname MizerParams-methods
#' @aliases MizerParams,data.frame,missing-method
setMethod('MizerParams', signature(object='data.frame', interaction='missing'),
          function(object, ...){
            #interaction <- matrix(1,nrow=nrow(object), ncol=nrow(object))
            interaction <- matrix(object$cannibalism,nrow=nrow(object), ncol=nrow(object)) #temporary
            diag(interaction) <- object$cannibalism # tweaking cannibalism
            res <- MizerParams(object,interaction, ...)
            return(res)
          })

# Check that the species_params dataset is OK
# internal only
check_species_params_dataframe <- function(species_params){
  # Check species_params dataframe (with a function) for essential cols
  # Essential columns: species (name) # wInf # wMat # h # gamma - search Volume #  ks # beta # z0
  essential_cols <- c("species","w_inf","w_mat","h","gamma","ks","beta","sigma", "z0")
  missing_cols <- !(essential_cols %in% colnames(species_params))
  if(any(missing_cols))
  {
    errors <- character()
    for (i in essential_cols[missing_cols])
      errors <- paste(errors, i, sep=" ")
    stop("You are missing these columns from the input dataframe:\n", errors)
  }
  return(TRUE)
}

