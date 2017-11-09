#' Functional Evenness, Richness and Divergence of Communities,
#' Species or Populations
#'
#' \code{REND} computes Functional Richness, Functional Evenness and
#'    Functional Divergence, the three primary components of functional
#'    diversity (Mason et al. 2005) for single or multiple traits. Although
#'    these components were originally intended to be calculated for communities,
#'    \code{REND} also allows to compute them for populations or species.
#'    In the case of communities, all the calculations are based on the TPDc of
#'    the considered communities; therefore results are independent of any
#'    underlying feature of the species that compose the communities.
#'
#' @param TPDc An object of class "TPDcomm", generated with the
#'    \code{\link{TPDc}} function, containing the TPDc of the considered
#'    communities.
#' @param TPDs An object of class "TPDsp", generated with the
#'    \code{\link{TPDs}} function, containing the TPDs of the considered
#'    populations or species.
#'
#' @return \code{REND} returns a list with an element for each of the
#'    provided parameters (ie. communities and/or populations/species). These
#'    lists contain in turn one element for the Functional Richness  of each
#'    unit, one for Functional Evenness, and one for Functional Divergence.
#'
#' @references Mason, NWH, Mouillot, D, Lee, WG and Wilson, JB (2005),
#'    Functional richness, functional evenness and functional divergence:
#'    the primary components of functional diversity. \emph{Oikos},
#'    111: 112--118.
#'
#'   @examples
#' # 1.  Compute the TPDs of five different species. SP3 is in the center of
#' #   the trait space, and the rest of species in the corners
#' set.seed(1)
#' species_ex <- c(rep("SP1",20), rep("SP2",20), rep("SP3",20), rep("SP4",20),
#'              rep("SP5",20))
#' traits_ex <- data.frame(trait1 = c(rnorm(20, 10, 1),
#'                                    rnorm(20, 10, 1),
#'                                    rnorm(20, 15, 1),
#'                                    rnorm(20, 20, 1),
#'                                    rnorm(20, 20, 1)),
#'                         trait2 = c(rnorm(20, 10, 1),
#'                                    rnorm(20, 20, 1),
#'                                    rnorm(20, 15, 1),
#'                                    rnorm(20, 10, 1),
#'                                    rnorm(20, 20, 1)))
#' species_TPDs <- TPDs (species = species_ex, traits = traits_ex)
#' #2. Five different communities with different abundances of each species
#' abundances_ex <- matrix(c(0.05, 0.05, 0.8,  0.05, 0.05, # 1. Low divergence
#'                    0.9,  0,    0,    0,    0.1,  # 2. High divergence
#'                    0,    0,    1,    0,    0,    # 3. Low Richness
#'                    0.2,  0.2,  0.2,  0.2,  0.2,  # 4. High Evenness
#'                    0.8,  0.05, 0.05, 0.05, 0.05), # 5. Low Evenness
#'                    ncol = 5, byrow = TRUE, dimnames = list(paste0("Comm.",1:5),
#'                    unique(species_ex)))
#'
#' example_TPDc <- TPDc (TPDs = species_TPDs, sampUnit = abundances_ex)
#' #3. Estimate functional richness, evenness and divergence
#' example_RicEveDiv <- REND (TPDc = example_TPDc)

#' @export
REND <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
	# 1. At least one of TPDc or TPDs must be supplied.
	if (is.null(TPDc) & is.null(TPDs)) {
		stop("At least one of 'TPDc' or 'TPDs' must be supplied")
	}
	if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
	  stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
	if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
	  stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
	}
	# Creation of lists to store results:
	results <- list()
	# 1. Functional Richness
	Calc_FRich <- function(x) {
		results_FR <- numeric()
		if (class(x) == "TPDcomm") {
			TPD <- x$TPDc$TPDc
			names_aux <- names(x$TPDc$TPDc)
			cell_volume <- x$data$cell_volume
		}
		if (class(x) == "TPDsp") {
			TPD <- x$TPDs
			names_aux <- names(x$TPDs)
			cell_volume <- x$data$cell_volume
		}
		for (i in 1:length(TPD)) {
			TPD_aux <- TPD[[i]]
			TPD_aux[TPD_aux > 0] <- cell_volume
			results_FR[i] <- sum(TPD_aux)
		}
		names(results_FR) <- names_aux
		return(results_FR)
	}
	# 2. Functional Evenness
	Calc_FEve <- function(x) {
		results_FE <- numeric()
		if (class(x) == "TPDcomm") {
			TPD <- x$TPDc$TPDc
			names_aux <- names(x$TPDc$TPDc)
			cell_volume <- x$data$cell_volume
		}
		if (class(x) == "TPDsp") {
			TPD <- x$TPDs
			names_aux <- names(x$TPDs)
			cell_volume <- x$data$cell_volume
		}
		for (i in 1:length(TPD)) {
			TPD_aux <- TPD[[i]][TPD[[i]] > 0]
			TPD_eve <- rep((1 / length(TPD_aux)), times = length(TPD_aux))
			results_FE[i] <- sum(pmin(TPD_aux, TPD_eve))
		}
		names(results_FE) <- names_aux
		return(results_FE)
	}
	# 3. Functional Divergence
	Calc_FDiv <- function(x) {
		results_FD <- numeric()
		if (class(x) == "TPDcomm") {
			TPD <- x$TPDc$TPDc
			evaluation_grid<-x$data$evaluation_grid
			names_aux <- names(x$TPDc$TPDc)
			cell_volume <- x$data$cell_volume
		}
		if (class(x) == "TPDsp") {
			TPD <- x$TPDs
			evaluation_grid<-x$data$evaluation_grid
			names_aux <- names(x$TPDs)
			cell_volume <- x$data$cell_volume
		}
		for (i in 1:length(TPD)) {
			functional_volume <- evaluation_grid[TPD[[i]]>0 , , drop=F]
			# Functional volume has to be standardised so that distances are
      # independent of the scale of the axes
      for (j in 1:ncol(functional_volume)){
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) /
          (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      TPD_aux <- TPD[[i]][TPD[[i]] > 0]
			# 1. Calculate the center of gravity
			COG <- colMeans(functional_volume, na.rm=T)
			# 2. Calculate the distance of each point in the functional volume to the
			#   COG:
			dist_COG <- function(x, COG) {
				result_aux<-stats::dist(rbind(x, COG))
				return(result_aux)
			}
			COGDist <- apply(functional_volume, 1, dist_COG, COG)
			# 3. Calculate the mean of the COGDist's
			meanCOGDist <- mean(COGDist)
			# 4. Calculate the sum of the abundance-weighted deviances for distaces
			#   from the COG (AWdistDeviances) and the absolute abundance-weighted
			#   deviances:
			distDeviances <- COGDist-meanCOGDist
			AWdistDeviances <- sum(TPD_aux * distDeviances)
			absdistDeviances <- abs(COGDist-meanCOGDist)
			AWabsdistDeviances <- sum(TPD_aux * absdistDeviances)
			#Finally, calculate FDiv:
			results_FD[i] <- (AWdistDeviances + meanCOGDist) /
			                 ( AWabsdistDeviances +  meanCOGDist)
		}
		names(results_FD) <- names_aux
		return(results_FD)
	}
	# IMPLEMNENTATION
	if (!is.null(TPDc)) {
		results$communities <- list()
		message("Calculating FRichness of communities")
		results$communities$FRichness <- Calc_FRich(TPDc)
		message("Calculating FEvenness of communities")
		results$communities$FEvenness <- Calc_FEve(TPDc)
		message("Calculating FDivergence of communities")
		results$communities$FDivergence <- Calc_FDiv(TPDc)
	}
	if (!is.null(TPDs)) {
	  if (TPDs$data$type == "One population_One species" |
	      TPDs$data$type == "One population_Multiple species") {
	    results$species <- list()
			message("Calculating FRichness of species")
			results$species$FRichness <- Calc_FRich(TPDs)
			message("Calculating FEvenness of species")
			results$species$FEvenness <- Calc_FEve(TPDs)
			message("Calculating FDivergence of species")
			results$species$FDivergence <- Calc_FDiv(TPDs)
		} else {
			results$populations <- list()
			message("Calculating FRichness of populations")
			results$populations$FRichness <- Calc_FRich(TPDs)
			message("Calculating FEvenness of populations")
			results$populations$FEvenness <- Calc_FEve(TPDs)
			message("Calculating FDivergence of populations")
			results$populations$FDivergence <- Calc_FDiv(TPDs)
		}
		if (TPDs$data$method == "mean") {
		  message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
		}
	}
	return(results)
}


