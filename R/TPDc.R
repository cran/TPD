#' Trait Probability Density of Communities
#'
#' \code{TPDc} computes the trait probability density functions (TPD) of communities, for single or multiple traits. A TPDc for each community is calculated based on the TPDs of the species (or populations) present in the community and their relative abundances. The TPDs of all species should have been calculated beforehand using the \code{\link{TPDs}} function.
#'
#' @param TPDs An object of class "TPDsp" calculated with the function \code{\link{TPDs}}, containing the TPDs of all the species or populations present in the communities whose TPDc is going to be computed.
#' @param sampUnit A matrix or data.frame containing the abundances of each species (in columns) in each community (in rows). The names of the species must be given in 'colnames(sampUnit)', whereas the names of the species must appear in 'rownames(sampUnit)'. Species names must match the names of the species used to build the TPDs object provided. In addition, in the cases in which there are different populations of the same species, the names of the communities must also match the names of the 'samples' argument provided to TPDs. In case that there is some combination of species x community with abundance greater than 0 that is not present in 'TPDs', the function will fail. NA values are not allowed in the matrix nor in the row or column names.
#'	@return \code{TPDc} returns an object of class "TPDcomm", which is a list containing the following components:
#'
#'	  \emph{data:} A list containing information used to perform the calculations, including the coordinates --in trait space-- in which the TPD function has been evaluated, the volume --in trait units-- of each cell of the grid, the length of each edges of the cells of the grid, the original trait data, the names of the species, the name of the populations in case sample is not NULL, the alpha level specified by the user, the traits of the individuals of each population, the type of TPDs calculates, which can be either "species" or "populations" depending on whether sample is or not NULL.
#'
#'    \emph{TPDc:} A list containing information related with the TPDc of each community, including the species present in each community, the abundance of those species in each community, the abundance-rescaled TPDs of each species in each community, and the TPDc of each community, which is the probability associated to each cell of the grid in in which the trait space has been divided.
#'
#'  @examples
#' # 1.  Compute the TPDs of three different species
#' traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
#' sp_iris <- iris$Species
#' example_TPDs <- TPDs(species = sp_iris, traits = traits_iris)
#'
#' #2. Three different communities with different abundances of each species
#' example_abundances <- matrix(c(c(0.5, 0.3, 0.2,
#'                                  0.1, 0.8, 0.1,
#'                                  0.5, 0,   0.5)), #I. virg. dominates; setosa absent
#'                          ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
#'                          unique(iris$Species)))
#' example_TPDc <- TPDc (TPDs = example_TPDs, sampUnit = example_abundances)

#' @export
TPDc <- function(TPDs, sampUnit){
  sampUnit <- as.matrix(sampUnit)
  # 1. species names:
  if (is.null(colnames(sampUnit)) | any(is.na(colnames(sampUnit)))) {
    stop("colnames(sampUnit), must contain the names of the species.
      NA values are not allowed")
  }
  # 2. communities names:
  if (is.null(rownames(sampUnit)) | any(is.na(rownames(sampUnit)))){
    stop("rownames(sampUnit), must contain the names of the sampling units.
      NA values are not allowed")
  }
  # 3. Data values will be inherithed from TPDs, which must be of class TPD
  if (class(TPDs) != "TPDsp") {
    stop("TPDs must be an object of class 'TPDsp', created with the function
      'TPDs'")
  }
  species <- samples <- abundances <- numeric()
  for (i in 1:nrow(sampUnit)){
    samples <- c(samples, rep(rownames(sampUnit)[i], ncol(sampUnit)))
    species <- c(species, colnames(sampUnit))
    abundances <- c(abundances, sampUnit[i, ])
  }
  nonZero <- which(abundances > 0)
  samples <- samples[nonZero]
  species <- species[nonZero]
  abundances <- abundances[nonZero]

	# Creation of lists to store results:
	results <- list()
	results$data <- TPDs$data
	results$data$sampUnit <- sampUnit
	type <- results$data$type
	# All the species or populations in 'species' must be in the species or
	#   populations of TPDs:
	if (type == "Multiple populations_One species" |
	    type == "Multiple populations_Multiple species") {
	  species_base <- paste(species, samples, sep = ".")
		if (!all(unique(species_base) %in% unique(results$data$populations))) {
			non_found_pops <- which(unique(species_base) %in%
			                        unique(results$data$populations) == 0)
			stop("All the populations TPDs must be present in 'TPDs'. Not present:\n",
			  paste(species_base[non_found_pops], collapse=" / "))
		}
	}
	if (type == "One population_One species" |
	    type == "One population_Multiple species") {
	  species_base <- species
		if (!all(unique(species_base) %in% unique(results$data$species))) {
			non_found_sps <- which(unique(species_base) %in%
			                       unique(results$data$species) == 0)
			stop("All the species TPDs must be present in 'TPDs'. Not present:\n",
			  paste(species_base[non_found_sps], collapse=" / "))
		}
	}
	# END OF INITIAL CHECKS
	# TPDc computation
	results$TPDc <- list()
	results$TPDc$species <- list()
	results$TPDc$abundances <- list()
	results$TPDc$RTPDs <- list()
	results$TPDc$TPDc <- list()

	for (samp in 1:length(unique(samples))) {
		selected_rows <- which(samples == unique(samples)[samp])
		species_aux <- species_base[selected_rows]
		abundances_aux <- abundances[selected_rows] / sum(abundances[selected_rows])
		results$TPDc$RTPDs[[samp]] <- matrix(NA, ncol=length(species_aux),
		                                  nrow=nrow((results$data$evaluation_grid)))
		TPDs_aux <- TPDs$TPDs[names(TPDs$TPDs) %in% species_aux]
		for (sp in 1:length(TPDs_aux)) {
		  selected_name <- which(names(TPDs_aux) == species_aux[sp])
			 results$TPDc$RTPDs[[samp]][,sp] <- TPDs_aux[[selected_name]] * abundances_aux[sp]
		}
		colnames(results$TPDc$RTPDs[[samp]]) <- species_aux
		results$TPDc$TPDc[[samp]] <- rowSums(results$TPDc$RTPDs[[samp]])
		results$TPDc$species[[samp]] <- species_aux
		results$TPDc$abundances[[samp]] <- abundances_aux
		names(results$TPDc$RTPDs)[samp] <- names(results$TPDc$TPDc)[samp] <-
		names(results$TPDc$species)[samp] <- names(results$TPDc$abundances)[samp] <-
		unique(samples)[samp]
	}
	class(results) <- "TPDcomm"
	return(results)
}


