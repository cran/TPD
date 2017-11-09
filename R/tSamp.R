#'  Trait Values Sampling
#'
#' \code{tSamp} samples (with replacement) trait values from populations, species or communities. The probability of sampling each trait --or combination of traits in multidimensional cases-- is proportional to the value of TPDs or TPDc for the corresponding cell (the trait space is divided in a grid composed of cells, see \code{\link{TPDs}} for further information).
#'
#' @param TPDc An object of class "TPDcomm", generated with the \code{\link{TPDc}} function, containing the TPDc of the considered communities.
#' @param TPDs An object of class "TPDsp", generated with the \code{\link{TPDs}} function, containing the TPDs of the considered populations or species.
#' @param size Non-negative integer giving the number of observations to choose. Defaults to 1.
#'
#' @return \code{tSamp} returns a list containing sampled trait values for each community of TPDc or species/populations from TPDs.
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
#'
#' #3. Sample 1,000 trait values from each species and community
#' example_sampling <- tSamp(TPDc = example_TPDc, TPDs = example_TPDs,
#'                                  size = 1000)

#' @export
tSamp <- function(TPDc = NULL, TPDs = NULL, size = 1){
	# INITIAL CHECKS:
	# 1. At least one of TPDc or TPDs must be supplied.
	if (is.null(TPDc) & is.null(TPDs)) {
		stop("At least one of 'TPDc' or 'TPDs' must be supplied")
	}
	results <- list()
	sample_traits <- function(x) {
		results_samp <- list()
		if (class(x)=="TPDcomm") {
			TPD <- x$TPDc$TPDc
			names_aux <- names(x$TPDc$TPDc)
			grid <- x$data$evaluation_grid
		}
		if (class(x)=="TPDsp") {
			TPD <- x$TPDs
			names_aux <- names(x$TPDs)
			grid <- x$data$evaluation_grid
		}
		for (i in 1:length(TPD)) {
			TPD_aux <- TPD[[i]]
			selected_rows <- sample(1:nrow(grid), size = size, replace = T,
			                        prob = TPD_aux)
			results_samp[[i]] <- grid[selected_rows, ]
		}
		names(results_samp) <- names_aux
		return(results_samp)
	}
	# IMPLEMENTATION
	if (!is.null(TPDc)) {
		results$communities_samples <- sample_traits(TPDc)
	}
	if (!is.null(TPDs)) {
	  if (TPDs$data$type == "One population_One species" |
	      TPDs$data$type == "One population_Multiple species") {
	  	results$species_samples <- sample_traits(TPDs)
		} else {
			results$populations_samples <- sample_traits(TPDs)
		}
	}
	return(results)
}


