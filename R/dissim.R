#' Overlap-Based Functional Dissimilarity and its Decomposition
#'
#' \code{dissim} calculates the  functional  dissimilarity between pairs of communities or populations, as well as its decomposition into shared and non-shared trait volume.
#'
#' @param TPDc An object of class "TPDcomm", generated with the \code{\link{TPDc}} function, containing the TPDc of the considered communities.
#' @param TPDs An object of class "TPDsp", generated with the \code{\link{TPDs}} function, containing the TPDs of the considered populations or species

#' @return \code{dissim} returns lists for populations, species or communities, depending on the input data. These lists contain the overlap-based functional dissimilarity between all pairs of populations/species/communities, along with the decomposition of dissimilarity between shared and non-shared trait volume.
#'
#' @examples
#' # 1.  Compute the TPDs of three different species:
#' traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
#' sp_iris <- iris$Species
#' TPDs_iris <- TPDs(species = sp_iris, traits_iris)
#'
#'  #2. Compute the TPDc of three different communities:
#' abundances_comm_iris <- matrix(c(c(0.9, 0.1, 0), #I. setosa dominates
#'                          c(0.0, 0.9,  0.1 ), #I. Versic. dominates; setosa absent
#'                          c(0.0, 0.1,  0.9 )), #I. virg. dominates; setosa absent
#'                          ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
#'                           unique(iris$Species)))
#' TPDc_iris <- TPDc(TPDs = TPDs_iris, sampUnit = abundances_comm_iris)
#'
#' #3. Estimate functional dissimilarity
#' example_dissimilarity <- dissim (TPDc = TPDc_iris, TPDs = TPDs_iris)

#' @export
dissim <- function(TPDc = NULL, TPDs = NULL) {
	# INITIAL CHECKS:
	# 1. At least one of TPDc or TPDs must be supplied.
	if (is.null(TPDc) & is.null(TPDs)) {
		stop("At least one of 'TPDc' or 'TPDs' must be supplied")
	}
	results <- list()
	Calc_dissim <- function(x) {
	  # 1. BetaO (functional dissimilarity)
		results_samp <- list()
		if (class(x) == "TPDcomm") {
			TPD <- x$TPDc$TPDc
			names_aux <- names(x$TPDc$TPDc)
		}
		if (class(x) == "TPDsp") {
			TPD <- x$TPDs
			names_aux <- names(x$TPDs)
		}
		results_samp$dissimilarity <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
			                                  dimnames = list(names_aux, names_aux))
		results_samp$P_shared <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
			                             dimnames = list(names_aux, names_aux))
		results_samp$P_non_shared <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
			                                 dimnames = list(names_aux, names_aux))
		for (i in 1: length(TPD)) {
			TPD_i <- TPD[[i]]
			for (j in 1:length(TPD)) {
			  if (i > j) {
					TPD_j <- TPD[[j]]
					O_aux <- sum(pmin(TPD_i, TPD_j))
					shared_aux <- which(TPD_i > 0 & TPD_j > 0)
					A_aux <- sum(pmax(TPD_i[shared_aux], TPD_j[shared_aux])) - O_aux
					only_in_i_aux <- which(TPD_i > 0 & TPD_j == 0)
					B_aux <- sum(TPD_i[only_in_i_aux])
					only_in_j_aux <- which(TPD_i == 0 & TPD_j > 0)
					C_aux <- sum(TPD_j[only_in_j_aux])

					results_samp$dissimilarity[i, j] <-
					results_samp$dissimilarity[j, i] <-	1 - O_aux
					if (results_samp$dissimilarity[j, i] == 0) {
					  results_samp$P_non_shared[i, j] <- NA
					  results_samp$P_non_shared[j, i] <- NA
					  results_samp$P_shared[i, j] <- NA
					  results_samp$P_shared[j, i] <- NA
					}	else {
					  results_samp$P_non_shared[i, j] <-
					  results_samp$P_non_shared[j, i] <-
					  (2 * min(B_aux, C_aux)) / (A_aux + 2 * min(B_aux, C_aux))
					  results_samp$P_shared[i, j] <- results_samp$P_shared[j, i] <-
					  1 - results_samp$P_non_shared[i, j]
					}
				}
				if (i == j) {
					results_samp$dissimilarity[i, j] <- 0
		    }
			}
		}
	return(results_samp)
	}
	# IMPLEMENTATION
	if (!is.null(TPDc)) {
	  message("Calculating dissimilarities between ", length(TPDc$TPDc$TPDc)," communities. It might take some time")
		results$communities <- Calc_dissim(TPDc)
	}
	if (!is.null(TPDs)) {
	  if (TPDs$data$type == "One population_One species" |
	      TPDs$data$type == "One population_Multiple species") {
	    message("Calculating dissimilarities between ", length(TPDs$TPDs)," species. It might take some time")
	    results$species <- Calc_dissim(TPDs)
		} else {
		  message("Calculating dissimilarities between ", length(TPDs$TPDs)," populations. It might take some time")
		  results$populations <- Calc_dissim(TPDs)
		}
	}
	class(results) <- "OverlapDiss"
	return(results)
}
