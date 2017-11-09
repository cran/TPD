#'  Functional Redundancy of Communities
#'
#' \code{redundancy} calculates the functional redundancy of communities, considering single or multiple traits. The functional volume (indicated by Functional Richness) occupied by a community with high functional redundancy should not decrease substantially when some species are lost, and vice versa.
#'
#' @param TPDc An object of class "TPDcomm", generated with the \code{\link{TPDc}} function, containing the TPDc of the considered communities.
#' @return \code{redundancy} returns a list containing the functional redundancy values of all the communities from TDPc, along with the number of species of each community.
#' @examples
#' #1. Compute the TPDs of three different species.
#' traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
#' sp_iris <- iris$Species
#' TPDs_iris <- TPDs(species = sp_iris, traits_iris)
#'
#' #2. Compute the TPDc of five different communities:
#' abundances_comm_iris <- matrix(c(c(0.9,  0.05, 0.05), #I. setosa dominates
#'                                  c(0.0,  0.5,  0.5 ), #I. setosa absent
#'                                  c(0.33, 0.33, 0.33), #Equal abundances
#'                                  c(0.1,  0.45, 0.45), #Versicolor and virginica dominate
#'                                  c(0.5,  0,    0.5)), #versicolor absent
#'                            ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:5),
#'                            unique(iris$Species)))
#' TPDc_iris <- TPDc( TPDs = TPDs_iris, sampUnit = abundances_comm_iris)
#'
#' #3. Estimate functional redundancy
#' FRed_iris <- redundancy(TPDc = TPDc_iris)

#' @export
redundancy <- function(TPDc = NULL) {
	if (class(TPDc) != "TPDcomm") {
		stop("TPDc must be an object of class TPDcomm generated with the TPDc
		    function")
	}
  x <- TPDc
	results <- list()
	results$redundancy <- numeric()
	for (i in 1: length(x$TPDc$TPDc)) {
		TPDc_aux <- x$TPDc$TPDc[[i]]
		RTPDs_aux <- x$TPDc$RTPDs[[i]]
		RTPDs_aux[RTPDs_aux > 0] <- 1
		M <- rowSums(RTPDs_aux)
		results$redundancy[i] <- sum(M * TPDc_aux) - 1
		results$richness[i] <- sum(x$TPDc$abundances[[i]] >0)
	}
	names(results$redundancy) <- names(results$richness) <- names(x$TPDc$TPDc)
	return(results)
}
