#' Functional Uniqueness of Ecological Units
#'
#' \code{uniqueness} estimates the functional uniqueness of species, communities by comparing the TPD of lower  levels (i.e. species), with that of higher levels (i.e. communities). TPD's are compared by means of overlap. High overlap means low uniqueness (i.e. the species traits are frequent in the community), whereas low overlap means high uniqueness. Uniqueness is then estimated as 1-overlap. The function is hence basically the same as 'dissim', with some slight modifications. Despite functional uniqueness can be estimated at any scale, current implementation is limited to species within communities (although communities can be easily created to represent regions, or regional pools of species).
#'
#' @param TPDs An object of class "TPDsp", generated with the \code{\link{TPDs}} function, containing the TPDs of the considered species
#' @param TPDc An object of class "TPDcomm", generated with the \code{\link{TPDc}} function, containing the TPDc of the considered communities.
#'
#' @return \code{uniqueness} returns a matrix, with the communities in rows and the species in columns. The values in the matrix represent the functional uniqueness of each species in each community. Very unique species will have values close to 1, whereas non-unique species will have values close to 0.
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
#'
#' #3. Calculate the uniqueness of each species in each community
#' example_uniqueness <- uniqueness (TPDs = example_TPDs, TPDc = example_TPDc)
#'
#' @export

uniqueness <- function(TPDs = NULL, TPDc = NULL) {
  # INITIAL CHECKS:
  # 1. Both TPDc and TPDs must be supplied.
  if (is.null(TPDs) | is.null(TPDc)) {
    stop("Both 'TPDs' and 'TPDc' must be supplied")
  }
  if (class(TPDs) != "TPDsp" | class(TPDc) != "TPDcomm") {
    stop("'TPDs' and 'TPDc' have to be of the classes 'TPDsp' and 'TPDcomm' respectively")
  }
  Calc_uniq <- function(x, y) { #where x is TPDs and y is TPDc
    TPDs_aux <- x$TPDs
    TPDc_aux <- y$TPDc$TPDc
    names_aux_com <- names(y$TPDc$TPDc)
    names_aux_sp <- names(x$TPDs)
    results <- matrix(NA, ncol= length(TPDs_aux), nrow= length(TPDc_aux),
                      dimnames = dimnames(y$data$sampUnit))
    for (i in 1: length(TPDc_aux)) {
      TPD_ci <- TPDc_aux[[i]]
      for (j in 1:length(TPDs_aux)) {
        TPD_sj <- TPDs_aux[[j]]
        O_aux <- sum(pmin(TPD_ci, TPD_sj))
        results[i, j] <- 1 - O_aux
      }
    }
    return(results)
  }
  # IMPLEMENTATION
  results <- Calc_uniq(x = TPDs, y = TPDc)
  return(results)
}
