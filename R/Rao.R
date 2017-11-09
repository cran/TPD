#' Rao's Quadratic Entropy and its Partition
#'
#' \code{Rao}
#'
#' @param diss An object of class "OverlapDiss", generated with the \code{\link{dissim}} function, containing the dissimilarity of the considered populations or species.
#' @param TPDc An object of class "TPDcomm", generated with the \code{\link{TPDc}} function, containing the containing the TPDc of all the communities whose functional diversity is going to be calculated. Species (or populations) identities and their relative abundance in each community will be extracted from this object.
#' @param regional Logical indicating if the correction by Villeger and Mouillot (2008) is applied or not. Defaults to TRUE.
#'
#' @return \code{Rao} returns a list containing functional diversity at different scales for the whole dataset and for pairs of communities.
#'
#'    Information for the whole dataset include: i) alpha functional diversity of each sampling unit expressed as raw rao values (alpha_rao) and in equivalent numbers alpha_eqv), ii) the average alpha functional diversity of the sampling units, calculated following de Bello et al. (2010) (mean_alpha_eqv), iii) gamma functional diversity for the whole dataset, expressed as raw rao values (gamma_rao) and in equivalent numbers (gamma_eqv), and iv) beta functional diversity for the whole dataset expressed in proportional terms (see de Bello et al. 2010) (beta_prop).
#'
#'    Information for pairs of communities (contained in the element pairwise) include the average alpha (expressed in equivalent numbers) of each pair of communities, gamma of each pair of communities and beta functional diversity for each pair of communities.
#'
#' @examples
#' # 1.  Compute the TPDs of three different species.
#' traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
#' sp_iris <- iris$Species
#' TPDs_iris <- TPDs(species = sp_iris, traits_iris)
#' #2. Compute the dissimilarity between the three species:
#' dissim_iris_sp <- dissim(TPDs = TPDs_iris)

#' #3. Compute the TPDc of five different communities:
#' abundances_comm_iris <- matrix(c(c(0.9,  0.1, 0), # setosa dominates
#'                                  c(0.4,  0.5, 0.1 ),
#'                                  c(0.15, 0.7, 0.15), #versicolor dominates
#'                                  c(0.1,  0.5, 0.4),
#'                                  c(0,    0.1, 0.9)), #virginica dominates
#'                            ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:5),
#'                              unique(iris$Species)))
#' TPDc_iris <- TPDc( TPDs = TPDs_iris, sampUnit = abundances_comm_iris)
#'
#' #4. Compute Rao:
#' Rao_iris <- Rao(diss = dissim_iris_sp, TPDc = TPDc_iris)

#' @export
Rao <- function(diss = NULL, TPDc = NULL, regional = TRUE ) {
	# INITIAL CHECKS:
	# 1. Both diss and TPDc must be supplied
	if (is.null(TPDc) | is.null(TPDs)) {
		stop("Both 'diss' and 'TPDc' must be supplied")
	}
  # 2. Both diss and TPDc must be of the correct class
  if (class(TPDc) != "TPDcomm") {
    stop("TPDc must be an object of class 'TPDcomm', created with the function
        'TPDc'")
  }
  if (class(diss) != "OverlapDiss") {
    stop("diss must be an object of class 'OverlapDiss', created with the
        function 'FunctionalDissimilarity'")
  }
  # 3. Determine if dissimilarities are calculated for species or populations,
  #   and check that all species or populations in the communities are present
  #   in the dissimilarity matrix.
  type <- TPDc$data$type
  if (type == "One population_One species" |
      type == "One population_Multiple species") {
    if (is.null(diss$species)) {
      stop("TPDc contains information at the ", type," level, whereas diss
          contains information at the population level (more than one population
          for some species). Recalculate either TPDc or diss so that their
          levels match")
    }
    sp_in_samples <- unique(TPDc$data$species)
    diss <- diss$species$dissimilarity
    sp_in_diss <- colnames(diss)
    if (!all(sp_in_samples %in% sp_in_diss)) {
      stop("There are some ", type, " in TPDc that are not present in diss")
    }
    samples_matrix <- matrix(0, nrow = length(TPDc$TPDc$species),
      ncol = length(sp_in_samples), dimnames = list(names(TPDc$TPDc$species),
      sp_in_samples))
    for (i in 1:length(TPDc$TPDc$species)){
      species_aux <- TPDc$TPDc$species[[i]]
      abundances_aux <- TPDc$TPDc$abundances[[i]]
      for (j in 1:length(species_aux)){
        samples_matrix[i, species_aux[j]] <- abundances_aux[j]
      }
    }
  }
  if (type == "Multiple populations_One species" |
      type == "Multiple populations_Multiple species") {
    if (is.null(diss$populations)) {
      stop("TPDc contains information at the ", type," level (more than one
          population for some species), whereas diss contains information at the
          species level (only one population per species). Recalculate either
          TPDc or diss so that their levels match")
    }
    pop_in_samples <- unique(TPDc$data$populations)
    diss <- diss$populations$dissimilarity
    pop_in_diss <- colnames(diss)
    if (!all(pop_in_samples %in% pop_in_diss)) {
      stop("There are some ", type, " in TPDc that are not present in diss")
    }
    samples_matrix <- matrix(0, nrow = length(TPDc$TPDc$species),
      ncol = length(pop_in_samples),
      dimnames = list(names(TPDc$TPDc$species), pop_in_samples))
    for (i in 1:length(TPDc$TPDc$species)) {
      species_aux <- TPDc$TPDc$species[[i]]
      abundances_aux <- TPDc$TPDc$abundances[[i]]
      for (j in 1:length(species_aux)){
        samples_matrix[i, species_aux[j]] <- abundances_aux[j]
      }
    }
  }

  results <- list()
  results$mean_alpha_eqv <- numeric()
  results$alpha_eqv <- results$alpha_rao <- rep(NA, nrow(samples_matrix))
  names(results$alpha_eqv) <- names(results$alpha_rao) <- rownames(samples_matrix)
  results$gamma_eqv <- results$gamma_rao <- numeric()
  results$beta_prop <- numeric()
  results$pairwise <- list()
  results$pairwise$gamma <- results$pairwise$beta_prop <-
    results$pairwise$mean_alpha <- matrix(NA, nrow = nrow(samples_matrix),
      ncol = nrow(samples_matrix),dimnames = list(rownames(samples_matrix),
      rownames(samples_matrix)))
  rao_raw_region <- rep(NA, nrow(samples_matrix))
  rao_raw <- rep(NA, nrow(samples_matrix))
  #alpha
  for (i in 1:nrow(samples_matrix)) {
    sp_aux <- colnames(samples_matrix)[samples_matrix[i, ]>0]
    abund_aux <- samples_matrix[i, sp_aux]
    if (length(sp_aux) == 1) {
      rao_raw_region[i] <- rao_raw[i] <- 0
    }
    if (length(sp_aux) > 1) {
      sample_dis<-diss[sp_aux, sp_aux]
      sample_weights <- outer(abund_aux , abund_aux)
      w <- 1 / nrow(samples_matrix)
      rao_raw_region[i] <- sum(sample_weights * sample_dis) * w
      rao_raw[i] <- sum(sample_weights * sample_dis)
    }
    Rao_aux <- ifelse(regional, rao_raw_region[i], rao_raw[i])
    results$alpha_rao[i] <- Rao_aux
    results$alpha_eqv[i] <- 1 / (1 - Rao_aux)
  }
  #gamma_pair
  for (i in 1:nrow(samples_matrix)) {
    for (j in 1:nrow(samples_matrix)) {
      if (j > i){
        samples_aux <- colMeans(samples_matrix[c(i,j), ])
        sp_aux <- names(samples_aux)[samples_aux>0]
        abund_aux <- samples_aux[sp_aux]
        if (length(sp_aux) == 1) {
          Rao_aux <- 0
        }
        if (length(sp_aux) > 1) {
          sample_dis<-diss[sp_aux, sp_aux]
          sample_weights <- outer(abund_aux , abund_aux)
          Rao_aux <- sum(sample_weights * sample_dis)
        }
        results$pairwise$mean_alpha[i, j] <- results$pairwise$mean_alpha[j, i] <-
          1 / (1 - mean(c(rao_raw[i], rao_raw[j])))
        results$pairwise$gamma[i, j] <- results$pairwise$gamma[j, i] <-
          1 / (1 - Rao_aux)
      }
    }
  }
  results$pairwise$beta_prop <- (results$pairwise$gamma -
      results$pairwise$mean_alpha) / (0.5 * results$pairwise$gamma)
  #gamma_regional
  samples_aux <- colMeans(samples_matrix)
  sp_aux <- names(samples_aux)[samples_aux>0]
  abund_aux <- samples_aux[sp_aux]
  if (length(sp_aux) == 1) {
    Rao_aux <- 0
  }
  if (length(sp_aux) > 1) {
    sample_dis<-diss[sp_aux, sp_aux]
    sample_weights <- outer(abund_aux , abund_aux)
    Rao_aux <- sum(sample_weights * sample_dis)
  }
  results$gamma_rao <- Rao_aux
  results$gamma_eqv <- 1 / (1 - Rao_aux)
  results$mean_alpha_eqv <- 1 / (1 - mean(rao_raw))
  results$beta_prop <- 100 * (results$gamma_eqv - results$mean_alpha_eqv) /
      results$gamma_eqv

	return(results)
}
