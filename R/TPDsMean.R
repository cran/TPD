#' Creating TPDs without individual observations
#'
#' \code{TPDsMean} estimates the TPDs of species using the mean trait values and covariance matrix of traits.
#'  It is most useful when there is no trait information at the individual level, but the mean and
#'  variance (and optionally covariance) of traits are known.
#'
#' @param species A vector containing the names of the species whose TPD is to be calculated.
#'  The length of 'species' must match the number of rows of 'traits', and the length of 'sigmas'.
#'  NA values are not allowed.
#' @param means A matrix or data.frame containing the average trait values of each  species (or population).
#'  Species are in rows and  traits in columns, with one column for each trait. NA values are not allowed.
#' @param sds A matrix or data.frame containing the standard deviation of trait values of each  species
#'  (or population). Species are in rows and  traits in columns, with one column for each trait. NA values are not allowed.
#' @param samples A vector containing the identity of the sampling unit in which each individual was present.
#'  Defaults to NULL, in which case, a TPDs is calculated for each species. If it is not NULL,
#'  a TPDs is calculated for each combination of 'species x sampling unit'  (ie, a TPDs for each population).
#'  NA values are not allowed.
#' @param covar Logical; if TRUE, the covariance between traits is calculated using the mean trait values and
#'  considered in the construction of the multivariate normal distributions. Defaults to FALSE.
#' @param alpha A number between 0 and 1, indicating the proportion of the probability density function of
#'  each population to include. A value of 1 includes the whole density function, but may be sensitive
#'  to the presence of outliers. Defaults to 0.95.
#' @param trait_ranges A vector or a list indicating the range of trait values that will be considered in
#'  the calculations. If a vector is provided, each	element should indicate the interval (in number of
#'  standard deviations) by which the range of	each trait should be expanded in each direction.
#'  In that case, an interval of n standard deviations is calculated around each species mean trait
#'  value, and the absolute maximum and minimum across all species are selected as the range for
#'  each trait.If a list is provided, it should contain the range (minimum and maximum) of trait
#'  values that will be considered.Each element of the vector or list corresponds with one trait.
#'  The order of the traits must be the same as the order of the columns in the 'traits' arguments.
#'  Defaults to	NULL, in which case, the ranges of all the traits are automatically calculated by
#'  expanding the range of the columns in 'means' by 5 times the values in 'sds' in each direction.
#'  Trait ranges that are too short may result in an inadequate characterization of TPDs.
#' @param n_divisions The number of equal-length parts in which each trait should be divided to calculate
#'  the grid in which calculations are based. Note	the number of cells composing the grid increases
#'  exponentially as	dimensionality increases, which can result in  high computation times. Defaults to NULL,
#'  in which case one trait is divided into 1000 parts, two traits are divided into 200 parts (40,000 cells),
#'  three traits are divided into 50 parts (125,000 cells), and 4 traits are divided into 25 parts
#'  (390,625 cells).
#' @param tolerance A number between 0 and 1, giving the admissible proportion of deviation from 1 in the
#'  integral of the TPDs of each population. Integrals can be lower than 1 when the extent of the evaluation
#'  grid is not enough to capture all the probability density function of the species. These problems
#'  are usually solved by increasing 'trait_ranges'. When the absolute deviation is greater than 'tolerance',
#'  a warning message is produced, but the TPDsMean function does not fail. Defaults to 0.05.
#'
#' @return \code{TPDsMean} returns an object of class "TPDsp", which is a list containing the following components:
#'
#'	\emph{data:} A list containing information used to perform the calculations, including the coordinates
#'   --in trait space-- in which the TPD function has been evaluated, the volume --in trait units--
#'   of each cell of the grid, the length of each edges of the cells of the grid, the original trait data
#'   (means and sds matrices), the names of the species, the alpha level specified by the user,
#'   and the type of TPDs calculated.
#'
#'  \emph{TPDs:} A list, with one element per species or population, containing the probability associated
#'  to each cell of the grid in which the trait space has been divided.
#'
#' @examples
#'
#' # 1.  Compute the TPDs of three different species (1 dimension)
#' sp_ex <- unique(iris$Species)
#' mt1 <- tapply(iris[, "Sepal.Length"], iris$Species, mean)
#' means_ex <- matrix(c(mt1), ncol=1)
#' st1 <- tapply(iris[, "Sepal.Length"], iris$Species, sd)
#' sds_ex <- matrix(c(st1), ncol=1)
#' TPDs_iris<- TPDsMean(species = sp_ex, means = means_ex, sds = sds_ex)
#'
#' # 2.  Compute the TPDs of three different species (2 dimensions)
#' sp_ex <- unique(iris$Species)
#' mt1 <- tapply(iris[, "Sepal.Length"], iris$Species, mean)
#' mt2 <- tapply(iris[, "Sepal.Width"], iris$Species, mean)
#' means_ex <- matrix(c(mt1, mt2), ncol=2)
#' st1 <- tapply(iris[, "Sepal.Length"], iris$Species, sd)
#' st2 <- tapply(iris[, "Sepal.Width"], iris$Species, sd)
#' sds_ex <- matrix(c(st1, st2), ncol=2)
#' TPDs_iris<- TPDsMean(species = sp_ex, means = means_ex, sds = sds_ex)
#'
#' # 3.  Two different populations of each species
#' samples_aux <- rep(c(rep(1, 25), rep(2, 25)), 3)
#' sp_ex <- rep(unique(iris$Species), each=2)
#' mt1 <- tapply(iris[, "Sepal.Length"], (paste0(iris$Species,samples_aux)), mean)
#' mt2 <- tapply(iris[, "Sepal.Width"], (paste0(iris$Species,samples_aux)), mean)
#' means_ex <- matrix(c(mt1, mt2), ncol=2)
#' st1 <- tapply(iris[, "Sepal.Length"], (paste0(iris$Species,samples_aux)), sd)
#' st2 <- tapply(iris[, "Sepal.Width"], (paste0(iris$Species,samples_aux)), sd)
#' sds_ex <- matrix(c(st1, st2), ncol=2)
#' samples_ex<- rep(c("Comm.1","Comm.2"),3)
#' TPDs_iris_pop <- TPDsMean (species = sp_ex, means = means_ex, sds = sds_ex,
#'    samples = samples_ex)
#'
#' @import mvtnorm
#'
#' @export
TPDsMean<- function(species, means, sds, covar = FALSE, alpha = 0.95, samples = NULL,
                 trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05) {

  # INITIAL CHECKS:
  #   1. Compute the number of dimensions (traits):
  means <- as.matrix(means)
  dimensions <- ncol(means)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the
         number of dimensions")
  }
  #   2. sds and means must have the same dimensions:
  sds <- as.matrix(sds)
  if (all(dim(means) != dim(sds))) {
    stop("'means' and 'sds' must have the same dimensions")
  }
  #   3. species and means must have the same "dimensions":
  if (length(species) != nrow(means)) {
    stop("The length of 'species' does not match the number of rows of 'means'
         and 'sds'")
  }
  #	4. NA's not allowed in means, sds & species:
  if (any(is.na(means)) | any(is.na(sds)) | any(is.na(species))) {
    stop("NA values are not allowed in 'means', 'sds' or 'species'")
  }
  #	5. Compute the species or populations upon which calculations will be done:
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1){
      type <- "One population_One species"
    } else{
      type <- "One population_Multiple species"
    }
  } else {
    if (length(samples) != nrow(means)) {
      stop("The length of 'samples' does not match the number of rows of 'means'
         and 'sds'")
    }
    if (any(is.na(samples))) {
        stop("NA values are not allowed in 'samples'")
    }
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1){
      type <- "Multiple populations_One species"
    } else{
      type <- "Multiple populations_Multiple species"
    }
  }

  #	6. Define trait ranges:
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (5, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      max_aux <- max(means[, dimens] + trait_ranges_aux[dimens] * sds[, dimens])
      min_aux <- min(means[, dimens] - trait_ranges_aux[dimens] * sds[, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  #	6. Create the grid of cells in which the density function is evaluated:
  if (is.null(n_divisions)) {
    n_divisions_choose<- c(1000, 200, 50, 25)
    n_divisions<- n_divisions_choose[dimensions]
  }
  grid_evaluate<-list()
  edge_length <- list()
  cell_volume<-1
  for (dimens in 1:dimensions){
    grid_evaluate[[dimens]] <- seq(from = trait_ranges[[dimens]][1],
                                   to = trait_ranges[[dimens]][2],
                                   length=n_divisions)
    edge_length[[dimens]] <- grid_evaluate[[dimens]][2] -
      grid_evaluate[[dimens]][1]
    cell_volume <- cell_volume * edge_length[[dimens]]
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  if (is.null(colnames(means))){
    names(evaluation_grid) <- paste0("Trait.",1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(means)
  }
  if (dimensions==1){
    evaluation_grid <- as.matrix(evaluation_grid)
  }
  # Creation of lists to store results:
  results <- list()
  # DATA: To store data and common features
  results$data <- list()
  results$data$evaluation_grid <- evaluation_grid
  results$data$cell_volume <- cell_volume
  results$data$edge_length <- edge_length
  results$data$species <- species
  results$data$means <- means
  results$data$sds <- sds
  if (is.null(samples)){
    results$data$populations <-  NA
  } else{
    results$data$populations <-  species_base
  }

  results$data$alpha <- alpha
  results$data$pop_means <- list()
  results$data$pop_sds <- list()
  results$data$pop_sigma <- list()
  results$data$dimensions <- dimensions
  results$data$type <- type
  results$data$method <- "mean"

  # TPDs: To store TPDs features of each species/population
  results$TPDs<-list()


  ########Multivariate normal density calculation
  for (spi in 1:length(unique(species_base))) {
    # Some information messages
    if (spi == 1) { message(paste0("-------Calculating densities for ", type, "-----------\n")) }
    #Data selection
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_means[[spi]] <- means[selected_rows, ]
    results$data$pop_sds[[spi]] <- sds[selected_rows, ]
    names(results$data$pop_means)[spi] <- names(results$data$pop_sds)[spi]<-
      unique(species_base)[spi]
    if (dimensions > 1) {
      results$data$pop_sigma[[spi]] <- diag(results$data$pop_sds[[spi]]^2)
      if (covar == TRUE) {
        for (tri in 1:dimensions){
          for (trj in 1:dimensions){
            if (tri != trj){
              results$data$pop_sigma[[spi]][tri, trj] <- results$data$pop_sigma[[spi]][trj, tri] <-
                results$data$pop_sds[[spi]][tri] * results$data$pop_sds[[spi]][trj] *
                stats::cor(means)[tri, trj]
            }
          }
        }
      }
      results$TPDs[[spi]] <- mvtnorm::dmvnorm(x = evaluation_grid,
                                              mean = results$data$pop_means[[spi]],
                                              sigma = results$data$pop_sigma[[spi]])
    }
    if (dimensions == 1) {
      results$data$pop_sigma[[spi]] <- results$data$pop_sds[[spi]]
      results$TPDs[[spi]] <- stats::dnorm(x = evaluation_grid,
                                   mean = results$data$pop_means[[spi]],
                                   sd = results$data$pop_sigma[[spi]])
    }
  }
  names(results$TPDs) <- unique(species_base)
  check_volume <- function(x) sum(x * cell_volume)
  volumes_checked <- sapply(results$TPDs, check_volume )
  if (any(abs(volumes_checked - 1) > tolerance)) {
    names_fail <- unique(species_base)[ which(abs(volumes_checked - 1) > tolerance)]
    message("Be careful, the integral of the pdf of some cases differ from 1.
        They have been reescaled, but you should consider increasing
        'trait_ranges' \n", paste(names_fail,collapse=" / "))
  }
  rescale_estimate <- function(x){
    x<- x / sum(x)
  }
  results$TPDs<- lapply(results$TPDs, rescale_estimate )
# Now, we extract the selected fraction of volume (alpha), if necessary
  extract_alpha <- function(x){
    # 1. Order the 'cells' according to their probabilities:
    alphaSpace_aux <- x[order(x, decreasing=T)]
    # 2. Find where does the accumulated sum of ordered probabilities becomes
    #   greater than the threshold (alpha):
    greater_prob <- alphaSpace_aux[which(cumsum(alphaSpace_aux ) > alpha) [1]]
    # 3. Substitute smaller probabilities with 0:
    x[x < greater_prob] <- 0
    # 5. Finally, reescale, so that the total density is finally 1:
    x <- x / sum(x)
    return(x)
  }
  if (alpha < 1) {
    results$TPDs<- lapply(results$TPDs, extract_alpha )
  }
  names(results$TPDs) <- unique(species_base)
  class(results) <- "TPDsp"
  return(results)
}
