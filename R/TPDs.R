#' Trait Probability Density of Populations
#'
#' \code{TPDs} computes the trait probability density functions (TPD) of populations. A TPD for each population is calculated using kernel density estimators around each trait value provided \code{TPDs} can be used for single or multiple traits (up to four traits at the present time).
#'
#' @param species A vector containing the names of the species whose TPD is to be calculated. It must have the same length than 'traits'; therefore, repeated names are allowed. NA values are not allowed.
#' @param traits A matrix or data.frame containing the trait values of each	individual of each species. Individuals are in rows and trait values in columns, with one column for each trait. NA values are not allowed.
#' @param samples A vector containing the identity of the sampling unit in which each individual was present. Defaults to NULL, in which case, a TPDs is calculated for each species. If it is not NULL, a TPDs is calculated for each combination of 'species x sampling unit'  (ie, a TPDs for each population). NA values are not allowed.
#' @param weight A vector containing the weights to apply to each observation. It can be useful when individuals differ in their biomass. Defaults to NULL, in which case, all individuals are given the same weight.
#' @param alpha A number between 0 and 1, indicating the proportion of the probability density function of each population to include. A value of 1 includes the whole density function, but may be sensitive to the presence of outliers. Defaults to 0.95.
#' @param trait_ranges A vector or a list indicating the range of trait values that will be considered in the calculations. If a vector is provided, each	element should indicate the percentage (0-Inf) by which the range of	each trait should be expanded in each direction. If a list is provided, it should contain the range (minimum and maximum) of trait values that will be considered. Each element of the vector or list corresponds with one trait. The order of the traits must be the same as the order of the columns in the 'traits' arguments. Defaults to	NULL, in which case, the ranges of all the traits are automatically calculated	by expanding the range of the columns in 'traits' by 15\% in each direction. Trait ranges that are too short may result in an inadequate characterization of TPDs.
#'	@param n_divisions The number of equal-length parts in which each trait should be divided to calculate the grid in which calculations are based. Note	the number of cells composing the grid increases exponentially as	dimensionality increases, which can result in  high computation times. Defaults to NULL, in which case one trait is divided into 1000 parts, two traits are divided into 200 parts (40,000 cells), three traits are divided into 50 parts (125,000 cells), and 4 traits are divided into 25 parts (390,625 cells).
#'	@param tolerance A number between 0 and 1, giving the admissible proportion of deviation from 1 in the integral of the TPDs of each population. Integrals can be lower than 1 when the extent of the evaluation grid is not enough to capture all the probability density function of the species. These problems  are usually solved by increasing 'trait_ranges'. When the absolute deviation is greater than 'tolerance', a warning message is produced, but the TPDs function does not fail. Defaults to 0.05.

#'	@return \code{TPDs} returns an object of class "TPDsp", which is a list containing the following components:
#'
#'	\emph{data:} A list containing information used to perform the calculations, including the coordinates --in trait space-- in which the TPD function has been evaluated, the volume --in trait units-- of each cell of the grid, the length of each edges of the cells of the grid, the original trait data, the names of the species, the name of the populations in case sample is not NULL, the alpha level specified by the user, the traits of the individuals of each population, and the type of TPDs calculated.
#'
#'  \emph{TPDs:} A list, with one element per species or population, containing the probability associated to each cell of the grid in which the trait space has been divided.
#'
#'  @examples
#'
#' # 1.  Compute the TPDs of three different species
#' traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
#' sp_iris <- iris$Species
#' example_sp <- TPDs(species = sp_iris, traits = traits_iris)
#'
#' # 2.  Two different populations of each species
#' samples_ex <- rep(c(rep(1, 25), rep(2, 25)), 3)
#' example_pop <- TPDs (species = sp_iris, traits = traits_iris,
#'      samples = samples_ex)
#'
#' @import ks
#'
#' @export
TPDs <- function(species, traits, samples = NULL, weight = NULL, alpha = 0.95,
  trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05){

  # INITIAL CHECKS:
  # 	1. Compute the number of dimensions (traits):
  traits <- as.matrix(traits)
  dimensions <- ncol(traits)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the
      number of dimensions")
  }
  #	2. NA's not allowed in traits & species:
  if (any(is.na(traits)) | any(is.na(species))) {
    stop("NA values are not allowed in 'traits' or 'species'")
  }
  #	3. Compute the species or populations upon which calculations will be done:
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1){
      type <- "One population_One species"
    } else{
      type <- "One population_Multiple species"
    }
    if (min(table(species_base)) <= dimensions ){
      non_good <- which(table(species_base) <= dimensions)
      stop("You must have more observations (individuals) than traits for all
        species. Consider removing species with too few observations.\n
        Check these species:",
        paste(names(non_good), collapse=" / "))
    }
  } else {
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1){
      type <- "Multiple populations_One species"
    } else{
      type <- "Multiple populations_Multiple species"
    }
  if (min(table(species_base)) <= dimensions ){
      non_good <- which(table(species_base) <= dimensions)
      stop("You must have more observations (individuals) than traits for all
        populations. \n Consider removing populations with too few observations
        or pooling populations of the same species together.\n
        Otherwise, consider using the TPDs2 function.
        Check these populations:",
        paste(names(non_good), collapse=" / "))
    }
  }
  #	4. Define trait ranges:
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (15, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      obs_range <- as.numeric(stats::dist(range(traits[, dimens])))
      min_aux <- min(traits[, dimens]) -
        (obs_range * trait_ranges_aux[dimens] / 100)
      max_aux <- max(traits[, dimens]) +
        (obs_range * trait_ranges_aux[dimens] / 100)
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  #	5. Create the grid of cells in which the density function is evaluated:
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
  if (is.null(colnames(traits))){
    names(evaluation_grid) <- paste0("Trait.",1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(traits)
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
  results$data$traits <- traits
  results$data$dimensions <- dimensions
  results$data$species <- species
  if (is.null(samples)){
    results$data$populations <-  NA
  } else{
    results$data$populations <-  species_base
  }
  if (is.null(weight)){
    results$data$weight <-  weight
  } else{
    results$data$weight <-  NA
  }
  results$data$alpha <- alpha
  results$data$pop_traits <- list()
  results$data$pop_weight <- list()
  results$data$type <- type
  results$data$method <- "kernel"

  # TPDs: To store TPDs features of each species/population
  results$TPDs<-list()

  ########KDE CALCULATIONS
  for (spi in 1:length(unique(species_base))) {
    # Some information messages
    if (spi == 1) { message(paste0("-------Calculating densities for ", type, "-----------\n")) }
    #Data selection
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_traits[[spi]] <- traits[selected_rows, ]
    results$data$pop_weight[[spi]] <- length(selected_rows) *
      weight[selected_rows] / sum(weight[selected_rows])
    names(results$data$pop_traits)[spi] <- unique(species_base)[spi]
  }
  if (is.null(weight)){
    results$TPDs <- lapply(results$data$pop_traits, ks::kde,
      eval.points=evaluation_grid)
  } else{
    for (spi in 1:length(unique(species_base))) {
      results$TPDs[[spi]] <- ks::kde(results$data$pop_traits[[spi]],
        eval.points=evaluation_grid, w = results$data$pop_weight[[spi]])
    }
  }
  check_volume <- function(x) sum(x$estimate * cell_volume)
  volumes_checked <- sapply(results$TPDs, check_volume )
  if (any(abs(volumes_checked - 1) > tolerance)) {
    names_fail <- unique(species_base)[ which(abs(volumes_checked - 1) > tolerance)]
    message("Be careful, the integral of the pdf of some cases differ from 1.
      They have been reescaled, but you should consider increasing
      'trait_ranges' \n", paste(names_fail,collapse=" / "))
  }
  rescale_estimate <- function(x){
    x$estimate<- x$estimate / sum(x$estimate)
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
