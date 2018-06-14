#' Plotting Trait Probability Distributions
#'
#' \code{plotTPD} plots a TPD object (created with either TPDs ot TPDc) of 1 or 2 dimensions. In the 1-dimension case, \code{plotTPD} displays the trait in the x-axis and the probability associated to each trait value in the y-axis. In the 2-dimensions case, \code{plotTPD} displays traits in the x- and y-axes, and probabilities are indicated by a gradient of colors. The function yields a panel for each TPD calculated (one for each species or population in the case of TPDs and one for each community in the case of TPDc).
#'
#'    In the 2-dimensions case, \code{plotTPD} requires the packages \code{ggplot2} and \code{gridExtra} in order to work.
#'
#' @param TPD An object of class "TPDsp" or TPDcomm", generated with the \code{\link{TPDs}} or \code{\link{TPDc}} functions, respectively, containing the TPDs of the considered populations or species or the TPDc of the considered communities.
#' @param whichPlot A vector indicating the identity of the species, populations or communities to plots. Defaults to NULL, in which case, all cases are plotted.
#' @param nRowCol A vector with two integers indicating the number of rows and columns of the layout. The product of the two numbers must be greater or equal than the length of whichPlot, so that all plots can be included. Defaults to NULL, in which case the layout is automatically selected.
#' @param color1 The color used to fill the TPD in the 1-dimension case. Defaults to "grey60".
#' @param leg Logical, indicating whether a legend with the name of the species, population or community should be added in the 1-dimension case.
#' @param leg.text Vector, containing the names to be used in the legend of each plot in the 1-dimension case, or in the plots title in the 2-dimensions case. If provided, it must have the same length as whichPlot. Defaults to NULL, in which case the names of the species/populations or communities are used.
#' @param leg.pos Character, indicating the location of the legend in the plot. Possible values are: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Defaults to "topright".
#' @param leg.cex Numeric, indicating the character expansion factor relative to current par("cex") for the text of the legend.
#'
#' @examples
#' # 1.  Compute the TPDs of five different species. SP3 is in the center of
#' #   the trait space, and the rest of species in the corners
#' set.seed(1)
#' nind <- 10
#' species_ex <- c(rep("SP1",nind), rep("SP2",nind), rep("SP3",nind),
#'  rep("SP4",nind), rep("SP5",nind))
#' traits_ex <- data.frame(trait1 = c(rnorm(nind, 10, 3),
#'                                    rnorm(nind, 10, 3),
#'                                    rnorm(nind, 15, 3),
#'                                    rnorm(nind, 20, 3),
#'                                    rnorm(nind, 20, 3)),
#'                         trait2 = c(rnorm(nind, 10, 3),
#'                                    rnorm(nind, 20, 3),
#'                                    rnorm(nind, 15, 3),
#'                                    rnorm(nind, 10, 3),
#'                                    rnorm(nind, 20, 3)))
#' \dontrun{
#' species_TPDs_2D <- TPDs (species = species_ex, traits = traits_ex)
#' # Plot all species
#' plotTPD(species_TPDs_2D)
#' # Plot only species 3
#' plotTPD(species_TPDs_2D, whichPlot = 3)
#' }
#' #1 dimension case:
#' species_TPDs_1D <- TPDs (species = species_ex, traits = traits_ex$trait1)
#' plotTPD(species_TPDs_1D)
#'
#' #Now, let us plot communities (TPDc)
#' #2. three different communities with different abundances of each species
#' abundances_ex <- matrix(c(c(0.05, 0.05, 0.8,  0.05, 0.05),
#'                           c(0.9,  0,    0,    0,    0.1 ),
#'                           c(0,    0,    1,    0,    0   )),
#'                    ncol = 5, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
#'                      unique(species_ex)))
#' \dontrun{
#' example_TPDc_2D <- TPDc (TPDs = species_TPDs_2D, sampUnit = abundances_ex)
#' plotTPD(example_TPDc_2D)
#' }
#' example_TPDc_1D <- TPDc (TPDs = species_TPDs_1D, sampUnit = abundances_ex)
#' plotTPD(example_TPDc_1D)
#'
#' @import ggplot2
#' @import gridExtra
#'
#' @export
plotTPD<- function(TPD, whichPlot = NULL, nRowCol = NULL, color1 = "grey60",
  leg = TRUE, leg.text = NULL, leg.pos = "topright", leg.cex = 1) {
  requireNamespace("ggplot2", quietly = T)
  requireNamespace("gridExtra", quietly = T)
  dimensions <- TPD$data$dimensions
  if (dimensions > 2){
    stop("plotTPD only works for 1 or 2 dimensional trait spaces.")
  }
  if (is.null(whichPlot)) { #if whichPlot is NULL, plot everything
    if (class(TPD) == "TPDsp")  whichPlot <- 1:length(TPD$TPDs)
    if (class(TPD) == "TPDcomm")  whichPlot <- 1:length(TPD$TPDc$TPDc)
  }
  if (class(TPD) == "TPDsp") {
    TPD_probs <- TPD$TPDs[whichPlot]
    TPD_grid <- TPD$data$evaluation_grid
    if(is.null(leg.text)){
      names_plots <- names(TPD$TPDs)[whichPlot]
    } else {
      if(length(leg.text) != length(whichPlot)){
        stop("leg.text and whichPlot must have the same length")
      }
      names_plots <- leg.text
    }
  }
  if (class(TPD) == "TPDcomm") {
    TPD_probs <- TPD$TPDc$TPDc[whichPlot]
    TPD_grid <- TPD$data$evaluation_grid
    if(is.null(leg.text)){
      names_plots <- names(TPD$TPDc$TPDc)[whichPlot]
    } else {
      if(length(leg.text) != length(whichPlot)){
        stop("leg.text and whichPlot must have the same length")
      }
      names_plots <- leg.text
    }  }
  if (is.null(nRowCol)) {
    nrowCol <- numeric()
    nRowCol[1] <- ceiling(sqrt(length(whichPlot)))
    nRowCol[2] <- ceiling(length(whichPlot) / nRowCol[1])
  }
  if (nRowCol[1]*nRowCol[2] < length(whichPlot)) {
    stop("You have to specify a higher number of columns and/or rows to plot
        all your communities")
  }
  graphics::par(mfrow = nRowCol, mar = c(0,0,0,0), oma = c(4,4,1,1))
  limSup <- 0
  limInf <- 100
  for (pli in 1:length(whichPlot)) {
    limSup <- max(limSup, max(TPD_probs[[pli]]))
    limInf <- min(limInf, min(TPD_probs[[pli]]))
  }
  if (dimensions == 1) {
    limx <- range(unlist(TPD_grid))
    for (pli in whichPlot) {
      axisy <- axisx <- FALSE
      if (pli %% (2*nRowCol[2]) == 1){
        axisy <- TRUE
      }
      Nbottom <- length(whichPlot) %% nRowCol[2]
      Nprebottom <- nRowCol[2] - Nbottom
      if (pli > (length(whichPlot) - (Nbottom + Nprebottom))) {
        axisx<-TRUE
      }
      #POLYGON
      graphics::plot(TPD_probs[[pli]], type = "n", xaxt = "n", yaxt = "n", xlab = "",
        ylab = "", ylim = c(0, limSup), lwd = 2, xlim = limx)
      graphics::polygon( x = c(unlist(TPD_grid), rev(unlist(TPD_grid))),
        y = c(TPD_probs[[pli]], rep(0, length(TPD_probs[[pli]]))),
        col = color1)
      if (leg) {
        graphics::legend(leg.pos, bty = "n", legend = names_plots[whichPlot[pli]], cex = leg.cex)
      }
      if (axisy) {
        graphics::axis(side = 2, at = c(0, round(limSup, 4)))
      }
      if (axisx) {
        graphics::axis(side = 1)
      }
      graphics::abline(h=0)
    }
  }
  if (dimensions == 2) {
    if (length(whichPlot) > 1) {
    message("Be patient, in the 2-dimensional case plots can take some time.\nIf you think it takes too long, consider reducing the number of plots using 'whichPlot'")
    }
    limx <- range(TPD_grid[, 1])
    limy <- range(TPD_grid[, 2])
    pl <- list()
    for (pli in 1:length(whichPlot)) {
      axisy <- axisx <- FALSE
      if (pli %% (2*nRowCol[2]) == 1) {
        axisy <- TRUE
      }
      if (pli >= length(whichPlot) -
          (nRowCol[2] -(length(whichPlot) %% nRowCol[2])) - 1) {
        axisx<-TRUE
      }
      mat <- cbind(TPD_grid, TPD_probs[[pli]])
      names(mat)<-c("T1", "T2", "prob")
      mat$prob[mat$prob == 0] <- NA
      pl[[pli]] <- with(mat, ggplot2::ggplot(mat, ggplot2::aes(T1, T2,  fill = prob),
        interpolate=T)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_gradient(na.value = "grey80") +
        ggplot2::xlim(limx[1], limx[2]) +
        ggplot2::ylim(limy[1], limy[2]) +
        ggplot2::ggtitle(names_plots[pli] ) +
        ggplot2::theme(
          axis.line = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(colour = "black", size = 8) ,
          axis.text.y = ggplot2::element_text(colour = "black", size = 8) ,
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          legend.position = "none",
          plot.title = ggplot2::element_text(size = 10),
          plot.background = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(size = 1, linetype = "solid" ,
            color = "black"),
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        )
    }
    do.call(gridExtra::grid.arrange, c(pl, ncol = nRowCol[2]))
  }
}
