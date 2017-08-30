#' @title Dimensionality reduction
#'
#' @description \code{reduce.dimensionality} performs an eigenanalysis of the given dissimilarity matrix and returns coordinates of the samples represented in an \code{ndim}-dimensional space.
#'
#' @usage
#' reduce.dimensionality(dist, ndim, rescale=TRUE)
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param ndim The number of dimensions in the new space.
#' @param rescale A logical indicating whether or not the returned space should be rescaled and centered.
#'
#' @return A matrix containing the coordinates of each sample, represented in an \code{ndim}-dimensional space.
#'
#' @seealso \code{\link{correlation.distance}}, \code{\link{euclidean.distance}}, \code{\link{rescale.and.center}}, \code{\link{draw.trajectory.plot}}
#'
#' @export
#'
#' @importFrom stats cmdscale
#'
#' @examples
#' ## Generate an example dataset
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#'
#' ## Reduce the dimensionality of this dataset
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#'
#' ## Visualise the dataset
#' draw.trajectory.plot(space, progression.group=dataset$sample.info$group.name)
reduce.dimensionality <- function(dist, ndim = 3, rescale = TRUE) {
  requireNamespace("stats")

  # input check
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  if (!is.finite(ndim) || round(ndim) != ndim || length(ndim) != 1 || ndim < 1 || ndim >= nrow(dist))
    stop(sQuote("ndim"), " must be a whole number and 1 <= ndim <= nrow(dist)-1")

  space <- stats::cmdscale(dist, k = ndim)
  if (rescale) space <- rescale.and.center(space)
  colnames(space) <- paste("Comp", seq_len(ncol(space)), sep="")
  space
}
