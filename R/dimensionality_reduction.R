#' @title Dimensionality reduction
#'
#' @description \code{reduce_dimensionality} performs an eigenanalysis of the given dissimilarity matrix and returns coordinates of the samples represented in an \code{ndim}-dimensional space.
#'
#' @usage
#' reduce_dimensionality(dist, ndim, rescale=TRUE)
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param ndim The number of dimensions in the new space.
#' @param rescale A logical indicating whether or not the returned space should be rescaled and centered.
#'
#' @return A matrix containing the coordinates of each sample, represented in an \code{ndim}-dimensional space.
#'
#' @seealso \code{\link{correlation_distance}}, \code{\link{scale_uniform}}, \code{\link{draw_trajectory_plot}}
#'
#' @export
#'
#' @importFrom stats cmdscale
#'
#' @examples
#' ## Generate an example dataset
#' dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)
#'
#' ## Reduce the dimensionality of this dataset
#' dist <- correlation_distance(dataset$expression)
#' space <- reduce_dimensionality(dist, ndim = 2)
#'
#' ## Visualise the dataset
#' draw_trajectory_plot(space, progression_group=dataset$sample_info$group_name)
reduce_dimensionality <- function(dist, ndim = 3, rescale = TRUE) {
  # input check
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  if (!is.finite(ndim) || round(ndim) != ndim || length(ndim) != 1 || ndim < 1 || ndim >= nrow(dist))
    stop(sQuote("ndim"), " must be a whole number and 1 <= ndim <= nrow(dist)-1")

  space <- stats::cmdscale(dist, k = ndim)
  if (rescale) space <- scale_uniform(space)
  colnames(space) <- paste("Comp", seq_len(ncol(space)), sep = "")
  space
}
