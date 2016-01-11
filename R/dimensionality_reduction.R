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
reduce.dimensionality <- function(dist, ndim, rescale=TRUE) {
  # input check
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  if (!is.finite(ndim) || round(ndim) != ndim || length(ndim) != 1 || ndim < 1 || ndim >= nrow(dist))
    stop(sQuote("ndim"), " must be a whole number and 1 <= ndim <= nrow(dist)-1")

  space <- cmdscale(dist, k = ndim)
  if (rescale) space <- rescale.and.center(space)
  colnames(space) <- paste("Comp", seq_len(ncol(space)), sep="")
  space
}

#' @title Scaling and centering of matrix-like objects
#'
#' @description \code{rescale.and.center} uniformily scales a given matrix such that the returned space is centered on \code{center}, and each column was scaled equally such that the range of each column is at most \code{max.range}.
#'
#' @usage
#' rescale.and.center(x, center=0, max.range=1)
#'
#' @param x A numeric matrix or data frame.
#' @param center The new center point of the data.
#' @param max.range The maximum range of each column.
#'
#' @return The centered, scaled matrix. ZThe numeric centering and scalings used are returned as attributes.
#'
#' @export
#'
#' @examples
#' ## Generate a matrix from a normal distribution with a large standard deviation, and approximately centered at c(5, 5)
#' x <- matrix(rnorm(200*2, sd = 10, mean = 5), ncol=2)
#'
#' ## Center the dataset at c(0, 0) with a minimum of c(-.5, -.5) and a maximum of c(.5, .5)
#' x.scaled <- rescale.and.center(x, center=0, max.range=1)
#'
#' ## Plot rescaled data
#' plot(x.scaled)
#'
#' ## Show ranges of each column
#' apply(x.scaled, 2, range)
rescale.and.center <- function(x, center=0, max.range=1) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.finite(center))
    stop(sQuote("center"), " must be numeric")
  if (length(center) != 1 && length(center) != ncol(x))
    stop("length of ", sQuote("center"), " must be either 1 or ncol(x)")
  if (!is.finite(max.range))
    stop(sQuote("max.range"), " must be numeric")
  if (length(max.range) != 1)
    stop("length of ", sQuote("max.range"), " must be exactly 1")

  # calculate the minimum values of each column
  mins <- apply(x, 2, min)

  # calculate the maximum values of each column
  maxs <- apply(x, 2, max)

  # calculate the old center point
  old.center <- (maxs + mins) / 2
  new.center <- old.center - center

  # calculate the scale
  scale <- max(maxs - mins)

  # calculate rescaled data
  for (i in seq_len(nrow(x))) {
    x[i,] <- (x[i,] - new.center) / scale
  }

  # attach scaling information to output
  attr(x, "center") <- new.center
  attr(x, "scale") <- scale

  # return output
  x
}
