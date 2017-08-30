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
#' @return The centered, scaled matrix. The numeric centering and scalings used are returned as attributes.
#'
#' @export
#'
#' @examples
#' ## Generate a matrix from a normal distribution
#' ## with a large standard deviation, centered at c(5, 5)
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
rescale.and.center <- function(x, center = 0, max.range = 1) {
  ranges <- apply(x, 2, range)

  new.scale <- max(ranges[2,] - ranges[1,]) / max.range
  new.center <- (ranges[1,] + ranges[2,]) / 2 - (center * new.scale / 2)

  # calculate rescaled data
  out_x <- apply.scale(x, new.center, new.scale)
}

#' Apply scaling factors
#'
#' @param x A numeric matrix or data frame.
#' @param center A centering vector for each column
#' @param scale A scaling vector for each column
#'
#' @return The centered, scaled matrix. The numeric centering and scalings used are returned as attributes.
#' @export
apply.scale <- function(x, center, scale) {
  y <- t(apply(x, 1, function(x) (x - center) / scale))
  attr(y, "center") <- center
  attr(y, "scale") <- scale
  y
}

#' Calculate and apply a quantile scale
#'
#' @param x A numeric matrix or data frame.
#' @param outlier.cutoff The quantile cutoff for outliers (default 0.05).
#'
#' @return The centered, scaled matrix. The numeric centering and scalings used are returned as attributes.
#'
#' @export
#'
#' @importFrom stats quantile
#'
#' @examples
#' ## Generate a matrix from a normal distribution
#' ## with a large standard deviation, centered at c(5, 5)
#' x <- matrix(rnorm(200*2, sd = 10, mean = 5), ncol=2)
#'
#' ## Center the dataset at c(0, 0) with a minimum of c(-.5, -.5) and a maximum of c(.5, .5)
#' x.scaled <- quant.scale(x)
#'
#' ## Plot rescaled data
#' plot(x.scaled)
#'
#' ## Show ranges of each column
#' apply(x.scaled, 2, range)
quant.scale <- function(x, outlier.cutoff = .05) {
  quants <- apply(x, 2, quantile, c(outlier.cutoff, 1 - outlier.cutoff), na.rm = T)

  center <- quants[1,]
  scale <- quants[2,] - quants[1,]
  scale[scale == 0] <- 1

  apply.quant.scale(x, center, scale)
}

#' Apply a quantile scale.
#'
#' Anything outside the range of [0, 1] will be set to 0 or 1.
#'
#' @param x A numeric matrix or data frame.
#' @param center A centering vector for each column
#' @param scale A scaling vector for each column
#'
#' @return The centered, scaled matrix. The numeric centering and scalings used are returned as attributes.
#' @export
apply.quant.scale <- function(x, center, scale) {
  y <- apply.scale(x, center, scale)
  y[y > 1] <- 1
  y[y < 0] <- 0
  y
}

