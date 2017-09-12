#' @title Scaling and centering of matrix-like objects
#'
#' @description \code{rescale_and_center} uniformily scales a given matrix such that
#' the returned space is centered on \code{center}, and each column was scaled equally
#' such that the range of each column is at most \code{max_range}.
#'
#' @usage
#' rescale_and_center(x, center=0, max_range=1)
#'
#' @param x A numeric matrix or data frame.
#' @param center The new center point of the data.
#' @param max_range The maximum range of each column.
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
#' x_scaled <- rescale_and_center(x, center=0, max_range=1)
#'
#' ## Plot rescaled data
#' plot(x_scaled)
#'
#' ## Show ranges of each column
#' apply(x_scaled, 2, range)
rescale_and_center <- function(x, center = 0, max_range = 1) {
  ranges <- apply(x, 2, range)

  new_scale <- max(ranges[2,] - ranges[1,]) / max_range
  new_center <- (ranges[1,] + ranges[2,]) / 2 - (center * new_scale / 2)

  # calculate rescaled data
  apply_scale(x, new_center, new_scale)
}

#' Apply a uniform scale
#'
#' @param x A numeric matrix or data frame.
#' @param center A centering vector for each column
#' @param scale A scaling vector for each column
#'
#' @return The centered, scaled matrix. The numeric centering and scalings used are returned as attributes.
#' @export
apply_scale <- function(x, center, scale) {
  y <- t(apply(x, 1, function(x) (x - center) / scale))
  attr(y, "center") <- center
  attr(y, "scale") <- scale
  y
}

#' Calculate and apply a quantile scale
#'
#' @param x A numeric matrix or data frame.
#' @param outlier_cutoff The quantile cutoff for outliers (default 0.05).
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
#' x_scaled <- quant_scale(x)
#'
#' ## Plot rescaled data
#' plot(x_scaled)
#'
#' ## Show ranges of each column
#' apply(x_scaled, 2, range)
quant_scale <- function(x, outlier_cutoff = .05) {
  quants <- apply(x, 2, quantile, c(outlier_cutoff, 1 - outlier_cutoff), na.rm = T)

  center <- quants[1,]
  scale <- quants[2,] - quants[1,]
  scale[scale == 0] <- 1

  apply_quant_scale(x, center, scale)
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
apply_quant_scale <- function(x, center, scale) {
  y <- apply_scale(x, center, scale)
  y[y > 1] <- 1
  y[y < 0] <- 0
  y
}

