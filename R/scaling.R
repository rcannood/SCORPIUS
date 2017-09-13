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
  cur_range <- max(apply(ranges, 2, diff))

  if (cur_range != 0) {
    multiplier <- max_range / cur_range
  } else {
    multiplier <- 1
  }

  addend <- center / multiplier - colMeans(ranges)

  # calculate rescaled data
  apply_scale(x, addend, multiplier)
}

#' Apply a uniform scale
#'
#' @param x A numeric matrix or data frame.
#' @param addend A centering vector for each column
#' @param multiplier A scaling vector for each column
#'
#' @return The centered, scaled matrix. The numeric centering and scalings used are returned as attributes.
#' @export
apply_scale <- function(x, addend, multiplier) {
  y <- x %>%
    sweep(2, addend, "+") %>%
    sweep(2, multiplier, "*")
  attr(y, "addend") <- addend
  attr(y, "multiplier") <- multiplier
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

  addend <- -quants[1,]
  divisor <- apply(quants, 2, diff)
  divisor[divisor == 0] <- 1

  apply_quant_scale(x, addend, 1 / divisor)
}

#' Apply a quantile scale.
#'
#' Anything outside the range of [0, 1] will be set to 0 or 1.
#'
#' @param x A numeric matrix or data frame.
#' @param addend A minimum vector for each column
#' @param multiplier A scaling vector for each column
#'
#' @return The scaled matrix. The numeric centering and scalings used are returned as attributes.
#' @export
apply_quant_scale <- function(x, addend, multiplier) {
  y <- apply_scale(x, addend, multiplier)
  y[y > 1] <- 1
  y[y < 0] <- 0
  y
}

