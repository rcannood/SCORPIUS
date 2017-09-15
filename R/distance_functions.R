#' @title (Pairwise) Euclidean distances between two sets of samples
#'
#' @description \code{euclidean_distance} calculates the (pairwise) Euclidean distances between one or two sets of samples.
#'
#' @usage
#' euclidean_distance(x, y)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param y \code{NULL} (default) or a numeric matrix or data frame with \emph{N} rows (one per sample) and \emph{P} columns (one per feature).
#'
#' @return An \emph{M}-by-\emph{M} (if \code{y} is \code{NULL}) or an \emph{M}-by-\emph{N} (otherwise) matrix containing the Euclidean distances between the given sets of samples.
#'
#' @importFrom stats dist
#' @useDynLib SCORPIUS
#'
#' @export
#'
#' @examples
#' ## generate two matrices with 50 and 100 samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- euclidean_distance(x, y)
#'
#' ## compare with the standard dist function
#' dist2 <- as.matrix(dist(rbind(x, y)))[1:50, 51:150]
#' plot(dist, dist2)
euclidean_distance <- function (x, y=NULL) {
  # input checks
  if (!any(c("matrix", "data.frame", "dgCMatrix") %in% class(x)))
    stop(sQuote("x"), " must be a numeric matrix, a data frame, or a Matrix::dgCMatrix")
  if (!is.null(y) && !any(c("matrix", "data.frame", "dgCMatrix") %in% class(y)))
    stop(sQuote("y"), " must be NULL, a numeric matrix, a data frame, or a Matrix::dgCMatrix")
  if (!is.null(y) && ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")

  # casting, just to make sure
  x <- as.matrix(x)
  if (!is.null(y)) y <- as.matrix(y)

  # if y is null, we can simply use the normal dist function
  if (is.null(y)) {
    return(as.matrix(stats::dist(x)))
  } else {
    euclidean_distance_rcpp(x, y)
  }
}

#' @title Correlation distance
#'
#' @description \code{correlation_distance} calculates the (pairwise) correlation distances between one or two sets of samples.
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param y \code{NULL} (default) or a numeric matrix or data frame with \emph{N} rows (one per sample) and \emph{P} columns (one per feature).
#' @param method A character string indicating which correlation coefficient (or covariance) is to be computed. One of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#' @param use See \code{\link[stats]{cor}}.
#'
#' @return An \emph{M}-by-\emph{M} (if \code{y} is \code{NULL}) or an \emph{M}-by-\emph{N} (otherwise) matrix containing the correlation distances between the given sets of samples.
#'
#' @importFrom stats cor
#'
#' @export
#'
#' @examples
#' ## Generate two matrices with 50 and 100 samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- correlation_distance(x, y, method="spearman")
#'
#' ## Compare with the standard correlation function
#' dist2 <- cor(t(x), t(y), method="spearman")
#' plot(dist, dist2)
correlation_distance <- function(x, y = NULL, method = c("spearman", "pearson", "kendall"), use = "everything") {
  # input checks
  if (!any(c("matrix", "data.frame", "dgCMatrix") %in% class(x)))
    stop(sQuote("x"), " must be a numeric matrix, a data frame, or a Matrix::dgCMatrix")
  if (!is.null(y) && !any(c("matrix", "data.frame", "dgCMatrix") %in% class(y)))
    stop(sQuote("y"), " must be NULL, a numeric matrix, a data frame, or a Matrix::dgCMatrix")
  if (!is.null(y) && ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"))
  if (is.na(na.method))
    stop("invalid 'use' argument")
  method <- match.arg(method)

  # transpose if necessary
  x <- t(as.matrix(x))
  if (!is.null(y)) y <- t(as.matrix(y))

  # calculate and return correlation distance
  1 - (stats::cor(x, y, method=method, use = use)+1)/2
}

