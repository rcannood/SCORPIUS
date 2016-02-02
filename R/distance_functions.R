#' @title (Pairwise) Euclidean distances between two sets of samples
#'
#' @description \code{euclidean.distance} calculates the (pairwise) Euclidean distances between one or two sets of samples.
#'
#' @usage
#' euclidean.distance(x, y)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param y \code{NULL} (default) or a numeric matrix or data frame with \emph{N} rows (one per sample) and \emph{P} columns (one per feature).
#'
#' @return An \emph{M}-by-\emph{M} (if \code{y} is \code{NULL}) or an \emph{M}-by-\emph{N} (otherwise) matrix containing the Euclidean distances between the given sets of samples.
#'
#' @export
#'
#' @examples
#' ## generate two matrices with 50 and 100 samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- euclidean.distance(x, y)
#'
#' ## compare with the standard dist function
#' dist2 <- as.matrix(dist(rbind(x, y)))[1:50, 51:150]
#' plot(dist, dist2)
euclidean.distance <- function (x, y=NULL) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")

  # if y is null, we can simply use the normal dist function
  if (is.null(y)) return(as.matrix(dist(x)))

  # more input checks
  if (!is.matrix(y) && !is.data.frame(y))
    stop(sQuote("y"), " must be a numeric matrix or data frame")
  if (ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")

  # initialise a matrix with NAs
  z <- matrix(NA, nrow = nrow(x), ncol = nrow(y), dimnames=list(rownames(x), rownames(y)))

  # fill matrix by column
  for (k in seq_len(nrow(y))) {
    z[,k] <- sqrt(colSums((t(x) - y[k,])^2))
  }

  # return distances
  return(z)
}

#' @title Correlation distance
#'
#' @description \code{correlation.distance} calculates the (pairwise) correlation distances between one or two sets of samples.
#'
#' @usage
#' correlation.distance(x, y=NULL, method="spearman")
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param y \code{NULL} (default) or a numeric matrix or data frame with \emph{N} rows (one per sample) and \emph{P} columns (one per feature).
#' @param method A character string indicating which correlation coefficient (or covariance) is to be computed. One of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#' @param use See \code{\link[stats]{cor}}.
#'
#' @return An \emph{M}-by-\emph{M} (if \code{y} is \code{NULL}) or an \emph{M}-by-\emph{N} (otherwise) matrix containing the correlation distances between the given sets of samples.
#'
#' @export
#'
#' @examples
#' ## Generate two matrices with 50 and 100 samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- correlation.distance(x, y, method="spearman")
#'
#' ## Compare with the standard correlation function
#' dist2 <- cor(t(x), t(y), method="spearman")
#' plot(dist, dist2)
correlation.distance <- function(x, y=NULL, method="spearman", use = "everything") {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.null(y) && !is.matrix(y) && !is.data.frame(y))
    stop(sQuote("y"), " must be NULL, a numeric matrix or a data frame")
  if (!is.null(y) && ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"))
  if (is.na(na.method))
    stop("invalid 'use' argument")
  method <- match.arg(method)

  # transpose if necessary
  x <- t(x)
  if (!is.null(y)) y <- t(y)

  # calculate and return correlation distance
  1 - (cor(x, y, method=method, use = use)+1)/2
}

#' @title k Nearest Neighbour distances
#'
#' @description \code{knn.distances} returns the distances of the \emph{k} nearest neighbours of each sample.
#'
#' @usage knn.distances(dist, k, self.loops=F)
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param k The maximum number of nearest neighbours to search.
#' @param self.loops \code{TRUE} if samples with the same index or name are allowed to be neighbours.
#'
#' @return A matrix containing the distances of the \emph{k} nearest neighbours of each sample.
#'
#' @export
#'
#' @examples
#' ## Calculate the kNN distances within a set of samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' dist <- dist(x)
#' knnd <- knn.distances(dist, 10)
#' plot(density(knnd))
#'
#' ## Calculate the kNN distances between two sets of samples
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- euclidean.distance(x, y)
#' knnd <- knn.distances(dist, 10)
#' plot(density(knnd))
knn.distances <- function(dist, k, self.loops=F) {
  knn(dist, k, self.loops = self.loops)$distances
}

#' @title k Nearest Neighbour indices and distances
#'
#' @description \code{knn} returns the indices and distances of the \emph{k} nearest neighbours of each sample.
#'
#' @usage knn(dist, k, self.loops=F)
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param k The maximum number of nearest neighbours to search.
#' @param self.loops \code{TRUE} if samples with the same index or name are allowed to be neighbours.
#'
#' @return A list containing two matrices \code{indices} and \code{distances}
#'
#' @export
#'
#' @examples
#' ## Calculate the kNN distances within a set of samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' dist <- dist(x)
#' knnd <- knn(dist, 10)
#' plot(density(knnd$distances))
#'
#' ## Calculate the kNN distances between two sets of samples
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- euclidean.distance(x, y)
#' knnd <- knn(dist, 10)
#' plot(density(knnd$distances))
knn <- function(dist, k, self.loops=F) {
  # input checks
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  if (nrow(dist) < 2)
    stop(sQuote("dist"), " needs to consist of at least 2 rows")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 0)
    stop(sQuote("k"), " must be a whole number and k >= 1")
  if (!is.logical(self.loops) || length(self.loops) != 1)
    stop(sQuote("self.loops"), " must be a logical value")

  # k can't be larger than nrow(dist)-1
  K <- min(k, nrow(dist)-1)

  # initialise matrices with NAs
  indices <- knndist <- matrix(
    NA,
    nrow = nrow(dist),
    ncol = K,
    dimnames=list(rownames(dist), if (K==0) c() else paste0("knn", seq_len(K)))
  )

  # use dimnames if possible
  row.ix <- if (!is.null(rownames(dist))) rownames(dist) else seq_len(nrow(dist))
  col.ix <- if (!is.null(colnames(dist))) colnames(dist) else seq_len(ncol(dist))

  # fill matrix by sample
  if (self.loops) {
    for (i in row.ix) {
      indices[i,] <- head(order(dist[i,]), K)
      knndist[i,] <- dist[i,indices[i,]]
    }
  } else {
    diag(dist) <- 0 # just to make sure
    for (i in row.ix) {
      indices[i,] <- head(order(dist[i,]), K+1)[-1]
      knndist[i,] <- dist[i,indices[i,]]
    }
  }

  # return KNN distances
  list(
    indices = indices,
    distances = knndist
  )
}

