

#' @title k Nearest Neighbour distances
#'
#' @description \code{knn_distances} returns the distances of the \emph{k} nearest neighbours of each sample.
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param k The maximum number of nearest neighbours to search.
#' @param self_loops \code{TRUE} if samples with the same index or name are allowed to be neighbours.
#'
#' @return A matrix containing the distances of the \emph{k} nearest neighbours of each sample.
#'
#' @export
#'
#' @examples
#' ## Calculate the kNN distances within a set of samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' dist <- dist(x)
#' knnd <- knn_distances(dist, 10)
#' plot(density(knnd))
#'
#' ## Calculate the kNN distances between two sets of samples
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- euclidean_distance(x, y)
#' knnd <- knn_distances(dist, 10)
#' plot(density(knnd))
knn_distances <- function(dist, k, self_loops=FALSE) {
  knn(dist, k, self_loops = self_loops)$distances
}

#' @title k Nearest Neighbour indices and distances
#'
#' @description \code{knn} returns the indices and distances of the \emph{k} nearest neighbours of each sample.
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param k The maximum number of nearest neighbours to search.
#' @param self_loops \code{TRUE} if samples with the same index or name are allowed to be neighbours.
#'
#' @return A list containing two matrices \code{indices} and \code{distances}
#'
#' @export
#'
#' @importFrom utils head
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
#' dist <- euclidean_distance(x, y)
#' knnd <- knn(dist, 10)
#' plot(density(knnd$distances))
knn <- function(dist, k, self_loops=FALSE) {
  requireNamespace("utils")

  # input checks
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  if (nrow(dist) < 2)
    stop(sQuote("dist"), " needs to consist of at least 2 rows")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 0)
    stop(sQuote("k"), " must be a whole number and k >= 1")
  if (!is.logical(self_loops) || length(self_loops) != 1)
    stop(sQuote("self_loops"), " must be a logical value")

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
  row_ix <- if (!is.null(rownames(dist))) rownames(dist) else seq_len(nrow(dist))
  col_ix <- if (!is.null(colnames(dist))) colnames(dist) else seq_len(ncol(dist))

  # fill matrix by sample
  if (self_loops) {
    for (i in row_ix) {
      indices[i,] <- utils::head(order(dist[i,]), K)
      knndist[i,] <- dist[i,indices[i,]]
    }
  } else {
    diag(dist) <- 0 # just to make sure
    for (i in row_ix) {
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
