

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
#' library(dynutils)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- calculate_distance(x, y, "euclidean")
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
#' library(dynutils)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- calculate_distance(x, y, "euclidean")
#' knnd <- knn(dist, 10)
#' plot(density(knnd$distances))
knn <- function(dist, k, self_loops = FALSE) {
  # input checks
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  check_numeric_matrix(dist, "dist", finite = TRUE)

  if (nrow(dist) < 2)
    stop(sQuote("dist"), " needs to consist of at least 2 rows")
  check_numeric_vector(k, "k", finite = TRUE, whole = TRUE, range = c(1, nrow(dist) - 1), length = 1)
  check_logical_vector(self_loops, "self_loops", length = 1)

  # initialise matrices with NAs
  indices <- knndist <- matrix(
    NA,
    nrow = nrow(dist),
    ncol = k,
    dimnames = list(
      rownames(dist),
      if (k == 0) c() else paste0("knn", seq_len(k))
    )
  )

  # use dimnames if possible
  row_ix <- if (!is.null(rownames(dist))) rownames(dist) else seq_len(nrow(dist))
  col_ix <- if (!is.null(colnames(dist))) colnames(dist) else seq_len(ncol(dist))

  # fill matrix by sample
  if (self_loops) {
    for (i in row_ix) {
      indices[i,] <- utils::head(order(dist[i,]), k)
      knndist[i,] <- dist[i,indices[i,]]
    }
  } else {
    diag(dist) <- 0 # just to make sure
    for (i in row_ix) {
      indices[i,] <- head(order(dist[i,]), k+1)[-1]
      knndist[i,] <- dist[i,indices[i,]]
    }
  }

  # return KNN distances
  list(
    indices = indices,
    distances = knndist
  )
}
