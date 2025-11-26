#' Dimensionality reduction
#'
#' \code{reduce_dimensionality} performs an eigenanalysis of the given dissimilarity matrix
#' and returns coordinates of the samples represented in an \code{ndim}-dimensional space.
#'
#' @param x a numeric matrix
#' @param dist the distance metric to be used; can be any of the metrics listed in [dynutils::calculate_distance()].
#' @param ndim the maximum dimension of the space which the data are to be represented in; must be in \eqn{[1, n - 1]}, with \eqn{n} the number of samples (rows) in \code{x}.
#' @param num_landmarks the number of landmarks to be selected.
#'
#' @return A matrix containing the coordinates of each sample, represented in an \code{ndim}-dimensional space.
#'
#' @seealso [SCORPIUS]
#'
#' @export
#'
#' @importFrom stats cmdscale
#' @importFrom lmds lmds
#'
#' @examples
#' ## Generate an example dataset
#' dataset <- generate_dataset(num_genes = 200, num_samples = 400, num_groups = 4)
#'
#' ## Reduce the dimensionality of this dataset
#' space <- reduce_dimensionality(dataset$expression, ndim = 2)
#'
#' ## Visualise the dataset
#' draw_trajectory_plot(space, progression_group = dataset$sample_info$group_name)
reduce_dimensionality <- function(
  x,
  dist = c("spearman", "pearson", "euclidean", "cosine", "manhattan"),
  ndim = 3,
  num_landmarks = 1000
) {
  # input check
  check_numeric_matrix(x, "x", finite = TRUE, sparse = TRUE)
  check_numeric_vector(ndim, "ndim", finite = TRUE, whole = TRUE, range = c(1, nrow(x)), length = 1)
  dist <- match.arg(dist)

  space <- lmds::lmds(
    x = x,
    distance_method = dist,
    ndim = ndim,
    num_landmarks = num_landmarks
  )

  colnames(space) <- paste0("Comp", seq_len(ncol(space)))

  space
}
