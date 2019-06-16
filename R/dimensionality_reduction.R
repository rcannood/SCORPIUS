#' Dimensionality reduction
#'
#' \code{reduce_dimensionality} performs an eigenanalysis of the given dissimilarity matrix
#' and returns coordinates of the samples represented in an \code{ndim}-dimensional space.
#'
#' @param x a numeric matrix
#' @param dist the distance metric to be used; can be any of the metrics listed in [dynutils::calculate_distance()].
#' @param ndim the maximum dimension of the space which the data are to be represented in; must be in {1, 2, \ldots, n-1}.
#' @param num_landmarks the number of landmarks to be selected.
#' @param rescale A logical indicating whether or not the returned space should be rescaled and centered.
#'
#' @return A matrix containing the coordinates of each sample, represented in an \code{ndim}-dimensional space.
#'
#' @seealso [SCORPIUS]
#'
#' @export
#'
#' @importFrom stats cmdscale
#' @importFrom dynutils calculate_distance list_distance_methods
#' @importFrom dyndimred dimred_landmark_mds
#'
#' @examples
#' ## Generate an example dataset
#' dataset <- generate_dataset(num_genes = 500, num_samples = 1000, num_groups = 4)
#'
#' ## Reduce the dimensionality of this dataset
#' space <- reduce_dimensionality(dataset$expression, ndim = 2)
#'
#' ## Visualise the dataset
#' draw_trajectory_plot(space, progression_group = dataset$sample_info$group_name)
reduce_dimensionality <- function(
  x,
  dist,
  ndim = 3,
  num_landmarks = 1000,
  rescale = TRUE
) {
  # input check
  check_numeric_matrix(x, "x", finite = TRUE, sparse = TRUE)
  check_numeric_vector(ndim, "ndim", finite = TRUE, whole = TRUE, range = c(1, nrow(x)), length = 1)

  space <- dyndimred::dimred_landmark_mds(
    x = x,
    distance_method = dist,
    ndim = ndim,
    num_landmarks = num_landmarks
  )

  colnames(space) <- paste0("Comp", seq_len(ncol(space)))

  space
}
formals(reduce_dimensionality)$dist <- dynutils::list_distance_methods()
