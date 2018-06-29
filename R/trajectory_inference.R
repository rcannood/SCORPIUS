#' @title Infer an initial trajectory through space
#'
#' @description \code{infer_initial_trajectory} infers an initial trajectory for
#' \code{\link{infer_trajectory}} by clustering the points and calculating
#' the shortest path through cluster centers. The shortest path takes into account
#' the euclidean distance between cluster centers, and the density between those two
#' points.
#'
#' @usage
#' infer_initial_trajectory(space, k)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param k The number of clusters to cluster the data into.
#'
#' @return the initial trajectory obtained by this method
#'
#' @export
#'
#' @importFrom TSP TSP insert_dummy solve_TSP
#' @importFrom stats kmeans dist
#'
#' @examples
#' ## Generate an example dataset and visualise it
#' dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)
#' space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim=2)
#' draw_trajectory_plot(space, progression_group = dataset$sample_info$group_name)
#'
#' ## Infer a trajectory through this space
#' init_traj <- infer_initial_trajectory(space, k = 4)
#'
#' ## Visualise the trajectory
#' draw_trajectory_plot(space, path = init_traj, progression_group = dataset$sample_info$group_name)
infer_initial_trajectory <- function(space, k) {
  # input checks
  if (is.data.frame(space))
    space <- as.matrix(space)
  if (!is.matrix(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 2)
    stop(sQuote("k"), " must be a whole number and k >= 2")

  # cluster space into k clusters
  kmeans_clust <- stats::kmeans(space, centers = k)
  centers <- kmeans_clust$centers

  # calculate the euclidean space between clusters
  eucl_dist <- as.matrix(stats::dist(centers))

  # calculate the densities along the straight lines between any two cluster centers
  density_dist <- sapply(seq_len(k), function(i) {
    sapply(seq_len(k), function(j) {
      if (i == j) {
        0
      } else {
        twocent <- centers[c(i,j),,drop=FALSE]
        segment_pts <- apply(twocent, 2, function(x) seq(x[[1]], x[[2]], length.out = 20))
        dists <- euclidean_distance(segment_pts, space)
        mean(knn_distances(dists, 10, self_loops=TRUE))
      }
    })
  })

  requireNamespace("TSP")
  # combine both distance matrices
  cluster_distances <- eucl_dist * density_dist

  # find the shortest path through all clusters
  tsp <- TSP::insert_dummy(TSP::TSP(cluster_distances))
  tour <- as.vector(TSP::solve_TSP(tsp))
  tour2 <- c(tour, tour)
  start <- min(which(tour2 == k+1))
  stop <- max(which(tour2 == k+1))
  best_ord <- tour2[(start+1):(stop-1)]

  # use this ordering as the initial curve
  init_traj <- centers[best_ord,,drop=FALSE]

  init_traj
}


#' @title Infer linear trajectory through space
#'
#' @description \code{infer_trajectory} infers a trajectory through samples in a given space in a four-step process:
#' \enumerate{
#'   \item Perform \emph{k}-means clustering
#'   \item Calculate distance matrix between cluster centers using a custom distance function
#'   \item Find the shortest path connecting all cluster centers using the custom distance matrix
#'   \item Iteratively fit a curve to the given data using principal curves
#' }
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param k The number of clusters to cluster the data into.
#' @param thresh \code{\link[princurve]{principal_curve}} parameter: convergence threshhold on shortest distances to the curve
#' @param maxit \code{\link[princurve]{principal_curve}} parameter: maximum number of iterations
#' @param stretch \code{\link[princurve]{principal_curve}} parameter: a factor by which the curve can be extrapolated when points are projected
#' @param smoother \code{\link[princurve]{principal_curve}} parameter: choice of smoother
#'
#' @return A list containing several objects:
#' \itemize{
#'   \item \code{path}: the trajectory obtained by principal curves.
#'   \item \code{time}: the time point of each sample along the inferred trajectory.
#' }
#'
#' @seealso \code{\link{reduce_dimensionality}}, \code{\link{draw_trajectory_plot}}
#'
#' @export
#'
#' @importFrom princurve principal_curve
#'
#' @examples
#' ## Generate an example dataset and visualise it
#' dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)
#' space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim=2)
#' draw_trajectory_plot(space, progression_group = dataset$sample_info$group_name)
#'
#' ## Infer a trajectory through this space
#' traj <- infer_trajectory(space)
#'
#' ## Visualise the trajectory
#' draw_trajectory_plot(space, path=traj$path, progression_group=dataset$sample_info$group_name)
infer_trajectory <- function(space, k = 4, thresh = .001, maxit = 10, stretch = 0, smoother = "smooth_spline") {
  # input checks
  if (is.data.frame(space))
    space <- as.matrix(space)
  if (!is.matrix(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")

  if (!is.null(k)) {
    # use a clustering and shortest path based approach to define an intiial trajectory
    init_traj <- infer_initial_trajectory(space, k)
  } else {
    init_traj <- NULL
  }

  # iteratively improve this curve using principal_curve
  fit <- princurve::principal_curve(
    space,
    start = init_traj,
    thresh = thresh,
    plot_iterations = FALSE,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother,
    trace = FALSE
  )

  # construct final trajectory
  path <- fit$s[fit$ord, , drop = FALSE]
  dimnames(path) <- list(NULL, paste0("Comp", seq_len(ncol(path))))

  # construct timeline values
  time <- fit$lambda
  time <- (time - min(time)) / (max(time) - min(time))
  names(time) <- rownames(space)

  # output result
  trajectory <- list(
    path = path,
    time = time
  )
  class(trajectory) <- c(class(trajectory), "SCORPIUS::trajectory")
  trajectory
}

#' @title Reverse a trajectory
#'
#' @description Since the direction of the trajectory is not specified, the ordering of a trajectory may be inverted using \code{reverse_trajectory}.
#'
#' @usage
#' reverse_trajectory(trajectory)
#'
#' @param trajectory A trajectory as returned by \code{\link{infer_trajectory}}.
#'
#' @return The same trajectory, but in the other direction.
#'
#' @seealso \code{\link{infer_trajectory}}
#'
#' @export
#'
#' @examples
#' ## Generate an example dataset and infer a trajectory through it
#' dataset <- generate_dataset(type="poly", num_genes=500, num_samples=1000, num_groups=4)
#' group_name <- dataset$sample_info$group_name
#' space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim=2)
#' traj <- infer_trajectory(space)
#'
#' ## Visualise the trajectory
#' draw_trajectory_plot(space, group_name, path=traj$path)
#'
#' ## Reverse the trajectory
#' reverse_traj <- reverse_trajectory(traj)
#' draw_trajectory_plot(space, group_name, path=reverse_traj$path)
#'
#' ## It's the same but reversed?!
#' plot(traj$time, reverse_traj$time, type="l")
reverse_trajectory <- function(trajectory) {
  if (! "SCORPIUS::trajectory" %in% class(trajectory))
    stop(sQuote("trajectory"), " needs to be an object returned by infer_trajectory")
  trajectory$time <- 1-trajectory$time
  trajectory$path <- trajectory$path[rev(seq_len(nrow(trajectory$path))),,drop=FALSE]
  trajectory
}
