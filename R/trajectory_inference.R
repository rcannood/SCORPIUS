#' @title Infer an initial trajectory through space
#'
#' @description \code{infer_initial_trajectory} infers an initial trajectory for
#' \code{\link{infer_trajectory}} by clustering the points and calculating
#' the shortest path through cluster centers. The shortest path takes into account
#' the euclidean distance between cluster centers, and the density between those two
#' points.
#'
#' @param space A numeric matrix or a data frame containing the coordinates of samples.
#' @param k The number of clusters
#'
#' @return the initial trajectory obtained by this method
#'
#' @importFrom RANN nn2
#' @importFrom TSP TSP insert_dummy solve_TSP
#' @importFrom stats kmeans dist
#' @importFrom dynutils scale_minmax add_class
infer_initial_trajectory <- function(space, k) {
  # input checks
  check_numeric_matrix(space, "space", finite = TRUE)
  check_numeric_vector(k, "k", whole = TRUE, finite = TRUE, range = c(1, nrow(space) - 1), length = 1)

  # cluster space into max k clusters
  fit <- stats::kmeans(space, k)
  centers <- fit$centers

  # calculate the euclidean space between clusters
  eucl_dist <- as.matrix(stats::dist(centers))

  # calculate the densities along the straight lines between any two cluster centers
  i <- j <- NULL # satisfy r cmd check
  pts <-
    crossing(
      i = seq_len(k),
      j = seq_len(k),
      pct = seq(0, 1, length.out = 21)
    ) %>%
    filter(i < j)
  pts_space <- (1 - pts$pct) * centers[pts$i, ] + pts$pct * centers[pts$j, ]

  pts$dist <- rowMeans(RANN::nn2(space, pts_space, k = 10)$nn.dist)
  dendis <- pts %>% group_by(i, j) %>% summarise(dist = mean(dist)) %>% ungroup()

  density_dist <- matrix(0, nrow = k, ncol = k)
  density_dist[cbind(dendis$i, dendis$j)] <- dendis$dist
  density_dist[cbind(dendis$j, dendis$i)] <- dendis$dist

  # combine both distance matrices
  cluster_distances <- eucl_dist * density_dist

  # find the shortest path through all clusters
  tsp <- TSP::insert_dummy(TSP::TSP(cluster_distances))
  tour <- as.vector(TSP::solve_TSP(tsp))
  tour2 <- c(tour, tour)
  start <- min(which(tour2 == k + 1))
  stop <- max(which(tour2 == k + 1))
  best_ord <- tour2[(start + 1):(stop - 1)]

  # use this ordering as the initial curve
  init_traj <- centers[best_ord, , drop = FALSE]

  init_traj
}


#' Infer linear trajectory through space
#'
#' \code{infer_trajectory} infers a trajectory through samples in a given space in a four-step process:
#' \enumerate{
#'   \item Perform \emph{k}-means clustering
#'   \item Calculate distance matrix between cluster centers using a custom distance function
#'   \item Find the shortest path connecting all cluster centers using the custom distance matrix
#'   \item Iteratively fit a curve to the given data using principal curves
#' }
#'
#' @inheritParams princurve::principal_curve
#' @param space A numeric matrix or a data frame containing the coordinates of samples.
#' @param k The number of clusters to cluster the data into.
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
#' @importFrom dynutils scale_minmax
#'
#' @examples
#' ## Generate an example dataset and visualise it
#' dataset <- generate_dataset(num_genes = 500, num_samples = 1000, num_groups = 4)
#' space <- reduce_dimensionality(dataset$expression, ndim = 2)
#' draw_trajectory_plot(space, progression_group = dataset$sample_info$group_name)
#'
#' ## Infer a trajectory through this space
#' traj <- infer_trajectory(space)
#'
#' ## Visualise the trajectory
#' draw_trajectory_plot(space, path=traj$path, progression_group = dataset$sample_info$group_name)
infer_trajectory <- function(
  space,
  k = 4,
  thresh = .001,
  maxit = 10,
  stretch = 0,
  smoother = "smooth_spline",
  approx_points = 100
) {
  # input checks
  check_numeric_matrix(space, "space", finite = TRUE)

  init_traj <-
    if (k <= 1) {
      NULL
    } else {
      infer_initial_trajectory(space, k = k)
    }

  # iteratively improve this curve using principal_curve
  fit <- princurve::principal_curve(
    as.matrix(space),
    start = init_traj,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother,
    approx_points = approx_points,
    trace = FALSE,
    plot_iterations = FALSE
  )

  # construct final trajectory
  path <- fit$s[fit$ord, , drop = FALSE]
  dimnames(path) <- list(NULL, paste0("Comp", seq_len(ncol(path))))

  # construct timeline values
  time <- dynutils::scale_minmax(fit$lambda)

  # output result
  list(
    path = path,
    time = time
  ) %>% dynutils::add_class(
    "SCORPIUS::trajectory"
  )
}

#' Reverse a trajectory
#'
#' Since the direction of the trajectory is not specified, the ordering of a trajectory may be inverted using \code{reverse_trajectory}.
#'
#' @param trajectory A trajectory as returned by \code{\link{infer_trajectory}}.
#'
#' @return The same trajectory, but in the other direction.
#'
#' @seealso \code{\link{infer_trajectory}}
#'
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' ## Generate an example dataset and infer a trajectory through it
#' dataset <- generate_dataset(num_genes = 500, num_samples = 1000, num_groups = 4)
#' group_name <- dataset$sample_info$group_name
#' space <- reduce_dimensionality(dataset$expression, ndim = 2)
#' traj <- infer_trajectory(space)
#'
#' ## Visualise the trajectory
#' draw_trajectory_plot(space, group_name, path = traj$path)
#'
#' ## Reverse the trajectory
#' reverse_traj <- reverse_trajectory(traj)
#' draw_trajectory_plot(space, group_name, path = reverse_traj$path)
#'
#' plot(traj$time, reverse_traj$time, type = "l")
reverse_trajectory <- function(trajectory) {
  if (!is(trajectory, "SCORPIUS::trajectory"))
    stop(sQuote("trajectory"), " needs to be an object returned by infer_trajectory")
  trajectory$time <- 1 - trajectory$time
  trajectory$path <- trajectory$path[rev(seq_len(nrow(trajectory$path))), , drop = FALSE]
  trajectory
}
