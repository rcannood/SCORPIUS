#' @title Infer an initial trajectory through space
#'
#' @description \code{infer.initial.trajectory} infers an initial trajectory for
#' \code{\link{infer.trajectory}} by clustering the points and calculating
#' the shortest path through cluster centers. The shortest path takes into account
#' the euclidean distance between cluster centers, and the density between those two
#' points.
#'
#' @usage
#' infer.initial.trajectory(space, k)
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
#' dataset <- generate.dataset(type = "poly", num.genes = 500, num.samples = 1000, num.groups = 4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim = 2)
#' draw.trajectory.plot(space, progression.group = dataset$sample.info$group.name)
#'
#' ## Infer a trajectory through this space
#' init.traj <- infer.initial.trajectory(space, k = 4)
#'
#' ## Visualise the trajectory
#' draw.trajectory.plot(space, path = init.traj, progression.group = dataset$sample.info$group.name)
infer.initial.trajectory <- function(space, k) {
  requireNamespace("TSP")
  requireNamespace("stats")

  # input checks
  if (is.data.frame(space))
    space <- as.matrix(space)
  if (!is.matrix(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 2)
    stop(sQuote("k"), " must be a whole number and k >= 2")

  # cluster space into k clusters
  kmeans.clust <- stats::kmeans(space, centers = k)
  centers <- kmeans.clust$centers

  # calculate the euclidean space between clusters
  eucl.dist <- as.matrix(stats::dist(centers))

  # calculate the densities along the straight lines between any two cluster centers
  density.dist <- sapply(seq_len(k), function(i) {
    sapply(seq_len(k), function(j) {
      if (i == j) {
        0
      } else {
        twocent <- centers[c(i,j),]
        segment.pts <- apply(twocent, 2, function(x) seq(x[[1]], x[[2]], length.out = 20))
        dists <- euclidean.distance(segment.pts, space)
        mean(knn.distances(dists, 10, self.loops=TRUE))
      }
    })
  })

  requireNamespace("TSP")
  # combine both distance matrices
  cluster.distances <- eucl.dist * density.dist

  # find the shortest path through all clusters
  tsp <- TSP::insert_dummy(TSP::TSP(cluster.distances))
  tour <- as.vector(TSP::solve_TSP(tsp))
  tour2 <- c(tour, tour)
  start <- min(which(tour2 == k+1))
  stop <- max(which(tour2 == k+1))
  best.ord <- tour2[(start+1):(stop-1)]

  # use this ordering as the initial curve
  init.traj <- centers[best.ord,,drop=F]

  init.traj
}


#' @title Infer linear trajectory through space
#'
#' @description \code{infer.trajectory} infers a trajectory through samples in a given space in a four-step process:
#' \enumerate{
#'   \item Perform \emph{k}-means clustering
#'   \item Calculate distance matrix between cluster centers using a custom distance function
#'   \item Find the shortest path connecting all cluster centers using the custom distance matrix
#'   \item Iteratively fit a curve to the given data using principal curves
#' }
#'
#' @usage
#' infer.trajectory(
#'   space,
#'   k = 4,
#'   thresh = .001,
#'   maxit = 10,
#'   stretch = 0,
#'   smoother = "smooth.spline"
#' )
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param k The number of clusters to cluster the data into.
#' @param thresh \code{\link[princurve]{principal.curve}} parameter: convergence threshhold on shortest distances to the curve
#' @param maxit \code{\link[princurve]{principal.curve}} parameter: maximum number of iterations
#' @param stretch \code{\link[princurve]{principal.curve}} parameter: a factor by which the curve can be extrapolated when points are projected
#' @param smoother \code{\link[princurve]{principal.curve}} parameter: choice of smoother
#'
#' @return A list containing several objects:
#' \itemize{
#'   \item \code{path}: the trajectory obtained by principal curves.
#'   \item \code{time}: the time point of each sample along the inferred trajectory.
#' }
#'
#' @seealso \code{\link{reduce.dimensionality}}, \code{\link{draw.trajectory.plot}}
#'
#' @export
#'
#' @importFrom princurve principal.curve
#'
#' @examples
#' ## Generate an example dataset and visualise it
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' draw.trajectory.plot(space, progression.group=dataset$sample.info$group.name)
#'
#' ## Infer a trajectory through this space
#' traj <- infer.trajectory(space)
#'
#' ## Visualise the trajectory
#' draw.trajectory.plot(space, path=traj$path, progression.group=dataset$sample.info$group.name)
infer.trajectory <- function(space, k = 4, thresh = .001, maxit = 10, stretch = 0, smoother = "smooth.spline") {
  requireNamespace("princurve")

  # input checks
  if (is.data.frame(space))
    space <- as.matrix(space)
  if (!is.matrix(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")

  if (!is.null(k)) {
    # use a clustering and shortest path based approach to define an intiial trajectory
    init.traj <- infer.initial.trajectory(space, k)
  } else {
    init.traj <- NULL
  }

  # iteratively improve this curve using principal.curve
  fit <- princurve::principal.curve(
    space,
    start = init.traj,
    thresh = thresh,
    plot.true = F,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother,
    trace = F
  )

  # construct final trajectory
  path <- fit$s[fit$tag,,drop=F]
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
  class(trajectory) <- "SCORPIUS::trajectory"
  trajectory
}

#' @title Reverse a trajectory
#'
#' @description Since the direction of the trajectory is not specified, the ordering of a trajectory may be inverted using \code{reverse.trajectory}.
#'
#' @usage
#' reverse.trajectory(trajectory)
#'
#' @param trajectory A trajectory as returned by \code{\link{infer.trajectory}}.
#'
#' @return The same trajectory, but in the other direction.
#'
#' @seealso \code{\link{infer.trajectory}}
#'
#' @export
#'
#' @examples
#' ## Generate an example dataset and infer a trajectory through it
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#' group.name <- dataset$sample.info$group.name
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space)
#'
#' ## Visualise the trajectory
#' draw.trajectory.plot(space, group.name, path=traj$path)
#'
#' ## Reverse the trajectory
#' reverse.traj <- reverse.trajectory(traj)
#' draw.trajectory.plot(space, group.name, path=reverse.traj$path)
#'
#' ## It's the same but reversed?!
#' plot(traj$time, reverse.traj$time, type="l")
reverse.trajectory <- function(trajectory) {
  if (class(trajectory) != "SCORPIUS::trajectory" && c("path", "time") %in% names(trajectory))
    stop(sQuote("trajectory"), " needs to be an object returned by infer.trajectory")
  trajectory$time <- 1-trajectory$time
  trajectory$path <- trajectory$path[rev(seq_len(nrow(trajectory$path))),,drop=F]
  trajectory
}
