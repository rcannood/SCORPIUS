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
  if (!is.matrix(space) && !is.data.frame(space))
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


  # requireNamespace("GA")
  # # define a fitness function for an ordering of the cluster centers according to the cluster.distances such that the path length will be minimised.
  # fitness <- function(ord, dist=cluster.distances) {
  #   -sum(mapply(ord[-1], ord[-length(ord)], FUN=function(i, j) dist[[i, j]]))
  # }
  #
  # # if k <= 7, it's easier to just check all permutations of seq_len(k),
  # # else a genetic algorithm is used.
  # if (k <= 7) {
  #   permutations <- function( x, prefix = c() ) {
  #     if(length(x) == 0 ) return(prefix)
  #     do.call(rbind, sapply(1:length(x), FUN = function(idx) permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
  #   }
  #   ords <- permutations(seq_len(k))
  #   fitnesses <- apply(ords, 1, fitness)
  #   best.ord <- ords[which.max(fitnesses),]
  # } else {
  #   ga.fit <- GA::ga(type="permutation", fitness=fitness, min=1, max=nrow(centers), monitor=function(...) { })
  #   best.ord <- ga.fit@solution[1,]
  # }

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
#' infer.trajectory(space, k)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param k The number of clusters to cluster the data into.
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
infer.trajectory <- function(space, k = 4) {
  requireNamespace("princurve")

  # input checks
  if (!is.matrix(space) && !is.data.frame(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")

  if (!is.null(k)) {
    # use a clustering and shortest path based approach to define an intiial trajectory
    init.traj <- infer.initial.trajectory(space, k)
  } else {
    init.traj <- NULL
  }

  # iteratively improve this curve using principal.curve
  fit <- princurve::principal.curve(space, start = init.traj, plot.true = F, trace = F, stretch = 0)

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

#' @title Infer a guided trajectory
#'
#' @description \code{infer.guided.trajectory} infers multiple
#' trajectories and uses the original space with the combined trajectory orderings
#' to infer a guided trajectory.
#'
#' @usage
#' infer.guided.trajectory(space, number = 10, ...)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param number The number of trajectories to infer
#' @param ... Parameters for \code{\link{infer.trajectory}}
#'
#' @return A list containing several objects:
#' \itemize{
#'   \item \code{path}: NULL
#'   \item \code{time}: the time point of each sample along the inferred trajectory.
#' }
#'
#' @seealso \code{\link{infer.trajectory}}, \code{\link{reduce.dimensionality}}, \code{\link{draw.trajectory.plot}}
#'
#' @importFrom stats cor
#'
#' @export
infer.guided.trajectory <- function(space, number = 10, ...) {
  requireNamespace("stats")

  times <- sapply(seq_len(number), function(zzz) {
    infer.trajectory(space, ...)$time
  })

  for (i in seq_len(number)) {
    if (stats::cor(times[,1], times[,i]) < 0) {
      times[,i] <- 1 - times[,i]
    }
  }

  space.times <- cbind(space, times)
  SCORPIUS::infer.trajectory(space.times)
}

#' @title Average orderings of multiple trajectories
#'
#' @description \code{infer.consensus.trajectory} infers multiple
#' trajectories with the \code{\link{infer.trajectory}} method, and aligns
#' and averages the orderings.
#'
#' @usage
#' infer.consensus.trajectory(space, number = 10, ...)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param number The number of trajectories to infer
#' @param ... Parameters for \code{\link{infer.trajectory}}
#'
#' @return A list containing several objects:
#' \itemize{
#'   \item \code{path}: NULL
#'   \item \code{time}: the time point of each sample along the inferred trajectory.
#' }
#'
#' @seealso \code{\link{infer.trajectory}}, \code{\link{reduce.dimensionality}}, \code{\link{draw.trajectory.plot}}
#'
#' @importFrom stats cor
#'
#' @export
infer.consensus.trajectory <- function(space, number = 10, ...) {
  requireNamespace("stats")

  times <- sapply(seq_len(number), function(zzz) {
    infer.trajectory(space, ...)$time
  })

  for (i in seq_len(number)) {
    if (stats::cor(times[,1], times[,i]) < 0) {
      times[,i] <- 1 - times[,i]
    }
  }

  time <- rowMeans(times)
  path <- space[order(time),,drop=F]

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

#' @title Extract modules of features
#'
#' @description \code{extract.modules} uses adaptive branch pruning to extract modules of features, which is typically done on the smoothed expression returned by \code{\link{gene.importances}}.
#'
#' @usage
#' extract.modules(x, ...)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param ... Extra parameters passed to Mclust
#'
#' @return A data frame containing meta-data for the features in \code{x}, namely the order in which to visualise the features in and which module they belong to.
#'
#' @seealso \code{\link{gene.importances}}
#'
#' @export
#'
#' @importFrom mclust Mclust
#' @importFrom stats as.dist hclust
#' @importFrom dplyr bind_rows
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' group.name <- dataset$sample.info$group.name
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space)
#' time <- traj$time
#' draw.trajectory.plot(space, path=traj$path, group.name)
#'
#' ## Select most important genes
#' gimp <- gene.importances(expression, traj$time, num.permutations = 0)
#' gene.sel <- gimp[1:50,]
#' expr.sel <- expression[,gene.sel$gene]
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract.modules(quant.scale(expr.sel))
#' draw.trajectory.heatmap(expr.sel, time, group.name, modules)
extract.modules <- function(x, ...) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")

  requireNamespace("mclust")
  requireNamespace("dplyr")
  requireNamespace("stats")

  feature.names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))

  # sigh.. mclust doesn't do well with requireNamespace
  mclustBIC <- mclust::mclustBIC

  # cluster with mclust
  labels <- mclust::Mclust(t(x), ...)$classification

  # hierarchically cluster the features
  dist <- correlation.distance(t(x))
  hcl <- stats::hclust(stats::as.dist(dist), method="average")

  # order features within one module according to a dimensionality reduction of the correlation distance
  modules <- dplyr::bind_rows(lapply(unique(labels), function(l) {
    ix <- which(labels==l)
    dimred <- reduce.dimensionality(dist[ix, ix, drop=F], ndim=1)
    data.frame(feature=feature.names[ix], index=ix, module=l, value=dimred[,1], stringsAsFactors = F, row.names = NULL)
  }))

  modules <- as.data.frame(modules[order(modules$module, modules$value),,drop=F])

  # return output
  modules[,c("feature", "index", "module")]
}

#' @title Calculate the importance of a feature
#'
#' @description Calculates the feature importance of each column in \code{x} in trying to predict the time ordering.
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory as returned by \code{\link{infer.trajectory}}.
#' @param num.permutations The number of permutations to test against for calculating the p-values (default: 10).
#' @param ntree The number of trees to grow (default: 10000).
#' @param mtry The number of variables randomly samples at each split (default: 1\% of features).
#' @param ... Extra parameters passed to randomForest.
#'
#' @return a data frame containing the importance of each feature for the given time line
#'
#' @importFrom randomForest randomForest
#' @importFrom pbapply pblapply
#' @export
#'
#' @examples
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' group.name <- dataset$sample.info$group.name
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space)
#' gene.importances(expression, traj$time, num.permutations = 0)
gene.importances <- function(x, time, num.permutations = 10, ntree = 10000, mtry = ncol(x) * .01, ...) {
  importance <- randomForest::randomForest(x, time, ntree = ntree, mtry = mtry)$importance[,1]
  if (num.permutations > 0) {
    perms <- unlist(pbapply::pblapply(seq_len(num.permutations), function(i) {
      randomForest::randomForest(x, sample(time), ntree = ntree, mtry = mtry)$importance[,1]
    }))
    pvalue <- sapply(importance, function(x) mean(x < perms))
  } else {
    pvalue <- rep(NA, length(importance))
  }
  df <- data.frame(gene = colnames(x), importance, pvalue, stringsAsFactors = F)
  df[order(df$importance, decreasing = T), , drop = F]
}
