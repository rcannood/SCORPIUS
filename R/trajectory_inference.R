
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
#'   \item \code{clustering}: the initial \emph{k}-means clustering.
#'   \item \code{initial.path}: the initial shortest path through the different cluster centers.
#'   \item \code{final.path}: the final trajectory obtained by principal curves.
#'   \item \code{time}: the time point of each sample along the inferred trajectory.
#' }
#'
#' @seealso \code{\link{reduce.dimensionality}}, \code{\link{draw.trajectory.plot}}
#'
#' @export
#'
#' @importFrom princurve principal.curve
#' @importFrom dplyr percent_rank
#' @importFrom GA ga
#'
#' @examples
#' ## Generate an example dataset and visualise it
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' draw.trajectory.plot(space, progression.group=dataset$sample.info$group.name)
#'
#' ## Infer a trajectory through this space
#' traj <- infer.trajectory(space, k=4)
#'
#' ## Visualise the trajectory
#' draw.trajectory.plot(space, path=traj$final.path, progression.group=dataset$sample.info$group.name)
infer.trajectory <- function(space, k) {
  requireNamespace("princurve")
  requireNamespace("dplyr")
  requireNamespace("GA")

  # input checks
  if (!is.matrix(space) && !is.data.frame(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 2)
    stop(sQuote("k"), " must be a whole number and k >= 2")

  # cluster space into k clusters
  kmeans.clust <- kmeans(space, centers = k)
  centers <- kmeans.clust$centers

  # calculate the euclidean space between clusters
  eucl.dist <- as.matrix(dist(centers))

  # calculate the densities along the straight lines between any two cluster centers
  density.dist <- sapply(seq_len(k), function(i) {
    sapply(seq_len(k), function(j) {
      if (i == j) {
        0
      } else {
        twocent <- centers[c(i,j),]
        unit <- twocent[2,] - twocent[1,]
        unit <- unit / sqrt(sum(unit^2)) * .01
        num.pts <- ceiling(mean((twocent[2,] - twocent[1,])/unit))
        segment.pts <- apply(twocent, 2, function(x) seq(x[[1]], x[[2]], length.out = num.pts))
        dists <- euclidean.distance(segment.pts, space)
        mean(knn.distances(dists, 10, self.loops=TRUE))
      }
    })
  })

  # combine both distance matrices
  cluster.distances <- eucl.dist * density.dist

  # define a fitness function for an ordering of the cluster centers according to the cluster.distances such that the path length will be minimised.
  fitness <- function(ord, dist=cluster.distances) {
    -sum(mapply(ord[-1], ord[-length(ord)], FUN=function(i, j) dist[[i, j]]))
  }

  # if k <= 7, it's easier to just check all permutations of seq_len(k),
  # else a genetic algorithm is used.
  if (k <= 7) {
    permutations <- function( x, prefix = c() ) {
      if(length(x) == 0 ) return(prefix)
      do.call(rbind, sapply(1:length(x), FUN = function(idx) permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
    }
    ords <- permutations(seq_len(k))
    fitnesses <- apply(ords, 1, fitness)
    best.ord <- ords[which.max(fitnesses),]
  } else {
    ga.fit <- GA::ga(type="permutation", fitness=fitness, min=1, max=nrow(centers), monitor=function(...) { })
    best.ord <- ga.fit@solution[1,]
  }

  # use this ordering as the initial curve
  initial.path <- centers[best.ord,,drop=F]

  # iteratively improve this curve using principal.curve
  fit <- princurve::principal.curve(space, start=initial.path, plot.true=F, trace=F, stretch = 0)

  # construct final trajectory
  final.path <- fit$s[fit$tag,,drop=F]
  colnames(final.path) <- paste0("Comp", seq_len(ncol(final.path)))
  rownames(final.path) <- NULL

  # construct timeline values
  time <- dplyr::percent_rank(order(fit$tag))
  names(time) <- rownames(space)

  # output result
  trajectory <- list(clustering=kmeans.clust, initial.path=initial.path, final.path=final.path, time=time)
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
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#'
#' ## Visualise the trajectory
#' draw.trajectory.plot(space, path=traj$final.path, progression.group=dataset$sample.info$group.name)
#'
#' ## Reverse the trajectory
#' reverse.traj <- reverse.trajectory(traj)
#' draw.trajectory.plot(space, path=reverse.traj$final.path, progression.group=dataset$sample.info$group.name)
#'
#' ## It's the same but reversed?!
#' plot(traj$time, reverse.traj$time, type="l")
reverse.trajectory <- function(trajectory) {
  if (class(trajectory) != "SCORPIUS::trajectory" && c("clustering", "initial.path", "final.path", "time") %in% names(trajectory))
    stop(sQuote("trajectory"), " needs to be an object returned by infer.trajectory")
  trajectory$time <- 1-trajectory$time
  trajectory$final.path <- trajectory$final.path[rev(seq_len(nrow(trajectory$final.path))),,drop=F]
  trajectory$initial.path <- trajectory$initial.path[rev(seq_len(nrow(trajectory$initial.path))),,drop=F]
  trajectory
}

#' @title Find trajectory-aligned features
#'
#' @description \code{find.trajectory.aligned.features} searches for features whose values . It evaluates how well splines are able model the values of a given feature over the timeline of a trajectory.
#'
#' @usage
#' find.trajectory.aligned.features(x, time, p.adjust.method="BH", q.value.cutoff=1e-10, df=8, parallel=F)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory as returned by \code{\link{infer.trajectory}}.
#' @param p.adjust.method The correction method used by \code{\link[stats]{p.adjust}}.
#' @param q.value.cutoff The cutoff used on the q-values.
#' @param df The degree of freedom used by \code{\link[splines]{ns}}.
#' @param parallel Must be \code{FALSE} (default), \code{TRUE} (will auto-detect number of cores), or a numeric value specifying the number of cores to use.
#'
#' @return Returns a list, containing the following items: \itemize{
#'   \item \code{tafs}: the names or indices of all trajectory-aligned features (q-value < \code{q.value.cutoff}),
#'   \item \code{p.values}: a data frame containing an adjusted p-value for each feature,
#'   \item \code{smooth.x}: the input matrix \code{x} smoothed by splines.
#' }
#'
#' @seealso \code{\link{infer.trajectory}}, \code{\link{draw.trajectory.heatmap}}
#'
#' @export
#'
#' @importFrom splines ns
#' @importFrom pbapply pblapply
#' @importFrom parallel mclapply detectCores
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#' draw.trajectory.plot(space, path=traj$final.path, progression.group=dataset$sample.info$group.name)
#'
#' ## Show which genes are most trajectory aligned
#' tafs <- find.trajectory.aligned.features(expression, traj$time)
#' head(tafs$p.values, 10)
#'
#' ## Visualise the expression of these genes in a heatmap
#' expr.tafs <- expression[,tafs$tafs]
#' draw.trajectory.heatmap(expr.tafs, time=traj$time, progression.group=dataset$sample.info$group.name)
find.trajectory.aligned.features <- function(x, time, p.adjust.method="BH", q.value.cutoff=1e-10, df=8, parallel=F) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.vector(time) || !is.numeric(time))
    stop(sQuote("time"), " must be a numeric vector")
  if (nrow(x) != length(time))
    stop(sQuote("time"), " must have one value for each row in ", sQuote("x"))
  if (!p.adjust.method %in% p.adjust.methods)
    stop(sQuote("p.adjust.method"), " must be one of: ", paste(p.adjust.methods, collapse=", "))
  if (!is.numeric(q.value.cutoff))
    stop(sQuote("q.value.cutoff"), " needs to be a numerical")
  if (!is.finite(df) || df < 1)
    stop(sQuote("df"), " needs to be a numerical and larger than 0")

  requireNamespace("splines")
  requireNamespace("pbapply")
  requireNamespace("parallel")

  # if 'parallel' is a number, use this as the number of cores
  # if 'parallel' is not a number nor a logical, throw an exception
  # if 'parallel' is TRUE, automatically detect the number of cores
  # if 'parallel' is FALSE, use pbapply
  lapply.fun <-
    if (is.numeric(parallel)) {
      function(...) parallel::mclapply(..., mc.cores=parallel)
    } else if (!is.logical(parallel)) {
      stop(sQuote("parallel"), " needs to be FALSE (default), TRUE (will auto-detect number of cores), or a numeric value specifying the number of cores to use")
    } else if (parallel) {
      function(...) parallel::mclapply(..., mc.cores=parallel::detectCores())
    } else {
      function(...) pbapply::pblapply(...)
    }

  # apply splines to time
  splt <- splines::ns(time, df)

  # perform anova test for each gene
  feature.results <- lapply.fun(seq_len(ncol(x)), function(i) {
    fit <- glm(x[,i]~splt, family=gaussian(), epsilon=1e-5)
    fit1 <- glm(x[,i]~1, family=gaussian(), epsilon=1e-5)
    test <- anova(fit1, fit, test = "F")
    pval <- test[6][2, 1]
    smooth <- predict(fit, splt)
    list(pval=pval, smooth=smooth)
  })

  # gather and adjust p-values
  pvals <- sapply(feature.results, function(x) x$pval)
  pvals[!is.finite(pvals)] <- 1
  qvals <- p.adjust(pvals, method=p.adjust.method, n=ncol(x))

  # gather smoothed values
  smooth.x <- sapply(feature.results, function(x) x$smooth)
  dimnames(smooth.x) <- dimnames(x)

  # categorise q-values
  qval.cat <- cut(
    qvals,
    breaks=c(1, .05, .001, 1e-5, 1e-10, 1e-20, 1e-40, -Inf),
    labels=c("p <= 1e-40", "1e-40 < p <= 1e-20", "1e-20 < p <= 1e-10", "1e-10 < p <= 1e-5", "1e-5 < p <= 0.001", "0.001 < p <= 0.05", "p > 0.05")
  )

  # select final set of tafs
  is.taf <- qvals < q.value.cutoff

  # construct output data frame and order by q.value
  feature.names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))
  p.values <- data.frame(feature=feature.names, p.value=pvals, q.value=qvals, is.taf=is.taf, category=qval.cat, stringsAsFactors = F)
  p.values <- p.values[order(p.values$q.value),,drop=F]

  # also return the names or indexes of the TAFs
  features <- p.values$feature[p.values$is.taf]

  # return all output
  list(tafs=features, p.values=p.values, smooth.x=smooth.x)
}

#' @title Extract modules of features
#'
#' @description \code{extract.modules} uses adaptive branch pruning to extract modules of features, which is typically done on the smoothed expression returned by \code{\link{find.trajectory.aligned.features}}.
#'
#' @usage
#' extract.modules(x)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#'
#' @return A data frame containing meta-data for the features in \code{x}, namely the order in which to visualise the features in and which module they belong to.
#'
#' @seealso \code{\link{find.trajectory.aligned.features}}
#'
#' @export
#'
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom WGCNA mergeCloseModules
#' @importFrom dplyr bind_rows
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#' draw.trajectory.plot(space, path=traj$final.path, progression.group=dataset$sample.info$group.name)
#'
#' ## Show which genes are most trajectory aligned
#' tafs <- find.trajectory.aligned.features(expression, traj$time)
#' expr.tafs <- expression[,tafs$tafs]
#' smooth.tafs <- tafs$smooth.x[,tafs$tafs]
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract.modules(smooth.tafs)
#' draw.trajectory.heatmap(expr.tafs, time=traj$time, progression.group=dataset$sample.info$group.name, modules=modules)
extract.modules <- function(x) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")

  requireNamespace("dynamicTreeCut")
  requireNamespace("WGCNA")
  requireNamespace("dplyr")

  feature.names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))

  # calculate the correlation distance between features
  dist <- correlation.distance(t(x))

  # hierarchically cluster the features
  hcl <- hclust(as.dist(dist), method="average")

  # perform adaptive branch pruning
  labels <- dynamicTreeCut::cutreeDynamic(hcl, distM=dist, cutHeight = 0.8, deepSplit=1, pamRespectsDendro = F, method="hybrid", verbose=0)

  # merge close modules
  labels <- WGCNA::mergeCloseModules(x, labels, cutHeight = 0.3, verbose=0)$colors
  labels <- match(labels, unique(labels))

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
