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
  if ((!is.matrix(x) && !is.data.frame(x)) || !is.numeric(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  
  # if y is null, we can simply use the normal dist function
  if (is.null(y)) return(as.matrix(dist(x)))
  
  # more input checks
  if ((!is.matrix(y) && !is.data.frame(y)) || !is.numeric(y))
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
#' correlation.distance(x, y)
#' 
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param y \code{NULL} (default) or a numeric matrix or data frame with \emph{N} rows (one per sample) and \emph{P} columns (one per feature).
#' @param method A character string indicating which correlation coefficient (or covariance) is to be computed. One of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
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
correlation.distance <- function(x, y=NULL, method="spearman") {
  # input checks
  if ((!is.matrix(x) && !is.data.frame(x)) || !is.numeric(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if ((!is.null(y) && !is.matrix(y) && !is.data.frame(y)) || ((is.matrix(y) || is.data.frame(y)) && !is.numeric(y)))
    stop(sQuote("y"), " must be NULL, a numeric matrix or a data frame")
  if (!is.null(y) && ncol(x) != ncol(y)) 
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  
  # transpose if necessary
  x <- t(x)
  if (!is.null(y)) y <- t(y)
  
  # calculate and return correlation distance
  1 - (cor(x, y, method=method)+1)/2
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
  # input checks
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (nrow(dist) < 2)
    stop(sQuote("dist"), " needs to consist of at least 2 rows")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 0)
    stop(sQuote("k"), " must be a whole number and k >= 1")
  if (!is.logical(self.loops) || length(self.loops) != 1)
    stop(sQuote("self.loops"), " must be a logical value")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  
  # k can't be larger than nrow(dist)-1
  K <- min(k, nrow(dist)-1)
  
  # initialise matrix with NAs
  z <- matrix(
    NA, 
    nrow = nrow(dist), 
    ncol = K, 
    dimnames=list(rownames(dist), if (K==0) c() else paste0("knn", seq_len(K)))
  )
  
  # use dimnames if possible
  row.ix <- if (!is.null(rownames(dist))) rownames(dist) else seq_len(nrow(dist))
  col.ix <- if (!is.null(colnames(dist))) colnames(dist) else seq_len(ncol(dist))
  
  # fill matrix by sample
  if (self.loops)
    for (i in row.ix)
      z[i,] <- head(sort(dist[i,]), K)
  else
    for (i in row.ix)
      z[i,] <- head(sort(dist[i,col.ix!=i]), K)
  
  # return KNN distances
  z
}

#' @title Outlierness
#' 
#' @description \code{outlierness} calculates the mean distance of each sample to its \emph{k} nearest neighbours.
#' 
#' @usage outlierness(dist, k=10)
#' 
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param k The maximum number of nearest neighbours to search.
#'
#' @return The outlierness values for each of the samples.
#' 
#' @export
#'
#' @examples
#' ## Generate example dataset
#' x <- matrix(rnorm(100*2, mean=0, sd=1), ncol=2)
#' dist <- dist(x)
#' outl <- outlierness(dist, 10)
#' # Visualise the outlierness scores for each of the points
#' plot(x, cex=outl, pch=20)
outlierness <- function(dist, k=10) {
  # input check
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 1)
      stop(sQuote("k"), " must be a whole number and k >= 1")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)  
  
  # calculate and return each sample outlierness'
  rowMeans(knn.distances(dist, k))
}

#' @title Outlier detection
#' 
#' @description \code{outlier.filter} calculates which samples are outliers by iteratively
#' removing the samples with the highest \emph{outlierness}' and fitting a normal distribution
#' to the remaining outlierness values. A selection of samples is made by picking the iteration
#' at which the log likelihood is maximised.
#'
#' @usage 
#' outlier.filter(dist)
#'
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#'
#' @return A boolean vector indicating whether samples are \emph{not} outliers.
#' 
#' @export
#'
#' @importFrom fitdistrplus fitdist
#' 
#' @examples
#' ## Generate normally distributed points, calculate their outliernesses and which points are outliers
#' x <- matrix(rnorm(200*2), ncol=2)
#' dist <- euclidean.distance(x)
#' filt <- outlier.filter(dist)
#' # plot points using their outlierness value as size and whether or not they were outliers as colours
#' plot(x, col=filt+2, cex=outlierness(dist)+1, pch=20)
#' # plot the score at each iteration of the removal process
#' likelihood.df <- attr(filt, "loglikelihood")
#' plot(likelihood.df$amount.removed, likelihood.df$log.likelihood, type="l")
#' 
#' ## Generate a random expression dataset
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=200, num.groups=4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' filt <- outlier.filter(dist)
#' # plot points using their outlierness value as size and whether or not they were outliers as colours
#' plot(space, col=filt+2, cex=outlierness(dist)+1, pch=20)
#' # plot the score at each iteration of the removal process
#' likelihood.df <- attr(filt, "loglikelihood")
#' plot(likelihood.df$amount.removed, likelihood.df$log.likelihood, type="l")
outlier.filter <- function(dist) {
  # input check
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  
  # make sure fitdistrplus is installed
  requireNamespace("fitdistrplus")
  
  # initialise data structures
  num.pts.removed <- nrow(dist)
  removed <- rep(NA, num.pts.removed+1)
  logliks <- rep(NA, num.pts.removed+1)
  ix <- seq_len(nrow(dist))
  
  # calculate log likelihood when no samples are removed
  filt <- rep(T, nrow(dist))
  outliernesses <- outlierness(dist)
  dist.fit <- fitdistrplus::fitdist(outliernesses, distr="norm")
  logliks[[1]] <- dist.fit$loglik
  
  # iteratively remove samples and calculate log likelihood of fits
  for (i in seq_len(num.pts.removed-3)) {
    removed.sample <- which(filt)[[which.max(outliernesses)]]
    removed[[i+1]] <- removed.sample
    filt <- !ix %in% removed
    outliernesses <- outlierness(dist[filt,filt,drop=F])
    
    tryCatch({
      dist.fit <- fitdistrplus::fitdist(outliernesses, distr="norm", keepdata=F)
      logliks[[i+1]] <- dist.fit$loglik
    }, error=function(e) {})
  }
  
  # finish up tail of execution
  remaining <- which(filt)[order(outliernesses, decreasing=T)]
  removed[seq(num.pts.removed-1, num.pts.removed+1)] <- remaining
  
  # return a vector indicating which samples are /not/ outliers
  filt <- !ix %in% removed[seq(1, which.max(logliks))]
  
  # attach log likelihood information as an attribute
  loglik.df <- data.frame(amount.removed=seq_along(removed)-1, sample.index=removed, log.likelihood=logliks)
  attr(filt, "loglikelihoods") <- loglik.df
  
  # return outlier output
  filt
}

#' @title Scaling and centering of matrix-like objects
#' 
#' @description \code{rescale.and.center} uniformily scales a given matrix such that the returned space is centered on \code{center}, and each column was scaled equally such that the range of each column is at most \code{max.range}.
#' 
#' @usage
#' rescale.and.center(x, center=0, max.range=1)
#' 
#' @param x A numeric matrix or data frame.
#' @param center The new center point of the data.
#' @param max.range The maximum range of each column.
#'
#' @return The centered, scaled matrix. ZThe numeric centering and scalings used are returned as attributes.
#' 
#' @export
#'
#' @examples
#' ## Generate a matrix from a normal distribution with a large standard deviation, and approximately centered at c(5, 5)
#' x <- matrix(rnorm(200*2, sd = 10, mean = 5), ncol=2)
#' 
#' ## Center the dataset at c(0, 0) with a minimum of c(-.5, -.5) and a maximum of c(.5, .5)
#' x.scaled <- rescale.and.center(x, center=0, max.range=1)
#' 
#' ## Plot rescaled data
#' plot(x.scaled)
#' 
#' ## Show ranges of each column
#' apply(x.scaled, 2, range) 
rescale.and.center <- function(x, center=0, max.range=1) {
  # input checks 
  if ((!is.matrix(x) && !is.data.frame(x)) || !is.numeric(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.finite(center))
    stop(sQuote("center"), " must be numeric")
  if (length(center) != 1 && length(center) != ncol(x)) 
    stop("length of ", sQuote("center"), " must be either 1 or ncol(x)")
  if (!is.finite(max.range))
    stop(sQuote("max.range"), " must be numeric")
  if (length(max.range) != 1) 
    stop("length of ", sQuote("max.range"), " must be exactly 1")
  
  # calculate the minimum values of each column
  mins <- apply(x, 2, min)
  
  # calculate the maximum values of each column
  maxs <- apply(x, 2, max)
  
  # calculate the old center point
  old.center <- (maxs + mins) / 2
  new.center <- old.center - center
  
  # calculate the scale
  scale <- max(maxs - mins)
  
  # calculate rescaled data
  for (i in seq_len(nrow(x))) {
    x[i,] <- (x[i,] - new.center) / scale
  }
    
  # attach scaling information to output
  attr(x, "center") <- new.center
  attr(x, "scale") <- scale
  
  # return output
  x
}

#' @title Dimensionality reduction
#' 
#' @description \code{reduce.dimensionality} performs an eigenanalysis of the given dissimilarity matrix and returns coordinates of the samples represented in an \code{ndim}-dimensional space.
#' 
#' @usage 
#' reduce.dimensionality(dist, ndim, rescale=T)
#' 
#' @param dist A numeric matrix, data frame or "\code{dist}" object. 
#' @param ndim The number of dimensions in the new space.
#' @param rescale A logical indicating whether or not the returned space should be rescaled and centered.
#'
#' @return A matrix containing the coordinates of each sample, represented in an \code{ndim}-dimensional space.
#' 
#' @seealso If \code{rescale} is TRUE, \code{\link{rescale.and.center}} is used to rescale and center the output space.
#' 
#' @export
#'
#' @examples
#' ## Generate an example dataset
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#' 
#' ## Reduce the dimensionality of this dataset
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' 
#' ## Visualise the dataset
#' plot.dimensionality.reduction(space, colour=dataset$sample.info$group.name)
reduce.dimensionality <- function(dist, ndim, rescale=T) {
  # input check
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (!is.finite(ndim) || round(ndim) != ndim || length(ndim) != 1 || ndim < 1 || ndim >= nrow(dist))
    stop(sQuote("ndim"), " must be a whole number and 1 <= ndim <= nrow(dist)-1")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  
  space <- cmdscale(dist, k = ndim)
  if (rescale) space <- rescale.and.center(space)
  colnames(space) <- paste("Comp", seq_len(ncol(space)), sep="")
  space
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
#' @param space A numeric matrix or data frame.
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
#' @export
#' 
#' @importFrom princurve principal.curve
#' @importFrom dplyr percent_rank
#'
#' @examples 
#' ## Generate an example dataset and visualise it
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' plot.dimensionality.reduction(space, colour=dataset$sample.info$group.name)
#' 
#' ## Infer a trajectory through this space
#' traj <- infer.trajectory(space, k=4)
#' 
#' ## Visualise the trajectory
#' plot.trajectory(space, traj$final.path, colour=dataset$sample.info$group.name)
infer.trajectory <- function(space, k) {
  requireNamespace("princurve")
  requireNamespace("dplyr")
  
  # input checks
  if ((!is.matrix(space) && !is.data.frame(space)) || !is.numeric(space))
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
        mean(knn.distances(dists, 10, self.loops=T))
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
#' @param trajectory A trajectory as returned by \code{infer.trajectory}.
#'
#' @return The same trajectory, but in the other direction.
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
#' plot.trajectory(space, traj$final.path, colour=dataset$sample.info$group.name)
#' 
#' ## Reverse the trajectory
#' reverse.traj <- reverse.trajectory(traj)
#' plot.trajectory(space, reverse.traj$final.path, colour=dataset$sample.info$group.name)
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
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory.
#' @param p.adjust.method The correction method used by \link{\code{p.adjust}}.
#' @param q.value.cutoff The cutoff used on the q-values.
#' @param df The degree of freedom used by \link{\code{ns}}.
#' @param parallel Must be \code{FALSE} (default), \code{TRUE} (will auto-detect number of cores), or a numeric value specifying the number of cores to use.
#'
#' @return Returns a list, containing the following items: \itemize{
#'   \item \code{tafs}: the names or indices of all trajectory-aligned features (q-value < \code{q.value.cutoff}),
#'   \item \code{p.values}: a data frame containing an adjusted p-value for each feature,
#'   \item \code{smooth.x}: the input matrix \code{x} smoothed by splines.
#' }
#' 
#' @export
#' 
#' @importFrom splines ns
#' @importFrom pbapply pblapply
#' @importFrom parallel mclapply detectcores
#' @importFrom dplyr arrange
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#' plot.trajectory(space, traj$final.path, colour=dataset$sample.info$group.name)
#' 
#' ## Show which genes are most trajectory aligned
#' tafs <- find.trajectory.aligned.features(expression, traj$time)
#' head(tafs$p.values, 10)
#' 
#' ## Visualise the expression of these genes in a heatmap
#' expr.tafs <- expression[,tafs$tafs]
#' plot.trajectory.heatmap(expr.tafs, time=traj$time, progression=dataset$sample.info$group.name)
find.trajectory.aligned.features <- function(x, time, p.adjust.method="BH", q.value.cutoff=1e-10, df=8, parallel=F) {
  # input checks
  if ((!is.matrix(x) && !is.data.frame(x)) || !is.numeric(x))
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
  p.values <- dplyr::arrange(p.values, q.value)
  
  # also return the names or indexes of the TAFs
  features <- p.values$feature[p.values$is.taf]
  
  # return all output
  list(tafs=features, p.values=p.values, smooth.x=smooth.x)
}

#' @title Extract modules of features
#' 
#' @description \code{extract.modules} uses adaptive branch pruning 
#' 
#' @usage 
#' extract.modules(x)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#'
#' @return A data frame containing meta-data for the features in \code{x}, namely the order in which to visualise the features in and which module they belong to.
#' 
#' @export
#' 
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom WGCNA mergeCloseModules
#' @importFrom dplyr bind_rows arrange
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#' plot.trajectory(space, traj$final.path, colour=dataset$sample.info$group.name)
#' 
#' ## Show which genes are most trajectory aligned
#' tafs <- find.trajectory.aligned.features(expression, traj$time)
#' expr.tafs <- expression[,tafs$tafs]
#' smooth.tafs <- tafs$smooth.x[,tafs$tafs]
#' 
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract.modules(smooth.tafs)
#' plot.trajectory.heatmap(expr.tafs, time=traj$time, progression=dataset$sample.info$group.name, modules=modules)
extract.modules <- function(x) {
  # input checks
  if ((!is.matrix(x) && !is.data.frame(x)) || !is.numeric(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  
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
  modules <- as.data.frame(dplyr::arrange(modules, module, value))
  
  # return output
  modules[,c("feature", "index", "module")]
}


#' Title
#'
#' @param space 
#' @param colour 
#' @param contour
#'
#' @return
#' @export
#' 
#' @import ggplot2
#' @importFrom MASS kde2d
#'
#' @examples
plot.dimensionality.reduction <- function(space, colour=NULL, contour=F) {
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2
  
  g <- ggplot2::ggplot() + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x="Component 1", y="Component 2", colour="Group", fill="Group") + 
    ggplot2::xlim(min-diff, max+diff) + 
    ggplot2::ylim(min-diff, max+diff) + 
    ggplot2::coord_equal(xlim=c(min-0.1*diff, max+0.1*diff), ylim=c(min-0.1*diff, max+0.1*diff))

  if (contour) {
    if (is.null(colour)) {
      g <- g + ggplot2::stat_density2d(ggplot2::aes(Comp1, Comp2), breaks=1, geom="polygon", data.frame(space), alpha=.1)
    } else {
      g <- g + ggplot2::stat_density2d(ggplot2::aes(Comp1, Comp2, fill=colour), breaks=1, geom="polygon", data.frame(space, colour), alpha=.1)
    }
  }
  
  if (is.null(colour)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(Comp1, Comp2), data.frame(space))
  } else {
    g <- g + ggplot2::geom_point(ggplot2::aes(Comp1, Comp2, colour=colour), data.frame(space, colour))
  }
  
  g
}

#' Title
#'
#' @param space 
#' @param path 
#' @param colour 
#'
#' @return
#' @export
#' 
#' @import ggplot2
#'
#' @examples
plot.trajectory <- function(space, path, colour=NULL, contour=F) {
  plot.dimensionality.reduction(space, colour=colour, contour=contour) + ggplot2::geom_path(ggplot2::aes(Comp1, Comp2), data.frame(path))
}

#' Title
#'
#' @param time 
#' @param progression 
#' @param colour 
#'
#' @return
#' @export
#' 
#' @import ggplot2
#'
#' @examples
plot.trajectory.density <- function(time, progression, colour=NULL) {
  ggplot2::ggplot() + ggplot2::geom_density(ggplot2::aes(time, colour=colour), data.frame(time, colour=factor(progression))) + ggplot2::theme_classic()
}

#' Title
#'
#' @param x 
#' @param time 
#' @param progression 
#' @param modules 
#' @param show.labels.row 
#' @param show.labels.col 
#' @param scale.features
#'
#' @return
#' @export
#'
#' @examples
plot.trajectory.heatmap <- function(x, time, progression=NULL, modules=NULL, show.labels.row=F, show.labels.col=F, scale.features=T, narrow.breaks=T) {
  col.ann <- data.frame(row.names = rownames(x), Time=time)
  
  x.part <- x[order(time),,drop=F]
  if (scale.features) x.part <- scale(x.part)
  x.part <- t(x.part)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  ann.col <- list(
    Time=c("#ca0020", "#f4a582", "#ffffff", "#bababa", "#404040")
  )
  
  if (!is.null(progression)) {
    if (!is.factor(progression)) progression <- factor(progression)
    col.ann$Progression <- progression
    ann.col$Progression <- setNames(gg_color_hue(length(levels(progression))), levels(progression))
  }
  
  if (!is.null(modules)) {
    x.part <- x.part[modules$index,]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  } else {
    gaps_row <- NULL
    cluster_rows <- T
  }
  
  labels_row <- if (!show.labels.row) rep("", nrow(x.part)) else NULL
  labels_col <- if (!show.labels.col) rep("", ncol(x.part)) else NULL
  
  if (scale.features && narrow.breaks) {
    break.cutoff <- quantile(abs(x.part), .9)
    breaks <- c(min(x.part), seq(-break.cutoff, break.cutoff, length.out=99), max(x.part))
  } else {
    breaks <- NA
  }
  
  pheatmap::pheatmap(
    x.part, 
    cluster_cols = F, 
    cluster_rows = cluster_rows,
    annotation_col = col.ann, 
    annotation_colors = ann.col, 
    gaps_row = gaps_row,
    breaks = breaks,
    labels_row = labels_row,
    labels_col = labels_col
  )
}


#' Title
#'
#' @param time 
#' @param progression 
#'
#' @return
#' @export
#'
#' @examples
evaluate.trajectory <- function(time, progression) {
  stime <- sort(time)
  
  # if some timepoints contain multiple cells, randomly order these
  diff <- stime[-1]-stime[-length(stime)]
  min.diff <- min(diff[diff!=0])
  noises <- runif(length(time), 0, .01) * min.diff
  noised.time <- time + noises
  rank <- rank(noised.time)
  
  consistencies <- t(sapply(seq_along(rank), function(i) {
    left <- rank[[i]] > rank
    right <- rank[[i]] < rank
    up <- progression[[i]] > progression
    down <- progression[[i]] < progression
    
    leftup <- sum(left & up)
    leftdown <- sum(left & down)
    rightup <- sum(right & up)
    rightdown <- sum(right & down)
    
    one <- (leftup + rightdown)
    two <- (rightup + leftdown)
    
    consistency.direction1 <- one / (one + two) * 2 - 1
    consistency.direction2 <- two / (one + two) * 2 - 1
    
    c(consistency.direction1, consistency.direction2)
  }))
  
  consistency <- colMeans(consistencies)
  
  max(consistency)
}


#' Title
#'
#' @param space 
#' @param progression 
#' @param k 
#'
#' @return
#' @export
#' 
#' @importFrom class knn
#'
#' @examples
evaluate.space <- function(space, progression, k=5) {
  pred <- sapply(seq_len(nrow(space)), function(i) {
    as.integer(as.character(class::knn(space[-i,,drop=F], space[i,,drop=F], progression[-i], k=k)))
  })
  mean(progression==pred)
}

#' @title Naive dataset generation
#' 
#' @description \code{generate.dataset} generates an expression dataset which can be used for visualisation purposes.
#' 
#' @usage
#' generate.dataset(type=c("splines", "polynomial"), num.samples=1000, num.genes=100, num.groups=4)
#'
#' @param type The type of function used in order to generate the expression data. Must be either \code{"splines"} (default) or \code{"polynomial"} (or abbreviations thereof).
#' @param num.samples The number of samples the dataset will contain.
#' @param num.genes The number of genes the dataset will contain.
#' @param num.groups The number of groups the samples will be split up in.
#'
#' @return A list containing the expression data and the meta data of the samples.
#' 
#' @export
#'
#' @examples
#' ## Generate a dataset
#' dataset <- generate.dataset(type="poly", num.genes=500, num.samples=1000, num.groups=4)
#' 
#' ## Reduce dimensionality and infer trajectory with SCORPIUS
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#' 
#' ## Visualise
#' plot.trajectory(space, traj$final.path, colour=dataset$sample.info$group.name)
generate.dataset <- function(type=c("splines", "polynomial"), num.samples=1000, num.genes=100, num.groups=4) {
  group.names <- paste0("Group ", seq_len(num.groups))
  gene.names <- paste0("Gene", seq_len(num.genes))
  sample.names <- paste0("Sample", seq_len(num.samples))
  
  x <- seq(-1, 1, length.out=num.samples)
  group <- cut(x, breaks=num.groups, labels = group.names)
  sample.info <- data.frame(name=sample.names, group.name=group)
  
  type <- match.arg(type, c("splines", "polynomial"))
  
  switch(type, polynomial={
    y <- poly(x, 2)
    sd <- .012 * sqrt(num.genes)
  }, splines={
    y <- splines::ns(x, df=3)
    sd <- .06 * sqrt(num.genes)
  })
  
  expression <- sapply(seq_len(num.genes), function(g) {
    scale <- rnorm(ncol(y), mean=0, sd=1)
    rowSums(sweep(y, 2, scale, "*")) + rnorm(length(x), sd=sd)
  })
  dimnames(expression) <- list(sample.names, gene.names)
  undetectable <- which(expression < 0)
  undetectable <- sample(undetectable, length(undetectable)*.5)
  throw.back <- -.4
  expression[expression < throw.back | seq_along(expression) %in% undetectable] <- throw.back
  expression <- (expression - throw.back) * 4
  
  list(expression=expression, sample.info=sample.info)
}
