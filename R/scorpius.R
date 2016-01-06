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
#' ## generate two matrices with 50 and 100 samples
#' x <- matrix(rnorm(50*10, mean=0, sd=1), ncol=10)
#' y <- matrix(rnorm(100*10, mean=1, sd=2), ncol=10)
#' dist <- correlation.distance(x, y, method="spearman")
#' 
#' ## compare with the standard correlation function
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
#' @usage knn.distances(dist, k)
#' 
#' @param dist A numeric matrix, data frame or "\code{dist}" object.
#' @param k The maximum number of nearest neighbours to search.
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
knn.distances <- function(dist, k) {
  # input checks
  if (!is.matrix(dist) && !is.data.frame(dist) && class(dist) != "dist")
    stop(sQuote("dist"), " must be a numeric matrix, data frame or a ", sQuote("dist"), " object")
  if (!is.finite(k) || round(k) != k || length(k) != 1 || k < 0)
    stop(sQuote("k"), " must be a whole number and k >= 1")
  if (class(dist) == "dist")
    dist <- as.matrix(dist)
  
  # k can't be larger than nrow(dist)-1
  K <- min(k, nrow(dist)-1)
  
  # initialise matrix with NAs
  z <- matrix(
    NA, 
    nrow = nrow(dist), 
    ncol = K, 
    dimnames=list(rownames(dist), paste0("knn", seq_len(K)))
  )
  
  # fill matrix by sample
  for (i in seq_len(nrow(z))) {
    z[i,] <- head(sort(dist[i,-i]), K)
  }
  
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
#' x <- matrix(rnorm(100*2, mean=0, sd=1), ncol=2)
#' dist <- dist(x)
#' outl <- outlierness(dist, 10)
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
#' ## generate uniformily distributed points, calculate their outliernesses and which points are outliers
#' x <- matrix(runif(200*2), ncol=2)
#' dist <- euclidean.distance(x)
#' filt <- outlier.filter(dist)
#' plot(x, col=filt+2, cex=outlierness(dist)+1, pch=20)
#' 
#' ## generate normally distributed points, calculate their outliernesses and which points are outliers
#' x <- matrix(rnorm(200*2), ncol=2)
#' dist <- euclidean.distance(x)
#' filt <- outlier.filter(dist)
#' plot(x, col=filt+2, cex=outlierness(dist)+.5, pch=20)
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
  for (i in seq_len(num.pts.removed)) {
    if (i < num.pts.removed) {
      removed.sample <- which(filt)[[which.max(outliernesses)]]
      removed[[i+1]] <- removed.sample
      filt <- !ix %in% removed
      outliernesses <- outlierness(dist[filt,filt,drop=F])
    } else {
      removed[[i+1]] <- which(filt)
    }
    
    # it doesn't make sense to fit a distribution to less than three values.
    if (i < (num.pts.removed-3)) {
      dist.fit <- fitdistrplus::fitdist(outliernesses, distr="norm")
      logliks[[i+1]] <- dist.fit$loglik
    } else {
      logliks[[i+1]] <- -Inf
    }
  }
  
  # return a vector indicating which samples are /not/ outliers
  filt <- !ix %in% removed[seq(2, which.max(logliks))]
  
  # attach log likelihood information as an attribute
  attr(filt, "loglikelihoods") <- data.frame(amount.removed=seq_along(removed)-1, sample.index=removed, log.likelihood=logliks)
  
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
#' x <- matrix(rnorm(200*2, sd = 10, mean = 5), ncol=2)
#' x.scaled <- rescale.and.center(x, center=0, max.range=1)
#' # x.scaled has a maximum range of 1 per column, and is centered at 0.
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
  rescaled <- t(apply(x, 1, function(row) (row-new.center)/scale))
    
  # attach scaling information to output
  attr(rescaled, "center") <- new.center
  attr(rescaled, "scale") <- scale
  
  # return output
  rescaled
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
#' x <- matrix(runif(1000*10), ncol=10)
#' dist <- euclidean.distance(x)
#' space <- reduce.dimensionality(dist, ndim=2, rescale = F)
#' plot(space)
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
#' ## generate example dataset and visualise it
#' x <- seq(-1, 1, by=.002)
#' group <- cut(x, breaks=4)
#' y <- rowSums(poly(x, 3))*4
#' space <- cbind(Comp1=x + rnorm(length(x), sd=.1), Comp2=y + rnorm(length(y), sd=.2))
#' space <- rescale.and.center(space)
#' plot.dimensionality.reduction(space, colour=group)
#' 
#' ## infer trajectory and visualise
#' traj <- infer.trajectory(space, k=4)
#' plot.trajectory(space, traj$final.path, colour=group)
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
        mean(knn.distances(dists, 10))
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

#' @title Infer random trajectory through space
#' 
#' @description \code{infer.random.trajectory} infers a random trajectory in the given space.
#' 
#' @usage 
#' infer.random.trajectory(space, k)
#' 
#' @param space A numeric matrix or data frame.
#' @param k The number of clusters to cluster the data into.
#'
#' @return A list of the same structure as \link{\code{infer.trajectory}}
#' 
#' @export
#' 
#' @importFrom princurve principal.curve
#' @importFrom dplyr percent_rank
#'
#' @examples 
#' ## generate example dataset and visualise it
#' x <- seq(-1, 1, by=.002)
#' group <- cut(x, breaks=4)
#' y <- rowSums(poly(x, 3))*4
#' space <- cbind(Comp1=x + rnorm(length(x), sd=.1), Comp2=y + rnorm(length(y), sd=.2))
#' space <- rescale.and.center(space)
#' plot.dimensionality.reduction(space, colour=group)
#' 
#' ## infer trajectory and visualise
#' traj <- infer.random.trajectory(space, k=10)
#' plot.trajectory(space, traj$final.path, colour=group)
infer.random.trajectory <- function(space, k) {
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
  
  best.ord <- sample(seq_len(k))
  
  # use this ordering as the initial curve
  initial.path <- centers[best.ord,,drop=F]
  
  # iteratively improve this curve using principal.curve
  fit <- princurve::principal.curve(space, start=initial.path, plot.true=F, trace=F, stretch = 0, maxit=0)
  
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
#' ## generate example dataset and visualise it
#' x <- seq(-1, 1, by=.002)
#' group <- cut(x, breaks=4)
#' y <- rowSums(poly(x, 3))*4
#' space <- cbind(Comp1=x + rnorm(length(x), sd=.1), Comp2=y + rnorm(length(y), sd=.2))
#' space <- rescale.and.center(space)
#' plot.dimensionality.reduction(space, colour=group)
#' 
#' ## infer trajectory and visualise
#' traj <- infer.trajectory(space, k=4)
#' plot.trajectory(space, traj$final.path, colour=group)
#' 
#' ## reverse this trajectory
#' reverse.traj <- reverse.trajectory(traj)
#' plot.trajectory(space, traj$final.path, colour=group)
#' 
#' ## it's the same but reversed?!
#' plot(traj$time, reverse.traj$time, type="l")
reverse.trajectory <- function(trajectory) {
  if (class(trajectory) != "SCORPIUS::trajectory" && c("clustering", "initial.path", "final.path", "time") %in% names(trajectory)) 
    stop(sQuote("trajectory"), " needs to be an object returned by infer.trajectory")
  trajectory$time <- 1-trajectory$time
  trajectory$final.path <- trajectory$final.path[rev(seq_len(nrow(trajectory$final.path))),,drop=F]
  trajectory$initial.path <- trajectory$initial.path[rev(seq_len(nrow(trajectory$initial.path))),,drop=F]
  trajectory
}

#' Title
#'
#' @param counts 
#' @param time 
#' @param degrees.of.freedom 
#' @param p.adjust.method
#'
#' @return
#' @export
#' 
#' @importFrom splines ns
#' @importFrom pbapply pblapply
#' @importFrom parallel mclapply
#' @importFrom dplyr arrange
#'
#' @examples
find.trajectory.aligned.genes <- function(expression, time, degrees.of.freedom=8, p.adjust.method="BH", q.value.cutoff=1e-10, mc.cores=1) {
  lapply.fun <- 
    if (mc.cores==1) {
      function(...) pbapply::pblapply(...)
    } else {
      function(...) parallel::mclapply(..., mc.cores=mc.cores)
    }
  runs <- lapply.fun(seq_len(ncol(expression)), function(i) {
    fit <- glm(expression[,i]~splines::ns(time, degrees.of.freedom), family=gaussian(), epsilon=1e-5)
    fit1 <- glm(expression[,i]~1, family=gaussian(), epsilon=1e-5)
    test <- anova(fit1, fit, test = "F")
    pval <- test[6][2, 1]
    
    smooth <- predict(fit, splines::ns(time, degrees.of.freedom))
    list(pval=pval, smooth=smooth)
  })
  pvals <- sapply(runs, function(x) x$pval)
  pvals[!is.finite(pvals)] <- 1
  qvals <- p.adjust(pvals, method=p.adjust.method, n=ncol(expression))
  smooth.expression <- sapply(runs, function(x) x$smooth)
  dimnames(smooth.expression) <- dimnames(expression)
  
  breaks <- c(1, .05, .001, 1e-5, 1e-10, 1e-20, 1e-40, -Inf)
  qval.cat <- cut(qvals, breaks=breaks)
  levels(qval.cat) <- c("p <= 1e-40", "1e-40 < p <= 1e-20", "1e-20 < p <= 1e-10", "1e-10 < p <= 1e-5", "1e-5 < p <= 0.001", "0.001 < p <= 0.05", "p > 0.05")
  
  is.tag <- qvals < q.value.cutoff
  
  p.values <- data.frame(gene=colnames(expression), p.value=pvals, q.value=qvals, is.tag=is.tag, category=qval.cat, stringsAsFactors = F)
  p.values <- dplyr::arrange(p.values, q.value)
  
  genes <- p.values$gene[p.values$is.tag]
  list(genes=genes, p.values=p.values, smooth.expression=smooth.expression)
}

#' Title
#'
#' @param counts 
#' @param tag.genes 
#' @param time 
#'
#' @return
#' @export
#' 
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom WGCNA mergeCloseModules
#' @importFrom dplyr bind_rows arrange
#'
#' @examples
find.modules <- function(smooth.expression, tag.genes) {
  smooth.expr <- smooth.expression[,tag.genes,drop=F]
  dissim <- correlation.distance(t(smooth.expr))
  hcl <- hclust(as.dist(dissim), method="average")
  
  labels <- dynamicTreeCut::cutreeDynamic(hcl, distM=dissim, cutHeight = 0.8, deepSplit=1, pamRespectsDendro = F, method="hybrid")
  labels2 <- WGCNA::mergeCloseModules(smooth.expr, labels, cutHeight = 0.3)$colors
  labels2 <- match(labels2, unique(labels2))
  # labels2 <- dynamicTreeCut::cutreeDynamic(hcl, distM=dissim, cutHeight = 0.8, deepSplit=1, pamRespectsDendro = F, method="hybrid")
  
  modules <- dplyr::bind_rows(lapply(unique(labels2), function(l) {
    ix <- which(labels2==l)
    dimred <- reduce.dimensionality(dissim[ix, ix, drop=F], ndim=2)
    value <- dimred[,1]
    data.frame(gene=tag.genes[ix], index=ix, module=l, value=value, stringsAsFactors = F)
  }))
  modules <- dplyr::arrange(modules, module, value)
  
  modules
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
#' @importFrom dplyr bind_rows
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
#' @param tags 
#' @param progression.str 
#' @param time 
#' @param modules 
#'
#' @return
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#'
#' @examples
plot.modules.heatmap <- function(smooth.expression, tag.genes, progression.str, time, modules) {
  smooth.expr <- smooth.expression[,tag.genes,drop=F]
  smooth.expr.ord <- smooth.expr[order(time), modules$index]
  smooth.expr.ord.scaled <- t(scale(smooth.expr.ord))
  gaps <- which(modules$module[-1] != modules$module[-length(modules$module)])
  col.ann <- data.frame(row.names = rownames(smooth.expression), Progression=progression.str, Time=time)
  row.ann <- data.frame(row.names = modules$gene, Module=paste0("Module ", modules$module))
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  num.modules <- length(unique(modules$module))
  Module.colours <- setNames(RColorBrewer::brewer.pal(max(num.modules, 3), "Set2"), paste0("Module ", seq_len(num.modules)))
  Prog.colours <- setNames(gg_color_hue(length(levels(progression.str))), levels(progression.str))
  ann.col <- list(
    Module=Module.colours,
    Time=Prog.colours,
    Progression=Prog.colours
  )
  
  pheatmap::pheatmap(
    smooth.expr.ord.scaled, 
    cluster_rows=F, 
    cluster_cols=F, 
    gaps_row=gaps, 
    annotation_row = row.ann, 
    annotation_col = col.ann, 
    annotation_colors = ann.col, 
    breaks=c(min(smooth.expr.ord.scaled), seq(-2, 2, length.out=99), max(smooth.expr.ord.scaled)),
    labels_row=rep("", nrow(smooth.expr.ord.scaled)),
    labels_col=rep("", ncol(smooth.expr.ord.scaled))
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

