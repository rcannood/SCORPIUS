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
#'
#' ## Visualise the outlierness scores for each of the points
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
#' @seealso \code{\link{correlation.distance}}, \code{\link{euclidean.distance}}, \code{\link{outlierness}}
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
  remaining <- which(filt)[order(outliernesses, decreasing=TRUE)]
  removed[seq(num.pts.removed-1, num.pts.removed+1)] <- remaining

  # return a vector indicating which samples are /not/ outliers
  filt <- !ix %in% removed[seq(1, which.max(logliks))]

  # attach log likelihood information as an attribute
  loglik.df <- data.frame(amount.removed=seq_along(removed)-1, sample.index=removed, log.likelihood=logliks)
  attr(filt, "loglikelihoods") <- loglik.df

  # return outlier output
  filt
}
