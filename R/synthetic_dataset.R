#' @title Generate a synthetic dataset
#'
#' @description \code{generate.dataset} generates an synthetic dataset which can be used for visualisation purposes.
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
#' @seealso \code{\link{correlation.distance}}, \code{\link{reduce.dimensionality}}, \code{\link{infer.trajectory}}, \code{\link{draw.trajectory.plot}}
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
#' draw.trajectory.plot(space, path=traj$final.path, progression.group=dataset$sample.info$group.name)
generate.dataset <- function(type=c("splines", "polynomial"), num.samples=1000, num.genes=100, num.groups=4) {
  # make names for each group, gene and sample
  group.names <- paste0("Group ", seq_len(num.groups))
  gene.names <- paste0("Gene", seq_len(num.genes))
  sample.names <- paste0("Sample", seq_len(num.samples))

  # match the type argument
  type <- match.arg(type, c("splines", "polynomial"))

  # construct the sample info
  x <- seq(-1, 1, length.out=num.samples)
  group <- cut(x, breaks=num.groups, labels = group.names)
  sample.info <- data.frame(row.names=sample.names, group.name=group)

  # apply function and determine noise sd
  switch(type, polynomial={
    y <- poly(x, 2)
    sd <- .012 * sqrt(num.genes)
  }, splines={
    y <- splines::ns(x, df=3)
    sd <- .06 * sqrt(num.genes)
  })

  # generate expression data
  expression <- sapply(seq_len(num.genes), function(g) {
    scale <- rnorm(ncol(y), mean=0, sd=1)
    noise <- rnorm(length(x), sd=sd)
    rowSums(sweep(y, 2, scale, "*")) + noise
  })
  dimnames(expression) <- list(sample.names, gene.names)

  # simulate genes that are not expressed
  weighted.random.sample <- function(data, weights, n){
    key <- runif(length(data)) ^ (1 / weights)
    data[order(key, decreasing=TRUE)][seq_len(n)]
  }

  undetectable <- which(expression < 0)
  undetectable <- weighted.random.sample(undetectable, -expression[undetectable], round(length(undetectable)*.5))

  # shift expression
  expression <- expression + .5

  # set everything below 0 or that was marked as undetectable to zero
  expression[expression < 0 | seq_along(expression) %in% undetectable] <- 0

  # rescale to a reasonable value
  expression <- expression / max(expression) * 20

  list(expression=expression, sample.info=sample.info)
}
