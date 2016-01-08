#' @title Evaluate the inferred timeline
#' 
#' @description \code{evaluate.trajectory} calculates the \emph{consistency} of 
#' the predicted time points of samples versus the known progression stages.
#' 
#' @usage
#' evaluate.trajectory(time, progression)
#' 
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory.
#' @param progression A factor or a numeric vector which represents the progression stages of each sample.
#'
#' @return The consistency value for the predicted timeline.
#' 
#' @export
#'
#' @examples
#' ## Generate a dataset
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space, k=4)
#' 
#' ## Evaluate the trajectory timeline
#' evaluate.trajectory(traj$time, dataset$sample.info$group.name)
evaluate.trajectory <- function(time, progression) {
  # input checks
  if (!is.vector(time) || !is.numeric(time)) 
    stop(sQuote("time"), " must be a numeric vector")
  if (!is.factor(progression) && (!is.vector(progression) || !is.numeric(progression)))
    stop(sQuote("progression"), " must be a numeric vector or a factor")
  if (length(time) != length(progression)) 
    stop(sQuote("time"), " and ", sQuote("progression"), " must have equal lengths.")
  
  ## Calculate the smallest distance between any two time values other than 0
  stime <- sort(time)
  diff <- stime[-1] - stime[-length(stime)]
  min.diff <- min(diff[diff != 0])
  
  ## Add small values to the time points. If there are time points with same values, samples will now be ordered randomly.
  noises <- runif(length(time), 0, 0.01) * min.diff
  noised.time <- time + noises
  
  ## Rank the time points
  rank <- rank(noised.time)
  
  ## If progression is a factor, convert it to an integer
  if (is.factor(progression)) progression <- as.integer(progression)
  
  ## Calculate whether or not pairs of samples are consistent in terms of its progression and rank
  comp <- expand.grid(i=seq_along(progression), j=seq_along(progression)) 
  comp$pi <- progression[comp$i]
  comp$pj <- progression[comp$j]
  comp$ri <- rank[comp$i]
  comp$rj <- rank[comp$j]
  comp <- comp[comp$pi != comp$pj,,drop=F]
  comp$consistent <- with(comp, (pi < pj) == (ri < rj))
  
  ## Calculate the mean consistency
  con <- mean(comp$consistent)
  
  ## Take into account undirectionality of the timeline
  con <- max(con, 1-con)
  
  ## Rescale and return
  (con-.5)*2
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