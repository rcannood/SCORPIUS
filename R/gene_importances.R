#' @title Calculate the importance of a feature
#'
#' @description Calculates the feature importance of each column in \code{x} in trying to predict the time ordering.
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory as returned by \code{\link{infer.trajectory}}.
#' @param num_permutations The number of permutations to test against for calculating the p-values (default: 0).
#' @param ntree The number of trees to grow (default: 10000).
#' @param ntree_perm The number of trees to grow for each of the permutations (default: ntree / 10).
#' @param mtry The number of variables randomly samples at each split (default: 1\% of features).
#' @param num_threads Number of threads. Default is 1.
#' @param ... Extra parameters passed to \code{\link[ranger]{ranger}}.
#'
#' @return a data frame containing the importance of each feature for the given time line
#'
#' @importFrom ranger ranger
#' @importFrom pbapply pblapply
#' @export
#'
#' @examples
#' dataset <- generate_dataset(type="s", num_genes=500, num_samples=300, num_groups=4)
#' expression <- dataset$expression
#' group_name <- dataset$sample_info$group_name
#' dist <- correlation_distance(expression)
#' space <- reduce_dimensionality(dist, ndim=2)
#' traj <- infer_trajectory(space)
#' # set ntree to at least 1000!
#' gene_importances(expression, traj$time, num_permutations = 0, ntree = 1000)
gene_importances <- function(
  x,
  time,
  num_permutations = 0,
  ntree = 10000,
  ntree_perm = ntree / 10,
  mtry = ncol(x) * .01,
  num_threads = 1,
  ...
) {
  data <- data.frame(x, XXXtimeXXX = time, check.names = F, stringsAsFactors = F)
  importance <- ranger::ranger(
    data = data,
    dependent.variable.name = "XXXtimeXXX",
    num.trees = ntree,
    mtry = mtry,
    importance = "impurity",
    num.threads = num_threads,
    ...
  )$variable.importance
  if (num_permutations > 0) {
    perms <- unlist(pbapply::pblapply(seq_len(num_permutations), function(i) {
      data$XXXtimeXXX <- sample(data$XXXtimeXXX)
      ranger::ranger(
        data = data,
        dependent.variable.name = "XXXtimeXXX",
        num.trees = ntree_perm,
        mtry = mtry,
        importance = "impurity",
        num.threads = num_threads,
        ...
      )$variable.importance
    }))
    pvalue <- sapply(importance, function(x) mean(x < perms))
  } else {
    pvalue <- rep(NA, length(importance))
  }
  data_frame(gene = colnames(x), importance, pvalue) %>% arrange(desc(importance))
}
