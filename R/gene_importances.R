#' @title Calculate the importance of a feature
#'
#' @description Calculates the feature importance of each column in \code{x} in trying to predict the time ordering.
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory as returned by \code{\link{infer.trajectory}}.
#' @param num.permutations The number of permutations to test against for calculating the p-values (default: 0).
#' @param ntree The number of trees to grow (default: 10000).
#' @param mtry The number of variables randomly samples at each split (default: 1\% of features).
#' @param num.threads Number of threads. Default is number of CPU cores available.
#' @param ... Extra parameters passed to ranger.
#'
#' @return a data frame containing the importance of each feature for the given time line
#'
#' @importFrom ranger ranger
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
gene.importances <- function(x, time, num.permutations = 0, ntree = 10000, mtry = ncol(x) * .01, num.threads = NULL, ...) {
  data <- data.frame(x, XXXtimeXXX = time, check.names = F, stringsAsFactors = F)
  importance <- ranger::ranger(data = data, dependent.variable.name = "XXXtimeXXX", num.trees = ntree, mtry = mtry, importance = "impurity", ...)$variable.importance
  if (num.permutations > 0) {
    perms <- unlist(pbapply::pblapply(seq_len(num.permutations), function(i) {
      data$time <- sample(data$time)
      ranger::ranger(
        data = data,
        dependent.variable.name = "XXXtimeXXX",
        num.trees = ntree,
        mtry = mtry,
        importance = "impurity",
        num.threads = num.threads,
        ...
      )$variable.importance
    }))
    pvalue <- sapply(importance, function(x) mean(x < perms))
  } else {
    pvalue <- rep(NA, length(importance))
  }
  data_frame(gene = colnames(x), importance, pvalue) %>% arrange(desc(importance))
}
