#' @title Extract modules of features
#'
#' @description \code{extract_modules} uses adaptive branch pruning to extract modules of features,
#'  which is typically done on the smoothed expression returned by \code{\link{gene_importances}}.
#'
#' @usage
#' extract_modules(x, ...)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param time (Optional) Order the modules according to a pseudotime
#' @param ... Extra parameters passed to \code{\link[mclust]{Mclust}}
#'
#' @return A data frame containing meta-data for the features in \code{x}, namely the order in which to visualise the features in and which module they belong to.
#'
#' @seealso \code{\link{gene_importances}}
#'
#' @export
#'
#' @importFrom mclust Mclust
#' @importFrom stats as.dist
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate_dataset(type="s", num_genes=500, num_samples=1000, num_groups=4)
#' expression <- dataset$expression
#' group_name <- dataset$sample_info$group_name
#' dist <- correlation_distance(expression)
#' space <- reduce_dimensionality(dist, ndim=2)
#' traj <- infer_trajectory(space)
#' time <- traj$time
#' draw_trajectory_plot(space, path=traj$path, group_name)
#'
#' ## Select most important genes
#' gimp <- gene_importances(expression, traj$time, num_permutations = 0)
#' gene_sel <- gimp[1:50,]
#' expr_sel <- expression[,gene_sel$gene]
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract_modules(quant_scale(expr_sel))
#' draw_trajectory_heatmap(expr_sel, time, group_name, modules)
extract_modules <- function(x, time = NULL, ...) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")

  feature_names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))

  # sigh.. mclust doesn't do well with requireNamespace
  mclustBIC <- mclust::mclustBIC

  # cluster with mclust
  labels <- mclust::Mclust(t(x), ...)$classification

  # determine mean module expression
  module_means <- do.call(cbind, tapply(feature_names, labels, function(fn) rowMeans(x[,fn])))

  # reorder modules
  if (!is.null(time)) {
    module_traj <- infer_trajectory(t(module_means), k = NULL)
    if (cor(time[apply(module_means, 2, which.max)], module_traj$time) < 0) {
      module_traj <- reverse_trajectory(module_traj)
    }
    labels <- order(order(module_traj$time))[labels]
  }

  # order features within one module according to a dimensionality reduction of the correlation distance
  modules <- bind_rows(lapply(sort(unique(labels)), function(l) {
    ix <- which(labels==l)

    if (length(ix) > 1) {
      value <- infer_trajectory(t(x[,ix]))$time

      if (!is.null(time) && cor(time[apply(x[,ix], 2, which.max)], value) < 0) {
        value <- 1 - value
      }

    } else {
      value <- 1
    }

    data_frame(
      feature = feature_names[ix],
      orig_index = ix,
      module = l,
      within_module_ordering = percent_rank(value)
    ) %>%
      arrange(within_module_ordering)
  }))

  # return output
  modules
}
