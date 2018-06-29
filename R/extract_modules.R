#' @title Extract modules of features
#'
#' @description \code{extract_modules} uses adaptive branch pruning to extract modules of features,
#'  which is typically done on the smoothed expression returned by \code{\link{gene_importances}}.
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param time (Optional) Order the modules according to a pseudotime
#' @param suppress_warnings Whether or not to suppress warnings when P > 1000
#' @param verbose Whether or not Mclust will print output or not
#' @param ... Extra parameters passed to \code{\link[mclust]{Mclust}}
#'
#' @return A data frame containing meta-data for the features in \code{x}, namely the order in which to visualise the features in and which module they belong to.
#'
#' @seealso \code{\link{gene_importances}}
#'
#' @export
#'
#' @importFrom mclust Mclust
#' @importFrom stats as.dist cor
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate_dataset(type="s", num_genes=500, num_samples=300, num_groups=4)
#' expression <- dataset$expression
#' group_name <- dataset$sample_info$group_name
#' space <- reduce_dimensionality(expression, correlation_distance, ndim=2)
#' traj <- infer_trajectory(space)
#' time <- traj$time
#' draw_trajectory_plot(space, path=traj$path, group_name)
#'
#' ## Select most important genes (set ntree to at least 10000!)
#' gimp <- gene_importances(expression, traj$time, num_permutations = 0, ntree = 1000)
#' gene_sel <- gimp[1:50,]
#' expr_sel <- expression[,gene_sel$gene]
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract_modules(scale_quantile(expr_sel))
#' draw_trajectory_heatmap(expr_sel, time, group_name, modules)
extract_modules <- function(x, time = NULL, suppress_warnings = FALSE, verbose = FALSE, ...) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")

  if (!suppress_warnings && ncol(x) > 1000) {
    warning(sQuote("x"), " has more than 1000 features. ",
            "This might take a while. Use suppress_warnings = TRUE ",
            "if you do not wish to get this message.")
  }

  feature_names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))

  # sigh.. mclust doesn't do well with requireNamespace
  mclustBIC <- mclust::mclustBIC

  # cluster with mclust
  labels <- mclust::Mclust(t(x), verbose = verbose, ...)$classification

  # determine mean module expression
  module_means <- do.call(cbind, tapply(feature_names, labels, function(fn) rowMeans(x[,fn,drop=FALSE])))

  order_data <- function(z) {
    if (ncol(z) <= 2) {
      pct <- seq(0, 1, length.out = ncol(z))
    } else if (ncol(z) == 3) {
      pct <- reduce_dimensionality(t(z), correlation_distance, ndim = 1)[,1]
      pct <- (pct - min(pct)) / (max(pct) - min(pct))
    } else {
      pct <- suppressWarnings(infer_trajectory(t(z), k = NULL)$time)
    }
    if (!is.null(time) && ncol(z) > 1 && stats::cor(stats::cor(time, z)[1,], pct) < 0) {
      pct <- -pct
    }
    pct
  }

  # reorder modules
  module_time <- order_data(module_means)
  labels <- order(order(module_time))[labels]

  # order features within one module according to a dimensionality reduction of the correlation distance
  modules <- bind_rows(lapply(sort(unique(labels)), function(l) {
    ix <- which(labels==l)

    value <- order_data(x[,ix,drop=FALSE])
    within_module_ordering <- percent_rank(value)

    data_frame(
      feature = feature_names[ix],
      orig_index = ix,
      module = l,
      within_module_ordering
    ) %>%
      arrange(within_module_ordering)
  }))

  # return output
  modules
}
