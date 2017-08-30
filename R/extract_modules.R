#' @title Extract modules of features
#'
#' @description \code{extract.modules} uses adaptive branch pruning to extract modules of features, which is typically done on the smoothed expression returned by \code{\link{gene.importances}}.
#'
#' @usage
#' extract.modules(x, ...)
#'
#' @param x A numeric matrix or data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
#' @param ... Extra parameters passed to Mclust
#'
#' @return A data frame containing meta-data for the features in \code{x}, namely the order in which to visualise the features in and which module they belong to.
#'
#' @seealso \code{\link{gene.importances}}
#'
#' @export
#'
#' @importFrom mclust Mclust
#' @importFrom stats as.dist hclust
#' @importFrom dplyr bind_rows
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' group.name <- dataset$sample.info$group.name
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' traj <- infer.trajectory(space)
#' time <- traj$time
#' draw.trajectory.plot(space, path=traj$path, group.name)
#'
#' ## Select most important genes
#' gimp <- gene.importances(expression, traj$time, num.permutations = 0)
#' gene.sel <- gimp[1:50,]
#' expr.sel <- expression[,gene.sel$gene]
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract.modules(quant.scale(expr.sel))
#' draw.trajectory.heatmap(expr.sel, time, group.name, modules)
extract.modules <- function(x, ...) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")

  feature.names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))

  # sigh.. mclust doesn't do well with requireNamespace
  mclustBIC <- mclust::mclustBIC

  # cluster with mclust
  labels <- mclust::Mclust(t(x), ...)$classification

  # hierarchically cluster the features
  dist <- correlation.distance(t(x))
  hcl <- stats::hclust(stats::as.dist(dist), method="average")

  # order features within one module according to a dimensionality reduction of the correlation distance
  modules <- dplyr::bind_rows(lapply(unique(labels), function(l) {
    ix <- which(labels==l)
    if (length(ix) > 1) {
      dimred <- reduce.dimensionality(dist[ix, ix, drop=F], ndim=1)
      data.frame(feature=feature.names[ix], index=ix, module=l, value=dimred[,1], stringsAsFactors = F, row.names = NULL)
    } else {
      data.frame(feature=feature.names[ix], index=ix, module=l, value=0, stringsAsFactors = F, row.names = NULL)
    }
  }))

  modules <- as.data.frame(modules[order(modules$module, modules$value),,drop=F])

  # return output
  modules[,c("feature", "index", "module")]
}
