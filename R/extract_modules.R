#' @title Extract modules of features
#'
#' @description \code{extract_modules} uses adaptive branch pruning to extract modules of features,
#'  which is typically done on the smoothed expression returned by \code{\link{gene_importances}}.
#'
#' @param x A numeric matrix or a data frame with \emph{M} rows (one per sample) and \emph{P} columns (one per feature).
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
#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats as.dist cor
#' @importFrom Matrix t rowMeans
#'
#' @examples
#' ## Generate a dataset and visualise
#' dataset <- generate_dataset(num_genes=300, num_samples=200, num_groups=4)
#' expression <- dataset$expression
#' group_name <- dataset$sample_info$group_name
#' space <- reduce_dimensionality(expression, ndim=2)
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
  check_numeric_matrix(x, "x", finite = TRUE, sparse = TRUE)
  check_numeric_vector(time, "time", finite = TRUE, length = nrow(x), is_nullable = TRUE)

  if (!suppress_warnings && ncol(x) > 1000) {
    warning(sQuote("x"), " has more than 1000 features. ",
            "This might take a while. Use suppress_warnings = TRUE ",
            "if you do not wish to get this message.")
  }

  feature_names <- if (!is.null(colnames(x))) colnames(x) else seq_len(ncol(x))

  # cluster with mclust
  labels <- mclust::Mclust(Matrix::t(x), verbose = verbose, ...)$classification

  # determine mean module expression
  module_means <- do.call(cbind, tapply(feature_names, labels, function(fn) Matrix::rowMeans(x[,fn,drop=FALSE])))

  # reorder modules
  module_time <- .extract_modules_order_data(module_means, time)
  labels <- order(order(module_time))[labels]

  # order features within one module according to a dimensionality reduction of the correlation distance
  modules <- map_df(sort(unique(labels)), function(l) {
    ix <- which(labels==l)

    value <- .extract_modules_order_data(x[,ix,drop=FALSE], time)
    within_module_ordering <- percent_rank(value)

    tibble(
      feature = feature_names[ix],
      orig_index = ix,
      module = l,
      within_module_ordering
    ) %>%
      arrange(within_module_ordering)
  })

  # return output
  modules
}

#' @importFrom Matrix t
.extract_modules_order_data <- function(z, time) {
  if (ncol(z) <= 2) {
    pct <- seq(0, 1, length.out = ncol(z))
  } else if (ncol(z) <= 5) {
    pct <- reduce_dimensionality(Matrix::t(z), "spearman", ndim = 1)[,1]
    pct <- (pct - min(pct)) / (max(pct) - min(pct))
  } else {
    tz <- Matrix::t(z)
    fi <- apply(tz, 2, function(x) length(unique(x))) > 3
    tz <- tz[, fi, drop = FALSE]
    pct <- suppressWarnings(infer_trajectory(as.matrix(tz), k = 0)$time)
  }
  if (!is.null(time) && ncol(z) > 1 && stats::cor(stats::cor(time, z)[1,], pct) < 0) {
    pct <- -pct
  }
  pct
}
