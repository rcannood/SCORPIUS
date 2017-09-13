#' Deprecated function(s) in the SCORPIUS package
#'
#' These functions are provided for compatibility with older version of
#' the SCORPIUS package.  They may eventually be completely removed.
#' @rdname SCORPIUS-deprecated
#' @name SCORPIUS-deprecated
#' @param ... Parameters to be passed to the new version of the function
#' @docType package
#' @export  generate.dataset reduce.dimensionality reduce.dimensionality.landmarked euclidean.distance correlation.distance knn.distances evaluate.dim.red evaluate.trajectory extract.modules gene.importances outlier.filter draw.trajectory.plot draw.trajectory.heatmap rescale.and.center apply.scale quant.scale apply.quant.scale infer.initial.trajectory infer.trajectory reverse.trajectory
#' @aliases generate.dataset reduce.dimensionality reduce.dimensionality.landmarked euclidean.distance correlation.distance knn.distances evaluate.dim.red evaluate.trajectory extract.modules gene.importances outlier.filter draw.trajectory.plot draw.trajectory.heatmap rescale.and.center apply.scale quant.scale apply.quant.scale infer.initial.trajectory infer.trajectory reverse.trajectory
#' @section Details:
#' \tabular{rl}{
#'   \code{generate.dataset} \tab now a synonym for \code{\link{generate_dataset}}\cr
#' }
#'
generate.dataset <- function(...) {
  .Defunct(new = "generate_dataset", package = "SCORPIUS")
}
reduce.dimensionality <- function(...) {
  .Defunct(new = "reduce_dimensionality", package = "SCORPIUS")
}
reduce.dimensionality.landmarked <- function(...) {
  .Defunct(new = "reduce_dimensionality_landmarked", package = "SCORPIUS")
}
euclidean.distance <- function(...) {
  .Defunct(new = "euclidean_distance", package = "SCORPIUS")
}
correlation.distance <- function(...) {
  .Defunct(new = "correlation_distance", package = "SCORPIUS")
}
knn.distances <- function(...) {
  .Defunct(new = "knn_distances", package = "SCORPIUS")
}
evaluate.dim.red <- function(...) {
  .Defunct(new = "evaluate_dim_red", package = "SCORPIUS")
}
evaluate.trajectory <- function(...) {
  .Defunct(new = "evaluate_trajectory", package = "SCORPIUS")
}
extract.modules <- function(...) {
  .Defunct(new = "extract_modules", package = "SCORPIUS")
}
gene.importances <- function(...) {
  .Defunct(new = "gene_importances", package = "SCORPIUS")
}
outlier.filter <- function(...) {
  .Defunct(new = "outlier_filter", package = "SCORPIUS")
}
draw.trajectory.plot <- function(...) {
  .Defunct(new = "draw_trajectory_plot", package = "SCORPIUS")
}
draw.trajectory.heatmap <- function(...) {
  .Defunct(new = "draw_trajectory_heatmap", package = "SCORPIUS")
}
rescale.and.center <- function(...) {
  .Defunct(new = "scale_uniform", package = "SCORPIUS")
}
apply.scale <- function(...) {
  .Defunct(new = "apply_uniform_scale", package = "SCORPIUS")
}
quant.scale <- function(...) {
  .Defunct(new = "scale_quantile", package = "SCORPIUS")
}
apply.quant.scale <- function(...) {
  .Defunct(new = "apply_quantile_Scale", package = "SCORPIUS")
}
infer.initial.trajectory <- function(...) {
  .Defunct(new = "infer_initial_trajectory", package = "SCORPIUS")
}
infer.trajectory <- function(...) {
  .Defunct(new = "infer_trajectory", package = "SCORPIUS")
}
reverse.trajectory <- function(...) {
  .Defunct(new = "reverse_trajectory", package = "SCORPIUS")
}


