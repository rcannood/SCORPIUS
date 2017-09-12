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
  .Deprecated(new = "generate_dataset", package = "SCORPIUS")
  generate_dataset(...)
}
reduce.dimensionality <- function(...) {
  .Deprecated(new = "reduce_dimensionality", package = "SCORPIUS")
  reduce_dimensionality(...)
}
reduce.dimensionality.landmarked <- function(...) {
  .Deprecated(new = "reduce_dimensionality_landmarked", package = "SCORPIUS")
  reduce_dimensionality_landmarked(...)
}
euclidean.distance <- function(...) {
  .Deprecated(new = "euclidean_distance", package = "SCORPIUS")
  euclidean_distance(...)
}
correlation.distance <- function(...) {
  .Deprecated(new = "correlation_distance", package = "SCORPIUS")
  correlation_distance(...)
}
knn.distances <- function(...) {
  .Deprecated(new = "knn_distances", package = "SCORPIUS")
  knn_distances(...)
}
evaluate.dim.red <- function(...) {
  .Deprecated(new = "evaluate_dim_red", package = "SCORPIUS")
  evaluate_dim_red(...)
}
evaluate.trajectory <- function(...) {
  .Deprecated(new = "evaluate_trajectory", package = "SCORPIUS")
  evaluate_trajectory(...)
}
extract.modules <- function(...) {
  .Deprecated(new = "extract_modules", package = "SCORPIUS")
  extract_modules(...)
}
gene.importances <- function(...) {
  .Deprecated(new = "gene_importances", package = "SCORPIUS")
  gene_importances(...)
}
outlier.filter <- function(...) {
  .Deprecated(new = "outlier_filter", package = "SCORPIUS")
  outlier_filter(...)
}
draw.trajectory.plot <- function(...) {
  .Deprecated(new = "draw_trajectory_plot", package = "SCORPIUS")
  draw_trajectory_plot(...)
}
draw.trajectory.heatmap <- function(...) {
  .Deprecated(new = "draw_trajectory_heatmap", package = "SCORPIUS")
  draw_trajectory_heatmap(...)
}
rescale.and.center <- function(...) {
  .Deprecated(new = "rescale_and_center", package = "SCORPIUS")
  rescale_and_center(...)
}
apply.scale <- function(...) {
  .Deprecated(new = "apply_scale", package = "SCORPIUS")
  apply_scale(...)
}
quant.scale <- function(...) {
  .Deprecated(new = "quant_scale", package = "SCORPIUS")
  quant_scale(...)
}
apply.quant.scale <- function(...) {
  .Deprecated(new = "apply_quant_scale", package = "SCORPIUS")
  apply_quant_scale(...)
}
infer.initial.trajectory <- function(...) {
  .Deprecated(new = "infer_initial_trajectory", package = "SCORPIUS")
  infer_initial_trajectory(...)
}
infer.trajectory <- function(...) {
  .Deprecated(new = "infer_trajectory", package = "SCORPIUS")
  infer_trajectory(...)
}
reverse.trajectory <- function(...) {
  .Deprecated(new = "reverse_trajectory", package = "SCORPIUS")
  reverse_trajectory(...)
}


