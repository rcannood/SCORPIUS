#' SCORPIUS: Trajectory inference from single-cell RNA sequencing data.
#'
#' SCORPIUS orders single cells with regard to an implicit timeline,
#' such as cellular development or progression over time.
#'
#' @section Outlier functions:
#' \code{\link{outlierness}} \code{\link{outlier.filter}}
#'
#' @section Distance functions:
#' \code{\link{correlation.distance}}, \code{\link{euclidean.distance}}, \code{\link{knn}}
#'
#' @section Dimensionality Reduction functions:
#' \code{\link{reduce.dimensionality}}, \code{\link{rescale.and.center}}
#'
#' @section Trajectory Inference functions:
#' \code{\link{infer.trajectory}}, \code{\link{reverse.trajectory}}, \code{\link{find.trajectory.aligned.features}}, \code{\link{extract.modules}}
#'
#' @section Visualisation functions:
#' \code{\link{draw.trajectory.plot}}, \code{\link{draw.trajectory.heatmap}}
#'
#' @section Datasets:
#' \code{\link{generate.dataset}}, \code{\link{ginhoux}}
#'
#' @section Evaluation functions:
#' \code{\link{evaluate.trajectory}}, \code{\link{evaluate.dim.red}}
#'
#' @docType package
#' @name SCORPIUS
NULL
