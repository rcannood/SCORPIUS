#' SCORPIUS: Trajectory inference from single-cell RNA sequencing data.
#'
#' SCORPIUS orders single cells with regard to an implicit timeline,
#' such as cellular development or progression over time.
#'
#' @section Distance functions:
#' \code{\link{correlation_distance}}, \code{\link{euclidean_distance}}
#'
#' @section Dimensionality Reduction functions:
#' \code{\link{reduce_dimensionality}}, \code{\link{scale_uniform}}, \code{\link{scale_quantile}}
#'
#' @section Trajectory Inference functions:
#' \code{\link{infer_trajectory}}, \code{\link{infer_initial_trajectory}}, \code{\link{reverse_trajectory}}, \code{\link{gene_importances}}, \code{\link{extract_modules}}
#'
#' @section Visualisation functions:
#' \code{\link{draw_trajectory_plot}}, \code{\link{draw_trajectory_heatmap}}
#'
#' @section Datasets:
#' \code{\link{generate_dataset}}, \code{\link{ginhoux}}
#'
#' @section Evaluation functions:
#' \code{\link{evaluate_trajectory}}, \code{\link{evaluate_dim_red}}
#'
#' @docType package
#' @name SCORPIUS
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr spread gather
#' @importFrom purrr %>% map map_df map_chr map_lgl map_int map_dbl keep
#' @importFrom magrittr %<>% %$%
NULL
