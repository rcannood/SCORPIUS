#' SCORPIUS: Trajectory inference from single-cell RNA sequencing data.
#'
#' SCORPIUS orders single cells with regard to an implicit timeline,
#' such as cellular development or progression over time.
#'
#' @docType package
#' @name SCORPIUS-package
#' @aliases SCORPIUS-package SCORPIUS
#'
#' @references
#' Cannoodt R. et al.,
#' \href{https://www.biorxiv.org/content/early/2016/10/07/079509}{SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development},
#' bioRxiv (Oct., 2016). DOI: \href{https://doi.org/10.1101/079509}{10.1101/079509}
#' (\href{https://www.biorxiv.org/content/early/2016/10/07/079509.full.pdf}{PDF}).
#'
#' @section Dimensionality Reduction functions:
#' \code{\link{reduce_dimensionality}}, \code{\link{correlation_distance}}, \code{\link{euclidean_distance}}
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
#' @section Scaling functions:
#' \code{\link{scale_uniform}}, \code{\link{scale_quantile}}, \code{\link{scale_minmax}}
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr spread gather crossing
#' @importFrom purrr %>% map map_df map_chr map_lgl map_int map_dbl keep
#' @importFrom magrittr %<>% %$%
NULL
