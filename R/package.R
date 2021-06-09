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
#' Cannoodt R. et al., SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development,
#' bioRxiv (Oct., 2016). \doi{10.1101/079509}
#' ([PDF](https://www.biorxiv.org/content/biorxiv/early/2016/10/07/079509.full.pdf)).
#'
#' @section Dimensionality Reduction functions:
#' \code{\link{reduce_dimensionality}}
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
#' @importFrom dplyr tibble as_tibble transmute arrange percent_rank desc slice filter group_by summarise ungroup select
#' @import ggplot2
#' @importFrom tidyr spread gather crossing
#' @importFrom purrr %>% map map_df map_chr map_lgl map_int map_dbl
#'
#' @examples
#' ## Load dataset from Schlitzer et al., 2015
#' data("ginhoux")
#'
#' ## Reduce dimensionality and infer trajectory with SCORPIUS
#' space <- reduce_dimensionality(ginhoux$expression, "spearman")
#' traj <- infer_trajectory(space)
#'
#' ## Visualise
#' draw_trajectory_plot(
#'   space,
#'   path = traj$path,
#'   progression_group = ginhoux$sample_info$group_name
#' )
NULL
