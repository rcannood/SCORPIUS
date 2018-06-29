#' @title Visualise SCORPIUS
#'
#' @description \code{draw_trajectory_plot} is used to plot samples after performing dimensionality reduction.
#' Additional arguments can be provided to colour the samples, plot the trajectory inferred by SCORPIUS,
#' and draw a contour around the samples.
#'
#' @usage
#' draw_trajectory_plot(space, progression_group=NULL, path=NULL, contour=FALSE)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param progression_group \code{NULL} or a vector (or factor) containing the groupings of the samples (default \code{NULL}).
#' @param path A numeric matrix or data frame containing the coordinates of the inferred path.
#' @param contour \code{TRUE} if contours are to be drawn around the samples.
#'
#' @return A ggplot2 plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom MASS kde2d
#' @importFrom reshape2 melt
#'
#' @examples
#' ## Generate a synthetic dataset
#' dataset <- generate_dataset(type="p", num_genes=500, num_samples=300, num_groups=4)
#' space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim=2)
#' groups <- dataset$sample_info$group_name
#'
#' ## Simply plot the samples
#' draw_trajectory_plot(space)
#'
#' ## Colour each sample according to its group
#' draw_trajectory_plot(space, progression_group=groups)
#'
#' ## Add contours to the plot
#' draw_trajectory_plot(space, progression_group=groups, contour=TRUE)
#'
#' ## Plot contours without colours
#' draw_trajectory_plot(space, contour=TRUE)
#'
#' ## Infer a trajectory and plot it
#' traj <- infer_trajectory(space)
#' draw_trajectory_plot(space, progression_group=groups, path=traj$path)
#' draw_trajectory_plot(space, progression_group=groups, path=traj$path, contour=TRUE)
draw_trajectory_plot <- function(space, progression_group = NULL, path = NULL, contour = FALSE) {
  # input checks
  if (!is.matrix(space) && !is.data.frame(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")
  if ((!is.null(progression_group) && !is.vector(progression_group) && !is.factor(progression_group)) || (!is.null(progression_group) && length(progression_group) != nrow(space)))
    stop(sQuote("progression_group"), " must be a vector or a factor of length nrow(space)")
  if (!is.null(path) && !is.matrix(path) && !is.data.frame(path))
    stop(sQuote("path"), " must be NULL, a numeric matrix or a data frame")
  if (!is.logical(contour))
    stop(sQuote("contour"), " must be a logical")

  # retrieve data about the range of the plot
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2

  # construct data frame
  space_df <- data.frame(space[,1:2], check.rows = F, check.names = F, stringsAsFactors = F)
  colnames(space_df) <- c("Comp1", "Comp2")

  # if the grouping colours are specified, add these to the data frame
  if (!is.null(progression_group))
    space_df$progression_group <- progression_group

  lim <- if (contour) c(min-.1*diff, max+.1*diff) else c(min, max)

  # construct base ggplot
  g <- ggplot() +
    theme_classic() +
    labs(x="Component 1", y="Component 2", colour="Group", fill="Group") +
    xlim(min-diff, max+diff) +
    ylim(min-diff, max+diff) +
    coord_equal(xlim=lim, ylim=lim)

  # if a contour is desirable, add the contour layer
  if (contour) {
    aes_contour <- aes_string("Comp1", "Comp2", z="density")
    if (!is.null(progression_group)) aes_contour$fill <- quote(progression_group)

    groupings <-
      if (is.null(progression_group)) {
        list(group=seq_len(nrow(space_df)))
      } else {
        unique_groups <- unique(progression_group)
        gr <- lapply(unique_groups, function(col) which(col==progression_group))
        names(gr) <- unique_groups
        gr
      }

    density_df <- as.data.frame(dplyr::bind_rows(lapply(names(groupings), FUN=function(group_name) {
      group_ix <- groupings[[group_name]]
      kde_out <- MASS::kde2d(space_df[group_ix,1], space_df[group_ix,2], lims=c(min-diff, max+diff, min-diff, max+diff))
      z_melt <- reshape2::melt(kde_out$z)
      df <- data.frame(group_name, kde_out$x[z_melt$Var1], kde_out$y[z_melt$Var2], z_melt$value, stringsAsFactors = F)
      colnames(df) <- c("progression_group", "Comp1", "Comp2", "density")
      df
    })))

    if (!is.null(progression_group) && is.factor(progression_group))
      density_df$progression_group <- factor(density_df$progression_group, levels = levels(progression_group))

    g <- g + stat_contour(geom="polygon", aes_contour, density_df, breaks=c(1), alpha=.2)
  }

  # add the point layer
  aes_point <- aes_string("Comp1", "Comp2")
  if (!is.null(progression_group))
    aes_point$colour <- quote(progression_group)
  g <- g + geom_point(aes_point, space_df)

  # if a path is desirable, add the path layer
  if (!is.null(path))
    g <- g + geom_path(aes_string("Comp1", "Comp2"), data.frame(path))

  # return the plot
  g
}

#' @title Draw time-series heatmap
#'
#' @description \code{draw_trajectory_heatmap} draws a heatmap in which the samples
#' are ranked according their position in an inferred trajectory. In addition, the progression groups and
#' feature modules can be passed along to further enhance the visualisation.
#'
#' @usage
#' draw_trajectory_heatmap(
#'   x,
#'   time,
#'   progression_group = NULL,
#'   modules = NULL,
#'   show_labels_row = FALSE,
#'   show_labels_col = FALSE,
#'   scale_features = TRUE,
#'   ...
#' )
#'
#' @param x A numeric matrix or data frame with one row per sample and one column per feature.
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory.
#' @param progression_group \code{NULL} or a vector (or factor) containing the groupings of the samples (default \code{NULL}).
#' @param modules \code{NULL} or a data frame as returned by \code{\link{extract_modules}}.
#' @param show_labels_row \code{TRUE} if the labels of the rows are to be plotted (default \code{FALSE}).
#' @param show_labels_col \code{TRUE} if the labels of the cols are to be plotted (default \code{FALSE}).
#' @param scale_features \code{TRUE} if the values of each feature is to be scaled (default \code{TRUE}).
#' @param ... Optional arguments to \code{\link[pheatmap]{pheatmap}}
#'
#' @return The output of the \code{\link[pheatmap]{pheatmap}} function.
#'
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices hcl
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' ## Generate a dataset
#' dataset <- generate_dataset(type="s", num_genes=500, num_samples=300, num_groups=4)
#' expression <- dataset$expression
#' space <- reduce_dimensionality(expression, correlation_distance, ndim=2)
#' groups <- dataset$sample_info$group_name
#' traj <- infer_trajectory(space)
#' time <- traj$time
#'
#' gimp <- gene_importances(expression, traj$time, num_permutations = 0, ntree = 10000)
#' gene_sel <- gimp[1:50,]
#' expr_sel <- expression[,gene_sel$gene]
#'
#' ## Draw a time series heatmap
#' draw_trajectory_heatmap(expr_sel, time)
#'
#' ## Also show the progression groupings
#' draw_trajectory_heatmap(expr_sel, time, progression=groups)
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract_modules(scale_quantile(expr_sel))
#' draw_trajectory_heatmap(expr_sel, time, progression_group=groups, modules=modules)
#' }
draw_trajectory_heatmap <- function(x, time, progression_group=NULL, modules=NULL, show_labels_row=FALSE, show_labels_col=FALSE, scale_features=TRUE, ...) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.vector(time) || !is.numeric(time))
    stop(sQuote("time"), " must be a numeric vector")
  if (nrow(x) != length(time))
    stop(sQuote("time"), " must have one value for each row in ", sQuote("x"))
  if ((!is.null(progression_group) && !is.vector(progression_group) && !is.factor(progression_group)) || (!is.null(progression_group) && length(progression_group) != nrow(x)))
    stop(sQuote("progression_group"), " must be a vector or a factor of length nrow(x)")

  if (is.null(rownames(x))) {
    rownames(x) <- paste("Row ", seq_len(nrow(x)))
  }

  col_ann <- data.frame(row.names = rownames(x), Time=time)

  x_part <- x[order(time),,drop=FALSE]
  if (scale_features) {
    x_part <- scale_quantile(x_part)
  }
  x_part <- t(x_part)

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    grDevices::hcl(h=hues, l=65, c=100)[1:n]
  }

  ann_col <- list(
    Time=RColorBrewer::brewer.pal(5, "RdGy")
  )

  if (!is.null(progression_group)) {
    if (!is.factor(progression_group)) progression_group <- factor(progression_group)
    col_ann$Progression <- progression_group
    num_progressions <- length(levels(progression_group))
    progression_cols <-
      if (num_progressions <= 9) {
        RColorBrewer::brewer.pal(num_progressions, "Set1")
      } else {
        gg_color_hue(num_progressions)
      }
    ann_col$Progression <- stats::setNames(progression_cols, levels(progression_group))
  }

  labels_row <- if (!show_labels_row) rep("", nrow(x_part)) else NULL
  labels_col <- if (!show_labels_col) rep("", ncol(x_part)) else NULL

  if (!is.null(modules)) {
    x_part <- x_part[modules$feature,]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  } else {
    gaps_row <- NULL
    cluster_rows <- T
  }

  pheatmap::pheatmap(
    x_part,
    cluster_cols = F,
    cluster_rows = cluster_rows,
    annotation_col = col_ann,
    annotation_colors = ann_col,
    gaps_row = gaps_row,
    labels_row = labels_row,
    labels_col = labels_col,
    ...
  )
}
