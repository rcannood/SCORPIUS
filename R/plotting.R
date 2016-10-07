#' @title Visualise SCORPIUS
#'
#' @description \code{draw.trajectory.plot} is used to plot samples after performing dimensionality reduction.
#' Additional arguments can be provided to colour the samples, plot the trajectory inferred by SCORPIUS,
#' and draw a contour around the samples.
#'
#' @usage
#' draw.trajectory.plot(space, progression.group=NULL, path=NULL, contour=FALSE)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param progression.group \code{NULL} or a vector (or factor) containing the groupings of the samples (default \code{NULL}).
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
#' dataset <- generate.dataset(type="p", num.genes=500, num.samples=1000, num.groups=4)
#' dist <- correlation.distance(dataset$expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' groups <- dataset$sample.info$group.name
#'
#' ## Simply plot the samples
#' draw.trajectory.plot(space)
#'
#' ## Colour each sample according to its group
#' draw.trajectory.plot(space, progression.group=groups)
#'
#' ## Add contours to the plot
#' draw.trajectory.plot(space, progression.group=groups, contour=TRUE)
#'
#' ## Plot contours without colours
#' draw.trajectory.plot(space, contour=TRUE)
#'
#' ## Infer a trajectory and plot it
#' traj <- infer.trajectory(space)
#' draw.trajectory.plot(space, progression.group=groups, path=traj$path)
#' draw.trajectory.plot(space, progression.group=groups, path=traj$path, contour=TRUE)
draw.trajectory.plot <- function(space, progression.group=NULL, path=NULL, contour=FALSE) {
  # input checks
  if (!is.matrix(space) && !is.data.frame(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")
  if ((!is.null(progression.group) && !is.vector(progression.group) && !is.factor(progression.group)) || (!is.null(progression.group) && length(progression.group) != nrow(space)))
    stop(sQuote("progression.group"), " must be a vector or a factor of length nrow(space)")
  if (!is.null(path) && !is.matrix(path) && !is.data.frame(path))
    stop(sQuote("path"), " must be NULL, a numeric matrix or a data frame")
  if (!is.logical(contour))
    stop(sQuote("contour"), " must be a logical")

  requireNamespace("ggplot2")
  requireNamespace("MASS")
  requireNamespace("reshape2")
  requireNamespace("stats")

  # retrieve data about the range of the plot
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2

  # construct data frame
  space.df <- data.frame(space[,1:2], check.rows = F, check.names = F, stringsAsFactors = F)
  colnames(space.df) <- c("Comp1", "Comp2")

  # if the grouping colours are specified, add these to the data frame
  if (!is.null(progression.group))
    space.df$progression.group <- progression.group

  lim <- if (contour) c(min-.1*diff, max+.1*diff) else c(min, max)

  # construct base ggplot
  g <- ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::labs(x="Component 1", y="Component 2", colour="Group", fill="Group") +
    ggplot2::xlim(min-diff, max+diff) +
    ggplot2::ylim(min-diff, max+diff) +
    ggplot2::coord_equal(xlim=lim, ylim=lim)

  # if a contour is desirable, add the contour layer
  if (contour) {
    aes_contour <- ggplot2::aes_string("Comp1", "Comp2", z="density")
    if (!is.null(progression.group)) aes_contour$fill <- quote(progression.group)

    groupings <-
      if (is.null(progression.group)) {
        list(group=seq_len(nrow(space.df)))
      } else {
        unique.groups <- unique(progression.group)
        gr <- lapply(unique.groups, function(col) which(col==progression.group))
        names(gr) <- unique.groups
        gr
      }

    density.df <- as.data.frame(dplyr::bind_rows(lapply(names(groupings), FUN=function(group.name) {
      group.ix <- groupings[[group.name]]
      kde.out <- MASS::kde2d(space.df[group.ix,1], space.df[group.ix,2], lims=c(min-diff, max+diff, min-diff, max+diff))
      z.melt <- reshape2::melt(kde.out$z)
      df <- data.frame(group.name, kde.out$x[z.melt$Var1], kde.out$y[z.melt$Var2], z.melt$value, stringsAsFactors = F)
      colnames(df) <- c("progression.group", "Comp1", "Comp2", "density")
      df
    })))

    if (!is.null(progression.group) && is.factor(progression.group))
      density.df$progression.group <- factor(density.df$progression.group, levels = levels(progression.group))

    g <- g + ggplot2::stat_contour(geom="polygon", aes_contour, density.df, breaks=c(1), alpha=.2)

    # ## ggplot2 pre-2.0
    # aes_contour <- ggplot2::aes_string("Comp1", "Comp2")
    # if (!is.null(progression.group)) aes_contour$fill <- quote(progression.group)
    # g <- g + ggplot2::stat_density2d(aes_contour, breaks=1, geom="polygon", space.df, alpha=.1)
  }

  # add the point layer
  aes_point <- ggplot2::aes_string("Comp1", "Comp2")
  if (!is.null(progression.group))
    aes_point$colour <- quote(progression.group)
  g <- g + ggplot2::geom_point(aes_point, space.df)

  # if a path is desirable, add the path layer
  if (!is.null(path))
    g <- g + ggplot2::geom_path(ggplot2::aes_string("Comp1", "Comp2"), data.frame(path))

  # return the plot
  g
}

#' @title Draw time-series heatmap
#'
#' @description \code{draw.trajectory.heatmap} draws a heatmap in which the samples
#' are ranked according their position in an inferred trajectory. In addition, the progression groups and
#' feature modules can be passed along to further enhance the visualisation.
#'
#' @usage
#' draw.trajectory.heatmap(
#'   x,
#'   time,
#'   progression.group=NULL,
#'   modules=NULL,
#'   show.labels.row = FALSE,
#'   show.labels.col = FALSE,
#'   scale.features = TRUE,
#'   ...
#' )
#'
#' @param x A numeric matrix or data frame with one row per sample and one column per feature.
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory.
#' @param progression.group \code{NULL} or a vector (or factor) containing the groupings of the samples (default \code{NULL}).
#' @param modules \code{NULL} or a data frame as returned by \code{\link{extract.modules}}.
#' @param show.labels.row \code{TRUE} if the labels of the rows are to be plotted (default \code{FALSE}).
#' @param show.labels.col \code{TRUE} if the labels of the cols are to be plotted (default \code{FALSE}).
#' @param scale.features \code{TRUE} if the values of each feature is to be scaled (default \code{TRUE}).
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
#' ## Generate a dataset
#' dataset <- generate.dataset(type="s", num.genes=500, num.samples=1000, num.groups=4)
#' expression <- dataset$expression
#' dist <- correlation.distance(expression)
#' space <- reduce.dimensionality(dist, ndim=2)
#' groups <- dataset$sample.info$group.name
#' traj <- infer.trajectory(space)
#' time <- traj$time
#'
#' ## Select most important genes
#' gimp <- gene.importances(expression, traj$time, num.permutations = 0)
#' gene.sel <- gimp[1:50,]
#' expr.sel <- expression[,gene.sel$gene]
#'
#' ## Draw a time series heatmap
#' draw.trajectory.heatmap(expr.sel, time)
#'
#' ## Also show the progression groupings
#' draw.trajectory.heatmap(expr.sel, time, progression=groups)
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract.modules(quant.scale(expr.sel))
#' draw.trajectory.heatmap(expr.sel, time, progression.group=groups, modules=modules)
draw.trajectory.heatmap <- function(x, time, progression.group=NULL, modules=NULL, show.labels.row=FALSE, show.labels.col=FALSE, scale.features=TRUE, ...) {
  # input checks
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.vector(time) || !is.numeric(time))
    stop(sQuote("time"), " must be a numeric vector")
  if (nrow(x) != length(time))
    stop(sQuote("time"), " must have one value for each row in ", sQuote("x"))
  if ((!is.null(progression.group) && !is.vector(progression.group) && !is.factor(progression.group)) || (!is.null(progression.group) && length(progression.group) != nrow(x)))
    stop(sQuote("progression.group"), " must be a vector or a factor of length nrow(x)")

  requireNamespace("pheatmap")
  requireNamespace("RColorBrewer")
  requireNamespace("stats")
  requireNamespace("grDevices")

  if (is.null(rownames(x))) {
    rownames(x) <- paste("Row ", seq_len(nrow(x)))
  }

  col.ann <- data.frame(row.names = rownames(x), Time=time)

  x.part <- x[order(time),,drop=FALSE]
  if (scale.features) {
    x.part <- quant.scale(x.part)
  }
  x.part <- t(x.part)

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    grDevices::hcl(h=hues, l=65, c=100)[1:n]
  }

  ann.col <- list(
    Time=RColorBrewer::brewer.pal(5, "RdGy")
  )

  if (!is.null(progression.group)) {
    if (!is.factor(progression.group)) progression.group <- factor(progression.group)
    col.ann$Progression <- progression.group
    num.progressions <- length(levels(progression.group))
    progression.cols <-
      if (num.progressions <= 9) {
        RColorBrewer::brewer.pal(num.progressions, "Set1")
      } else {
        gg_color_hue(num.progressions)
      }
    ann.col$Progression <- stats::setNames(progression.cols, levels(progression.group))
  }

  labels_row <- if (!show.labels.row) rep("", nrow(x.part)) else NULL
  labels_col <- if (!show.labels.col) rep("", ncol(x.part)) else NULL

  if (!is.null(modules)) {
    x.part <- x.part[modules$index,]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  } else {
    gaps_row <- NULL
    cluster_rows <- T
  }

  pheatmap::pheatmap(
    x.part,
    cluster_cols = F,
    cluster_rows = cluster_rows,
    annotation_col = col.ann,
    annotation_colors = ann.col,
    gaps_row = gaps_row,
    labels_row = labels_row,
    labels_col = labels_col,
    ...
  )
}
