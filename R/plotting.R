#' @title Visualise SCORPIUS 
#' 
#' @description \code{plot.scorpius} is used to plot samples after performing dimensionality reduction. 
#' Additional arguments can be provided to colour the samples, plot the trajectory inferred by SCORPIUS,
#' and draw a contour around the samples.
#' 
#' @usage
#' plot.scorpius(space, colour=NULL, path=NULL, contour=F)
#'
#' @param space A numeric matrix or data frame containing the coordinates of samples.
#' @param colour A vector containing the groupings of the samples (default \code{NULL}).
#' @param path A numeric matrix or data frame containing the coordinates of the inferred path.
#' @param contour \code{TRUE} if contours are to be drawn around the samples. This will only work with ggplot2 (≥ 2.0).
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
#' plot.scorpius(space)
#' 
#' ## Colour each sample according to its group
#' plot.scorpius(space, colour=groups)
#' 
#' ## Add contours to the plot
#' plot.scorpius(space, colour=groups, contour=T)
#' 
#' ## Plot contours without colours
#' plot.scorpius(space, contour=T)
#' 
#' ## Infer a trajectory and plot it
#' traj <- infer.trajectory(space, k=4)
#' plot.scorpius(space, colour=groups, path=traj$final.path)
#' plot.scorpius(space, colour=groups, path=traj$final.path, contour=T)
plot.scorpius <- function(space, colour=NULL, path=NULL, contour=F) {
  # input checks
  if ((!is.matrix(space) && !is.data.frame(space)) || !is.numeric(space))
    stop(sQuote("space"), " must be a numeric matrix or data frame")
  if ((!is.null(colour) && !is.vector(colour) && !is.factor(colour)) || (!is.null(colour) && length(colour) != nrow(space)))
    stop(sQuote("colour"), " must be a vector or a factor of length nrow(space)")
  if ((!is.null(path) && !is.matrix(path) && !is.data.frame(path)) || (!is.null(path) && !is.numeric(path)))
    stop(sQuote("path"), " must be NULL, a numeric matrix or a data frame")
  if (!is.logical(contour))
    stop(sQuote("contour"), " must be a logical")
  
  requireNamespace("ggplot2")
  
  # retrieve data about the range of the plot
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2
  
  # construct data frame
  space.df <- data.frame(space)
  
  # if the grouping colours are specified, add these to the data frame
  if (!is.null(colour))
    space.df$colour <- colour
  
  # construct base ggplot
  g <- ggplot2::ggplot() + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x="Component 1", y="Component 2", colour="Group", fill="Group") + 
    ggplot2::xlim(min-diff, max+diff) + 
    ggplot2::ylim(min-diff, max+diff) + 
    ggplot2::coord_equal(xlim=c(min-0.1*diff, max+0.1*diff), ylim=c(min-0.1*diff, max+0.1*diff))
  
  # if a contour is desirable, add the contour layer
  if (contour) {
    aes_contour <- ggplot2::aes(Comp1, Comp2, z=density) 
    if (!is.null(colour)) aes_contour$fill <- quote(colour)
    
    groupings <-
      if (is.null(colour)) {
        list(group=seq_len(nrow(space.df))) 
      } else {
        unique.groups <- unique(colour)
        gr <- lapply(unique.groups, function(col) which(col==colour))
        names(gr) <- unique.groups
        gr
      }
    
    density.df <- as.data.frame(dplyr::bind_rows(lapply(names(groupings), FUN=function(group.name) {
      group.ix <- groupings[[group.name]]
      kde.out <- MASS::kde2d(space[group.ix,1], space[group.ix,2], lims=c(min-diff, max+diff, min-diff, max+diff))
      z.melt <- reshape2::melt(kde.out$z)
      df <- data.frame(group.name, kde.out$x[z.melt$Var1], kde.out$y[z.melt$Var2], z.melt$value, stringsAsFactors = F)
      colnames(df) <- c("colour", colnames(space)[1:2], "density")
      df
    })))
    
    g <- g + ggplot2::stat_contour(geom="polygon", aes_contour, density.df, breaks=c(1), alpha=.2)
    
    # ## ggplot2 pre-2.0
    # aes_contour <- ggplot2::aes(Comp1, Comp2) 
    # if (!is.null(colour)) aes_contour$fill <- quote(colour)
    # g <- g + ggplot2::stat_density2d(aes_contour, breaks=1, geom="polygon", space.df, alpha=.1)
  }
  
  # add the point layer
  aes_point <- ggplot2::aes(Comp1, Comp2) 
  if (!is.null(colour)) 
    aes_point$colour <- quote(colour)
  g <- g + ggplot2::geom_point(aes_point, space.df)
  
  # if a path is desirable, add the path layer
  if (!is.null(path)) 
    g <- g + ggplot2::geom_path(ggplot2::aes(Comp1, Comp2), data.frame(path))
  
  # return the plot
  g
}

#' Title
#'
#' @param x 
#' @param time 
#' @param progression 
#' @param modules 
#' @param show.labels.row 
#' @param show.labels.col 
#' @param scale.features
#'
#' @return
#' @export
#' 
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
plot.trajectory.heatmap <- function(x, time, progression=NULL, modules=NULL, show.labels.row=F, show.labels.col=F, scale.features=T, narrow.breaks=T) {
  col.ann <- data.frame(row.names = rownames(x), Time=time)
  
  x.part <- x[order(time),,drop=F]
  if (scale.features) x.part <- scale(x.part)
  x.part <- t(x.part)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  ann.col <- list(
    Time=RColorBrewer::brewer.pal(5, "RdGy")
  )
  
  if (!is.null(progression)) {
    if (!is.factor(progression)) progression <- factor(progression)
    col.ann$Progression <- progression
    num.progressions <- length(levels(progression))
    progression.cols <-
      if (num.progressions <= 9) {
        RColorBrewer::brewer.pal(num.progressions, "Set1") 
      } else {
        gg_color_hue(num.progressions)
      }
    ann.col$Progression <- setNames(progression.cols, levels(progression))
  }
  
  if (!is.null(modules)) {
    x.part <- x.part[modules$index,]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  } else {
    gaps_row <- NULL
    cluster_rows <- T
  }
  
  labels_row <- if (!show.labels.row) rep("", nrow(x.part)) else NULL
  labels_col <- if (!show.labels.col) rep("", ncol(x.part)) else NULL
  
  if (scale.features && narrow.breaks) {
    break.cutoff <- quantile(abs(x.part), .9)
    breaks <- c(min(x.part), seq(-break.cutoff, break.cutoff, length.out=99), max(x.part))
  } else {
    breaks <- NA
  }
  
  pheatmap::pheatmap(
    x.part, 
    cluster_cols = F, 
    cluster_rows = cluster_rows,
    annotation_col = col.ann, 
    annotation_colors = ann.col, 
    gaps_row = gaps_row,
    breaks = breaks,
    labels_row = labels_row,
    labels_col = labels_col
  )
}