
#' Title
#'
#' @param space 
#' @param colour 
#' @param contour
#'
#' @return
#' @export
#' 
#' @import ggplot2
#' @importFrom MASS kde2d
#'
#' @examples
plot.dimensionality.reduction <- function(space, colour=NULL, contour=F) {
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2
  
  g <- ggplot2::ggplot() + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x="Component 1", y="Component 2", colour="Group", fill="Group") + 
    ggplot2::xlim(min-diff, max+diff) + 
    ggplot2::ylim(min-diff, max+diff) + 
    ggplot2::coord_equal(xlim=c(min-0.1*diff, max+0.1*diff), ylim=c(min-0.1*diff, max+0.1*diff))
  
  if (contour) {
    if (is.null(colour)) {
      g <- g + ggplot2::stat_density2d(ggplot2::aes(Comp1, Comp2), breaks=1, geom="polygon", data.frame(space), alpha=.1)
    } else {
      g <- g + ggplot2::stat_density2d(ggplot2::aes(Comp1, Comp2, fill=colour), breaks=1, geom="polygon", data.frame(space, colour), alpha=.1)
    }
  }
  
  if (is.null(colour)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(Comp1, Comp2), data.frame(space))
  } else {
    g <- g + ggplot2::geom_point(ggplot2::aes(Comp1, Comp2, colour=colour), data.frame(space, colour))
  }
  
  g
}

#' Title
#'
#' @param space 
#' @param path 
#' @param colour 
#'
#' @return
#' @export
#' 
#' @import ggplot2
#'
#' @examples
plot.trajectory <- function(space, path, colour=NULL, contour=F) {
  plot.dimensionality.reduction(space, colour=colour, contour=contour) + ggplot2::geom_path(ggplot2::aes(Comp1, Comp2), data.frame(path))
}

#' Title
#'
#' @param time 
#' @param progression 
#' @param colour 
#'
#' @return
#' @export
#' 
#' @import ggplot2
#'
#' @examples
plot.trajectory.density <- function(time, progression, colour=NULL) {
  ggplot2::ggplot() + ggplot2::geom_density(ggplot2::aes(time, colour=colour), data.frame(time, colour=factor(progression))) + ggplot2::theme_classic()
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
    Time=c("#ca0020", "#f4a582", "#ffffff", "#bababa", "#404040")
  )
  
  if (!is.null(progression)) {
    if (!is.factor(progression)) progression <- factor(progression)
    col.ann$Progression <- progression
    ann.col$Progression <- setNames(gg_color_hue(length(levels(progression))), levels(progression))
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