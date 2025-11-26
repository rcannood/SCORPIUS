#' Visualise SCORPIUS
#'
#' \code{draw_trajectory_plot} is used to plot samples after performing dimensionality reduction.
#' Additional arguments can be provided to colour the samples, plot the trajectory inferred by SCORPIUS,
#' and draw a contour around the samples.
#'
#' @param space A numeric matrix or a data frame containing the coordinates of samples.
#' @param progression_group \code{NULL} or a vector (or factor) containing the groupings of the samples (default \code{NULL}).
#' @param path A numeric matrix or a data frame containing the coordinates of the inferred path.
#' @param contour \code{TRUE} if contours are to be drawn around the samples.
#' @param progression_group_palette A named vector palette for the progression group.
#' @param point_size The size of the points.
#' @param point_alpha The alpha of the points.
#' @param path_size The size of the path (if any).
#' @param path_alpha The alpha of the path (if any).
#' @param contour_alpha The alpha of the contour (if any).
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
#' dataset <- generate_dataset(num_genes = 500, num_samples = 300, num_groups = 4)
#' space <- reduce_dimensionality(dataset$expression, ndim = 2)
#' groups <- dataset$sample_info$group_name
#'
#' ## Simply plot the samples
#' draw_trajectory_plot(space)
#'
#' ## Colour each sample according to its group
#' draw_trajectory_plot(space, progression_group = groups)
#'
#' ## Add contours to the plot
#' draw_trajectory_plot(space, progression_group = groups, contour = TRUE)
#'
#' ## Plot contours without colours
#' draw_trajectory_plot(space, contour = TRUE)
#'
#' ## Infer a trajectory and plot it
#' traj <- infer_trajectory(space)
#' draw_trajectory_plot(space, progression_group = groups, path = traj$path)
#' draw_trajectory_plot(space, progression_group = groups, path = traj$path, contour = TRUE)
#'
#' ## Visualise gene expression
#' draw_trajectory_plot(space, progression_group = dataset$expression[,1])
draw_trajectory_plot <- function(
  space,
  progression_group = NULL,
  path = NULL,
  contour = FALSE,
  progression_group_palette = NULL,
  point_size = 2,
  point_alpha = 1,
  path_size = .5,
  path_alpha = 1,
  contour_alpha = .2
) {
  # input checks
  check_numeric_matrix(space, "space", finite = TRUE)
  check_numeric_matrix(path, "path", finite = TRUE, is_nullable = TRUE)
  check_logical_vector(contour, "contour", length = 1)
  check_numeric_vector(progression_group, "progression_group", is_nullable = TRUE, finite = TRUE, length = nrow(space), factor = TRUE)

  # retrieve data about the range of the plot
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2

  # construct data frame
  space_df <- data.frame(space[,1:2], check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(space_df) <- c("Comp1", "Comp2")

  # if the grouping colours are specified, add these to the data frame
  if (!is.null(progression_group))
    space_df$progression_group <- progression_group

  lim <-
    if (contour) {
      c(min - .1*diff, max + .1*diff)
    } else {
      c(min, max)
    }

  # construct base ggplot
  g <- ggplot() +
    theme_classic() +
    labs(x = "Component 1", y = "Component 2", colour = "Group", fill = "Group") +
    xlim(min - diff, max + diff) +
    ylim(min - diff, max + diff) +
    coord_equal(xlim = lim, ylim = lim)

  # if a contour is desirable, add the contour layer
  if (contour) {
    if (!is.null(progression_group) && is.numeric(progression_group)) {
      stop("If contour is TRUE, the progression group must be a factor or a character.")
    }
    aes_contour <- aes_string("Comp1", "Comp2", z="density")
    if (!is.null(progression_group)) aes_contour$fill <- quote(progression_group)

    groupings <-
      if (is.null(progression_group)) {
        list(group = seq_len(nrow(space_df)))
      } else {
        unique_groups <- unique(progression_group)
        gr <- lapply(unique_groups, function(col) which(col == progression_group))
        names(gr) <- unique_groups
        gr
      }

    density_df <- map_df(names(groupings), function(group_name) {
      group_ix <- groupings[[group_name]]

      kde_out <- MASS::kde2d(
        space_df[group_ix, 1],
        space_df[group_ix, 2],
        lims = c(min - diff, max + diff, min - diff, max + diff)
      )

      rownames(kde_out$z) <- names(kde_out$x) <- paste0("row", seq_along(kde_out$x))
      colnames(kde_out$z) <- names(kde_out$y) <- paste0("col", seq_along(kde_out$y))
      names(dimnames(kde_out$z)) <- c("x", "y")

      kde_out$z %>%
        reshape2::melt(kde_out$z, varnames = c("x", "y"), value.name = "density") %>%
        as_tibble() %>%
        transmute(
          progression_group = group_name,
          Comp1 = kde_out$x[.data$x],
          Comp2 = kde_out$y[.data$y],
          density = .data$density
        )
    })

    if (!is.null(progression_group) && is.factor(progression_group))
      density_df$progression_group <- factor(density_df$progression_group, levels = levels(progression_group))

    g <- g + stat_contour(geom = "polygon", aes_contour, density_df, breaks = 1, alpha = contour_alpha)
  }

  # add the point layer
  aes_point <- aes_string("Comp1", "Comp2")
  if (!is.null(progression_group))
    aes_point$colour <- quote(progression_group)
  g <- g + geom_point(aes_point, space_df, size = point_size, alpha = point_alpha)

  # if a path is desirable, add the path layer
  if (!is.null(path))
    g <- g + geom_path(aes_string("Comp1", "Comp2"), data.frame(path), size = path_size, alpha = path_alpha)

  palette <-
    if (!is.null(progression_group_palette)) {
      progression_group_palette
    } else if (is.character(progression_group) || is.factor(progression_group)) {
      .default_discrete_palette(progression_group)
    } else if (is.numeric(progression_group)) {
      .default_continuous_palette()
    }

  if (is.character(progression_group) || is.factor(progression_group)) {
    if (!is.null(progression_group_palette)) {
      if (
        is.null(names(progression_group_palette)) ||
        !setequal(names(progression_group_palette), progression_group)
      ) {
        stop(
          "progression_group_palette must be a named vector of colours\n",
          "where the names correspond to the unique values in progression_group"
        )
      }
    }

    g <- g + scale_color_manual(values = palette)
    if (contour) {
      g <- g + scale_fill_manual(values = palette)
    }
  } else if (is.numeric(progression_group)) {
    g <- g + scale_color_gradientn(colours = palette)
  }

  # return the plot
  g
}




#' Draw time-series heatmap
#'
#' \code{draw_trajectory_heatmap} draws a heatmap in which the samples
#' are ranked according their position in an inferred trajectory. In addition, the progression groups and
#' feature modules can be passed along to further enhance the visualisation.
#'
#' @param x A numeric matrix or a data frame with one row per sample and one column per feature.
#' @param time A numeric vector containing the inferred time points of each sample along a trajectory.
#' @param progression_group \code{NULL} or a vector (or factor) containing the groupings of the samples (default \code{NULL}).
#' @param modules \code{NULL} or a data frame as returned by \code{\link{extract_modules}}.
#' @param show_labels_row \code{TRUE} if the labels of the rows are to be plotted (default \code{FALSE}).
#' @param show_labels_col \code{TRUE} if the labels of the cols are to be plotted (default \code{FALSE}).
#' @param scale_features \code{TRUE} if the values of each feature is to be scaled (default \code{TRUE}).
#' @param progression_group_palette A named vector palette for the progression group.
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
#' \donttest{
#' ## Generate a dataset
#' dataset <- generate_dataset(num_genes=500, num_samples=300, num_groups=4)
#' expression <- dataset$expression
#' space <- reduce_dimensionality(expression, ndim=2)
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
#' draw_trajectory_heatmap(expr_sel, time, progression_group=groups)
#'
#' ## Use a different palette
#' draw_trajectory_heatmap(
#'   expr_sel, time, progression_group=groups,
#'   progression_group_palette = setNames(RColorBrewer::brewer.pal(4, "Set2"), paste0("Group ", 1:4))
#' )
#'
#' ## Group the genes into modules and visualise the modules in a heatmap
#' modules <- extract_modules(scale_quantile(expr_sel))
#' draw_trajectory_heatmap(expr_sel, time, progression_group=groups, modules=modules)
#' }
draw_trajectory_heatmap <- function(
  x,
  time,
  progression_group = NULL,
  modules = NULL,
  show_labels_row = FALSE,
  show_labels_col = FALSE,
  scale_features = TRUE,
  progression_group_palette = NULL,
  ...
) {
  # remove any irrelevant parameters from time
  attributes(time) <- attributes(time)[intersect(names(attributes(time)), "names")]

  # input checks
  check_numeric_matrix(x, "x")
  check_numeric_vector(time, "time", length = nrow(x))
  check_numeric_vector(progression_group, "progression_group", is_nullable = TRUE, factor = TRUE, length = nrow(x))

  if (is.null(rownames(x))) {
    rownames(x) <- paste("Row ", seq_len(nrow(x)))
  }

  # process ... params
  params <- list(...)

  # create column annotation
  col_ann <- data.frame(row.names = rownames(x), Time = time)

  # process expression data
  x_part <- x[order(time),,drop=FALSE]
  if (scale_features) {
    x_part <- scale_quantile(x_part)
  }
  x_part <- t(x_part)

  # get palettes
  ann_col <-
    if (!is.null(params$annotation_colors)) {
      params$annotation_colors
    } else {
      list()
    }

  # add palette for time
  ann_col$Time <- RColorBrewer::brewer.pal(5, "RdGy")

  # add palette for progression
  if (!is.null(progression_group)) {
    if (is.numeric(progression_group)) {
      ann_col$Progression <-
        if (!is.null(progression_group_palette)) {
          progression_group_palette
        } else {
          .default_continuous_palette()
        }
    } else {
      if (!is.factor(progression_group)) progression_group <- factor(progression_group)

      ann_col$Progression <-
        if (!is.null(progression_group_palette)) {
          progression_group_palette
        } else {
          .default_discrete_palette(progression_group)
        }
    }

    col_ann$Progression <- progression_group
  }

  # whether or not to show the labels
  labels_row <- if (!show_labels_row) rep("", nrow(x_part)) else NULL
  labels_col <- if (!show_labels_col) rep("", ncol(x_part)) else NULL

  # whether or not to cluster by modules
  if (!is.null(modules)) {
    x_part <- x_part[modules$feature,]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- FALSE
  } else {
    gaps_row <- NULL
    cluster_rows <- TRUE
  }

  # pass parameters to pheatmap
  params$mat <- x_part
  params$cluster_cols <- FALSE
  params$cluster_rows <- cluster_rows
  params$annotation_col <- col_ann
  params$annotation_colors <- ann_col
  params$gaps_row <- gaps_row
  params$labels_row <- labels_row
  params$labels_col <- labels_col

  do.call(pheatmap::pheatmap, params)
}


.gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  grDevices::hcl(h=hues, l=65, c=100)[1:n]
}

.default_discrete_palette <- function(progression_group) {
  num_progressions <- length(levels(progression_group))
  progression_cols <-
    if (num_progressions <= 9) {
      RColorBrewer::brewer.pal(num_progressions, "Set1")
    } else {
      .gg_color_hue(num_progressions)
    }
  stats::setNames(progression_cols, levels(progression_group))
}


.default_continuous_palette <- function() {
  rev(RColorBrewer::brewer.pal(9, "RdYlBu"))
}



#' Progression Group to Factor
#'
#' This function converts the Progression Group Column to a Factor (if not already), as expected by the draw_trajectory_plot()
#'
#' @param column_name A progression group column from metadata.
#'
#' @return Progression Group as a factor.

prog_col_to_factor <- function(column_name) {
  # Convert the column to a vector
  column_vector <- c(column_name)
  
  # Convert the vector to a factor
  column_factor <- as.factor(column_vector)
  
  return(column_factor)
}
