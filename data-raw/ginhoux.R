library(tidyverse)
library(SCORPIUS)

# load raw data
ginhoux <- readRDS("data-raw/ginhoux_orig.rds")
expression <- ginhoux$expression
sample_info <- ginhoux$sample_info

# check before normalisation
space <- reduce_dimensionality(expression)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, sample_info$group_name, path = traj$path)
gimp <- gene_importances(expression, traj$time)
draw_trajectory_heatmap(expression[, gimp$gene[1:150]], traj$time, sample_info$group_name)

# normalise
norm <- dynnormaliser::normalise_filter_counts(counts = round(2^expression) - 1)
expression <- norm$expression
sample_info <- sample_info[rownames(expression), , drop = FALSE]

# check after normalisation
space <- reduce_dimensionality(expression)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, sample_info$group_name, path = traj$path)
gimp <- gene_importances(expression, traj$time)
draw_trajectory_heatmap(expression[, gimp$gene[1:150]], traj$time, sample_info$group_name)

# filter genes for CRAN
hvg_genes <- names(sort(apply(expression, 2, var), decreasing = TRUE)[1:2000])
expression <- expression[, hvg_genes, drop = FALSE]

# check after gene filtering
space <- reduce_dimensionality(expression)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, sample_info$group_name, path = traj$path)
gimp <- gene_importances(expression, traj$time)
draw_trajectory_heatmap(expression[, gimp$gene[1:150]], traj$time, sample_info$group_name)

# save data
ginhoux <- lst(
  expression,
  sample_info
)
devtools::use_data(ginhoux, overwrite = TRUE)
