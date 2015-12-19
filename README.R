#' ---
#' title: "SCORPIUS"
#' output:
#'   md_document: 
#'      variant: markdown_github
#' ---

library(SCORPIUS)
library(dplyr)
library(tidyr)

#' ## Loading in the data 
#' Load in the dataset from Schlitzer et al.
data(ginhoux.dataset)

#' This matrix contains the expression data
counts <- ginhoux.dataset$counts
dim(counts)

#' This is info pertaining the cells
sample.info <- ginhoux.dataset$sample.info
head(sample.info)

celltypes <- ginhoux.dataset$celltypes %>% mutate(i=seq_along(parent)) %>% select(-parent)
celltypes

#' We construct a vector of the progression of cells as integers and as factors
progression <- sample.info %>% left_join(celltypes, by=colnames(celltypes)[-ncol(celltypes)]) %>% .$i

levels <- celltypes %>% select(-i) %>% apply(1, paste, collapse="_")
progression.str <- factor(sample.info[,colnames(celltypes)[-ncol(celltypes)],drop=F] %>% apply(1, paste, collapse="_"), levels)

names(progression.str) <- names(progression) <- rownames(counts)

#' ## SCORPIUS
#' Calculate distance
dist <- calculate.distance(counts)

#' Calculate outliers
out <- calculate.outlier(dist)
filter <- !out$is.outlier

#' Filter away outliers
counts <- counts[filter,,drop=F]
dist <- dist[filter, filter]
sample.info <- sample.info[filter,,drop=F]
progression.str <- progression.str[filter]
progression <- progression[filter]

#' Reduce dimensionality
space <- reduce.dimensionality(dist, ndim=3)
evaluate.space(space, progression)

#' Different plotting options
plot.dimensionality.reduction(space, colour = progression.str)
plot.dimensionality.reduction(space, colour = progression.str, contour=T)
plot.dimensionality.reduction(space, colour = progression.str, contour=T)
plot.dimensionality.reduction(space, contour=T)

#' Infer trajectory
trajectory <- infer.trajectory(space, k=4)
evaluate.trajectory(trajectory$time, progression)

#' Different plotting options
plot.trajectory(space, trajectory$initial.path, colour=progression.str, contour=T)
plot.trajectory(space, trajectory$final.path, colour=progression.str, contour=T)
plot.trajectory.density(trajectory$time, progression.str)

#' Reverse the trajectory if it is the wrong way around
prog.means <- tapply(trajectory$time, progression.str, mean)
if (prog.means[[last(levels)]] < prog.means[[first(levels)]]) {
  trajectory <- reverse.trajectory(trajectory)
}
time <- trajectory$time
plot.trajectory.density(trajectory$time, progression.str)

#' Find trajectory aligned genes (TAGs)
tags <- find.trajectory.aligned.genes(counts, trajectory$time, q.value.cutoff = 1e-10, mc.cores = 8)
head(tags$p.values, 20)

#' Group genes into modules
modules <- find.modules(tags$smooth.expression, tags$genes)
plot.modules.heatmap(tags$smooth.expression, tags$genes, progression.str, time, modules)
plot.modules.heatmap(counts, tags$genes, progression.str, time, modules)
