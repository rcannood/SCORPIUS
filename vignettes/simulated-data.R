## ---- echo=F-------------------------------------------------------------
set.seed(4)

## ----message=FALSE-------------------------------------------------------
library(SCORPIUS)
dataset <- generate.dataset(type="poly", num.genes=500, num.samples=384, num.groups=4)

## ------------------------------------------------------------------------
dataset$expression[1:6, 1:6]

## ------------------------------------------------------------------------
head(dataset$sample.info)

## ------------------------------------------------------------------------
expression <- dataset$expression
group.name <- dataset$sample.info$group.name
dist <- correlation.distance(expression)

## ------------------------------------------------------------------------
dim(dist)
plot(density(dist))

## ------------------------------------------------------------------------
space <- reduce.dimensionality(dist)

## ------------------------------------------------------------------------
draw_trajectory_plot(space)

## ------------------------------------------------------------------------
draw_trajectory_plot(space, progression_group = group.name)

## ------------------------------------------------------------------------
traj <- infer.trajectory(space)

## ------------------------------------------------------------------------
draw.trajectory.plot(space, progression.group = group.name, path = traj$path)

## ----find tafs-----------------------------------------------------------
gimp <- gene.importances(expression, traj$time, num.permutations = 0)
gene.sel <- gimp[1:50,]
expr.sel <- expression[,gene.sel$gene]

## ----visualise tafs------------------------------------------------------
draw.trajectory.heatmap(expr.sel, traj$time, group.name)

## ----moduled tafs--------------------------------------------------------
modules <- extract.modules(quant.scale(expr.sel))
draw.trajectory.heatmap(expr.sel, traj$time, group.name, modules)

