## ---- echo=F-------------------------------------------------------------
set.seed(4)

## ----message=FALSE-------------------------------------------------------
library(SCORPIUS)
data(ginhoux)

## ------------------------------------------------------------------------
ginhoux$expression[1:6, 1:6]

## ------------------------------------------------------------------------
head(ginhoux$sample.info)

## ------------------------------------------------------------------------
expression <- ginhoux$expression
group.name <- ginhoux$sample.info$group.name
dist <- correlation.distance(expression)

## ------------------------------------------------------------------------
dim(dist)
plot(density(dist))

## ------------------------------------------------------------------------
space <- reduce.dimensionality(dist)

## ------------------------------------------------------------------------
draw.trajectory.plot(space)

## ------------------------------------------------------------------------
draw.trajectory.plot(space, progression.group = group.name)

## ------------------------------------------------------------------------
library(ggplot2)
draw.trajectory.plot(space[, c(1, 3)]) + labs(y="Component 3")

## ------------------------------------------------------------------------
filt <- outlier.filter(dist)
expression <- expression[filt, ]
group.name <- group.name[filt]
dist <- dist[filt, filt]
space <- reduce.dimensionality(dist)

## ------------------------------------------------------------------------
draw.trajectory.plot(space[, c(1, 2)])
draw.trajectory.plot(space[, c(1, 3)]) + labs(y = "Component 3")
draw.trajectory.plot(space[, c(2, 3)]) + labs(x = "Component 2", y = "Component 3")

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

