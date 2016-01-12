## ----message=FALSE-------------------------------------------------------
library(SCORPIUS)
dataset <- generate.dataset(type="poly", num.genes=500, num.samples=384, num.groups=4)

## ------------------------------------------------------------------------
dataset$expression[1:6, 1:6]

## ------------------------------------------------------------------------
head(dataset$sample.info)

## ------------------------------------------------------------------------
dist <- correlation.distance(dataset$expression)

## ------------------------------------------------------------------------
dim(dist)
plot(density(dist))

## ------------------------------------------------------------------------
space <- reduce.dimensionality(dist, ndim=3)

## ------------------------------------------------------------------------
draw.trajectory.plot(space)

## ------------------------------------------------------------------------
group.name <- dataset$sample.info$group.name
draw.trajectory.plot(space, progression.group = group.name)

## ------------------------------------------------------------------------
traj <- infer.trajectory(space)

## ------------------------------------------------------------------------
draw.trajectory.plot(space, progression.group = group.name, path = traj$final.path)

## ----find tafs unparallel, eval=F----------------------------------------
#  tafs <- find.trajectory.aligned.features(dataset$expression, traj$time)
#  expr.tafs <- tafs$smooth.x[,tafs$tafs]
#  modules <- extract.modules(expr.tafs)
#  draw.trajectory.heatmap(expr.tafs, traj$time, group.name, modules)

## ----find tafs, echo=F, message=F----------------------------------------
tafs <- find.trajectory.aligned.features(dataset$expression, traj$time, parallel=8)
expr.tafs <- tafs$smooth.x[,tafs$tafs]
modules <- extract.modules(expr.tafs)
draw.trajectory.heatmap(expr.tafs, traj$time, group.name, modules)

