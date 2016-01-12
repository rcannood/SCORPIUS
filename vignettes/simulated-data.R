## ------------------------------------------------------------------------
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
group <- dataset$sample.info$group.name
draw.trajectory.plot(space, progression.group = group)

## ------------------------------------------------------------------------
traj <- infer.trajectory(space, k=4)

## ------------------------------------------------------------------------
draw.trajectory.plot(space, progression.group = group, path = traj$final.path)

