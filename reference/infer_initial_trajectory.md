# Infer an initial trajectory through space

`infer_initial_trajectory` infers an initial trajectory for
[`infer_trajectory`](rcannood.github.io/SCORPIUS/reference/infer_trajectory.md)
by clustering the points and calculating the shortest path through
cluster centers. The shortest path takes into account the euclidean
distance between cluster centers, and the density between those two
points.

## Usage

``` r
infer_initial_trajectory(space, k)
```

## Arguments

- space:

  A numeric matrix or a data frame containing the coordinates of
  samples.

- k:

  The number of clusters

## Value

the initial trajectory obtained by this method
