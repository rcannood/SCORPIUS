# Infer linear trajectory through space

`infer_trajectory` infers a trajectory through samples in a given space
in a four-step process:

1.  Perform *k*-means clustering

2.  Calculate distance matrix between cluster centers using a custom
    distance function

3.  Find the shortest path connecting all cluster centers using the
    custom distance matrix

4.  Iteratively fit a curve to the given data using principal curves

## Usage

``` r
infer_trajectory(
  space,
  k = 4,
  thresh = 0.001,
  maxit = 10,
  stretch = 0,
  smoother = "smooth_spline",
  approx_points = 100
)
```

## Arguments

- space:

  A numeric matrix or a data frame containing the coordinates of
  samples.

- k:

  The number of clusters to cluster the data into.

- thresh:

  convergence threshold on shortest distances to the curve.

- maxit:

  maximum number of iterations.

- stretch:

  A stretch factor for the endpoints of the curve, allowing the curve to
  grow to avoid bunching at the end. Must be a numeric value between 0
  and 2.

- smoother:

  choice of smoother. The default is `"smooth_spline"`, and other
  choices are `"lowess"` and `"periodic_lowess"`. The latter allows one
  to fit closed curves. Beware, you may want to use `iter = 0` with
  [`lowess()`](https://rdrr.io/r/stats/lowess.html).

- approx_points:

  Approximate curve after smoothing to reduce computational time. If
  `FALSE`, no approximation of the curve occurs. Otherwise,
  `approx_points` must be equal to the number of points the curve gets
  approximated to; preferably about 100.

## Value

A list containing several objects:

- `path`: the trajectory obtained by principal curves.

- `time`: the time point of each sample along the inferred trajectory.

## See also

[`reduce_dimensionality`](rcannood.github.io/SCORPIUS/reference/reduce_dimensionality.md),
[`draw_trajectory_plot`](rcannood.github.io/SCORPIUS/reference/draw_trajectory_plot.md)

## Examples

``` r
## Generate an example dataset and visualise it
dataset <- generate_dataset(num_genes = 200, num_samples = 400, num_groups = 4)
space <- reduce_dimensionality(dataset$expression, ndim = 2)
draw_trajectory_plot(space, progression_group = dataset$sample_info$group_name)
#> Ignoring unknown labels:
#> • fill : "Group"


## Infer a trajectory through this space
traj <- infer_trajectory(space)

## Visualise the trajectory
draw_trajectory_plot(space, path=traj$path, progression_group = dataset$sample_info$group_name)
#> Ignoring unknown labels:
#> • fill : "Group"
```
