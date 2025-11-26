# Infer a trajectory using SCORPIUS

Pass this object to
[`dynwrap::infer_trajectory()`](https://rdrr.io/pkg/dynwrap/man/infer_trajectories.html).

## Usage

``` r
ti_scorpius(
  distance_method = "spearman",
  ndim = 3L,
  k = 4L,
  thresh = 0.001,
  maxit = 10L,
  stretch = 0,
  smoother = "smooth_spline"
)
```

## Arguments

- distance_method:

  A character string indicating which correlationcoefficient (or
  covariance) is to be computed. One of "pearson", "spearman" (default),
  or "cosine". Domain: {spearman, pearson, cosine}. Default: spearman.
  Format: character.

- ndim:

  The number of dimensions in the new space. Domain: U(2, 20).
  Default: 3. Format: integer.

- k:

  The number of clusters to cluster the data into to construct the
  initial trajectory. Domain: U(1, 20). Default: 4. Format: integer.

- thresh:

  `principal_curve` parameter; convergence threshhold on shortest
  distances to the curve. Domain: e^U(-11.51, 11.51). Default: 0.001.
  Format: numeric.

- maxit:

  `principal_curve` parameter; maximum number of iterations. Domain:
  U(0, 50). Default: 10. Format: integer.

- stretch:

  `principal_curve` parameter; a factor by which the curve can be
  extrapolated when points are projected. Domain: U(0, 5). Default: 0.
  Format: numeric.

- smoother:

  `principal_curve` parameter; choice of smoother. Domain:
  {smooth_spline, lowess, periodic_lowess}. Default: smooth_spline.
  Format: character.

## Value

A dynwrap TI method.
