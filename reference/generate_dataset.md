# Generate a synthetic dataset

`generate_dataset` generates an synthetic dataset which can be used for
visualisation purposes.

## Usage

``` r
generate_dataset(
  num_samples = 400,
  num_genes = 500,
  num_groups = 4
)
```

## Arguments

- num_samples:

  The number of samples the dataset will contain.

- num_genes:

  The number of genes the dataset will contain.

- num_groups:

  The number of groups the samples will be split up in.

## Value

A list containing the expression data and the meta data of the samples.

## See also

[SCORPIUS](rcannood.github.io/SCORPIUS/reference/SCORPIUS-package.md)

## Examples

``` r
## Generate a dataset
dataset <- generate_dataset(num_genes = 200, num_samples = 400, num_groups = 4)

## Reduce dimensionality and infer trajectory with SCORPIUS
space <- reduce_dimensionality(dataset$expression, ndim = 2)
traj <- infer_trajectory(space)

## Visualise
draw_trajectory_plot(space, path=traj$path, progression_group=dataset$sample_info$group_name)
#> Ignoring unknown labels:
#> • fill : "Group"
```
