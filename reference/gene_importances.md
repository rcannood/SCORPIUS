# Calculate the importance of a feature

Calculates the feature importance of each column in `x` in trying to
predict the time ordering.

## Usage

``` r
gene_importances(
  x,
  time,
  num_permutations = 0,
  ntree = 10000,
  ntree_perm = ntree/10,
  mtry = ncol(x) * 0.01,
  num_threads = 1,
  ...
)
```

## Arguments

- x:

  A numeric matrix or a data frame with *M* rows (one per sample) and
  *P* columns (one per feature).

- time:

  A numeric vector containing the inferred time points of each sample
  along a trajectory as returned by
  [`infer_trajectory`](rcannood.github.io/SCORPIUS/reference/infer_trajectory.md).

- num_permutations:

  The number of permutations to test against for calculating the
  p-values (default: 0).

- ntree:

  The number of trees to grow (default: 10000).

- ntree_perm:

  The number of trees to grow for each of the permutations (default:
  ntree / 10).

- mtry:

  The number of variables randomly samples at each split (default: 1% of
  features).

- num_threads:

  Number of threads. Default is 1.

- ...:

  Extra parameters passed to
  [`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md).

## Value

a data frame containing the importance of each feature for the given
time line

## Examples

``` r
dataset <- generate_dataset(num_genes=500, num_samples=300, num_groups=4)
expression <- dataset$expression
group_name <- dataset$sample_info$group_name
space <- reduce_dimensionality(expression, ndim=2)
traj <- infer_trajectory(space)
# set ntree to at least 1000!
gene_importances(expression, traj$time, num_permutations = 0, ntree = 1000)
#> # A tibble: 500 × 3
#>    gene    importance pvalue
#>    <chr>        <dbl> <lgl> 
#>  1 Gene261      0.198 NA    
#>  2 Gene160      0.188 NA    
#>  3 Gene251      0.177 NA    
#>  4 Gene431      0.164 NA    
#>  5 Gene276      0.153 NA    
#>  6 Gene463      0.152 NA    
#>  7 Gene38       0.144 NA    
#>  8 Gene329      0.137 NA    
#>  9 Gene314      0.134 NA    
#> 10 Gene122      0.127 NA    
#> # ℹ 490 more rows
```
