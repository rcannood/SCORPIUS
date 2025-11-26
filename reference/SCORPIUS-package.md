# SCORPIUS: Trajectory inference from single-cell RNA sequencing data.

SCORPIUS orders single cells with regard to an implicit timeline, such
as cellular development or progression over time.

## Dimensionality Reduction functions

[`reduce_dimensionality`](rcannood.github.io/SCORPIUS/reference/reduce_dimensionality.md)

## Trajectory Inference functions

[`infer_trajectory`](rcannood.github.io/SCORPIUS/reference/infer_trajectory.md),
[`infer_initial_trajectory`](rcannood.github.io/SCORPIUS/reference/infer_initial_trajectory.md),
[`reverse_trajectory`](rcannood.github.io/SCORPIUS/reference/reverse_trajectory.md),
[`gene_importances`](rcannood.github.io/SCORPIUS/reference/gene_importances.md),
[`extract_modules`](rcannood.github.io/SCORPIUS/reference/extract_modules.md)

## Visualisation functions

[`draw_trajectory_plot`](rcannood.github.io/SCORPIUS/reference/draw_trajectory_plot.md),
[`draw_trajectory_heatmap`](rcannood.github.io/SCORPIUS/reference/draw_trajectory_heatmap.md)

## Datasets

[`generate_dataset`](rcannood.github.io/SCORPIUS/reference/generate_dataset.md),
[`ginhoux`](rcannood.github.io/SCORPIUS/reference/ginhoux.md)

## References

Cannoodt R. et al., SCORPIUS improves trajectory inference and
identifies novel modules in dendritic cell development, bioRxiv (Oct.,
2016). [doi:10.1101/079509](https://doi.org/10.1101/079509)
([PDF](https://www.biorxiv.org/content/biorxiv/early/2016/10/07/079509.full.pdf)).

## See also

Useful links:

- <https://github.com/rcannood/SCORPIUS>

- <http://rcannood.github.io/SCORPIUS/>

- Report bugs at <https://github.com/rcannood/SCORPIUS/issues>

## Author

**Maintainer**: Robrecht Cannoodt <rcannood@gmail.com>
([ORCID](https://orcid.org/0000-0003-3641-729X))

Other contributors:

- Wouter Saelens <wouter.saelens@ugent.be>
  ([ORCID](https://orcid.org/0000-0002-7114-6248)) \[contributor\]

## Examples

``` r
## Load dataset from Schlitzer et al., 2015
data("ginhoux")

## Reduce dimensionality and infer trajectory with SCORPIUS
space <- reduce_dimensionality(ginhoux$expression, "spearman")
traj <- infer_trajectory(space)

## Visualise
draw_trajectory_plot(
  space,
  path = traj$path,
  progression_group = ginhoux$sample_info$group_name
)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the SCORPIUS package.
#>   Please report the issue at <https://github.com/rcannood/SCORPIUS/issues>.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the SCORPIUS package.
#>   Please report the issue at <https://github.com/rcannood/SCORPIUS/issues>.
#> Ignoring unknown labels:
#> • fill : "Group"
```
