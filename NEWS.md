# SCORPIUS 1.0.8

* MINOR CHANGE: Allow adding row annotations to `draw_trajectory_heatmap()`.

# SCORPIUS 1.0.7 (2020-05-11)

Fix ahead of dplyr 1.0 release.

* MINOR CHANGE: substitute as.tbl_cube for reshape2::melt.

# SCORPIUS 1.0.6 (2020-03-16)

## Minor change

 * Resubmission of SCORPIUS. SCORPIUS was removed from CRAN because 
   dynwrap was removed from CRAN.
   
 * Added a vignette for using SCORPIUS to analyse Seurat data.

# SCORPIUS 1.0.5 (08-12-2019)

## Major change
 * Added a `ti_scorpius()` wrapper to SCORPIUS.
   
## Minor change
 * Use `RANN::nn2()` instead of own nearest neighbour functions. 
 
 * Remove deprecated functions.
 
 * Use `lmds` instead of `dyndimred`.

# SCORPIUS 1.0.4 (07-08-2019)

## Minor changes
 
 * Added extra customisation parameters to `draw_trajectory_plot()` and `draw_trajectory_heatmap()`.
 
## Optimisation

 * Fixed internal function `check_numeric_matrix()` such that it does not run for ages when applied to 
   a large sparse matrix.
   
 * Minor improvement in `infer_initial_trajectory()` when calculating the distance from points to 
   along candidate segments.

# SCORPIUS 1.0.3 (27-05-2019)

## Optimisation

 * `infer_trajectory()`: Use princurve's `approx_points` parameter, which greatly speeds up
   for trajectory inference for large number of samples.
   
## Major changes

 * Use dynutils' `calculate_distance()` instead of `correlation_distance()` and `euclidean_distance()`.
   
## Documentation

 * Vignettes were updated.

 * Added `cran-comments.md`.
 
 * Added recent news (`inst/NEWS.md`).
 
 * Added citation information (`inst/CITATION`).
 
 * Added support for sparsity in `extract_modules()` and `dimensionality_reduction()`.
 
## Minor changes

 * Use scaling functions from dynutils (`scale_minmax()`, `scale_quantile()`, `scale_uniform()`).
 
 * Expanded unit tests.
 
 * Renormalise the original ginhoux data using dynnormaliser and rerun all vignettes. 
 
 * Moved `cmdscale_withlandmarks()` to dyndimred.
 
## Bug fixes
 
 * `extract_modules()`: `smooth.spline()` now requires at least 4 unique values.
 
## Deprecation

 * Deprecated unused functions `evaluate_trajectory()` and `evaluate_dim_red()`.
   Use `dyneval::evaluate_ti_method()` instead.


# SCORPIUS 1.0.2 (2018-06-29)

 * MINOR CHANGE: Depend on dynutils for distance functions.
 
 * MAJOR CHANGE `reduce_dimensionality()`: Merge reduce_dimensionality_landmarked
   and reduce_dimensionality functions.

 * REMOVAL: Removed `outlier_filter()`; there are much better scRNA-seq preprocessing
   pipelines in existance by now.

# SCORPIUS 1.0.1 (2018-06-21)

 * MINOR CHANGE: Update for princurve 2.0.2.

 * REMOVAL: Removed deprecated functions `infer.trajectory()`, `reduce.dimensionality`, and so on.

 * DOCUMENTATION: Various fixes.

 * TESTING: Calculate code coverage on travis.
 
# SCORPIUS 1.0.0 (2017-09-15)

 * Initial release on CRAN.
