# SCORPIUS 1.0.3 (unreleased)

 * MINOR CHANGES: Cleaned up code style; mostly spacing issues, but now piping is used more.
 
 * MINOR CHANGE: Use scaling functions from dynutils (`scale_minmax()`, `scale_quantile()`, `scale_uniform()`)
 
 * OPTIMISATION `infer_trajectory()`: Added `approx_points` parameter, which greatly speeds up
   for trajectory inference for large sample numbers.

 * DOCUMENTATION: Vignettes were updated and regenerated.
 
 * DOCUMENTATION: Added `cran-comments.md` and `inst/NEWS.md`.
 
 * BUG FIX `cmdscale_withlandmarks()`: Fix colnames bug when some eigenvalues are equal to 0.
 
 * TESTING: Expanded unit tests.
 
 * DOCUMENTATION: Added citation information to package.
 
 * MINOR CHANGE: Renormalise the original ginhoux data using dynnormaliser and rerun all vignettes. 
 
 * MINOR CHANGE: Clean up `cmdscale_withlandmarks()` code.
 
 * DEPRECATION: Deprecated unused functions `evaluate_trajectory()` and `evaluate_dim_red()`.
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
