# SCORPIUS 1.0.3 (unreleased)


## Optimisation

 * `infer_trajectory()`: Use princurve's `approx_points` parameter, which greatly speeds up
   for trajectory inference for large sample numbers.
   
## Documentation

 * Vignettes were updated.

 * Added `cran-comments.md`.
 
 * Added recent news (`inst/NEWS.md`).
 
 * Added citation information (`inst/CITATION`).
 
## Minor changes

 * Use dynutils' `calculate_distance()`.
 
 * Use scaling functions from dynutils (`scale_minmax()`, `scale_quantile()`, `scale_uniform()`).
 
 * Expanded unit tests.
 
 * Renormalise the original ginhoux data using dynnormaliser and rerun all vignettes. 
 
 * Moved `cmdscale_withlandmarks()` to dyndimred.
 
## Bug fixes
 
 * BUG FIX `cmdscale_withlandmarks()`: Fix colnames bug when some eigenvalues are equal to 0.
 
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
