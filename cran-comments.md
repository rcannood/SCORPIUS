Submitting a new version of SCORPIUS due to breaking changes in R 3.6.0.

I've reuploaded the package because I forgot to remove the 'Date' field from the description, earlier.

# Changelog

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

# Checks
## Test environments
* local Fedora 28 install, R 3.6.0
* OS X (on travis-ci), R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder: no, because one of the dependencies is not built for windows yet.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.
