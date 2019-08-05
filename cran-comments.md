Added extra visualisation parameters and fixed a few minor bugs.

# Changelog

## Minor changes
 
 * Added extra customisation parameters to `draw_trajectory_plot()` and `draw_trajectory_heatmap()`.
 
## Optimisation

 * Fixed internal function `check_numeric_matrix()` such that it does not run for ages when applied to 
   a large sparse matrix.
   
 * Minor improvement in `infer_initial_trajectory()` when calculating the distance from points to 
   along candidate segments.

# Checks
## Test environments
* local Fedora 30 install, R 3.6.0
* OS X (on travis-ci), R 3.6.0
* ubuntu 16.04 (on travis-ci), R 3.6.0
* win-builder: release and devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.
