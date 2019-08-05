context("Sparsity support")

test_that("Large sparse matrix works", {
  expr <- Matrix::rsparsematrix(10000, 10000, .01)
  start <- Sys.time()
  space <- reduce_dimensionality(expr, dist = "spearman", ndim = 3)
  traj <- infer_trajectory(space)
  end <- Sys.time()
  expect_lte(as.numeric(difftime(end, start, units = "secs")), 20)
})


# library(SCORPIUS)
# profvis::profvis({
#   space <- reduce_dimensionality(expr, dist = "spearman", ndim = 3)
#   traj <- infer_trajectory(space)
# })
