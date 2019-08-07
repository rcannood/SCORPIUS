context("Sparsity support")

skip_on_cran()

test_that("Large sparse matrix works", {
  expr <- Matrix::rsparsematrix(10000, 10000, .01)
  start <- Sys.time()
  timing <- system.time({
    space <- reduce_dimensionality(expr, dist = "spearman", ndim = 3)
    traj <- infer_trajectory(space)
  })
  expect_lte(timing[["sys.self"]], 20)
})


# library(SCORPIUS)
# profvis::profvis({
#   space <- reduce_dimensionality(expr, dist = "spearman", ndim = 3)
#   traj <- infer_trajectory(space)
# })
