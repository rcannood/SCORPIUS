context("Scaling")

# generate some random data
num_samples <- 40
num_dims <- 10
data <- matrix(runif(num_samples * num_dims), nrow = num_samples)

test_that("Testing rescale_and_center", {
  data_sc <- rescale_and_center(data, center = 0, max_range = 1)
  ranges <- apply(data_sc, 2, range)

  expect_is( data_sc, "matrix" )
  expect_equal(max(ranges[2,] - ranges[1,]), 1)
  expect_equal(ranges[1,] + ranges[2,], rep(0, num_dims))

  # try with a different center and max_range
  data_sc <- rescale_and_center(data, center = 10, max_range = 1000)
  ranges <- apply(data_sc, 2, range)

  expect_is( data_sc, "matrix" )
  expect_equal(max(ranges[2,] - ranges[1,]), 1000)
  expect_equal(ranges[1,] + ranges[2,], rep(10, num_dims))
})
