context("Scaling")

# generate some random data
num.samples <- 40
num.dims <- 10
data <- matrix(runif(num.samples * num.dims), nrow = num.samples)

test_that("Testing rescale.and.center", {
  data.sc <- rescale.and.center(data, center = 0, max.range = 1)
  ranges <- apply(data.sc, 2, range)

  expect_true( is.matrix(data.sc) )
  expect_equal(max(ranges[2,] - ranges[1,]), 1)
  expect_equal(ranges[1,] + ranges[2,], rep(0, num.dims))

  # try with a different center and max.range
  data.sc <- rescale.and.center(data, center = 10, max.range = 1000)
  ranges <- apply(data.sc, 2, range)

  expect_true( is.matrix(data.sc) )
  expect_equal(max(ranges[2,] - ranges[1,]), 1000)
  expect_equal(ranges[1,] + ranges[2,], rep(10, num.dims))
})
