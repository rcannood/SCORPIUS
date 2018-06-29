context("Dimensionality reduction")

test_that("With generated data", {
  dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)

  # testing for 1 dimension
  space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim = 1, rescale = TRUE)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dataset$expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(1)) )

  ranges <- apply(space, 2, range)
  expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
  expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )

  # testing for 2 dimensions
  space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim = 2, rescale = TRUE)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dataset$expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(2)) )

  ranges <- apply(space, 2, range)
  expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
  expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )

  # testing for 5 dimensions without rescaling
  space <- reduce_dimensionality(dataset$expression, correlation_distance, ndim = 5, rescale = TRUE)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dataset$expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(5)) )
})
