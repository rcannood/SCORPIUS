context("Dimensionality reduction")

test_that("With generated data", {
  dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)

  dist <- correlation_distance(dataset$expression)

  # testing for 1 dimension
  space <- reduce_dimensionality(dist, ndim = 1, rescale = T)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dist) )
  expect_equal( colnames(space), paste0("Comp", seq_len(1)) )

  ranges <- apply(space, 2, range)
  expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
  expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )

  # testing for 2 dimensions
  space <- reduce_dimensionality(dist, ndim = 2, rescale = T)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dist) )
  expect_equal( colnames(space), paste0("Comp", seq_len(2)) )

  ranges <- apply(space, 2, range)
  expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
  expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )

  # testing for 5 dimensions without rescaling
  space <- reduce_dimensionality(dist, ndim = 5, rescale = F)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dist) )
  expect_equal( colnames(space), paste0("Comp", seq_len(5)) )
})
