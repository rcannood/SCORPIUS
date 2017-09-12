context("Distance functions")

test_that("Correlation", {
  dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)

  same_dist1 <- correlation_distance(dataset$expression)
  expect_is( same_dist1, "matrix" )
  expect_true( diag(same_dist1) %>% map_lgl(~ all.equal(0, .)) %>% all )
  expect_gt( median(same_dist1), 0 )
  expect_equal(rownames(same_dist1), rownames(dataset$expression) )
  expect_equal(colnames(same_dist1), rownames(dataset$expression) )

  same_dist2 <- correlation_distance(dataset$expression, dataset$expression)

  expect_true( all.equal(1, cor(as.vector(same_dist1), as.vector(same_dist2))) )
})

test_that("Euclidean", {
  dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)

  same_dist1 <- euclidean_distance(dataset$expression)
  expect_is( same_dist1, "matrix" )
  expect_true( diag(same_dist1) %>% map_lgl(~ all.equal(0, .)) %>% all )
  expect_gt( median(same_dist1), 0 )
  expect_equal(rownames(same_dist1), rownames(dataset$expression) )
  expect_equal(colnames(same_dist1), rownames(dataset$expression) )

  same_dist2 <- euclidean_distance(dataset$expression, dataset$expression)

  expect_true( all.equal(1, cor(as.vector(same_dist1), as.vector(same_dist2))) )
})

