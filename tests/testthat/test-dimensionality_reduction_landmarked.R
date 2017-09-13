context("Dimensionality reduction landmarked")

test_that("With generated data", {
  dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 10000, num_groups = 4)

  space <- reduce_dimensionality_landmarked(dataset$expression, correlation_distance,
                                            num_landmarks = 200, ndim = 2, rescale = T)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dataset$expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(2)) )

  ranges <- apply(space, 2, range)
  expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
  expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )

  space <- reduce_dimensionality_landmarked(dataset$expression, correlation_distance,
                                            num_landmarks = 100, ndim = 3, rescale = F)

  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(dataset$expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(3)) )
})

test_that("Compare to normal dimred", {
  dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 2000, num_groups = 4)

  dist <- correlation_distance(dataset$expression)
  space <- reduce_dimensionality(dist)

  dist_normal <- as.vector(as.matrix(dist(space)))

  for (i in 1:10) {
    space_lm <- reduce_dimensionality_landmarked(dataset$expression, correlation_distance,
                                                 num_landmarks = 500)
    dist_lm <- as.vector(as.matrix(dist(space_lm)))
    expect_gt( cor(dist_normal, dist_lm), .3 )
  }
})
