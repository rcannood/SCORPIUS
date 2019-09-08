context("Testing dimensionality_reduction.R")

dataset <- generate_dataset(num_genes = 500, num_samples = 1000, num_groups = 4)

check_space <- function(space, expression, ndim) {
  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(ndim)) )

  ranges <- apply(space, 2, range)
  expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
  expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )
}

test_that("dimred with 1 dimension", {
  space <- reduce_dimensionality(dataset$expression, "spearman", ndim = 1)

  check_space(space, dataset$expression, ndim = 1)
})

test_that("dimred with 2 dimensions", {
  space <- reduce_dimensionality(dataset$expression, "spearman", ndim = 2)

  check_space(space, dataset$expression, ndim = 2)
})

test_that("dimred with 5 dimensions", {
  space <- reduce_dimensionality(dataset$expression, "spearman", ndim = 5)

  check_space(space, dataset$expression, ndim = 5)
})

test_that("with landmarks", {
  space <- reduce_dimensionality(dataset$expression, "spearman", num_landmarks = 200, ndim = 2)

  check_space(space, dataset$expression, ndim = 2)

  space <- reduce_dimensionality(dataset$expression, "spearman", num_landmarks = 100, ndim = 5)

  check_space(space, dataset$expression, ndim = 5)

  space <- reduce_dimensionality(dataset$expression, "spearman", num_landmarks = 2000, ndim = 3)

  check_space(space, dataset$expression, ndim = 3)
})

test_that("compare landmarked with normal dimred", {
  space <- reduce_dimensionality(dataset$expression, "spearman")

  dist_normal <- as.vector(as.matrix(dist(space)))

  for (i in 1:10) {
    space_lm <- reduce_dimensionality(dataset$expression, "spearman", num_landmarks = 500)
    dist_lm <- as.vector(as.matrix(dist(space_lm)))
    expect_gt( cor(dist_normal, dist_lm), .3 )
  }
})

test_that("fail gracefully", {
  expect_error(reduce_dimensionality(list(), "spearman"), "must be a numeric matrix")

  expect_error(reduce_dimensionality(dataset$expression, "spearman", ndim = 0), "finite whole number")
  expect_error(reduce_dimensionality(dataset$expression, "spearman", ndim = 1.5), "finite whole number")
  expect_error(reduce_dimensionality(dataset$expression, "spearman", ndim = nrow(dataset$expression) + 10), "finite whole number")
})
