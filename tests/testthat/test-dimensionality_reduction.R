context("Testing dimensionality_reduction.R")

dataset <- generate_dataset(type = "poly", num_genes = 500, num_samples = 1000, num_groups = 4)

check_space <- function(space, expression, ndim, rescale) {
  expect_is( space, c("data.frame", "matrix") )
  expect_equal( rownames(space), rownames(expression) )
  expect_equal( colnames(space), paste0("Comp", seq_len(ndim)) )

  if (rescale) {
    ranges <- apply(space, 2, range)
    expect_true( all.equal(1, max(apply(ranges, 2, diff))) )
    expect_true( all(apply(ranges, 2, function(x) all.equal(0, mean(x)))) )
  }
}

test_that("dimred with 1 dimension", {
  space <- reduce_dimensionality(dataset$expression, correlation_distance, landmark_method = "none", ndim = 1, rescale = TRUE)

  check_space(space, dataset$expression, ndim = 1, rescale = TRUE)
})

test_that("dimred with 2 dimensions", {
  space <- reduce_dimensionality(dataset$expression, correlation_distance, landmark_method = "none", ndim = 2, rescale = TRUE)

  check_space(space, dataset$expression, ndim = 2, rescale = TRUE)
})

test_that("dimred with 5 dimensions", {
  space <- reduce_dimensionality(dataset$expression, correlation_distance, landmark_method = "none", ndim = 5, rescale = FALSE)

  check_space(space, dataset$expression, ndim = 5, rescale = FALSE)
})

test_that("with landmarks", {
  space <- reduce_dimensionality(dataset$expression, correlation_distance, num_landmarks = 200, ndim = 2, rescale = TRUE)

  check_space(space, dataset$expression, ndim = 2, rescale = TRUE)

  space <- reduce_dimensionality(dataset$expression, correlation_distance, num_landmarks = 100, ndim = 5, rescale = FALSE)

  check_space(space, dataset$expression, ndim = 5, rescale = FALSE)

  space <- reduce_dimensionality(dataset$expression, correlation_distance, num_landmarks = 2000, ndim = 3, rescale = FALSE)

  check_space(space, dataset$expression, ndim = 3, rescale = FALSE)
})

test_that("compare landmarked with normal dimred", {
  space <- reduce_dimensionality(dataset$expression, correlation_distance, landmark_method = "none")

  dist_normal <- as.vector(as.matrix(dist(space)))

  for (i in 1:10) {
    space_lm <- reduce_dimensionality(dataset$expression, correlation_distance, num_landmarks = 500)
    dist_lm <- as.vector(as.matrix(dist(space_lm)))
    expect_gt( cor(dist_normal, dist_lm), .3 )
  }
})

test_that("fail gracefully", {
  expect_error(reduce_dimensionality(list(), correlation_distance), "must be a numeric matrix")

  expect_error(reduce_dimensionality(dataset$expression, 1), "must be a function")

  expect_error(reduce_dimensionality(dataset$expression, correlation_distance, ndim = 0), "must be a whole number")
  expect_error(reduce_dimensionality(dataset$expression, correlation_distance, ndim = 1.5), "must be a whole number")
  expect_error(reduce_dimensionality(dataset$expression, correlation_distance, ndim = nrow(dataset$expression) + 10), "must be a whole number")
})

test_that("landmark_selection works as expected", {
  landmarks <- landmark_selection(dataset$expression, correlation_distance, landmark_method = "naive", num_landmarks = 77)
  dist_normal <- correlation_distance(dataset$expression)

  expect_true(all(landmarks$ix_lm %in% seq_len(nrow(dataset$expression))))
  expect_gte(cor(as.vector(landmarks$dist_lm), as.vector(dist_normal[landmarks$ix_lm, landmarks$ix_lm, drop = FALSE])), .99)
  expect_gte(cor(as.vector(landmarks$dist_2lm), as.vector(dist_normal[landmarks$ix_lm, , drop = FALSE])), .99)

  landmarks <- landmark_selection(dataset$expression, euclidean_distance, landmark_method = "naive", num_landmarks = 77)
  dist_normal <- euclidean_distance(dataset$expression)

  expect_true(all(landmarks$ix_lm %in% seq_len(nrow(dataset$expression))))
  expect_gte(cor(as.vector(landmarks$dist_lm), as.vector(dist_normal[landmarks$ix_lm, landmarks$ix_lm, drop = FALSE])), .99)
  expect_gte(cor(as.vector(landmarks$dist_2lm), as.vector(dist_normal[landmarks$ix_lm, , drop = FALSE])), .99)

  expect_error(landmark_selection(dataset$expression, correlation_distance, landmark_method = "nioewnioew", num_landmarks = 100), "must be one of")
})
