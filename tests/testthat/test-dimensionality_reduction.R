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

  expect_error(reduce_dimensionality(dataset$expression, correlation_distance, ndim = 0), "finite whole number")
  expect_error(reduce_dimensionality(dataset$expression, correlation_distance, ndim = 1.5), "finite whole number")
  expect_error(reduce_dimensionality(dataset$expression, correlation_distance, ndim = nrow(dataset$expression) + 10), "finite whole number")
})

test_that("landmark_selection works as expected", {
  landmarks <- landmark_selection(dataset$expression, correlation_distance, landmark_method = "naive", num_landmarks = 77)
  dist_normal <- correlation_distance(dataset$expression)

  expect_equal(length(landmarks$ix_lm), 77)
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

test_that("cmdscale_withlandmarks works as expected", {
  landmarks <- landmark_selection(dataset$expression, correlation_distance, landmark_method = "naive", num_landmarks = 77)
  num_samples <- 200
  num_landmarks <- 50
  sample_names <- paste0("Sample", seq_len(num_samples))
  landmark_names <- sample(sample_names, num_landmarks)

  dist_lm <- matrix(runif(num_landmarks * num_landmarks), nrow = num_landmarks, dimnames = list(landmark_names, landmark_names))
  dist_2lm <- matrix(runif(num_landmarks * num_samples), nrow = num_landmarks, dimnames = list(landmark_names, sample_names))

  cmdout <- cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = 4, rescale = TRUE)

  check_space(cmdout, t(dist_2lm), ndim = 4, rescale = TRUE)

  cmdout <- cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = 2, rescale = FALSE)

  check_space(cmdout, t(dist_2lm), ndim = 2, rescale = FALSE)

  dist_lm[1, 1] <- NA
  expect_error(cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = 2, rescale = FALSE), "finite numeric")
  dist_lm[1, 1] <- 0

  expect_error(cmdscale_withlandmarks(dist_lm[-1,], dist_2lm, ndim = 2, rescale = FALSE), "square matrix")

  expect_error(cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = 0), "finite whole number")
  expect_error(cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = 1.5), "finite whole number")
  expect_error(cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = nrow(dataset$expression) + 10), "finite whole number")

  num_samples <- 10
  sample_names <- paste0("Sample", seq_len(num_samples))
  dist_lm <- dist_2lm <- as.matrix(dist(matrix(rep(seq_len(num_samples), 2), ncol = 2)))
  expect_warning(cmdscale_withlandmarks(dist_lm, dist_2lm, ndim = num_samples, rescale = TRUE), "of the first [0-9]* eigenvalues are")
})
