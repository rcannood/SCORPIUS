context("KNN")

test_that("With generated data", {
  dataset <- generate_dataset(num_genes = 400, num_samples = 101, num_groups = 4)
  expression <- dataset$expression %>% scale_quantile(0)

  time <- seq(-1, 1, length.out = nrow(expression))

  dist <- dynutils::calculate_distance(expression)

  # testing without self loops allowed
  knn_out <- knn(dist, k = 5, self_loops = F)
  expect_equal(rownames(knn_out$indices), rownames(expression))
  expect_equal(rownames(knn_out$distances), rownames(expression))
  expect_equal(ncol(knn_out$indices), 5)
  expect_equal(ncol(knn_out$distances), 5)

  own_dists <- t(sapply(seq_len(nrow(dist)), function(i) dist[i, knn_out$indices[i,]]))
  expect_true(all(mapply(all.equal, as.vector(knn_out$distances), as.vector(own_dists))))

  expect_true(all(sapply(seq_len(nrow(dist)), function(i) {
    all(setdiff(order(dist[i,]), i)[1:5] == knn_out$indices[i,])
  })))

  # testing for different value of k
  knn_out <- knn(dist, k = 21, self_loops = F)
  expect_equal(rownames(knn_out$indices), rownames(expression))
  expect_equal(rownames(knn_out$distances), rownames(expression))
  expect_equal(ncol(knn_out$indices), 21)
  expect_equal(ncol(knn_out$distances), 21)

  own_dists <- t(sapply(seq_len(nrow(dist)), function(i) dist[i, knn_out$indices[i,]]))
  expect_true(all(mapply(all.equal, as.vector(knn_out$distances), as.vector(own_dists))))

  expect_true(all(sapply(seq_len(nrow(dist)), function(i) {
    all(setdiff(order(dist[i,]), i)[1:21] == knn_out$indices[i,])
  })))

  # testing with self loops allowed
  knn_out <- knn(dist, k = 5, self_loops = T)
  expect_equal(rownames(knn_out$indices), rownames(expression))
  expect_equal(rownames(knn_out$distances), rownames(expression))
  expect_equal(ncol(knn_out$indices), 5)
  expect_equal(ncol(knn_out$distances), 5)

  own_dists <- t(sapply(seq_len(nrow(dist)), function(i) dist[i, knn_out$indices[i,]]))
  expect_true(all(mapply(all.equal, as.vector(knn_out$distances), as.vector(own_dists))))

  expect_true(all(sapply(seq_len(nrow(dist)), function(i) {
    all(order(dist[i,])[1:5] == knn_out$indices[i,])
  })))

})
