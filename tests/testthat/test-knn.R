context("Gene importances")

test_that("With generated data", {
  dataset <- generate_dataset(type = "poly", num_genes = 400, num_samples = 101, num_groups = 4)
  expression <- dataset$expression %>% quant_scale(0)

  time <- seq(-1, 1, length.out = nrow(expression))

  dist <- correlation_distance(expression)

  knn_out <- knn(dist, k = 5, self_loops = F)
  expect_equal(rownames(knn_out$indices), rownames(expression))
  expect_equal(rownames(knn_out$distances), rownames(expression))
  expect_equal(ncol(knn_out$indices), 5)
  expect_equal(ncol(knn_out$distances), 5)

  own_dists <- t(sapply(seq_len(nrow(dist)), function(i) dist[i, knn_out$indices[i,]]))
  expect_true(all(mapply(all.equal, as.vector(knn_out$distances), as.vector(own_dists))))

  sapply(seq_len(nrow(dist)), function(i) {
    order(dist[i,])[1:5]

  })

})
