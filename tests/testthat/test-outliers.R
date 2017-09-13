context("Outliers")

test_that("With generated data", {
  dataset <- generate_dataset(type = "polynomial", num_genes = 200, num_samples = 201, num_groups = 4)
  expression <- dataset$expression %>% scale_quantile(0)

  # setting 10 samples to random noise
  outl_ix <- sort(sample.int(nrow(expression), 10))
  expression[outl_ix, ] <- runif(10 * ncol(expression))
  dist <- correlation_distance(expression)
  dist[outl_ix,] <- dist[outl_ix,] * 1.2
  dist[,outl_ix] <- dist[,outl_ix] * 1.2
  dist[dist > 1] <- 1

  # testing outlierness with k = 5
  outl5 <- outlierness(dist, 5)
  expect_equal(length(outl5), nrow(expression))
  expect_equal(names(outl5), rownames(expression))
  expect_gt(quantile(outl5[outl_ix], .2), quantile(outl5[-outl_ix], .95))

  # testing outlierness with k = 10
  outl10 <- outlierness(dist, 10)
  expect_equal(length(outl10), nrow(expression))
  expect_equal(names(outl10), rownames(expression))
  expect_gt(quantile(outl10[outl_ix], .2), quantile(outl10[-outl_ix], .95))

  expect_true(all(outl5 < outl10))

  # testing whether the outlier filter works
  filt <- outlier_filter(dist)
  expect_equal(length(filt), nrow(expression))
  expect_equal(which(!filt), outl_ix)
})
